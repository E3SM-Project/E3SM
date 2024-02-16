   subroutine ocn_setup_nmlrec_run_modes(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      character (len=StrKIND) :: config_ocean_run_mode = 'forward'

      namelist /run_modes/ &
         config_ocean_run_mode
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, run_modes, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_char(dminfo, config_ocean_run_mode)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record run_modes')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_ocean_run_mode = '//mpas_log_escape_dollars(config_ocean_run_mode))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record run_modes.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_ocean_run_mode', config_ocean_run_mode)

   end subroutine ocn_setup_nmlrec_run_modes


   subroutine ocn_setup_nmlrec_time_management(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_do_restart = .false.
      character (len=StrKIND) :: config_restart_timestamp_name = 'Restart_timestamp'
      character (len=StrKIND) :: config_start_time = '0001-01-01_00:00:00'
      character (len=StrKIND) :: config_stop_time = 'none'
      character (len=StrKIND) :: config_run_duration = '0010_00:00:00'
      character (len=StrKIND) :: config_calendar_type = 'noleap'
      character (len=StrKIND) :: config_output_reference_time = '0001-01-01_00:00:00'

      namelist /time_management/ &
         config_do_restart, &
         config_restart_timestamp_name, &
         config_start_time, &
         config_stop_time, &
         config_run_duration, &
         config_calendar_type, &
         config_output_reference_time
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, time_management, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_do_restart)
         call mpas_dmpar_bcast_char(dminfo, config_restart_timestamp_name)
         call mpas_dmpar_bcast_char(dminfo, config_start_time)
         call mpas_dmpar_bcast_char(dminfo, config_stop_time)
         call mpas_dmpar_bcast_char(dminfo, config_run_duration)
         call mpas_dmpar_bcast_char(dminfo, config_calendar_type)
         call mpas_dmpar_bcast_char(dminfo, config_output_reference_time)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record time_management')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_do_restart = $l', logicArgs=(/config_do_restart/))
            call mpas_log_write('        config_restart_timestamp_name = '//mpas_log_escape_dollars(config_restart_timestamp_name))
            call mpas_log_write('        config_start_time = '//mpas_log_escape_dollars(config_start_time))
            call mpas_log_write('        config_stop_time = '//mpas_log_escape_dollars(config_stop_time))
            call mpas_log_write('        config_run_duration = '//mpas_log_escape_dollars(config_run_duration))
            call mpas_log_write('        config_calendar_type = '//mpas_log_escape_dollars(config_calendar_type))
            call mpas_log_write('        config_output_reference_time = '//mpas_log_escape_dollars(config_output_reference_time))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record time_management.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_do_restart', config_do_restart)
      call mpas_pool_add_config(configPool, 'config_restart_timestamp_name', config_restart_timestamp_name)
      call mpas_pool_add_config(configPool, 'config_start_time', config_start_time)
      call mpas_pool_add_config(configPool, 'config_stop_time', config_stop_time)
      call mpas_pool_add_config(configPool, 'config_run_duration', config_run_duration)
      call mpas_pool_add_config(configPool, 'config_calendar_type', config_calendar_type)
      call mpas_pool_add_config(configPool, 'config_output_reference_time', config_output_reference_time)

   end subroutine ocn_setup_nmlrec_time_management


   subroutine ocn_setup_nmlrec_io(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_write_output_on_startup = .true.
      integer :: config_pio_num_iotasks = 0
      integer :: config_pio_stride = 1

      namelist /io/ &
         config_write_output_on_startup, &
         config_pio_num_iotasks, &
         config_pio_stride
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, io, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_write_output_on_startup)
         call mpas_dmpar_bcast_int(dminfo, config_pio_num_iotasks)
         call mpas_dmpar_bcast_int(dminfo, config_pio_stride)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record io')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_write_output_on_startup = $l', logicArgs=(/config_write_output_on_startup/))
            call mpas_log_write('        config_pio_num_iotasks = $i', intArgs=(/config_pio_num_iotasks/))
            call mpas_log_write('        config_pio_stride = $i', intArgs=(/config_pio_stride/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record io.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_write_output_on_startup', config_write_output_on_startup)
      call mpas_pool_add_config(configPool, 'config_pio_num_iotasks', config_pio_num_iotasks)
      call mpas_pool_add_config(configPool, 'config_pio_stride', config_pio_stride)

   end subroutine ocn_setup_nmlrec_io


   subroutine ocn_setup_nmlrec_decomposition(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      integer :: config_num_halos = 3
      character (len=StrKIND) :: config_block_decomp_file_prefix = 'graph.info.part.'
      integer :: config_number_of_blocks = 0
      logical :: config_explicit_proc_decomp = .false.
      character (len=StrKIND) :: config_proc_decomp_file_prefix = 'graph.info.part.'

      namelist /decomposition/ &
         config_num_halos, &
         config_block_decomp_file_prefix, &
         config_number_of_blocks, &
         config_explicit_proc_decomp, &
         config_proc_decomp_file_prefix
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, decomposition, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_int(dminfo, config_num_halos)
         call mpas_dmpar_bcast_char(dminfo, config_block_decomp_file_prefix)
         call mpas_dmpar_bcast_int(dminfo, config_number_of_blocks)
         call mpas_dmpar_bcast_logical(dminfo, config_explicit_proc_decomp)
         call mpas_dmpar_bcast_char(dminfo, config_proc_decomp_file_prefix)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record decomposition')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_num_halos = $i', intArgs=(/config_num_halos/))
            call mpas_log_write('        config_block_decomp_file_prefix = '&
//''//mpas_log_escape_dollars(config_block_decomp_file_prefix))
            call mpas_log_write('        config_number_of_blocks = $i', intArgs=(/config_number_of_blocks/))
            call mpas_log_write('        config_explicit_proc_decomp = $l', logicArgs=(/config_explicit_proc_decomp/))
            call mpas_log_write('        config_proc_decomp_file_prefix = '&
//''//mpas_log_escape_dollars(config_proc_decomp_file_prefix))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record decomposition.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_num_halos', config_num_halos)
      call mpas_pool_add_config(configPool, 'config_block_decomp_file_prefix', config_block_decomp_file_prefix)
      call mpas_pool_add_config(configPool, 'config_number_of_blocks', config_number_of_blocks)
      call mpas_pool_add_config(configPool, 'config_explicit_proc_decomp', config_explicit_proc_decomp)
      call mpas_pool_add_config(configPool, 'config_proc_decomp_file_prefix', config_proc_decomp_file_prefix)

   end subroutine ocn_setup_nmlrec_decomposition


   subroutine ocn_setup_nmlrec_time_integration(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      character (len=StrKIND) :: config_dt = '00:05:00'
      character (len=StrKIND) :: config_time_integrator = 'split_explicit_ab2'
      integer :: config_number_of_time_levels = 2

      namelist /time_integration/ &
         config_dt, &
         config_time_integrator, &
         config_number_of_time_levels
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, time_integration, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_char(dminfo, config_dt)
         call mpas_dmpar_bcast_char(dminfo, config_time_integrator)
         call mpas_dmpar_bcast_int(dminfo, config_number_of_time_levels)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record time_integration')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_dt = '//mpas_log_escape_dollars(config_dt))
            call mpas_log_write('        config_time_integrator = '//mpas_log_escape_dollars(config_time_integrator))
            call mpas_log_write('        config_number_of_time_levels = $i', intArgs=(/config_number_of_time_levels/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record time_integration.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_dt', config_dt)
      call mpas_pool_add_config(configPool, 'config_time_integrator', config_time_integrator)
      call mpas_pool_add_config(configPool, 'config_number_of_time_levels', config_number_of_time_levels)

   end subroutine ocn_setup_nmlrec_time_integration


   subroutine ocn_setup_nmlrec_hmix(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_hmix_scaleWithMesh = .false.
      real (kind=RKIND) :: config_maxMeshDensity = -1.000000
      logical :: config_hmix_use_ref_cell_width = .false.
      real (kind=RKIND) :: config_hmix_ref_cell_width = 30000.000000
      real (kind=RKIND) :: config_apvm_scale_factor = 0.000000

      namelist /hmix/ &
         config_hmix_scaleWithMesh, &
         config_maxMeshDensity, &
         config_hmix_use_ref_cell_width, &
         config_hmix_ref_cell_width, &
         config_apvm_scale_factor
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, hmix, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_hmix_scaleWithMesh)
         call mpas_dmpar_bcast_real(dminfo, config_maxMeshDensity)
         call mpas_dmpar_bcast_logical(dminfo, config_hmix_use_ref_cell_width)
         call mpas_dmpar_bcast_real(dminfo, config_hmix_ref_cell_width)
         call mpas_dmpar_bcast_real(dminfo, config_apvm_scale_factor)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record hmix')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_hmix_scaleWithMesh = $l', logicArgs=(/config_hmix_scaleWithMesh/))
            call mpas_log_write('        config_maxMeshDensity = $r', realArgs=(/config_maxMeshDensity/))
            call mpas_log_write('        config_hmix_use_ref_cell_width = $l', logicArgs=(/config_hmix_use_ref_cell_width/))
            call mpas_log_write('        config_hmix_ref_cell_width = $r', realArgs=(/config_hmix_ref_cell_width/))
            call mpas_log_write('        config_apvm_scale_factor = $r', realArgs=(/config_apvm_scale_factor/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record hmix.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_hmix_scaleWithMesh', config_hmix_scaleWithMesh)
      call mpas_pool_add_config(configPool, 'config_maxMeshDensity', config_maxMeshDensity)
      call mpas_pool_add_config(configPool, 'config_hmix_use_ref_cell_width', config_hmix_use_ref_cell_width)
      call mpas_pool_add_config(configPool, 'config_hmix_ref_cell_width', config_hmix_ref_cell_width)
      call mpas_pool_add_config(configPool, 'config_apvm_scale_factor', config_apvm_scale_factor)

   end subroutine ocn_setup_nmlrec_hmix


   subroutine ocn_setup_nmlrec_hmix_del2(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_use_mom_del2 = .false.
      real (kind=RKIND) :: config_mom_del2 = 1000.000000
      logical :: config_use_tracer_del2 = .false.
      real (kind=RKIND) :: config_tracer_del2 = 10.000000

      namelist /hmix_del2/ &
         config_use_mom_del2, &
         config_mom_del2, &
         config_use_tracer_del2, &
         config_tracer_del2
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, hmix_del2, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_use_mom_del2)
         call mpas_dmpar_bcast_real(dminfo, config_mom_del2)
         call mpas_dmpar_bcast_logical(dminfo, config_use_tracer_del2)
         call mpas_dmpar_bcast_real(dminfo, config_tracer_del2)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record hmix_del2')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_use_mom_del2 = $l', logicArgs=(/config_use_mom_del2/))
            call mpas_log_write('        config_mom_del2 = $r', realArgs=(/config_mom_del2/))
            call mpas_log_write('        config_use_tracer_del2 = $l', logicArgs=(/config_use_tracer_del2/))
            call mpas_log_write('        config_tracer_del2 = $r', realArgs=(/config_tracer_del2/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record hmix_del2.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_use_mom_del2', config_use_mom_del2)
      call mpas_pool_add_config(configPool, 'config_mom_del2', config_mom_del2)
      call mpas_pool_add_config(configPool, 'config_use_tracer_del2', config_use_tracer_del2)
      call mpas_pool_add_config(configPool, 'config_tracer_del2', config_tracer_del2)

   end subroutine ocn_setup_nmlrec_hmix_del2


   subroutine ocn_setup_nmlrec_hmix_del4(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_use_mom_del4 = .false.
      real (kind=RKIND) :: config_mom_del4 = 120000000000.000000
      real (kind=RKIND) :: config_mom_del4_div_factor = 1.000000
      logical :: config_use_tracer_del4 = .false.
      real (kind=RKIND) :: config_tracer_del4 = 0.000000

      namelist /hmix_del4/ &
         config_use_mom_del4, &
         config_mom_del4, &
         config_mom_del4_div_factor, &
         config_use_tracer_del4, &
         config_tracer_del4
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, hmix_del4, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_use_mom_del4)
         call mpas_dmpar_bcast_real(dminfo, config_mom_del4)
         call mpas_dmpar_bcast_real(dminfo, config_mom_del4_div_factor)
         call mpas_dmpar_bcast_logical(dminfo, config_use_tracer_del4)
         call mpas_dmpar_bcast_real(dminfo, config_tracer_del4)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record hmix_del4')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_use_mom_del4 = $l', logicArgs=(/config_use_mom_del4/))
            call mpas_log_write('        config_mom_del4 = $r', realArgs=(/config_mom_del4/))
            call mpas_log_write('        config_mom_del4_div_factor = $r', realArgs=(/config_mom_del4_div_factor/))
            call mpas_log_write('        config_use_tracer_del4 = $l', logicArgs=(/config_use_tracer_del4/))
            call mpas_log_write('        config_tracer_del4 = $r', realArgs=(/config_tracer_del4/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record hmix_del4.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_use_mom_del4', config_use_mom_del4)
      call mpas_pool_add_config(configPool, 'config_mom_del4', config_mom_del4)
      call mpas_pool_add_config(configPool, 'config_mom_del4_div_factor', config_mom_del4_div_factor)
      call mpas_pool_add_config(configPool, 'config_use_tracer_del4', config_use_tracer_del4)
      call mpas_pool_add_config(configPool, 'config_tracer_del4', config_tracer_del4)

   end subroutine ocn_setup_nmlrec_hmix_del4


   subroutine ocn_setup_nmlrec_hmix_Leith(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_use_Leith_del2 = .false.
      real (kind=RKIND) :: config_Leith_parameter = 1.000000
      real (kind=RKIND) :: config_Leith_dx = 15000.000000
      real (kind=RKIND) :: config_Leith_visc2_max = 2500.000000

      namelist /hmix_Leith/ &
         config_use_Leith_del2, &
         config_Leith_parameter, &
         config_Leith_dx, &
         config_Leith_visc2_max
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, hmix_Leith, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_use_Leith_del2)
         call mpas_dmpar_bcast_real(dminfo, config_Leith_parameter)
         call mpas_dmpar_bcast_real(dminfo, config_Leith_dx)
         call mpas_dmpar_bcast_real(dminfo, config_Leith_visc2_max)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record hmix_Leith')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_use_Leith_del2 = $l', logicArgs=(/config_use_Leith_del2/))
            call mpas_log_write('        config_Leith_parameter = $r', realArgs=(/config_Leith_parameter/))
            call mpas_log_write('        config_Leith_dx = $r', realArgs=(/config_Leith_dx/))
            call mpas_log_write('        config_Leith_visc2_max = $r', realArgs=(/config_Leith_visc2_max/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record hmix_Leith.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_use_Leith_del2', config_use_Leith_del2)
      call mpas_pool_add_config(configPool, 'config_Leith_parameter', config_Leith_parameter)
      call mpas_pool_add_config(configPool, 'config_Leith_dx', config_Leith_dx)
      call mpas_pool_add_config(configPool, 'config_Leith_visc2_max', config_Leith_visc2_max)

   end subroutine ocn_setup_nmlrec_hmix_Leith


   subroutine ocn_setup_nmlrec_Redi_isopycnal_mixing(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_use_Redi = .false.
      character (len=StrKIND) :: config_Redi_closure = 'constant'
      real (kind=RKIND) :: config_Redi_constant_kappa = 600.000000
      real (kind=RKIND) :: config_Redi_maximum_slope = 0.300000
      logical :: config_Redi_use_slope_taper = .true.
      logical :: config_Redi_use_surface_taper = .true.
      logical :: config_Redi_limit_term1 = .true.
      logical :: config_Redi_use_quasi_monotone_limiter = .true.
      real (kind=RKIND) :: config_Redi_quasi_monotone_safety_factor = 0.900000
      integer :: config_Redi_min_layers_diag_terms = 6
      character (len=StrKIND) :: config_Redi_horizontal_taper = 'ramp'
      real (kind=RKIND) :: config_Redi_horizontal_ramp_min = 20000.000000
      real (kind=RKIND) :: config_Redi_horizontal_ramp_max = 30000.000000

      namelist /Redi_isopycnal_mixing/ &
         config_use_Redi, &
         config_Redi_closure, &
         config_Redi_constant_kappa, &
         config_Redi_maximum_slope, &
         config_Redi_use_slope_taper, &
         config_Redi_use_surface_taper, &
         config_Redi_limit_term1, &
         config_Redi_use_quasi_monotone_limiter, &
         config_Redi_quasi_monotone_safety_factor, &
         config_Redi_min_layers_diag_terms, &
         config_Redi_horizontal_taper, &
         config_Redi_horizontal_ramp_min, &
         config_Redi_horizontal_ramp_max
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, Redi_isopycnal_mixing, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_use_Redi)
         call mpas_dmpar_bcast_char(dminfo, config_Redi_closure)
         call mpas_dmpar_bcast_real(dminfo, config_Redi_constant_kappa)
         call mpas_dmpar_bcast_real(dminfo, config_Redi_maximum_slope)
         call mpas_dmpar_bcast_logical(dminfo, config_Redi_use_slope_taper)
         call mpas_dmpar_bcast_logical(dminfo, config_Redi_use_surface_taper)
         call mpas_dmpar_bcast_logical(dminfo, config_Redi_limit_term1)
         call mpas_dmpar_bcast_logical(dminfo, config_Redi_use_quasi_monotone_limiter)
         call mpas_dmpar_bcast_real(dminfo, config_Redi_quasi_monotone_safety_factor)
         call mpas_dmpar_bcast_int(dminfo, config_Redi_min_layers_diag_terms)
         call mpas_dmpar_bcast_char(dminfo, config_Redi_horizontal_taper)
         call mpas_dmpar_bcast_real(dminfo, config_Redi_horizontal_ramp_min)
         call mpas_dmpar_bcast_real(dminfo, config_Redi_horizontal_ramp_max)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record Redi_isopycnal_mixing')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_use_Redi = $l', logicArgs=(/config_use_Redi/))
            call mpas_log_write('        config_Redi_closure = '//mpas_log_escape_dollars(config_Redi_closure))
            call mpas_log_write('        config_Redi_constant_kappa = $r', realArgs=(/config_Redi_constant_kappa/))
            call mpas_log_write('        config_Redi_maximum_slope = $r', realArgs=(/config_Redi_maximum_slope/))
            call mpas_log_write('        config_Redi_use_slope_taper = $l', logicArgs=(/config_Redi_use_slope_taper/))
            call mpas_log_write('        config_Redi_use_surface_taper = $l', logicArgs=(/config_Redi_use_surface_taper/))
            call mpas_log_write('        config_Redi_limit_term1 = $l', logicArgs=(/config_Redi_limit_term1/))
            call mpas_log_write('        config_Redi_use_quasi_monotone_limiter = $l', &
logicArgs=(/config_Redi_use_quasi_monotone_limiter/))
            call mpas_log_write('        config_Redi_quasi_monotone_safety_factor = $r', &
realArgs=(/config_Redi_quasi_monotone_safety_factor/))
            call mpas_log_write('        config_Redi_min_layers_diag_terms = $i', intArgs=(/config_Redi_min_layers_diag_terms/))
            call mpas_log_write('        config_Redi_horizontal_taper = '//mpas_log_escape_dollars(config_Redi_horizontal_taper))
            call mpas_log_write('        config_Redi_horizontal_ramp_min = $r', realArgs=(/config_Redi_horizontal_ramp_min/))
            call mpas_log_write('        config_Redi_horizontal_ramp_max = $r', realArgs=(/config_Redi_horizontal_ramp_max/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record Redi_isopycnal_mixing.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_use_Redi', config_use_Redi)
      call mpas_pool_add_config(configPool, 'config_Redi_closure', config_Redi_closure)
      call mpas_pool_add_config(configPool, 'config_Redi_constant_kappa', config_Redi_constant_kappa)
      call mpas_pool_add_config(configPool, 'config_Redi_maximum_slope', config_Redi_maximum_slope)
      call mpas_pool_add_config(configPool, 'config_Redi_use_slope_taper', config_Redi_use_slope_taper)
      call mpas_pool_add_config(configPool, 'config_Redi_use_surface_taper', config_Redi_use_surface_taper)
      call mpas_pool_add_config(configPool, 'config_Redi_limit_term1', config_Redi_limit_term1)
      call mpas_pool_add_config(configPool, 'config_Redi_use_quasi_monotone_limiter', config_Redi_use_quasi_monotone_limiter)
      call mpas_pool_add_config(configPool, 'config_Redi_quasi_monotone_safety_factor', config_Redi_quasi_monotone_safety_factor)
      call mpas_pool_add_config(configPool, 'config_Redi_min_layers_diag_terms', config_Redi_min_layers_diag_terms)
      call mpas_pool_add_config(configPool, 'config_Redi_horizontal_taper', config_Redi_horizontal_taper)
      call mpas_pool_add_config(configPool, 'config_Redi_horizontal_ramp_min', config_Redi_horizontal_ramp_min)
      call mpas_pool_add_config(configPool, 'config_Redi_horizontal_ramp_max', config_Redi_horizontal_ramp_max)

   end subroutine ocn_setup_nmlrec_Redi_isopycnal_mixing


   subroutine ocn_setup_nmlrec_submesoscale_eddy_parameterization(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_submesoscale_enable = .false.
      real (kind=RKIND) :: config_submesoscale_tau = 172800.000000
      real (kind=RKIND) :: config_submesoscale_Ce = 0.060000
      real (kind=RKIND) :: config_submesoscale_Lfmin = 1000.000000
      real (kind=RKIND) :: config_submesoscale_ds_max = 100000.000000

      namelist /submesoscale_eddy_parameterization/ &
         config_submesoscale_enable, &
         config_submesoscale_tau, &
         config_submesoscale_Ce, &
         config_submesoscale_Lfmin, &
         config_submesoscale_ds_max
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, submesoscale_eddy_parameterization, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_submesoscale_enable)
         call mpas_dmpar_bcast_real(dminfo, config_submesoscale_tau)
         call mpas_dmpar_bcast_real(dminfo, config_submesoscale_Ce)
         call mpas_dmpar_bcast_real(dminfo, config_submesoscale_Lfmin)
         call mpas_dmpar_bcast_real(dminfo, config_submesoscale_ds_max)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record '&
//'submesoscale_eddy_parameterization')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_submesoscale_enable = $l', logicArgs=(/config_submesoscale_enable/))
            call mpas_log_write('        config_submesoscale_tau = $r', realArgs=(/config_submesoscale_tau/))
            call mpas_log_write('        config_submesoscale_Ce = $r', realArgs=(/config_submesoscale_Ce/))
            call mpas_log_write('        config_submesoscale_Lfmin = $r', realArgs=(/config_submesoscale_Lfmin/))
            call mpas_log_write('        config_submesoscale_ds_max = $r', realArgs=(/config_submesoscale_ds_max/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record submesoscale_eddy_parameterization.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_submesoscale_enable', config_submesoscale_enable)
      call mpas_pool_add_config(configPool, 'config_submesoscale_tau', config_submesoscale_tau)
      call mpas_pool_add_config(configPool, 'config_submesoscale_Ce', config_submesoscale_Ce)
      call mpas_pool_add_config(configPool, 'config_submesoscale_Lfmin', config_submesoscale_Lfmin)
      call mpas_pool_add_config(configPool, 'config_submesoscale_ds_max', config_submesoscale_ds_max)

   end subroutine ocn_setup_nmlrec_submesoscale_eddy_parameterization


   subroutine ocn_setup_nmlrec_GM_eddy_parameterization(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_use_GM = .false.
      character (len=StrKIND) :: config_GM_closure = 'EdenGreatbatch'
      real (kind=RKIND) :: config_GM_constant_kappa = 600.000000
      real (kind=RKIND) :: config_GM_constant_bclModeSpeed = 0.300000
      character (len=StrKIND) :: config_GM_minBclModeSpeed_method = 'constant'
      real (kind=RKIND) :: config_GM_spatially_variable_min_kappa = 300.000000
      real (kind=RKIND) :: config_GM_spatially_variable_max_kappa = 1800.000000
      real (kind=RKIND) :: config_GM_spatially_variable_baroclinic_mode = 1.000000
      real (kind=RKIND) :: config_GM_Visbeck_alpha = 0.005000
      real (kind=RKIND) :: config_GM_Visbeck_max_depth = 1000.000000
      real (kind=RKIND) :: config_GM_EG_riMin = 200.000000
      real (kind=RKIND) :: config_GM_EG_kappa_factor = 3.000000
      real (kind=RKIND) :: config_GM_EG_Rossby_factor = 2.000000
      real (kind=RKIND) :: config_GM_EG_Rhines_factor = 0.300000
      character (len=StrKIND) :: config_GM_horizontal_taper = 'ramp'
      real (kind=RKIND) :: config_GM_horizontal_ramp_min = 20000.000000
      real (kind=RKIND) :: config_GM_horizontal_ramp_max = 30000.000000
      real (kind=RKIND) :: config_GMRedi_Rossby_ramp_min = 0.500000
      real (kind=RKIND) :: config_GMRedi_Rossby_ramp_max = 3.000000

      namelist /GM_eddy_parameterization/ &
         config_use_GM, &
         config_GM_closure, &
         config_GM_constant_kappa, &
         config_GM_constant_bclModeSpeed, &
         config_GM_minBclModeSpeed_method, &
         config_GM_spatially_variable_min_kappa, &
         config_GM_spatially_variable_max_kappa, &
         config_GM_spatially_variable_baroclinic_mode, &
         config_GM_Visbeck_alpha, &
         config_GM_Visbeck_max_depth, &
         config_GM_EG_riMin, &
         config_GM_EG_kappa_factor, &
         config_GM_EG_Rossby_factor, &
         config_GM_EG_Rhines_factor, &
         config_GM_horizontal_taper, &
         config_GM_horizontal_ramp_min, &
         config_GM_horizontal_ramp_max, &
         config_GMRedi_Rossby_ramp_min, &
         config_GMRedi_Rossby_ramp_max
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, GM_eddy_parameterization, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_use_GM)
         call mpas_dmpar_bcast_char(dminfo, config_GM_closure)
         call mpas_dmpar_bcast_real(dminfo, config_GM_constant_kappa)
         call mpas_dmpar_bcast_real(dminfo, config_GM_constant_bclModeSpeed)
         call mpas_dmpar_bcast_char(dminfo, config_GM_minBclModeSpeed_method)
         call mpas_dmpar_bcast_real(dminfo, config_GM_spatially_variable_min_kappa)
         call mpas_dmpar_bcast_real(dminfo, config_GM_spatially_variable_max_kappa)
         call mpas_dmpar_bcast_real(dminfo, config_GM_spatially_variable_baroclinic_mode)
         call mpas_dmpar_bcast_real(dminfo, config_GM_Visbeck_alpha)
         call mpas_dmpar_bcast_real(dminfo, config_GM_Visbeck_max_depth)
         call mpas_dmpar_bcast_real(dminfo, config_GM_EG_riMin)
         call mpas_dmpar_bcast_real(dminfo, config_GM_EG_kappa_factor)
         call mpas_dmpar_bcast_real(dminfo, config_GM_EG_Rossby_factor)
         call mpas_dmpar_bcast_real(dminfo, config_GM_EG_Rhines_factor)
         call mpas_dmpar_bcast_char(dminfo, config_GM_horizontal_taper)
         call mpas_dmpar_bcast_real(dminfo, config_GM_horizontal_ramp_min)
         call mpas_dmpar_bcast_real(dminfo, config_GM_horizontal_ramp_max)
         call mpas_dmpar_bcast_real(dminfo, config_GMRedi_Rossby_ramp_min)
         call mpas_dmpar_bcast_real(dminfo, config_GMRedi_Rossby_ramp_max)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record GM_eddy_parameterization')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_use_GM = $l', logicArgs=(/config_use_GM/))
            call mpas_log_write('        config_GM_closure = '//mpas_log_escape_dollars(config_GM_closure))
            call mpas_log_write('        config_GM_constant_kappa = $r', realArgs=(/config_GM_constant_kappa/))
            call mpas_log_write('        config_GM_constant_bclModeSpeed = $r', realArgs=(/config_GM_constant_bclModeSpeed/))
            call mpas_log_write('        config_GM_minBclModeSpeed_method = '&
//''//mpas_log_escape_dollars(config_GM_minBclModeSpeed_method))
            call mpas_log_write('        config_GM_spatially_variable_min_kappa = $r', &
realArgs=(/config_GM_spatially_variable_min_kappa/))
            call mpas_log_write('        config_GM_spatially_variable_max_kappa = $r', &
realArgs=(/config_GM_spatially_variable_max_kappa/))
            call mpas_log_write('        config_GM_spatially_variable_baroclinic_mode = $r', &
realArgs=(/config_GM_spatially_variable_baroclinic_mode/))
            call mpas_log_write('        config_GM_Visbeck_alpha = $r', realArgs=(/config_GM_Visbeck_alpha/))
            call mpas_log_write('        config_GM_Visbeck_max_depth = $r', realArgs=(/config_GM_Visbeck_max_depth/))
            call mpas_log_write('        config_GM_EG_riMin = $r', realArgs=(/config_GM_EG_riMin/))
            call mpas_log_write('        config_GM_EG_kappa_factor = $r', realArgs=(/config_GM_EG_kappa_factor/))
            call mpas_log_write('        config_GM_EG_Rossby_factor = $r', realArgs=(/config_GM_EG_Rossby_factor/))
            call mpas_log_write('        config_GM_EG_Rhines_factor = $r', realArgs=(/config_GM_EG_Rhines_factor/))
            call mpas_log_write('        config_GM_horizontal_taper = '//mpas_log_escape_dollars(config_GM_horizontal_taper))
            call mpas_log_write('        config_GM_horizontal_ramp_min = $r', realArgs=(/config_GM_horizontal_ramp_min/))
            call mpas_log_write('        config_GM_horizontal_ramp_max = $r', realArgs=(/config_GM_horizontal_ramp_max/))
            call mpas_log_write('        config_GMRedi_Rossby_ramp_min = $r', realArgs=(/config_GMRedi_Rossby_ramp_min/))
            call mpas_log_write('        config_GMRedi_Rossby_ramp_max = $r', realArgs=(/config_GMRedi_Rossby_ramp_max/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record GM_eddy_parameterization.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_use_GM', config_use_GM)
      call mpas_pool_add_config(configPool, 'config_GM_closure', config_GM_closure)
      call mpas_pool_add_config(configPool, 'config_GM_constant_kappa', config_GM_constant_kappa)
      call mpas_pool_add_config(configPool, 'config_GM_constant_bclModeSpeed', config_GM_constant_bclModeSpeed)
      call mpas_pool_add_config(configPool, 'config_GM_minBclModeSpeed_method', config_GM_minBclModeSpeed_method)
      call mpas_pool_add_config(configPool, 'config_GM_spatially_variable_min_kappa', config_GM_spatially_variable_min_kappa)
      call mpas_pool_add_config(configPool, 'config_GM_spatially_variable_max_kappa', config_GM_spatially_variable_max_kappa)
      call mpas_pool_add_config(configPool, 'config_GM_spatially_variable_baroclinic_mode', &
config_GM_spatially_variable_baroclinic_mode)
      call mpas_pool_add_config(configPool, 'config_GM_Visbeck_alpha', config_GM_Visbeck_alpha)
      call mpas_pool_add_config(configPool, 'config_GM_Visbeck_max_depth', config_GM_Visbeck_max_depth)
      call mpas_pool_add_config(configPool, 'config_GM_EG_riMin', config_GM_EG_riMin)
      call mpas_pool_add_config(configPool, 'config_GM_EG_kappa_factor', config_GM_EG_kappa_factor)
      call mpas_pool_add_config(configPool, 'config_GM_EG_Rossby_factor', config_GM_EG_Rossby_factor)
      call mpas_pool_add_config(configPool, 'config_GM_EG_Rhines_factor', config_GM_EG_Rhines_factor)
      call mpas_pool_add_config(configPool, 'config_GM_horizontal_taper', config_GM_horizontal_taper)
      call mpas_pool_add_config(configPool, 'config_GM_horizontal_ramp_min', config_GM_horizontal_ramp_min)
      call mpas_pool_add_config(configPool, 'config_GM_horizontal_ramp_max', config_GM_horizontal_ramp_max)
      call mpas_pool_add_config(configPool, 'config_GMRedi_Rossby_ramp_min', config_GMRedi_Rossby_ramp_min)
      call mpas_pool_add_config(configPool, 'config_GMRedi_Rossby_ramp_max', config_GMRedi_Rossby_ramp_max)

   end subroutine ocn_setup_nmlrec_GM_eddy_parameterization


   subroutine ocn_setup_nmlrec_eddy_parameterization(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      real (kind=RKIND) :: config_eddyMLD_dens_threshold = 0.030000
      real (kind=RKIND) :: config_eddyMLD_reference_depth = 10.000000
      real (kind=RKIND) :: config_eddyMLD_reference_pressure = 100000.000000
      logical :: config_eddyMLD_use_old = .true.

      namelist /eddy_parameterization/ &
         config_eddyMLD_dens_threshold, &
         config_eddyMLD_reference_depth, &
         config_eddyMLD_reference_pressure, &
         config_eddyMLD_use_old
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, eddy_parameterization, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_real(dminfo, config_eddyMLD_dens_threshold)
         call mpas_dmpar_bcast_real(dminfo, config_eddyMLD_reference_depth)
         call mpas_dmpar_bcast_real(dminfo, config_eddyMLD_reference_pressure)
         call mpas_dmpar_bcast_logical(dminfo, config_eddyMLD_use_old)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record eddy_parameterization')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_eddyMLD_dens_threshold = $r', realArgs=(/config_eddyMLD_dens_threshold/))
            call mpas_log_write('        config_eddyMLD_reference_depth = $r', realArgs=(/config_eddyMLD_reference_depth/))
            call mpas_log_write('        config_eddyMLD_reference_pressure = $r', realArgs=(/config_eddyMLD_reference_pressure/))
            call mpas_log_write('        config_eddyMLD_use_old = $l', logicArgs=(/config_eddyMLD_use_old/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record eddy_parameterization.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_eddyMLD_dens_threshold', config_eddyMLD_dens_threshold)
      call mpas_pool_add_config(configPool, 'config_eddyMLD_reference_depth', config_eddyMLD_reference_depth)
      call mpas_pool_add_config(configPool, 'config_eddyMLD_reference_pressure', config_eddyMLD_reference_pressure)
      call mpas_pool_add_config(configPool, 'config_eddyMLD_use_old', config_eddyMLD_use_old)

   end subroutine ocn_setup_nmlrec_eddy_parameterization


   subroutine ocn_setup_nmlrec_cvmix(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_use_cvmix = .true.
      real (kind=RKIND) :: config_cvmix_prandtl_number = 1.000000
      character (len=StrKIND) :: config_cvmix_background_scheme = 'constant'
      real (kind=RKIND) :: config_cvmix_background_diffusion = 0.000010
      real (kind=RKIND) :: config_cvmix_background_viscosity = 0.000100
      real (kind=RKIND) :: config_cvmix_BryanLewis_bl1 = 0.000080
      real (kind=RKIND) :: config_cvmix_BryanLewis_bl2 = 0.000105
      real (kind=RKIND) :: config_cvmix_BryanLewis_transitionDepth = 2500.000000
      real (kind=RKIND) :: config_cvmix_BryanLewis_transitionWidth = 222.000000
      logical :: config_use_cvmix_convection = .false.
      real (kind=RKIND) :: config_cvmix_convective_diffusion = 1.000000
      real (kind=RKIND) :: config_cvmix_convective_viscosity = 1.000000
      logical :: config_cvmix_convective_basedOnBVF = .true.
      real (kind=RKIND) :: config_cvmix_convective_triggerBVF = 0.000000
      logical :: config_use_cvmix_shear = .false.
      integer :: config_cvmix_num_ri_smooth_loops = 2
      logical :: config_cvmix_use_BLD_smoothing = .true.
      character (len=StrKIND) :: config_cvmix_shear_mixing_scheme = 'PP'
      real (kind=RKIND) :: config_cvmix_shear_PP_nu_zero = 0.005000
      real (kind=RKIND) :: config_cvmix_shear_PP_alpha = 5.000000
      real (kind=RKIND) :: config_cvmix_shear_PP_exp = 2.000000
      real (kind=RKIND) :: config_cvmix_shear_KPP_nu_zero = 0.005000
      real (kind=RKIND) :: config_cvmix_shear_KPP_Ri_zero = 0.700000
      real (kind=RKIND) :: config_cvmix_shear_KPP_exp = 3.000000
      logical :: config_use_cvmix_tidal_mixing = .false.
      logical :: config_use_cvmix_double_diffusion = .false.
      logical :: config_use_cvmix_kpp = .false.
      logical :: config_use_cvmix_fixed_boundary_layer = .false.
      real (kind=RKIND) :: config_cvmix_kpp_boundary_layer_depth = 30.000000
      real (kind=RKIND) :: config_cvmix_kpp_criticalBulkRichardsonNumber = 0.250000
      character (len=StrKIND) :: config_cvmix_kpp_matching = 'SimpleShapes'
      logical :: config_cvmix_kpp_EkmanOBL = .false.
      logical :: config_cvmix_kpp_MonObOBL = .false.
      character (len=StrKIND) :: config_cvmix_kpp_interpolationOMLType = 'quadratic'
      real (kind=RKIND) :: config_cvmix_kpp_surface_layer_extent = 0.100000
      real (kind=RKIND) :: config_cvmix_kpp_surface_layer_averaging = 5.000000
      real (kind=RKIND) :: configure_cvmix_kpp_minimum_OBL_under_sea_ice = 10.000000
      real (kind=RKIND) :: config_cvmix_kpp_stop_OBL_search = 100.000000
      logical :: config_cvmix_kpp_use_enhanced_diff = .true.
      logical :: config_cvmix_kpp_nonlocal_with_implicit_mix = .false.
      logical :: config_cvmix_kpp_use_theory_wave = .false.
      character (len=StrKIND) :: config_cvmix_kpp_langmuir_mixing_opt = 'NONE'
      character (len=StrKIND) :: config_cvmix_kpp_langmuir_entrainment_opt = 'NONE'
      logical :: config_cvmix_kpp_use_active_wave = .false.

      namelist /cvmix/ &
         config_use_cvmix, &
         config_cvmix_prandtl_number, &
         config_cvmix_background_scheme, &
         config_cvmix_background_diffusion, &
         config_cvmix_background_viscosity, &
         config_cvmix_BryanLewis_bl1, &
         config_cvmix_BryanLewis_bl2, &
         config_cvmix_BryanLewis_transitionDepth, &
         config_cvmix_BryanLewis_transitionWidth, &
         config_use_cvmix_convection, &
         config_cvmix_convective_diffusion, &
         config_cvmix_convective_viscosity, &
         config_cvmix_convective_basedOnBVF, &
         config_cvmix_convective_triggerBVF, &
         config_use_cvmix_shear, &
         config_cvmix_num_ri_smooth_loops, &
         config_cvmix_use_BLD_smoothing, &
         config_cvmix_shear_mixing_scheme, &
         config_cvmix_shear_PP_nu_zero, &
         config_cvmix_shear_PP_alpha, &
         config_cvmix_shear_PP_exp, &
         config_cvmix_shear_KPP_nu_zero, &
         config_cvmix_shear_KPP_Ri_zero, &
         config_cvmix_shear_KPP_exp, &
         config_use_cvmix_tidal_mixing, &
         config_use_cvmix_double_diffusion, &
         config_use_cvmix_kpp, &
         config_use_cvmix_fixed_boundary_layer, &
         config_cvmix_kpp_boundary_layer_depth, &
         config_cvmix_kpp_criticalBulkRichardsonNumber, &
         config_cvmix_kpp_matching, &
         config_cvmix_kpp_EkmanOBL, &
         config_cvmix_kpp_MonObOBL, &
         config_cvmix_kpp_interpolationOMLType, &
         config_cvmix_kpp_surface_layer_extent, &
         config_cvmix_kpp_surface_layer_averaging, &
         configure_cvmix_kpp_minimum_OBL_under_sea_ice, &
         config_cvmix_kpp_stop_OBL_search, &
         config_cvmix_kpp_use_enhanced_diff, &
         config_cvmix_kpp_nonlocal_with_implicit_mix, &
         config_cvmix_kpp_use_theory_wave, &
         config_cvmix_kpp_langmuir_mixing_opt, &
         config_cvmix_kpp_langmuir_entrainment_opt, &
         config_cvmix_kpp_use_active_wave
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, cvmix, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_use_cvmix)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_prandtl_number)
         call mpas_dmpar_bcast_char(dminfo, config_cvmix_background_scheme)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_background_diffusion)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_background_viscosity)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_BryanLewis_bl1)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_BryanLewis_bl2)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_BryanLewis_transitionDepth)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_BryanLewis_transitionWidth)
         call mpas_dmpar_bcast_logical(dminfo, config_use_cvmix_convection)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_convective_diffusion)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_convective_viscosity)
         call mpas_dmpar_bcast_logical(dminfo, config_cvmix_convective_basedOnBVF)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_convective_triggerBVF)
         call mpas_dmpar_bcast_logical(dminfo, config_use_cvmix_shear)
         call mpas_dmpar_bcast_int(dminfo, config_cvmix_num_ri_smooth_loops)
         call mpas_dmpar_bcast_logical(dminfo, config_cvmix_use_BLD_smoothing)
         call mpas_dmpar_bcast_char(dminfo, config_cvmix_shear_mixing_scheme)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_shear_PP_nu_zero)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_shear_PP_alpha)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_shear_PP_exp)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_shear_KPP_nu_zero)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_shear_KPP_Ri_zero)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_shear_KPP_exp)
         call mpas_dmpar_bcast_logical(dminfo, config_use_cvmix_tidal_mixing)
         call mpas_dmpar_bcast_logical(dminfo, config_use_cvmix_double_diffusion)
         call mpas_dmpar_bcast_logical(dminfo, config_use_cvmix_kpp)
         call mpas_dmpar_bcast_logical(dminfo, config_use_cvmix_fixed_boundary_layer)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_kpp_boundary_layer_depth)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_kpp_criticalBulkRichardsonNumber)
         call mpas_dmpar_bcast_char(dminfo, config_cvmix_kpp_matching)
         call mpas_dmpar_bcast_logical(dminfo, config_cvmix_kpp_EkmanOBL)
         call mpas_dmpar_bcast_logical(dminfo, config_cvmix_kpp_MonObOBL)
         call mpas_dmpar_bcast_char(dminfo, config_cvmix_kpp_interpolationOMLType)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_kpp_surface_layer_extent)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_kpp_surface_layer_averaging)
         call mpas_dmpar_bcast_real(dminfo, configure_cvmix_kpp_minimum_OBL_under_sea_ice)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_kpp_stop_OBL_search)
         call mpas_dmpar_bcast_logical(dminfo, config_cvmix_kpp_use_enhanced_diff)
         call mpas_dmpar_bcast_logical(dminfo, config_cvmix_kpp_nonlocal_with_implicit_mix)
         call mpas_dmpar_bcast_logical(dminfo, config_cvmix_kpp_use_theory_wave)
         call mpas_dmpar_bcast_char(dminfo, config_cvmix_kpp_langmuir_mixing_opt)
         call mpas_dmpar_bcast_char(dminfo, config_cvmix_kpp_langmuir_entrainment_opt)
         call mpas_dmpar_bcast_logical(dminfo, config_cvmix_kpp_use_active_wave)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record cvmix')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_use_cvmix = $l', logicArgs=(/config_use_cvmix/))
            call mpas_log_write('        config_cvmix_prandtl_number = $r', realArgs=(/config_cvmix_prandtl_number/))
            call mpas_log_write('        config_cvmix_background_scheme = '&
//''//mpas_log_escape_dollars(config_cvmix_background_scheme))
            call mpas_log_write('        config_cvmix_background_diffusion = $r', realArgs=(/config_cvmix_background_diffusion/))
            call mpas_log_write('        config_cvmix_background_viscosity = $r', realArgs=(/config_cvmix_background_viscosity/))
            call mpas_log_write('        config_cvmix_BryanLewis_bl1 = $r', realArgs=(/config_cvmix_BryanLewis_bl1/))
            call mpas_log_write('        config_cvmix_BryanLewis_bl2 = $r', realArgs=(/config_cvmix_BryanLewis_bl2/))
            call mpas_log_write('        config_cvmix_BryanLewis_transitionDepth = $r', &
realArgs=(/config_cvmix_BryanLewis_transitionDepth/))
            call mpas_log_write('        config_cvmix_BryanLewis_transitionWidth = $r', &
realArgs=(/config_cvmix_BryanLewis_transitionWidth/))
            call mpas_log_write('        config_use_cvmix_convection = $l', logicArgs=(/config_use_cvmix_convection/))
            call mpas_log_write('        config_cvmix_convective_diffusion = $r', realArgs=(/config_cvmix_convective_diffusion/))
            call mpas_log_write('        config_cvmix_convective_viscosity = $r', realArgs=(/config_cvmix_convective_viscosity/))
            call mpas_log_write('        config_cvmix_convective_basedOnBVF = $l', &
logicArgs=(/config_cvmix_convective_basedOnBVF/))
            call mpas_log_write('        config_cvmix_convective_triggerBVF = $r', realArgs=(/config_cvmix_convective_triggerBVF/))
            call mpas_log_write('        config_use_cvmix_shear = $l', logicArgs=(/config_use_cvmix_shear/))
            call mpas_log_write('        config_cvmix_num_ri_smooth_loops = $i', intArgs=(/config_cvmix_num_ri_smooth_loops/))
            call mpas_log_write('        config_cvmix_use_BLD_smoothing = $l', logicArgs=(/config_cvmix_use_BLD_smoothing/))
            call mpas_log_write('        config_cvmix_shear_mixing_scheme = '&
//''//mpas_log_escape_dollars(config_cvmix_shear_mixing_scheme))
            call mpas_log_write('        config_cvmix_shear_PP_nu_zero = $r', realArgs=(/config_cvmix_shear_PP_nu_zero/))
            call mpas_log_write('        config_cvmix_shear_PP_alpha = $r', realArgs=(/config_cvmix_shear_PP_alpha/))
            call mpas_log_write('        config_cvmix_shear_PP_exp = $r', realArgs=(/config_cvmix_shear_PP_exp/))
            call mpas_log_write('        config_cvmix_shear_KPP_nu_zero = $r', realArgs=(/config_cvmix_shear_KPP_nu_zero/))
            call mpas_log_write('        config_cvmix_shear_KPP_Ri_zero = $r', realArgs=(/config_cvmix_shear_KPP_Ri_zero/))
            call mpas_log_write('        config_cvmix_shear_KPP_exp = $r', realArgs=(/config_cvmix_shear_KPP_exp/))
            call mpas_log_write('        config_use_cvmix_tidal_mixing = $l', logicArgs=(/config_use_cvmix_tidal_mixing/))
            call mpas_log_write('        config_use_cvmix_double_diffusion = $l', logicArgs=(/config_use_cvmix_double_diffusion/))
            call mpas_log_write('        config_use_cvmix_kpp = $l', logicArgs=(/config_use_cvmix_kpp/))
            call mpas_log_write('        config_use_cvmix_fixed_boundary_layer = $l', &
logicArgs=(/config_use_cvmix_fixed_boundary_layer/))
            call mpas_log_write('        config_cvmix_kpp_boundary_layer_depth = $r', &
realArgs=(/config_cvmix_kpp_boundary_layer_depth/))
            call mpas_log_write('        config_cvmix_kpp_criticalBulkRichardsonNumber = $r', &
realArgs=(/config_cvmix_kpp_criticalBulkRichardsonNumber/))
            call mpas_log_write('        config_cvmix_kpp_matching = '//mpas_log_escape_dollars(config_cvmix_kpp_matching))
            call mpas_log_write('        config_cvmix_kpp_EkmanOBL = $l', logicArgs=(/config_cvmix_kpp_EkmanOBL/))
            call mpas_log_write('        config_cvmix_kpp_MonObOBL = $l', logicArgs=(/config_cvmix_kpp_MonObOBL/))
            call mpas_log_write('        config_cvmix_kpp_interpolationOMLType = '&
//''//mpas_log_escape_dollars(config_cvmix_kpp_interpolationOMLType))
            call mpas_log_write('        config_cvmix_kpp_surface_layer_extent = $r', &
realArgs=(/config_cvmix_kpp_surface_layer_extent/))
            call mpas_log_write('        config_cvmix_kpp_surface_layer_averaging = $r', &
realArgs=(/config_cvmix_kpp_surface_layer_averaging/))
            call mpas_log_write('        configure_cvmix_kpp_minimum_OBL_under_sea_ice = $r', &
realArgs=(/configure_cvmix_kpp_minimum_OBL_under_sea_ice/))
            call mpas_log_write('        config_cvmix_kpp_stop_OBL_search = $r', realArgs=(/config_cvmix_kpp_stop_OBL_search/))
            call mpas_log_write('        config_cvmix_kpp_use_enhanced_diff = $l', &
logicArgs=(/config_cvmix_kpp_use_enhanced_diff/))
            call mpas_log_write('        config_cvmix_kpp_nonlocal_with_implicit_mix = $l', &
logicArgs=(/config_cvmix_kpp_nonlocal_with_implicit_mix/))
            call mpas_log_write('        config_cvmix_kpp_use_theory_wave = $l', logicArgs=(/config_cvmix_kpp_use_theory_wave/))
            call mpas_log_write('        config_cvmix_kpp_langmuir_mixing_opt = '&
//''//mpas_log_escape_dollars(config_cvmix_kpp_langmuir_mixing_opt))
            call mpas_log_write('        config_cvmix_kpp_langmuir_entrainment_opt = '&
//''//mpas_log_escape_dollars(config_cvmix_kpp_langmuir_entrainment_opt))
            call mpas_log_write('        config_cvmix_kpp_use_active_wave = $l', logicArgs=(/config_cvmix_kpp_use_active_wave/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record cvmix.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_use_cvmix', config_use_cvmix)
      call mpas_pool_add_config(configPool, 'config_cvmix_prandtl_number', config_cvmix_prandtl_number)
      call mpas_pool_add_config(configPool, 'config_cvmix_background_scheme', config_cvmix_background_scheme)
      call mpas_pool_add_config(configPool, 'config_cvmix_background_diffusion', config_cvmix_background_diffusion)
      call mpas_pool_add_config(configPool, 'config_cvmix_background_viscosity', config_cvmix_background_viscosity)
      call mpas_pool_add_config(configPool, 'config_cvmix_BryanLewis_bl1', config_cvmix_BryanLewis_bl1)
      call mpas_pool_add_config(configPool, 'config_cvmix_BryanLewis_bl2', config_cvmix_BryanLewis_bl2)
      call mpas_pool_add_config(configPool, 'config_cvmix_BryanLewis_transitionDepth', config_cvmix_BryanLewis_transitionDepth)
      call mpas_pool_add_config(configPool, 'config_cvmix_BryanLewis_transitionWidth', config_cvmix_BryanLewis_transitionWidth)
      call mpas_pool_add_config(configPool, 'config_use_cvmix_convection', config_use_cvmix_convection)
      call mpas_pool_add_config(configPool, 'config_cvmix_convective_diffusion', config_cvmix_convective_diffusion)
      call mpas_pool_add_config(configPool, 'config_cvmix_convective_viscosity', config_cvmix_convective_viscosity)
      call mpas_pool_add_config(configPool, 'config_cvmix_convective_basedOnBVF', config_cvmix_convective_basedOnBVF)
      call mpas_pool_add_config(configPool, 'config_cvmix_convective_triggerBVF', config_cvmix_convective_triggerBVF)
      call mpas_pool_add_config(configPool, 'config_use_cvmix_shear', config_use_cvmix_shear)
      call mpas_pool_add_config(configPool, 'config_cvmix_num_ri_smooth_loops', config_cvmix_num_ri_smooth_loops)
      call mpas_pool_add_config(configPool, 'config_cvmix_use_BLD_smoothing', config_cvmix_use_BLD_smoothing)
      call mpas_pool_add_config(configPool, 'config_cvmix_shear_mixing_scheme', config_cvmix_shear_mixing_scheme)
      call mpas_pool_add_config(configPool, 'config_cvmix_shear_PP_nu_zero', config_cvmix_shear_PP_nu_zero)
      call mpas_pool_add_config(configPool, 'config_cvmix_shear_PP_alpha', config_cvmix_shear_PP_alpha)
      call mpas_pool_add_config(configPool, 'config_cvmix_shear_PP_exp', config_cvmix_shear_PP_exp)
      call mpas_pool_add_config(configPool, 'config_cvmix_shear_KPP_nu_zero', config_cvmix_shear_KPP_nu_zero)
      call mpas_pool_add_config(configPool, 'config_cvmix_shear_KPP_Ri_zero', config_cvmix_shear_KPP_Ri_zero)
      call mpas_pool_add_config(configPool, 'config_cvmix_shear_KPP_exp', config_cvmix_shear_KPP_exp)
      call mpas_pool_add_config(configPool, 'config_use_cvmix_tidal_mixing', config_use_cvmix_tidal_mixing)
      call mpas_pool_add_config(configPool, 'config_use_cvmix_double_diffusion', config_use_cvmix_double_diffusion)
      call mpas_pool_add_config(configPool, 'config_use_cvmix_kpp', config_use_cvmix_kpp)
      call mpas_pool_add_config(configPool, 'config_use_cvmix_fixed_boundary_layer', config_use_cvmix_fixed_boundary_layer)
      call mpas_pool_add_config(configPool, 'config_cvmix_kpp_boundary_layer_depth', config_cvmix_kpp_boundary_layer_depth)
      call mpas_pool_add_config(configPool, 'config_cvmix_kpp_criticalBulkRichardsonNumber', &
config_cvmix_kpp_criticalBulkRichardsonNumber)
      call mpas_pool_add_config(configPool, 'config_cvmix_kpp_matching', config_cvmix_kpp_matching)
      call mpas_pool_add_config(configPool, 'config_cvmix_kpp_EkmanOBL', config_cvmix_kpp_EkmanOBL)
      call mpas_pool_add_config(configPool, 'config_cvmix_kpp_MonObOBL', config_cvmix_kpp_MonObOBL)
      call mpas_pool_add_config(configPool, 'config_cvmix_kpp_interpolationOMLType', config_cvmix_kpp_interpolationOMLType)
      call mpas_pool_add_config(configPool, 'config_cvmix_kpp_surface_layer_extent', config_cvmix_kpp_surface_layer_extent)
      call mpas_pool_add_config(configPool, 'config_cvmix_kpp_surface_layer_averaging', config_cvmix_kpp_surface_layer_averaging)
      call mpas_pool_add_config(configPool, 'configure_cvmix_kpp_minimum_OBL_under_sea_ice', &
configure_cvmix_kpp_minimum_OBL_under_sea_ice)
      call mpas_pool_add_config(configPool, 'config_cvmix_kpp_stop_OBL_search', config_cvmix_kpp_stop_OBL_search)
      call mpas_pool_add_config(configPool, 'config_cvmix_kpp_use_enhanced_diff', config_cvmix_kpp_use_enhanced_diff)
      call mpas_pool_add_config(configPool, 'config_cvmix_kpp_nonlocal_with_implicit_mix', &
config_cvmix_kpp_nonlocal_with_implicit_mix)
      call mpas_pool_add_config(configPool, 'config_cvmix_kpp_use_theory_wave', config_cvmix_kpp_use_theory_wave)
      call mpas_pool_add_config(configPool, 'config_cvmix_kpp_langmuir_mixing_opt', config_cvmix_kpp_langmuir_mixing_opt)
      call mpas_pool_add_config(configPool, 'config_cvmix_kpp_langmuir_entrainment_opt', config_cvmix_kpp_langmuir_entrainment_opt)
      call mpas_pool_add_config(configPool, 'config_cvmix_kpp_use_active_wave', config_cvmix_kpp_use_active_wave)

   end subroutine ocn_setup_nmlrec_cvmix


   subroutine ocn_setup_nmlrec_wave_coupling(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_use_active_wave = .false.
      integer :: config_n_stokes_drift_wavenumber_partitions = 6

      namelist /wave_coupling/ &
         config_use_active_wave, &
         config_n_stokes_drift_wavenumber_partitions
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, wave_coupling, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_use_active_wave)
         call mpas_dmpar_bcast_int(dminfo, config_n_stokes_drift_wavenumber_partitions)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record wave_coupling')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_use_active_wave = $l', logicArgs=(/config_use_active_wave/))
            call mpas_log_write('        config_n_stokes_drift_wavenumber_partitions = $i', &
intArgs=(/config_n_stokes_drift_wavenumber_partitions/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record wave_coupling.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_use_active_wave', config_use_active_wave)
      call mpas_pool_add_config(configPool, 'config_n_stokes_drift_wavenumber_partitions', &
config_n_stokes_drift_wavenumber_partitions)

   end subroutine ocn_setup_nmlrec_wave_coupling


   subroutine ocn_setup_nmlrec_gotm(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_use_gotm = .false.
      character (len=StrKIND) :: config_gotm_namelist_file = 'gotmturb.nml'
      real (kind=RKIND) :: config_gotm_constant_surface_roughness_length = 0.020000
      real (kind=RKIND) :: config_gotm_constant_bottom_roughness_length = 0.001500
      real (kind=RKIND) :: config_gotm_constant_bottom_drag_coeff = 0.001000

      namelist /gotm/ &
         config_use_gotm, &
         config_gotm_namelist_file, &
         config_gotm_constant_surface_roughness_length, &
         config_gotm_constant_bottom_roughness_length, &
         config_gotm_constant_bottom_drag_coeff
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, gotm, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_use_gotm)
         call mpas_dmpar_bcast_char(dminfo, config_gotm_namelist_file)
         call mpas_dmpar_bcast_real(dminfo, config_gotm_constant_surface_roughness_length)
         call mpas_dmpar_bcast_real(dminfo, config_gotm_constant_bottom_roughness_length)
         call mpas_dmpar_bcast_real(dminfo, config_gotm_constant_bottom_drag_coeff)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record gotm')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_use_gotm = $l', logicArgs=(/config_use_gotm/))
            call mpas_log_write('        config_gotm_namelist_file = '//mpas_log_escape_dollars(config_gotm_namelist_file))
            call mpas_log_write('        config_gotm_constant_surface_roughness_length = $r', &
realArgs=(/config_gotm_constant_surface_roughness_length/))
            call mpas_log_write('        config_gotm_constant_bottom_roughness_length = $r', &
realArgs=(/config_gotm_constant_bottom_roughness_length/))
            call mpas_log_write('        config_gotm_constant_bottom_drag_coeff = $r', &
realArgs=(/config_gotm_constant_bottom_drag_coeff/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record gotm.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_use_gotm', config_use_gotm)
      call mpas_pool_add_config(configPool, 'config_gotm_namelist_file', config_gotm_namelist_file)
      call mpas_pool_add_config(configPool, 'config_gotm_constant_surface_roughness_length', &
config_gotm_constant_surface_roughness_length)
      call mpas_pool_add_config(configPool, 'config_gotm_constant_bottom_roughness_length', &
config_gotm_constant_bottom_roughness_length)
      call mpas_pool_add_config(configPool, 'config_gotm_constant_bottom_drag_coeff', config_gotm_constant_bottom_drag_coeff)

   end subroutine ocn_setup_nmlrec_gotm


   subroutine ocn_setup_nmlrec_forcing(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_use_variable_drag = .false.
      logical :: config_use_bulk_wind_stress = .false.
      logical :: config_use_bulk_thickness_flux = .false.
      real (kind=RKIND) :: config_flux_attenuation_coefficient = 0.001000
      real (kind=RKIND) :: config_flux_attenuation_coefficient_runoff = 0.001000

      namelist /forcing/ &
         config_use_variable_drag, &
         config_use_bulk_wind_stress, &
         config_use_bulk_thickness_flux, &
         config_flux_attenuation_coefficient, &
         config_flux_attenuation_coefficient_runoff
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, forcing, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_use_variable_drag)
         call mpas_dmpar_bcast_logical(dminfo, config_use_bulk_wind_stress)
         call mpas_dmpar_bcast_logical(dminfo, config_use_bulk_thickness_flux)
         call mpas_dmpar_bcast_real(dminfo, config_flux_attenuation_coefficient)
         call mpas_dmpar_bcast_real(dminfo, config_flux_attenuation_coefficient_runoff)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record forcing')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_use_variable_drag = $l', logicArgs=(/config_use_variable_drag/))
            call mpas_log_write('        config_use_bulk_wind_stress = $l', logicArgs=(/config_use_bulk_wind_stress/))
            call mpas_log_write('        config_use_bulk_thickness_flux = $l', logicArgs=(/config_use_bulk_thickness_flux/))
            call mpas_log_write('        config_flux_attenuation_coefficient = $r', &
realArgs=(/config_flux_attenuation_coefficient/))
            call mpas_log_write('        config_flux_attenuation_coefficient_runoff = $r', &
realArgs=(/config_flux_attenuation_coefficient_runoff/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record forcing.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_use_variable_drag', config_use_variable_drag)
      call mpas_pool_add_config(configPool, 'config_use_bulk_wind_stress', config_use_bulk_wind_stress)
      call mpas_pool_add_config(configPool, 'config_use_bulk_thickness_flux', config_use_bulk_thickness_flux)
      call mpas_pool_add_config(configPool, 'config_flux_attenuation_coefficient', config_flux_attenuation_coefficient)
      call mpas_pool_add_config(configPool, 'config_flux_attenuation_coefficient_runoff', &
config_flux_attenuation_coefficient_runoff)

   end subroutine ocn_setup_nmlrec_forcing


   subroutine ocn_setup_nmlrec_time_varying_forcing(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_use_time_varying_atmospheric_forcing = .false.
      character (len=StrKIND) :: config_time_varying_atmospheric_forcing_type = 'WINDPRES'
      character (len=StrKIND) :: config_time_varying_atmospheric_forcing_start_time = '0001-01-01_00:00:00'
      character (len=StrKIND) :: config_time_varying_atmospheric_forcing_reference_time = '0001-01-01_00:00:00'
      character (len=StrKIND) :: config_time_varying_atmospheric_forcing_cycle_start = '0001-01-01_00:00:00'
      character (len=StrKIND) :: config_time_varying_atmospheric_forcing_cycle_duration = '2-00-00_00:00:00'
      character (len=StrKIND) :: config_time_varying_atmospheric_forcing_interval = '01:00:00'
      real (kind=RKIND) :: config_time_varying_atmospheric_forcing_ramp = 10.000000
      real (kind=RKIND) :: config_time_varying_atmospheric_forcing_ramp_delay = 0.000000
      logical :: config_use_time_varying_land_ice_forcing = .false.
      character (len=StrKIND) :: config_time_varying_land_ice_forcing_start_time = '0001-01-01_00:00:00'
      character (len=StrKIND) :: config_time_varying_land_ice_forcing_reference_time = '0001-01-01_00:00:00'
      character (len=StrKIND) :: config_time_varying_land_ice_forcing_cycle_start = '0001-01-01_00:00:00'
      character (len=StrKIND) :: config_time_varying_land_ice_forcing_cycle_duration = '2-00-00_00:00:00'
      character (len=StrKIND) :: config_time_varying_land_ice_forcing_interval = '01:00:00'

      namelist /time_varying_forcing/ &
         config_use_time_varying_atmospheric_forcing, &
         config_time_varying_atmospheric_forcing_type, &
         config_time_varying_atmospheric_forcing_start_time, &
         config_time_varying_atmospheric_forcing_reference_time, &
         config_time_varying_atmospheric_forcing_cycle_start, &
         config_time_varying_atmospheric_forcing_cycle_duration, &
         config_time_varying_atmospheric_forcing_interval, &
         config_time_varying_atmospheric_forcing_ramp, &
         config_time_varying_atmospheric_forcing_ramp_delay, &
         config_use_time_varying_land_ice_forcing, &
         config_time_varying_land_ice_forcing_start_time, &
         config_time_varying_land_ice_forcing_reference_time, &
         config_time_varying_land_ice_forcing_cycle_start, &
         config_time_varying_land_ice_forcing_cycle_duration, &
         config_time_varying_land_ice_forcing_interval
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, time_varying_forcing, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_use_time_varying_atmospheric_forcing)
         call mpas_dmpar_bcast_char(dminfo, config_time_varying_atmospheric_forcing_type)
         call mpas_dmpar_bcast_char(dminfo, config_time_varying_atmospheric_forcing_start_time)
         call mpas_dmpar_bcast_char(dminfo, config_time_varying_atmospheric_forcing_reference_time)
         call mpas_dmpar_bcast_char(dminfo, config_time_varying_atmospheric_forcing_cycle_start)
         call mpas_dmpar_bcast_char(dminfo, config_time_varying_atmospheric_forcing_cycle_duration)
         call mpas_dmpar_bcast_char(dminfo, config_time_varying_atmospheric_forcing_interval)
         call mpas_dmpar_bcast_real(dminfo, config_time_varying_atmospheric_forcing_ramp)
         call mpas_dmpar_bcast_real(dminfo, config_time_varying_atmospheric_forcing_ramp_delay)
         call mpas_dmpar_bcast_logical(dminfo, config_use_time_varying_land_ice_forcing)
         call mpas_dmpar_bcast_char(dminfo, config_time_varying_land_ice_forcing_start_time)
         call mpas_dmpar_bcast_char(dminfo, config_time_varying_land_ice_forcing_reference_time)
         call mpas_dmpar_bcast_char(dminfo, config_time_varying_land_ice_forcing_cycle_start)
         call mpas_dmpar_bcast_char(dminfo, config_time_varying_land_ice_forcing_cycle_duration)
         call mpas_dmpar_bcast_char(dminfo, config_time_varying_land_ice_forcing_interval)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record time_varying_forcing')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_use_time_varying_atmospheric_forcing = $l', &
logicArgs=(/config_use_time_varying_atmospheric_forcing/))
            call mpas_log_write('        config_time_varying_atmospheric_forcing_type = '&
//''//mpas_log_escape_dollars(config_time_varying_atmospheric_forcing_type))
            call mpas_log_write('        config_time_varying_atmospheric_forcing_start_time = '&
//''//mpas_log_escape_dollars(config_time_varying_atmospheric_forcing_start_time))
            call mpas_log_write('        config_time_varying_atmospheric_forcing_reference_time = '&
//''//mpas_log_escape_dollars(config_time_varying_atmospheric_forcing_reference_time))
            call mpas_log_write('        config_time_varying_atmospheric_forcing_cycle_start = '&
//''//mpas_log_escape_dollars(config_time_varying_atmospheric_forcing_cycle_start))
            call mpas_log_write('        config_time_varying_atmospheric_forcing_cycle_duration = '&
//''//mpas_log_escape_dollars(config_time_varying_atmospheric_forcing_cycle_duration))
            call mpas_log_write('        config_time_varying_atmospheric_forcing_interval = '&
//''//mpas_log_escape_dollars(config_time_varying_atmospheric_forcing_interval))
            call mpas_log_write('        config_time_varying_atmospheric_forcing_ramp = $r', &
realArgs=(/config_time_varying_atmospheric_forcing_ramp/))
            call mpas_log_write('        config_time_varying_atmospheric_forcing_ramp_delay = $r', &
realArgs=(/config_time_varying_atmospheric_forcing_ramp_delay/))
            call mpas_log_write('        config_use_time_varying_land_ice_forcing = $l', &
logicArgs=(/config_use_time_varying_land_ice_forcing/))
            call mpas_log_write('        config_time_varying_land_ice_forcing_start_time = '&
//''//mpas_log_escape_dollars(config_time_varying_land_ice_forcing_start_time))
            call mpas_log_write('        config_time_varying_land_ice_forcing_reference_time = '&
//''//mpas_log_escape_dollars(config_time_varying_land_ice_forcing_reference_time))
            call mpas_log_write('        config_time_varying_land_ice_forcing_cycle_start = '&
//''//mpas_log_escape_dollars(config_time_varying_land_ice_forcing_cycle_start))
            call mpas_log_write('        config_time_varying_land_ice_forcing_cycle_duration = '&
//''//mpas_log_escape_dollars(config_time_varying_land_ice_forcing_cycle_duration))
            call mpas_log_write('        config_time_varying_land_ice_forcing_interval = '&
//''//mpas_log_escape_dollars(config_time_varying_land_ice_forcing_interval))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record time_varying_forcing.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_use_time_varying_atmospheric_forcing', &
config_use_time_varying_atmospheric_forcing)
      call mpas_pool_add_config(configPool, 'config_time_varying_atmospheric_forcing_type', &
config_time_varying_atmospheric_forcing_type)
      call mpas_pool_add_config(configPool, 'config_time_varying_atmospheric_forcing_start_time', &
config_time_varying_atmospheric_forcing_start_time)
      call mpas_pool_add_config(configPool, 'config_time_varying_atmospheric_forcing_reference_time', &
config_time_varying_atmospheric_forcing_reference_time)
      call mpas_pool_add_config(configPool, 'config_time_varying_atmospheric_forcing_cycle_start', &
config_time_varying_atmospheric_forcing_cycle_start)
      call mpas_pool_add_config(configPool, 'config_time_varying_atmospheric_forcing_cycle_duration', &
config_time_varying_atmospheric_forcing_cycle_duration)
      call mpas_pool_add_config(configPool, 'config_time_varying_atmospheric_forcing_interval', &
config_time_varying_atmospheric_forcing_interval)
      call mpas_pool_add_config(configPool, 'config_time_varying_atmospheric_forcing_ramp', &
config_time_varying_atmospheric_forcing_ramp)
      call mpas_pool_add_config(configPool, 'config_time_varying_atmospheric_forcing_ramp_delay', &
config_time_varying_atmospheric_forcing_ramp_delay)
      call mpas_pool_add_config(configPool, 'config_use_time_varying_land_ice_forcing', config_use_time_varying_land_ice_forcing)
      call mpas_pool_add_config(configPool, 'config_time_varying_land_ice_forcing_start_time', &
config_time_varying_land_ice_forcing_start_time)
      call mpas_pool_add_config(configPool, 'config_time_varying_land_ice_forcing_reference_time', &
config_time_varying_land_ice_forcing_reference_time)
      call mpas_pool_add_config(configPool, 'config_time_varying_land_ice_forcing_cycle_start', &
config_time_varying_land_ice_forcing_cycle_start)
      call mpas_pool_add_config(configPool, 'config_time_varying_land_ice_forcing_cycle_duration', &
config_time_varying_land_ice_forcing_cycle_duration)
      call mpas_pool_add_config(configPool, 'config_time_varying_land_ice_forcing_interval', &
config_time_varying_land_ice_forcing_interval)

   end subroutine ocn_setup_nmlrec_time_varying_forcing


   subroutine ocn_setup_nmlrec_coupling(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      real (kind=RKIND) :: config_ssh_grad_relax_timescale = 0.000000
      logical :: config_remove_AIS_coupler_runoff = .false.

      namelist /coupling/ &
         config_ssh_grad_relax_timescale, &
         config_remove_AIS_coupler_runoff
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, coupling, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_real(dminfo, config_ssh_grad_relax_timescale)
         call mpas_dmpar_bcast_logical(dminfo, config_remove_AIS_coupler_runoff)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record coupling')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_ssh_grad_relax_timescale = $r', realArgs=(/config_ssh_grad_relax_timescale/))
            call mpas_log_write('        config_remove_AIS_coupler_runoff = $l', logicArgs=(/config_remove_AIS_coupler_runoff/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record coupling.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_ssh_grad_relax_timescale', config_ssh_grad_relax_timescale)
      call mpas_pool_add_config(configPool, 'config_remove_AIS_coupler_runoff', config_remove_AIS_coupler_runoff)

   end subroutine ocn_setup_nmlrec_coupling


   subroutine ocn_setup_nmlrec_shortwaveRadiation(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      character (len=StrKIND) :: config_sw_absorption_type = 'none'
      integer :: config_jerlov_water_type = 3
      real (kind=RKIND) :: config_surface_buoyancy_depth = 1.000000
      logical :: config_enable_shortwave_energy_fixer = .false.

      namelist /shortwaveRadiation/ &
         config_sw_absorption_type, &
         config_jerlov_water_type, &
         config_surface_buoyancy_depth, &
         config_enable_shortwave_energy_fixer
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, shortwaveRadiation, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_char(dminfo, config_sw_absorption_type)
         call mpas_dmpar_bcast_int(dminfo, config_jerlov_water_type)
         call mpas_dmpar_bcast_real(dminfo, config_surface_buoyancy_depth)
         call mpas_dmpar_bcast_logical(dminfo, config_enable_shortwave_energy_fixer)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record shortwaveRadiation')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_sw_absorption_type = '//mpas_log_escape_dollars(config_sw_absorption_type))
            call mpas_log_write('        config_jerlov_water_type = $i', intArgs=(/config_jerlov_water_type/))
            call mpas_log_write('        config_surface_buoyancy_depth = $r', realArgs=(/config_surface_buoyancy_depth/))
            call mpas_log_write('        config_enable_shortwave_energy_fixer = $l', &
logicArgs=(/config_enable_shortwave_energy_fixer/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record shortwaveRadiation.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_sw_absorption_type', config_sw_absorption_type)
      call mpas_pool_add_config(configPool, 'config_jerlov_water_type', config_jerlov_water_type)
      call mpas_pool_add_config(configPool, 'config_surface_buoyancy_depth', config_surface_buoyancy_depth)
      call mpas_pool_add_config(configPool, 'config_enable_shortwave_energy_fixer', config_enable_shortwave_energy_fixer)

   end subroutine ocn_setup_nmlrec_shortwaveRadiation


   subroutine ocn_setup_nmlrec_tidal_forcing(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_use_tidal_forcing = .false.
      real (kind=RKIND) :: config_use_tidal_forcing_tau = 10000.000000
      character (len=StrKIND) :: config_tidal_forcing_type = 'off'
      character (len=StrKIND) :: config_tidal_forcing_model = 'off'
      real (kind=RKIND) :: config_tidal_forcing_monochromatic_amp = 2.000000
      real (kind=RKIND) :: config_tidal_forcing_monochromatic_period = 0.500000
      real (kind=RKIND) :: config_tidal_forcing_monochromatic_phaseLag = 0.000000
      real (kind=RKIND) :: config_tidal_forcing_monochromatic_baseline = 0.000000

      namelist /tidal_forcing/ &
         config_use_tidal_forcing, &
         config_use_tidal_forcing_tau, &
         config_tidal_forcing_type, &
         config_tidal_forcing_model, &
         config_tidal_forcing_monochromatic_amp, &
         config_tidal_forcing_monochromatic_period, &
         config_tidal_forcing_monochromatic_phaseLag, &
         config_tidal_forcing_monochromatic_baseline
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, tidal_forcing, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_use_tidal_forcing)
         call mpas_dmpar_bcast_real(dminfo, config_use_tidal_forcing_tau)
         call mpas_dmpar_bcast_char(dminfo, config_tidal_forcing_type)
         call mpas_dmpar_bcast_char(dminfo, config_tidal_forcing_model)
         call mpas_dmpar_bcast_real(dminfo, config_tidal_forcing_monochromatic_amp)
         call mpas_dmpar_bcast_real(dminfo, config_tidal_forcing_monochromatic_period)
         call mpas_dmpar_bcast_real(dminfo, config_tidal_forcing_monochromatic_phaseLag)
         call mpas_dmpar_bcast_real(dminfo, config_tidal_forcing_monochromatic_baseline)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record tidal_forcing')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_use_tidal_forcing = $l', logicArgs=(/config_use_tidal_forcing/))
            call mpas_log_write('        config_use_tidal_forcing_tau = $r', realArgs=(/config_use_tidal_forcing_tau/))
            call mpas_log_write('        config_tidal_forcing_type = '//mpas_log_escape_dollars(config_tidal_forcing_type))
            call mpas_log_write('        config_tidal_forcing_model = '//mpas_log_escape_dollars(config_tidal_forcing_model))
            call mpas_log_write('        config_tidal_forcing_monochromatic_amp = $r', &
realArgs=(/config_tidal_forcing_monochromatic_amp/))
            call mpas_log_write('        config_tidal_forcing_monochromatic_period = $r', &
realArgs=(/config_tidal_forcing_monochromatic_period/))
            call mpas_log_write('        config_tidal_forcing_monochromatic_phaseLag = $r', &
realArgs=(/config_tidal_forcing_monochromatic_phaseLag/))
            call mpas_log_write('        config_tidal_forcing_monochromatic_baseline = $r', &
realArgs=(/config_tidal_forcing_monochromatic_baseline/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record tidal_forcing.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_use_tidal_forcing', config_use_tidal_forcing)
      call mpas_pool_add_config(configPool, 'config_use_tidal_forcing_tau', config_use_tidal_forcing_tau)
      call mpas_pool_add_config(configPool, 'config_tidal_forcing_type', config_tidal_forcing_type)
      call mpas_pool_add_config(configPool, 'config_tidal_forcing_model', config_tidal_forcing_model)
      call mpas_pool_add_config(configPool, 'config_tidal_forcing_monochromatic_amp', config_tidal_forcing_monochromatic_amp)
      call mpas_pool_add_config(configPool, 'config_tidal_forcing_monochromatic_period', config_tidal_forcing_monochromatic_period)
      call mpas_pool_add_config(configPool, 'config_tidal_forcing_monochromatic_phaseLag', &
config_tidal_forcing_monochromatic_phaseLag)
      call mpas_pool_add_config(configPool, 'config_tidal_forcing_monochromatic_baseline', &
config_tidal_forcing_monochromatic_baseline)

   end subroutine ocn_setup_nmlrec_tidal_forcing


   subroutine ocn_setup_nmlrec_self_attraction_loading(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_use_self_attraction_loading = .false.
      real (kind=RKIND) :: config_self_attraction_loading_smoothing_width = 1.000000
      character (len=StrKIND) :: config_mpas_to_grid_weights_file = 'mpas_to_grid.nc'
      character (len=StrKIND) :: config_grid_to_mpas_weights_file = 'grid_to_mpas.nc'
      character (len=StrKIND) :: config_self_attraction_loading_compute_interval = '0000-00-00_00:30:00'
      integer :: config_nLatitude = 128
      integer :: config_nLongitude = 256
      logical :: config_use_parallel_self_attraction_loading = .false.
      integer :: config_parallel_self_attraction_loading_order = 10
      integer :: config_parallel_self_attraction_loading_n_cells_per_block = 600
      logical :: config_parallel_self_attraction_loading_bfb = .false.

      namelist /self_attraction_loading/ &
         config_use_self_attraction_loading, &
         config_self_attraction_loading_smoothing_width, &
         config_mpas_to_grid_weights_file, &
         config_grid_to_mpas_weights_file, &
         config_self_attraction_loading_compute_interval, &
         config_nLatitude, &
         config_nLongitude, &
         config_use_parallel_self_attraction_loading, &
         config_parallel_self_attraction_loading_order, &
         config_parallel_self_attraction_loading_n_cells_per_block, &
         config_parallel_self_attraction_loading_bfb
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, self_attraction_loading, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_use_self_attraction_loading)
         call mpas_dmpar_bcast_real(dminfo, config_self_attraction_loading_smoothing_width)
         call mpas_dmpar_bcast_char(dminfo, config_mpas_to_grid_weights_file)
         call mpas_dmpar_bcast_char(dminfo, config_grid_to_mpas_weights_file)
         call mpas_dmpar_bcast_char(dminfo, config_self_attraction_loading_compute_interval)
         call mpas_dmpar_bcast_int(dminfo, config_nLatitude)
         call mpas_dmpar_bcast_int(dminfo, config_nLongitude)
         call mpas_dmpar_bcast_logical(dminfo, config_use_parallel_self_attraction_loading)
         call mpas_dmpar_bcast_int(dminfo, config_parallel_self_attraction_loading_order)
         call mpas_dmpar_bcast_int(dminfo, config_parallel_self_attraction_loading_n_cells_per_block)
         call mpas_dmpar_bcast_logical(dminfo, config_parallel_self_attraction_loading_bfb)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record self_attraction_loading')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_use_self_attraction_loading = $l', &
logicArgs=(/config_use_self_attraction_loading/))
            call mpas_log_write('        config_self_attraction_loading_smoothing_width = $r', &
realArgs=(/config_self_attraction_loading_smoothing_width/))
            call mpas_log_write('        config_mpas_to_grid_weights_file = '&
//''//mpas_log_escape_dollars(config_mpas_to_grid_weights_file))
            call mpas_log_write('        config_grid_to_mpas_weights_file = '&
//''//mpas_log_escape_dollars(config_grid_to_mpas_weights_file))
            call mpas_log_write('        config_self_attraction_loading_compute_interval = '&
//''//mpas_log_escape_dollars(config_self_attraction_loading_compute_interval))
            call mpas_log_write('        config_nLatitude = $i', intArgs=(/config_nLatitude/))
            call mpas_log_write('        config_nLongitude = $i', intArgs=(/config_nLongitude/))
            call mpas_log_write('        config_use_parallel_self_attraction_loading = $l', &
logicArgs=(/config_use_parallel_self_attraction_loading/))
            call mpas_log_write('        config_parallel_self_attraction_loading_order = $i', &
intArgs=(/config_parallel_self_attraction_loading_order/))
            call mpas_log_write('        config_parallel_self_attraction_loading_n_cells_per_block = $i', &
intArgs=(/config_parallel_self_attraction_loading_n_cells_per_block/))
            call mpas_log_write('        config_parallel_self_attraction_loading_bfb = $l', &
logicArgs=(/config_parallel_self_attraction_loading_bfb/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record self_attraction_loading.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_use_self_attraction_loading', config_use_self_attraction_loading)
      call mpas_pool_add_config(configPool, 'config_self_attraction_loading_smoothing_width', &
config_self_attraction_loading_smoothing_width)
      call mpas_pool_add_config(configPool, 'config_mpas_to_grid_weights_file', config_mpas_to_grid_weights_file)
      call mpas_pool_add_config(configPool, 'config_grid_to_mpas_weights_file', config_grid_to_mpas_weights_file)
      call mpas_pool_add_config(configPool, 'config_self_attraction_loading_compute_interval', &
config_self_attraction_loading_compute_interval)
      call mpas_pool_add_config(configPool, 'config_nLatitude', config_nLatitude)
      call mpas_pool_add_config(configPool, 'config_nLongitude', config_nLongitude)
      call mpas_pool_add_config(configPool, 'config_use_parallel_self_attraction_loading', &
config_use_parallel_self_attraction_loading)
      call mpas_pool_add_config(configPool, 'config_parallel_self_attraction_loading_order', &
config_parallel_self_attraction_loading_order)
      call mpas_pool_add_config(configPool, 'config_parallel_self_attraction_loading_n_cells_per_block', &
config_parallel_self_attraction_loading_n_cells_per_block)
      call mpas_pool_add_config(configPool, 'config_parallel_self_attraction_loading_bfb', &
config_parallel_self_attraction_loading_bfb)

   end subroutine ocn_setup_nmlrec_self_attraction_loading


   subroutine ocn_setup_nmlrec_tidal_potential_forcing(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_use_tidal_potential_forcing = .false.
      character (len=StrKIND) :: config_tidal_potential_reference_time = '0001-01-01_00:00:00'
      logical :: config_use_tidal_potential_forcing_M2 = .true.
      logical :: config_use_tidal_potential_forcing_S2 = .true.
      logical :: config_use_tidal_potential_forcing_N2 = .true.
      logical :: config_use_tidal_potential_forcing_K2 = .true.
      logical :: config_use_tidal_potential_forcing_K1 = .true.
      logical :: config_use_tidal_potential_forcing_O1 = .true.
      logical :: config_use_tidal_potential_forcing_Q1 = .true.
      logical :: config_use_tidal_potential_forcing_P1 = .true.
      real (kind=RKIND) :: config_tidal_potential_ramp = 10.000000
      real (kind=RKIND) :: config_self_attraction_and_loading_beta = 0.090000

      namelist /tidal_potential_forcing/ &
         config_use_tidal_potential_forcing, &
         config_tidal_potential_reference_time, &
         config_use_tidal_potential_forcing_M2, &
         config_use_tidal_potential_forcing_S2, &
         config_use_tidal_potential_forcing_N2, &
         config_use_tidal_potential_forcing_K2, &
         config_use_tidal_potential_forcing_K1, &
         config_use_tidal_potential_forcing_O1, &
         config_use_tidal_potential_forcing_Q1, &
         config_use_tidal_potential_forcing_P1, &
         config_tidal_potential_ramp, &
         config_self_attraction_and_loading_beta
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, tidal_potential_forcing, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_use_tidal_potential_forcing)
         call mpas_dmpar_bcast_char(dminfo, config_tidal_potential_reference_time)
         call mpas_dmpar_bcast_logical(dminfo, config_use_tidal_potential_forcing_M2)
         call mpas_dmpar_bcast_logical(dminfo, config_use_tidal_potential_forcing_S2)
         call mpas_dmpar_bcast_logical(dminfo, config_use_tidal_potential_forcing_N2)
         call mpas_dmpar_bcast_logical(dminfo, config_use_tidal_potential_forcing_K2)
         call mpas_dmpar_bcast_logical(dminfo, config_use_tidal_potential_forcing_K1)
         call mpas_dmpar_bcast_logical(dminfo, config_use_tidal_potential_forcing_O1)
         call mpas_dmpar_bcast_logical(dminfo, config_use_tidal_potential_forcing_Q1)
         call mpas_dmpar_bcast_logical(dminfo, config_use_tidal_potential_forcing_P1)
         call mpas_dmpar_bcast_real(dminfo, config_tidal_potential_ramp)
         call mpas_dmpar_bcast_real(dminfo, config_self_attraction_and_loading_beta)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record tidal_potential_forcing')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_use_tidal_potential_forcing = $l', &
logicArgs=(/config_use_tidal_potential_forcing/))
            call mpas_log_write('        config_tidal_potential_reference_time = '&
//''//mpas_log_escape_dollars(config_tidal_potential_reference_time))
            call mpas_log_write('        config_use_tidal_potential_forcing_M2 = $l', &
logicArgs=(/config_use_tidal_potential_forcing_M2/))
            call mpas_log_write('        config_use_tidal_potential_forcing_S2 = $l', &
logicArgs=(/config_use_tidal_potential_forcing_S2/))
            call mpas_log_write('        config_use_tidal_potential_forcing_N2 = $l', &
logicArgs=(/config_use_tidal_potential_forcing_N2/))
            call mpas_log_write('        config_use_tidal_potential_forcing_K2 = $l', &
logicArgs=(/config_use_tidal_potential_forcing_K2/))
            call mpas_log_write('        config_use_tidal_potential_forcing_K1 = $l', &
logicArgs=(/config_use_tidal_potential_forcing_K1/))
            call mpas_log_write('        config_use_tidal_potential_forcing_O1 = $l', &
logicArgs=(/config_use_tidal_potential_forcing_O1/))
            call mpas_log_write('        config_use_tidal_potential_forcing_Q1 = $l', &
logicArgs=(/config_use_tidal_potential_forcing_Q1/))
            call mpas_log_write('        config_use_tidal_potential_forcing_P1 = $l', &
logicArgs=(/config_use_tidal_potential_forcing_P1/))
            call mpas_log_write('        config_tidal_potential_ramp = $r', realArgs=(/config_tidal_potential_ramp/))
            call mpas_log_write('        config_self_attraction_and_loading_beta = $r', &
realArgs=(/config_self_attraction_and_loading_beta/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record tidal_potential_forcing.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_use_tidal_potential_forcing', config_use_tidal_potential_forcing)
      call mpas_pool_add_config(configPool, 'config_tidal_potential_reference_time', config_tidal_potential_reference_time)
      call mpas_pool_add_config(configPool, 'config_use_tidal_potential_forcing_M2', config_use_tidal_potential_forcing_M2)
      call mpas_pool_add_config(configPool, 'config_use_tidal_potential_forcing_S2', config_use_tidal_potential_forcing_S2)
      call mpas_pool_add_config(configPool, 'config_use_tidal_potential_forcing_N2', config_use_tidal_potential_forcing_N2)
      call mpas_pool_add_config(configPool, 'config_use_tidal_potential_forcing_K2', config_use_tidal_potential_forcing_K2)
      call mpas_pool_add_config(configPool, 'config_use_tidal_potential_forcing_K1', config_use_tidal_potential_forcing_K1)
      call mpas_pool_add_config(configPool, 'config_use_tidal_potential_forcing_O1', config_use_tidal_potential_forcing_O1)
      call mpas_pool_add_config(configPool, 'config_use_tidal_potential_forcing_Q1', config_use_tidal_potential_forcing_Q1)
      call mpas_pool_add_config(configPool, 'config_use_tidal_potential_forcing_P1', config_use_tidal_potential_forcing_P1)
      call mpas_pool_add_config(configPool, 'config_tidal_potential_ramp', config_tidal_potential_ramp)
      call mpas_pool_add_config(configPool, 'config_self_attraction_and_loading_beta', config_self_attraction_and_loading_beta)

   end subroutine ocn_setup_nmlrec_tidal_potential_forcing


   subroutine ocn_setup_nmlrec_frazil_ice(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_use_frazil_ice_formation = .false.
      logical :: config_frazil_in_open_ocean = .true.
      logical :: config_frazil_under_land_ice = .true.
      real (kind=RKIND) :: config_frazil_heat_of_fusion = 334000.000000
      real (kind=RKIND) :: config_frazil_ice_density = 1000.000000
      real (kind=RKIND) :: config_frazil_fractional_thickness_limit = 0.100000
      real (kind=RKIND) :: config_specific_heat_sea_water = 3985.000000
      real (kind=RKIND) :: config_frazil_maximum_depth = 100.000000
      real (kind=RKIND) :: config_frazil_sea_ice_reference_salinity = 4.000000
      real (kind=RKIND) :: config_frazil_land_ice_reference_salinity = 0.000000
      real (kind=RKIND) :: config_frazil_maximum_freezing_temperature = 0.000000
      logical :: config_frazil_use_surface_pressure = .true.

      namelist /frazil_ice/ &
         config_use_frazil_ice_formation, &
         config_frazil_in_open_ocean, &
         config_frazil_under_land_ice, &
         config_frazil_heat_of_fusion, &
         config_frazil_ice_density, &
         config_frazil_fractional_thickness_limit, &
         config_specific_heat_sea_water, &
         config_frazil_maximum_depth, &
         config_frazil_sea_ice_reference_salinity, &
         config_frazil_land_ice_reference_salinity, &
         config_frazil_maximum_freezing_temperature, &
         config_frazil_use_surface_pressure
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, frazil_ice, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_use_frazil_ice_formation)
         call mpas_dmpar_bcast_logical(dminfo, config_frazil_in_open_ocean)
         call mpas_dmpar_bcast_logical(dminfo, config_frazil_under_land_ice)
         call mpas_dmpar_bcast_real(dminfo, config_frazil_heat_of_fusion)
         call mpas_dmpar_bcast_real(dminfo, config_frazil_ice_density)
         call mpas_dmpar_bcast_real(dminfo, config_frazil_fractional_thickness_limit)
         call mpas_dmpar_bcast_real(dminfo, config_specific_heat_sea_water)
         call mpas_dmpar_bcast_real(dminfo, config_frazil_maximum_depth)
         call mpas_dmpar_bcast_real(dminfo, config_frazil_sea_ice_reference_salinity)
         call mpas_dmpar_bcast_real(dminfo, config_frazil_land_ice_reference_salinity)
         call mpas_dmpar_bcast_real(dminfo, config_frazil_maximum_freezing_temperature)
         call mpas_dmpar_bcast_logical(dminfo, config_frazil_use_surface_pressure)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record frazil_ice')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_use_frazil_ice_formation = $l', logicArgs=(/config_use_frazil_ice_formation/))
            call mpas_log_write('        config_frazil_in_open_ocean = $l', logicArgs=(/config_frazil_in_open_ocean/))
            call mpas_log_write('        config_frazil_under_land_ice = $l', logicArgs=(/config_frazil_under_land_ice/))
            call mpas_log_write('        config_frazil_heat_of_fusion = $r', realArgs=(/config_frazil_heat_of_fusion/))
            call mpas_log_write('        config_frazil_ice_density = $r', realArgs=(/config_frazil_ice_density/))
            call mpas_log_write('        config_frazil_fractional_thickness_limit = $r', &
realArgs=(/config_frazil_fractional_thickness_limit/))
            call mpas_log_write('        config_specific_heat_sea_water = $r', realArgs=(/config_specific_heat_sea_water/))
            call mpas_log_write('        config_frazil_maximum_depth = $r', realArgs=(/config_frazil_maximum_depth/))
            call mpas_log_write('        config_frazil_sea_ice_reference_salinity = $r', &
realArgs=(/config_frazil_sea_ice_reference_salinity/))
            call mpas_log_write('        config_frazil_land_ice_reference_salinity = $r', &
realArgs=(/config_frazil_land_ice_reference_salinity/))
            call mpas_log_write('        config_frazil_maximum_freezing_temperature = $r', &
realArgs=(/config_frazil_maximum_freezing_temperature/))
            call mpas_log_write('        config_frazil_use_surface_pressure = $l', &
logicArgs=(/config_frazil_use_surface_pressure/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record frazil_ice.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_use_frazil_ice_formation', config_use_frazil_ice_formation)
      call mpas_pool_add_config(configPool, 'config_frazil_in_open_ocean', config_frazil_in_open_ocean)
      call mpas_pool_add_config(configPool, 'config_frazil_under_land_ice', config_frazil_under_land_ice)
      call mpas_pool_add_config(configPool, 'config_frazil_heat_of_fusion', config_frazil_heat_of_fusion)
      call mpas_pool_add_config(configPool, 'config_frazil_ice_density', config_frazil_ice_density)
      call mpas_pool_add_config(configPool, 'config_frazil_fractional_thickness_limit', config_frazil_fractional_thickness_limit)
      call mpas_pool_add_config(configPool, 'config_specific_heat_sea_water', config_specific_heat_sea_water)
      call mpas_pool_add_config(configPool, 'config_frazil_maximum_depth', config_frazil_maximum_depth)
      call mpas_pool_add_config(configPool, 'config_frazil_sea_ice_reference_salinity', config_frazil_sea_ice_reference_salinity)
      call mpas_pool_add_config(configPool, 'config_frazil_land_ice_reference_salinity', config_frazil_land_ice_reference_salinity)
      call mpas_pool_add_config(configPool, 'config_frazil_maximum_freezing_temperature', &
config_frazil_maximum_freezing_temperature)
      call mpas_pool_add_config(configPool, 'config_frazil_use_surface_pressure', config_frazil_use_surface_pressure)

   end subroutine ocn_setup_nmlrec_frazil_ice


   subroutine ocn_setup_nmlrec_land_ice_fluxes(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      character (len=StrKIND) :: config_land_ice_flux_mode = 'off'
      character (len=StrKIND) :: config_land_ice_flux_formulation = 'Jenkins'
      logical :: config_land_ice_flux_useHollandJenkinsAdvDiff = .false.
      real (kind=RKIND) :: config_land_ice_flux_attenuation_coefficient = 10.000000
      real (kind=RKIND) :: config_land_ice_flux_boundaryLayerThickness = 10.000000
      real (kind=RKIND) :: config_land_ice_flux_boundaryLayerNeighborWeight = 0.000000
      real (kind=RKIND) :: config_land_ice_flux_cp_ice = 2009.000000
      real (kind=RKIND) :: config_land_ice_flux_rho_ice = 918.000000
      real (kind=RKIND) :: config_land_ice_flux_explicit_topDragCoeff = 0.002500
      real (kind=RKIND) :: config_land_ice_flux_ISOMIP_gammaT = 0.000100
      real (kind=RKIND) :: config_land_ice_flux_rms_tidal_velocity = 0.050000
      real (kind=RKIND) :: config_land_ice_flux_jenkins_heat_transfer_coefficient = 0.011000
      real (kind=RKIND) :: config_land_ice_flux_jenkins_salt_transfer_coefficient = 0.000310

      namelist /land_ice_fluxes/ &
         config_land_ice_flux_mode, &
         config_land_ice_flux_formulation, &
         config_land_ice_flux_useHollandJenkinsAdvDiff, &
         config_land_ice_flux_attenuation_coefficient, &
         config_land_ice_flux_boundaryLayerThickness, &
         config_land_ice_flux_boundaryLayerNeighborWeight, &
         config_land_ice_flux_cp_ice, &
         config_land_ice_flux_rho_ice, &
         config_land_ice_flux_explicit_topDragCoeff, &
         config_land_ice_flux_ISOMIP_gammaT, &
         config_land_ice_flux_rms_tidal_velocity, &
         config_land_ice_flux_jenkins_heat_transfer_coefficient, &
         config_land_ice_flux_jenkins_salt_transfer_coefficient
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, land_ice_fluxes, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_char(dminfo, config_land_ice_flux_mode)
         call mpas_dmpar_bcast_char(dminfo, config_land_ice_flux_formulation)
         call mpas_dmpar_bcast_logical(dminfo, config_land_ice_flux_useHollandJenkinsAdvDiff)
         call mpas_dmpar_bcast_real(dminfo, config_land_ice_flux_attenuation_coefficient)
         call mpas_dmpar_bcast_real(dminfo, config_land_ice_flux_boundaryLayerThickness)
         call mpas_dmpar_bcast_real(dminfo, config_land_ice_flux_boundaryLayerNeighborWeight)
         call mpas_dmpar_bcast_real(dminfo, config_land_ice_flux_cp_ice)
         call mpas_dmpar_bcast_real(dminfo, config_land_ice_flux_rho_ice)
         call mpas_dmpar_bcast_real(dminfo, config_land_ice_flux_explicit_topDragCoeff)
         call mpas_dmpar_bcast_real(dminfo, config_land_ice_flux_ISOMIP_gammaT)
         call mpas_dmpar_bcast_real(dminfo, config_land_ice_flux_rms_tidal_velocity)
         call mpas_dmpar_bcast_real(dminfo, config_land_ice_flux_jenkins_heat_transfer_coefficient)
         call mpas_dmpar_bcast_real(dminfo, config_land_ice_flux_jenkins_salt_transfer_coefficient)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record land_ice_fluxes')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_land_ice_flux_mode = '//mpas_log_escape_dollars(config_land_ice_flux_mode))
            call mpas_log_write('        config_land_ice_flux_formulation = '&
//''//mpas_log_escape_dollars(config_land_ice_flux_formulation))
            call mpas_log_write('        config_land_ice_flux_useHollandJenkinsAdvDiff = $l', &
logicArgs=(/config_land_ice_flux_useHollandJenkinsAdvDiff/))
            call mpas_log_write('        config_land_ice_flux_attenuation_coefficient = $r', &
realArgs=(/config_land_ice_flux_attenuation_coefficient/))
            call mpas_log_write('        config_land_ice_flux_boundaryLayerThickness = $r', &
realArgs=(/config_land_ice_flux_boundaryLayerThickness/))
            call mpas_log_write('        config_land_ice_flux_boundaryLayerNeighborWeight = $r', &
realArgs=(/config_land_ice_flux_boundaryLayerNeighborWeight/))
            call mpas_log_write('        config_land_ice_flux_cp_ice = $r', realArgs=(/config_land_ice_flux_cp_ice/))
            call mpas_log_write('        config_land_ice_flux_rho_ice = $r', realArgs=(/config_land_ice_flux_rho_ice/))
            call mpas_log_write('        config_land_ice_flux_explicit_topDragCoeff = $r', &
realArgs=(/config_land_ice_flux_explicit_topDragCoeff/))
            call mpas_log_write('        config_land_ice_flux_ISOMIP_gammaT = $r', realArgs=(/config_land_ice_flux_ISOMIP_gammaT/))
            call mpas_log_write('        config_land_ice_flux_rms_tidal_velocity = $r', &
realArgs=(/config_land_ice_flux_rms_tidal_velocity/))
            call mpas_log_write('        config_land_ice_flux_jenkins_heat_transfer_coefficient = $r', &
realArgs=(/config_land_ice_flux_jenkins_heat_transfer_coefficient/))
            call mpas_log_write('        config_land_ice_flux_jenkins_salt_transfer_coefficient = $r', &
realArgs=(/config_land_ice_flux_jenkins_salt_transfer_coefficient/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record land_ice_fluxes.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_land_ice_flux_mode', config_land_ice_flux_mode)
      call mpas_pool_add_config(configPool, 'config_land_ice_flux_formulation', config_land_ice_flux_formulation)
      call mpas_pool_add_config(configPool, 'config_land_ice_flux_useHollandJenkinsAdvDiff', &
config_land_ice_flux_useHollandJenkinsAdvDiff)
      call mpas_pool_add_config(configPool, 'config_land_ice_flux_attenuation_coefficient', &
config_land_ice_flux_attenuation_coefficient)
      call mpas_pool_add_config(configPool, 'config_land_ice_flux_boundaryLayerThickness', &
config_land_ice_flux_boundaryLayerThickness)
      call mpas_pool_add_config(configPool, 'config_land_ice_flux_boundaryLayerNeighborWeight', &
config_land_ice_flux_boundaryLayerNeighborWeight)
      call mpas_pool_add_config(configPool, 'config_land_ice_flux_cp_ice', config_land_ice_flux_cp_ice)
      call mpas_pool_add_config(configPool, 'config_land_ice_flux_rho_ice', config_land_ice_flux_rho_ice)
      call mpas_pool_add_config(configPool, 'config_land_ice_flux_explicit_topDragCoeff', &
config_land_ice_flux_explicit_topDragCoeff)
      call mpas_pool_add_config(configPool, 'config_land_ice_flux_ISOMIP_gammaT', config_land_ice_flux_ISOMIP_gammaT)
      call mpas_pool_add_config(configPool, 'config_land_ice_flux_rms_tidal_velocity', config_land_ice_flux_rms_tidal_velocity)
      call mpas_pool_add_config(configPool, 'config_land_ice_flux_jenkins_heat_transfer_coefficient', &
config_land_ice_flux_jenkins_heat_transfer_coefficient)
      call mpas_pool_add_config(configPool, 'config_land_ice_flux_jenkins_salt_transfer_coefficient', &
config_land_ice_flux_jenkins_salt_transfer_coefficient)

   end subroutine ocn_setup_nmlrec_land_ice_fluxes


   subroutine ocn_setup_nmlrec_advection(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      character (len=StrKIND) :: config_vert_advection_method = 'flux-form'
      integer :: config_vert_remap_order = 3
      integer :: config_vert_remap_interval = 0
      integer :: config_vert_tracer_adv_flux_order = 3
      integer :: config_horiz_tracer_adv_order = 3
      real (kind=RKIND) :: config_coef_3rd_order = 0.250000
      character (len=StrKIND) :: config_flux_limiter = 'monotonic'
      character (len=StrKIND) :: config_remap_limiter = 'monotonic'
      character (len=StrKIND) :: config_thickness_flux_type = 'centered'

      namelist /advection/ &
         config_vert_advection_method, &
         config_vert_remap_order, &
         config_vert_remap_interval, &
         config_vert_tracer_adv_flux_order, &
         config_horiz_tracer_adv_order, &
         config_coef_3rd_order, &
         config_flux_limiter, &
         config_remap_limiter, &
         config_thickness_flux_type
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, advection, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_char(dminfo, config_vert_advection_method)
         call mpas_dmpar_bcast_int(dminfo, config_vert_remap_order)
         call mpas_dmpar_bcast_int(dminfo, config_vert_remap_interval)
         call mpas_dmpar_bcast_int(dminfo, config_vert_tracer_adv_flux_order)
         call mpas_dmpar_bcast_int(dminfo, config_horiz_tracer_adv_order)
         call mpas_dmpar_bcast_real(dminfo, config_coef_3rd_order)
         call mpas_dmpar_bcast_char(dminfo, config_flux_limiter)
         call mpas_dmpar_bcast_char(dminfo, config_remap_limiter)
         call mpas_dmpar_bcast_char(dminfo, config_thickness_flux_type)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record advection')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_vert_advection_method = '//mpas_log_escape_dollars(config_vert_advection_method))
            call mpas_log_write('        config_vert_remap_order = $i', intArgs=(/config_vert_remap_order/))
            call mpas_log_write('        config_vert_remap_interval = $i', intArgs=(/config_vert_remap_interval/))
            call mpas_log_write('        config_vert_tracer_adv_flux_order = $i', intArgs=(/config_vert_tracer_adv_flux_order/))
            call mpas_log_write('        config_horiz_tracer_adv_order = $i', intArgs=(/config_horiz_tracer_adv_order/))
            call mpas_log_write('        config_coef_3rd_order = $r', realArgs=(/config_coef_3rd_order/))
            call mpas_log_write('        config_flux_limiter = '//mpas_log_escape_dollars(config_flux_limiter))
            call mpas_log_write('        config_remap_limiter = '//mpas_log_escape_dollars(config_remap_limiter))
            call mpas_log_write('        config_thickness_flux_type = '//mpas_log_escape_dollars(config_thickness_flux_type))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record advection.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_vert_advection_method', config_vert_advection_method)
      call mpas_pool_add_config(configPool, 'config_vert_remap_order', config_vert_remap_order)
      call mpas_pool_add_config(configPool, 'config_vert_remap_interval', config_vert_remap_interval)
      call mpas_pool_add_config(configPool, 'config_vert_tracer_adv_flux_order', config_vert_tracer_adv_flux_order)
      call mpas_pool_add_config(configPool, 'config_horiz_tracer_adv_order', config_horiz_tracer_adv_order)
      call mpas_pool_add_config(configPool, 'config_coef_3rd_order', config_coef_3rd_order)
      call mpas_pool_add_config(configPool, 'config_flux_limiter', config_flux_limiter)
      call mpas_pool_add_config(configPool, 'config_remap_limiter', config_remap_limiter)
      call mpas_pool_add_config(configPool, 'config_thickness_flux_type', config_thickness_flux_type)

   end subroutine ocn_setup_nmlrec_advection


   subroutine ocn_setup_nmlrec_bottom_drag(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      character (len=StrKIND) :: config_bottom_drag_mode = 'implicit'
      character (len=StrKIND) :: config_implicit_bottom_drag_type = 'constant'
      real (kind=RKIND) :: config_implicit_constant_bottom_drag_coeff = 0.001000
      logical :: config_use_implicit_top_drag = .false.
      real (kind=RKIND) :: config_implicit_top_drag_coeff = 0.001000
      real (kind=RKIND) :: config_loglaw_bottom_roughness = 0.001000
      real (kind=RKIND) :: config_loglaw_layer_depth_max = 10.000000
      real (kind=RKIND) :: config_loglaw_bottom_drag_min = 0.002500
      real (kind=RKIND) :: config_loglaw_bottom_drag_max = 0.100000
      real (kind=RKIND) :: config_explicit_bottom_drag_coeff = 0.001000
      logical :: config_use_topographic_wave_drag = .false.
      real (kind=RKIND) :: config_topographic_wave_drag_coeff = 0.000500
      character (len=StrKIND) :: config_thickness_drag_type = 'centered'

      namelist /bottom_drag/ &
         config_bottom_drag_mode, &
         config_implicit_bottom_drag_type, &
         config_implicit_constant_bottom_drag_coeff, &
         config_use_implicit_top_drag, &
         config_implicit_top_drag_coeff, &
         config_loglaw_bottom_roughness, &
         config_loglaw_layer_depth_max, &
         config_loglaw_bottom_drag_min, &
         config_loglaw_bottom_drag_max, &
         config_explicit_bottom_drag_coeff, &
         config_use_topographic_wave_drag, &
         config_topographic_wave_drag_coeff, &
         config_thickness_drag_type
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, bottom_drag, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_char(dminfo, config_bottom_drag_mode)
         call mpas_dmpar_bcast_char(dminfo, config_implicit_bottom_drag_type)
         call mpas_dmpar_bcast_real(dminfo, config_implicit_constant_bottom_drag_coeff)
         call mpas_dmpar_bcast_logical(dminfo, config_use_implicit_top_drag)
         call mpas_dmpar_bcast_real(dminfo, config_implicit_top_drag_coeff)
         call mpas_dmpar_bcast_real(dminfo, config_loglaw_bottom_roughness)
         call mpas_dmpar_bcast_real(dminfo, config_loglaw_layer_depth_max)
         call mpas_dmpar_bcast_real(dminfo, config_loglaw_bottom_drag_min)
         call mpas_dmpar_bcast_real(dminfo, config_loglaw_bottom_drag_max)
         call mpas_dmpar_bcast_real(dminfo, config_explicit_bottom_drag_coeff)
         call mpas_dmpar_bcast_logical(dminfo, config_use_topographic_wave_drag)
         call mpas_dmpar_bcast_real(dminfo, config_topographic_wave_drag_coeff)
         call mpas_dmpar_bcast_char(dminfo, config_thickness_drag_type)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record bottom_drag')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_bottom_drag_mode = '//mpas_log_escape_dollars(config_bottom_drag_mode))
            call mpas_log_write('        config_implicit_bottom_drag_type = '&
//''//mpas_log_escape_dollars(config_implicit_bottom_drag_type))
            call mpas_log_write('        config_implicit_constant_bottom_drag_coeff = $r', &
realArgs=(/config_implicit_constant_bottom_drag_coeff/))
            call mpas_log_write('        config_use_implicit_top_drag = $l', logicArgs=(/config_use_implicit_top_drag/))
            call mpas_log_write('        config_implicit_top_drag_coeff = $r', realArgs=(/config_implicit_top_drag_coeff/))
            call mpas_log_write('        config_loglaw_bottom_roughness = $r', realArgs=(/config_loglaw_bottom_roughness/))
            call mpas_log_write('        config_loglaw_layer_depth_max = $r', realArgs=(/config_loglaw_layer_depth_max/))
            call mpas_log_write('        config_loglaw_bottom_drag_min = $r', realArgs=(/config_loglaw_bottom_drag_min/))
            call mpas_log_write('        config_loglaw_bottom_drag_max = $r', realArgs=(/config_loglaw_bottom_drag_max/))
            call mpas_log_write('        config_explicit_bottom_drag_coeff = $r', realArgs=(/config_explicit_bottom_drag_coeff/))
            call mpas_log_write('        config_use_topographic_wave_drag = $l', logicArgs=(/config_use_topographic_wave_drag/))
            call mpas_log_write('        config_topographic_wave_drag_coeff = $r', realArgs=(/config_topographic_wave_drag_coeff/))
            call mpas_log_write('        config_thickness_drag_type = '//mpas_log_escape_dollars(config_thickness_drag_type))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record bottom_drag.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_bottom_drag_mode', config_bottom_drag_mode)
      call mpas_pool_add_config(configPool, 'config_implicit_bottom_drag_type', config_implicit_bottom_drag_type)
      call mpas_pool_add_config(configPool, 'config_implicit_constant_bottom_drag_coeff', &
config_implicit_constant_bottom_drag_coeff)
      call mpas_pool_add_config(configPool, 'config_use_implicit_top_drag', config_use_implicit_top_drag)
      call mpas_pool_add_config(configPool, 'config_implicit_top_drag_coeff', config_implicit_top_drag_coeff)
      call mpas_pool_add_config(configPool, 'config_loglaw_bottom_roughness', config_loglaw_bottom_roughness)
      call mpas_pool_add_config(configPool, 'config_loglaw_layer_depth_max', config_loglaw_layer_depth_max)
      call mpas_pool_add_config(configPool, 'config_loglaw_bottom_drag_min', config_loglaw_bottom_drag_min)
      call mpas_pool_add_config(configPool, 'config_loglaw_bottom_drag_max', config_loglaw_bottom_drag_max)
      call mpas_pool_add_config(configPool, 'config_explicit_bottom_drag_coeff', config_explicit_bottom_drag_coeff)
      call mpas_pool_add_config(configPool, 'config_use_topographic_wave_drag', config_use_topographic_wave_drag)
      call mpas_pool_add_config(configPool, 'config_topographic_wave_drag_coeff', config_topographic_wave_drag_coeff)
      call mpas_pool_add_config(configPool, 'config_thickness_drag_type', config_thickness_drag_type)

   end subroutine ocn_setup_nmlrec_bottom_drag


   subroutine ocn_setup_nmlrec_Rayleigh_damping(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      real (kind=RKIND) :: config_Rayleigh_damping_coeff = 0.000100
      logical :: config_Rayleigh_damping_depth_variable = .false.
      logical :: config_Rayleigh_bottom_friction = .false.
      real (kind=RKIND) :: config_Rayleigh_bottom_damping_coeff = 0.000100

      namelist /Rayleigh_damping/ &
         config_Rayleigh_damping_coeff, &
         config_Rayleigh_damping_depth_variable, &
         config_Rayleigh_bottom_friction, &
         config_Rayleigh_bottom_damping_coeff
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, Rayleigh_damping, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_real(dminfo, config_Rayleigh_damping_coeff)
         call mpas_dmpar_bcast_logical(dminfo, config_Rayleigh_damping_depth_variable)
         call mpas_dmpar_bcast_logical(dminfo, config_Rayleigh_bottom_friction)
         call mpas_dmpar_bcast_real(dminfo, config_Rayleigh_bottom_damping_coeff)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record Rayleigh_damping')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_Rayleigh_damping_coeff = $r', realArgs=(/config_Rayleigh_damping_coeff/))
            call mpas_log_write('        config_Rayleigh_damping_depth_variable = $l', &
logicArgs=(/config_Rayleigh_damping_depth_variable/))
            call mpas_log_write('        config_Rayleigh_bottom_friction = $l', logicArgs=(/config_Rayleigh_bottom_friction/))
            call mpas_log_write('        config_Rayleigh_bottom_damping_coeff = $r', &
realArgs=(/config_Rayleigh_bottom_damping_coeff/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record Rayleigh_damping.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_Rayleigh_damping_coeff', config_Rayleigh_damping_coeff)
      call mpas_pool_add_config(configPool, 'config_Rayleigh_damping_depth_variable', config_Rayleigh_damping_depth_variable)
      call mpas_pool_add_config(configPool, 'config_Rayleigh_bottom_friction', config_Rayleigh_bottom_friction)
      call mpas_pool_add_config(configPool, 'config_Rayleigh_bottom_damping_coeff', config_Rayleigh_bottom_damping_coeff)

   end subroutine ocn_setup_nmlrec_Rayleigh_damping


   subroutine ocn_setup_nmlrec_vegetation_drag(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_use_vegetation_drag = .false.
      real (kind=RKIND) :: config_vegetation_drag_coefficient = 1.090000

      namelist /vegetation_drag/ &
         config_use_vegetation_drag, &
         config_vegetation_drag_coefficient
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, vegetation_drag, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_use_vegetation_drag)
         call mpas_dmpar_bcast_real(dminfo, config_vegetation_drag_coefficient)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record vegetation_drag')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_use_vegetation_drag = $l', logicArgs=(/config_use_vegetation_drag/))
            call mpas_log_write('        config_vegetation_drag_coefficient = $r', realArgs=(/config_vegetation_drag_coefficient/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record vegetation_drag.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_use_vegetation_drag', config_use_vegetation_drag)
      call mpas_pool_add_config(configPool, 'config_vegetation_drag_coefficient', config_vegetation_drag_coefficient)

   end subroutine ocn_setup_nmlrec_vegetation_drag


   subroutine ocn_setup_nmlrec_wetting_drying(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_use_wetting_drying = .false.
      logical :: config_prevent_drying = .false.
      real (kind=RKIND) :: config_drying_min_cell_height = 0.001000
      logical :: config_zero_drying_velocity = .false.
      logical :: config_zero_drying_velocity_ramp = .false.
      real (kind=RKIND) :: config_zero_drying_velocity_ramp_hmin = 0.001000
      real (kind=RKIND) :: config_zero_drying_velocity_ramp_hmax = 0.002000
      logical :: config_verify_not_dry = .false.
      real (kind=RKIND) :: config_drying_safety_height = 0.000000

      namelist /wetting_drying/ &
         config_use_wetting_drying, &
         config_prevent_drying, &
         config_drying_min_cell_height, &
         config_zero_drying_velocity, &
         config_zero_drying_velocity_ramp, &
         config_zero_drying_velocity_ramp_hmin, &
         config_zero_drying_velocity_ramp_hmax, &
         config_verify_not_dry, &
         config_drying_safety_height
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, wetting_drying, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_use_wetting_drying)
         call mpas_dmpar_bcast_logical(dminfo, config_prevent_drying)
         call mpas_dmpar_bcast_real(dminfo, config_drying_min_cell_height)
         call mpas_dmpar_bcast_logical(dminfo, config_zero_drying_velocity)
         call mpas_dmpar_bcast_logical(dminfo, config_zero_drying_velocity_ramp)
         call mpas_dmpar_bcast_real(dminfo, config_zero_drying_velocity_ramp_hmin)
         call mpas_dmpar_bcast_real(dminfo, config_zero_drying_velocity_ramp_hmax)
         call mpas_dmpar_bcast_logical(dminfo, config_verify_not_dry)
         call mpas_dmpar_bcast_real(dminfo, config_drying_safety_height)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record wetting_drying')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_use_wetting_drying = $l', logicArgs=(/config_use_wetting_drying/))
            call mpas_log_write('        config_prevent_drying = $l', logicArgs=(/config_prevent_drying/))
            call mpas_log_write('        config_drying_min_cell_height = $r', realArgs=(/config_drying_min_cell_height/))
            call mpas_log_write('        config_zero_drying_velocity = $l', logicArgs=(/config_zero_drying_velocity/))
            call mpas_log_write('        config_zero_drying_velocity_ramp = $l', logicArgs=(/config_zero_drying_velocity_ramp/))
            call mpas_log_write('        config_zero_drying_velocity_ramp_hmin = $r', &
realArgs=(/config_zero_drying_velocity_ramp_hmin/))
            call mpas_log_write('        config_zero_drying_velocity_ramp_hmax = $r', &
realArgs=(/config_zero_drying_velocity_ramp_hmax/))
            call mpas_log_write('        config_verify_not_dry = $l', logicArgs=(/config_verify_not_dry/))
            call mpas_log_write('        config_drying_safety_height = $r', realArgs=(/config_drying_safety_height/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record wetting_drying.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_use_wetting_drying', config_use_wetting_drying)
      call mpas_pool_add_config(configPool, 'config_prevent_drying', config_prevent_drying)
      call mpas_pool_add_config(configPool, 'config_drying_min_cell_height', config_drying_min_cell_height)
      call mpas_pool_add_config(configPool, 'config_zero_drying_velocity', config_zero_drying_velocity)
      call mpas_pool_add_config(configPool, 'config_zero_drying_velocity_ramp', config_zero_drying_velocity_ramp)
      call mpas_pool_add_config(configPool, 'config_zero_drying_velocity_ramp_hmin', config_zero_drying_velocity_ramp_hmin)
      call mpas_pool_add_config(configPool, 'config_zero_drying_velocity_ramp_hmax', config_zero_drying_velocity_ramp_hmax)
      call mpas_pool_add_config(configPool, 'config_verify_not_dry', config_verify_not_dry)
      call mpas_pool_add_config(configPool, 'config_drying_safety_height', config_drying_safety_height)

   end subroutine ocn_setup_nmlrec_wetting_drying


   subroutine ocn_setup_nmlrec_ocean_constants(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      real (kind=RKIND) :: config_density0 = 1026.000000

      namelist /ocean_constants/ &
         config_density0
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, ocean_constants, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_real(dminfo, config_density0)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record ocean_constants')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_density0 = $r', realArgs=(/config_density0/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record ocean_constants.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_density0', config_density0)

   end subroutine ocn_setup_nmlrec_ocean_constants


   subroutine ocn_setup_nmlrec_lts(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      integer :: config_dt_scaling_LTS = 1

      namelist /lts/ &
         config_dt_scaling_LTS
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, lts, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_int(dminfo, config_dt_scaling_LTS)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record lts')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_dt_scaling_LTS = $i', intArgs=(/config_dt_scaling_LTS/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record lts.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_dt_scaling_LTS', config_dt_scaling_LTS)

   end subroutine ocn_setup_nmlrec_lts


   subroutine ocn_setup_nmlrec_pressure_gradient(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      character (len=StrKIND) :: config_pressure_gradient_type = 'pressure_and_zmid'
      real (kind=RKIND) :: config_common_level_weight = 0.500000
      real (kind=RKIND) :: config_zonal_ssh_grad = 0.000000
      real (kind=RKIND) :: config_meridional_ssh_grad = 0.000000

      namelist /pressure_gradient/ &
         config_pressure_gradient_type, &
         config_common_level_weight, &
         config_zonal_ssh_grad, &
         config_meridional_ssh_grad
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, pressure_gradient, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_char(dminfo, config_pressure_gradient_type)
         call mpas_dmpar_bcast_real(dminfo, config_common_level_weight)
         call mpas_dmpar_bcast_real(dminfo, config_zonal_ssh_grad)
         call mpas_dmpar_bcast_real(dminfo, config_meridional_ssh_grad)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record pressure_gradient')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_pressure_gradient_type = '//mpas_log_escape_dollars(config_pressure_gradient_type))
            call mpas_log_write('        config_common_level_weight = $r', realArgs=(/config_common_level_weight/))
            call mpas_log_write('        config_zonal_ssh_grad = $r', realArgs=(/config_zonal_ssh_grad/))
            call mpas_log_write('        config_meridional_ssh_grad = $r', realArgs=(/config_meridional_ssh_grad/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record pressure_gradient.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_pressure_gradient_type', config_pressure_gradient_type)
      call mpas_pool_add_config(configPool, 'config_common_level_weight', config_common_level_weight)
      call mpas_pool_add_config(configPool, 'config_zonal_ssh_grad', config_zonal_ssh_grad)
      call mpas_pool_add_config(configPool, 'config_meridional_ssh_grad', config_meridional_ssh_grad)

   end subroutine ocn_setup_nmlrec_pressure_gradient


   subroutine ocn_setup_nmlrec_eos(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      character (len=StrKIND) :: config_eos_type = 'linear'
      real (kind=RKIND) :: config_open_ocean_freezing_temperature_coeff_0 = 0.000000
      real (kind=RKIND) :: config_open_ocean_freezing_temperature_coeff_S = 0.000000
      real (kind=RKIND) :: config_open_ocean_freezing_temperature_coeff_p = 0.000000
      real (kind=RKIND) :: config_open_ocean_freezing_temperature_coeff_pS = 0.000000
      real (kind=RKIND) :: config_open_ocean_freezing_temperature_coeff_mushy_az1_liq = -18.480000
      real (kind=RKIND) :: config_land_ice_cavity_freezing_temperature_coeff_0 = 0.062200
      real (kind=RKIND) :: config_land_ice_cavity_freezing_temperature_coeff_S = -0.056300
      real (kind=RKIND) :: config_land_ice_cavity_freezing_temperature_coeff_p = -0.000000
      real (kind=RKIND) :: config_land_ice_cavity_freezing_temperature_coeff_pS = -0.000000

      namelist /eos/ &
         config_eos_type, &
         config_open_ocean_freezing_temperature_coeff_0, &
         config_open_ocean_freezing_temperature_coeff_S, &
         config_open_ocean_freezing_temperature_coeff_p, &
         config_open_ocean_freezing_temperature_coeff_pS, &
         config_open_ocean_freezing_temperature_coeff_mushy_az1_liq, &
         config_land_ice_cavity_freezing_temperature_coeff_0, &
         config_land_ice_cavity_freezing_temperature_coeff_S, &
         config_land_ice_cavity_freezing_temperature_coeff_p, &
         config_land_ice_cavity_freezing_temperature_coeff_pS
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, eos, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_char(dminfo, config_eos_type)
         call mpas_dmpar_bcast_real(dminfo, config_open_ocean_freezing_temperature_coeff_0)
         call mpas_dmpar_bcast_real(dminfo, config_open_ocean_freezing_temperature_coeff_S)
         call mpas_dmpar_bcast_real(dminfo, config_open_ocean_freezing_temperature_coeff_p)
         call mpas_dmpar_bcast_real(dminfo, config_open_ocean_freezing_temperature_coeff_pS)
         call mpas_dmpar_bcast_real(dminfo, config_open_ocean_freezing_temperature_coeff_mushy_az1_liq)
         call mpas_dmpar_bcast_real(dminfo, config_land_ice_cavity_freezing_temperature_coeff_0)
         call mpas_dmpar_bcast_real(dminfo, config_land_ice_cavity_freezing_temperature_coeff_S)
         call mpas_dmpar_bcast_real(dminfo, config_land_ice_cavity_freezing_temperature_coeff_p)
         call mpas_dmpar_bcast_real(dminfo, config_land_ice_cavity_freezing_temperature_coeff_pS)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record eos')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_eos_type = '//mpas_log_escape_dollars(config_eos_type))
            call mpas_log_write('        config_open_ocean_freezing_temperature_coeff_0 = $r', &
realArgs=(/config_open_ocean_freezing_temperature_coeff_0/))
            call mpas_log_write('        config_open_ocean_freezing_temperature_coeff_S = $r', &
realArgs=(/config_open_ocean_freezing_temperature_coeff_S/))
            call mpas_log_write('        config_open_ocean_freezing_temperature_coeff_p = $r', &
realArgs=(/config_open_ocean_freezing_temperature_coeff_p/))
            call mpas_log_write('        config_open_ocean_freezing_temperature_coeff_pS = $r', &
realArgs=(/config_open_ocean_freezing_temperature_coeff_pS/))
            call mpas_log_write('        config_open_ocean_freezing_temperature_coeff_mushy_az1_liq = $r', &
realArgs=(/config_open_ocean_freezing_temperature_coeff_mushy_az1_liq/))
            call mpas_log_write('        config_land_ice_cavity_freezing_temperature_coeff_0 = $r', &
realArgs=(/config_land_ice_cavity_freezing_temperature_coeff_0/))
            call mpas_log_write('        config_land_ice_cavity_freezing_temperature_coeff_S = $r', &
realArgs=(/config_land_ice_cavity_freezing_temperature_coeff_S/))
            call mpas_log_write('        config_land_ice_cavity_freezing_temperature_coeff_p = $r', &
realArgs=(/config_land_ice_cavity_freezing_temperature_coeff_p/))
            call mpas_log_write('        config_land_ice_cavity_freezing_temperature_coeff_pS = $r', &
realArgs=(/config_land_ice_cavity_freezing_temperature_coeff_pS/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record eos.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_eos_type', config_eos_type)
      call mpas_pool_add_config(configPool, 'config_open_ocean_freezing_temperature_coeff_0', &
config_open_ocean_freezing_temperature_coeff_0)
      call mpas_pool_add_config(configPool, 'config_open_ocean_freezing_temperature_coeff_S', &
config_open_ocean_freezing_temperature_coeff_S)
      call mpas_pool_add_config(configPool, 'config_open_ocean_freezing_temperature_coeff_p', &
config_open_ocean_freezing_temperature_coeff_p)
      call mpas_pool_add_config(configPool, 'config_open_ocean_freezing_temperature_coeff_pS', &
config_open_ocean_freezing_temperature_coeff_pS)
      call mpas_pool_add_config(configPool, 'config_open_ocean_freezing_temperature_coeff_mushy_az1_liq', &
config_open_ocean_freezing_temperature_coeff_mushy_az1_liq)
      call mpas_pool_add_config(configPool, 'config_land_ice_cavity_freezing_temperature_coeff_0', &
config_land_ice_cavity_freezing_temperature_coeff_0)
      call mpas_pool_add_config(configPool, 'config_land_ice_cavity_freezing_temperature_coeff_S', &
config_land_ice_cavity_freezing_temperature_coeff_S)
      call mpas_pool_add_config(configPool, 'config_land_ice_cavity_freezing_temperature_coeff_p', &
config_land_ice_cavity_freezing_temperature_coeff_p)
      call mpas_pool_add_config(configPool, 'config_land_ice_cavity_freezing_temperature_coeff_pS', &
config_land_ice_cavity_freezing_temperature_coeff_pS)

   end subroutine ocn_setup_nmlrec_eos


   subroutine ocn_setup_nmlrec_eos_linear(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      real (kind=RKIND) :: config_eos_linear_alpha = 0.200000
      real (kind=RKIND) :: config_eos_linear_beta = 0.800000
      real (kind=RKIND) :: config_eos_linear_Tref = 5.000000
      real (kind=RKIND) :: config_eos_linear_Sref = 35.000000
      real (kind=RKIND) :: config_eos_linear_densityref = 1000.000000

      namelist /eos_linear/ &
         config_eos_linear_alpha, &
         config_eos_linear_beta, &
         config_eos_linear_Tref, &
         config_eos_linear_Sref, &
         config_eos_linear_densityref
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, eos_linear, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_real(dminfo, config_eos_linear_alpha)
         call mpas_dmpar_bcast_real(dminfo, config_eos_linear_beta)
         call mpas_dmpar_bcast_real(dminfo, config_eos_linear_Tref)
         call mpas_dmpar_bcast_real(dminfo, config_eos_linear_Sref)
         call mpas_dmpar_bcast_real(dminfo, config_eos_linear_densityref)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record eos_linear')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_eos_linear_alpha = $r', realArgs=(/config_eos_linear_alpha/))
            call mpas_log_write('        config_eos_linear_beta = $r', realArgs=(/config_eos_linear_beta/))
            call mpas_log_write('        config_eos_linear_Tref = $r', realArgs=(/config_eos_linear_Tref/))
            call mpas_log_write('        config_eos_linear_Sref = $r', realArgs=(/config_eos_linear_Sref/))
            call mpas_log_write('        config_eos_linear_densityref = $r', realArgs=(/config_eos_linear_densityref/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record eos_linear.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_eos_linear_alpha', config_eos_linear_alpha)
      call mpas_pool_add_config(configPool, 'config_eos_linear_beta', config_eos_linear_beta)
      call mpas_pool_add_config(configPool, 'config_eos_linear_Tref', config_eos_linear_Tref)
      call mpas_pool_add_config(configPool, 'config_eos_linear_Sref', config_eos_linear_Sref)
      call mpas_pool_add_config(configPool, 'config_eos_linear_densityref', config_eos_linear_densityref)

   end subroutine ocn_setup_nmlrec_eos_linear


   subroutine ocn_setup_nmlrec_eos_wright(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      real (kind=RKIND) :: config_eos_wright_ref_pressure = 0.000000

      namelist /eos_wright/ &
         config_eos_wright_ref_pressure
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, eos_wright, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_real(dminfo, config_eos_wright_ref_pressure)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record eos_wright')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_eos_wright_ref_pressure = $r', realArgs=(/config_eos_wright_ref_pressure/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record eos_wright.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_eos_wright_ref_pressure', config_eos_wright_ref_pressure)

   end subroutine ocn_setup_nmlrec_eos_wright


   subroutine ocn_setup_nmlrec_split_timestep_share(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      integer :: config_n_ts_iter = 2
      integer :: config_n_bcl_iter_beg = 1
      integer :: config_n_bcl_iter_mid = 2
      integer :: config_n_bcl_iter_end = 2

      namelist /split_timestep_share/ &
         config_n_ts_iter, &
         config_n_bcl_iter_beg, &
         config_n_bcl_iter_mid, &
         config_n_bcl_iter_end
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, split_timestep_share, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_int(dminfo, config_n_ts_iter)
         call mpas_dmpar_bcast_int(dminfo, config_n_bcl_iter_beg)
         call mpas_dmpar_bcast_int(dminfo, config_n_bcl_iter_mid)
         call mpas_dmpar_bcast_int(dminfo, config_n_bcl_iter_end)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record split_timestep_share')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_n_ts_iter = $i', intArgs=(/config_n_ts_iter/))
            call mpas_log_write('        config_n_bcl_iter_beg = $i', intArgs=(/config_n_bcl_iter_beg/))
            call mpas_log_write('        config_n_bcl_iter_mid = $i', intArgs=(/config_n_bcl_iter_mid/))
            call mpas_log_write('        config_n_bcl_iter_end = $i', intArgs=(/config_n_bcl_iter_end/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record split_timestep_share.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_n_ts_iter', config_n_ts_iter)
      call mpas_pool_add_config(configPool, 'config_n_bcl_iter_beg', config_n_bcl_iter_beg)
      call mpas_pool_add_config(configPool, 'config_n_bcl_iter_mid', config_n_bcl_iter_mid)
      call mpas_pool_add_config(configPool, 'config_n_bcl_iter_end', config_n_bcl_iter_end)

   end subroutine ocn_setup_nmlrec_split_timestep_share


   subroutine ocn_setup_nmlrec_split_explicit_ts(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      character (len=StrKIND) :: config_btr_dt = '0000_00:00:15'
      integer :: config_n_btr_cor_iter = 2
      logical :: config_vel_correction = .true.
      integer :: config_btr_subcycle_loop_factor = 2
      real (kind=RKIND) :: config_btr_gam1_velWt1 = 0.533300
      real (kind=RKIND) :: config_btr_gam2_SSHWt1 = 0.533300
      real (kind=RKIND) :: config_btr_gam3_velWt2 = 1.000000
      logical :: config_btr_solve_SSH2 = .false.

      namelist /split_explicit_ts/ &
         config_btr_dt, &
         config_n_btr_cor_iter, &
         config_vel_correction, &
         config_btr_subcycle_loop_factor, &
         config_btr_gam1_velWt1, &
         config_btr_gam2_SSHWt1, &
         config_btr_gam3_velWt2, &
         config_btr_solve_SSH2
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, split_explicit_ts, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_char(dminfo, config_btr_dt)
         call mpas_dmpar_bcast_int(dminfo, config_n_btr_cor_iter)
         call mpas_dmpar_bcast_logical(dminfo, config_vel_correction)
         call mpas_dmpar_bcast_int(dminfo, config_btr_subcycle_loop_factor)
         call mpas_dmpar_bcast_real(dminfo, config_btr_gam1_velWt1)
         call mpas_dmpar_bcast_real(dminfo, config_btr_gam2_SSHWt1)
         call mpas_dmpar_bcast_real(dminfo, config_btr_gam3_velWt2)
         call mpas_dmpar_bcast_logical(dminfo, config_btr_solve_SSH2)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record split_explicit_ts')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_btr_dt = '//mpas_log_escape_dollars(config_btr_dt))
            call mpas_log_write('        config_n_btr_cor_iter = $i', intArgs=(/config_n_btr_cor_iter/))
            call mpas_log_write('        config_vel_correction = $l', logicArgs=(/config_vel_correction/))
            call mpas_log_write('        config_btr_subcycle_loop_factor = $i', intArgs=(/config_btr_subcycle_loop_factor/))
            call mpas_log_write('        config_btr_gam1_velWt1 = $r', realArgs=(/config_btr_gam1_velWt1/))
            call mpas_log_write('        config_btr_gam2_SSHWt1 = $r', realArgs=(/config_btr_gam2_SSHWt1/))
            call mpas_log_write('        config_btr_gam3_velWt2 = $r', realArgs=(/config_btr_gam3_velWt2/))
            call mpas_log_write('        config_btr_solve_SSH2 = $l', logicArgs=(/config_btr_solve_SSH2/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record split_explicit_ts.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_btr_dt', config_btr_dt)
      call mpas_pool_add_config(configPool, 'config_n_btr_cor_iter', config_n_btr_cor_iter)
      call mpas_pool_add_config(configPool, 'config_vel_correction', config_vel_correction)
      call mpas_pool_add_config(configPool, 'config_btr_subcycle_loop_factor', config_btr_subcycle_loop_factor)
      call mpas_pool_add_config(configPool, 'config_btr_gam1_velWt1', config_btr_gam1_velWt1)
      call mpas_pool_add_config(configPool, 'config_btr_gam2_SSHWt1', config_btr_gam2_SSHWt1)
      call mpas_pool_add_config(configPool, 'config_btr_gam3_velWt2', config_btr_gam3_velWt2)
      call mpas_pool_add_config(configPool, 'config_btr_solve_SSH2', config_btr_solve_SSH2)

   end subroutine ocn_setup_nmlrec_split_explicit_ts


   subroutine ocn_setup_nmlrec_split_implicit_ts(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      character (len=StrKIND) :: config_btr_si_preconditioner = 'ras'
      real (kind=RKIND) :: config_btr_si_tolerance = 0.000000
      integer :: config_n_btr_si_large_iter = 1
      logical :: config_btr_si_partition_match_mode = .false.

      namelist /split_implicit_ts/ &
         config_btr_si_preconditioner, &
         config_btr_si_tolerance, &
         config_n_btr_si_large_iter, &
         config_btr_si_partition_match_mode
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, split_implicit_ts, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_char(dminfo, config_btr_si_preconditioner)
         call mpas_dmpar_bcast_real(dminfo, config_btr_si_tolerance)
         call mpas_dmpar_bcast_int(dminfo, config_n_btr_si_large_iter)
         call mpas_dmpar_bcast_logical(dminfo, config_btr_si_partition_match_mode)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record split_implicit_ts')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_btr_si_preconditioner = '//mpas_log_escape_dollars(config_btr_si_preconditioner))
            call mpas_log_write('        config_btr_si_tolerance = $r', realArgs=(/config_btr_si_tolerance/))
            call mpas_log_write('        config_n_btr_si_large_iter = $i', intArgs=(/config_n_btr_si_large_iter/))
            call mpas_log_write('        config_btr_si_partition_match_mode = $l', &
logicArgs=(/config_btr_si_partition_match_mode/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record split_implicit_ts.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_btr_si_preconditioner', config_btr_si_preconditioner)
      call mpas_pool_add_config(configPool, 'config_btr_si_tolerance', config_btr_si_tolerance)
      call mpas_pool_add_config(configPool, 'config_n_btr_si_large_iter', config_n_btr_si_large_iter)
      call mpas_pool_add_config(configPool, 'config_btr_si_partition_match_mode', config_btr_si_partition_match_mode)

   end subroutine ocn_setup_nmlrec_split_implicit_ts


   subroutine ocn_setup_nmlrec_ALE_vertical_grid(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      character (len=StrKIND) :: config_vert_coord_movement = 'uniform_stretching'
      character (len=StrKIND) :: config_ALE_thickness_proportionality = 'restingThickness_times_weights'
      real (kind=RKIND) :: config_vert_taper_weight_depth_1 = 250.000000
      real (kind=RKIND) :: config_vert_taper_weight_depth_2 = 500.000000
      logical :: config_use_min_max_thickness = .false.
      real (kind=RKIND) :: config_min_thickness = 1.000000
      real (kind=RKIND) :: config_max_thickness_factor = 6.000000
      logical :: config_dzdk_positive = .false.

      namelist /ALE_vertical_grid/ &
         config_vert_coord_movement, &
         config_ALE_thickness_proportionality, &
         config_vert_taper_weight_depth_1, &
         config_vert_taper_weight_depth_2, &
         config_use_min_max_thickness, &
         config_min_thickness, &
         config_max_thickness_factor, &
         config_dzdk_positive
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, ALE_vertical_grid, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_char(dminfo, config_vert_coord_movement)
         call mpas_dmpar_bcast_char(dminfo, config_ALE_thickness_proportionality)
         call mpas_dmpar_bcast_real(dminfo, config_vert_taper_weight_depth_1)
         call mpas_dmpar_bcast_real(dminfo, config_vert_taper_weight_depth_2)
         call mpas_dmpar_bcast_logical(dminfo, config_use_min_max_thickness)
         call mpas_dmpar_bcast_real(dminfo, config_min_thickness)
         call mpas_dmpar_bcast_real(dminfo, config_max_thickness_factor)
         call mpas_dmpar_bcast_logical(dminfo, config_dzdk_positive)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record ALE_vertical_grid')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_vert_coord_movement = '//mpas_log_escape_dollars(config_vert_coord_movement))
            call mpas_log_write('        config_ALE_thickness_proportionality = '&
//''//mpas_log_escape_dollars(config_ALE_thickness_proportionality))
            call mpas_log_write('        config_vert_taper_weight_depth_1 = $r', realArgs=(/config_vert_taper_weight_depth_1/))
            call mpas_log_write('        config_vert_taper_weight_depth_2 = $r', realArgs=(/config_vert_taper_weight_depth_2/))
            call mpas_log_write('        config_use_min_max_thickness = $l', logicArgs=(/config_use_min_max_thickness/))
            call mpas_log_write('        config_min_thickness = $r', realArgs=(/config_min_thickness/))
            call mpas_log_write('        config_max_thickness_factor = $r', realArgs=(/config_max_thickness_factor/))
            call mpas_log_write('        config_dzdk_positive = $l', logicArgs=(/config_dzdk_positive/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record ALE_vertical_grid.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_vert_coord_movement', config_vert_coord_movement)
      call mpas_pool_add_config(configPool, 'config_ALE_thickness_proportionality', config_ALE_thickness_proportionality)
      call mpas_pool_add_config(configPool, 'config_vert_taper_weight_depth_1', config_vert_taper_weight_depth_1)
      call mpas_pool_add_config(configPool, 'config_vert_taper_weight_depth_2', config_vert_taper_weight_depth_2)
      call mpas_pool_add_config(configPool, 'config_use_min_max_thickness', config_use_min_max_thickness)
      call mpas_pool_add_config(configPool, 'config_min_thickness', config_min_thickness)
      call mpas_pool_add_config(configPool, 'config_max_thickness_factor', config_max_thickness_factor)
      call mpas_pool_add_config(configPool, 'config_dzdk_positive', config_dzdk_positive)

   end subroutine ocn_setup_nmlrec_ALE_vertical_grid


   subroutine ocn_setup_nmlrec_ALE_frequency_filtered_thickness(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_use_freq_filtered_thickness = .false.
      real (kind=RKIND) :: config_thickness_filter_timescale = 5.000000
      logical :: config_use_highFreqThick_restore = .false.
      real (kind=RKIND) :: config_highFreqThick_restore_time = 30.000000
      logical :: config_use_highFreqThick_del2 = .false.
      real (kind=RKIND) :: config_highFreqThick_del2 = 100.000000

      namelist /ALE_frequency_filtered_thickness/ &
         config_use_freq_filtered_thickness, &
         config_thickness_filter_timescale, &
         config_use_highFreqThick_restore, &
         config_highFreqThick_restore_time, &
         config_use_highFreqThick_del2, &
         config_highFreqThick_del2
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, ALE_frequency_filtered_thickness, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_use_freq_filtered_thickness)
         call mpas_dmpar_bcast_real(dminfo, config_thickness_filter_timescale)
         call mpas_dmpar_bcast_logical(dminfo, config_use_highFreqThick_restore)
         call mpas_dmpar_bcast_real(dminfo, config_highFreqThick_restore_time)
         call mpas_dmpar_bcast_logical(dminfo, config_use_highFreqThick_del2)
         call mpas_dmpar_bcast_real(dminfo, config_highFreqThick_del2)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record '&
//'ALE_frequency_filtered_thickness')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_use_freq_filtered_thickness = $l', &
logicArgs=(/config_use_freq_filtered_thickness/))
            call mpas_log_write('        config_thickness_filter_timescale = $r', realArgs=(/config_thickness_filter_timescale/))
            call mpas_log_write('        config_use_highFreqThick_restore = $l', logicArgs=(/config_use_highFreqThick_restore/))
            call mpas_log_write('        config_highFreqThick_restore_time = $r', realArgs=(/config_highFreqThick_restore_time/))
            call mpas_log_write('        config_use_highFreqThick_del2 = $l', logicArgs=(/config_use_highFreqThick_del2/))
            call mpas_log_write('        config_highFreqThick_del2 = $r', realArgs=(/config_highFreqThick_del2/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record ALE_frequency_filtered_thickness.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_use_freq_filtered_thickness', config_use_freq_filtered_thickness)
      call mpas_pool_add_config(configPool, 'config_thickness_filter_timescale', config_thickness_filter_timescale)
      call mpas_pool_add_config(configPool, 'config_use_highFreqThick_restore', config_use_highFreqThick_restore)
      call mpas_pool_add_config(configPool, 'config_highFreqThick_restore_time', config_highFreqThick_restore_time)
      call mpas_pool_add_config(configPool, 'config_use_highFreqThick_del2', config_use_highFreqThick_del2)
      call mpas_pool_add_config(configPool, 'config_highFreqThick_del2', config_highFreqThick_del2)

   end subroutine ocn_setup_nmlrec_ALE_frequency_filtered_thickness


   subroutine ocn_setup_nmlrec_debug(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_check_zlevel_consistency = .false.
      logical :: config_check_ssh_consistency = .true.
      logical :: config_filter_btr_mode = .false.
      logical :: config_prescribe_velocity = .false.
      logical :: config_prescribe_thickness = .false.
      logical :: config_include_KE_vertex = .false.
      logical :: config_check_tracer_monotonicity = .false.
      logical :: config_compute_active_tracer_budgets = .true.
      logical :: config_disable_thick_all_tend = .false.
      logical :: config_disable_thick_hadv = .false.
      logical :: config_disable_thick_vadv = .false.
      logical :: config_disable_thick_sflux = .false.
      logical :: config_disable_vel_all_tend = .false.
      logical :: config_disable_vel_hadv = .false.
      logical :: config_disable_vel_coriolis = .false.
      logical :: config_disable_vel_pgrad = .false.
      logical :: config_disable_vel_hmix = .false.
      logical :: config_disable_vel_surface_stress = .false.
      logical :: config_disable_vel_topographic_wave_drag = .false.
      logical :: config_disable_vel_explicit_bottom_drag = .false.
      logical :: config_disable_vel_vmix = .false.
      logical :: config_disable_vel_vadv = .false.
      logical :: config_disable_tr_all_tend = .false.
      logical :: config_disable_tr_adv = .false.
      logical :: config_disable_tr_hmix = .false.
      logical :: config_disable_tr_vmix = .false.
      logical :: config_disable_tr_sflux = .false.
      logical :: config_disable_tr_nonlocalflux = .false.
      logical :: config_disable_redi_k33 = .false.
      logical :: config_read_nearest_restart = .false.

      namelist /debug/ &
         config_check_zlevel_consistency, &
         config_check_ssh_consistency, &
         config_filter_btr_mode, &
         config_prescribe_velocity, &
         config_prescribe_thickness, &
         config_include_KE_vertex, &
         config_check_tracer_monotonicity, &
         config_compute_active_tracer_budgets, &
         config_disable_thick_all_tend, &
         config_disable_thick_hadv, &
         config_disable_thick_vadv, &
         config_disable_thick_sflux, &
         config_disable_vel_all_tend, &
         config_disable_vel_hadv, &
         config_disable_vel_coriolis, &
         config_disable_vel_pgrad, &
         config_disable_vel_hmix, &
         config_disable_vel_surface_stress, &
         config_disable_vel_topographic_wave_drag, &
         config_disable_vel_explicit_bottom_drag, &
         config_disable_vel_vmix, &
         config_disable_vel_vadv, &
         config_disable_tr_all_tend, &
         config_disable_tr_adv, &
         config_disable_tr_hmix, &
         config_disable_tr_vmix, &
         config_disable_tr_sflux, &
         config_disable_tr_nonlocalflux, &
         config_disable_redi_k33, &
         config_read_nearest_restart
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, debug, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_check_zlevel_consistency)
         call mpas_dmpar_bcast_logical(dminfo, config_check_ssh_consistency)
         call mpas_dmpar_bcast_logical(dminfo, config_filter_btr_mode)
         call mpas_dmpar_bcast_logical(dminfo, config_prescribe_velocity)
         call mpas_dmpar_bcast_logical(dminfo, config_prescribe_thickness)
         call mpas_dmpar_bcast_logical(dminfo, config_include_KE_vertex)
         call mpas_dmpar_bcast_logical(dminfo, config_check_tracer_monotonicity)
         call mpas_dmpar_bcast_logical(dminfo, config_compute_active_tracer_budgets)
         call mpas_dmpar_bcast_logical(dminfo, config_disable_thick_all_tend)
         call mpas_dmpar_bcast_logical(dminfo, config_disable_thick_hadv)
         call mpas_dmpar_bcast_logical(dminfo, config_disable_thick_vadv)
         call mpas_dmpar_bcast_logical(dminfo, config_disable_thick_sflux)
         call mpas_dmpar_bcast_logical(dminfo, config_disable_vel_all_tend)
         call mpas_dmpar_bcast_logical(dminfo, config_disable_vel_hadv)
         call mpas_dmpar_bcast_logical(dminfo, config_disable_vel_coriolis)
         call mpas_dmpar_bcast_logical(dminfo, config_disable_vel_pgrad)
         call mpas_dmpar_bcast_logical(dminfo, config_disable_vel_hmix)
         call mpas_dmpar_bcast_logical(dminfo, config_disable_vel_surface_stress)
         call mpas_dmpar_bcast_logical(dminfo, config_disable_vel_topographic_wave_drag)
         call mpas_dmpar_bcast_logical(dminfo, config_disable_vel_explicit_bottom_drag)
         call mpas_dmpar_bcast_logical(dminfo, config_disable_vel_vmix)
         call mpas_dmpar_bcast_logical(dminfo, config_disable_vel_vadv)
         call mpas_dmpar_bcast_logical(dminfo, config_disable_tr_all_tend)
         call mpas_dmpar_bcast_logical(dminfo, config_disable_tr_adv)
         call mpas_dmpar_bcast_logical(dminfo, config_disable_tr_hmix)
         call mpas_dmpar_bcast_logical(dminfo, config_disable_tr_vmix)
         call mpas_dmpar_bcast_logical(dminfo, config_disable_tr_sflux)
         call mpas_dmpar_bcast_logical(dminfo, config_disable_tr_nonlocalflux)
         call mpas_dmpar_bcast_logical(dminfo, config_disable_redi_k33)
         call mpas_dmpar_bcast_logical(dminfo, config_read_nearest_restart)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record debug')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_check_zlevel_consistency = $l', logicArgs=(/config_check_zlevel_consistency/))
            call mpas_log_write('        config_check_ssh_consistency = $l', logicArgs=(/config_check_ssh_consistency/))
            call mpas_log_write('        config_filter_btr_mode = $l', logicArgs=(/config_filter_btr_mode/))
            call mpas_log_write('        config_prescribe_velocity = $l', logicArgs=(/config_prescribe_velocity/))
            call mpas_log_write('        config_prescribe_thickness = $l', logicArgs=(/config_prescribe_thickness/))
            call mpas_log_write('        config_include_KE_vertex = $l', logicArgs=(/config_include_KE_vertex/))
            call mpas_log_write('        config_check_tracer_monotonicity = $l', logicArgs=(/config_check_tracer_monotonicity/))
            call mpas_log_write('        config_compute_active_tracer_budgets = $l', &
logicArgs=(/config_compute_active_tracer_budgets/))
            call mpas_log_write('        config_disable_thick_all_tend = $l', logicArgs=(/config_disable_thick_all_tend/))
            call mpas_log_write('        config_disable_thick_hadv = $l', logicArgs=(/config_disable_thick_hadv/))
            call mpas_log_write('        config_disable_thick_vadv = $l', logicArgs=(/config_disable_thick_vadv/))
            call mpas_log_write('        config_disable_thick_sflux = $l', logicArgs=(/config_disable_thick_sflux/))
            call mpas_log_write('        config_disable_vel_all_tend = $l', logicArgs=(/config_disable_vel_all_tend/))
            call mpas_log_write('        config_disable_vel_hadv = $l', logicArgs=(/config_disable_vel_hadv/))
            call mpas_log_write('        config_disable_vel_coriolis = $l', logicArgs=(/config_disable_vel_coriolis/))
            call mpas_log_write('        config_disable_vel_pgrad = $l', logicArgs=(/config_disable_vel_pgrad/))
            call mpas_log_write('        config_disable_vel_hmix = $l', logicArgs=(/config_disable_vel_hmix/))
            call mpas_log_write('        config_disable_vel_surface_stress = $l', logicArgs=(/config_disable_vel_surface_stress/))
            call mpas_log_write('        config_disable_vel_topographic_wave_drag = $l', &
logicArgs=(/config_disable_vel_topographic_wave_drag/))
            call mpas_log_write('        config_disable_vel_explicit_bottom_drag = $l', &
logicArgs=(/config_disable_vel_explicit_bottom_drag/))
            call mpas_log_write('        config_disable_vel_vmix = $l', logicArgs=(/config_disable_vel_vmix/))
            call mpas_log_write('        config_disable_vel_vadv = $l', logicArgs=(/config_disable_vel_vadv/))
            call mpas_log_write('        config_disable_tr_all_tend = $l', logicArgs=(/config_disable_tr_all_tend/))
            call mpas_log_write('        config_disable_tr_adv = $l', logicArgs=(/config_disable_tr_adv/))
            call mpas_log_write('        config_disable_tr_hmix = $l', logicArgs=(/config_disable_tr_hmix/))
            call mpas_log_write('        config_disable_tr_vmix = $l', logicArgs=(/config_disable_tr_vmix/))
            call mpas_log_write('        config_disable_tr_sflux = $l', logicArgs=(/config_disable_tr_sflux/))
            call mpas_log_write('        config_disable_tr_nonlocalflux = $l', logicArgs=(/config_disable_tr_nonlocalflux/))
            call mpas_log_write('        config_disable_redi_k33 = $l', logicArgs=(/config_disable_redi_k33/))
            call mpas_log_write('        config_read_nearest_restart = $l', logicArgs=(/config_read_nearest_restart/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record debug.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_check_zlevel_consistency', config_check_zlevel_consistency)
      call mpas_pool_add_config(configPool, 'config_check_ssh_consistency', config_check_ssh_consistency)
      call mpas_pool_add_config(configPool, 'config_filter_btr_mode', config_filter_btr_mode)
      call mpas_pool_add_config(configPool, 'config_prescribe_velocity', config_prescribe_velocity)
      call mpas_pool_add_config(configPool, 'config_prescribe_thickness', config_prescribe_thickness)
      call mpas_pool_add_config(configPool, 'config_include_KE_vertex', config_include_KE_vertex)
      call mpas_pool_add_config(configPool, 'config_check_tracer_monotonicity', config_check_tracer_monotonicity)
      call mpas_pool_add_config(configPool, 'config_compute_active_tracer_budgets', config_compute_active_tracer_budgets)
      call mpas_pool_add_config(configPool, 'config_disable_thick_all_tend', config_disable_thick_all_tend)
      call mpas_pool_add_config(configPool, 'config_disable_thick_hadv', config_disable_thick_hadv)
      call mpas_pool_add_config(configPool, 'config_disable_thick_vadv', config_disable_thick_vadv)
      call mpas_pool_add_config(configPool, 'config_disable_thick_sflux', config_disable_thick_sflux)
      call mpas_pool_add_config(configPool, 'config_disable_vel_all_tend', config_disable_vel_all_tend)
      call mpas_pool_add_config(configPool, 'config_disable_vel_hadv', config_disable_vel_hadv)
      call mpas_pool_add_config(configPool, 'config_disable_vel_coriolis', config_disable_vel_coriolis)
      call mpas_pool_add_config(configPool, 'config_disable_vel_pgrad', config_disable_vel_pgrad)
      call mpas_pool_add_config(configPool, 'config_disable_vel_hmix', config_disable_vel_hmix)
      call mpas_pool_add_config(configPool, 'config_disable_vel_surface_stress', config_disable_vel_surface_stress)
      call mpas_pool_add_config(configPool, 'config_disable_vel_topographic_wave_drag', config_disable_vel_topographic_wave_drag)
      call mpas_pool_add_config(configPool, 'config_disable_vel_explicit_bottom_drag', config_disable_vel_explicit_bottom_drag)
      call mpas_pool_add_config(configPool, 'config_disable_vel_vmix', config_disable_vel_vmix)
      call mpas_pool_add_config(configPool, 'config_disable_vel_vadv', config_disable_vel_vadv)
      call mpas_pool_add_config(configPool, 'config_disable_tr_all_tend', config_disable_tr_all_tend)
      call mpas_pool_add_config(configPool, 'config_disable_tr_adv', config_disable_tr_adv)
      call mpas_pool_add_config(configPool, 'config_disable_tr_hmix', config_disable_tr_hmix)
      call mpas_pool_add_config(configPool, 'config_disable_tr_vmix', config_disable_tr_vmix)
      call mpas_pool_add_config(configPool, 'config_disable_tr_sflux', config_disable_tr_sflux)
      call mpas_pool_add_config(configPool, 'config_disable_tr_nonlocalflux', config_disable_tr_nonlocalflux)
      call mpas_pool_add_config(configPool, 'config_disable_redi_k33', config_disable_redi_k33)
      call mpas_pool_add_config(configPool, 'config_read_nearest_restart', config_read_nearest_restart)

   end subroutine ocn_setup_nmlrec_debug


   subroutine ocn_setup_nmlrec_testing(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_conduct_tests = .false.
      logical :: config_test_tensors = .false.
      character (len=StrKIND) :: config_tensor_test_function = 'sph_uCosCos'

      namelist /testing/ &
         config_conduct_tests, &
         config_test_tensors, &
         config_tensor_test_function
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, testing, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_conduct_tests)
         call mpas_dmpar_bcast_logical(dminfo, config_test_tensors)
         call mpas_dmpar_bcast_char(dminfo, config_tensor_test_function)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record testing')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_conduct_tests = $l', logicArgs=(/config_conduct_tests/))
            call mpas_log_write('        config_test_tensors = $l', logicArgs=(/config_test_tensors/))
            call mpas_log_write('        config_tensor_test_function = '//mpas_log_escape_dollars(config_tensor_test_function))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record testing.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_conduct_tests', config_conduct_tests)
      call mpas_pool_add_config(configPool, 'config_test_tensors', config_test_tensors)
      call mpas_pool_add_config(configPool, 'config_tensor_test_function', config_tensor_test_function)

   end subroutine ocn_setup_nmlrec_testing


   subroutine ocn_setup_nmlrec_transport_tests(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      integer :: config_transport_tests_vert_levels = 3
      real (kind=RKIND) :: config_transport_tests_temperature = 15.000000
      real (kind=RKIND) :: config_transport_tests_salinity = 35.000000
      integer :: config_transport_tests_flow_id = 0

      namelist /transport_tests/ &
         config_transport_tests_vert_levels, &
         config_transport_tests_temperature, &
         config_transport_tests_salinity, &
         config_transport_tests_flow_id
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, transport_tests, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_int(dminfo, config_transport_tests_vert_levels)
         call mpas_dmpar_bcast_real(dminfo, config_transport_tests_temperature)
         call mpas_dmpar_bcast_real(dminfo, config_transport_tests_salinity)
         call mpas_dmpar_bcast_int(dminfo, config_transport_tests_flow_id)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record transport_tests')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_transport_tests_vert_levels = $i', intArgs=(/config_transport_tests_vert_levels/))
            call mpas_log_write('        config_transport_tests_temperature = $r', realArgs=(/config_transport_tests_temperature/))
            call mpas_log_write('        config_transport_tests_salinity = $r', realArgs=(/config_transport_tests_salinity/))
            call mpas_log_write('        config_transport_tests_flow_id = $i', intArgs=(/config_transport_tests_flow_id/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record transport_tests.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_transport_tests_vert_levels', config_transport_tests_vert_levels)
      call mpas_pool_add_config(configPool, 'config_transport_tests_temperature', config_transport_tests_temperature)
      call mpas_pool_add_config(configPool, 'config_transport_tests_salinity', config_transport_tests_salinity)
      call mpas_pool_add_config(configPool, 'config_transport_tests_flow_id', config_transport_tests_flow_id)

   end subroutine ocn_setup_nmlrec_transport_tests


   subroutine ocn_setup_nmlrec_init_mode_vert_levels(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      integer :: config_vert_levels = -1

      namelist /init_mode_vert_levels/ &
         config_vert_levels
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, init_mode_vert_levels, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_int(dminfo, config_vert_levels)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record init_mode_vert_levels')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_vert_levels = $i', intArgs=(/config_vert_levels/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record init_mode_vert_levels.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_vert_levels', config_vert_levels)

   end subroutine ocn_setup_nmlrec_init_mode_vert_levels


   subroutine ocn_setup_nmlrec_manufactured_solution(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_use_manufactured_solution = .false.
      real (kind=RKIND) :: config_manufactured_solution_wavelength_x = 2000000.000000
      real (kind=RKIND) :: config_manufactured_solution_wavelength_y = 2000000.000000
      real (kind=RKIND) :: config_manufactured_solution_amplitude = 1.000000

      namelist /manufactured_solution/ &
         config_use_manufactured_solution, &
         config_manufactured_solution_wavelength_x, &
         config_manufactured_solution_wavelength_y, &
         config_manufactured_solution_amplitude
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, manufactured_solution, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_use_manufactured_solution)
         call mpas_dmpar_bcast_real(dminfo, config_manufactured_solution_wavelength_x)
         call mpas_dmpar_bcast_real(dminfo, config_manufactured_solution_wavelength_y)
         call mpas_dmpar_bcast_real(dminfo, config_manufactured_solution_amplitude)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record manufactured_solution')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_use_manufactured_solution = $l', logicArgs=(/config_use_manufactured_solution/))
            call mpas_log_write('        config_manufactured_solution_wavelength_x = $r', &
realArgs=(/config_manufactured_solution_wavelength_x/))
            call mpas_log_write('        config_manufactured_solution_wavelength_y = $r', &
realArgs=(/config_manufactured_solution_wavelength_y/))
            call mpas_log_write('        config_manufactured_solution_amplitude = $r', &
realArgs=(/config_manufactured_solution_amplitude/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record manufactured_solution.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_use_manufactured_solution', config_use_manufactured_solution)
      call mpas_pool_add_config(configPool, 'config_manufactured_solution_wavelength_x', config_manufactured_solution_wavelength_x)
      call mpas_pool_add_config(configPool, 'config_manufactured_solution_wavelength_y', config_manufactured_solution_wavelength_y)
      call mpas_pool_add_config(configPool, 'config_manufactured_solution_amplitude', config_manufactured_solution_amplitude)

   end subroutine ocn_setup_nmlrec_manufactured_solution


   subroutine ocn_setup_nmlrec_tracer_forcing_activeTracers(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_use_activeTracers = .true.
      logical :: config_use_activeTracers_surface_bulk_forcing = .false.
      logical :: config_use_activeTracers_surface_restoring = .false.
      logical :: config_use_activeTracers_interior_restoring = .false.
      logical :: config_use_activeTracers_exponential_decay = .false.
      logical :: config_use_activeTracers_idealAge_forcing = .false.
      logical :: config_use_activeTracers_ttd_forcing = .false.
      logical :: config_use_surface_salinity_monthly_restoring = .false.
      character (len=StrKIND) :: config_surface_salinity_monthly_restoring_compute_interval = '0000-00-01_00:00:00'
      real (kind=RKIND) :: config_salinity_restoring_constant_piston_velocity = 0.000016
      real (kind=RKIND) :: config_salinity_restoring_max_difference = 100.000000
      logical :: config_salinity_restoring_under_sea_ice = .false.

      namelist /tracer_forcing_activeTracers/ &
         config_use_activeTracers, &
         config_use_activeTracers_surface_bulk_forcing, &
         config_use_activeTracers_surface_restoring, &
         config_use_activeTracers_interior_restoring, &
         config_use_activeTracers_exponential_decay, &
         config_use_activeTracers_idealAge_forcing, &
         config_use_activeTracers_ttd_forcing, &
         config_use_surface_salinity_monthly_restoring, &
         config_surface_salinity_monthly_restoring_compute_interval, &
         config_salinity_restoring_constant_piston_velocity, &
         config_salinity_restoring_max_difference, &
         config_salinity_restoring_under_sea_ice
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, tracer_forcing_activeTracers, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_use_activeTracers)
         call mpas_dmpar_bcast_logical(dminfo, config_use_activeTracers_surface_bulk_forcing)
         call mpas_dmpar_bcast_logical(dminfo, config_use_activeTracers_surface_restoring)
         call mpas_dmpar_bcast_logical(dminfo, config_use_activeTracers_interior_restoring)
         call mpas_dmpar_bcast_logical(dminfo, config_use_activeTracers_exponential_decay)
         call mpas_dmpar_bcast_logical(dminfo, config_use_activeTracers_idealAge_forcing)
         call mpas_dmpar_bcast_logical(dminfo, config_use_activeTracers_ttd_forcing)
         call mpas_dmpar_bcast_logical(dminfo, config_use_surface_salinity_monthly_restoring)
         call mpas_dmpar_bcast_char(dminfo, config_surface_salinity_monthly_restoring_compute_interval)
         call mpas_dmpar_bcast_real(dminfo, config_salinity_restoring_constant_piston_velocity)
         call mpas_dmpar_bcast_real(dminfo, config_salinity_restoring_max_difference)
         call mpas_dmpar_bcast_logical(dminfo, config_salinity_restoring_under_sea_ice)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record tracer_forcing_activeTracers')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_use_activeTracers = $l', logicArgs=(/config_use_activeTracers/))
            call mpas_log_write('        config_use_activeTracers_surface_bulk_forcing = $l', &
logicArgs=(/config_use_activeTracers_surface_bulk_forcing/))
            call mpas_log_write('        config_use_activeTracers_surface_restoring = $l', &
logicArgs=(/config_use_activeTracers_surface_restoring/))
            call mpas_log_write('        config_use_activeTracers_interior_restoring = $l', &
logicArgs=(/config_use_activeTracers_interior_restoring/))
            call mpas_log_write('        config_use_activeTracers_exponential_decay = $l', &
logicArgs=(/config_use_activeTracers_exponential_decay/))
            call mpas_log_write('        config_use_activeTracers_idealAge_forcing = $l', &
logicArgs=(/config_use_activeTracers_idealAge_forcing/))
            call mpas_log_write('        config_use_activeTracers_ttd_forcing = $l', &
logicArgs=(/config_use_activeTracers_ttd_forcing/))
            call mpas_log_write('        config_use_surface_salinity_monthly_restoring = $l', &
logicArgs=(/config_use_surface_salinity_monthly_restoring/))
            call mpas_log_write('        config_surface_salinity_monthly_restoring_compute_interval = '&
//''//mpas_log_escape_dollars(config_surface_salinity_monthly_restoring_compute_interval))
            call mpas_log_write('        config_salinity_restoring_constant_piston_velocity = $r', &
realArgs=(/config_salinity_restoring_constant_piston_velocity/))
            call mpas_log_write('        config_salinity_restoring_max_difference = $r', &
realArgs=(/config_salinity_restoring_max_difference/))
            call mpas_log_write('        config_salinity_restoring_under_sea_ice = $l', &
logicArgs=(/config_salinity_restoring_under_sea_ice/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record tracer_forcing_activeTracers.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_use_activeTracers', config_use_activeTracers)
      call mpas_pool_add_config(configPool, 'config_use_activeTracers_surface_bulk_forcing', &
config_use_activeTracers_surface_bulk_forcing)
      call mpas_pool_add_config(configPool, 'config_use_activeTracers_surface_restoring', &
config_use_activeTracers_surface_restoring)
      call mpas_pool_add_config(configPool, 'config_use_activeTracers_interior_restoring', &
config_use_activeTracers_interior_restoring)
      call mpas_pool_add_config(configPool, 'config_use_activeTracers_exponential_decay', &
config_use_activeTracers_exponential_decay)
      call mpas_pool_add_config(configPool, 'config_use_activeTracers_idealAge_forcing', config_use_activeTracers_idealAge_forcing)
      call mpas_pool_add_config(configPool, 'config_use_activeTracers_ttd_forcing', config_use_activeTracers_ttd_forcing)
      call mpas_pool_add_config(configPool, 'config_use_surface_salinity_monthly_restoring', &
config_use_surface_salinity_monthly_restoring)
      call mpas_pool_add_config(configPool, 'config_surface_salinity_monthly_restoring_compute_interval', &
config_surface_salinity_monthly_restoring_compute_interval)
      call mpas_pool_add_config(configPool, 'config_salinity_restoring_constant_piston_velocity', &
config_salinity_restoring_constant_piston_velocity)
      call mpas_pool_add_config(configPool, 'config_salinity_restoring_max_difference', config_salinity_restoring_max_difference)
      call mpas_pool_add_config(configPool, 'config_salinity_restoring_under_sea_ice', config_salinity_restoring_under_sea_ice)

   end subroutine ocn_setup_nmlrec_tracer_forcing_activeTracers


   subroutine ocn_setup_nmlrec_tracer_forcing_debugTracers(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_use_debugTracers = .false.
      logical :: config_reset_debugTracers_near_surface = .false.
      integer :: config_reset_debugTracers_top_nLayers = 20
      logical :: config_use_debugTracers_surface_bulk_forcing = .false.
      logical :: config_use_debugTracers_surface_restoring = .false.
      logical :: config_use_debugTracers_interior_restoring = .false.
      logical :: config_use_debugTracers_exponential_decay = .false.
      logical :: config_use_debugTracers_idealAge_forcing = .false.
      logical :: config_use_debugTracers_ttd_forcing = .false.

      namelist /tracer_forcing_debugTracers/ &
         config_use_debugTracers, &
         config_reset_debugTracers_near_surface, &
         config_reset_debugTracers_top_nLayers, &
         config_use_debugTracers_surface_bulk_forcing, &
         config_use_debugTracers_surface_restoring, &
         config_use_debugTracers_interior_restoring, &
         config_use_debugTracers_exponential_decay, &
         config_use_debugTracers_idealAge_forcing, &
         config_use_debugTracers_ttd_forcing
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, tracer_forcing_debugTracers, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_use_debugTracers)
         call mpas_dmpar_bcast_logical(dminfo, config_reset_debugTracers_near_surface)
         call mpas_dmpar_bcast_int(dminfo, config_reset_debugTracers_top_nLayers)
         call mpas_dmpar_bcast_logical(dminfo, config_use_debugTracers_surface_bulk_forcing)
         call mpas_dmpar_bcast_logical(dminfo, config_use_debugTracers_surface_restoring)
         call mpas_dmpar_bcast_logical(dminfo, config_use_debugTracers_interior_restoring)
         call mpas_dmpar_bcast_logical(dminfo, config_use_debugTracers_exponential_decay)
         call mpas_dmpar_bcast_logical(dminfo, config_use_debugTracers_idealAge_forcing)
         call mpas_dmpar_bcast_logical(dminfo, config_use_debugTracers_ttd_forcing)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record tracer_forcing_debugTracers')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_use_debugTracers = $l', logicArgs=(/config_use_debugTracers/))
            call mpas_log_write('        config_reset_debugTracers_near_surface = $l', &
logicArgs=(/config_reset_debugTracers_near_surface/))
            call mpas_log_write('        config_reset_debugTracers_top_nLayers = $i', &
intArgs=(/config_reset_debugTracers_top_nLayers/))
            call mpas_log_write('        config_use_debugTracers_surface_bulk_forcing = $l', &
logicArgs=(/config_use_debugTracers_surface_bulk_forcing/))
            call mpas_log_write('        config_use_debugTracers_surface_restoring = $l', &
logicArgs=(/config_use_debugTracers_surface_restoring/))
            call mpas_log_write('        config_use_debugTracers_interior_restoring = $l', &
logicArgs=(/config_use_debugTracers_interior_restoring/))
            call mpas_log_write('        config_use_debugTracers_exponential_decay = $l', &
logicArgs=(/config_use_debugTracers_exponential_decay/))
            call mpas_log_write('        config_use_debugTracers_idealAge_forcing = $l', &
logicArgs=(/config_use_debugTracers_idealAge_forcing/))
            call mpas_log_write('        config_use_debugTracers_ttd_forcing = $l', &
logicArgs=(/config_use_debugTracers_ttd_forcing/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record tracer_forcing_debugTracers.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_use_debugTracers', config_use_debugTracers)
      call mpas_pool_add_config(configPool, 'config_reset_debugTracers_near_surface', config_reset_debugTracers_near_surface)
      call mpas_pool_add_config(configPool, 'config_reset_debugTracers_top_nLayers', config_reset_debugTracers_top_nLayers)
      call mpas_pool_add_config(configPool, 'config_use_debugTracers_surface_bulk_forcing', &
config_use_debugTracers_surface_bulk_forcing)
      call mpas_pool_add_config(configPool, 'config_use_debugTracers_surface_restoring', config_use_debugTracers_surface_restoring)
      call mpas_pool_add_config(configPool, 'config_use_debugTracers_interior_restoring', &
config_use_debugTracers_interior_restoring)
      call mpas_pool_add_config(configPool, 'config_use_debugTracers_exponential_decay', config_use_debugTracers_exponential_decay)
      call mpas_pool_add_config(configPool, 'config_use_debugTracers_idealAge_forcing', config_use_debugTracers_idealAge_forcing)
      call mpas_pool_add_config(configPool, 'config_use_debugTracers_ttd_forcing', config_use_debugTracers_ttd_forcing)

   end subroutine ocn_setup_nmlrec_tracer_forcing_debugTracers


   subroutine ocn_setup_nmlrec_tracer_forcing_ecosysTracers(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_use_ecosysTracers = .false.
      character (len=StrKIND) :: config_ecosys_atm_co2_option = 'none'
      character (len=StrKIND) :: config_ecosys_atm_alt_co2_option = 'none'
      logical :: config_ecosys_atm_alt_co2_use_eco = .true.
      real (kind=RKIND) :: config_ecosys_atm_co2_constant_value = 379.000000
      logical :: config_use_ecosysTracers_surface_bulk_forcing = .false.
      logical :: config_use_ecosysTracers_surface_restoring = .false.
      logical :: config_use_ecosysTracers_interior_restoring = .false.
      logical :: config_use_ecosysTracers_exponential_decay = .false.
      logical :: config_use_ecosysTracers_idealAge_forcing = .false.
      logical :: config_use_ecosysTracers_ttd_forcing = .false.
      logical :: config_use_ecosysTracers_surface_value = .false.
      logical :: config_use_ecosysTracers_river_inputs_from_coupler = .false.
      logical :: config_use_ecosysTracers_sea_ice_coupling = .false.
      logical :: config_ecosysTracers_diagnostic_fields_level1 = .false.
      logical :: config_ecosysTracers_diagnostic_fields_level2 = .false.
      logical :: config_ecosysTracers_diagnostic_fields_level3 = .false.
      logical :: config_ecosysTracers_diagnostic_fields_level4 = .false.
      logical :: config_ecosysTracers_diagnostic_fields_level5 = .false.

      namelist /tracer_forcing_ecosysTracers/ &
         config_use_ecosysTracers, &
         config_ecosys_atm_co2_option, &
         config_ecosys_atm_alt_co2_option, &
         config_ecosys_atm_alt_co2_use_eco, &
         config_ecosys_atm_co2_constant_value, &
         config_use_ecosysTracers_surface_bulk_forcing, &
         config_use_ecosysTracers_surface_restoring, &
         config_use_ecosysTracers_interior_restoring, &
         config_use_ecosysTracers_exponential_decay, &
         config_use_ecosysTracers_idealAge_forcing, &
         config_use_ecosysTracers_ttd_forcing, &
         config_use_ecosysTracers_surface_value, &
         config_use_ecosysTracers_river_inputs_from_coupler, &
         config_use_ecosysTracers_sea_ice_coupling, &
         config_ecosysTracers_diagnostic_fields_level1, &
         config_ecosysTracers_diagnostic_fields_level2, &
         config_ecosysTracers_diagnostic_fields_level3, &
         config_ecosysTracers_diagnostic_fields_level4, &
         config_ecosysTracers_diagnostic_fields_level5
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, tracer_forcing_ecosysTracers, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_use_ecosysTracers)
         call mpas_dmpar_bcast_char(dminfo, config_ecosys_atm_co2_option)
         call mpas_dmpar_bcast_char(dminfo, config_ecosys_atm_alt_co2_option)
         call mpas_dmpar_bcast_logical(dminfo, config_ecosys_atm_alt_co2_use_eco)
         call mpas_dmpar_bcast_real(dminfo, config_ecosys_atm_co2_constant_value)
         call mpas_dmpar_bcast_logical(dminfo, config_use_ecosysTracers_surface_bulk_forcing)
         call mpas_dmpar_bcast_logical(dminfo, config_use_ecosysTracers_surface_restoring)
         call mpas_dmpar_bcast_logical(dminfo, config_use_ecosysTracers_interior_restoring)
         call mpas_dmpar_bcast_logical(dminfo, config_use_ecosysTracers_exponential_decay)
         call mpas_dmpar_bcast_logical(dminfo, config_use_ecosysTracers_idealAge_forcing)
         call mpas_dmpar_bcast_logical(dminfo, config_use_ecosysTracers_ttd_forcing)
         call mpas_dmpar_bcast_logical(dminfo, config_use_ecosysTracers_surface_value)
         call mpas_dmpar_bcast_logical(dminfo, config_use_ecosysTracers_river_inputs_from_coupler)
         call mpas_dmpar_bcast_logical(dminfo, config_use_ecosysTracers_sea_ice_coupling)
         call mpas_dmpar_bcast_logical(dminfo, config_ecosysTracers_diagnostic_fields_level1)
         call mpas_dmpar_bcast_logical(dminfo, config_ecosysTracers_diagnostic_fields_level2)
         call mpas_dmpar_bcast_logical(dminfo, config_ecosysTracers_diagnostic_fields_level3)
         call mpas_dmpar_bcast_logical(dminfo, config_ecosysTracers_diagnostic_fields_level4)
         call mpas_dmpar_bcast_logical(dminfo, config_ecosysTracers_diagnostic_fields_level5)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record tracer_forcing_ecosysTracers')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_use_ecosysTracers = $l', logicArgs=(/config_use_ecosysTracers/))
            call mpas_log_write('        config_ecosys_atm_co2_option = '//mpas_log_escape_dollars(config_ecosys_atm_co2_option))
            call mpas_log_write('        config_ecosys_atm_alt_co2_option = '&
//''//mpas_log_escape_dollars(config_ecosys_atm_alt_co2_option))
            call mpas_log_write('        config_ecosys_atm_alt_co2_use_eco = $l', logicArgs=(/config_ecosys_atm_alt_co2_use_eco/))
            call mpas_log_write('        config_ecosys_atm_co2_constant_value = $r', &
realArgs=(/config_ecosys_atm_co2_constant_value/))
            call mpas_log_write('        config_use_ecosysTracers_surface_bulk_forcing = $l', &
logicArgs=(/config_use_ecosysTracers_surface_bulk_forcing/))
            call mpas_log_write('        config_use_ecosysTracers_surface_restoring = $l', &
logicArgs=(/config_use_ecosysTracers_surface_restoring/))
            call mpas_log_write('        config_use_ecosysTracers_interior_restoring = $l', &
logicArgs=(/config_use_ecosysTracers_interior_restoring/))
            call mpas_log_write('        config_use_ecosysTracers_exponential_decay = $l', &
logicArgs=(/config_use_ecosysTracers_exponential_decay/))
            call mpas_log_write('        config_use_ecosysTracers_idealAge_forcing = $l', &
logicArgs=(/config_use_ecosysTracers_idealAge_forcing/))
            call mpas_log_write('        config_use_ecosysTracers_ttd_forcing = $l', &
logicArgs=(/config_use_ecosysTracers_ttd_forcing/))
            call mpas_log_write('        config_use_ecosysTracers_surface_value = $l', &
logicArgs=(/config_use_ecosysTracers_surface_value/))
            call mpas_log_write('        config_use_ecosysTracers_river_inputs_from_coupler = $l', &
logicArgs=(/config_use_ecosysTracers_river_inputs_from_coupler/))
            call mpas_log_write('        config_use_ecosysTracers_sea_ice_coupling = $l', &
logicArgs=(/config_use_ecosysTracers_sea_ice_coupling/))
            call mpas_log_write('        config_ecosysTracers_diagnostic_fields_level1 = $l', &
logicArgs=(/config_ecosysTracers_diagnostic_fields_level1/))
            call mpas_log_write('        config_ecosysTracers_diagnostic_fields_level2 = $l', &
logicArgs=(/config_ecosysTracers_diagnostic_fields_level2/))
            call mpas_log_write('        config_ecosysTracers_diagnostic_fields_level3 = $l', &
logicArgs=(/config_ecosysTracers_diagnostic_fields_level3/))
            call mpas_log_write('        config_ecosysTracers_diagnostic_fields_level4 = $l', &
logicArgs=(/config_ecosysTracers_diagnostic_fields_level4/))
            call mpas_log_write('        config_ecosysTracers_diagnostic_fields_level5 = $l', &
logicArgs=(/config_ecosysTracers_diagnostic_fields_level5/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record tracer_forcing_ecosysTracers.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_use_ecosysTracers', config_use_ecosysTracers)
      call mpas_pool_add_config(configPool, 'config_ecosys_atm_co2_option', config_ecosys_atm_co2_option)
      call mpas_pool_add_config(configPool, 'config_ecosys_atm_alt_co2_option', config_ecosys_atm_alt_co2_option)
      call mpas_pool_add_config(configPool, 'config_ecosys_atm_alt_co2_use_eco', config_ecosys_atm_alt_co2_use_eco)
      call mpas_pool_add_config(configPool, 'config_ecosys_atm_co2_constant_value', config_ecosys_atm_co2_constant_value)
      call mpas_pool_add_config(configPool, 'config_use_ecosysTracers_surface_bulk_forcing', &
config_use_ecosysTracers_surface_bulk_forcing)
      call mpas_pool_add_config(configPool, 'config_use_ecosysTracers_surface_restoring', &
config_use_ecosysTracers_surface_restoring)
      call mpas_pool_add_config(configPool, 'config_use_ecosysTracers_interior_restoring', &
config_use_ecosysTracers_interior_restoring)
      call mpas_pool_add_config(configPool, 'config_use_ecosysTracers_exponential_decay', &
config_use_ecosysTracers_exponential_decay)
      call mpas_pool_add_config(configPool, 'config_use_ecosysTracers_idealAge_forcing', config_use_ecosysTracers_idealAge_forcing)
      call mpas_pool_add_config(configPool, 'config_use_ecosysTracers_ttd_forcing', config_use_ecosysTracers_ttd_forcing)
      call mpas_pool_add_config(configPool, 'config_use_ecosysTracers_surface_value', config_use_ecosysTracers_surface_value)
      call mpas_pool_add_config(configPool, 'config_use_ecosysTracers_river_inputs_from_coupler', &
config_use_ecosysTracers_river_inputs_from_coupler)
      call mpas_pool_add_config(configPool, 'config_use_ecosysTracers_sea_ice_coupling', config_use_ecosysTracers_sea_ice_coupling)
      call mpas_pool_add_config(configPool, 'config_ecosysTracers_diagnostic_fields_level1', &
config_ecosysTracers_diagnostic_fields_level1)
      call mpas_pool_add_config(configPool, 'config_ecosysTracers_diagnostic_fields_level2', &
config_ecosysTracers_diagnostic_fields_level2)
      call mpas_pool_add_config(configPool, 'config_ecosysTracers_diagnostic_fields_level3', &
config_ecosysTracers_diagnostic_fields_level3)
      call mpas_pool_add_config(configPool, 'config_ecosysTracers_diagnostic_fields_level4', &
config_ecosysTracers_diagnostic_fields_level4)
      call mpas_pool_add_config(configPool, 'config_ecosysTracers_diagnostic_fields_level5', &
config_ecosysTracers_diagnostic_fields_level5)

   end subroutine ocn_setup_nmlrec_tracer_forcing_ecosysTracers


   subroutine ocn_setup_nmlrec_tracer_forcing_DMSTracers(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_use_DMSTracers = .false.
      logical :: config_use_DMSTracers_surface_bulk_forcing = .false.
      logical :: config_use_DMSTracers_surface_restoring = .false.
      logical :: config_use_DMSTracers_interior_restoring = .false.
      logical :: config_use_DMSTracers_exponential_decay = .false.
      logical :: config_use_DMSTracers_idealAge_forcing = .false.
      logical :: config_use_DMSTracers_ttd_forcing = .false.
      logical :: config_use_DMSTracers_surface_value = .false.
      logical :: config_use_DMSTracers_sea_ice_coupling = .false.

      namelist /tracer_forcing_DMSTracers/ &
         config_use_DMSTracers, &
         config_use_DMSTracers_surface_bulk_forcing, &
         config_use_DMSTracers_surface_restoring, &
         config_use_DMSTracers_interior_restoring, &
         config_use_DMSTracers_exponential_decay, &
         config_use_DMSTracers_idealAge_forcing, &
         config_use_DMSTracers_ttd_forcing, &
         config_use_DMSTracers_surface_value, &
         config_use_DMSTracers_sea_ice_coupling
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, tracer_forcing_DMSTracers, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_use_DMSTracers)
         call mpas_dmpar_bcast_logical(dminfo, config_use_DMSTracers_surface_bulk_forcing)
         call mpas_dmpar_bcast_logical(dminfo, config_use_DMSTracers_surface_restoring)
         call mpas_dmpar_bcast_logical(dminfo, config_use_DMSTracers_interior_restoring)
         call mpas_dmpar_bcast_logical(dminfo, config_use_DMSTracers_exponential_decay)
         call mpas_dmpar_bcast_logical(dminfo, config_use_DMSTracers_idealAge_forcing)
         call mpas_dmpar_bcast_logical(dminfo, config_use_DMSTracers_ttd_forcing)
         call mpas_dmpar_bcast_logical(dminfo, config_use_DMSTracers_surface_value)
         call mpas_dmpar_bcast_logical(dminfo, config_use_DMSTracers_sea_ice_coupling)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record tracer_forcing_DMSTracers')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_use_DMSTracers = $l', logicArgs=(/config_use_DMSTracers/))
            call mpas_log_write('        config_use_DMSTracers_surface_bulk_forcing = $l', &
logicArgs=(/config_use_DMSTracers_surface_bulk_forcing/))
            call mpas_log_write('        config_use_DMSTracers_surface_restoring = $l', &
logicArgs=(/config_use_DMSTracers_surface_restoring/))
            call mpas_log_write('        config_use_DMSTracers_interior_restoring = $l', &
logicArgs=(/config_use_DMSTracers_interior_restoring/))
            call mpas_log_write('        config_use_DMSTracers_exponential_decay = $l', &
logicArgs=(/config_use_DMSTracers_exponential_decay/))
            call mpas_log_write('        config_use_DMSTracers_idealAge_forcing = $l', &
logicArgs=(/config_use_DMSTracers_idealAge_forcing/))
            call mpas_log_write('        config_use_DMSTracers_ttd_forcing = $l', logicArgs=(/config_use_DMSTracers_ttd_forcing/))
            call mpas_log_write('        config_use_DMSTracers_surface_value = $l', &
logicArgs=(/config_use_DMSTracers_surface_value/))
            call mpas_log_write('        config_use_DMSTracers_sea_ice_coupling = $l', &
logicArgs=(/config_use_DMSTracers_sea_ice_coupling/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record tracer_forcing_DMSTracers.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_use_DMSTracers', config_use_DMSTracers)
      call mpas_pool_add_config(configPool, 'config_use_DMSTracers_surface_bulk_forcing', &
config_use_DMSTracers_surface_bulk_forcing)
      call mpas_pool_add_config(configPool, 'config_use_DMSTracers_surface_restoring', config_use_DMSTracers_surface_restoring)
      call mpas_pool_add_config(configPool, 'config_use_DMSTracers_interior_restoring', config_use_DMSTracers_interior_restoring)
      call mpas_pool_add_config(configPool, 'config_use_DMSTracers_exponential_decay', config_use_DMSTracers_exponential_decay)
      call mpas_pool_add_config(configPool, 'config_use_DMSTracers_idealAge_forcing', config_use_DMSTracers_idealAge_forcing)
      call mpas_pool_add_config(configPool, 'config_use_DMSTracers_ttd_forcing', config_use_DMSTracers_ttd_forcing)
      call mpas_pool_add_config(configPool, 'config_use_DMSTracers_surface_value', config_use_DMSTracers_surface_value)
      call mpas_pool_add_config(configPool, 'config_use_DMSTracers_sea_ice_coupling', config_use_DMSTracers_sea_ice_coupling)

   end subroutine ocn_setup_nmlrec_tracer_forcing_DMSTracers


   subroutine ocn_setup_nmlrec_tracer_forcing_MacroMoleculesTracers(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_use_MacroMoleculesTracers = .false.
      logical :: config_use_MacroMoleculesTracers_surface_bulk_forcing = .false.
      logical :: config_use_MacroMoleculesTracers_surface_restoring = .false.
      logical :: config_use_MacroMoleculesTracers_interior_restoring = .false.
      logical :: config_use_MacroMoleculesTracers_exponential_decay = .false.
      logical :: config_use_MacroMoleculesTracers_idealAge_forcing = .false.
      logical :: config_use_MacroMoleculesTracers_ttd_forcing = .false.
      logical :: config_use_MacroMoleculesTracers_surface_value = .false.
      logical :: config_use_MacroMoleculesTracers_sea_ice_coupling = .false.

      namelist /tracer_forcing_MacroMoleculesTracers/ &
         config_use_MacroMoleculesTracers, &
         config_use_MacroMoleculesTracers_surface_bulk_forcing, &
         config_use_MacroMoleculesTracers_surface_restoring, &
         config_use_MacroMoleculesTracers_interior_restoring, &
         config_use_MacroMoleculesTracers_exponential_decay, &
         config_use_MacroMoleculesTracers_idealAge_forcing, &
         config_use_MacroMoleculesTracers_ttd_forcing, &
         config_use_MacroMoleculesTracers_surface_value, &
         config_use_MacroMoleculesTracers_sea_ice_coupling
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, tracer_forcing_MacroMoleculesTracers, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_use_MacroMoleculesTracers)
         call mpas_dmpar_bcast_logical(dminfo, config_use_MacroMoleculesTracers_surface_bulk_forcing)
         call mpas_dmpar_bcast_logical(dminfo, config_use_MacroMoleculesTracers_surface_restoring)
         call mpas_dmpar_bcast_logical(dminfo, config_use_MacroMoleculesTracers_interior_restoring)
         call mpas_dmpar_bcast_logical(dminfo, config_use_MacroMoleculesTracers_exponential_decay)
         call mpas_dmpar_bcast_logical(dminfo, config_use_MacroMoleculesTracers_idealAge_forcing)
         call mpas_dmpar_bcast_logical(dminfo, config_use_MacroMoleculesTracers_ttd_forcing)
         call mpas_dmpar_bcast_logical(dminfo, config_use_MacroMoleculesTracers_surface_value)
         call mpas_dmpar_bcast_logical(dminfo, config_use_MacroMoleculesTracers_sea_ice_coupling)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record '&
//'tracer_forcing_MacroMoleculesTracers')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_use_MacroMoleculesTracers = $l', logicArgs=(/config_use_MacroMoleculesTracers/))
            call mpas_log_write('        config_use_MacroMoleculesTracers_surface_bulk_forcing = $l', &
logicArgs=(/config_use_MacroMoleculesTracers_surface_bulk_forcing/))
            call mpas_log_write('        config_use_MacroMoleculesTracers_surface_restoring = $l', &
logicArgs=(/config_use_MacroMoleculesTracers_surface_restoring/))
            call mpas_log_write('        config_use_MacroMoleculesTracers_interior_restoring = $l', &
logicArgs=(/config_use_MacroMoleculesTracers_interior_restoring/))
            call mpas_log_write('        config_use_MacroMoleculesTracers_exponential_decay = $l', &
logicArgs=(/config_use_MacroMoleculesTracers_exponential_decay/))
            call mpas_log_write('        config_use_MacroMoleculesTracers_idealAge_forcing = $l', &
logicArgs=(/config_use_MacroMoleculesTracers_idealAge_forcing/))
            call mpas_log_write('        config_use_MacroMoleculesTracers_ttd_forcing = $l', &
logicArgs=(/config_use_MacroMoleculesTracers_ttd_forcing/))
            call mpas_log_write('        config_use_MacroMoleculesTracers_surface_value = $l', &
logicArgs=(/config_use_MacroMoleculesTracers_surface_value/))
            call mpas_log_write('        config_use_MacroMoleculesTracers_sea_ice_coupling = $l', &
logicArgs=(/config_use_MacroMoleculesTracers_sea_ice_coupling/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record tracer_forcing_MacroMoleculesTracers.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_use_MacroMoleculesTracers', config_use_MacroMoleculesTracers)
      call mpas_pool_add_config(configPool, 'config_use_MacroMoleculesTracers_surface_bulk_forcing', &
config_use_MacroMoleculesTracers_surface_bulk_forcing)
      call mpas_pool_add_config(configPool, 'config_use_MacroMoleculesTracers_surface_restoring', &
config_use_MacroMoleculesTracers_surface_restoring)
      call mpas_pool_add_config(configPool, 'config_use_MacroMoleculesTracers_interior_restoring', &
config_use_MacroMoleculesTracers_interior_restoring)
      call mpas_pool_add_config(configPool, 'config_use_MacroMoleculesTracers_exponential_decay', &
config_use_MacroMoleculesTracers_exponential_decay)
      call mpas_pool_add_config(configPool, 'config_use_MacroMoleculesTracers_idealAge_forcing', &
config_use_MacroMoleculesTracers_idealAge_forcing)
      call mpas_pool_add_config(configPool, 'config_use_MacroMoleculesTracers_ttd_forcing', &
config_use_MacroMoleculesTracers_ttd_forcing)
      call mpas_pool_add_config(configPool, 'config_use_MacroMoleculesTracers_surface_value', &
config_use_MacroMoleculesTracers_surface_value)
      call mpas_pool_add_config(configPool, 'config_use_MacroMoleculesTracers_sea_ice_coupling', &
config_use_MacroMoleculesTracers_sea_ice_coupling)

   end subroutine ocn_setup_nmlrec_tracer_forcing_MacroMoleculesTracers


   subroutine ocn_setup_nmlrec_tracer_forcing_idealAgeTracers(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_use_idealAgeTracers = .false.
      logical :: config_use_idealAgeTracers_surface_bulk_forcing = .false.
      logical :: config_use_idealAgeTracers_surface_restoring = .false.
      logical :: config_use_idealAgeTracers_interior_restoring = .false.
      logical :: config_use_idealAgeTracers_exponential_decay = .false.
      logical :: config_use_idealAgeTracers_idealAge_forcing = .true.
      logical :: config_use_idealAgeTracers_ttd_forcing = .false.

      namelist /tracer_forcing_idealAgeTracers/ &
         config_use_idealAgeTracers, &
         config_use_idealAgeTracers_surface_bulk_forcing, &
         config_use_idealAgeTracers_surface_restoring, &
         config_use_idealAgeTracers_interior_restoring, &
         config_use_idealAgeTracers_exponential_decay, &
         config_use_idealAgeTracers_idealAge_forcing, &
         config_use_idealAgeTracers_ttd_forcing
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, tracer_forcing_idealAgeTracers, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_use_idealAgeTracers)
         call mpas_dmpar_bcast_logical(dminfo, config_use_idealAgeTracers_surface_bulk_forcing)
         call mpas_dmpar_bcast_logical(dminfo, config_use_idealAgeTracers_surface_restoring)
         call mpas_dmpar_bcast_logical(dminfo, config_use_idealAgeTracers_interior_restoring)
         call mpas_dmpar_bcast_logical(dminfo, config_use_idealAgeTracers_exponential_decay)
         call mpas_dmpar_bcast_logical(dminfo, config_use_idealAgeTracers_idealAge_forcing)
         call mpas_dmpar_bcast_logical(dminfo, config_use_idealAgeTracers_ttd_forcing)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record tracer_forcing_idealAgeTracers')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_use_idealAgeTracers = $l', logicArgs=(/config_use_idealAgeTracers/))
            call mpas_log_write('        config_use_idealAgeTracers_surface_bulk_forcing = $l', &
logicArgs=(/config_use_idealAgeTracers_surface_bulk_forcing/))
            call mpas_log_write('        config_use_idealAgeTracers_surface_restoring = $l', &
logicArgs=(/config_use_idealAgeTracers_surface_restoring/))
            call mpas_log_write('        config_use_idealAgeTracers_interior_restoring = $l', &
logicArgs=(/config_use_idealAgeTracers_interior_restoring/))
            call mpas_log_write('        config_use_idealAgeTracers_exponential_decay = $l', &
logicArgs=(/config_use_idealAgeTracers_exponential_decay/))
            call mpas_log_write('        config_use_idealAgeTracers_idealAge_forcing = $l', &
logicArgs=(/config_use_idealAgeTracers_idealAge_forcing/))
            call mpas_log_write('        config_use_idealAgeTracers_ttd_forcing = $l', &
logicArgs=(/config_use_idealAgeTracers_ttd_forcing/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record tracer_forcing_idealAgeTracers.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_use_idealAgeTracers', config_use_idealAgeTracers)
      call mpas_pool_add_config(configPool, 'config_use_idealAgeTracers_surface_bulk_forcing', &
config_use_idealAgeTracers_surface_bulk_forcing)
      call mpas_pool_add_config(configPool, 'config_use_idealAgeTracers_surface_restoring', &
config_use_idealAgeTracers_surface_restoring)
      call mpas_pool_add_config(configPool, 'config_use_idealAgeTracers_interior_restoring', &
config_use_idealAgeTracers_interior_restoring)
      call mpas_pool_add_config(configPool, 'config_use_idealAgeTracers_exponential_decay', &
config_use_idealAgeTracers_exponential_decay)
      call mpas_pool_add_config(configPool, 'config_use_idealAgeTracers_idealAge_forcing', &
config_use_idealAgeTracers_idealAge_forcing)
      call mpas_pool_add_config(configPool, 'config_use_idealAgeTracers_ttd_forcing', config_use_idealAgeTracers_ttd_forcing)

   end subroutine ocn_setup_nmlrec_tracer_forcing_idealAgeTracers


   subroutine ocn_setup_nmlrec_tracer_forcing_CFCTracers(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_use_CFCTracers = .false.
      logical :: config_use_CFCTracers_surface_bulk_forcing = .false.
      logical :: config_use_CFCTracers_surface_restoring = .false.
      logical :: config_use_CFCTracers_interior_restoring = .false.
      logical :: config_use_CFCTracers_exponential_decay = .false.
      logical :: config_use_CFCTracers_idealAge_forcing = .false.
      logical :: config_use_CFCTracers_ttd_forcing = .false.
      logical :: config_use_CFC11 = .true.
      logical :: config_use_CFC12 = .true.

      namelist /tracer_forcing_CFCTracers/ &
         config_use_CFCTracers, &
         config_use_CFCTracers_surface_bulk_forcing, &
         config_use_CFCTracers_surface_restoring, &
         config_use_CFCTracers_interior_restoring, &
         config_use_CFCTracers_exponential_decay, &
         config_use_CFCTracers_idealAge_forcing, &
         config_use_CFCTracers_ttd_forcing, &
         config_use_CFC11, &
         config_use_CFC12
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, tracer_forcing_CFCTracers, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_use_CFCTracers)
         call mpas_dmpar_bcast_logical(dminfo, config_use_CFCTracers_surface_bulk_forcing)
         call mpas_dmpar_bcast_logical(dminfo, config_use_CFCTracers_surface_restoring)
         call mpas_dmpar_bcast_logical(dminfo, config_use_CFCTracers_interior_restoring)
         call mpas_dmpar_bcast_logical(dminfo, config_use_CFCTracers_exponential_decay)
         call mpas_dmpar_bcast_logical(dminfo, config_use_CFCTracers_idealAge_forcing)
         call mpas_dmpar_bcast_logical(dminfo, config_use_CFCTracers_ttd_forcing)
         call mpas_dmpar_bcast_logical(dminfo, config_use_CFC11)
         call mpas_dmpar_bcast_logical(dminfo, config_use_CFC12)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record tracer_forcing_CFCTracers')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_use_CFCTracers = $l', logicArgs=(/config_use_CFCTracers/))
            call mpas_log_write('        config_use_CFCTracers_surface_bulk_forcing = $l', &
logicArgs=(/config_use_CFCTracers_surface_bulk_forcing/))
            call mpas_log_write('        config_use_CFCTracers_surface_restoring = $l', &
logicArgs=(/config_use_CFCTracers_surface_restoring/))
            call mpas_log_write('        config_use_CFCTracers_interior_restoring = $l', &
logicArgs=(/config_use_CFCTracers_interior_restoring/))
            call mpas_log_write('        config_use_CFCTracers_exponential_decay = $l', &
logicArgs=(/config_use_CFCTracers_exponential_decay/))
            call mpas_log_write('        config_use_CFCTracers_idealAge_forcing = $l', &
logicArgs=(/config_use_CFCTracers_idealAge_forcing/))
            call mpas_log_write('        config_use_CFCTracers_ttd_forcing = $l', logicArgs=(/config_use_CFCTracers_ttd_forcing/))
            call mpas_log_write('        config_use_CFC11 = $l', logicArgs=(/config_use_CFC11/))
            call mpas_log_write('        config_use_CFC12 = $l', logicArgs=(/config_use_CFC12/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record tracer_forcing_CFCTracers.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_use_CFCTracers', config_use_CFCTracers)
      call mpas_pool_add_config(configPool, 'config_use_CFCTracers_surface_bulk_forcing', &
config_use_CFCTracers_surface_bulk_forcing)
      call mpas_pool_add_config(configPool, 'config_use_CFCTracers_surface_restoring', config_use_CFCTracers_surface_restoring)
      call mpas_pool_add_config(configPool, 'config_use_CFCTracers_interior_restoring', config_use_CFCTracers_interior_restoring)
      call mpas_pool_add_config(configPool, 'config_use_CFCTracers_exponential_decay', config_use_CFCTracers_exponential_decay)
      call mpas_pool_add_config(configPool, 'config_use_CFCTracers_idealAge_forcing', config_use_CFCTracers_idealAge_forcing)
      call mpas_pool_add_config(configPool, 'config_use_CFCTracers_ttd_forcing', config_use_CFCTracers_ttd_forcing)
      call mpas_pool_add_config(configPool, 'config_use_CFC11', config_use_CFC11)
      call mpas_pool_add_config(configPool, 'config_use_CFC12', config_use_CFC12)

   end subroutine ocn_setup_nmlrec_tracer_forcing_CFCTracers


   subroutine ocn_setup_nmlrec_AM_globalStats(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_AM_globalStats_enable = .false.
      character (len=StrKIND) :: config_AM_globalStats_compute_interval = 'output_interval'
      logical :: config_AM_globalStats_compute_on_startup = .false.
      logical :: config_AM_globalStats_write_on_startup = .false.
      logical :: config_AM_globalStats_text_file = .false.
      character (len=StrKIND) :: config_AM_globalStats_directory = 'analysis_members'
      character (len=StrKIND) :: config_AM_globalStats_output_stream = 'globalStatsOutput'

      namelist /AM_globalStats/ &
         config_AM_globalStats_enable, &
         config_AM_globalStats_compute_interval, &
         config_AM_globalStats_compute_on_startup, &
         config_AM_globalStats_write_on_startup, &
         config_AM_globalStats_text_file, &
         config_AM_globalStats_directory, &
         config_AM_globalStats_output_stream
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, AM_globalStats, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_AM_globalStats_enable)
         call mpas_dmpar_bcast_char(dminfo, config_AM_globalStats_compute_interval)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_globalStats_compute_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_globalStats_write_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_globalStats_text_file)
         call mpas_dmpar_bcast_char(dminfo, config_AM_globalStats_directory)
         call mpas_dmpar_bcast_char(dminfo, config_AM_globalStats_output_stream)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record AM_globalStats')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_AM_globalStats_enable = $l', logicArgs=(/config_AM_globalStats_enable/))
            call mpas_log_write('        config_AM_globalStats_compute_interval = '&
//''//mpas_log_escape_dollars(config_AM_globalStats_compute_interval))
            call mpas_log_write('        config_AM_globalStats_compute_on_startup = $l', &
logicArgs=(/config_AM_globalStats_compute_on_startup/))
            call mpas_log_write('        config_AM_globalStats_write_on_startup = $l', &
logicArgs=(/config_AM_globalStats_write_on_startup/))
            call mpas_log_write('        config_AM_globalStats_text_file = $l', logicArgs=(/config_AM_globalStats_text_file/))
            call mpas_log_write('        config_AM_globalStats_directory = '&
//''//mpas_log_escape_dollars(config_AM_globalStats_directory))
            call mpas_log_write('        config_AM_globalStats_output_stream = '&
//''//mpas_log_escape_dollars(config_AM_globalStats_output_stream))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record AM_globalStats.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_AM_globalStats_enable', config_AM_globalStats_enable)
      call mpas_pool_add_config(configPool, 'config_AM_globalStats_compute_interval', config_AM_globalStats_compute_interval)
      call mpas_pool_add_config(configPool, 'config_AM_globalStats_compute_on_startup', config_AM_globalStats_compute_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_globalStats_write_on_startup', config_AM_globalStats_write_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_globalStats_text_file', config_AM_globalStats_text_file)
      call mpas_pool_add_config(configPool, 'config_AM_globalStats_directory', config_AM_globalStats_directory)
      call mpas_pool_add_config(configPool, 'config_AM_globalStats_output_stream', config_AM_globalStats_output_stream)

   end subroutine ocn_setup_nmlrec_AM_globalStats


   subroutine ocn_setup_nmlrec_AM_surfaceAreaWeightedAverages(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_AM_surfaceAreaWeightedAverages_enable = .false.
      logical :: config_AM_surfaceAreaWeightedAverages_compute_on_startup = .true.
      logical :: config_AM_surfaceAreaWeightedAverages_write_on_startup = .true.
      character (len=StrKIND) :: config_AM_surfaceAreaWeightedAverages_compute_interval = 'output_interval'
      character (len=StrKIND) :: config_AM_surfaceAreaWeightedAverages_output_stream = 'surfaceAreaWeightedAveragesOutput'

      namelist /AM_surfaceAreaWeightedAverages/ &
         config_AM_surfaceAreaWeightedAverages_enable, &
         config_AM_surfaceAreaWeightedAverages_compute_on_startup, &
         config_AM_surfaceAreaWeightedAverages_write_on_startup, &
         config_AM_surfaceAreaWeightedAverages_compute_interval, &
         config_AM_surfaceAreaWeightedAverages_output_stream
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, AM_surfaceAreaWeightedAverages, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_AM_surfaceAreaWeightedAverages_enable)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_surfaceAreaWeightedAverages_compute_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_surfaceAreaWeightedAverages_write_on_startup)
         call mpas_dmpar_bcast_char(dminfo, config_AM_surfaceAreaWeightedAverages_compute_interval)
         call mpas_dmpar_bcast_char(dminfo, config_AM_surfaceAreaWeightedAverages_output_stream)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record AM_surfaceAreaWeightedAverages')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_AM_surfaceAreaWeightedAverages_enable = $l', &
logicArgs=(/config_AM_surfaceAreaWeightedAverages_enable/))
            call mpas_log_write('        config_AM_surfaceAreaWeightedAverages_compute_on_startup = $l', &
logicArgs=(/config_AM_surfaceAreaWeightedAverages_compute_on_startup/))
            call mpas_log_write('        config_AM_surfaceAreaWeightedAverages_write_on_startup = $l', &
logicArgs=(/config_AM_surfaceAreaWeightedAverages_write_on_startup/))
            call mpas_log_write('        config_AM_surfaceAreaWeightedAverages_compute_interval = '&
//''//mpas_log_escape_dollars(config_AM_surfaceAreaWeightedAverages_compute_interval))
            call mpas_log_write('        config_AM_surfaceAreaWeightedAverages_output_stream = '&
//''//mpas_log_escape_dollars(config_AM_surfaceAreaWeightedAverages_output_stream))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record AM_surfaceAreaWeightedAverages.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_AM_surfaceAreaWeightedAverages_enable', &
config_AM_surfaceAreaWeightedAverages_enable)
      call mpas_pool_add_config(configPool, 'config_AM_surfaceAreaWeightedAverages_compute_on_startup', &
config_AM_surfaceAreaWeightedAverages_compute_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_surfaceAreaWeightedAverages_write_on_startup', &
config_AM_surfaceAreaWeightedAverages_write_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_surfaceAreaWeightedAverages_compute_interval', &
config_AM_surfaceAreaWeightedAverages_compute_interval)
      call mpas_pool_add_config(configPool, 'config_AM_surfaceAreaWeightedAverages_output_stream', &
config_AM_surfaceAreaWeightedAverages_output_stream)

   end subroutine ocn_setup_nmlrec_AM_surfaceAreaWeightedAverages


   subroutine ocn_setup_nmlrec_AM_waterMassCensus(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_AM_waterMassCensus_enable = .false.
      character (len=StrKIND) :: config_AM_waterMassCensus_compute_interval = 'output_interval'
      character (len=StrKIND) :: config_AM_waterMassCensus_output_stream = 'waterMassCensusOutput'
      logical :: config_AM_waterMassCensus_compute_on_startup = .false.
      logical :: config_AM_waterMassCensus_write_on_startup = .false.
      real (kind=RKIND) :: config_AM_waterMassCensus_minTemperature = -2.000000
      real (kind=RKIND) :: config_AM_waterMassCensus_maxTemperature = 30.000000
      real (kind=RKIND) :: config_AM_waterMassCensus_minSalinity = 32.000000
      real (kind=RKIND) :: config_AM_waterMassCensus_maxSalinity = 37.000000
      logical :: config_AM_waterMassCensus_compute_predefined_regions = .true.
      character (len=StrKIND) :: config_AM_waterMassCensus_region_group = ''

      namelist /AM_waterMassCensus/ &
         config_AM_waterMassCensus_enable, &
         config_AM_waterMassCensus_compute_interval, &
         config_AM_waterMassCensus_output_stream, &
         config_AM_waterMassCensus_compute_on_startup, &
         config_AM_waterMassCensus_write_on_startup, &
         config_AM_waterMassCensus_minTemperature, &
         config_AM_waterMassCensus_maxTemperature, &
         config_AM_waterMassCensus_minSalinity, &
         config_AM_waterMassCensus_maxSalinity, &
         config_AM_waterMassCensus_compute_predefined_regions, &
         config_AM_waterMassCensus_region_group
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, AM_waterMassCensus, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_AM_waterMassCensus_enable)
         call mpas_dmpar_bcast_char(dminfo, config_AM_waterMassCensus_compute_interval)
         call mpas_dmpar_bcast_char(dminfo, config_AM_waterMassCensus_output_stream)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_waterMassCensus_compute_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_waterMassCensus_write_on_startup)
         call mpas_dmpar_bcast_real(dminfo, config_AM_waterMassCensus_minTemperature)
         call mpas_dmpar_bcast_real(dminfo, config_AM_waterMassCensus_maxTemperature)
         call mpas_dmpar_bcast_real(dminfo, config_AM_waterMassCensus_minSalinity)
         call mpas_dmpar_bcast_real(dminfo, config_AM_waterMassCensus_maxSalinity)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_waterMassCensus_compute_predefined_regions)
         call mpas_dmpar_bcast_char(dminfo, config_AM_waterMassCensus_region_group)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record AM_waterMassCensus')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_AM_waterMassCensus_enable = $l', logicArgs=(/config_AM_waterMassCensus_enable/))
            call mpas_log_write('        config_AM_waterMassCensus_compute_interval = '&
//''//mpas_log_escape_dollars(config_AM_waterMassCensus_compute_interval))
            call mpas_log_write('        config_AM_waterMassCensus_output_stream = '&
//''//mpas_log_escape_dollars(config_AM_waterMassCensus_output_stream))
            call mpas_log_write('        config_AM_waterMassCensus_compute_on_startup = $l', &
logicArgs=(/config_AM_waterMassCensus_compute_on_startup/))
            call mpas_log_write('        config_AM_waterMassCensus_write_on_startup = $l', &
logicArgs=(/config_AM_waterMassCensus_write_on_startup/))
            call mpas_log_write('        config_AM_waterMassCensus_minTemperature = $r', &
realArgs=(/config_AM_waterMassCensus_minTemperature/))
            call mpas_log_write('        config_AM_waterMassCensus_maxTemperature = $r', &
realArgs=(/config_AM_waterMassCensus_maxTemperature/))
            call mpas_log_write('        config_AM_waterMassCensus_minSalinity = $r', &
realArgs=(/config_AM_waterMassCensus_minSalinity/))
            call mpas_log_write('        config_AM_waterMassCensus_maxSalinity = $r', &
realArgs=(/config_AM_waterMassCensus_maxSalinity/))
            call mpas_log_write('        config_AM_waterMassCensus_compute_predefined_regions = $l', &
logicArgs=(/config_AM_waterMassCensus_compute_predefined_regions/))
            call mpas_log_write('        config_AM_waterMassCensus_region_group = '&
//''//mpas_log_escape_dollars(config_AM_waterMassCensus_region_group))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record AM_waterMassCensus.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_AM_waterMassCensus_enable', config_AM_waterMassCensus_enable)
      call mpas_pool_add_config(configPool, 'config_AM_waterMassCensus_compute_interval', &
config_AM_waterMassCensus_compute_interval)
      call mpas_pool_add_config(configPool, 'config_AM_waterMassCensus_output_stream', config_AM_waterMassCensus_output_stream)
      call mpas_pool_add_config(configPool, 'config_AM_waterMassCensus_compute_on_startup', &
config_AM_waterMassCensus_compute_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_waterMassCensus_write_on_startup', &
config_AM_waterMassCensus_write_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_waterMassCensus_minTemperature', config_AM_waterMassCensus_minTemperature)
      call mpas_pool_add_config(configPool, 'config_AM_waterMassCensus_maxTemperature', config_AM_waterMassCensus_maxTemperature)
      call mpas_pool_add_config(configPool, 'config_AM_waterMassCensus_minSalinity', config_AM_waterMassCensus_minSalinity)
      call mpas_pool_add_config(configPool, 'config_AM_waterMassCensus_maxSalinity', config_AM_waterMassCensus_maxSalinity)
      call mpas_pool_add_config(configPool, 'config_AM_waterMassCensus_compute_predefined_regions', &
config_AM_waterMassCensus_compute_predefined_regions)
      call mpas_pool_add_config(configPool, 'config_AM_waterMassCensus_region_group', config_AM_waterMassCensus_region_group)

   end subroutine ocn_setup_nmlrec_AM_waterMassCensus


   subroutine ocn_setup_nmlrec_AM_layerVolumeWeightedAverage(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_AM_layerVolumeWeightedAverage_enable = .false.
      character (len=StrKIND) :: config_AM_layerVolumeWeightedAverage_compute_interval = 'output_interval'
      logical :: config_AM_layerVolumeWeightedAverage_compute_on_startup = .false.
      logical :: config_AM_layerVolumeWeightedAverage_write_on_startup = .false.
      character (len=StrKIND) :: config_AM_layerVolumeWeightedAverage_output_stream = 'layerVolumeWeightedAverageOutput'

      namelist /AM_layerVolumeWeightedAverage/ &
         config_AM_layerVolumeWeightedAverage_enable, &
         config_AM_layerVolumeWeightedAverage_compute_interval, &
         config_AM_layerVolumeWeightedAverage_compute_on_startup, &
         config_AM_layerVolumeWeightedAverage_write_on_startup, &
         config_AM_layerVolumeWeightedAverage_output_stream
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, AM_layerVolumeWeightedAverage, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_AM_layerVolumeWeightedAverage_enable)
         call mpas_dmpar_bcast_char(dminfo, config_AM_layerVolumeWeightedAverage_compute_interval)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_layerVolumeWeightedAverage_compute_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_layerVolumeWeightedAverage_write_on_startup)
         call mpas_dmpar_bcast_char(dminfo, config_AM_layerVolumeWeightedAverage_output_stream)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record AM_layerVolumeWeightedAverage')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_AM_layerVolumeWeightedAverage_enable = $l', &
logicArgs=(/config_AM_layerVolumeWeightedAverage_enable/))
            call mpas_log_write('        config_AM_layerVolumeWeightedAverage_compute_interval = '&
//''//mpas_log_escape_dollars(config_AM_layerVolumeWeightedAverage_compute_interval))
            call mpas_log_write('        config_AM_layerVolumeWeightedAverage_compute_on_startup = $l', &
logicArgs=(/config_AM_layerVolumeWeightedAverage_compute_on_startup/))
            call mpas_log_write('        config_AM_layerVolumeWeightedAverage_write_on_startup = $l', &
logicArgs=(/config_AM_layerVolumeWeightedAverage_write_on_startup/))
            call mpas_log_write('        config_AM_layerVolumeWeightedAverage_output_stream = '&
//''//mpas_log_escape_dollars(config_AM_layerVolumeWeightedAverage_output_stream))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record AM_layerVolumeWeightedAverage.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_AM_layerVolumeWeightedAverage_enable', &
config_AM_layerVolumeWeightedAverage_enable)
      call mpas_pool_add_config(configPool, 'config_AM_layerVolumeWeightedAverage_compute_interval', &
config_AM_layerVolumeWeightedAverage_compute_interval)
      call mpas_pool_add_config(configPool, 'config_AM_layerVolumeWeightedAverage_compute_on_startup', &
config_AM_layerVolumeWeightedAverage_compute_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_layerVolumeWeightedAverage_write_on_startup', &
config_AM_layerVolumeWeightedAverage_write_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_layerVolumeWeightedAverage_output_stream', &
config_AM_layerVolumeWeightedAverage_output_stream)

   end subroutine ocn_setup_nmlrec_AM_layerVolumeWeightedAverage


   subroutine ocn_setup_nmlrec_AM_zonalMean(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_AM_zonalMean_enable = .false.
      logical :: config_AM_zonalMean_compute_on_startup = .true.
      logical :: config_AM_zonalMean_write_on_startup = .true.
      character (len=StrKIND) :: config_AM_zonalMean_compute_interval = 'output_interval'
      character (len=StrKIND) :: config_AM_zonalMean_output_stream = 'zonalMeanOutput'
      integer :: config_AM_zonalMean_num_bins = 180
      real (kind=RKIND) :: config_AM_zonalMean_min_bin = -9999999999999999455752309870428160.000000
      real (kind=RKIND) :: config_AM_zonalMean_max_bin = -9999999999999999455752309870428160.000000

      namelist /AM_zonalMean/ &
         config_AM_zonalMean_enable, &
         config_AM_zonalMean_compute_on_startup, &
         config_AM_zonalMean_write_on_startup, &
         config_AM_zonalMean_compute_interval, &
         config_AM_zonalMean_output_stream, &
         config_AM_zonalMean_num_bins, &
         config_AM_zonalMean_min_bin, &
         config_AM_zonalMean_max_bin
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, AM_zonalMean, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_AM_zonalMean_enable)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_zonalMean_compute_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_zonalMean_write_on_startup)
         call mpas_dmpar_bcast_char(dminfo, config_AM_zonalMean_compute_interval)
         call mpas_dmpar_bcast_char(dminfo, config_AM_zonalMean_output_stream)
         call mpas_dmpar_bcast_int(dminfo, config_AM_zonalMean_num_bins)
         call mpas_dmpar_bcast_real(dminfo, config_AM_zonalMean_min_bin)
         call mpas_dmpar_bcast_real(dminfo, config_AM_zonalMean_max_bin)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record AM_zonalMean')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_AM_zonalMean_enable = $l', logicArgs=(/config_AM_zonalMean_enable/))
            call mpas_log_write('        config_AM_zonalMean_compute_on_startup = $l', &
logicArgs=(/config_AM_zonalMean_compute_on_startup/))
            call mpas_log_write('        config_AM_zonalMean_write_on_startup = $l', &
logicArgs=(/config_AM_zonalMean_write_on_startup/))
            call mpas_log_write('        config_AM_zonalMean_compute_interval = '&
//''//mpas_log_escape_dollars(config_AM_zonalMean_compute_interval))
            call mpas_log_write('        config_AM_zonalMean_output_stream = '&
//''//mpas_log_escape_dollars(config_AM_zonalMean_output_stream))
            call mpas_log_write('        config_AM_zonalMean_num_bins = $i', intArgs=(/config_AM_zonalMean_num_bins/))
            call mpas_log_write('        config_AM_zonalMean_min_bin = $r', realArgs=(/config_AM_zonalMean_min_bin/))
            call mpas_log_write('        config_AM_zonalMean_max_bin = $r', realArgs=(/config_AM_zonalMean_max_bin/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record AM_zonalMean.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_AM_zonalMean_enable', config_AM_zonalMean_enable)
      call mpas_pool_add_config(configPool, 'config_AM_zonalMean_compute_on_startup', config_AM_zonalMean_compute_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_zonalMean_write_on_startup', config_AM_zonalMean_write_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_zonalMean_compute_interval', config_AM_zonalMean_compute_interval)
      call mpas_pool_add_config(configPool, 'config_AM_zonalMean_output_stream', config_AM_zonalMean_output_stream)
      call mpas_pool_add_config(configPool, 'config_AM_zonalMean_num_bins', config_AM_zonalMean_num_bins)
      call mpas_pool_add_config(configPool, 'config_AM_zonalMean_min_bin', config_AM_zonalMean_min_bin)
      call mpas_pool_add_config(configPool, 'config_AM_zonalMean_max_bin', config_AM_zonalMean_max_bin)

   end subroutine ocn_setup_nmlrec_AM_zonalMean


   subroutine ocn_setup_nmlrec_AM_okuboWeiss(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_AM_okuboWeiss_enable = .false.
      logical :: config_AM_okuboWeiss_compute_on_startup = .true.
      logical :: config_AM_okuboWeiss_write_on_startup = .true.
      character (len=StrKIND) :: config_AM_okuboWeiss_compute_interval = 'output_interval'
      character (len=StrKIND) :: config_AM_okuboWeiss_output_stream = 'okuboWeissOutput'
      character (len=StrKIND) :: config_AM_okuboWeiss_directory = 'analysis_members'
      real (kind=RKIND) :: config_AM_okuboWeiss_threshold_value = -0.200000
      real (kind=RKIND) :: config_AM_okuboWeiss_normalization = 0.000000
      real (kind=RKIND) :: config_AM_okuboWeiss_lambda2_normalization = 0.000000
      logical :: config_AM_okuboWeiss_use_lat_lon_coords = .true.
      logical :: config_AM_okuboWeiss_compute_eddy_census = .true.
      integer :: config_AM_okuboWeiss_eddy_min_cells = 20

      namelist /AM_okuboWeiss/ &
         config_AM_okuboWeiss_enable, &
         config_AM_okuboWeiss_compute_on_startup, &
         config_AM_okuboWeiss_write_on_startup, &
         config_AM_okuboWeiss_compute_interval, &
         config_AM_okuboWeiss_output_stream, &
         config_AM_okuboWeiss_directory, &
         config_AM_okuboWeiss_threshold_value, &
         config_AM_okuboWeiss_normalization, &
         config_AM_okuboWeiss_lambda2_normalization, &
         config_AM_okuboWeiss_use_lat_lon_coords, &
         config_AM_okuboWeiss_compute_eddy_census, &
         config_AM_okuboWeiss_eddy_min_cells
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, AM_okuboWeiss, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_AM_okuboWeiss_enable)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_okuboWeiss_compute_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_okuboWeiss_write_on_startup)
         call mpas_dmpar_bcast_char(dminfo, config_AM_okuboWeiss_compute_interval)
         call mpas_dmpar_bcast_char(dminfo, config_AM_okuboWeiss_output_stream)
         call mpas_dmpar_bcast_char(dminfo, config_AM_okuboWeiss_directory)
         call mpas_dmpar_bcast_real(dminfo, config_AM_okuboWeiss_threshold_value)
         call mpas_dmpar_bcast_real(dminfo, config_AM_okuboWeiss_normalization)
         call mpas_dmpar_bcast_real(dminfo, config_AM_okuboWeiss_lambda2_normalization)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_okuboWeiss_use_lat_lon_coords)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_okuboWeiss_compute_eddy_census)
         call mpas_dmpar_bcast_int(dminfo, config_AM_okuboWeiss_eddy_min_cells)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record AM_okuboWeiss')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_AM_okuboWeiss_enable = $l', logicArgs=(/config_AM_okuboWeiss_enable/))
            call mpas_log_write('        config_AM_okuboWeiss_compute_on_startup = $l', &
logicArgs=(/config_AM_okuboWeiss_compute_on_startup/))
            call mpas_log_write('        config_AM_okuboWeiss_write_on_startup = $l', &
logicArgs=(/config_AM_okuboWeiss_write_on_startup/))
            call mpas_log_write('        config_AM_okuboWeiss_compute_interval = '&
//''//mpas_log_escape_dollars(config_AM_okuboWeiss_compute_interval))
            call mpas_log_write('        config_AM_okuboWeiss_output_stream = '&
//''//mpas_log_escape_dollars(config_AM_okuboWeiss_output_stream))
            call mpas_log_write('        config_AM_okuboWeiss_directory = '&
//''//mpas_log_escape_dollars(config_AM_okuboWeiss_directory))
            call mpas_log_write('        config_AM_okuboWeiss_threshold_value = $r', &
realArgs=(/config_AM_okuboWeiss_threshold_value/))
            call mpas_log_write('        config_AM_okuboWeiss_normalization = $r', realArgs=(/config_AM_okuboWeiss_normalization/))
            call mpas_log_write('        config_AM_okuboWeiss_lambda2_normalization = $r', &
realArgs=(/config_AM_okuboWeiss_lambda2_normalization/))
            call mpas_log_write('        config_AM_okuboWeiss_use_lat_lon_coords = $l', &
logicArgs=(/config_AM_okuboWeiss_use_lat_lon_coords/))
            call mpas_log_write('        config_AM_okuboWeiss_compute_eddy_census = $l', &
logicArgs=(/config_AM_okuboWeiss_compute_eddy_census/))
            call mpas_log_write('        config_AM_okuboWeiss_eddy_min_cells = $i', &
intArgs=(/config_AM_okuboWeiss_eddy_min_cells/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record AM_okuboWeiss.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_AM_okuboWeiss_enable', config_AM_okuboWeiss_enable)
      call mpas_pool_add_config(configPool, 'config_AM_okuboWeiss_compute_on_startup', config_AM_okuboWeiss_compute_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_okuboWeiss_write_on_startup', config_AM_okuboWeiss_write_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_okuboWeiss_compute_interval', config_AM_okuboWeiss_compute_interval)
      call mpas_pool_add_config(configPool, 'config_AM_okuboWeiss_output_stream', config_AM_okuboWeiss_output_stream)
      call mpas_pool_add_config(configPool, 'config_AM_okuboWeiss_directory', config_AM_okuboWeiss_directory)
      call mpas_pool_add_config(configPool, 'config_AM_okuboWeiss_threshold_value', config_AM_okuboWeiss_threshold_value)
      call mpas_pool_add_config(configPool, 'config_AM_okuboWeiss_normalization', config_AM_okuboWeiss_normalization)
      call mpas_pool_add_config(configPool, 'config_AM_okuboWeiss_lambda2_normalization', &
config_AM_okuboWeiss_lambda2_normalization)
      call mpas_pool_add_config(configPool, 'config_AM_okuboWeiss_use_lat_lon_coords', config_AM_okuboWeiss_use_lat_lon_coords)
      call mpas_pool_add_config(configPool, 'config_AM_okuboWeiss_compute_eddy_census', config_AM_okuboWeiss_compute_eddy_census)
      call mpas_pool_add_config(configPool, 'config_AM_okuboWeiss_eddy_min_cells', config_AM_okuboWeiss_eddy_min_cells)

   end subroutine ocn_setup_nmlrec_AM_okuboWeiss


   subroutine ocn_setup_nmlrec_AM_meridionalHeatTransport(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_AM_meridionalHeatTransport_enable = .false.
      character (len=StrKIND) :: config_AM_meridionalHeatTransport_compute_interval = 'output_interval'
      logical :: config_AM_meridionalHeatTransport_compute_on_startup = .true.
      logical :: config_AM_meridionalHeatTransport_write_on_startup = .true.
      character (len=StrKIND) :: config_AM_meridionalHeatTransport_output_stream = 'meridionalHeatTransportOutput'
      integer :: config_AM_meridionalHeatTransport_num_bins = 180
      real (kind=RKIND) :: config_AM_meridionalHeatTransport_min_bin = -9999999999999999455752309870428160.000000
      real (kind=RKIND) :: config_AM_meridionalHeatTransport_max_bin = -9999999999999999455752309870428160.000000
      character (len=StrKIND) :: config_AM_meridionalHeatTransport_region_group = ''

      namelist /AM_meridionalHeatTransport/ &
         config_AM_meridionalHeatTransport_enable, &
         config_AM_meridionalHeatTransport_compute_interval, &
         config_AM_meridionalHeatTransport_compute_on_startup, &
         config_AM_meridionalHeatTransport_write_on_startup, &
         config_AM_meridionalHeatTransport_output_stream, &
         config_AM_meridionalHeatTransport_num_bins, &
         config_AM_meridionalHeatTransport_min_bin, &
         config_AM_meridionalHeatTransport_max_bin, &
         config_AM_meridionalHeatTransport_region_group
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, AM_meridionalHeatTransport, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_AM_meridionalHeatTransport_enable)
         call mpas_dmpar_bcast_char(dminfo, config_AM_meridionalHeatTransport_compute_interval)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_meridionalHeatTransport_compute_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_meridionalHeatTransport_write_on_startup)
         call mpas_dmpar_bcast_char(dminfo, config_AM_meridionalHeatTransport_output_stream)
         call mpas_dmpar_bcast_int(dminfo, config_AM_meridionalHeatTransport_num_bins)
         call mpas_dmpar_bcast_real(dminfo, config_AM_meridionalHeatTransport_min_bin)
         call mpas_dmpar_bcast_real(dminfo, config_AM_meridionalHeatTransport_max_bin)
         call mpas_dmpar_bcast_char(dminfo, config_AM_meridionalHeatTransport_region_group)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record AM_meridionalHeatTransport')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_AM_meridionalHeatTransport_enable = $l', &
logicArgs=(/config_AM_meridionalHeatTransport_enable/))
            call mpas_log_write('        config_AM_meridionalHeatTransport_compute_interval = '&
//''//mpas_log_escape_dollars(config_AM_meridionalHeatTransport_compute_interval))
            call mpas_log_write('        config_AM_meridionalHeatTransport_compute_on_startup = $l', &
logicArgs=(/config_AM_meridionalHeatTransport_compute_on_startup/))
            call mpas_log_write('        config_AM_meridionalHeatTransport_write_on_startup = $l', &
logicArgs=(/config_AM_meridionalHeatTransport_write_on_startup/))
            call mpas_log_write('        config_AM_meridionalHeatTransport_output_stream = '&
//''//mpas_log_escape_dollars(config_AM_meridionalHeatTransport_output_stream))
            call mpas_log_write('        config_AM_meridionalHeatTransport_num_bins = $i', &
intArgs=(/config_AM_meridionalHeatTransport_num_bins/))
            call mpas_log_write('        config_AM_meridionalHeatTransport_min_bin = $r', &
realArgs=(/config_AM_meridionalHeatTransport_min_bin/))
            call mpas_log_write('        config_AM_meridionalHeatTransport_max_bin = $r', &
realArgs=(/config_AM_meridionalHeatTransport_max_bin/))
            call mpas_log_write('        config_AM_meridionalHeatTransport_region_group = '&
//''//mpas_log_escape_dollars(config_AM_meridionalHeatTransport_region_group))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record AM_meridionalHeatTransport.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_AM_meridionalHeatTransport_enable', config_AM_meridionalHeatTransport_enable)
      call mpas_pool_add_config(configPool, 'config_AM_meridionalHeatTransport_compute_interval', &
config_AM_meridionalHeatTransport_compute_interval)
      call mpas_pool_add_config(configPool, 'config_AM_meridionalHeatTransport_compute_on_startup', &
config_AM_meridionalHeatTransport_compute_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_meridionalHeatTransport_write_on_startup', &
config_AM_meridionalHeatTransport_write_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_meridionalHeatTransport_output_stream', &
config_AM_meridionalHeatTransport_output_stream)
      call mpas_pool_add_config(configPool, 'config_AM_meridionalHeatTransport_num_bins', &
config_AM_meridionalHeatTransport_num_bins)
      call mpas_pool_add_config(configPool, 'config_AM_meridionalHeatTransport_min_bin', config_AM_meridionalHeatTransport_min_bin)
      call mpas_pool_add_config(configPool, 'config_AM_meridionalHeatTransport_max_bin', config_AM_meridionalHeatTransport_max_bin)
      call mpas_pool_add_config(configPool, 'config_AM_meridionalHeatTransport_region_group', &
config_AM_meridionalHeatTransport_region_group)

   end subroutine ocn_setup_nmlrec_AM_meridionalHeatTransport


   subroutine ocn_setup_nmlrec_AM_testComputeInterval(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_AM_testComputeInterval_enable = .false.
      character (len=StrKIND) :: config_AM_testComputeInterval_compute_interval = '00-00-01_00:00:00'
      logical :: config_AM_testComputeInterval_compute_on_startup = .true.
      logical :: config_AM_testComputeInterval_write_on_startup = .true.
      character (len=StrKIND) :: config_AM_testComputeInterval_output_stream = 'testComputeIntervalOutput'

      namelist /AM_testComputeInterval/ &
         config_AM_testComputeInterval_enable, &
         config_AM_testComputeInterval_compute_interval, &
         config_AM_testComputeInterval_compute_on_startup, &
         config_AM_testComputeInterval_write_on_startup, &
         config_AM_testComputeInterval_output_stream
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, AM_testComputeInterval, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_AM_testComputeInterval_enable)
         call mpas_dmpar_bcast_char(dminfo, config_AM_testComputeInterval_compute_interval)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_testComputeInterval_compute_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_testComputeInterval_write_on_startup)
         call mpas_dmpar_bcast_char(dminfo, config_AM_testComputeInterval_output_stream)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record AM_testComputeInterval')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_AM_testComputeInterval_enable = $l', &
logicArgs=(/config_AM_testComputeInterval_enable/))
            call mpas_log_write('        config_AM_testComputeInterval_compute_interval = '&
//''//mpas_log_escape_dollars(config_AM_testComputeInterval_compute_interval))
            call mpas_log_write('        config_AM_testComputeInterval_compute_on_startup = $l', &
logicArgs=(/config_AM_testComputeInterval_compute_on_startup/))
            call mpas_log_write('        config_AM_testComputeInterval_write_on_startup = $l', &
logicArgs=(/config_AM_testComputeInterval_write_on_startup/))
            call mpas_log_write('        config_AM_testComputeInterval_output_stream = '&
//''//mpas_log_escape_dollars(config_AM_testComputeInterval_output_stream))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record AM_testComputeInterval.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_AM_testComputeInterval_enable', config_AM_testComputeInterval_enable)
      call mpas_pool_add_config(configPool, 'config_AM_testComputeInterval_compute_interval', &
config_AM_testComputeInterval_compute_interval)
      call mpas_pool_add_config(configPool, 'config_AM_testComputeInterval_compute_on_startup', &
config_AM_testComputeInterval_compute_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_testComputeInterval_write_on_startup', &
config_AM_testComputeInterval_write_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_testComputeInterval_output_stream', &
config_AM_testComputeInterval_output_stream)

   end subroutine ocn_setup_nmlrec_AM_testComputeInterval


   subroutine ocn_setup_nmlrec_AM_highFrequencyOutput(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_AM_highFrequencyOutput_enable = .false.
      character (len=StrKIND) :: config_AM_highFrequencyOutput_compute_interval = 'output_interval'
      character (len=StrKIND) :: config_AM_highFrequencyOutput_output_stream = 'highFrequencyOutput'
      logical :: config_AM_highFrequencyOutput_compute_on_startup = .true.
      logical :: config_AM_highFrequencyOutput_write_on_startup = .true.

      namelist /AM_highFrequencyOutput/ &
         config_AM_highFrequencyOutput_enable, &
         config_AM_highFrequencyOutput_compute_interval, &
         config_AM_highFrequencyOutput_output_stream, &
         config_AM_highFrequencyOutput_compute_on_startup, &
         config_AM_highFrequencyOutput_write_on_startup
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, AM_highFrequencyOutput, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_AM_highFrequencyOutput_enable)
         call mpas_dmpar_bcast_char(dminfo, config_AM_highFrequencyOutput_compute_interval)
         call mpas_dmpar_bcast_char(dminfo, config_AM_highFrequencyOutput_output_stream)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_highFrequencyOutput_compute_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_highFrequencyOutput_write_on_startup)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record AM_highFrequencyOutput')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_AM_highFrequencyOutput_enable = $l', &
logicArgs=(/config_AM_highFrequencyOutput_enable/))
            call mpas_log_write('        config_AM_highFrequencyOutput_compute_interval = '&
//''//mpas_log_escape_dollars(config_AM_highFrequencyOutput_compute_interval))
            call mpas_log_write('        config_AM_highFrequencyOutput_output_stream = '&
//''//mpas_log_escape_dollars(config_AM_highFrequencyOutput_output_stream))
            call mpas_log_write('        config_AM_highFrequencyOutput_compute_on_startup = $l', &
logicArgs=(/config_AM_highFrequencyOutput_compute_on_startup/))
            call mpas_log_write('        config_AM_highFrequencyOutput_write_on_startup = $l', &
logicArgs=(/config_AM_highFrequencyOutput_write_on_startup/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record AM_highFrequencyOutput.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_AM_highFrequencyOutput_enable', config_AM_highFrequencyOutput_enable)
      call mpas_pool_add_config(configPool, 'config_AM_highFrequencyOutput_compute_interval', &
config_AM_highFrequencyOutput_compute_interval)
      call mpas_pool_add_config(configPool, 'config_AM_highFrequencyOutput_output_stream', &
config_AM_highFrequencyOutput_output_stream)
      call mpas_pool_add_config(configPool, 'config_AM_highFrequencyOutput_compute_on_startup', &
config_AM_highFrequencyOutput_compute_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_highFrequencyOutput_write_on_startup', &
config_AM_highFrequencyOutput_write_on_startup)

   end subroutine ocn_setup_nmlrec_AM_highFrequencyOutput


   subroutine ocn_setup_nmlrec_AM_timeFilters(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_AM_timeFilters_enable = .false.
      character (len=StrKIND) :: config_AM_timeFilters_compute_interval = 'dt'
      character (len=StrKIND) :: config_AM_timeFilters_output_stream = 'timeFiltersOutput'
      character (len=StrKIND) :: config_AM_timeFilters_restart_stream = 'timeFiltersRestart'
      logical :: config_AM_timeFilters_compute_on_startup = .true.
      logical :: config_AM_timeFilters_write_on_startup = .true.
      logical :: config_AM_timeFilters_initialize_filters = .true.
      character (len=StrKIND) :: config_AM_timeFilters_tau = '90_00:00:00'
      logical :: config_AM_timeFilters_compute_cell_centered_values = .true.

      namelist /AM_timeFilters/ &
         config_AM_timeFilters_enable, &
         config_AM_timeFilters_compute_interval, &
         config_AM_timeFilters_output_stream, &
         config_AM_timeFilters_restart_stream, &
         config_AM_timeFilters_compute_on_startup, &
         config_AM_timeFilters_write_on_startup, &
         config_AM_timeFilters_initialize_filters, &
         config_AM_timeFilters_tau, &
         config_AM_timeFilters_compute_cell_centered_values
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, AM_timeFilters, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_AM_timeFilters_enable)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeFilters_compute_interval)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeFilters_output_stream)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeFilters_restart_stream)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_timeFilters_compute_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_timeFilters_write_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_timeFilters_initialize_filters)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeFilters_tau)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_timeFilters_compute_cell_centered_values)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record AM_timeFilters')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_AM_timeFilters_enable = $l', logicArgs=(/config_AM_timeFilters_enable/))
            call mpas_log_write('        config_AM_timeFilters_compute_interval = '&
//''//mpas_log_escape_dollars(config_AM_timeFilters_compute_interval))
            call mpas_log_write('        config_AM_timeFilters_output_stream = '&
//''//mpas_log_escape_dollars(config_AM_timeFilters_output_stream))
            call mpas_log_write('        config_AM_timeFilters_restart_stream = '&
//''//mpas_log_escape_dollars(config_AM_timeFilters_restart_stream))
            call mpas_log_write('        config_AM_timeFilters_compute_on_startup = $l', &
logicArgs=(/config_AM_timeFilters_compute_on_startup/))
            call mpas_log_write('        config_AM_timeFilters_write_on_startup = $l', &
logicArgs=(/config_AM_timeFilters_write_on_startup/))
            call mpas_log_write('        config_AM_timeFilters_initialize_filters = $l', &
logicArgs=(/config_AM_timeFilters_initialize_filters/))
            call mpas_log_write('        config_AM_timeFilters_tau = '//mpas_log_escape_dollars(config_AM_timeFilters_tau))
            call mpas_log_write('        config_AM_timeFilters_compute_cell_centered_values = $l', &
logicArgs=(/config_AM_timeFilters_compute_cell_centered_values/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record AM_timeFilters.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_AM_timeFilters_enable', config_AM_timeFilters_enable)
      call mpas_pool_add_config(configPool, 'config_AM_timeFilters_compute_interval', config_AM_timeFilters_compute_interval)
      call mpas_pool_add_config(configPool, 'config_AM_timeFilters_output_stream', config_AM_timeFilters_output_stream)
      call mpas_pool_add_config(configPool, 'config_AM_timeFilters_restart_stream', config_AM_timeFilters_restart_stream)
      call mpas_pool_add_config(configPool, 'config_AM_timeFilters_compute_on_startup', config_AM_timeFilters_compute_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_timeFilters_write_on_startup', config_AM_timeFilters_write_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_timeFilters_initialize_filters', config_AM_timeFilters_initialize_filters)
      call mpas_pool_add_config(configPool, 'config_AM_timeFilters_tau', config_AM_timeFilters_tau)
      call mpas_pool_add_config(configPool, 'config_AM_timeFilters_compute_cell_centered_values', &
config_AM_timeFilters_compute_cell_centered_values)

   end subroutine ocn_setup_nmlrec_AM_timeFilters


   subroutine ocn_setup_nmlrec_AM_lagrPartTrack(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_AM_lagrPartTrack_enable = .false.
      character (len=StrKIND) :: config_AM_lagrPartTrack_compute_interval = 'dt'
      logical :: config_AM_lagrPartTrack_compute_on_startup = .false.
      character (len=StrKIND) :: config_AM_lagrPartTrack_output_stream = 'lagrPartTrackOutput'
      character (len=StrKIND) :: config_AM_lagrPartTrack_restart_stream = 'lagrPartTrackRestart'
      character (len=StrKIND) :: config_AM_lagrPartTrack_input_stream = 'lagrPartTrackInput'
      logical :: config_AM_lagrPartTrack_write_on_startup = .true.
      integer :: config_AM_lagrPartTrack_filter_number = 0
      integer :: config_AM_lagrPartTrack_timeIntegration = 2
      character (len=StrKIND) :: config_AM_lagrPartTrack_reset_criteria = 'none'
      character (len=StrKIND) :: config_AM_lagrPartTrack_reset_global_timestamp = '0000_00:00:00'
      character (len=StrKIND) :: config_AM_lagrPartTrack_region_stream = 'lagrPartTrackRegions'
      logical :: config_AM_lagrPartTrack_reset_if_outside_region = .false.
      logical :: config_AM_lagrPartTrack_reset_if_inside_region = .false.
      logical :: config_AM_lagrPartTrack_sample_horizontal_interp = .true.
      logical :: config_AM_lagrPartTrack_sample_temperature = .true.
      logical :: config_AM_lagrPartTrack_sample_salinity = .true.
      logical :: config_AM_lagrPartTrack_sample_DIC = .false.
      logical :: config_AM_lagrPartTrack_sample_ALK = .false.
      logical :: config_AM_lagrPartTrack_sample_PO4 = .false.
      logical :: config_AM_lagrPartTrack_sample_NO3 = .false.
      logical :: config_AM_lagrPartTrack_sample_SiO3 = .false.
      logical :: config_AM_lagrPartTrack_sample_NH4 = .false.
      logical :: config_AM_lagrPartTrack_sample_Fe = .false.
      logical :: config_AM_lagrPartTrack_sample_O2 = .false.

      namelist /AM_lagrPartTrack/ &
         config_AM_lagrPartTrack_enable, &
         config_AM_lagrPartTrack_compute_interval, &
         config_AM_lagrPartTrack_compute_on_startup, &
         config_AM_lagrPartTrack_output_stream, &
         config_AM_lagrPartTrack_restart_stream, &
         config_AM_lagrPartTrack_input_stream, &
         config_AM_lagrPartTrack_write_on_startup, &
         config_AM_lagrPartTrack_filter_number, &
         config_AM_lagrPartTrack_timeIntegration, &
         config_AM_lagrPartTrack_reset_criteria, &
         config_AM_lagrPartTrack_reset_global_timestamp, &
         config_AM_lagrPartTrack_region_stream, &
         config_AM_lagrPartTrack_reset_if_outside_region, &
         config_AM_lagrPartTrack_reset_if_inside_region, &
         config_AM_lagrPartTrack_sample_horizontal_interp, &
         config_AM_lagrPartTrack_sample_temperature, &
         config_AM_lagrPartTrack_sample_salinity, &
         config_AM_lagrPartTrack_sample_DIC, &
         config_AM_lagrPartTrack_sample_ALK, &
         config_AM_lagrPartTrack_sample_PO4, &
         config_AM_lagrPartTrack_sample_NO3, &
         config_AM_lagrPartTrack_sample_SiO3, &
         config_AM_lagrPartTrack_sample_NH4, &
         config_AM_lagrPartTrack_sample_Fe, &
         config_AM_lagrPartTrack_sample_O2
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, AM_lagrPartTrack, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_AM_lagrPartTrack_enable)
         call mpas_dmpar_bcast_char(dminfo, config_AM_lagrPartTrack_compute_interval)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_lagrPartTrack_compute_on_startup)
         call mpas_dmpar_bcast_char(dminfo, config_AM_lagrPartTrack_output_stream)
         call mpas_dmpar_bcast_char(dminfo, config_AM_lagrPartTrack_restart_stream)
         call mpas_dmpar_bcast_char(dminfo, config_AM_lagrPartTrack_input_stream)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_lagrPartTrack_write_on_startup)
         call mpas_dmpar_bcast_int(dminfo, config_AM_lagrPartTrack_filter_number)
         call mpas_dmpar_bcast_int(dminfo, config_AM_lagrPartTrack_timeIntegration)
         call mpas_dmpar_bcast_char(dminfo, config_AM_lagrPartTrack_reset_criteria)
         call mpas_dmpar_bcast_char(dminfo, config_AM_lagrPartTrack_reset_global_timestamp)
         call mpas_dmpar_bcast_char(dminfo, config_AM_lagrPartTrack_region_stream)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_lagrPartTrack_reset_if_outside_region)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_lagrPartTrack_reset_if_inside_region)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_lagrPartTrack_sample_horizontal_interp)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_lagrPartTrack_sample_temperature)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_lagrPartTrack_sample_salinity)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_lagrPartTrack_sample_DIC)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_lagrPartTrack_sample_ALK)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_lagrPartTrack_sample_PO4)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_lagrPartTrack_sample_NO3)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_lagrPartTrack_sample_SiO3)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_lagrPartTrack_sample_NH4)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_lagrPartTrack_sample_Fe)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_lagrPartTrack_sample_O2)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record AM_lagrPartTrack')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_AM_lagrPartTrack_enable = $l', logicArgs=(/config_AM_lagrPartTrack_enable/))
            call mpas_log_write('        config_AM_lagrPartTrack_compute_interval = '&
//''//mpas_log_escape_dollars(config_AM_lagrPartTrack_compute_interval))
            call mpas_log_write('        config_AM_lagrPartTrack_compute_on_startup = $l', &
logicArgs=(/config_AM_lagrPartTrack_compute_on_startup/))
            call mpas_log_write('        config_AM_lagrPartTrack_output_stream = '&
//''//mpas_log_escape_dollars(config_AM_lagrPartTrack_output_stream))
            call mpas_log_write('        config_AM_lagrPartTrack_restart_stream = '&
//''//mpas_log_escape_dollars(config_AM_lagrPartTrack_restart_stream))
            call mpas_log_write('        config_AM_lagrPartTrack_input_stream = '&
//''//mpas_log_escape_dollars(config_AM_lagrPartTrack_input_stream))
            call mpas_log_write('        config_AM_lagrPartTrack_write_on_startup = $l', &
logicArgs=(/config_AM_lagrPartTrack_write_on_startup/))
            call mpas_log_write('        config_AM_lagrPartTrack_filter_number = $i', &
intArgs=(/config_AM_lagrPartTrack_filter_number/))
            call mpas_log_write('        config_AM_lagrPartTrack_timeIntegration = $i', &
intArgs=(/config_AM_lagrPartTrack_timeIntegration/))
            call mpas_log_write('        config_AM_lagrPartTrack_reset_criteria = '&
//''//mpas_log_escape_dollars(config_AM_lagrPartTrack_reset_criteria))
            call mpas_log_write('        config_AM_lagrPartTrack_reset_global_timestamp = '&
//''//mpas_log_escape_dollars(config_AM_lagrPartTrack_reset_global_timestamp))
            call mpas_log_write('        config_AM_lagrPartTrack_region_stream = '&
//''//mpas_log_escape_dollars(config_AM_lagrPartTrack_region_stream))
            call mpas_log_write('        config_AM_lagrPartTrack_reset_if_outside_region = $l', &
logicArgs=(/config_AM_lagrPartTrack_reset_if_outside_region/))
            call mpas_log_write('        config_AM_lagrPartTrack_reset_if_inside_region = $l', &
logicArgs=(/config_AM_lagrPartTrack_reset_if_inside_region/))
            call mpas_log_write('        config_AM_lagrPartTrack_sample_horizontal_interp = $l', &
logicArgs=(/config_AM_lagrPartTrack_sample_horizontal_interp/))
            call mpas_log_write('        config_AM_lagrPartTrack_sample_temperature = $l', &
logicArgs=(/config_AM_lagrPartTrack_sample_temperature/))
            call mpas_log_write('        config_AM_lagrPartTrack_sample_salinity = $l', &
logicArgs=(/config_AM_lagrPartTrack_sample_salinity/))
            call mpas_log_write('        config_AM_lagrPartTrack_sample_DIC = $l', &
logicArgs=(/config_AM_lagrPartTrack_sample_DIC/))
            call mpas_log_write('        config_AM_lagrPartTrack_sample_ALK = $l', &
logicArgs=(/config_AM_lagrPartTrack_sample_ALK/))
            call mpas_log_write('        config_AM_lagrPartTrack_sample_PO4 = $l', &
logicArgs=(/config_AM_lagrPartTrack_sample_PO4/))
            call mpas_log_write('        config_AM_lagrPartTrack_sample_NO3 = $l', &
logicArgs=(/config_AM_lagrPartTrack_sample_NO3/))
            call mpas_log_write('        config_AM_lagrPartTrack_sample_SiO3 = $l', &
logicArgs=(/config_AM_lagrPartTrack_sample_SiO3/))
            call mpas_log_write('        config_AM_lagrPartTrack_sample_NH4 = $l', &
logicArgs=(/config_AM_lagrPartTrack_sample_NH4/))
            call mpas_log_write('        config_AM_lagrPartTrack_sample_Fe = $l', logicArgs=(/config_AM_lagrPartTrack_sample_Fe/))
            call mpas_log_write('        config_AM_lagrPartTrack_sample_O2 = $l', logicArgs=(/config_AM_lagrPartTrack_sample_O2/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record AM_lagrPartTrack.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_AM_lagrPartTrack_enable', config_AM_lagrPartTrack_enable)
      call mpas_pool_add_config(configPool, 'config_AM_lagrPartTrack_compute_interval', config_AM_lagrPartTrack_compute_interval)
      call mpas_pool_add_config(configPool, 'config_AM_lagrPartTrack_compute_on_startup', &
config_AM_lagrPartTrack_compute_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_lagrPartTrack_output_stream', config_AM_lagrPartTrack_output_stream)
      call mpas_pool_add_config(configPool, 'config_AM_lagrPartTrack_restart_stream', config_AM_lagrPartTrack_restart_stream)
      call mpas_pool_add_config(configPool, 'config_AM_lagrPartTrack_input_stream', config_AM_lagrPartTrack_input_stream)
      call mpas_pool_add_config(configPool, 'config_AM_lagrPartTrack_write_on_startup', config_AM_lagrPartTrack_write_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_lagrPartTrack_filter_number', config_AM_lagrPartTrack_filter_number)
      call mpas_pool_add_config(configPool, 'config_AM_lagrPartTrack_timeIntegration', config_AM_lagrPartTrack_timeIntegration)
      call mpas_pool_add_config(configPool, 'config_AM_lagrPartTrack_reset_criteria', config_AM_lagrPartTrack_reset_criteria)
      call mpas_pool_add_config(configPool, 'config_AM_lagrPartTrack_reset_global_timestamp', &
config_AM_lagrPartTrack_reset_global_timestamp)
      call mpas_pool_add_config(configPool, 'config_AM_lagrPartTrack_region_stream', config_AM_lagrPartTrack_region_stream)
      call mpas_pool_add_config(configPool, 'config_AM_lagrPartTrack_reset_if_outside_region', &
config_AM_lagrPartTrack_reset_if_outside_region)
      call mpas_pool_add_config(configPool, 'config_AM_lagrPartTrack_reset_if_inside_region', &
config_AM_lagrPartTrack_reset_if_inside_region)
      call mpas_pool_add_config(configPool, 'config_AM_lagrPartTrack_sample_horizontal_interp', &
config_AM_lagrPartTrack_sample_horizontal_interp)
      call mpas_pool_add_config(configPool, 'config_AM_lagrPartTrack_sample_temperature', &
config_AM_lagrPartTrack_sample_temperature)
      call mpas_pool_add_config(configPool, 'config_AM_lagrPartTrack_sample_salinity', config_AM_lagrPartTrack_sample_salinity)
      call mpas_pool_add_config(configPool, 'config_AM_lagrPartTrack_sample_DIC', config_AM_lagrPartTrack_sample_DIC)
      call mpas_pool_add_config(configPool, 'config_AM_lagrPartTrack_sample_ALK', config_AM_lagrPartTrack_sample_ALK)
      call mpas_pool_add_config(configPool, 'config_AM_lagrPartTrack_sample_PO4', config_AM_lagrPartTrack_sample_PO4)
      call mpas_pool_add_config(configPool, 'config_AM_lagrPartTrack_sample_NO3', config_AM_lagrPartTrack_sample_NO3)
      call mpas_pool_add_config(configPool, 'config_AM_lagrPartTrack_sample_SiO3', config_AM_lagrPartTrack_sample_SiO3)
      call mpas_pool_add_config(configPool, 'config_AM_lagrPartTrack_sample_NH4', config_AM_lagrPartTrack_sample_NH4)
      call mpas_pool_add_config(configPool, 'config_AM_lagrPartTrack_sample_Fe', config_AM_lagrPartTrack_sample_Fe)
      call mpas_pool_add_config(configPool, 'config_AM_lagrPartTrack_sample_O2', config_AM_lagrPartTrack_sample_O2)

   end subroutine ocn_setup_nmlrec_AM_lagrPartTrack


   subroutine ocn_setup_nmlrec_AM_eliassenPalm(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_AM_eliassenPalm_enable = .false.
      character (len=StrKIND) :: config_AM_eliassenPalm_compute_interval = 'output_interval'
      character (len=StrKIND) :: config_AM_eliassenPalm_output_stream = 'eliassenPalmOutput'
      character (len=StrKIND) :: config_AM_eliassenPalm_restart_stream = 'eliassenPalmRestart'
      logical :: config_AM_eliassenPalm_compute_on_startup = .true.
      logical :: config_AM_eliassenPalm_write_on_startup = .true.
      logical :: config_AM_eliassenPalm_debug = .false.
      integer :: config_AM_eliassenPalm_nBuoyancyLayers = 45
      real (kind=RKIND) :: config_AM_eliassenPalm_rhomin_buoycoor = 900.000000
      real (kind=RKIND) :: config_AM_eliassenPalm_rhomax_buoycoor = 1080.000000

      namelist /AM_eliassenPalm/ &
         config_AM_eliassenPalm_enable, &
         config_AM_eliassenPalm_compute_interval, &
         config_AM_eliassenPalm_output_stream, &
         config_AM_eliassenPalm_restart_stream, &
         config_AM_eliassenPalm_compute_on_startup, &
         config_AM_eliassenPalm_write_on_startup, &
         config_AM_eliassenPalm_debug, &
         config_AM_eliassenPalm_nBuoyancyLayers, &
         config_AM_eliassenPalm_rhomin_buoycoor, &
         config_AM_eliassenPalm_rhomax_buoycoor
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, AM_eliassenPalm, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_AM_eliassenPalm_enable)
         call mpas_dmpar_bcast_char(dminfo, config_AM_eliassenPalm_compute_interval)
         call mpas_dmpar_bcast_char(dminfo, config_AM_eliassenPalm_output_stream)
         call mpas_dmpar_bcast_char(dminfo, config_AM_eliassenPalm_restart_stream)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_eliassenPalm_compute_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_eliassenPalm_write_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_eliassenPalm_debug)
         call mpas_dmpar_bcast_int(dminfo, config_AM_eliassenPalm_nBuoyancyLayers)
         call mpas_dmpar_bcast_real(dminfo, config_AM_eliassenPalm_rhomin_buoycoor)
         call mpas_dmpar_bcast_real(dminfo, config_AM_eliassenPalm_rhomax_buoycoor)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record AM_eliassenPalm')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_AM_eliassenPalm_enable = $l', logicArgs=(/config_AM_eliassenPalm_enable/))
            call mpas_log_write('        config_AM_eliassenPalm_compute_interval = '&
//''//mpas_log_escape_dollars(config_AM_eliassenPalm_compute_interval))
            call mpas_log_write('        config_AM_eliassenPalm_output_stream = '&
//''//mpas_log_escape_dollars(config_AM_eliassenPalm_output_stream))
            call mpas_log_write('        config_AM_eliassenPalm_restart_stream = '&
//''//mpas_log_escape_dollars(config_AM_eliassenPalm_restart_stream))
            call mpas_log_write('        config_AM_eliassenPalm_compute_on_startup = $l', &
logicArgs=(/config_AM_eliassenPalm_compute_on_startup/))
            call mpas_log_write('        config_AM_eliassenPalm_write_on_startup = $l', &
logicArgs=(/config_AM_eliassenPalm_write_on_startup/))
            call mpas_log_write('        config_AM_eliassenPalm_debug = $l', logicArgs=(/config_AM_eliassenPalm_debug/))
            call mpas_log_write('        config_AM_eliassenPalm_nBuoyancyLayers = $i', &
intArgs=(/config_AM_eliassenPalm_nBuoyancyLayers/))
            call mpas_log_write('        config_AM_eliassenPalm_rhomin_buoycoor = $r', &
realArgs=(/config_AM_eliassenPalm_rhomin_buoycoor/))
            call mpas_log_write('        config_AM_eliassenPalm_rhomax_buoycoor = $r', &
realArgs=(/config_AM_eliassenPalm_rhomax_buoycoor/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record AM_eliassenPalm.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_AM_eliassenPalm_enable', config_AM_eliassenPalm_enable)
      call mpas_pool_add_config(configPool, 'config_AM_eliassenPalm_compute_interval', config_AM_eliassenPalm_compute_interval)
      call mpas_pool_add_config(configPool, 'config_AM_eliassenPalm_output_stream', config_AM_eliassenPalm_output_stream)
      call mpas_pool_add_config(configPool, 'config_AM_eliassenPalm_restart_stream', config_AM_eliassenPalm_restart_stream)
      call mpas_pool_add_config(configPool, 'config_AM_eliassenPalm_compute_on_startup', config_AM_eliassenPalm_compute_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_eliassenPalm_write_on_startup', config_AM_eliassenPalm_write_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_eliassenPalm_debug', config_AM_eliassenPalm_debug)
      call mpas_pool_add_config(configPool, 'config_AM_eliassenPalm_nBuoyancyLayers', config_AM_eliassenPalm_nBuoyancyLayers)
      call mpas_pool_add_config(configPool, 'config_AM_eliassenPalm_rhomin_buoycoor', config_AM_eliassenPalm_rhomin_buoycoor)
      call mpas_pool_add_config(configPool, 'config_AM_eliassenPalm_rhomax_buoycoor', config_AM_eliassenPalm_rhomax_buoycoor)

   end subroutine ocn_setup_nmlrec_AM_eliassenPalm


   subroutine ocn_setup_nmlrec_AM_mixedLayerDepths(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_AM_mixedLayerDepths_enable = .false.
      character (len=StrKIND) :: config_AM_mixedLayerDepths_compute_interval = 'output_interval'
      character (len=StrKIND) :: config_AM_mixedLayerDepths_output_stream = 'mixedLayerDepthsOutput'
      logical :: config_AM_mixedLayerDepths_write_on_startup = .true.
      logical :: config_AM_mixedLayerDepths_compute_on_startup = .true.
      logical :: config_AM_mixedLayerDepths_Tthreshold = .true.
      real (kind=RKIND) :: config_AM_mixedLayerDepths_crit_temp_threshold = 0.200000
      real (kind=RKIND) :: config_AM_mixedLayerDepths_reference_pressure = 100000.000000
      logical :: config_AM_mixedLayerDepths_Tgradient = .false.
      logical :: config_AM_mixedLayerDepths_Dgradient = .false.
      real (kind=RKIND) :: config_AM_mixedLayerDepths_temp_gradient_threshold = 0.000000
      real (kind=RKIND) :: config_AM_mixedLayerDepths_den_gradient_threshold = 0.000000
      integer :: config_AM_mixedLayerDepths_interp_method = 1

      namelist /AM_mixedLayerDepths/ &
         config_AM_mixedLayerDepths_enable, &
         config_AM_mixedLayerDepths_compute_interval, &
         config_AM_mixedLayerDepths_output_stream, &
         config_AM_mixedLayerDepths_write_on_startup, &
         config_AM_mixedLayerDepths_compute_on_startup, &
         config_AM_mixedLayerDepths_Tthreshold, &
         config_AM_mixedLayerDepths_crit_temp_threshold, &
         config_AM_mixedLayerDepths_reference_pressure, &
         config_AM_mixedLayerDepths_Tgradient, &
         config_AM_mixedLayerDepths_Dgradient, &
         config_AM_mixedLayerDepths_temp_gradient_threshold, &
         config_AM_mixedLayerDepths_den_gradient_threshold, &
         config_AM_mixedLayerDepths_interp_method
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, AM_mixedLayerDepths, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_AM_mixedLayerDepths_enable)
         call mpas_dmpar_bcast_char(dminfo, config_AM_mixedLayerDepths_compute_interval)
         call mpas_dmpar_bcast_char(dminfo, config_AM_mixedLayerDepths_output_stream)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_mixedLayerDepths_write_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_mixedLayerDepths_compute_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_mixedLayerDepths_Tthreshold)
         call mpas_dmpar_bcast_real(dminfo, config_AM_mixedLayerDepths_crit_temp_threshold)
         call mpas_dmpar_bcast_real(dminfo, config_AM_mixedLayerDepths_reference_pressure)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_mixedLayerDepths_Tgradient)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_mixedLayerDepths_Dgradient)
         call mpas_dmpar_bcast_real(dminfo, config_AM_mixedLayerDepths_temp_gradient_threshold)
         call mpas_dmpar_bcast_real(dminfo, config_AM_mixedLayerDepths_den_gradient_threshold)
         call mpas_dmpar_bcast_int(dminfo, config_AM_mixedLayerDepths_interp_method)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record AM_mixedLayerDepths')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_AM_mixedLayerDepths_enable = $l', logicArgs=(/config_AM_mixedLayerDepths_enable/))
            call mpas_log_write('        config_AM_mixedLayerDepths_compute_interval = '&
//''//mpas_log_escape_dollars(config_AM_mixedLayerDepths_compute_interval))
            call mpas_log_write('        config_AM_mixedLayerDepths_output_stream = '&
//''//mpas_log_escape_dollars(config_AM_mixedLayerDepths_output_stream))
            call mpas_log_write('        config_AM_mixedLayerDepths_write_on_startup = $l', &
logicArgs=(/config_AM_mixedLayerDepths_write_on_startup/))
            call mpas_log_write('        config_AM_mixedLayerDepths_compute_on_startup = $l', &
logicArgs=(/config_AM_mixedLayerDepths_compute_on_startup/))
            call mpas_log_write('        config_AM_mixedLayerDepths_Tthreshold = $l', &
logicArgs=(/config_AM_mixedLayerDepths_Tthreshold/))
            call mpas_log_write('        config_AM_mixedLayerDepths_crit_temp_threshold = $r', &
realArgs=(/config_AM_mixedLayerDepths_crit_temp_threshold/))
            call mpas_log_write('        config_AM_mixedLayerDepths_reference_pressure = $r', &
realArgs=(/config_AM_mixedLayerDepths_reference_pressure/))
            call mpas_log_write('        config_AM_mixedLayerDepths_Tgradient = $l', &
logicArgs=(/config_AM_mixedLayerDepths_Tgradient/))
            call mpas_log_write('        config_AM_mixedLayerDepths_Dgradient = $l', &
logicArgs=(/config_AM_mixedLayerDepths_Dgradient/))
            call mpas_log_write('        config_AM_mixedLayerDepths_temp_gradient_threshold = $r', &
realArgs=(/config_AM_mixedLayerDepths_temp_gradient_threshold/))
            call mpas_log_write('        config_AM_mixedLayerDepths_den_gradient_threshold = $r', &
realArgs=(/config_AM_mixedLayerDepths_den_gradient_threshold/))
            call mpas_log_write('        config_AM_mixedLayerDepths_interp_method = $i', &
intArgs=(/config_AM_mixedLayerDepths_interp_method/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record AM_mixedLayerDepths.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_AM_mixedLayerDepths_enable', config_AM_mixedLayerDepths_enable)
      call mpas_pool_add_config(configPool, 'config_AM_mixedLayerDepths_compute_interval', &
config_AM_mixedLayerDepths_compute_interval)
      call mpas_pool_add_config(configPool, 'config_AM_mixedLayerDepths_output_stream', config_AM_mixedLayerDepths_output_stream)
      call mpas_pool_add_config(configPool, 'config_AM_mixedLayerDepths_write_on_startup', &
config_AM_mixedLayerDepths_write_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_mixedLayerDepths_compute_on_startup', &
config_AM_mixedLayerDepths_compute_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_mixedLayerDepths_Tthreshold', config_AM_mixedLayerDepths_Tthreshold)
      call mpas_pool_add_config(configPool, 'config_AM_mixedLayerDepths_crit_temp_threshold', &
config_AM_mixedLayerDepths_crit_temp_threshold)
      call mpas_pool_add_config(configPool, 'config_AM_mixedLayerDepths_reference_pressure', &
config_AM_mixedLayerDepths_reference_pressure)
      call mpas_pool_add_config(configPool, 'config_AM_mixedLayerDepths_Tgradient', config_AM_mixedLayerDepths_Tgradient)
      call mpas_pool_add_config(configPool, 'config_AM_mixedLayerDepths_Dgradient', config_AM_mixedLayerDepths_Dgradient)
      call mpas_pool_add_config(configPool, 'config_AM_mixedLayerDepths_temp_gradient_threshold', &
config_AM_mixedLayerDepths_temp_gradient_threshold)
      call mpas_pool_add_config(configPool, 'config_AM_mixedLayerDepths_den_gradient_threshold', &
config_AM_mixedLayerDepths_den_gradient_threshold)
      call mpas_pool_add_config(configPool, 'config_AM_mixedLayerDepths_interp_method', config_AM_mixedLayerDepths_interp_method)

   end subroutine ocn_setup_nmlrec_AM_mixedLayerDepths


   subroutine ocn_setup_nmlrec_AM_regionalStatsDaily(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_AM_regionalStatsDaily_enable = .false.
      logical :: config_AM_regionalStatsDaily_compute_on_startup = .false.
      logical :: config_AM_regionalStatsDaily_write_on_startup = .false.
      character (len=StrKIND) :: config_AM_regionalStatsDaily_compute_interval = 'output_interval'
      character (len=StrKIND) :: config_AM_regionalStatsDaily_output_stream = 'regionalStatsDailyOutput'
      character (len=StrKIND) :: config_AM_regionalStatsDaily_restart_stream = 'regionalMasksInput'
      character (len=StrKIND) :: config_AM_regionalStatsDaily_input_stream = 'regionalMasksInput'
      character (len=StrKIND) :: config_AM_regionalStatsDaily_operation = 'avg'
      character (len=StrKIND) :: config_AM_regionalStatsDaily_region_type = 'cell'
      character (len=StrKIND) :: config_AM_regionalStatsDaily_region_group = 'all'
      character (len=StrKIND) :: config_AM_regionalStatsDaily_1d_weighting_function = 'mul'
      character (len=StrKIND) :: config_AM_regionalStatsDaily_2d_weighting_function = 'mul'
      character (len=StrKIND) :: config_AM_regionalStatsDaily_1d_weighting_field = 'areaCell'
      character (len=StrKIND) :: config_AM_regionalStatsDaily_2d_weighting_field = 'volumeCell'
      character (len=StrKIND) :: config_AM_regionalStatsDaily_vertical_mask = 'cellMask'
      character (len=StrKIND) :: config_AM_regionalStatsDaily_vertical_dimension = 'nVertLevels'

      namelist /AM_regionalStatsDaily/ &
         config_AM_regionalStatsDaily_enable, &
         config_AM_regionalStatsDaily_compute_on_startup, &
         config_AM_regionalStatsDaily_write_on_startup, &
         config_AM_regionalStatsDaily_compute_interval, &
         config_AM_regionalStatsDaily_output_stream, &
         config_AM_regionalStatsDaily_restart_stream, &
         config_AM_regionalStatsDaily_input_stream, &
         config_AM_regionalStatsDaily_operation, &
         config_AM_regionalStatsDaily_region_type, &
         config_AM_regionalStatsDaily_region_group, &
         config_AM_regionalStatsDaily_1d_weighting_function, &
         config_AM_regionalStatsDaily_2d_weighting_function, &
         config_AM_regionalStatsDaily_1d_weighting_field, &
         config_AM_regionalStatsDaily_2d_weighting_field, &
         config_AM_regionalStatsDaily_vertical_mask, &
         config_AM_regionalStatsDaily_vertical_dimension
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, AM_regionalStatsDaily, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_AM_regionalStatsDaily_enable)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_regionalStatsDaily_compute_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_regionalStatsDaily_write_on_startup)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsDaily_compute_interval)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsDaily_output_stream)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsDaily_restart_stream)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsDaily_input_stream)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsDaily_operation)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsDaily_region_type)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsDaily_region_group)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsDaily_1d_weighting_function)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsDaily_2d_weighting_function)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsDaily_1d_weighting_field)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsDaily_2d_weighting_field)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsDaily_vertical_mask)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsDaily_vertical_dimension)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record AM_regionalStatsDaily')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_AM_regionalStatsDaily_enable = $l', &
logicArgs=(/config_AM_regionalStatsDaily_enable/))
            call mpas_log_write('        config_AM_regionalStatsDaily_compute_on_startup = $l', &
logicArgs=(/config_AM_regionalStatsDaily_compute_on_startup/))
            call mpas_log_write('        config_AM_regionalStatsDaily_write_on_startup = $l', &
logicArgs=(/config_AM_regionalStatsDaily_write_on_startup/))
            call mpas_log_write('        config_AM_regionalStatsDaily_compute_interval = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsDaily_compute_interval))
            call mpas_log_write('        config_AM_regionalStatsDaily_output_stream = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsDaily_output_stream))
            call mpas_log_write('        config_AM_regionalStatsDaily_restart_stream = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsDaily_restart_stream))
            call mpas_log_write('        config_AM_regionalStatsDaily_input_stream = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsDaily_input_stream))
            call mpas_log_write('        config_AM_regionalStatsDaily_operation = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsDaily_operation))
            call mpas_log_write('        config_AM_regionalStatsDaily_region_type = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsDaily_region_type))
            call mpas_log_write('        config_AM_regionalStatsDaily_region_group = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsDaily_region_group))
            call mpas_log_write('        config_AM_regionalStatsDaily_1d_weighting_function = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsDaily_1d_weighting_function))
            call mpas_log_write('        config_AM_regionalStatsDaily_2d_weighting_function = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsDaily_2d_weighting_function))
            call mpas_log_write('        config_AM_regionalStatsDaily_1d_weighting_field = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsDaily_1d_weighting_field))
            call mpas_log_write('        config_AM_regionalStatsDaily_2d_weighting_field = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsDaily_2d_weighting_field))
            call mpas_log_write('        config_AM_regionalStatsDaily_vertical_mask = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsDaily_vertical_mask))
            call mpas_log_write('        config_AM_regionalStatsDaily_vertical_dimension = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsDaily_vertical_dimension))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record AM_regionalStatsDaily.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsDaily_enable', config_AM_regionalStatsDaily_enable)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsDaily_compute_on_startup', &
config_AM_regionalStatsDaily_compute_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsDaily_write_on_startup', &
config_AM_regionalStatsDaily_write_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsDaily_compute_interval', &
config_AM_regionalStatsDaily_compute_interval)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsDaily_output_stream', &
config_AM_regionalStatsDaily_output_stream)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsDaily_restart_stream', &
config_AM_regionalStatsDaily_restart_stream)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsDaily_input_stream', config_AM_regionalStatsDaily_input_stream)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsDaily_operation', config_AM_regionalStatsDaily_operation)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsDaily_region_type', config_AM_regionalStatsDaily_region_type)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsDaily_region_group', config_AM_regionalStatsDaily_region_group)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsDaily_1d_weighting_function', &
config_AM_regionalStatsDaily_1d_weighting_function)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsDaily_2d_weighting_function', &
config_AM_regionalStatsDaily_2d_weighting_function)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsDaily_1d_weighting_field', &
config_AM_regionalStatsDaily_1d_weighting_field)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsDaily_2d_weighting_field', &
config_AM_regionalStatsDaily_2d_weighting_field)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsDaily_vertical_mask', &
config_AM_regionalStatsDaily_vertical_mask)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsDaily_vertical_dimension', &
config_AM_regionalStatsDaily_vertical_dimension)

   end subroutine ocn_setup_nmlrec_AM_regionalStatsDaily


   subroutine ocn_setup_nmlrec_AM_regionalStatsWeekly(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_AM_regionalStatsWeekly_enable = .false.
      logical :: config_AM_regionalStatsWeekly_compute_on_startup = .false.
      logical :: config_AM_regionalStatsWeekly_write_on_startup = .false.
      character (len=StrKIND) :: config_AM_regionalStatsWeekly_compute_interval = 'output_interval'
      character (len=StrKIND) :: config_AM_regionalStatsWeekly_output_stream = 'regionalStatsWeeklyOutput'
      character (len=StrKIND) :: config_AM_regionalStatsWeekly_restart_stream = 'regionalMasksInput'
      character (len=StrKIND) :: config_AM_regionalStatsWeekly_input_stream = 'regionalMasksInput'
      character (len=StrKIND) :: config_AM_regionalStatsWeekly_operation = 'avg'
      character (len=StrKIND) :: config_AM_regionalStatsWeekly_region_type = 'cell'
      character (len=StrKIND) :: config_AM_regionalStatsWeekly_region_group = 'all'
      character (len=StrKIND) :: config_AM_regionalStatsWeekly_1d_weighting_function = 'mul'
      character (len=StrKIND) :: config_AM_regionalStatsWeekly_2d_weighting_function = 'mul'
      character (len=StrKIND) :: config_AM_regionalStatsWeekly_1d_weighting_field = 'areaCell'
      character (len=StrKIND) :: config_AM_regionalStatsWeekly_2d_weighting_field = 'volumeCell'
      character (len=StrKIND) :: config_AM_regionalStatsWeekly_vertical_mask = 'cellMask'
      character (len=StrKIND) :: config_AM_regionalStatsWeekly_vertical_dimension = 'nVertLevels'

      namelist /AM_regionalStatsWeekly/ &
         config_AM_regionalStatsWeekly_enable, &
         config_AM_regionalStatsWeekly_compute_on_startup, &
         config_AM_regionalStatsWeekly_write_on_startup, &
         config_AM_regionalStatsWeekly_compute_interval, &
         config_AM_regionalStatsWeekly_output_stream, &
         config_AM_regionalStatsWeekly_restart_stream, &
         config_AM_regionalStatsWeekly_input_stream, &
         config_AM_regionalStatsWeekly_operation, &
         config_AM_regionalStatsWeekly_region_type, &
         config_AM_regionalStatsWeekly_region_group, &
         config_AM_regionalStatsWeekly_1d_weighting_function, &
         config_AM_regionalStatsWeekly_2d_weighting_function, &
         config_AM_regionalStatsWeekly_1d_weighting_field, &
         config_AM_regionalStatsWeekly_2d_weighting_field, &
         config_AM_regionalStatsWeekly_vertical_mask, &
         config_AM_regionalStatsWeekly_vertical_dimension
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, AM_regionalStatsWeekly, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_AM_regionalStatsWeekly_enable)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_regionalStatsWeekly_compute_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_regionalStatsWeekly_write_on_startup)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsWeekly_compute_interval)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsWeekly_output_stream)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsWeekly_restart_stream)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsWeekly_input_stream)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsWeekly_operation)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsWeekly_region_type)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsWeekly_region_group)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsWeekly_1d_weighting_function)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsWeekly_2d_weighting_function)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsWeekly_1d_weighting_field)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsWeekly_2d_weighting_field)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsWeekly_vertical_mask)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsWeekly_vertical_dimension)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record AM_regionalStatsWeekly')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_AM_regionalStatsWeekly_enable = $l', &
logicArgs=(/config_AM_regionalStatsWeekly_enable/))
            call mpas_log_write('        config_AM_regionalStatsWeekly_compute_on_startup = $l', &
logicArgs=(/config_AM_regionalStatsWeekly_compute_on_startup/))
            call mpas_log_write('        config_AM_regionalStatsWeekly_write_on_startup = $l', &
logicArgs=(/config_AM_regionalStatsWeekly_write_on_startup/))
            call mpas_log_write('        config_AM_regionalStatsWeekly_compute_interval = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsWeekly_compute_interval))
            call mpas_log_write('        config_AM_regionalStatsWeekly_output_stream = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsWeekly_output_stream))
            call mpas_log_write('        config_AM_regionalStatsWeekly_restart_stream = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsWeekly_restart_stream))
            call mpas_log_write('        config_AM_regionalStatsWeekly_input_stream = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsWeekly_input_stream))
            call mpas_log_write('        config_AM_regionalStatsWeekly_operation = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsWeekly_operation))
            call mpas_log_write('        config_AM_regionalStatsWeekly_region_type = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsWeekly_region_type))
            call mpas_log_write('        config_AM_regionalStatsWeekly_region_group = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsWeekly_region_group))
            call mpas_log_write('        config_AM_regionalStatsWeekly_1d_weighting_function = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsWeekly_1d_weighting_function))
            call mpas_log_write('        config_AM_regionalStatsWeekly_2d_weighting_function = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsWeekly_2d_weighting_function))
            call mpas_log_write('        config_AM_regionalStatsWeekly_1d_weighting_field = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsWeekly_1d_weighting_field))
            call mpas_log_write('        config_AM_regionalStatsWeekly_2d_weighting_field = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsWeekly_2d_weighting_field))
            call mpas_log_write('        config_AM_regionalStatsWeekly_vertical_mask = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsWeekly_vertical_mask))
            call mpas_log_write('        config_AM_regionalStatsWeekly_vertical_dimension = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsWeekly_vertical_dimension))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record AM_regionalStatsWeekly.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsWeekly_enable', config_AM_regionalStatsWeekly_enable)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsWeekly_compute_on_startup', &
config_AM_regionalStatsWeekly_compute_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsWeekly_write_on_startup', &
config_AM_regionalStatsWeekly_write_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsWeekly_compute_interval', &
config_AM_regionalStatsWeekly_compute_interval)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsWeekly_output_stream', &
config_AM_regionalStatsWeekly_output_stream)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsWeekly_restart_stream', &
config_AM_regionalStatsWeekly_restart_stream)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsWeekly_input_stream', &
config_AM_regionalStatsWeekly_input_stream)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsWeekly_operation', config_AM_regionalStatsWeekly_operation)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsWeekly_region_type', config_AM_regionalStatsWeekly_region_type)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsWeekly_region_group', &
config_AM_regionalStatsWeekly_region_group)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsWeekly_1d_weighting_function', &
config_AM_regionalStatsWeekly_1d_weighting_function)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsWeekly_2d_weighting_function', &
config_AM_regionalStatsWeekly_2d_weighting_function)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsWeekly_1d_weighting_field', &
config_AM_regionalStatsWeekly_1d_weighting_field)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsWeekly_2d_weighting_field', &
config_AM_regionalStatsWeekly_2d_weighting_field)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsWeekly_vertical_mask', &
config_AM_regionalStatsWeekly_vertical_mask)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsWeekly_vertical_dimension', &
config_AM_regionalStatsWeekly_vertical_dimension)

   end subroutine ocn_setup_nmlrec_AM_regionalStatsWeekly


   subroutine ocn_setup_nmlrec_AM_regionalStatsMonthly(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_AM_regionalStatsMonthly_enable = .false.
      logical :: config_AM_regionalStatsMonthly_compute_on_startup = .false.
      logical :: config_AM_regionalStatsMonthly_write_on_startup = .false.
      character (len=StrKIND) :: config_AM_regionalStatsMonthly_compute_interval = 'output_interval'
      character (len=StrKIND) :: config_AM_regionalStatsMonthly_output_stream = 'regionalStatsMonthlyOutput'
      character (len=StrKIND) :: config_AM_regionalStatsMonthly_restart_stream = 'regionalMasksInput'
      character (len=StrKIND) :: config_AM_regionalStatsMonthly_input_stream = 'regionalMasksInput'
      character (len=StrKIND) :: config_AM_regionalStatsMonthly_operation = 'avg'
      character (len=StrKIND) :: config_AM_regionalStatsMonthly_region_type = 'cell'
      character (len=StrKIND) :: config_AM_regionalStatsMonthly_region_group = 'all'
      character (len=StrKIND) :: config_AM_regionalStatsMonthly_1d_weighting_function = 'mul'
      character (len=StrKIND) :: config_AM_regionalStatsMonthly_2d_weighting_function = 'mul'
      character (len=StrKIND) :: config_AM_regionalStatsMonthly_1d_weighting_field = 'areaCell'
      character (len=StrKIND) :: config_AM_regionalStatsMonthly_2d_weighting_field = 'volumeCell'
      character (len=StrKIND) :: config_AM_regionalStatsMonthly_vertical_mask = 'cellMask'
      character (len=StrKIND) :: config_AM_regionalStatsMonthly_vertical_dimension = 'nVertLevels'

      namelist /AM_regionalStatsMonthly/ &
         config_AM_regionalStatsMonthly_enable, &
         config_AM_regionalStatsMonthly_compute_on_startup, &
         config_AM_regionalStatsMonthly_write_on_startup, &
         config_AM_regionalStatsMonthly_compute_interval, &
         config_AM_regionalStatsMonthly_output_stream, &
         config_AM_regionalStatsMonthly_restart_stream, &
         config_AM_regionalStatsMonthly_input_stream, &
         config_AM_regionalStatsMonthly_operation, &
         config_AM_regionalStatsMonthly_region_type, &
         config_AM_regionalStatsMonthly_region_group, &
         config_AM_regionalStatsMonthly_1d_weighting_function, &
         config_AM_regionalStatsMonthly_2d_weighting_function, &
         config_AM_regionalStatsMonthly_1d_weighting_field, &
         config_AM_regionalStatsMonthly_2d_weighting_field, &
         config_AM_regionalStatsMonthly_vertical_mask, &
         config_AM_regionalStatsMonthly_vertical_dimension
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, AM_regionalStatsMonthly, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_AM_regionalStatsMonthly_enable)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_regionalStatsMonthly_compute_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_regionalStatsMonthly_write_on_startup)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsMonthly_compute_interval)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsMonthly_output_stream)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsMonthly_restart_stream)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsMonthly_input_stream)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsMonthly_operation)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsMonthly_region_type)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsMonthly_region_group)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsMonthly_1d_weighting_function)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsMonthly_2d_weighting_function)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsMonthly_1d_weighting_field)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsMonthly_2d_weighting_field)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsMonthly_vertical_mask)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsMonthly_vertical_dimension)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record AM_regionalStatsMonthly')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_AM_regionalStatsMonthly_enable = $l', &
logicArgs=(/config_AM_regionalStatsMonthly_enable/))
            call mpas_log_write('        config_AM_regionalStatsMonthly_compute_on_startup = $l', &
logicArgs=(/config_AM_regionalStatsMonthly_compute_on_startup/))
            call mpas_log_write('        config_AM_regionalStatsMonthly_write_on_startup = $l', &
logicArgs=(/config_AM_regionalStatsMonthly_write_on_startup/))
            call mpas_log_write('        config_AM_regionalStatsMonthly_compute_interval = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsMonthly_compute_interval))
            call mpas_log_write('        config_AM_regionalStatsMonthly_output_stream = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsMonthly_output_stream))
            call mpas_log_write('        config_AM_regionalStatsMonthly_restart_stream = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsMonthly_restart_stream))
            call mpas_log_write('        config_AM_regionalStatsMonthly_input_stream = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsMonthly_input_stream))
            call mpas_log_write('        config_AM_regionalStatsMonthly_operation = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsMonthly_operation))
            call mpas_log_write('        config_AM_regionalStatsMonthly_region_type = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsMonthly_region_type))
            call mpas_log_write('        config_AM_regionalStatsMonthly_region_group = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsMonthly_region_group))
            call mpas_log_write('        config_AM_regionalStatsMonthly_1d_weighting_function = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsMonthly_1d_weighting_function))
            call mpas_log_write('        config_AM_regionalStatsMonthly_2d_weighting_function = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsMonthly_2d_weighting_function))
            call mpas_log_write('        config_AM_regionalStatsMonthly_1d_weighting_field = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsMonthly_1d_weighting_field))
            call mpas_log_write('        config_AM_regionalStatsMonthly_2d_weighting_field = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsMonthly_2d_weighting_field))
            call mpas_log_write('        config_AM_regionalStatsMonthly_vertical_mask = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsMonthly_vertical_mask))
            call mpas_log_write('        config_AM_regionalStatsMonthly_vertical_dimension = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsMonthly_vertical_dimension))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record AM_regionalStatsMonthly.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsMonthly_enable', config_AM_regionalStatsMonthly_enable)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsMonthly_compute_on_startup', &
config_AM_regionalStatsMonthly_compute_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsMonthly_write_on_startup', &
config_AM_regionalStatsMonthly_write_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsMonthly_compute_interval', &
config_AM_regionalStatsMonthly_compute_interval)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsMonthly_output_stream', &
config_AM_regionalStatsMonthly_output_stream)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsMonthly_restart_stream', &
config_AM_regionalStatsMonthly_restart_stream)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsMonthly_input_stream', &
config_AM_regionalStatsMonthly_input_stream)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsMonthly_operation', config_AM_regionalStatsMonthly_operation)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsMonthly_region_type', &
config_AM_regionalStatsMonthly_region_type)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsMonthly_region_group', &
config_AM_regionalStatsMonthly_region_group)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsMonthly_1d_weighting_function', &
config_AM_regionalStatsMonthly_1d_weighting_function)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsMonthly_2d_weighting_function', &
config_AM_regionalStatsMonthly_2d_weighting_function)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsMonthly_1d_weighting_field', &
config_AM_regionalStatsMonthly_1d_weighting_field)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsMonthly_2d_weighting_field', &
config_AM_regionalStatsMonthly_2d_weighting_field)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsMonthly_vertical_mask', &
config_AM_regionalStatsMonthly_vertical_mask)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsMonthly_vertical_dimension', &
config_AM_regionalStatsMonthly_vertical_dimension)

   end subroutine ocn_setup_nmlrec_AM_regionalStatsMonthly


   subroutine ocn_setup_nmlrec_AM_regionalStatsCustom(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_AM_regionalStatsCustom_enable = .false.
      logical :: config_AM_regionalStatsCustom_compute_on_startup = .false.
      logical :: config_AM_regionalStatsCustom_write_on_startup = .false.
      character (len=StrKIND) :: config_AM_regionalStatsCustom_compute_interval = 'output_interval'
      character (len=StrKIND) :: config_AM_regionalStatsCustom_output_stream = 'regionalStatsCustomOutput'
      character (len=StrKIND) :: config_AM_regionalStatsCustom_restart_stream = 'regionalMasksInput'
      character (len=StrKIND) :: config_AM_regionalStatsCustom_input_stream = 'regionalMasksInput'
      character (len=StrKIND) :: config_AM_regionalStatsCustom_operation = 'avg'
      character (len=StrKIND) :: config_AM_regionalStatsCustom_region_type = 'cell'
      character (len=StrKIND) :: config_AM_regionalStatsCustom_region_group = 'all'
      character (len=StrKIND) :: config_AM_regionalStatsCustom_1d_weighting_function = 'mul'
      character (len=StrKIND) :: config_AM_regionalStatsCustom_2d_weighting_function = 'mul'
      character (len=StrKIND) :: config_AM_regionalStatsCustom_1d_weighting_field = 'areaCell'
      character (len=StrKIND) :: config_AM_regionalStatsCustom_2d_weighting_field = 'volumeCell'
      character (len=StrKIND) :: config_AM_regionalStatsCustom_vertical_mask = 'cellMask'
      character (len=StrKIND) :: config_AM_regionalStatsCustom_vertical_dimension = 'nVertLevels'

      namelist /AM_regionalStatsCustom/ &
         config_AM_regionalStatsCustom_enable, &
         config_AM_regionalStatsCustom_compute_on_startup, &
         config_AM_regionalStatsCustom_write_on_startup, &
         config_AM_regionalStatsCustom_compute_interval, &
         config_AM_regionalStatsCustom_output_stream, &
         config_AM_regionalStatsCustom_restart_stream, &
         config_AM_regionalStatsCustom_input_stream, &
         config_AM_regionalStatsCustom_operation, &
         config_AM_regionalStatsCustom_region_type, &
         config_AM_regionalStatsCustom_region_group, &
         config_AM_regionalStatsCustom_1d_weighting_function, &
         config_AM_regionalStatsCustom_2d_weighting_function, &
         config_AM_regionalStatsCustom_1d_weighting_field, &
         config_AM_regionalStatsCustom_2d_weighting_field, &
         config_AM_regionalStatsCustom_vertical_mask, &
         config_AM_regionalStatsCustom_vertical_dimension
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, AM_regionalStatsCustom, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_AM_regionalStatsCustom_enable)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_regionalStatsCustom_compute_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_regionalStatsCustom_write_on_startup)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsCustom_compute_interval)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsCustom_output_stream)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsCustom_restart_stream)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsCustom_input_stream)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsCustom_operation)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsCustom_region_type)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsCustom_region_group)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsCustom_1d_weighting_function)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsCustom_2d_weighting_function)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsCustom_1d_weighting_field)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsCustom_2d_weighting_field)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsCustom_vertical_mask)
         call mpas_dmpar_bcast_char(dminfo, config_AM_regionalStatsCustom_vertical_dimension)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record AM_regionalStatsCustom')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_AM_regionalStatsCustom_enable = $l', &
logicArgs=(/config_AM_regionalStatsCustom_enable/))
            call mpas_log_write('        config_AM_regionalStatsCustom_compute_on_startup = $l', &
logicArgs=(/config_AM_regionalStatsCustom_compute_on_startup/))
            call mpas_log_write('        config_AM_regionalStatsCustom_write_on_startup = $l', &
logicArgs=(/config_AM_regionalStatsCustom_write_on_startup/))
            call mpas_log_write('        config_AM_regionalStatsCustom_compute_interval = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsCustom_compute_interval))
            call mpas_log_write('        config_AM_regionalStatsCustom_output_stream = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsCustom_output_stream))
            call mpas_log_write('        config_AM_regionalStatsCustom_restart_stream = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsCustom_restart_stream))
            call mpas_log_write('        config_AM_regionalStatsCustom_input_stream = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsCustom_input_stream))
            call mpas_log_write('        config_AM_regionalStatsCustom_operation = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsCustom_operation))
            call mpas_log_write('        config_AM_regionalStatsCustom_region_type = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsCustom_region_type))
            call mpas_log_write('        config_AM_regionalStatsCustom_region_group = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsCustom_region_group))
            call mpas_log_write('        config_AM_regionalStatsCustom_1d_weighting_function = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsCustom_1d_weighting_function))
            call mpas_log_write('        config_AM_regionalStatsCustom_2d_weighting_function = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsCustom_2d_weighting_function))
            call mpas_log_write('        config_AM_regionalStatsCustom_1d_weighting_field = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsCustom_1d_weighting_field))
            call mpas_log_write('        config_AM_regionalStatsCustom_2d_weighting_field = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsCustom_2d_weighting_field))
            call mpas_log_write('        config_AM_regionalStatsCustom_vertical_mask = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsCustom_vertical_mask))
            call mpas_log_write('        config_AM_regionalStatsCustom_vertical_dimension = '&
//''//mpas_log_escape_dollars(config_AM_regionalStatsCustom_vertical_dimension))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record AM_regionalStatsCustom.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsCustom_enable', config_AM_regionalStatsCustom_enable)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsCustom_compute_on_startup', &
config_AM_regionalStatsCustom_compute_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsCustom_write_on_startup', &
config_AM_regionalStatsCustom_write_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsCustom_compute_interval', &
config_AM_regionalStatsCustom_compute_interval)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsCustom_output_stream', &
config_AM_regionalStatsCustom_output_stream)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsCustom_restart_stream', &
config_AM_regionalStatsCustom_restart_stream)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsCustom_input_stream', &
config_AM_regionalStatsCustom_input_stream)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsCustom_operation', config_AM_regionalStatsCustom_operation)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsCustom_region_type', config_AM_regionalStatsCustom_region_type)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsCustom_region_group', &
config_AM_regionalStatsCustom_region_group)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsCustom_1d_weighting_function', &
config_AM_regionalStatsCustom_1d_weighting_function)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsCustom_2d_weighting_function', &
config_AM_regionalStatsCustom_2d_weighting_function)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsCustom_1d_weighting_field', &
config_AM_regionalStatsCustom_1d_weighting_field)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsCustom_2d_weighting_field', &
config_AM_regionalStatsCustom_2d_weighting_field)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsCustom_vertical_mask', &
config_AM_regionalStatsCustom_vertical_mask)
      call mpas_pool_add_config(configPool, 'config_AM_regionalStatsCustom_vertical_dimension', &
config_AM_regionalStatsCustom_vertical_dimension)

   end subroutine ocn_setup_nmlrec_AM_regionalStatsCustom


   subroutine ocn_setup_nmlrec_AM_timeSeriesStatsDaily(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_AM_timeSeriesStatsDaily_enable = .false.
      logical :: config_AM_timeSeriesStatsDaily_compute_on_startup = .false.
      logical :: config_AM_timeSeriesStatsDaily_write_on_startup = .false.
      character (len=StrKIND) :: config_AM_timeSeriesStatsDaily_compute_interval = '00-00-00_01:00:00'
      character (len=StrKIND) :: config_AM_timeSeriesStatsDaily_output_stream = 'timeSeriesStatsDailyOutput'
      character (len=StrKIND) :: config_AM_timeSeriesStatsDaily_restart_stream = 'timeSeriesStatsDailyRestart'
      character (len=StrKIND) :: config_AM_timeSeriesStatsDaily_operation = 'avg'
      character (len=StrKIND) :: config_AM_timeSeriesStatsDaily_reference_times = 'initial_time'
      character (len=StrKIND) :: config_AM_timeSeriesStatsDaily_duration_intervals = 'repeat_interval'
      character (len=StrKIND) :: config_AM_timeSeriesStatsDaily_repeat_intervals = 'reset_interval'
      character (len=StrKIND) :: config_AM_timeSeriesStatsDaily_reset_intervals = '00-00-01_00:00:00'
      character (len=StrKIND) :: config_AM_timeSeriesStatsDaily_backward_output_offset = '00-00-01_00:00:00'

      namelist /AM_timeSeriesStatsDaily/ &
         config_AM_timeSeriesStatsDaily_enable, &
         config_AM_timeSeriesStatsDaily_compute_on_startup, &
         config_AM_timeSeriesStatsDaily_write_on_startup, &
         config_AM_timeSeriesStatsDaily_compute_interval, &
         config_AM_timeSeriesStatsDaily_output_stream, &
         config_AM_timeSeriesStatsDaily_restart_stream, &
         config_AM_timeSeriesStatsDaily_operation, &
         config_AM_timeSeriesStatsDaily_reference_times, &
         config_AM_timeSeriesStatsDaily_duration_intervals, &
         config_AM_timeSeriesStatsDaily_repeat_intervals, &
         config_AM_timeSeriesStatsDaily_reset_intervals, &
         config_AM_timeSeriesStatsDaily_backward_output_offset
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, AM_timeSeriesStatsDaily, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_AM_timeSeriesStatsDaily_enable)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_timeSeriesStatsDaily_compute_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_timeSeriesStatsDaily_write_on_startup)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsDaily_compute_interval)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsDaily_output_stream)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsDaily_restart_stream)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsDaily_operation)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsDaily_reference_times)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsDaily_duration_intervals)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsDaily_repeat_intervals)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsDaily_reset_intervals)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsDaily_backward_output_offset)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record AM_timeSeriesStatsDaily')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_AM_timeSeriesStatsDaily_enable = $l', &
logicArgs=(/config_AM_timeSeriesStatsDaily_enable/))
            call mpas_log_write('        config_AM_timeSeriesStatsDaily_compute_on_startup = $l', &
logicArgs=(/config_AM_timeSeriesStatsDaily_compute_on_startup/))
            call mpas_log_write('        config_AM_timeSeriesStatsDaily_write_on_startup = $l', &
logicArgs=(/config_AM_timeSeriesStatsDaily_write_on_startup/))
            call mpas_log_write('        config_AM_timeSeriesStatsDaily_compute_interval = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsDaily_compute_interval))
            call mpas_log_write('        config_AM_timeSeriesStatsDaily_output_stream = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsDaily_output_stream))
            call mpas_log_write('        config_AM_timeSeriesStatsDaily_restart_stream = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsDaily_restart_stream))
            call mpas_log_write('        config_AM_timeSeriesStatsDaily_operation = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsDaily_operation))
            call mpas_log_write('        config_AM_timeSeriesStatsDaily_reference_times = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsDaily_reference_times))
            call mpas_log_write('        config_AM_timeSeriesStatsDaily_duration_intervals = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsDaily_duration_intervals))
            call mpas_log_write('        config_AM_timeSeriesStatsDaily_repeat_intervals = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsDaily_repeat_intervals))
            call mpas_log_write('        config_AM_timeSeriesStatsDaily_reset_intervals = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsDaily_reset_intervals))
            call mpas_log_write('        config_AM_timeSeriesStatsDaily_backward_output_offset = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsDaily_backward_output_offset))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record AM_timeSeriesStatsDaily.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsDaily_enable', config_AM_timeSeriesStatsDaily_enable)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsDaily_compute_on_startup', &
config_AM_timeSeriesStatsDaily_compute_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsDaily_write_on_startup', &
config_AM_timeSeriesStatsDaily_write_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsDaily_compute_interval', &
config_AM_timeSeriesStatsDaily_compute_interval)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsDaily_output_stream', &
config_AM_timeSeriesStatsDaily_output_stream)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsDaily_restart_stream', &
config_AM_timeSeriesStatsDaily_restart_stream)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsDaily_operation', config_AM_timeSeriesStatsDaily_operation)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsDaily_reference_times', &
config_AM_timeSeriesStatsDaily_reference_times)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsDaily_duration_intervals', &
config_AM_timeSeriesStatsDaily_duration_intervals)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsDaily_repeat_intervals', &
config_AM_timeSeriesStatsDaily_repeat_intervals)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsDaily_reset_intervals', &
config_AM_timeSeriesStatsDaily_reset_intervals)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsDaily_backward_output_offset', &
config_AM_timeSeriesStatsDaily_backward_output_offset)

   end subroutine ocn_setup_nmlrec_AM_timeSeriesStatsDaily


   subroutine ocn_setup_nmlrec_AM_timeSeriesStatsMonthly(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_AM_timeSeriesStatsMonthly_enable = .false.
      logical :: config_AM_timeSeriesStatsMonthly_compute_on_startup = .false.
      logical :: config_AM_timeSeriesStatsMonthly_write_on_startup = .false.
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthly_compute_interval = '00-00-00_01:00:00'
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthly_output_stream = 'timeSeriesStatsMonthlyOutput'
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthly_restart_stream = 'timeSeriesStatsMonthlyRestart'
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthly_operation = 'avg'
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthly_reference_times = 'initial_time'
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthly_duration_intervals = 'repeat_interval'
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthly_repeat_intervals = 'reset_interval'
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthly_reset_intervals = '00-01-00_00:00:00'
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthly_backward_output_offset = '00-01-00_00:00:00'

      namelist /AM_timeSeriesStatsMonthly/ &
         config_AM_timeSeriesStatsMonthly_enable, &
         config_AM_timeSeriesStatsMonthly_compute_on_startup, &
         config_AM_timeSeriesStatsMonthly_write_on_startup, &
         config_AM_timeSeriesStatsMonthly_compute_interval, &
         config_AM_timeSeriesStatsMonthly_output_stream, &
         config_AM_timeSeriesStatsMonthly_restart_stream, &
         config_AM_timeSeriesStatsMonthly_operation, &
         config_AM_timeSeriesStatsMonthly_reference_times, &
         config_AM_timeSeriesStatsMonthly_duration_intervals, &
         config_AM_timeSeriesStatsMonthly_repeat_intervals, &
         config_AM_timeSeriesStatsMonthly_reset_intervals, &
         config_AM_timeSeriesStatsMonthly_backward_output_offset
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, AM_timeSeriesStatsMonthly, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_AM_timeSeriesStatsMonthly_enable)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_timeSeriesStatsMonthly_compute_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_timeSeriesStatsMonthly_write_on_startup)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsMonthly_compute_interval)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsMonthly_output_stream)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsMonthly_restart_stream)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsMonthly_operation)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsMonthly_reference_times)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsMonthly_duration_intervals)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsMonthly_repeat_intervals)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsMonthly_reset_intervals)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsMonthly_backward_output_offset)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record AM_timeSeriesStatsMonthly')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_AM_timeSeriesStatsMonthly_enable = $l', &
logicArgs=(/config_AM_timeSeriesStatsMonthly_enable/))
            call mpas_log_write('        config_AM_timeSeriesStatsMonthly_compute_on_startup = $l', &
logicArgs=(/config_AM_timeSeriesStatsMonthly_compute_on_startup/))
            call mpas_log_write('        config_AM_timeSeriesStatsMonthly_write_on_startup = $l', &
logicArgs=(/config_AM_timeSeriesStatsMonthly_write_on_startup/))
            call mpas_log_write('        config_AM_timeSeriesStatsMonthly_compute_interval = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsMonthly_compute_interval))
            call mpas_log_write('        config_AM_timeSeriesStatsMonthly_output_stream = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsMonthly_output_stream))
            call mpas_log_write('        config_AM_timeSeriesStatsMonthly_restart_stream = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsMonthly_restart_stream))
            call mpas_log_write('        config_AM_timeSeriesStatsMonthly_operation = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsMonthly_operation))
            call mpas_log_write('        config_AM_timeSeriesStatsMonthly_reference_times = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsMonthly_reference_times))
            call mpas_log_write('        config_AM_timeSeriesStatsMonthly_duration_intervals = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsMonthly_duration_intervals))
            call mpas_log_write('        config_AM_timeSeriesStatsMonthly_repeat_intervals = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsMonthly_repeat_intervals))
            call mpas_log_write('        config_AM_timeSeriesStatsMonthly_reset_intervals = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsMonthly_reset_intervals))
            call mpas_log_write('        config_AM_timeSeriesStatsMonthly_backward_output_offset = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsMonthly_backward_output_offset))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record AM_timeSeriesStatsMonthly.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsMonthly_enable', config_AM_timeSeriesStatsMonthly_enable)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsMonthly_compute_on_startup', &
config_AM_timeSeriesStatsMonthly_compute_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsMonthly_write_on_startup', &
config_AM_timeSeriesStatsMonthly_write_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsMonthly_compute_interval', &
config_AM_timeSeriesStatsMonthly_compute_interval)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsMonthly_output_stream', &
config_AM_timeSeriesStatsMonthly_output_stream)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsMonthly_restart_stream', &
config_AM_timeSeriesStatsMonthly_restart_stream)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsMonthly_operation', &
config_AM_timeSeriesStatsMonthly_operation)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsMonthly_reference_times', &
config_AM_timeSeriesStatsMonthly_reference_times)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsMonthly_duration_intervals', &
config_AM_timeSeriesStatsMonthly_duration_intervals)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsMonthly_repeat_intervals', &
config_AM_timeSeriesStatsMonthly_repeat_intervals)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsMonthly_reset_intervals', &
config_AM_timeSeriesStatsMonthly_reset_intervals)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsMonthly_backward_output_offset', &
config_AM_timeSeriesStatsMonthly_backward_output_offset)

   end subroutine ocn_setup_nmlrec_AM_timeSeriesStatsMonthly


   subroutine ocn_setup_nmlrec_AM_timeSeriesStatsClimatology(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_AM_timeSeriesStatsClimatology_enable = .false.
      logical :: config_AM_timeSeriesStatsClimatology_compute_on_startup = .false.
      logical :: config_AM_timeSeriesStatsClimatology_write_on_startup = .false.
      character (len=StrKIND) :: config_AM_timeSeriesStatsClimatology_compute_interval = '00-00-00_01:00:00'
      character (len=StrKIND) :: config_AM_timeSeriesStatsClimatology_output_stream = 'timeSeriesStatsClimatologyOutput'
      character (len=StrKIND) :: config_AM_timeSeriesStatsClimatology_restart_stream = 'timeSeriesStatsClimatologyRestart'
      character (len=StrKIND) :: config_AM_timeSeriesStatsClimatology_operation = 'avg'
      character (len=StrKIND) :: config_AM_timeSeriesStatsClimatology_reference_times = &
'00-03-01_00:00:00;00-06-01_00:00:00;00-09-01_00:00:00;00-12-01_00:00:00'
      character (len=StrKIND) :: config_AM_timeSeriesStatsClimatology_duration_intervals = &
'00-03-00_00:00:00;00-03-00_00:00:00;00-03-00_00:00:00;00-03-00_00:00:00'
      character (len=StrKIND) :: config_AM_timeSeriesStatsClimatology_repeat_intervals = &
'01-00-00_00:00:00;01-00-00_00:00:00;01-00-00_00:00:00;01-00-00_00:00:00'
      character (len=StrKIND) :: config_AM_timeSeriesStatsClimatology_reset_intervals = &
'1000-00-00_00:00:00;1000-00-00_00:00:00;1000-00-00_00:00:00;1000-00-00_00:00:00'
      character (len=StrKIND) :: config_AM_timeSeriesStatsClimatology_backward_output_offset = '00-03-00_00:00:00'

      namelist /AM_timeSeriesStatsClimatology/ &
         config_AM_timeSeriesStatsClimatology_enable, &
         config_AM_timeSeriesStatsClimatology_compute_on_startup, &
         config_AM_timeSeriesStatsClimatology_write_on_startup, &
         config_AM_timeSeriesStatsClimatology_compute_interval, &
         config_AM_timeSeriesStatsClimatology_output_stream, &
         config_AM_timeSeriesStatsClimatology_restart_stream, &
         config_AM_timeSeriesStatsClimatology_operation, &
         config_AM_timeSeriesStatsClimatology_reference_times, &
         config_AM_timeSeriesStatsClimatology_duration_intervals, &
         config_AM_timeSeriesStatsClimatology_repeat_intervals, &
         config_AM_timeSeriesStatsClimatology_reset_intervals, &
         config_AM_timeSeriesStatsClimatology_backward_output_offset
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, AM_timeSeriesStatsClimatology, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_AM_timeSeriesStatsClimatology_enable)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_timeSeriesStatsClimatology_compute_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_timeSeriesStatsClimatology_write_on_startup)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsClimatology_compute_interval)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsClimatology_output_stream)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsClimatology_restart_stream)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsClimatology_operation)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsClimatology_reference_times)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsClimatology_duration_intervals)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsClimatology_repeat_intervals)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsClimatology_reset_intervals)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsClimatology_backward_output_offset)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record AM_timeSeriesStatsClimatology')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_AM_timeSeriesStatsClimatology_enable = $l', &
logicArgs=(/config_AM_timeSeriesStatsClimatology_enable/))
            call mpas_log_write('        config_AM_timeSeriesStatsClimatology_compute_on_startup = $l', &
logicArgs=(/config_AM_timeSeriesStatsClimatology_compute_on_startup/))
            call mpas_log_write('        config_AM_timeSeriesStatsClimatology_write_on_startup = $l', &
logicArgs=(/config_AM_timeSeriesStatsClimatology_write_on_startup/))
            call mpas_log_write('        config_AM_timeSeriesStatsClimatology_compute_interval = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsClimatology_compute_interval))
            call mpas_log_write('        config_AM_timeSeriesStatsClimatology_output_stream = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsClimatology_output_stream))
            call mpas_log_write('        config_AM_timeSeriesStatsClimatology_restart_stream = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsClimatology_restart_stream))
            call mpas_log_write('        config_AM_timeSeriesStatsClimatology_operation = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsClimatology_operation))
            call mpas_log_write('        config_AM_timeSeriesStatsClimatology_reference_times = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsClimatology_reference_times))
            call mpas_log_write('        config_AM_timeSeriesStatsClimatology_duration_intervals = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsClimatology_duration_intervals))
            call mpas_log_write('        config_AM_timeSeriesStatsClimatology_repeat_intervals = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsClimatology_repeat_intervals))
            call mpas_log_write('        config_AM_timeSeriesStatsClimatology_reset_intervals = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsClimatology_reset_intervals))
            call mpas_log_write('        config_AM_timeSeriesStatsClimatology_backward_output_offset = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsClimatology_backward_output_offset))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record AM_timeSeriesStatsClimatology.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsClimatology_enable', &
config_AM_timeSeriesStatsClimatology_enable)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsClimatology_compute_on_startup', &
config_AM_timeSeriesStatsClimatology_compute_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsClimatology_write_on_startup', &
config_AM_timeSeriesStatsClimatology_write_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsClimatology_compute_interval', &
config_AM_timeSeriesStatsClimatology_compute_interval)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsClimatology_output_stream', &
config_AM_timeSeriesStatsClimatology_output_stream)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsClimatology_restart_stream', &
config_AM_timeSeriesStatsClimatology_restart_stream)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsClimatology_operation', &
config_AM_timeSeriesStatsClimatology_operation)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsClimatology_reference_times', &
config_AM_timeSeriesStatsClimatology_reference_times)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsClimatology_duration_intervals', &
config_AM_timeSeriesStatsClimatology_duration_intervals)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsClimatology_repeat_intervals', &
config_AM_timeSeriesStatsClimatology_repeat_intervals)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsClimatology_reset_intervals', &
config_AM_timeSeriesStatsClimatology_reset_intervals)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsClimatology_backward_output_offset', &
config_AM_timeSeriesStatsClimatology_backward_output_offset)

   end subroutine ocn_setup_nmlrec_AM_timeSeriesStatsClimatology


   subroutine ocn_setup_nmlrec_AM_timeSeriesStatsMonthlyMax(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_AM_timeSeriesStatsMonthlyMax_enable = .false.
      logical :: config_AM_timeSeriesStatsMonthlyMax_compute_on_startup = .false.
      logical :: config_AM_timeSeriesStatsMonthlyMax_write_on_startup = .false.
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthlyMax_compute_interval = '00-00-00_01:00:00'
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthlyMax_output_stream = 'timeSeriesStatsMonthlyMaxOutput'
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthlyMax_restart_stream = 'timeSeriesStatsMonthlyMaxRestart'
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthlyMax_operation = 'max'
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthlyMax_reference_times = 'initial_time'
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthlyMax_duration_intervals = 'repeat_interval'
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthlyMax_repeat_intervals = 'reset_interval'
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthlyMax_reset_intervals = '00-01-00_00:00:00'
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthlyMax_backward_output_offset = '00-01-00_00:00:00'

      namelist /AM_timeSeriesStatsMonthlyMax/ &
         config_AM_timeSeriesStatsMonthlyMax_enable, &
         config_AM_timeSeriesStatsMonthlyMax_compute_on_startup, &
         config_AM_timeSeriesStatsMonthlyMax_write_on_startup, &
         config_AM_timeSeriesStatsMonthlyMax_compute_interval, &
         config_AM_timeSeriesStatsMonthlyMax_output_stream, &
         config_AM_timeSeriesStatsMonthlyMax_restart_stream, &
         config_AM_timeSeriesStatsMonthlyMax_operation, &
         config_AM_timeSeriesStatsMonthlyMax_reference_times, &
         config_AM_timeSeriesStatsMonthlyMax_duration_intervals, &
         config_AM_timeSeriesStatsMonthlyMax_repeat_intervals, &
         config_AM_timeSeriesStatsMonthlyMax_reset_intervals, &
         config_AM_timeSeriesStatsMonthlyMax_backward_output_offset
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, AM_timeSeriesStatsMonthlyMax, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_AM_timeSeriesStatsMonthlyMax_enable)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_timeSeriesStatsMonthlyMax_compute_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_timeSeriesStatsMonthlyMax_write_on_startup)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsMonthlyMax_compute_interval)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsMonthlyMax_output_stream)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsMonthlyMax_restart_stream)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsMonthlyMax_operation)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsMonthlyMax_reference_times)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsMonthlyMax_duration_intervals)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsMonthlyMax_repeat_intervals)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsMonthlyMax_reset_intervals)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsMonthlyMax_backward_output_offset)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record AM_timeSeriesStatsMonthlyMax')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_AM_timeSeriesStatsMonthlyMax_enable = $l', &
logicArgs=(/config_AM_timeSeriesStatsMonthlyMax_enable/))
            call mpas_log_write('        config_AM_timeSeriesStatsMonthlyMax_compute_on_startup = $l', &
logicArgs=(/config_AM_timeSeriesStatsMonthlyMax_compute_on_startup/))
            call mpas_log_write('        config_AM_timeSeriesStatsMonthlyMax_write_on_startup = $l', &
logicArgs=(/config_AM_timeSeriesStatsMonthlyMax_write_on_startup/))
            call mpas_log_write('        config_AM_timeSeriesStatsMonthlyMax_compute_interval = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsMonthlyMax_compute_interval))
            call mpas_log_write('        config_AM_timeSeriesStatsMonthlyMax_output_stream = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsMonthlyMax_output_stream))
            call mpas_log_write('        config_AM_timeSeriesStatsMonthlyMax_restart_stream = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsMonthlyMax_restart_stream))
            call mpas_log_write('        config_AM_timeSeriesStatsMonthlyMax_operation = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsMonthlyMax_operation))
            call mpas_log_write('        config_AM_timeSeriesStatsMonthlyMax_reference_times = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsMonthlyMax_reference_times))
            call mpas_log_write('        config_AM_timeSeriesStatsMonthlyMax_duration_intervals = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsMonthlyMax_duration_intervals))
            call mpas_log_write('        config_AM_timeSeriesStatsMonthlyMax_repeat_intervals = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsMonthlyMax_repeat_intervals))
            call mpas_log_write('        config_AM_timeSeriesStatsMonthlyMax_reset_intervals = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsMonthlyMax_reset_intervals))
            call mpas_log_write('        config_AM_timeSeriesStatsMonthlyMax_backward_output_offset = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsMonthlyMax_backward_output_offset))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record AM_timeSeriesStatsMonthlyMax.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsMonthlyMax_enable', &
config_AM_timeSeriesStatsMonthlyMax_enable)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsMonthlyMax_compute_on_startup', &
config_AM_timeSeriesStatsMonthlyMax_compute_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsMonthlyMax_write_on_startup', &
config_AM_timeSeriesStatsMonthlyMax_write_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsMonthlyMax_compute_interval', &
config_AM_timeSeriesStatsMonthlyMax_compute_interval)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsMonthlyMax_output_stream', &
config_AM_timeSeriesStatsMonthlyMax_output_stream)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsMonthlyMax_restart_stream', &
config_AM_timeSeriesStatsMonthlyMax_restart_stream)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsMonthlyMax_operation', &
config_AM_timeSeriesStatsMonthlyMax_operation)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsMonthlyMax_reference_times', &
config_AM_timeSeriesStatsMonthlyMax_reference_times)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsMonthlyMax_duration_intervals', &
config_AM_timeSeriesStatsMonthlyMax_duration_intervals)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsMonthlyMax_repeat_intervals', &
config_AM_timeSeriesStatsMonthlyMax_repeat_intervals)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsMonthlyMax_reset_intervals', &
config_AM_timeSeriesStatsMonthlyMax_reset_intervals)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsMonthlyMax_backward_output_offset', &
config_AM_timeSeriesStatsMonthlyMax_backward_output_offset)

   end subroutine ocn_setup_nmlrec_AM_timeSeriesStatsMonthlyMax


   subroutine ocn_setup_nmlrec_AM_timeSeriesStatsMonthlyMin(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_AM_timeSeriesStatsMonthlyMin_enable = .false.
      logical :: config_AM_timeSeriesStatsMonthlyMin_compute_on_startup = .false.
      logical :: config_AM_timeSeriesStatsMonthlyMin_write_on_startup = .false.
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthlyMin_compute_interval = '00-00-00_01:00:00'
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthlyMin_output_stream = 'timeSeriesStatsMonthlyMinOutput'
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthlyMin_restart_stream = 'timeSeriesStatsMonthlyMinRestart'
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthlyMin_operation = 'min'
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthlyMin_reference_times = 'initial_time'
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthlyMin_duration_intervals = 'repeat_interval'
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthlyMin_repeat_intervals = 'reset_interval'
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthlyMin_reset_intervals = '00-01-00_00:00:00'
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthlyMin_backward_output_offset = '00-01-00_00:00:00'

      namelist /AM_timeSeriesStatsMonthlyMin/ &
         config_AM_timeSeriesStatsMonthlyMin_enable, &
         config_AM_timeSeriesStatsMonthlyMin_compute_on_startup, &
         config_AM_timeSeriesStatsMonthlyMin_write_on_startup, &
         config_AM_timeSeriesStatsMonthlyMin_compute_interval, &
         config_AM_timeSeriesStatsMonthlyMin_output_stream, &
         config_AM_timeSeriesStatsMonthlyMin_restart_stream, &
         config_AM_timeSeriesStatsMonthlyMin_operation, &
         config_AM_timeSeriesStatsMonthlyMin_reference_times, &
         config_AM_timeSeriesStatsMonthlyMin_duration_intervals, &
         config_AM_timeSeriesStatsMonthlyMin_repeat_intervals, &
         config_AM_timeSeriesStatsMonthlyMin_reset_intervals, &
         config_AM_timeSeriesStatsMonthlyMin_backward_output_offset
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, AM_timeSeriesStatsMonthlyMin, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_AM_timeSeriesStatsMonthlyMin_enable)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_timeSeriesStatsMonthlyMin_compute_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_timeSeriesStatsMonthlyMin_write_on_startup)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsMonthlyMin_compute_interval)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsMonthlyMin_output_stream)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsMonthlyMin_restart_stream)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsMonthlyMin_operation)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsMonthlyMin_reference_times)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsMonthlyMin_duration_intervals)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsMonthlyMin_repeat_intervals)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsMonthlyMin_reset_intervals)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsMonthlyMin_backward_output_offset)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record AM_timeSeriesStatsMonthlyMin')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_AM_timeSeriesStatsMonthlyMin_enable = $l', &
logicArgs=(/config_AM_timeSeriesStatsMonthlyMin_enable/))
            call mpas_log_write('        config_AM_timeSeriesStatsMonthlyMin_compute_on_startup = $l', &
logicArgs=(/config_AM_timeSeriesStatsMonthlyMin_compute_on_startup/))
            call mpas_log_write('        config_AM_timeSeriesStatsMonthlyMin_write_on_startup = $l', &
logicArgs=(/config_AM_timeSeriesStatsMonthlyMin_write_on_startup/))
            call mpas_log_write('        config_AM_timeSeriesStatsMonthlyMin_compute_interval = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsMonthlyMin_compute_interval))
            call mpas_log_write('        config_AM_timeSeriesStatsMonthlyMin_output_stream = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsMonthlyMin_output_stream))
            call mpas_log_write('        config_AM_timeSeriesStatsMonthlyMin_restart_stream = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsMonthlyMin_restart_stream))
            call mpas_log_write('        config_AM_timeSeriesStatsMonthlyMin_operation = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsMonthlyMin_operation))
            call mpas_log_write('        config_AM_timeSeriesStatsMonthlyMin_reference_times = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsMonthlyMin_reference_times))
            call mpas_log_write('        config_AM_timeSeriesStatsMonthlyMin_duration_intervals = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsMonthlyMin_duration_intervals))
            call mpas_log_write('        config_AM_timeSeriesStatsMonthlyMin_repeat_intervals = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsMonthlyMin_repeat_intervals))
            call mpas_log_write('        config_AM_timeSeriesStatsMonthlyMin_reset_intervals = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsMonthlyMin_reset_intervals))
            call mpas_log_write('        config_AM_timeSeriesStatsMonthlyMin_backward_output_offset = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsMonthlyMin_backward_output_offset))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record AM_timeSeriesStatsMonthlyMin.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsMonthlyMin_enable', &
config_AM_timeSeriesStatsMonthlyMin_enable)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsMonthlyMin_compute_on_startup', &
config_AM_timeSeriesStatsMonthlyMin_compute_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsMonthlyMin_write_on_startup', &
config_AM_timeSeriesStatsMonthlyMin_write_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsMonthlyMin_compute_interval', &
config_AM_timeSeriesStatsMonthlyMin_compute_interval)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsMonthlyMin_output_stream', &
config_AM_timeSeriesStatsMonthlyMin_output_stream)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsMonthlyMin_restart_stream', &
config_AM_timeSeriesStatsMonthlyMin_restart_stream)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsMonthlyMin_operation', &
config_AM_timeSeriesStatsMonthlyMin_operation)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsMonthlyMin_reference_times', &
config_AM_timeSeriesStatsMonthlyMin_reference_times)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsMonthlyMin_duration_intervals', &
config_AM_timeSeriesStatsMonthlyMin_duration_intervals)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsMonthlyMin_repeat_intervals', &
config_AM_timeSeriesStatsMonthlyMin_repeat_intervals)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsMonthlyMin_reset_intervals', &
config_AM_timeSeriesStatsMonthlyMin_reset_intervals)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsMonthlyMin_backward_output_offset', &
config_AM_timeSeriesStatsMonthlyMin_backward_output_offset)

   end subroutine ocn_setup_nmlrec_AM_timeSeriesStatsMonthlyMin


   subroutine ocn_setup_nmlrec_AM_timeSeriesStatsCustom(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_AM_timeSeriesStatsCustom_enable = .false.
      logical :: config_AM_timeSeriesStatsCustom_compute_on_startup = .false.
      logical :: config_AM_timeSeriesStatsCustom_write_on_startup = .false.
      character (len=StrKIND) :: config_AM_timeSeriesStatsCustom_compute_interval = '00-00-00_01:00:00'
      character (len=StrKIND) :: config_AM_timeSeriesStatsCustom_output_stream = 'timeSeriesStatsCustomOutput'
      character (len=StrKIND) :: config_AM_timeSeriesStatsCustom_restart_stream = 'timeSeriesStatsCustomRestart'
      character (len=StrKIND) :: config_AM_timeSeriesStatsCustom_operation = 'avg'
      character (len=StrKIND) :: config_AM_timeSeriesStatsCustom_reference_times = 'initial_time'
      character (len=StrKIND) :: config_AM_timeSeriesStatsCustom_duration_intervals = 'repeat_interval'
      character (len=StrKIND) :: config_AM_timeSeriesStatsCustom_repeat_intervals = 'reset_interval'
      character (len=StrKIND) :: config_AM_timeSeriesStatsCustom_reset_intervals = '00-00-07_00:00:00'
      character (len=StrKIND) :: config_AM_timeSeriesStatsCustom_backward_output_offset = '00-00-01_00:00:00'

      namelist /AM_timeSeriesStatsCustom/ &
         config_AM_timeSeriesStatsCustom_enable, &
         config_AM_timeSeriesStatsCustom_compute_on_startup, &
         config_AM_timeSeriesStatsCustom_write_on_startup, &
         config_AM_timeSeriesStatsCustom_compute_interval, &
         config_AM_timeSeriesStatsCustom_output_stream, &
         config_AM_timeSeriesStatsCustom_restart_stream, &
         config_AM_timeSeriesStatsCustom_operation, &
         config_AM_timeSeriesStatsCustom_reference_times, &
         config_AM_timeSeriesStatsCustom_duration_intervals, &
         config_AM_timeSeriesStatsCustom_repeat_intervals, &
         config_AM_timeSeriesStatsCustom_reset_intervals, &
         config_AM_timeSeriesStatsCustom_backward_output_offset
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, AM_timeSeriesStatsCustom, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_AM_timeSeriesStatsCustom_enable)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_timeSeriesStatsCustom_compute_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_timeSeriesStatsCustom_write_on_startup)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsCustom_compute_interval)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsCustom_output_stream)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsCustom_restart_stream)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsCustom_operation)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsCustom_reference_times)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsCustom_duration_intervals)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsCustom_repeat_intervals)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsCustom_reset_intervals)
         call mpas_dmpar_bcast_char(dminfo, config_AM_timeSeriesStatsCustom_backward_output_offset)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record AM_timeSeriesStatsCustom')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_AM_timeSeriesStatsCustom_enable = $l', &
logicArgs=(/config_AM_timeSeriesStatsCustom_enable/))
            call mpas_log_write('        config_AM_timeSeriesStatsCustom_compute_on_startup = $l', &
logicArgs=(/config_AM_timeSeriesStatsCustom_compute_on_startup/))
            call mpas_log_write('        config_AM_timeSeriesStatsCustom_write_on_startup = $l', &
logicArgs=(/config_AM_timeSeriesStatsCustom_write_on_startup/))
            call mpas_log_write('        config_AM_timeSeriesStatsCustom_compute_interval = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsCustom_compute_interval))
            call mpas_log_write('        config_AM_timeSeriesStatsCustom_output_stream = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsCustom_output_stream))
            call mpas_log_write('        config_AM_timeSeriesStatsCustom_restart_stream = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsCustom_restart_stream))
            call mpas_log_write('        config_AM_timeSeriesStatsCustom_operation = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsCustom_operation))
            call mpas_log_write('        config_AM_timeSeriesStatsCustom_reference_times = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsCustom_reference_times))
            call mpas_log_write('        config_AM_timeSeriesStatsCustom_duration_intervals = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsCustom_duration_intervals))
            call mpas_log_write('        config_AM_timeSeriesStatsCustom_repeat_intervals = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsCustom_repeat_intervals))
            call mpas_log_write('        config_AM_timeSeriesStatsCustom_reset_intervals = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsCustom_reset_intervals))
            call mpas_log_write('        config_AM_timeSeriesStatsCustom_backward_output_offset = '&
//''//mpas_log_escape_dollars(config_AM_timeSeriesStatsCustom_backward_output_offset))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record AM_timeSeriesStatsCustom.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsCustom_enable', config_AM_timeSeriesStatsCustom_enable)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsCustom_compute_on_startup', &
config_AM_timeSeriesStatsCustom_compute_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsCustom_write_on_startup', &
config_AM_timeSeriesStatsCustom_write_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsCustom_compute_interval', &
config_AM_timeSeriesStatsCustom_compute_interval)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsCustom_output_stream', &
config_AM_timeSeriesStatsCustom_output_stream)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsCustom_restart_stream', &
config_AM_timeSeriesStatsCustom_restart_stream)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsCustom_operation', config_AM_timeSeriesStatsCustom_operation)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsCustom_reference_times', &
config_AM_timeSeriesStatsCustom_reference_times)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsCustom_duration_intervals', &
config_AM_timeSeriesStatsCustom_duration_intervals)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsCustom_repeat_intervals', &
config_AM_timeSeriesStatsCustom_repeat_intervals)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsCustom_reset_intervals', &
config_AM_timeSeriesStatsCustom_reset_intervals)
      call mpas_pool_add_config(configPool, 'config_AM_timeSeriesStatsCustom_backward_output_offset', &
config_AM_timeSeriesStatsCustom_backward_output_offset)

   end subroutine ocn_setup_nmlrec_AM_timeSeriesStatsCustom


   subroutine ocn_setup_nmlrec_AM_pointwiseStats(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_AM_pointwiseStats_enable = .false.
      character (len=StrKIND) :: config_AM_pointwiseStats_compute_interval = 'output_interval'
      character (len=StrKIND) :: config_AM_pointwiseStats_output_stream = 'pointwiseStatsOutput'
      logical :: config_AM_pointwiseStats_compute_on_startup = .true.
      logical :: config_AM_pointwiseStats_write_on_startup = .true.

      namelist /AM_pointwiseStats/ &
         config_AM_pointwiseStats_enable, &
         config_AM_pointwiseStats_compute_interval, &
         config_AM_pointwiseStats_output_stream, &
         config_AM_pointwiseStats_compute_on_startup, &
         config_AM_pointwiseStats_write_on_startup
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, AM_pointwiseStats, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_AM_pointwiseStats_enable)
         call mpas_dmpar_bcast_char(dminfo, config_AM_pointwiseStats_compute_interval)
         call mpas_dmpar_bcast_char(dminfo, config_AM_pointwiseStats_output_stream)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_pointwiseStats_compute_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_pointwiseStats_write_on_startup)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record AM_pointwiseStats')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_AM_pointwiseStats_enable = $l', logicArgs=(/config_AM_pointwiseStats_enable/))
            call mpas_log_write('        config_AM_pointwiseStats_compute_interval = '&
//''//mpas_log_escape_dollars(config_AM_pointwiseStats_compute_interval))
            call mpas_log_write('        config_AM_pointwiseStats_output_stream = '&
//''//mpas_log_escape_dollars(config_AM_pointwiseStats_output_stream))
            call mpas_log_write('        config_AM_pointwiseStats_compute_on_startup = $l', &
logicArgs=(/config_AM_pointwiseStats_compute_on_startup/))
            call mpas_log_write('        config_AM_pointwiseStats_write_on_startup = $l', &
logicArgs=(/config_AM_pointwiseStats_write_on_startup/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record AM_pointwiseStats.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_AM_pointwiseStats_enable', config_AM_pointwiseStats_enable)
      call mpas_pool_add_config(configPool, 'config_AM_pointwiseStats_compute_interval', config_AM_pointwiseStats_compute_interval)
      call mpas_pool_add_config(configPool, 'config_AM_pointwiseStats_output_stream', config_AM_pointwiseStats_output_stream)
      call mpas_pool_add_config(configPool, 'config_AM_pointwiseStats_compute_on_startup', &
config_AM_pointwiseStats_compute_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_pointwiseStats_write_on_startup', config_AM_pointwiseStats_write_on_startup)

   end subroutine ocn_setup_nmlrec_AM_pointwiseStats


   subroutine ocn_setup_nmlrec_AM_debugDiagnostics(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_AM_debugDiagnostics_enable = .false.
      character (len=StrKIND) :: config_AM_debugDiagnostics_compute_interval = 'output_interval'
      character (len=StrKIND) :: config_AM_debugDiagnostics_output_stream = 'debugDiagnosticsOutput'
      logical :: config_AM_debugDiagnostics_compute_on_startup = .true.
      logical :: config_AM_debugDiagnostics_write_on_startup = .true.
      logical :: config_AM_debugDiagnostics_check_state = .false.

      namelist /AM_debugDiagnostics/ &
         config_AM_debugDiagnostics_enable, &
         config_AM_debugDiagnostics_compute_interval, &
         config_AM_debugDiagnostics_output_stream, &
         config_AM_debugDiagnostics_compute_on_startup, &
         config_AM_debugDiagnostics_write_on_startup, &
         config_AM_debugDiagnostics_check_state
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, AM_debugDiagnostics, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_AM_debugDiagnostics_enable)
         call mpas_dmpar_bcast_char(dminfo, config_AM_debugDiagnostics_compute_interval)
         call mpas_dmpar_bcast_char(dminfo, config_AM_debugDiagnostics_output_stream)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_debugDiagnostics_compute_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_debugDiagnostics_write_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_debugDiagnostics_check_state)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record AM_debugDiagnostics')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_AM_debugDiagnostics_enable = $l', logicArgs=(/config_AM_debugDiagnostics_enable/))
            call mpas_log_write('        config_AM_debugDiagnostics_compute_interval = '&
//''//mpas_log_escape_dollars(config_AM_debugDiagnostics_compute_interval))
            call mpas_log_write('        config_AM_debugDiagnostics_output_stream = '&
//''//mpas_log_escape_dollars(config_AM_debugDiagnostics_output_stream))
            call mpas_log_write('        config_AM_debugDiagnostics_compute_on_startup = $l', &
logicArgs=(/config_AM_debugDiagnostics_compute_on_startup/))
            call mpas_log_write('        config_AM_debugDiagnostics_write_on_startup = $l', &
logicArgs=(/config_AM_debugDiagnostics_write_on_startup/))
            call mpas_log_write('        config_AM_debugDiagnostics_check_state = $l', &
logicArgs=(/config_AM_debugDiagnostics_check_state/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record AM_debugDiagnostics.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_AM_debugDiagnostics_enable', config_AM_debugDiagnostics_enable)
      call mpas_pool_add_config(configPool, 'config_AM_debugDiagnostics_compute_interval', &
config_AM_debugDiagnostics_compute_interval)
      call mpas_pool_add_config(configPool, 'config_AM_debugDiagnostics_output_stream', config_AM_debugDiagnostics_output_stream)
      call mpas_pool_add_config(configPool, 'config_AM_debugDiagnostics_compute_on_startup', &
config_AM_debugDiagnostics_compute_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_debugDiagnostics_write_on_startup', &
config_AM_debugDiagnostics_write_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_debugDiagnostics_check_state', config_AM_debugDiagnostics_check_state)

   end subroutine ocn_setup_nmlrec_AM_debugDiagnostics


   subroutine ocn_setup_nmlrec_AM_rpnCalculator(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_AM_rpnCalculator_enable = .false.
      logical :: config_AM_rpnCalculator_compute_on_startup = .true.
      logical :: config_AM_rpnCalculator_write_on_startup = .false.
      character (len=StrKIND) :: config_AM_rpnCalculator_compute_interval = '0010-00-00_00:00:00'
      character (len=StrKIND) :: config_AM_rpnCalculator_output_stream = 'none'
      character (len=StrKIND) :: config_AM_rpnCalculator_variable_a = 'layerThickness'
      character (len=StrKIND) :: config_AM_rpnCalculator_variable_b = 'areaCell'
      character (len=StrKIND) :: config_AM_rpnCalculator_variable_c = 'none'
      character (len=StrKIND) :: config_AM_rpnCalculator_variable_d = 'none'
      character (len=StrKIND) :: config_AM_rpnCalculator_variable_e = 'none'
      character (len=StrKIND) :: config_AM_rpnCalculator_variable_f = 'none'
      character (len=StrKIND) :: config_AM_rpnCalculator_variable_g = 'none'
      character (len=StrKIND) :: config_AM_rpnCalculator_variable_h = 'none'
      character (len=StrKIND) :: config_AM_rpnCalculator_expression_1 = 'a b *'
      character (len=StrKIND) :: config_AM_rpnCalculator_expression_2 = 'none'
      character (len=StrKIND) :: config_AM_rpnCalculator_expression_3 = 'none'
      character (len=StrKIND) :: config_AM_rpnCalculator_expression_4 = 'none'
      character (len=StrKIND) :: config_AM_rpnCalculator_output_name_1 = 'volumeCell'
      character (len=StrKIND) :: config_AM_rpnCalculator_output_name_2 = 'none'
      character (len=StrKIND) :: config_AM_rpnCalculator_output_name_3 = 'none'
      character (len=StrKIND) :: config_AM_rpnCalculator_output_name_4 = 'none'

      namelist /AM_rpnCalculator/ &
         config_AM_rpnCalculator_enable, &
         config_AM_rpnCalculator_compute_on_startup, &
         config_AM_rpnCalculator_write_on_startup, &
         config_AM_rpnCalculator_compute_interval, &
         config_AM_rpnCalculator_output_stream, &
         config_AM_rpnCalculator_variable_a, &
         config_AM_rpnCalculator_variable_b, &
         config_AM_rpnCalculator_variable_c, &
         config_AM_rpnCalculator_variable_d, &
         config_AM_rpnCalculator_variable_e, &
         config_AM_rpnCalculator_variable_f, &
         config_AM_rpnCalculator_variable_g, &
         config_AM_rpnCalculator_variable_h, &
         config_AM_rpnCalculator_expression_1, &
         config_AM_rpnCalculator_expression_2, &
         config_AM_rpnCalculator_expression_3, &
         config_AM_rpnCalculator_expression_4, &
         config_AM_rpnCalculator_output_name_1, &
         config_AM_rpnCalculator_output_name_2, &
         config_AM_rpnCalculator_output_name_3, &
         config_AM_rpnCalculator_output_name_4
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, AM_rpnCalculator, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_AM_rpnCalculator_enable)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_rpnCalculator_compute_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_rpnCalculator_write_on_startup)
         call mpas_dmpar_bcast_char(dminfo, config_AM_rpnCalculator_compute_interval)
         call mpas_dmpar_bcast_char(dminfo, config_AM_rpnCalculator_output_stream)
         call mpas_dmpar_bcast_char(dminfo, config_AM_rpnCalculator_variable_a)
         call mpas_dmpar_bcast_char(dminfo, config_AM_rpnCalculator_variable_b)
         call mpas_dmpar_bcast_char(dminfo, config_AM_rpnCalculator_variable_c)
         call mpas_dmpar_bcast_char(dminfo, config_AM_rpnCalculator_variable_d)
         call mpas_dmpar_bcast_char(dminfo, config_AM_rpnCalculator_variable_e)
         call mpas_dmpar_bcast_char(dminfo, config_AM_rpnCalculator_variable_f)
         call mpas_dmpar_bcast_char(dminfo, config_AM_rpnCalculator_variable_g)
         call mpas_dmpar_bcast_char(dminfo, config_AM_rpnCalculator_variable_h)
         call mpas_dmpar_bcast_char(dminfo, config_AM_rpnCalculator_expression_1)
         call mpas_dmpar_bcast_char(dminfo, config_AM_rpnCalculator_expression_2)
         call mpas_dmpar_bcast_char(dminfo, config_AM_rpnCalculator_expression_3)
         call mpas_dmpar_bcast_char(dminfo, config_AM_rpnCalculator_expression_4)
         call mpas_dmpar_bcast_char(dminfo, config_AM_rpnCalculator_output_name_1)
         call mpas_dmpar_bcast_char(dminfo, config_AM_rpnCalculator_output_name_2)
         call mpas_dmpar_bcast_char(dminfo, config_AM_rpnCalculator_output_name_3)
         call mpas_dmpar_bcast_char(dminfo, config_AM_rpnCalculator_output_name_4)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record AM_rpnCalculator')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_AM_rpnCalculator_enable = $l', logicArgs=(/config_AM_rpnCalculator_enable/))
            call mpas_log_write('        config_AM_rpnCalculator_compute_on_startup = $l', &
logicArgs=(/config_AM_rpnCalculator_compute_on_startup/))
            call mpas_log_write('        config_AM_rpnCalculator_write_on_startup = $l', &
logicArgs=(/config_AM_rpnCalculator_write_on_startup/))
            call mpas_log_write('        config_AM_rpnCalculator_compute_interval = '&
//''//mpas_log_escape_dollars(config_AM_rpnCalculator_compute_interval))
            call mpas_log_write('        config_AM_rpnCalculator_output_stream = '&
//''//mpas_log_escape_dollars(config_AM_rpnCalculator_output_stream))
            call mpas_log_write('        config_AM_rpnCalculator_variable_a = '&
//''//mpas_log_escape_dollars(config_AM_rpnCalculator_variable_a))
            call mpas_log_write('        config_AM_rpnCalculator_variable_b = '&
//''//mpas_log_escape_dollars(config_AM_rpnCalculator_variable_b))
            call mpas_log_write('        config_AM_rpnCalculator_variable_c = '&
//''//mpas_log_escape_dollars(config_AM_rpnCalculator_variable_c))
            call mpas_log_write('        config_AM_rpnCalculator_variable_d = '&
//''//mpas_log_escape_dollars(config_AM_rpnCalculator_variable_d))
            call mpas_log_write('        config_AM_rpnCalculator_variable_e = '&
//''//mpas_log_escape_dollars(config_AM_rpnCalculator_variable_e))
            call mpas_log_write('        config_AM_rpnCalculator_variable_f = '&
//''//mpas_log_escape_dollars(config_AM_rpnCalculator_variable_f))
            call mpas_log_write('        config_AM_rpnCalculator_variable_g = '&
//''//mpas_log_escape_dollars(config_AM_rpnCalculator_variable_g))
            call mpas_log_write('        config_AM_rpnCalculator_variable_h = '&
//''//mpas_log_escape_dollars(config_AM_rpnCalculator_variable_h))
            call mpas_log_write('        config_AM_rpnCalculator_expression_1 = '&
//''//mpas_log_escape_dollars(config_AM_rpnCalculator_expression_1))
            call mpas_log_write('        config_AM_rpnCalculator_expression_2 = '&
//''//mpas_log_escape_dollars(config_AM_rpnCalculator_expression_2))
            call mpas_log_write('        config_AM_rpnCalculator_expression_3 = '&
//''//mpas_log_escape_dollars(config_AM_rpnCalculator_expression_3))
            call mpas_log_write('        config_AM_rpnCalculator_expression_4 = '&
//''//mpas_log_escape_dollars(config_AM_rpnCalculator_expression_4))
            call mpas_log_write('        config_AM_rpnCalculator_output_name_1 = '&
//''//mpas_log_escape_dollars(config_AM_rpnCalculator_output_name_1))
            call mpas_log_write('        config_AM_rpnCalculator_output_name_2 = '&
//''//mpas_log_escape_dollars(config_AM_rpnCalculator_output_name_2))
            call mpas_log_write('        config_AM_rpnCalculator_output_name_3 = '&
//''//mpas_log_escape_dollars(config_AM_rpnCalculator_output_name_3))
            call mpas_log_write('        config_AM_rpnCalculator_output_name_4 = '&
//''//mpas_log_escape_dollars(config_AM_rpnCalculator_output_name_4))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record AM_rpnCalculator.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_AM_rpnCalculator_enable', config_AM_rpnCalculator_enable)
      call mpas_pool_add_config(configPool, 'config_AM_rpnCalculator_compute_on_startup', &
config_AM_rpnCalculator_compute_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_rpnCalculator_write_on_startup', config_AM_rpnCalculator_write_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_rpnCalculator_compute_interval', config_AM_rpnCalculator_compute_interval)
      call mpas_pool_add_config(configPool, 'config_AM_rpnCalculator_output_stream', config_AM_rpnCalculator_output_stream)
      call mpas_pool_add_config(configPool, 'config_AM_rpnCalculator_variable_a', config_AM_rpnCalculator_variable_a)
      call mpas_pool_add_config(configPool, 'config_AM_rpnCalculator_variable_b', config_AM_rpnCalculator_variable_b)
      call mpas_pool_add_config(configPool, 'config_AM_rpnCalculator_variable_c', config_AM_rpnCalculator_variable_c)
      call mpas_pool_add_config(configPool, 'config_AM_rpnCalculator_variable_d', config_AM_rpnCalculator_variable_d)
      call mpas_pool_add_config(configPool, 'config_AM_rpnCalculator_variable_e', config_AM_rpnCalculator_variable_e)
      call mpas_pool_add_config(configPool, 'config_AM_rpnCalculator_variable_f', config_AM_rpnCalculator_variable_f)
      call mpas_pool_add_config(configPool, 'config_AM_rpnCalculator_variable_g', config_AM_rpnCalculator_variable_g)
      call mpas_pool_add_config(configPool, 'config_AM_rpnCalculator_variable_h', config_AM_rpnCalculator_variable_h)
      call mpas_pool_add_config(configPool, 'config_AM_rpnCalculator_expression_1', config_AM_rpnCalculator_expression_1)
      call mpas_pool_add_config(configPool, 'config_AM_rpnCalculator_expression_2', config_AM_rpnCalculator_expression_2)
      call mpas_pool_add_config(configPool, 'config_AM_rpnCalculator_expression_3', config_AM_rpnCalculator_expression_3)
      call mpas_pool_add_config(configPool, 'config_AM_rpnCalculator_expression_4', config_AM_rpnCalculator_expression_4)
      call mpas_pool_add_config(configPool, 'config_AM_rpnCalculator_output_name_1', config_AM_rpnCalculator_output_name_1)
      call mpas_pool_add_config(configPool, 'config_AM_rpnCalculator_output_name_2', config_AM_rpnCalculator_output_name_2)
      call mpas_pool_add_config(configPool, 'config_AM_rpnCalculator_output_name_3', config_AM_rpnCalculator_output_name_3)
      call mpas_pool_add_config(configPool, 'config_AM_rpnCalculator_output_name_4', config_AM_rpnCalculator_output_name_4)

   end subroutine ocn_setup_nmlrec_AM_rpnCalculator


   subroutine ocn_setup_nmlrec_AM_transectTransport(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_AM_transectTransport_enable = .false.
      character (len=StrKIND) :: config_AM_transectTransport_compute_interval = 'output_interval'
      character (len=StrKIND) :: config_AM_transectTransport_output_stream = 'transectTransportOutput'
      logical :: config_AM_transectTransport_compute_on_startup = .true.
      logical :: config_AM_transectTransport_write_on_startup = .true.
      character (len=StrKIND) :: config_AM_transectTransport_transect_group = 'all'

      namelist /AM_transectTransport/ &
         config_AM_transectTransport_enable, &
         config_AM_transectTransport_compute_interval, &
         config_AM_transectTransport_output_stream, &
         config_AM_transectTransport_compute_on_startup, &
         config_AM_transectTransport_write_on_startup, &
         config_AM_transectTransport_transect_group
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, AM_transectTransport, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_AM_transectTransport_enable)
         call mpas_dmpar_bcast_char(dminfo, config_AM_transectTransport_compute_interval)
         call mpas_dmpar_bcast_char(dminfo, config_AM_transectTransport_output_stream)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_transectTransport_compute_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_transectTransport_write_on_startup)
         call mpas_dmpar_bcast_char(dminfo, config_AM_transectTransport_transect_group)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record AM_transectTransport')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_AM_transectTransport_enable = $l', &
logicArgs=(/config_AM_transectTransport_enable/))
            call mpas_log_write('        config_AM_transectTransport_compute_interval = '&
//''//mpas_log_escape_dollars(config_AM_transectTransport_compute_interval))
            call mpas_log_write('        config_AM_transectTransport_output_stream = '&
//''//mpas_log_escape_dollars(config_AM_transectTransport_output_stream))
            call mpas_log_write('        config_AM_transectTransport_compute_on_startup = $l', &
logicArgs=(/config_AM_transectTransport_compute_on_startup/))
            call mpas_log_write('        config_AM_transectTransport_write_on_startup = $l', &
logicArgs=(/config_AM_transectTransport_write_on_startup/))
            call mpas_log_write('        config_AM_transectTransport_transect_group = '&
//''//mpas_log_escape_dollars(config_AM_transectTransport_transect_group))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record AM_transectTransport.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_AM_transectTransport_enable', config_AM_transectTransport_enable)
      call mpas_pool_add_config(configPool, 'config_AM_transectTransport_compute_interval', &
config_AM_transectTransport_compute_interval)
      call mpas_pool_add_config(configPool, 'config_AM_transectTransport_output_stream', config_AM_transectTransport_output_stream)
      call mpas_pool_add_config(configPool, 'config_AM_transectTransport_compute_on_startup', &
config_AM_transectTransport_compute_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_transectTransport_write_on_startup', &
config_AM_transectTransport_write_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_transectTransport_transect_group', &
config_AM_transectTransport_transect_group)

   end subroutine ocn_setup_nmlrec_AM_transectTransport


   subroutine ocn_setup_nmlrec_AM_eddyProductVariables(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_AM_eddyProductVariables_enable = .false.
      character (len=StrKIND) :: config_AM_eddyProductVariables_compute_interval = 'dt'
      character (len=StrKIND) :: config_AM_eddyProductVariables_output_stream = 'eddyProductVariablesOutput'
      logical :: config_AM_eddyProductVariables_compute_on_startup = .true.
      logical :: config_AM_eddyProductVariables_write_on_startup = .false.

      namelist /AM_eddyProductVariables/ &
         config_AM_eddyProductVariables_enable, &
         config_AM_eddyProductVariables_compute_interval, &
         config_AM_eddyProductVariables_output_stream, &
         config_AM_eddyProductVariables_compute_on_startup, &
         config_AM_eddyProductVariables_write_on_startup
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, AM_eddyProductVariables, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_AM_eddyProductVariables_enable)
         call mpas_dmpar_bcast_char(dminfo, config_AM_eddyProductVariables_compute_interval)
         call mpas_dmpar_bcast_char(dminfo, config_AM_eddyProductVariables_output_stream)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_eddyProductVariables_compute_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_eddyProductVariables_write_on_startup)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record AM_eddyProductVariables')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_AM_eddyProductVariables_enable = $l', &
logicArgs=(/config_AM_eddyProductVariables_enable/))
            call mpas_log_write('        config_AM_eddyProductVariables_compute_interval = '&
//''//mpas_log_escape_dollars(config_AM_eddyProductVariables_compute_interval))
            call mpas_log_write('        config_AM_eddyProductVariables_output_stream = '&
//''//mpas_log_escape_dollars(config_AM_eddyProductVariables_output_stream))
            call mpas_log_write('        config_AM_eddyProductVariables_compute_on_startup = $l', &
logicArgs=(/config_AM_eddyProductVariables_compute_on_startup/))
            call mpas_log_write('        config_AM_eddyProductVariables_write_on_startup = $l', &
logicArgs=(/config_AM_eddyProductVariables_write_on_startup/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record AM_eddyProductVariables.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_AM_eddyProductVariables_enable', config_AM_eddyProductVariables_enable)
      call mpas_pool_add_config(configPool, 'config_AM_eddyProductVariables_compute_interval', &
config_AM_eddyProductVariables_compute_interval)
      call mpas_pool_add_config(configPool, 'config_AM_eddyProductVariables_output_stream', &
config_AM_eddyProductVariables_output_stream)
      call mpas_pool_add_config(configPool, 'config_AM_eddyProductVariables_compute_on_startup', &
config_AM_eddyProductVariables_compute_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_eddyProductVariables_write_on_startup', &
config_AM_eddyProductVariables_write_on_startup)

   end subroutine ocn_setup_nmlrec_AM_eddyProductVariables


   subroutine ocn_setup_nmlrec_AM_mocStreamfunction(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_AM_mocStreamfunction_enable = .false.
      character (len=StrKIND) :: config_AM_mocStreamfunction_compute_interval = 'output_interval'
      character (len=StrKIND) :: config_AM_mocStreamfunction_output_stream = 'mocStreamfunctionOutput'
      logical :: config_AM_mocStreamfunction_compute_on_startup = .true.
      logical :: config_AM_mocStreamfunction_write_on_startup = .true.
      real (kind=RKIND) :: config_AM_mocStreamfunction_min_bin = -9999999999999999455752309870428160.000000
      real (kind=RKIND) :: config_AM_mocStreamfunction_max_bin = -9999999999999999455752309870428160.000000
      integer :: config_AM_mocStreamfunction_num_bins = 180
      character (len=StrKIND) :: config_AM_mocStreamfunction_region_group = 'all'
      character (len=StrKIND) :: config_AM_mocStreamfunction_transect_group = 'all'

      namelist /AM_mocStreamfunction/ &
         config_AM_mocStreamfunction_enable, &
         config_AM_mocStreamfunction_compute_interval, &
         config_AM_mocStreamfunction_output_stream, &
         config_AM_mocStreamfunction_compute_on_startup, &
         config_AM_mocStreamfunction_write_on_startup, &
         config_AM_mocStreamfunction_min_bin, &
         config_AM_mocStreamfunction_max_bin, &
         config_AM_mocStreamfunction_num_bins, &
         config_AM_mocStreamfunction_region_group, &
         config_AM_mocStreamfunction_transect_group
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, AM_mocStreamfunction, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_AM_mocStreamfunction_enable)
         call mpas_dmpar_bcast_char(dminfo, config_AM_mocStreamfunction_compute_interval)
         call mpas_dmpar_bcast_char(dminfo, config_AM_mocStreamfunction_output_stream)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_mocStreamfunction_compute_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_mocStreamfunction_write_on_startup)
         call mpas_dmpar_bcast_real(dminfo, config_AM_mocStreamfunction_min_bin)
         call mpas_dmpar_bcast_real(dminfo, config_AM_mocStreamfunction_max_bin)
         call mpas_dmpar_bcast_int(dminfo, config_AM_mocStreamfunction_num_bins)
         call mpas_dmpar_bcast_char(dminfo, config_AM_mocStreamfunction_region_group)
         call mpas_dmpar_bcast_char(dminfo, config_AM_mocStreamfunction_transect_group)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record AM_mocStreamfunction')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_AM_mocStreamfunction_enable = $l', &
logicArgs=(/config_AM_mocStreamfunction_enable/))
            call mpas_log_write('        config_AM_mocStreamfunction_compute_interval = '&
//''//mpas_log_escape_dollars(config_AM_mocStreamfunction_compute_interval))
            call mpas_log_write('        config_AM_mocStreamfunction_output_stream = '&
//''//mpas_log_escape_dollars(config_AM_mocStreamfunction_output_stream))
            call mpas_log_write('        config_AM_mocStreamfunction_compute_on_startup = $l', &
logicArgs=(/config_AM_mocStreamfunction_compute_on_startup/))
            call mpas_log_write('        config_AM_mocStreamfunction_write_on_startup = $l', &
logicArgs=(/config_AM_mocStreamfunction_write_on_startup/))
            call mpas_log_write('        config_AM_mocStreamfunction_min_bin = $r', &
realArgs=(/config_AM_mocStreamfunction_min_bin/))
            call mpas_log_write('        config_AM_mocStreamfunction_max_bin = $r', &
realArgs=(/config_AM_mocStreamfunction_max_bin/))
            call mpas_log_write('        config_AM_mocStreamfunction_num_bins = $i', &
intArgs=(/config_AM_mocStreamfunction_num_bins/))
            call mpas_log_write('        config_AM_mocStreamfunction_region_group = '&
//''//mpas_log_escape_dollars(config_AM_mocStreamfunction_region_group))
            call mpas_log_write('        config_AM_mocStreamfunction_transect_group = '&
//''//mpas_log_escape_dollars(config_AM_mocStreamfunction_transect_group))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record AM_mocStreamfunction.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_AM_mocStreamfunction_enable', config_AM_mocStreamfunction_enable)
      call mpas_pool_add_config(configPool, 'config_AM_mocStreamfunction_compute_interval', &
config_AM_mocStreamfunction_compute_interval)
      call mpas_pool_add_config(configPool, 'config_AM_mocStreamfunction_output_stream', config_AM_mocStreamfunction_output_stream)
      call mpas_pool_add_config(configPool, 'config_AM_mocStreamfunction_compute_on_startup', &
config_AM_mocStreamfunction_compute_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_mocStreamfunction_write_on_startup', &
config_AM_mocStreamfunction_write_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_mocStreamfunction_min_bin', config_AM_mocStreamfunction_min_bin)
      call mpas_pool_add_config(configPool, 'config_AM_mocStreamfunction_max_bin', config_AM_mocStreamfunction_max_bin)
      call mpas_pool_add_config(configPool, 'config_AM_mocStreamfunction_num_bins', config_AM_mocStreamfunction_num_bins)
      call mpas_pool_add_config(configPool, 'config_AM_mocStreamfunction_region_group', config_AM_mocStreamfunction_region_group)
      call mpas_pool_add_config(configPool, 'config_AM_mocStreamfunction_transect_group', &
config_AM_mocStreamfunction_transect_group)

   end subroutine ocn_setup_nmlrec_AM_mocStreamfunction


   subroutine ocn_setup_nmlrec_AM_oceanHeatContent(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_AM_oceanHeatContent_enable = .false.
      character (len=StrKIND) :: config_AM_oceanHeatContent_compute_interval = 'output_interval'
      character (len=StrKIND) :: config_AM_oceanHeatContent_output_stream = 'oceanHeatContentOutput'
      logical :: config_AM_oceanHeatContent_compute_on_startup = .true.
      logical :: config_AM_oceanHeatContent_write_on_startup = .true.

      namelist /AM_oceanHeatContent/ &
         config_AM_oceanHeatContent_enable, &
         config_AM_oceanHeatContent_compute_interval, &
         config_AM_oceanHeatContent_output_stream, &
         config_AM_oceanHeatContent_compute_on_startup, &
         config_AM_oceanHeatContent_write_on_startup
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, AM_oceanHeatContent, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_AM_oceanHeatContent_enable)
         call mpas_dmpar_bcast_char(dminfo, config_AM_oceanHeatContent_compute_interval)
         call mpas_dmpar_bcast_char(dminfo, config_AM_oceanHeatContent_output_stream)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_oceanHeatContent_compute_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_oceanHeatContent_write_on_startup)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record AM_oceanHeatContent')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_AM_oceanHeatContent_enable = $l', logicArgs=(/config_AM_oceanHeatContent_enable/))
            call mpas_log_write('        config_AM_oceanHeatContent_compute_interval = '&
//''//mpas_log_escape_dollars(config_AM_oceanHeatContent_compute_interval))
            call mpas_log_write('        config_AM_oceanHeatContent_output_stream = '&
//''//mpas_log_escape_dollars(config_AM_oceanHeatContent_output_stream))
            call mpas_log_write('        config_AM_oceanHeatContent_compute_on_startup = $l', &
logicArgs=(/config_AM_oceanHeatContent_compute_on_startup/))
            call mpas_log_write('        config_AM_oceanHeatContent_write_on_startup = $l', &
logicArgs=(/config_AM_oceanHeatContent_write_on_startup/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record AM_oceanHeatContent.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_AM_oceanHeatContent_enable', config_AM_oceanHeatContent_enable)
      call mpas_pool_add_config(configPool, 'config_AM_oceanHeatContent_compute_interval', &
config_AM_oceanHeatContent_compute_interval)
      call mpas_pool_add_config(configPool, 'config_AM_oceanHeatContent_output_stream', config_AM_oceanHeatContent_output_stream)
      call mpas_pool_add_config(configPool, 'config_AM_oceanHeatContent_compute_on_startup', &
config_AM_oceanHeatContent_compute_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_oceanHeatContent_write_on_startup', &
config_AM_oceanHeatContent_write_on_startup)

   end subroutine ocn_setup_nmlrec_AM_oceanHeatContent


   subroutine ocn_setup_nmlrec_AM_mixedLayerHeatBudget(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_AM_mixedLayerHeatBudget_enable = .false.
      character (len=StrKIND) :: config_AM_mixedLayerHeatBudget_compute_interval = 'output_interval'
      character (len=StrKIND) :: config_AM_mixedLayerHeatBudget_output_stream = 'mixedLayerHeatBudgetOutput'
      logical :: config_AM_mixedLayerHeatBudget_compute_on_startup = .true.
      logical :: config_AM_mixedLayerHeatBudget_write_on_startup = .true.

      namelist /AM_mixedLayerHeatBudget/ &
         config_AM_mixedLayerHeatBudget_enable, &
         config_AM_mixedLayerHeatBudget_compute_interval, &
         config_AM_mixedLayerHeatBudget_output_stream, &
         config_AM_mixedLayerHeatBudget_compute_on_startup, &
         config_AM_mixedLayerHeatBudget_write_on_startup
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, AM_mixedLayerHeatBudget, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_AM_mixedLayerHeatBudget_enable)
         call mpas_dmpar_bcast_char(dminfo, config_AM_mixedLayerHeatBudget_compute_interval)
         call mpas_dmpar_bcast_char(dminfo, config_AM_mixedLayerHeatBudget_output_stream)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_mixedLayerHeatBudget_compute_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_mixedLayerHeatBudget_write_on_startup)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record AM_mixedLayerHeatBudget')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_AM_mixedLayerHeatBudget_enable = $l', &
logicArgs=(/config_AM_mixedLayerHeatBudget_enable/))
            call mpas_log_write('        config_AM_mixedLayerHeatBudget_compute_interval = '&
//''//mpas_log_escape_dollars(config_AM_mixedLayerHeatBudget_compute_interval))
            call mpas_log_write('        config_AM_mixedLayerHeatBudget_output_stream = '&
//''//mpas_log_escape_dollars(config_AM_mixedLayerHeatBudget_output_stream))
            call mpas_log_write('        config_AM_mixedLayerHeatBudget_compute_on_startup = $l', &
logicArgs=(/config_AM_mixedLayerHeatBudget_compute_on_startup/))
            call mpas_log_write('        config_AM_mixedLayerHeatBudget_write_on_startup = $l', &
logicArgs=(/config_AM_mixedLayerHeatBudget_write_on_startup/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record AM_mixedLayerHeatBudget.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_AM_mixedLayerHeatBudget_enable', config_AM_mixedLayerHeatBudget_enable)
      call mpas_pool_add_config(configPool, 'config_AM_mixedLayerHeatBudget_compute_interval', &
config_AM_mixedLayerHeatBudget_compute_interval)
      call mpas_pool_add_config(configPool, 'config_AM_mixedLayerHeatBudget_output_stream', &
config_AM_mixedLayerHeatBudget_output_stream)
      call mpas_pool_add_config(configPool, 'config_AM_mixedLayerHeatBudget_compute_on_startup', &
config_AM_mixedLayerHeatBudget_compute_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_mixedLayerHeatBudget_write_on_startup', &
config_AM_mixedLayerHeatBudget_write_on_startup)

   end subroutine ocn_setup_nmlrec_AM_mixedLayerHeatBudget


   subroutine ocn_setup_nmlrec_AM_sedimentFluxIndex(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_AM_sedimentFluxIndex_enable = .false.
      logical :: config_AM_sedimentFluxIndex_compute_on_startup = .true.
      logical :: config_AM_sedimentFluxIndex_write_on_startup = .true.
      character (len=StrKIND) :: config_AM_sedimentFluxIndex_compute_interval = 'output_interval'
      character (len=StrKIND) :: config_AM_sedimentFluxIndex_output_stream = 'sedimentFluxIndexOutput'
      character (len=StrKIND) :: config_AM_sedimentFluxIndex_directory = 'analysis_members'
      logical :: config_AM_sedimentFluxIndex_use_lat_lon_coords = .true.

      namelist /AM_sedimentFluxIndex/ &
         config_AM_sedimentFluxIndex_enable, &
         config_AM_sedimentFluxIndex_compute_on_startup, &
         config_AM_sedimentFluxIndex_write_on_startup, &
         config_AM_sedimentFluxIndex_compute_interval, &
         config_AM_sedimentFluxIndex_output_stream, &
         config_AM_sedimentFluxIndex_directory, &
         config_AM_sedimentFluxIndex_use_lat_lon_coords
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, AM_sedimentFluxIndex, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_AM_sedimentFluxIndex_enable)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_sedimentFluxIndex_compute_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_sedimentFluxIndex_write_on_startup)
         call mpas_dmpar_bcast_char(dminfo, config_AM_sedimentFluxIndex_compute_interval)
         call mpas_dmpar_bcast_char(dminfo, config_AM_sedimentFluxIndex_output_stream)
         call mpas_dmpar_bcast_char(dminfo, config_AM_sedimentFluxIndex_directory)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_sedimentFluxIndex_use_lat_lon_coords)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record AM_sedimentFluxIndex')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_AM_sedimentFluxIndex_enable = $l', &
logicArgs=(/config_AM_sedimentFluxIndex_enable/))
            call mpas_log_write('        config_AM_sedimentFluxIndex_compute_on_startup = $l', &
logicArgs=(/config_AM_sedimentFluxIndex_compute_on_startup/))
            call mpas_log_write('        config_AM_sedimentFluxIndex_write_on_startup = $l', &
logicArgs=(/config_AM_sedimentFluxIndex_write_on_startup/))
            call mpas_log_write('        config_AM_sedimentFluxIndex_compute_interval = '&
//''//mpas_log_escape_dollars(config_AM_sedimentFluxIndex_compute_interval))
            call mpas_log_write('        config_AM_sedimentFluxIndex_output_stream = '&
//''//mpas_log_escape_dollars(config_AM_sedimentFluxIndex_output_stream))
            call mpas_log_write('        config_AM_sedimentFluxIndex_directory = '&
//''//mpas_log_escape_dollars(config_AM_sedimentFluxIndex_directory))
            call mpas_log_write('        config_AM_sedimentFluxIndex_use_lat_lon_coords = $l', &
logicArgs=(/config_AM_sedimentFluxIndex_use_lat_lon_coords/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record AM_sedimentFluxIndex.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_AM_sedimentFluxIndex_enable', config_AM_sedimentFluxIndex_enable)
      call mpas_pool_add_config(configPool, 'config_AM_sedimentFluxIndex_compute_on_startup', &
config_AM_sedimentFluxIndex_compute_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_sedimentFluxIndex_write_on_startup', &
config_AM_sedimentFluxIndex_write_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_sedimentFluxIndex_compute_interval', &
config_AM_sedimentFluxIndex_compute_interval)
      call mpas_pool_add_config(configPool, 'config_AM_sedimentFluxIndex_output_stream', config_AM_sedimentFluxIndex_output_stream)
      call mpas_pool_add_config(configPool, 'config_AM_sedimentFluxIndex_directory', config_AM_sedimentFluxIndex_directory)
      call mpas_pool_add_config(configPool, 'config_AM_sedimentFluxIndex_use_lat_lon_coords', &
config_AM_sedimentFluxIndex_use_lat_lon_coords)

   end subroutine ocn_setup_nmlrec_AM_sedimentFluxIndex


   subroutine ocn_setup_nmlrec_AM_sedimentTransport(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_AM_sedimentTransport_enable = .false.
      logical :: config_AM_sedimentTransport_compute_on_startup = .true.
      logical :: config_AM_sedimentTransport_write_on_startup = .true.
      character (len=StrKIND) :: config_AM_sedimentTransport_compute_interval = 'output_interval'
      character (len=StrKIND) :: config_AM_sedimentTransport_output_stream = 'sedimentTransportOutput'
      character (len=StrKIND) :: config_AM_sedimentTransport_directory = 'analysis_members'
      real (kind=RKIND) :: config_AM_sedimentTransport_grain_size = 0.000250
      character (len=StrKIND) :: config_AM_sedimentTransport_ws_formula = 'VanRijn1993'
      character (len=StrKIND) :: config_AM_sedimentTransport_bedld_formula = 'Soulsby-Damgaard'
      character (len=StrKIND) :: config_AM_sedimentTransport_SSC_ref_formula = 'Lee2004'
      real (kind=RKIND) :: config_AM_sedimentTransport_drag_coefficient = 0.002500
      real (kind=RKIND) :: config_AM_sedimentTransport_erate = 0.000500
      real (kind=RKIND) :: config_AM_sedimentTransport_tau_ce = 0.100000
      real (kind=RKIND) :: config_AM_sedimentTransport_tau_cd = 0.100000
      real (kind=RKIND) :: config_AM_sedimentTransport_Manning_coef = 0.022000
      real (kind=RKIND) :: config_AM_sedimentTransport_grain_porosity = 0.500000
      real (kind=RKIND) :: config_AM_sedimentTransport_water_density = 1020.000000
      real (kind=RKIND) :: config_AM_sedimentTransport_grain_density = 2650.000000
      real (kind=RKIND) :: config_AM_sedimentTransport_alpha = 0.010000
      real (kind=RKIND) :: config_AM_sedimentTransport_kinematic_viscosity = 0.000001
      real (kind=RKIND) :: config_AM_sedimentTransport_vertical_diffusion_coefficient = 0.010000
      logical :: config_AM_sedimentTransport_bedload = .true.
      logical :: config_AM_sedimentTransport_suspended = .true.
      logical :: config_AM_sedimentTransport_use_lat_lon_coords = .true.

      namelist /AM_sedimentTransport/ &
         config_AM_sedimentTransport_enable, &
         config_AM_sedimentTransport_compute_on_startup, &
         config_AM_sedimentTransport_write_on_startup, &
         config_AM_sedimentTransport_compute_interval, &
         config_AM_sedimentTransport_output_stream, &
         config_AM_sedimentTransport_directory, &
         config_AM_sedimentTransport_grain_size, &
         config_AM_sedimentTransport_ws_formula, &
         config_AM_sedimentTransport_bedld_formula, &
         config_AM_sedimentTransport_SSC_ref_formula, &
         config_AM_sedimentTransport_drag_coefficient, &
         config_AM_sedimentTransport_erate, &
         config_AM_sedimentTransport_tau_ce, &
         config_AM_sedimentTransport_tau_cd, &
         config_AM_sedimentTransport_Manning_coef, &
         config_AM_sedimentTransport_grain_porosity, &
         config_AM_sedimentTransport_water_density, &
         config_AM_sedimentTransport_grain_density, &
         config_AM_sedimentTransport_alpha, &
         config_AM_sedimentTransport_kinematic_viscosity, &
         config_AM_sedimentTransport_vertical_diffusion_coefficient, &
         config_AM_sedimentTransport_bedload, &
         config_AM_sedimentTransport_suspended, &
         config_AM_sedimentTransport_use_lat_lon_coords
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, AM_sedimentTransport, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_AM_sedimentTransport_enable)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_sedimentTransport_compute_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_sedimentTransport_write_on_startup)
         call mpas_dmpar_bcast_char(dminfo, config_AM_sedimentTransport_compute_interval)
         call mpas_dmpar_bcast_char(dminfo, config_AM_sedimentTransport_output_stream)
         call mpas_dmpar_bcast_char(dminfo, config_AM_sedimentTransport_directory)
         call mpas_dmpar_bcast_real(dminfo, config_AM_sedimentTransport_grain_size)
         call mpas_dmpar_bcast_char(dminfo, config_AM_sedimentTransport_ws_formula)
         call mpas_dmpar_bcast_char(dminfo, config_AM_sedimentTransport_bedld_formula)
         call mpas_dmpar_bcast_char(dminfo, config_AM_sedimentTransport_SSC_ref_formula)
         call mpas_dmpar_bcast_real(dminfo, config_AM_sedimentTransport_drag_coefficient)
         call mpas_dmpar_bcast_real(dminfo, config_AM_sedimentTransport_erate)
         call mpas_dmpar_bcast_real(dminfo, config_AM_sedimentTransport_tau_ce)
         call mpas_dmpar_bcast_real(dminfo, config_AM_sedimentTransport_tau_cd)
         call mpas_dmpar_bcast_real(dminfo, config_AM_sedimentTransport_Manning_coef)
         call mpas_dmpar_bcast_real(dminfo, config_AM_sedimentTransport_grain_porosity)
         call mpas_dmpar_bcast_real(dminfo, config_AM_sedimentTransport_water_density)
         call mpas_dmpar_bcast_real(dminfo, config_AM_sedimentTransport_grain_density)
         call mpas_dmpar_bcast_real(dminfo, config_AM_sedimentTransport_alpha)
         call mpas_dmpar_bcast_real(dminfo, config_AM_sedimentTransport_kinematic_viscosity)
         call mpas_dmpar_bcast_real(dminfo, config_AM_sedimentTransport_vertical_diffusion_coefficient)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_sedimentTransport_bedload)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_sedimentTransport_suspended)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_sedimentTransport_use_lat_lon_coords)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record AM_sedimentTransport')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_AM_sedimentTransport_enable = $l', &
logicArgs=(/config_AM_sedimentTransport_enable/))
            call mpas_log_write('        config_AM_sedimentTransport_compute_on_startup = $l', &
logicArgs=(/config_AM_sedimentTransport_compute_on_startup/))
            call mpas_log_write('        config_AM_sedimentTransport_write_on_startup = $l', &
logicArgs=(/config_AM_sedimentTransport_write_on_startup/))
            call mpas_log_write('        config_AM_sedimentTransport_compute_interval = '&
//''//mpas_log_escape_dollars(config_AM_sedimentTransport_compute_interval))
            call mpas_log_write('        config_AM_sedimentTransport_output_stream = '&
//''//mpas_log_escape_dollars(config_AM_sedimentTransport_output_stream))
            call mpas_log_write('        config_AM_sedimentTransport_directory = '&
//''//mpas_log_escape_dollars(config_AM_sedimentTransport_directory))
            call mpas_log_write('        config_AM_sedimentTransport_grain_size = $r', &
realArgs=(/config_AM_sedimentTransport_grain_size/))
            call mpas_log_write('        config_AM_sedimentTransport_ws_formula = '&
//''//mpas_log_escape_dollars(config_AM_sedimentTransport_ws_formula))
            call mpas_log_write('        config_AM_sedimentTransport_bedld_formula = '&
//''//mpas_log_escape_dollars(config_AM_sedimentTransport_bedld_formula))
            call mpas_log_write('        config_AM_sedimentTransport_SSC_ref_formula = '&
//''//mpas_log_escape_dollars(config_AM_sedimentTransport_SSC_ref_formula))
            call mpas_log_write('        config_AM_sedimentTransport_drag_coefficient = $r', &
realArgs=(/config_AM_sedimentTransport_drag_coefficient/))
            call mpas_log_write('        config_AM_sedimentTransport_erate = $r', realArgs=(/config_AM_sedimentTransport_erate/))
            call mpas_log_write('        config_AM_sedimentTransport_tau_ce = $r', realArgs=(/config_AM_sedimentTransport_tau_ce/))
            call mpas_log_write('        config_AM_sedimentTransport_tau_cd = $r', realArgs=(/config_AM_sedimentTransport_tau_cd/))
            call mpas_log_write('        config_AM_sedimentTransport_Manning_coef = $r', &
realArgs=(/config_AM_sedimentTransport_Manning_coef/))
            call mpas_log_write('        config_AM_sedimentTransport_grain_porosity = $r', &
realArgs=(/config_AM_sedimentTransport_grain_porosity/))
            call mpas_log_write('        config_AM_sedimentTransport_water_density = $r', &
realArgs=(/config_AM_sedimentTransport_water_density/))
            call mpas_log_write('        config_AM_sedimentTransport_grain_density = $r', &
realArgs=(/config_AM_sedimentTransport_grain_density/))
            call mpas_log_write('        config_AM_sedimentTransport_alpha = $r', realArgs=(/config_AM_sedimentTransport_alpha/))
            call mpas_log_write('        config_AM_sedimentTransport_kinematic_viscosity = $r', &
realArgs=(/config_AM_sedimentTransport_kinematic_viscosity/))
            call mpas_log_write('        config_AM_sedimentTransport_vertical_diffusion_coefficient = $r', &
realArgs=(/config_AM_sedimentTransport_vertical_diffusion_coefficient/))
            call mpas_log_write('        config_AM_sedimentTransport_bedload = $l', &
logicArgs=(/config_AM_sedimentTransport_bedload/))
            call mpas_log_write('        config_AM_sedimentTransport_suspended = $l', &
logicArgs=(/config_AM_sedimentTransport_suspended/))
            call mpas_log_write('        config_AM_sedimentTransport_use_lat_lon_coords = $l', &
logicArgs=(/config_AM_sedimentTransport_use_lat_lon_coords/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record AM_sedimentTransport.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_AM_sedimentTransport_enable', config_AM_sedimentTransport_enable)
      call mpas_pool_add_config(configPool, 'config_AM_sedimentTransport_compute_on_startup', &
config_AM_sedimentTransport_compute_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_sedimentTransport_write_on_startup', &
config_AM_sedimentTransport_write_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_sedimentTransport_compute_interval', &
config_AM_sedimentTransport_compute_interval)
      call mpas_pool_add_config(configPool, 'config_AM_sedimentTransport_output_stream', config_AM_sedimentTransport_output_stream)
      call mpas_pool_add_config(configPool, 'config_AM_sedimentTransport_directory', config_AM_sedimentTransport_directory)
      call mpas_pool_add_config(configPool, 'config_AM_sedimentTransport_grain_size', config_AM_sedimentTransport_grain_size)
      call mpas_pool_add_config(configPool, 'config_AM_sedimentTransport_ws_formula', config_AM_sedimentTransport_ws_formula)
      call mpas_pool_add_config(configPool, 'config_AM_sedimentTransport_bedld_formula', config_AM_sedimentTransport_bedld_formula)
      call mpas_pool_add_config(configPool, 'config_AM_sedimentTransport_SSC_ref_formula', &
config_AM_sedimentTransport_SSC_ref_formula)
      call mpas_pool_add_config(configPool, 'config_AM_sedimentTransport_drag_coefficient', &
config_AM_sedimentTransport_drag_coefficient)
      call mpas_pool_add_config(configPool, 'config_AM_sedimentTransport_erate', config_AM_sedimentTransport_erate)
      call mpas_pool_add_config(configPool, 'config_AM_sedimentTransport_tau_ce', config_AM_sedimentTransport_tau_ce)
      call mpas_pool_add_config(configPool, 'config_AM_sedimentTransport_tau_cd', config_AM_sedimentTransport_tau_cd)
      call mpas_pool_add_config(configPool, 'config_AM_sedimentTransport_Manning_coef', config_AM_sedimentTransport_Manning_coef)
      call mpas_pool_add_config(configPool, 'config_AM_sedimentTransport_grain_porosity', &
config_AM_sedimentTransport_grain_porosity)
      call mpas_pool_add_config(configPool, 'config_AM_sedimentTransport_water_density', config_AM_sedimentTransport_water_density)
      call mpas_pool_add_config(configPool, 'config_AM_sedimentTransport_grain_density', config_AM_sedimentTransport_grain_density)
      call mpas_pool_add_config(configPool, 'config_AM_sedimentTransport_alpha', config_AM_sedimentTransport_alpha)
      call mpas_pool_add_config(configPool, 'config_AM_sedimentTransport_kinematic_viscosity', &
config_AM_sedimentTransport_kinematic_viscosity)
      call mpas_pool_add_config(configPool, 'config_AM_sedimentTransport_vertical_diffusion_coefficient', &
config_AM_sedimentTransport_vertical_diffusion_coefficient)
      call mpas_pool_add_config(configPool, 'config_AM_sedimentTransport_bedload', config_AM_sedimentTransport_bedload)
      call mpas_pool_add_config(configPool, 'config_AM_sedimentTransport_suspended', config_AM_sedimentTransport_suspended)
      call mpas_pool_add_config(configPool, 'config_AM_sedimentTransport_use_lat_lon_coords', &
config_AM_sedimentTransport_use_lat_lon_coords)

   end subroutine ocn_setup_nmlrec_AM_sedimentTransport


   subroutine ocn_setup_nmlrec_AM_harmonicAnalysis(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_AM_harmonicAnalysis_enable = .false.
      character (len=StrKIND) :: config_AM_harmonicAnalysis_compute_interval = 'output_interval'
      real (kind=RKIND) :: config_AM_harmonicAnalysis_start_delay = 20.000000
      real (kind=RKIND) :: config_AM_harmonicAnalysis_duration = 90.000000
      character (len=StrKIND) :: config_AM_harmonicAnalysis_output_stream = 'harmonicAnalysisOutput'
      character (len=StrKIND) :: config_AM_harmonicAnalysis_restart_stream = 'harmonicAnalysisRestart'
      logical :: config_AM_harmonicAnalysis_compute_on_startup = .false.
      logical :: config_AM_harmonicAnalysis_write_on_startup = .false.
      logical :: config_AM_harmonicAnalysis_use_M2 = .true.
      logical :: config_AM_harmonicAnalysis_use_S2 = .true.
      logical :: config_AM_harmonicAnalysis_use_N2 = .true.
      logical :: config_AM_harmonicAnalysis_use_K2 = .true.
      logical :: config_AM_harmonicAnalysis_use_K1 = .true.
      logical :: config_AM_harmonicAnalysis_use_O1 = .true.
      logical :: config_AM_harmonicAnalysis_use_Q1 = .true.
      logical :: config_AM_harmonicAnalysis_use_P1 = .true.

      namelist /AM_harmonicAnalysis/ &
         config_AM_harmonicAnalysis_enable, &
         config_AM_harmonicAnalysis_compute_interval, &
         config_AM_harmonicAnalysis_start_delay, &
         config_AM_harmonicAnalysis_duration, &
         config_AM_harmonicAnalysis_output_stream, &
         config_AM_harmonicAnalysis_restart_stream, &
         config_AM_harmonicAnalysis_compute_on_startup, &
         config_AM_harmonicAnalysis_write_on_startup, &
         config_AM_harmonicAnalysis_use_M2, &
         config_AM_harmonicAnalysis_use_S2, &
         config_AM_harmonicAnalysis_use_N2, &
         config_AM_harmonicAnalysis_use_K2, &
         config_AM_harmonicAnalysis_use_K1, &
         config_AM_harmonicAnalysis_use_O1, &
         config_AM_harmonicAnalysis_use_Q1, &
         config_AM_harmonicAnalysis_use_P1
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, AM_harmonicAnalysis, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_AM_harmonicAnalysis_enable)
         call mpas_dmpar_bcast_char(dminfo, config_AM_harmonicAnalysis_compute_interval)
         call mpas_dmpar_bcast_real(dminfo, config_AM_harmonicAnalysis_start_delay)
         call mpas_dmpar_bcast_real(dminfo, config_AM_harmonicAnalysis_duration)
         call mpas_dmpar_bcast_char(dminfo, config_AM_harmonicAnalysis_output_stream)
         call mpas_dmpar_bcast_char(dminfo, config_AM_harmonicAnalysis_restart_stream)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_harmonicAnalysis_compute_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_harmonicAnalysis_write_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_harmonicAnalysis_use_M2)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_harmonicAnalysis_use_S2)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_harmonicAnalysis_use_N2)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_harmonicAnalysis_use_K2)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_harmonicAnalysis_use_K1)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_harmonicAnalysis_use_O1)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_harmonicAnalysis_use_Q1)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_harmonicAnalysis_use_P1)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record AM_harmonicAnalysis')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_AM_harmonicAnalysis_enable = $l', logicArgs=(/config_AM_harmonicAnalysis_enable/))
            call mpas_log_write('        config_AM_harmonicAnalysis_compute_interval = '&
//''//mpas_log_escape_dollars(config_AM_harmonicAnalysis_compute_interval))
            call mpas_log_write('        config_AM_harmonicAnalysis_start_delay = $r', &
realArgs=(/config_AM_harmonicAnalysis_start_delay/))
            call mpas_log_write('        config_AM_harmonicAnalysis_duration = $r', &
realArgs=(/config_AM_harmonicAnalysis_duration/))
            call mpas_log_write('        config_AM_harmonicAnalysis_output_stream = '&
//''//mpas_log_escape_dollars(config_AM_harmonicAnalysis_output_stream))
            call mpas_log_write('        config_AM_harmonicAnalysis_restart_stream = '&
//''//mpas_log_escape_dollars(config_AM_harmonicAnalysis_restart_stream))
            call mpas_log_write('        config_AM_harmonicAnalysis_compute_on_startup = $l', &
logicArgs=(/config_AM_harmonicAnalysis_compute_on_startup/))
            call mpas_log_write('        config_AM_harmonicAnalysis_write_on_startup = $l', &
logicArgs=(/config_AM_harmonicAnalysis_write_on_startup/))
            call mpas_log_write('        config_AM_harmonicAnalysis_use_M2 = $l', logicArgs=(/config_AM_harmonicAnalysis_use_M2/))
            call mpas_log_write('        config_AM_harmonicAnalysis_use_S2 = $l', logicArgs=(/config_AM_harmonicAnalysis_use_S2/))
            call mpas_log_write('        config_AM_harmonicAnalysis_use_N2 = $l', logicArgs=(/config_AM_harmonicAnalysis_use_N2/))
            call mpas_log_write('        config_AM_harmonicAnalysis_use_K2 = $l', logicArgs=(/config_AM_harmonicAnalysis_use_K2/))
            call mpas_log_write('        config_AM_harmonicAnalysis_use_K1 = $l', logicArgs=(/config_AM_harmonicAnalysis_use_K1/))
            call mpas_log_write('        config_AM_harmonicAnalysis_use_O1 = $l', logicArgs=(/config_AM_harmonicAnalysis_use_O1/))
            call mpas_log_write('        config_AM_harmonicAnalysis_use_Q1 = $l', logicArgs=(/config_AM_harmonicAnalysis_use_Q1/))
            call mpas_log_write('        config_AM_harmonicAnalysis_use_P1 = $l', logicArgs=(/config_AM_harmonicAnalysis_use_P1/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record AM_harmonicAnalysis.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_AM_harmonicAnalysis_enable', config_AM_harmonicAnalysis_enable)
      call mpas_pool_add_config(configPool, 'config_AM_harmonicAnalysis_compute_interval', &
config_AM_harmonicAnalysis_compute_interval)
      call mpas_pool_add_config(configPool, 'config_AM_harmonicAnalysis_start_delay', config_AM_harmonicAnalysis_start_delay)
      call mpas_pool_add_config(configPool, 'config_AM_harmonicAnalysis_duration', config_AM_harmonicAnalysis_duration)
      call mpas_pool_add_config(configPool, 'config_AM_harmonicAnalysis_output_stream', config_AM_harmonicAnalysis_output_stream)
      call mpas_pool_add_config(configPool, 'config_AM_harmonicAnalysis_restart_stream', config_AM_harmonicAnalysis_restart_stream)
      call mpas_pool_add_config(configPool, 'config_AM_harmonicAnalysis_compute_on_startup', &
config_AM_harmonicAnalysis_compute_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_harmonicAnalysis_write_on_startup', &
config_AM_harmonicAnalysis_write_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_harmonicAnalysis_use_M2', config_AM_harmonicAnalysis_use_M2)
      call mpas_pool_add_config(configPool, 'config_AM_harmonicAnalysis_use_S2', config_AM_harmonicAnalysis_use_S2)
      call mpas_pool_add_config(configPool, 'config_AM_harmonicAnalysis_use_N2', config_AM_harmonicAnalysis_use_N2)
      call mpas_pool_add_config(configPool, 'config_AM_harmonicAnalysis_use_K2', config_AM_harmonicAnalysis_use_K2)
      call mpas_pool_add_config(configPool, 'config_AM_harmonicAnalysis_use_K1', config_AM_harmonicAnalysis_use_K1)
      call mpas_pool_add_config(configPool, 'config_AM_harmonicAnalysis_use_O1', config_AM_harmonicAnalysis_use_O1)
      call mpas_pool_add_config(configPool, 'config_AM_harmonicAnalysis_use_Q1', config_AM_harmonicAnalysis_use_Q1)
      call mpas_pool_add_config(configPool, 'config_AM_harmonicAnalysis_use_P1', config_AM_harmonicAnalysis_use_P1)

   end subroutine ocn_setup_nmlrec_AM_harmonicAnalysis


   subroutine ocn_setup_nmlrec_AM_conservationCheck(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_AM_conservationCheck_enable = .false.
      character (len=StrKIND) :: config_AM_conservationCheck_compute_interval = 'dt'
      character (len=StrKIND) :: config_AM_conservationCheck_output_stream = 'conservationCheckOutput'
      logical :: config_AM_conservationCheck_compute_on_startup = .false.
      logical :: config_AM_conservationCheck_write_on_startup = .false.
      logical :: config_AM_conservationCheck_write_to_logfile = .true.
      character (len=StrKIND) :: config_AM_conservationCheck_restart_stream = 'conservationCheckRestart'

      namelist /AM_conservationCheck/ &
         config_AM_conservationCheck_enable, &
         config_AM_conservationCheck_compute_interval, &
         config_AM_conservationCheck_output_stream, &
         config_AM_conservationCheck_compute_on_startup, &
         config_AM_conservationCheck_write_on_startup, &
         config_AM_conservationCheck_write_to_logfile, &
         config_AM_conservationCheck_restart_stream
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, AM_conservationCheck, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_AM_conservationCheck_enable)
         call mpas_dmpar_bcast_char(dminfo, config_AM_conservationCheck_compute_interval)
         call mpas_dmpar_bcast_char(dminfo, config_AM_conservationCheck_output_stream)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_conservationCheck_compute_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_conservationCheck_write_on_startup)
         call mpas_dmpar_bcast_logical(dminfo, config_AM_conservationCheck_write_to_logfile)
         call mpas_dmpar_bcast_char(dminfo, config_AM_conservationCheck_restart_stream)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record AM_conservationCheck')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_AM_conservationCheck_enable = $l', &
logicArgs=(/config_AM_conservationCheck_enable/))
            call mpas_log_write('        config_AM_conservationCheck_compute_interval = '&
//''//mpas_log_escape_dollars(config_AM_conservationCheck_compute_interval))
            call mpas_log_write('        config_AM_conservationCheck_output_stream = '&
//''//mpas_log_escape_dollars(config_AM_conservationCheck_output_stream))
            call mpas_log_write('        config_AM_conservationCheck_compute_on_startup = $l', &
logicArgs=(/config_AM_conservationCheck_compute_on_startup/))
            call mpas_log_write('        config_AM_conservationCheck_write_on_startup = $l', &
logicArgs=(/config_AM_conservationCheck_write_on_startup/))
            call mpas_log_write('        config_AM_conservationCheck_write_to_logfile = $l', &
logicArgs=(/config_AM_conservationCheck_write_to_logfile/))
            call mpas_log_write('        config_AM_conservationCheck_restart_stream = '&
//''//mpas_log_escape_dollars(config_AM_conservationCheck_restart_stream))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record AM_conservationCheck.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_AM_conservationCheck_enable', config_AM_conservationCheck_enable)
      call mpas_pool_add_config(configPool, 'config_AM_conservationCheck_compute_interval', &
config_AM_conservationCheck_compute_interval)
      call mpas_pool_add_config(configPool, 'config_AM_conservationCheck_output_stream', config_AM_conservationCheck_output_stream)
      call mpas_pool_add_config(configPool, 'config_AM_conservationCheck_compute_on_startup', &
config_AM_conservationCheck_compute_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_conservationCheck_write_on_startup', &
config_AM_conservationCheck_write_on_startup)
      call mpas_pool_add_config(configPool, 'config_AM_conservationCheck_write_to_logfile', &
config_AM_conservationCheck_write_to_logfile)
      call mpas_pool_add_config(configPool, 'config_AM_conservationCheck_restart_stream', &
config_AM_conservationCheck_restart_stream)

   end subroutine ocn_setup_nmlrec_AM_conservationCheck


   subroutine ocn_setup_nmlrec_baroclinic_channel(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      integer :: config_baroclinic_channel_vert_levels = 20
      logical :: config_baroclinic_channel_use_distances = .false.
      real (kind=RKIND) :: config_baroclinic_channel_surface_temperature = 13.100000
      real (kind=RKIND) :: config_baroclinic_channel_bottom_temperature = 10.100000
      real (kind=RKIND) :: config_baroclinic_channel_temperature_difference = 1.200000
      real (kind=RKIND) :: config_baroclinic_channel_gradient_width_frac = 0.080000
      real (kind=RKIND) :: config_baroclinic_channel_gradient_width_dist = 40000.000000
      real (kind=RKIND) :: config_baroclinic_channel_bottom_depth = 1000.000000
      real (kind=RKIND) :: config_baroclinic_channel_salinity = 35.000000
      real (kind=RKIND) :: config_baroclinic_channel_coriolis_parameter = -0.000120

      namelist /baroclinic_channel/ &
         config_baroclinic_channel_vert_levels, &
         config_baroclinic_channel_use_distances, &
         config_baroclinic_channel_surface_temperature, &
         config_baroclinic_channel_bottom_temperature, &
         config_baroclinic_channel_temperature_difference, &
         config_baroclinic_channel_gradient_width_frac, &
         config_baroclinic_channel_gradient_width_dist, &
         config_baroclinic_channel_bottom_depth, &
         config_baroclinic_channel_salinity, &
         config_baroclinic_channel_coriolis_parameter
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, baroclinic_channel, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_int(dminfo, config_baroclinic_channel_vert_levels)
         call mpas_dmpar_bcast_logical(dminfo, config_baroclinic_channel_use_distances)
         call mpas_dmpar_bcast_real(dminfo, config_baroclinic_channel_surface_temperature)
         call mpas_dmpar_bcast_real(dminfo, config_baroclinic_channel_bottom_temperature)
         call mpas_dmpar_bcast_real(dminfo, config_baroclinic_channel_temperature_difference)
         call mpas_dmpar_bcast_real(dminfo, config_baroclinic_channel_gradient_width_frac)
         call mpas_dmpar_bcast_real(dminfo, config_baroclinic_channel_gradient_width_dist)
         call mpas_dmpar_bcast_real(dminfo, config_baroclinic_channel_bottom_depth)
         call mpas_dmpar_bcast_real(dminfo, config_baroclinic_channel_salinity)
         call mpas_dmpar_bcast_real(dminfo, config_baroclinic_channel_coriolis_parameter)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record baroclinic_channel')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_baroclinic_channel_vert_levels = $i', &
intArgs=(/config_baroclinic_channel_vert_levels/))
            call mpas_log_write('        config_baroclinic_channel_use_distances = $l', &
logicArgs=(/config_baroclinic_channel_use_distances/))
            call mpas_log_write('        config_baroclinic_channel_surface_temperature = $r', &
realArgs=(/config_baroclinic_channel_surface_temperature/))
            call mpas_log_write('        config_baroclinic_channel_bottom_temperature = $r', &
realArgs=(/config_baroclinic_channel_bottom_temperature/))
            call mpas_log_write('        config_baroclinic_channel_temperature_difference = $r', &
realArgs=(/config_baroclinic_channel_temperature_difference/))
            call mpas_log_write('        config_baroclinic_channel_gradient_width_frac = $r', &
realArgs=(/config_baroclinic_channel_gradient_width_frac/))
            call mpas_log_write('        config_baroclinic_channel_gradient_width_dist = $r', &
realArgs=(/config_baroclinic_channel_gradient_width_dist/))
            call mpas_log_write('        config_baroclinic_channel_bottom_depth = $r', &
realArgs=(/config_baroclinic_channel_bottom_depth/))
            call mpas_log_write('        config_baroclinic_channel_salinity = $r', realArgs=(/config_baroclinic_channel_salinity/))
            call mpas_log_write('        config_baroclinic_channel_coriolis_parameter = $r', &
realArgs=(/config_baroclinic_channel_coriolis_parameter/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record baroclinic_channel.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_baroclinic_channel_vert_levels', config_baroclinic_channel_vert_levels)
      call mpas_pool_add_config(configPool, 'config_baroclinic_channel_use_distances', config_baroclinic_channel_use_distances)
      call mpas_pool_add_config(configPool, 'config_baroclinic_channel_surface_temperature', &
config_baroclinic_channel_surface_temperature)
      call mpas_pool_add_config(configPool, 'config_baroclinic_channel_bottom_temperature', &
config_baroclinic_channel_bottom_temperature)
      call mpas_pool_add_config(configPool, 'config_baroclinic_channel_temperature_difference', &
config_baroclinic_channel_temperature_difference)
      call mpas_pool_add_config(configPool, 'config_baroclinic_channel_gradient_width_frac', &
config_baroclinic_channel_gradient_width_frac)
      call mpas_pool_add_config(configPool, 'config_baroclinic_channel_gradient_width_dist', &
config_baroclinic_channel_gradient_width_dist)
      call mpas_pool_add_config(configPool, 'config_baroclinic_channel_bottom_depth', config_baroclinic_channel_bottom_depth)
      call mpas_pool_add_config(configPool, 'config_baroclinic_channel_salinity', config_baroclinic_channel_salinity)
      call mpas_pool_add_config(configPool, 'config_baroclinic_channel_coriolis_parameter', &
config_baroclinic_channel_coriolis_parameter)

   end subroutine ocn_setup_nmlrec_baroclinic_channel


   subroutine ocn_setup_nmlrec_lock_exchange(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      integer :: config_lock_exchange_vert_levels = 20
      real (kind=RKIND) :: config_lock_exchange_bottom_depth = 20.000000
      real (kind=RKIND) :: config_lock_exchange_cold_temperature = 5.000000
      real (kind=RKIND) :: config_lock_exchange_warm_temperature = 30.000000
      character (len=StrKIND) :: config_lock_exchange_direction = 'y'
      real (kind=RKIND) :: config_lock_exchange_salinity = 35.000000
      character (len=StrKIND) :: config_lock_exchange_layer_type = 'z-level'
      real (kind=RKIND) :: config_lock_exchange_isopycnal_min_thickness = 0.010000

      namelist /lock_exchange/ &
         config_lock_exchange_vert_levels, &
         config_lock_exchange_bottom_depth, &
         config_lock_exchange_cold_temperature, &
         config_lock_exchange_warm_temperature, &
         config_lock_exchange_direction, &
         config_lock_exchange_salinity, &
         config_lock_exchange_layer_type, &
         config_lock_exchange_isopycnal_min_thickness
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, lock_exchange, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_int(dminfo, config_lock_exchange_vert_levels)
         call mpas_dmpar_bcast_real(dminfo, config_lock_exchange_bottom_depth)
         call mpas_dmpar_bcast_real(dminfo, config_lock_exchange_cold_temperature)
         call mpas_dmpar_bcast_real(dminfo, config_lock_exchange_warm_temperature)
         call mpas_dmpar_bcast_char(dminfo, config_lock_exchange_direction)
         call mpas_dmpar_bcast_real(dminfo, config_lock_exchange_salinity)
         call mpas_dmpar_bcast_char(dminfo, config_lock_exchange_layer_type)
         call mpas_dmpar_bcast_real(dminfo, config_lock_exchange_isopycnal_min_thickness)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record lock_exchange')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_lock_exchange_vert_levels = $i', intArgs=(/config_lock_exchange_vert_levels/))
            call mpas_log_write('        config_lock_exchange_bottom_depth = $r', realArgs=(/config_lock_exchange_bottom_depth/))
            call mpas_log_write('        config_lock_exchange_cold_temperature = $r', &
realArgs=(/config_lock_exchange_cold_temperature/))
            call mpas_log_write('        config_lock_exchange_warm_temperature = $r', &
realArgs=(/config_lock_exchange_warm_temperature/))
            call mpas_log_write('        config_lock_exchange_direction = '&
//''//mpas_log_escape_dollars(config_lock_exchange_direction))
            call mpas_log_write('        config_lock_exchange_salinity = $r', realArgs=(/config_lock_exchange_salinity/))
            call mpas_log_write('        config_lock_exchange_layer_type = '&
//''//mpas_log_escape_dollars(config_lock_exchange_layer_type))
            call mpas_log_write('        config_lock_exchange_isopycnal_min_thickness = $r', &
realArgs=(/config_lock_exchange_isopycnal_min_thickness/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record lock_exchange.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_lock_exchange_vert_levels', config_lock_exchange_vert_levels)
      call mpas_pool_add_config(configPool, 'config_lock_exchange_bottom_depth', config_lock_exchange_bottom_depth)
      call mpas_pool_add_config(configPool, 'config_lock_exchange_cold_temperature', config_lock_exchange_cold_temperature)
      call mpas_pool_add_config(configPool, 'config_lock_exchange_warm_temperature', config_lock_exchange_warm_temperature)
      call mpas_pool_add_config(configPool, 'config_lock_exchange_direction', config_lock_exchange_direction)
      call mpas_pool_add_config(configPool, 'config_lock_exchange_salinity', config_lock_exchange_salinity)
      call mpas_pool_add_config(configPool, 'config_lock_exchange_layer_type', config_lock_exchange_layer_type)
      call mpas_pool_add_config(configPool, 'config_lock_exchange_isopycnal_min_thickness', &
config_lock_exchange_isopycnal_min_thickness)

   end subroutine ocn_setup_nmlrec_lock_exchange


   subroutine ocn_setup_nmlrec_internal_waves(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      integer :: config_internal_waves_vert_levels = 20
      logical :: config_internal_waves_use_distances = .false.
      real (kind=RKIND) :: config_internal_waves_surface_temperature = 20.100000
      real (kind=RKIND) :: config_internal_waves_bottom_temperature = 10.100000
      real (kind=RKIND) :: config_internal_waves_temperature_difference = 2.000000
      real (kind=RKIND) :: config_internal_waves_amplitude_width_frac = 0.330000
      real (kind=RKIND) :: config_internal_waves_amplitude_width_dist = 50000.000000
      real (kind=RKIND) :: config_internal_waves_bottom_depth = 500.000000
      real (kind=RKIND) :: config_internal_waves_salinity = 35.000000
      character (len=StrKIND) :: config_internal_waves_layer_type = 'z-level'
      real (kind=RKIND) :: config_internal_waves_isopycnal_displacement = 125.000000

      namelist /internal_waves/ &
         config_internal_waves_vert_levels, &
         config_internal_waves_use_distances, &
         config_internal_waves_surface_temperature, &
         config_internal_waves_bottom_temperature, &
         config_internal_waves_temperature_difference, &
         config_internal_waves_amplitude_width_frac, &
         config_internal_waves_amplitude_width_dist, &
         config_internal_waves_bottom_depth, &
         config_internal_waves_salinity, &
         config_internal_waves_layer_type, &
         config_internal_waves_isopycnal_displacement
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, internal_waves, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_int(dminfo, config_internal_waves_vert_levels)
         call mpas_dmpar_bcast_logical(dminfo, config_internal_waves_use_distances)
         call mpas_dmpar_bcast_real(dminfo, config_internal_waves_surface_temperature)
         call mpas_dmpar_bcast_real(dminfo, config_internal_waves_bottom_temperature)
         call mpas_dmpar_bcast_real(dminfo, config_internal_waves_temperature_difference)
         call mpas_dmpar_bcast_real(dminfo, config_internal_waves_amplitude_width_frac)
         call mpas_dmpar_bcast_real(dminfo, config_internal_waves_amplitude_width_dist)
         call mpas_dmpar_bcast_real(dminfo, config_internal_waves_bottom_depth)
         call mpas_dmpar_bcast_real(dminfo, config_internal_waves_salinity)
         call mpas_dmpar_bcast_char(dminfo, config_internal_waves_layer_type)
         call mpas_dmpar_bcast_real(dminfo, config_internal_waves_isopycnal_displacement)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record internal_waves')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_internal_waves_vert_levels = $i', intArgs=(/config_internal_waves_vert_levels/))
            call mpas_log_write('        config_internal_waves_use_distances = $l', &
logicArgs=(/config_internal_waves_use_distances/))
            call mpas_log_write('        config_internal_waves_surface_temperature = $r', &
realArgs=(/config_internal_waves_surface_temperature/))
            call mpas_log_write('        config_internal_waves_bottom_temperature = $r', &
realArgs=(/config_internal_waves_bottom_temperature/))
            call mpas_log_write('        config_internal_waves_temperature_difference = $r', &
realArgs=(/config_internal_waves_temperature_difference/))
            call mpas_log_write('        config_internal_waves_amplitude_width_frac = $r', &
realArgs=(/config_internal_waves_amplitude_width_frac/))
            call mpas_log_write('        config_internal_waves_amplitude_width_dist = $r', &
realArgs=(/config_internal_waves_amplitude_width_dist/))
            call mpas_log_write('        config_internal_waves_bottom_depth = $r', realArgs=(/config_internal_waves_bottom_depth/))
            call mpas_log_write('        config_internal_waves_salinity = $r', realArgs=(/config_internal_waves_salinity/))
            call mpas_log_write('        config_internal_waves_layer_type = '&
//''//mpas_log_escape_dollars(config_internal_waves_layer_type))
            call mpas_log_write('        config_internal_waves_isopycnal_displacement = $r', &
realArgs=(/config_internal_waves_isopycnal_displacement/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record internal_waves.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_internal_waves_vert_levels', config_internal_waves_vert_levels)
      call mpas_pool_add_config(configPool, 'config_internal_waves_use_distances', config_internal_waves_use_distances)
      call mpas_pool_add_config(configPool, 'config_internal_waves_surface_temperature', config_internal_waves_surface_temperature)
      call mpas_pool_add_config(configPool, 'config_internal_waves_bottom_temperature', config_internal_waves_bottom_temperature)
      call mpas_pool_add_config(configPool, 'config_internal_waves_temperature_difference', &
config_internal_waves_temperature_difference)
      call mpas_pool_add_config(configPool, 'config_internal_waves_amplitude_width_frac', &
config_internal_waves_amplitude_width_frac)
      call mpas_pool_add_config(configPool, 'config_internal_waves_amplitude_width_dist', &
config_internal_waves_amplitude_width_dist)
      call mpas_pool_add_config(configPool, 'config_internal_waves_bottom_depth', config_internal_waves_bottom_depth)
      call mpas_pool_add_config(configPool, 'config_internal_waves_salinity', config_internal_waves_salinity)
      call mpas_pool_add_config(configPool, 'config_internal_waves_layer_type', config_internal_waves_layer_type)
      call mpas_pool_add_config(configPool, 'config_internal_waves_isopycnal_displacement', &
config_internal_waves_isopycnal_displacement)

   end subroutine ocn_setup_nmlrec_internal_waves


   subroutine ocn_setup_nmlrec_overflow(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      integer :: config_overflow_vert_levels = 100
      logical :: config_overflow_use_distances = .false.
      real (kind=RKIND) :: config_overflow_bottom_depth = 2000.000000
      real (kind=RKIND) :: config_overflow_ridge_depth = 500.000000
      real (kind=RKIND) :: config_overflow_plug_temperature = 10.000000
      real (kind=RKIND) :: config_overflow_domain_temperature = 20.000000
      real (kind=RKIND) :: config_overflow_salinity = 35.000000
      real (kind=RKIND) :: config_overflow_plug_width_frac = 0.100000
      real (kind=RKIND) :: config_overflow_slope_center_frac = 0.200000
      real (kind=RKIND) :: config_overflow_slope_width_frac = 0.050000
      real (kind=RKIND) :: config_overflow_plug_width_dist = 20000.000000
      real (kind=RKIND) :: config_overflow_slope_center_dist = 40000.000000
      real (kind=RKIND) :: config_overflow_slope_width_dist = 7000.000000
      character (len=StrKIND) :: config_overflow_layer_type = 'z-level'
      real (kind=RKIND) :: config_overflow_isopycnal_min_thickness = 0.010000

      namelist /overflow/ &
         config_overflow_vert_levels, &
         config_overflow_use_distances, &
         config_overflow_bottom_depth, &
         config_overflow_ridge_depth, &
         config_overflow_plug_temperature, &
         config_overflow_domain_temperature, &
         config_overflow_salinity, &
         config_overflow_plug_width_frac, &
         config_overflow_slope_center_frac, &
         config_overflow_slope_width_frac, &
         config_overflow_plug_width_dist, &
         config_overflow_slope_center_dist, &
         config_overflow_slope_width_dist, &
         config_overflow_layer_type, &
         config_overflow_isopycnal_min_thickness
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, overflow, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_int(dminfo, config_overflow_vert_levels)
         call mpas_dmpar_bcast_logical(dminfo, config_overflow_use_distances)
         call mpas_dmpar_bcast_real(dminfo, config_overflow_bottom_depth)
         call mpas_dmpar_bcast_real(dminfo, config_overflow_ridge_depth)
         call mpas_dmpar_bcast_real(dminfo, config_overflow_plug_temperature)
         call mpas_dmpar_bcast_real(dminfo, config_overflow_domain_temperature)
         call mpas_dmpar_bcast_real(dminfo, config_overflow_salinity)
         call mpas_dmpar_bcast_real(dminfo, config_overflow_plug_width_frac)
         call mpas_dmpar_bcast_real(dminfo, config_overflow_slope_center_frac)
         call mpas_dmpar_bcast_real(dminfo, config_overflow_slope_width_frac)
         call mpas_dmpar_bcast_real(dminfo, config_overflow_plug_width_dist)
         call mpas_dmpar_bcast_real(dminfo, config_overflow_slope_center_dist)
         call mpas_dmpar_bcast_real(dminfo, config_overflow_slope_width_dist)
         call mpas_dmpar_bcast_char(dminfo, config_overflow_layer_type)
         call mpas_dmpar_bcast_real(dminfo, config_overflow_isopycnal_min_thickness)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record overflow')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_overflow_vert_levels = $i', intArgs=(/config_overflow_vert_levels/))
            call mpas_log_write('        config_overflow_use_distances = $l', logicArgs=(/config_overflow_use_distances/))
            call mpas_log_write('        config_overflow_bottom_depth = $r', realArgs=(/config_overflow_bottom_depth/))
            call mpas_log_write('        config_overflow_ridge_depth = $r', realArgs=(/config_overflow_ridge_depth/))
            call mpas_log_write('        config_overflow_plug_temperature = $r', realArgs=(/config_overflow_plug_temperature/))
            call mpas_log_write('        config_overflow_domain_temperature = $r', realArgs=(/config_overflow_domain_temperature/))
            call mpas_log_write('        config_overflow_salinity = $r', realArgs=(/config_overflow_salinity/))
            call mpas_log_write('        config_overflow_plug_width_frac = $r', realArgs=(/config_overflow_plug_width_frac/))
            call mpas_log_write('        config_overflow_slope_center_frac = $r', realArgs=(/config_overflow_slope_center_frac/))
            call mpas_log_write('        config_overflow_slope_width_frac = $r', realArgs=(/config_overflow_slope_width_frac/))
            call mpas_log_write('        config_overflow_plug_width_dist = $r', realArgs=(/config_overflow_plug_width_dist/))
            call mpas_log_write('        config_overflow_slope_center_dist = $r', realArgs=(/config_overflow_slope_center_dist/))
            call mpas_log_write('        config_overflow_slope_width_dist = $r', realArgs=(/config_overflow_slope_width_dist/))
            call mpas_log_write('        config_overflow_layer_type = '//mpas_log_escape_dollars(config_overflow_layer_type))
            call mpas_log_write('        config_overflow_isopycnal_min_thickness = $r', &
realArgs=(/config_overflow_isopycnal_min_thickness/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record overflow.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_overflow_vert_levels', config_overflow_vert_levels)
      call mpas_pool_add_config(configPool, 'config_overflow_use_distances', config_overflow_use_distances)
      call mpas_pool_add_config(configPool, 'config_overflow_bottom_depth', config_overflow_bottom_depth)
      call mpas_pool_add_config(configPool, 'config_overflow_ridge_depth', config_overflow_ridge_depth)
      call mpas_pool_add_config(configPool, 'config_overflow_plug_temperature', config_overflow_plug_temperature)
      call mpas_pool_add_config(configPool, 'config_overflow_domain_temperature', config_overflow_domain_temperature)
      call mpas_pool_add_config(configPool, 'config_overflow_salinity', config_overflow_salinity)
      call mpas_pool_add_config(configPool, 'config_overflow_plug_width_frac', config_overflow_plug_width_frac)
      call mpas_pool_add_config(configPool, 'config_overflow_slope_center_frac', config_overflow_slope_center_frac)
      call mpas_pool_add_config(configPool, 'config_overflow_slope_width_frac', config_overflow_slope_width_frac)
      call mpas_pool_add_config(configPool, 'config_overflow_plug_width_dist', config_overflow_plug_width_dist)
      call mpas_pool_add_config(configPool, 'config_overflow_slope_center_dist', config_overflow_slope_center_dist)
      call mpas_pool_add_config(configPool, 'config_overflow_slope_width_dist', config_overflow_slope_width_dist)
      call mpas_pool_add_config(configPool, 'config_overflow_layer_type', config_overflow_layer_type)
      call mpas_pool_add_config(configPool, 'config_overflow_isopycnal_min_thickness', config_overflow_isopycnal_min_thickness)

   end subroutine ocn_setup_nmlrec_overflow


   subroutine ocn_setup_nmlrec_dam_break(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      integer :: config_dam_break_vert_levels = 1
      real (kind=RKIND) :: config_dam_break_eta0 = 0.600000
      real (kind=RKIND) :: config_dam_break_dc = 0.040000
      real (kind=RKIND) :: config_dam_break_R0 = 24.200000
      real (kind=RKIND) :: config_dam_break_Xl = 1.000000
      real (kind=RKIND) :: config_dam_break_Yl = 2.000000
      real (kind=RKIND) :: config_dam_break_Inlet = 0.400000

      namelist /dam_break/ &
         config_dam_break_vert_levels, &
         config_dam_break_eta0, &
         config_dam_break_dc, &
         config_dam_break_R0, &
         config_dam_break_Xl, &
         config_dam_break_Yl, &
         config_dam_break_Inlet
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, dam_break, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_int(dminfo, config_dam_break_vert_levels)
         call mpas_dmpar_bcast_real(dminfo, config_dam_break_eta0)
         call mpas_dmpar_bcast_real(dminfo, config_dam_break_dc)
         call mpas_dmpar_bcast_real(dminfo, config_dam_break_R0)
         call mpas_dmpar_bcast_real(dminfo, config_dam_break_Xl)
         call mpas_dmpar_bcast_real(dminfo, config_dam_break_Yl)
         call mpas_dmpar_bcast_real(dminfo, config_dam_break_Inlet)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record dam_break')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_dam_break_vert_levels = $i', intArgs=(/config_dam_break_vert_levels/))
            call mpas_log_write('        config_dam_break_eta0 = $r', realArgs=(/config_dam_break_eta0/))
            call mpas_log_write('        config_dam_break_dc = $r', realArgs=(/config_dam_break_dc/))
            call mpas_log_write('        config_dam_break_R0 = $r', realArgs=(/config_dam_break_R0/))
            call mpas_log_write('        config_dam_break_Xl = $r', realArgs=(/config_dam_break_Xl/))
            call mpas_log_write('        config_dam_break_Yl = $r', realArgs=(/config_dam_break_Yl/))
            call mpas_log_write('        config_dam_break_Inlet = $r', realArgs=(/config_dam_break_Inlet/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record dam_break.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_dam_break_vert_levels', config_dam_break_vert_levels)
      call mpas_pool_add_config(configPool, 'config_dam_break_eta0', config_dam_break_eta0)
      call mpas_pool_add_config(configPool, 'config_dam_break_dc', config_dam_break_dc)
      call mpas_pool_add_config(configPool, 'config_dam_break_R0', config_dam_break_R0)
      call mpas_pool_add_config(configPool, 'config_dam_break_Xl', config_dam_break_Xl)
      call mpas_pool_add_config(configPool, 'config_dam_break_Yl', config_dam_break_Yl)
      call mpas_pool_add_config(configPool, 'config_dam_break_Inlet', config_dam_break_Inlet)

   end subroutine ocn_setup_nmlrec_dam_break


   subroutine ocn_setup_nmlrec_global_ocean(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      real (kind=RKIND) :: config_global_ocean_minimum_depth = 15.000000
      character (len=StrKIND) :: config_global_ocean_depth_file = 'vertical_grid.nc'
      character (len=StrKIND) :: config_global_ocean_depth_dimname = 'nVertLevels'
      character (len=StrKIND) :: config_global_ocean_depth_varname = 'refMidDepth'
      real (kind=RKIND) :: config_global_ocean_depth_conversion_factor = 1.000000
      character (len=StrKIND) :: config_global_ocean_temperature_file = 'none'
      character (len=StrKIND) :: config_global_ocean_salinity_file = 'none'
      character (len=StrKIND) :: config_global_ocean_tracer_nlat_dimname = 'none'
      character (len=StrKIND) :: config_global_ocean_tracer_nlon_dimname = 'none'
      character (len=StrKIND) :: config_global_ocean_tracer_ndepth_dimname = 'none'
      real (kind=RKIND) :: config_global_ocean_tracer_depth_conversion_factor = 1.000000
      integer :: config_global_ocean_tracer_vert_levels = -1
      character (len=StrKIND) :: config_global_ocean_temperature_varname = 'none'
      character (len=StrKIND) :: config_global_ocean_salinity_varname = 'none'
      logical :: config_global_ocean_tracer_latlon_degrees = .true.
      character (len=StrKIND) :: config_global_ocean_tracer_lat_varname = 'none'
      character (len=StrKIND) :: config_global_ocean_tracer_lon_varname = 'none'
      character (len=StrKIND) :: config_global_ocean_tracer_depth_varname = 'none'
      character (len=StrKIND) :: config_global_ocean_tracer_method = 'bilinear_interpolation'
      integer :: config_global_ocean_smooth_TS_iterations = 0
      character (len=StrKIND) :: config_global_ocean_swData_file = 'none'
      character (len=StrKIND) :: config_global_ocean_swData_nlat_dimname = 'none'
      character (len=StrKIND) :: config_global_ocean_swData_nlon_dimname = 'none'
      character (len=StrKIND) :: config_global_ocean_swData_lat_varname = 'none'
      character (len=StrKIND) :: config_global_ocean_swData_lon_varname = 'none'
      logical :: config_global_ocean_swData_latlon_degrees = .true.
      character (len=StrKIND) :: config_global_ocean_swData_method = 'bilinear_interpolation'
      character (len=StrKIND) :: config_global_ocean_chlorophyll_varname = 'none'
      character (len=StrKIND) :: config_global_ocean_zenithAngle_varname = 'none'
      character (len=StrKIND) :: config_global_ocean_clearSky_varname = 'none'
      real (kind=RKIND) :: config_global_ocean_piston_velocity = 0.000050
      real (kind=RKIND) :: config_global_ocean_interior_restore_rate = 0.000000
      character (len=StrKIND) :: config_global_ocean_topography_source = 'latlon_file'
      character (len=StrKIND) :: config_global_ocean_topography_file = 'none'
      character (len=StrKIND) :: config_global_ocean_topography_nlat_dimname = 'none'
      character (len=StrKIND) :: config_global_ocean_topography_nlon_dimname = 'none'
      logical :: config_global_ocean_topography_latlon_degrees = .true.
      character (len=StrKIND) :: config_global_ocean_topography_lat_varname = 'none'
      character (len=StrKIND) :: config_global_ocean_topography_lon_varname = 'none'
      character (len=StrKIND) :: config_global_ocean_topography_varname = 'none'
      logical :: config_global_ocean_topography_has_ocean_frac = .false.
      character (len=StrKIND) :: config_global_ocean_topography_ocean_frac_varname = 'none'
      character (len=StrKIND) :: config_global_ocean_topography_method = 'bilinear_interpolation'
      logical :: config_global_ocean_fill_bathymetry_holes = .true.
      integer :: config_global_ocean_topography_smooth_iterations = 0
      real (kind=RKIND) :: config_global_ocean_topography_smooth_weight = 0.900000
      logical :: config_global_ocean_deepen_critical_passages = .true.
      logical :: config_global_ocean_depress_by_land_ice = .false.
      character (len=StrKIND) :: config_global_ocean_land_ice_topo_file = 'none'
      character (len=StrKIND) :: config_global_ocean_land_ice_topo_nlat_dimname = 'none'
      character (len=StrKIND) :: config_global_ocean_land_ice_topo_nlon_dimname = 'none'
      logical :: config_global_ocean_land_ice_topo_latlon_degrees = .true.
      character (len=StrKIND) :: config_global_ocean_land_ice_topo_lat_varname = 'none'
      character (len=StrKIND) :: config_global_ocean_land_ice_topo_lon_varname = 'none'
      character (len=StrKIND) :: config_global_ocean_land_ice_topo_thickness_varname = 'none'
      character (len=StrKIND) :: config_global_ocean_land_ice_topo_draft_varname = 'none'
      character (len=StrKIND) :: config_global_ocean_land_ice_topo_ice_frac_varname = 'none'
      character (len=StrKIND) :: config_global_ocean_land_ice_topo_grounded_frac_varname = 'none'
      logical :: config_global_ocean_use_constant_land_ice_cavity_temperature = .false.
      real (kind=RKIND) :: config_global_ocean_constant_land_ice_cavity_temperature = -1.800000
      logical :: config_global_ocean_cull_inland_seas = .true.
      character (len=StrKIND) :: config_global_ocean_windstress_file = 'none'
      character (len=StrKIND) :: config_global_ocean_windstress_nlat_dimname = 'none'
      character (len=StrKIND) :: config_global_ocean_windstress_nlon_dimname = 'none'
      logical :: config_global_ocean_windstress_latlon_degrees = .true.
      character (len=StrKIND) :: config_global_ocean_windstress_lat_varname = 'none'
      character (len=StrKIND) :: config_global_ocean_windstress_lon_varname = 'none'
      character (len=StrKIND) :: config_global_ocean_windstress_zonal_varname = 'none'
      character (len=StrKIND) :: config_global_ocean_windstress_meridional_varname = 'none'
      character (len=StrKIND) :: config_global_ocean_windstress_method = 'bilinear_interpolation'
      real (kind=RKIND) :: config_global_ocean_windstress_conversion_factor = 1.000000
      character (len=StrKIND) :: config_global_ocean_ecosys_file = 'unknown'
      character (len=StrKIND) :: config_global_ocean_ecosys_forcing_file = 'unknown'
      character (len=StrKIND) :: config_global_ocean_ecosys_nlat_dimname = 'none'
      character (len=StrKIND) :: config_global_ocean_ecosys_nlon_dimname = 'none'
      character (len=StrKIND) :: config_global_ocean_ecosys_ndepth_dimname = 'none'
      real (kind=RKIND) :: config_global_ocean_ecosys_depth_conversion_factor = 1.000000
      integer :: config_global_ocean_ecosys_vert_levels = -1
      character (len=StrKIND) :: config_global_ocean_ecosys_lat_varname = 'none'
      character (len=StrKIND) :: config_global_ocean_ecosys_lon_varname = 'none'
      character (len=StrKIND) :: config_global_ocean_ecosys_depth_varname = 'none'
      logical :: config_global_ocean_ecosys_latlon_degrees = .true.
      character (len=StrKIND) :: config_global_ocean_ecosys_method = 'bilinear_interpolation'
      character (len=StrKIND) :: config_global_ocean_ecosys_forcing_time_dimname = 'none'
      integer :: config_global_ocean_smooth_ecosys_iterations = 0

      namelist /global_ocean/ &
         config_global_ocean_minimum_depth, &
         config_global_ocean_depth_file, &
         config_global_ocean_depth_dimname, &
         config_global_ocean_depth_varname, &
         config_global_ocean_depth_conversion_factor, &
         config_global_ocean_temperature_file, &
         config_global_ocean_salinity_file, &
         config_global_ocean_tracer_nlat_dimname, &
         config_global_ocean_tracer_nlon_dimname, &
         config_global_ocean_tracer_ndepth_dimname, &
         config_global_ocean_tracer_depth_conversion_factor, &
         config_global_ocean_tracer_vert_levels, &
         config_global_ocean_temperature_varname, &
         config_global_ocean_salinity_varname, &
         config_global_ocean_tracer_latlon_degrees, &
         config_global_ocean_tracer_lat_varname, &
         config_global_ocean_tracer_lon_varname, &
         config_global_ocean_tracer_depth_varname, &
         config_global_ocean_tracer_method, &
         config_global_ocean_smooth_TS_iterations, &
         config_global_ocean_swData_file, &
         config_global_ocean_swData_nlat_dimname, &
         config_global_ocean_swData_nlon_dimname, &
         config_global_ocean_swData_lat_varname, &
         config_global_ocean_swData_lon_varname, &
         config_global_ocean_swData_latlon_degrees, &
         config_global_ocean_swData_method, &
         config_global_ocean_chlorophyll_varname, &
         config_global_ocean_zenithAngle_varname, &
         config_global_ocean_clearSky_varname, &
         config_global_ocean_piston_velocity, &
         config_global_ocean_interior_restore_rate, &
         config_global_ocean_topography_source, &
         config_global_ocean_topography_file, &
         config_global_ocean_topography_nlat_dimname, &
         config_global_ocean_topography_nlon_dimname, &
         config_global_ocean_topography_latlon_degrees, &
         config_global_ocean_topography_lat_varname, &
         config_global_ocean_topography_lon_varname, &
         config_global_ocean_topography_varname, &
         config_global_ocean_topography_has_ocean_frac, &
         config_global_ocean_topography_ocean_frac_varname, &
         config_global_ocean_topography_method, &
         config_global_ocean_fill_bathymetry_holes, &
         config_global_ocean_topography_smooth_iterations, &
         config_global_ocean_topography_smooth_weight, &
         config_global_ocean_deepen_critical_passages, &
         config_global_ocean_depress_by_land_ice, &
         config_global_ocean_land_ice_topo_file, &
         config_global_ocean_land_ice_topo_nlat_dimname, &
         config_global_ocean_land_ice_topo_nlon_dimname, &
         config_global_ocean_land_ice_topo_latlon_degrees, &
         config_global_ocean_land_ice_topo_lat_varname, &
         config_global_ocean_land_ice_topo_lon_varname, &
         config_global_ocean_land_ice_topo_thickness_varname, &
         config_global_ocean_land_ice_topo_draft_varname, &
         config_global_ocean_land_ice_topo_ice_frac_varname, &
         config_global_ocean_land_ice_topo_grounded_frac_varname, &
         config_global_ocean_use_constant_land_ice_cavity_temperature, &
         config_global_ocean_constant_land_ice_cavity_temperature, &
         config_global_ocean_cull_inland_seas, &
         config_global_ocean_windstress_file, &
         config_global_ocean_windstress_nlat_dimname, &
         config_global_ocean_windstress_nlon_dimname, &
         config_global_ocean_windstress_latlon_degrees, &
         config_global_ocean_windstress_lat_varname, &
         config_global_ocean_windstress_lon_varname, &
         config_global_ocean_windstress_zonal_varname, &
         config_global_ocean_windstress_meridional_varname, &
         config_global_ocean_windstress_method, &
         config_global_ocean_windstress_conversion_factor, &
         config_global_ocean_ecosys_file, &
         config_global_ocean_ecosys_forcing_file, &
         config_global_ocean_ecosys_nlat_dimname, &
         config_global_ocean_ecosys_nlon_dimname, &
         config_global_ocean_ecosys_ndepth_dimname, &
         config_global_ocean_ecosys_depth_conversion_factor, &
         config_global_ocean_ecosys_vert_levels, &
         config_global_ocean_ecosys_lat_varname, &
         config_global_ocean_ecosys_lon_varname, &
         config_global_ocean_ecosys_depth_varname, &
         config_global_ocean_ecosys_latlon_degrees, &
         config_global_ocean_ecosys_method, &
         config_global_ocean_ecosys_forcing_time_dimname, &
         config_global_ocean_smooth_ecosys_iterations
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, global_ocean, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_real(dminfo, config_global_ocean_minimum_depth)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_depth_file)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_depth_dimname)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_depth_varname)
         call mpas_dmpar_bcast_real(dminfo, config_global_ocean_depth_conversion_factor)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_temperature_file)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_salinity_file)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_tracer_nlat_dimname)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_tracer_nlon_dimname)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_tracer_ndepth_dimname)
         call mpas_dmpar_bcast_real(dminfo, config_global_ocean_tracer_depth_conversion_factor)
         call mpas_dmpar_bcast_int(dminfo, config_global_ocean_tracer_vert_levels)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_temperature_varname)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_salinity_varname)
         call mpas_dmpar_bcast_logical(dminfo, config_global_ocean_tracer_latlon_degrees)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_tracer_lat_varname)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_tracer_lon_varname)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_tracer_depth_varname)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_tracer_method)
         call mpas_dmpar_bcast_int(dminfo, config_global_ocean_smooth_TS_iterations)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_swData_file)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_swData_nlat_dimname)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_swData_nlon_dimname)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_swData_lat_varname)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_swData_lon_varname)
         call mpas_dmpar_bcast_logical(dminfo, config_global_ocean_swData_latlon_degrees)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_swData_method)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_chlorophyll_varname)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_zenithAngle_varname)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_clearSky_varname)
         call mpas_dmpar_bcast_real(dminfo, config_global_ocean_piston_velocity)
         call mpas_dmpar_bcast_real(dminfo, config_global_ocean_interior_restore_rate)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_topography_source)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_topography_file)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_topography_nlat_dimname)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_topography_nlon_dimname)
         call mpas_dmpar_bcast_logical(dminfo, config_global_ocean_topography_latlon_degrees)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_topography_lat_varname)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_topography_lon_varname)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_topography_varname)
         call mpas_dmpar_bcast_logical(dminfo, config_global_ocean_topography_has_ocean_frac)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_topography_ocean_frac_varname)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_topography_method)
         call mpas_dmpar_bcast_logical(dminfo, config_global_ocean_fill_bathymetry_holes)
         call mpas_dmpar_bcast_int(dminfo, config_global_ocean_topography_smooth_iterations)
         call mpas_dmpar_bcast_real(dminfo, config_global_ocean_topography_smooth_weight)
         call mpas_dmpar_bcast_logical(dminfo, config_global_ocean_deepen_critical_passages)
         call mpas_dmpar_bcast_logical(dminfo, config_global_ocean_depress_by_land_ice)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_land_ice_topo_file)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_land_ice_topo_nlat_dimname)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_land_ice_topo_nlon_dimname)
         call mpas_dmpar_bcast_logical(dminfo, config_global_ocean_land_ice_topo_latlon_degrees)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_land_ice_topo_lat_varname)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_land_ice_topo_lon_varname)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_land_ice_topo_thickness_varname)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_land_ice_topo_draft_varname)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_land_ice_topo_ice_frac_varname)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_land_ice_topo_grounded_frac_varname)
         call mpas_dmpar_bcast_logical(dminfo, config_global_ocean_use_constant_land_ice_cavity_temperature)
         call mpas_dmpar_bcast_real(dminfo, config_global_ocean_constant_land_ice_cavity_temperature)
         call mpas_dmpar_bcast_logical(dminfo, config_global_ocean_cull_inland_seas)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_windstress_file)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_windstress_nlat_dimname)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_windstress_nlon_dimname)
         call mpas_dmpar_bcast_logical(dminfo, config_global_ocean_windstress_latlon_degrees)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_windstress_lat_varname)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_windstress_lon_varname)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_windstress_zonal_varname)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_windstress_meridional_varname)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_windstress_method)
         call mpas_dmpar_bcast_real(dminfo, config_global_ocean_windstress_conversion_factor)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_ecosys_file)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_ecosys_forcing_file)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_ecosys_nlat_dimname)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_ecosys_nlon_dimname)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_ecosys_ndepth_dimname)
         call mpas_dmpar_bcast_real(dminfo, config_global_ocean_ecosys_depth_conversion_factor)
         call mpas_dmpar_bcast_int(dminfo, config_global_ocean_ecosys_vert_levels)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_ecosys_lat_varname)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_ecosys_lon_varname)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_ecosys_depth_varname)
         call mpas_dmpar_bcast_logical(dminfo, config_global_ocean_ecosys_latlon_degrees)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_ecosys_method)
         call mpas_dmpar_bcast_char(dminfo, config_global_ocean_ecosys_forcing_time_dimname)
         call mpas_dmpar_bcast_int(dminfo, config_global_ocean_smooth_ecosys_iterations)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record global_ocean')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_global_ocean_minimum_depth = $r', realArgs=(/config_global_ocean_minimum_depth/))
            call mpas_log_write('        config_global_ocean_depth_file = '&
//''//mpas_log_escape_dollars(config_global_ocean_depth_file))
            call mpas_log_write('        config_global_ocean_depth_dimname = '&
//''//mpas_log_escape_dollars(config_global_ocean_depth_dimname))
            call mpas_log_write('        config_global_ocean_depth_varname = '&
//''//mpas_log_escape_dollars(config_global_ocean_depth_varname))
            call mpas_log_write('        config_global_ocean_depth_conversion_factor = $r', &
realArgs=(/config_global_ocean_depth_conversion_factor/))
            call mpas_log_write('        config_global_ocean_temperature_file = '&
//''//mpas_log_escape_dollars(config_global_ocean_temperature_file))
            call mpas_log_write('        config_global_ocean_salinity_file = '&
//''//mpas_log_escape_dollars(config_global_ocean_salinity_file))
            call mpas_log_write('        config_global_ocean_tracer_nlat_dimname = '&
//''//mpas_log_escape_dollars(config_global_ocean_tracer_nlat_dimname))
            call mpas_log_write('        config_global_ocean_tracer_nlon_dimname = '&
//''//mpas_log_escape_dollars(config_global_ocean_tracer_nlon_dimname))
            call mpas_log_write('        config_global_ocean_tracer_ndepth_dimname = '&
//''//mpas_log_escape_dollars(config_global_ocean_tracer_ndepth_dimname))
            call mpas_log_write('        config_global_ocean_tracer_depth_conversion_factor = $r', &
realArgs=(/config_global_ocean_tracer_depth_conversion_factor/))
            call mpas_log_write('        config_global_ocean_tracer_vert_levels = $i', &
intArgs=(/config_global_ocean_tracer_vert_levels/))
            call mpas_log_write('        config_global_ocean_temperature_varname = '&
//''//mpas_log_escape_dollars(config_global_ocean_temperature_varname))
            call mpas_log_write('        config_global_ocean_salinity_varname = '&
//''//mpas_log_escape_dollars(config_global_ocean_salinity_varname))
            call mpas_log_write('        config_global_ocean_tracer_latlon_degrees = $l', &
logicArgs=(/config_global_ocean_tracer_latlon_degrees/))
            call mpas_log_write('        config_global_ocean_tracer_lat_varname = '&
//''//mpas_log_escape_dollars(config_global_ocean_tracer_lat_varname))
            call mpas_log_write('        config_global_ocean_tracer_lon_varname = '&
//''//mpas_log_escape_dollars(config_global_ocean_tracer_lon_varname))
            call mpas_log_write('        config_global_ocean_tracer_depth_varname = '&
//''//mpas_log_escape_dollars(config_global_ocean_tracer_depth_varname))
            call mpas_log_write('        config_global_ocean_tracer_method = '&
//''//mpas_log_escape_dollars(config_global_ocean_tracer_method))
            call mpas_log_write('        config_global_ocean_smooth_TS_iterations = $i', &
intArgs=(/config_global_ocean_smooth_TS_iterations/))
            call mpas_log_write('        config_global_ocean_swData_file = '&
//''//mpas_log_escape_dollars(config_global_ocean_swData_file))
            call mpas_log_write('        config_global_ocean_swData_nlat_dimname = '&
//''//mpas_log_escape_dollars(config_global_ocean_swData_nlat_dimname))
            call mpas_log_write('        config_global_ocean_swData_nlon_dimname = '&
//''//mpas_log_escape_dollars(config_global_ocean_swData_nlon_dimname))
            call mpas_log_write('        config_global_ocean_swData_lat_varname = '&
//''//mpas_log_escape_dollars(config_global_ocean_swData_lat_varname))
            call mpas_log_write('        config_global_ocean_swData_lon_varname = '&
//''//mpas_log_escape_dollars(config_global_ocean_swData_lon_varname))
            call mpas_log_write('        config_global_ocean_swData_latlon_degrees = $l', &
logicArgs=(/config_global_ocean_swData_latlon_degrees/))
            call mpas_log_write('        config_global_ocean_swData_method = '&
//''//mpas_log_escape_dollars(config_global_ocean_swData_method))
            call mpas_log_write('        config_global_ocean_chlorophyll_varname = '&
//''//mpas_log_escape_dollars(config_global_ocean_chlorophyll_varname))
            call mpas_log_write('        config_global_ocean_zenithAngle_varname = '&
//''//mpas_log_escape_dollars(config_global_ocean_zenithAngle_varname))
            call mpas_log_write('        config_global_ocean_clearSky_varname = '&
//''//mpas_log_escape_dollars(config_global_ocean_clearSky_varname))
            call mpas_log_write('        config_global_ocean_piston_velocity = $r', &
realArgs=(/config_global_ocean_piston_velocity/))
            call mpas_log_write('        config_global_ocean_interior_restore_rate = $r', &
realArgs=(/config_global_ocean_interior_restore_rate/))
            call mpas_log_write('        config_global_ocean_topography_source = '&
//''//mpas_log_escape_dollars(config_global_ocean_topography_source))
            call mpas_log_write('        config_global_ocean_topography_file = '&
//''//mpas_log_escape_dollars(config_global_ocean_topography_file))
            call mpas_log_write('        config_global_ocean_topography_nlat_dimname = '&
//''//mpas_log_escape_dollars(config_global_ocean_topography_nlat_dimname))
            call mpas_log_write('        config_global_ocean_topography_nlon_dimname = '&
//''//mpas_log_escape_dollars(config_global_ocean_topography_nlon_dimname))
            call mpas_log_write('        config_global_ocean_topography_latlon_degrees = $l', &
logicArgs=(/config_global_ocean_topography_latlon_degrees/))
            call mpas_log_write('        config_global_ocean_topography_lat_varname = '&
//''//mpas_log_escape_dollars(config_global_ocean_topography_lat_varname))
            call mpas_log_write('        config_global_ocean_topography_lon_varname = '&
//''//mpas_log_escape_dollars(config_global_ocean_topography_lon_varname))
            call mpas_log_write('        config_global_ocean_topography_varname = '&
//''//mpas_log_escape_dollars(config_global_ocean_topography_varname))
            call mpas_log_write('        config_global_ocean_topography_has_ocean_frac = $l', &
logicArgs=(/config_global_ocean_topography_has_ocean_frac/))
            call mpas_log_write('        config_global_ocean_topography_ocean_frac_varname = '&
//''//mpas_log_escape_dollars(config_global_ocean_topography_ocean_frac_varname))
            call mpas_log_write('        config_global_ocean_topography_method = '&
//''//mpas_log_escape_dollars(config_global_ocean_topography_method))
            call mpas_log_write('        config_global_ocean_fill_bathymetry_holes = $l', &
logicArgs=(/config_global_ocean_fill_bathymetry_holes/))
            call mpas_log_write('        config_global_ocean_topography_smooth_iterations = $i', &
intArgs=(/config_global_ocean_topography_smooth_iterations/))
            call mpas_log_write('        config_global_ocean_topography_smooth_weight = $r', &
realArgs=(/config_global_ocean_topography_smooth_weight/))
            call mpas_log_write('        config_global_ocean_deepen_critical_passages = $l', &
logicArgs=(/config_global_ocean_deepen_critical_passages/))
            call mpas_log_write('        config_global_ocean_depress_by_land_ice = $l', &
logicArgs=(/config_global_ocean_depress_by_land_ice/))
            call mpas_log_write('        config_global_ocean_land_ice_topo_file = '&
//''//mpas_log_escape_dollars(config_global_ocean_land_ice_topo_file))
            call mpas_log_write('        config_global_ocean_land_ice_topo_nlat_dimname = '&
//''//mpas_log_escape_dollars(config_global_ocean_land_ice_topo_nlat_dimname))
            call mpas_log_write('        config_global_ocean_land_ice_topo_nlon_dimname = '&
//''//mpas_log_escape_dollars(config_global_ocean_land_ice_topo_nlon_dimname))
            call mpas_log_write('        config_global_ocean_land_ice_topo_latlon_degrees = $l', &
logicArgs=(/config_global_ocean_land_ice_topo_latlon_degrees/))
            call mpas_log_write('        config_global_ocean_land_ice_topo_lat_varname = '&
//''//mpas_log_escape_dollars(config_global_ocean_land_ice_topo_lat_varname))
            call mpas_log_write('        config_global_ocean_land_ice_topo_lon_varname = '&
//''//mpas_log_escape_dollars(config_global_ocean_land_ice_topo_lon_varname))
            call mpas_log_write('        config_global_ocean_land_ice_topo_thickness_varname = '&
//''//mpas_log_escape_dollars(config_global_ocean_land_ice_topo_thickness_varname))
            call mpas_log_write('        config_global_ocean_land_ice_topo_draft_varname = '&
//''//mpas_log_escape_dollars(config_global_ocean_land_ice_topo_draft_varname))
            call mpas_log_write('        config_global_ocean_land_ice_topo_ice_frac_varname = '&
//''//mpas_log_escape_dollars(config_global_ocean_land_ice_topo_ice_frac_varname))
            call mpas_log_write('        config_global_ocean_land_ice_topo_grounded_frac_varname = '&
//''//mpas_log_escape_dollars(config_global_ocean_land_ice_topo_grounded_frac_varname))
            call mpas_log_write('        config_global_ocean_use_constant_land_ice_cavity_temperature = $l', &
logicArgs=(/config_global_ocean_use_constant_land_ice_cavity_temperature/))
            call mpas_log_write('        config_global_ocean_constant_land_ice_cavity_temperature = $r', &
realArgs=(/config_global_ocean_constant_land_ice_cavity_temperature/))
            call mpas_log_write('        config_global_ocean_cull_inland_seas = $l', &
logicArgs=(/config_global_ocean_cull_inland_seas/))
            call mpas_log_write('        config_global_ocean_windstress_file = '&
//''//mpas_log_escape_dollars(config_global_ocean_windstress_file))
            call mpas_log_write('        config_global_ocean_windstress_nlat_dimname = '&
//''//mpas_log_escape_dollars(config_global_ocean_windstress_nlat_dimname))
            call mpas_log_write('        config_global_ocean_windstress_nlon_dimname = '&
//''//mpas_log_escape_dollars(config_global_ocean_windstress_nlon_dimname))
            call mpas_log_write('        config_global_ocean_windstress_latlon_degrees = $l', &
logicArgs=(/config_global_ocean_windstress_latlon_degrees/))
            call mpas_log_write('        config_global_ocean_windstress_lat_varname = '&
//''//mpas_log_escape_dollars(config_global_ocean_windstress_lat_varname))
            call mpas_log_write('        config_global_ocean_windstress_lon_varname = '&
//''//mpas_log_escape_dollars(config_global_ocean_windstress_lon_varname))
            call mpas_log_write('        config_global_ocean_windstress_zonal_varname = '&
//''//mpas_log_escape_dollars(config_global_ocean_windstress_zonal_varname))
            call mpas_log_write('        config_global_ocean_windstress_meridional_varname = '&
//''//mpas_log_escape_dollars(config_global_ocean_windstress_meridional_varname))
            call mpas_log_write('        config_global_ocean_windstress_method = '&
//''//mpas_log_escape_dollars(config_global_ocean_windstress_method))
            call mpas_log_write('        config_global_ocean_windstress_conversion_factor = $r', &
realArgs=(/config_global_ocean_windstress_conversion_factor/))
            call mpas_log_write('        config_global_ocean_ecosys_file = '&
//''//mpas_log_escape_dollars(config_global_ocean_ecosys_file))
            call mpas_log_write('        config_global_ocean_ecosys_forcing_file = '&
//''//mpas_log_escape_dollars(config_global_ocean_ecosys_forcing_file))
            call mpas_log_write('        config_global_ocean_ecosys_nlat_dimname = '&
//''//mpas_log_escape_dollars(config_global_ocean_ecosys_nlat_dimname))
            call mpas_log_write('        config_global_ocean_ecosys_nlon_dimname = '&
//''//mpas_log_escape_dollars(config_global_ocean_ecosys_nlon_dimname))
            call mpas_log_write('        config_global_ocean_ecosys_ndepth_dimname = '&
//''//mpas_log_escape_dollars(config_global_ocean_ecosys_ndepth_dimname))
            call mpas_log_write('        config_global_ocean_ecosys_depth_conversion_factor = $r', &
realArgs=(/config_global_ocean_ecosys_depth_conversion_factor/))
            call mpas_log_write('        config_global_ocean_ecosys_vert_levels = $i', &
intArgs=(/config_global_ocean_ecosys_vert_levels/))
            call mpas_log_write('        config_global_ocean_ecosys_lat_varname = '&
//''//mpas_log_escape_dollars(config_global_ocean_ecosys_lat_varname))
            call mpas_log_write('        config_global_ocean_ecosys_lon_varname = '&
//''//mpas_log_escape_dollars(config_global_ocean_ecosys_lon_varname))
            call mpas_log_write('        config_global_ocean_ecosys_depth_varname = '&
//''//mpas_log_escape_dollars(config_global_ocean_ecosys_depth_varname))
            call mpas_log_write('        config_global_ocean_ecosys_latlon_degrees = $l', &
logicArgs=(/config_global_ocean_ecosys_latlon_degrees/))
            call mpas_log_write('        config_global_ocean_ecosys_method = '&
//''//mpas_log_escape_dollars(config_global_ocean_ecosys_method))
            call mpas_log_write('        config_global_ocean_ecosys_forcing_time_dimname = '&
//''//mpas_log_escape_dollars(config_global_ocean_ecosys_forcing_time_dimname))
            call mpas_log_write('        config_global_ocean_smooth_ecosys_iterations = $i', &
intArgs=(/config_global_ocean_smooth_ecosys_iterations/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record global_ocean.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_global_ocean_minimum_depth', config_global_ocean_minimum_depth)
      call mpas_pool_add_config(configPool, 'config_global_ocean_depth_file', config_global_ocean_depth_file)
      call mpas_pool_add_config(configPool, 'config_global_ocean_depth_dimname', config_global_ocean_depth_dimname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_depth_varname', config_global_ocean_depth_varname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_depth_conversion_factor', &
config_global_ocean_depth_conversion_factor)
      call mpas_pool_add_config(configPool, 'config_global_ocean_temperature_file', config_global_ocean_temperature_file)
      call mpas_pool_add_config(configPool, 'config_global_ocean_salinity_file', config_global_ocean_salinity_file)
      call mpas_pool_add_config(configPool, 'config_global_ocean_tracer_nlat_dimname', config_global_ocean_tracer_nlat_dimname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_tracer_nlon_dimname', config_global_ocean_tracer_nlon_dimname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_tracer_ndepth_dimname', config_global_ocean_tracer_ndepth_dimname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_tracer_depth_conversion_factor', &
config_global_ocean_tracer_depth_conversion_factor)
      call mpas_pool_add_config(configPool, 'config_global_ocean_tracer_vert_levels', config_global_ocean_tracer_vert_levels)
      call mpas_pool_add_config(configPool, 'config_global_ocean_temperature_varname', config_global_ocean_temperature_varname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_salinity_varname', config_global_ocean_salinity_varname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_tracer_latlon_degrees', config_global_ocean_tracer_latlon_degrees)
      call mpas_pool_add_config(configPool, 'config_global_ocean_tracer_lat_varname', config_global_ocean_tracer_lat_varname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_tracer_lon_varname', config_global_ocean_tracer_lon_varname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_tracer_depth_varname', config_global_ocean_tracer_depth_varname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_tracer_method', config_global_ocean_tracer_method)
      call mpas_pool_add_config(configPool, 'config_global_ocean_smooth_TS_iterations', config_global_ocean_smooth_TS_iterations)
      call mpas_pool_add_config(configPool, 'config_global_ocean_swData_file', config_global_ocean_swData_file)
      call mpas_pool_add_config(configPool, 'config_global_ocean_swData_nlat_dimname', config_global_ocean_swData_nlat_dimname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_swData_nlon_dimname', config_global_ocean_swData_nlon_dimname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_swData_lat_varname', config_global_ocean_swData_lat_varname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_swData_lon_varname', config_global_ocean_swData_lon_varname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_swData_latlon_degrees', config_global_ocean_swData_latlon_degrees)
      call mpas_pool_add_config(configPool, 'config_global_ocean_swData_method', config_global_ocean_swData_method)
      call mpas_pool_add_config(configPool, 'config_global_ocean_chlorophyll_varname', config_global_ocean_chlorophyll_varname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_zenithAngle_varname', config_global_ocean_zenithAngle_varname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_clearSky_varname', config_global_ocean_clearSky_varname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_piston_velocity', config_global_ocean_piston_velocity)
      call mpas_pool_add_config(configPool, 'config_global_ocean_interior_restore_rate', config_global_ocean_interior_restore_rate)
      call mpas_pool_add_config(configPool, 'config_global_ocean_topography_source', config_global_ocean_topography_source)
      call mpas_pool_add_config(configPool, 'config_global_ocean_topography_file', config_global_ocean_topography_file)
      call mpas_pool_add_config(configPool, 'config_global_ocean_topography_nlat_dimname', &
config_global_ocean_topography_nlat_dimname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_topography_nlon_dimname', &
config_global_ocean_topography_nlon_dimname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_topography_latlon_degrees', &
config_global_ocean_topography_latlon_degrees)
      call mpas_pool_add_config(configPool, 'config_global_ocean_topography_lat_varname', &
config_global_ocean_topography_lat_varname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_topography_lon_varname', &
config_global_ocean_topography_lon_varname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_topography_varname', config_global_ocean_topography_varname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_topography_has_ocean_frac', &
config_global_ocean_topography_has_ocean_frac)
      call mpas_pool_add_config(configPool, 'config_global_ocean_topography_ocean_frac_varname', &
config_global_ocean_topography_ocean_frac_varname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_topography_method', config_global_ocean_topography_method)
      call mpas_pool_add_config(configPool, 'config_global_ocean_fill_bathymetry_holes', config_global_ocean_fill_bathymetry_holes)
      call mpas_pool_add_config(configPool, 'config_global_ocean_topography_smooth_iterations', &
config_global_ocean_topography_smooth_iterations)
      call mpas_pool_add_config(configPool, 'config_global_ocean_topography_smooth_weight', &
config_global_ocean_topography_smooth_weight)
      call mpas_pool_add_config(configPool, 'config_global_ocean_deepen_critical_passages', &
config_global_ocean_deepen_critical_passages)
      call mpas_pool_add_config(configPool, 'config_global_ocean_depress_by_land_ice', config_global_ocean_depress_by_land_ice)
      call mpas_pool_add_config(configPool, 'config_global_ocean_land_ice_topo_file', config_global_ocean_land_ice_topo_file)
      call mpas_pool_add_config(configPool, 'config_global_ocean_land_ice_topo_nlat_dimname', &
config_global_ocean_land_ice_topo_nlat_dimname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_land_ice_topo_nlon_dimname', &
config_global_ocean_land_ice_topo_nlon_dimname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_land_ice_topo_latlon_degrees', &
config_global_ocean_land_ice_topo_latlon_degrees)
      call mpas_pool_add_config(configPool, 'config_global_ocean_land_ice_topo_lat_varname', &
config_global_ocean_land_ice_topo_lat_varname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_land_ice_topo_lon_varname', &
config_global_ocean_land_ice_topo_lon_varname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_land_ice_topo_thickness_varname', &
config_global_ocean_land_ice_topo_thickness_varname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_land_ice_topo_draft_varname', &
config_global_ocean_land_ice_topo_draft_varname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_land_ice_topo_ice_frac_varname', &
config_global_ocean_land_ice_topo_ice_frac_varname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_land_ice_topo_grounded_frac_varname', &
config_global_ocean_land_ice_topo_grounded_frac_varname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_use_constant_land_ice_cavity_temperature', &
config_global_ocean_use_constant_land_ice_cavity_temperature)
      call mpas_pool_add_config(configPool, 'config_global_ocean_constant_land_ice_cavity_temperature', &
config_global_ocean_constant_land_ice_cavity_temperature)
      call mpas_pool_add_config(configPool, 'config_global_ocean_cull_inland_seas', config_global_ocean_cull_inland_seas)
      call mpas_pool_add_config(configPool, 'config_global_ocean_windstress_file', config_global_ocean_windstress_file)
      call mpas_pool_add_config(configPool, 'config_global_ocean_windstress_nlat_dimname', &
config_global_ocean_windstress_nlat_dimname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_windstress_nlon_dimname', &
config_global_ocean_windstress_nlon_dimname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_windstress_latlon_degrees', &
config_global_ocean_windstress_latlon_degrees)
      call mpas_pool_add_config(configPool, 'config_global_ocean_windstress_lat_varname', &
config_global_ocean_windstress_lat_varname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_windstress_lon_varname', &
config_global_ocean_windstress_lon_varname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_windstress_zonal_varname', &
config_global_ocean_windstress_zonal_varname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_windstress_meridional_varname', &
config_global_ocean_windstress_meridional_varname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_windstress_method', config_global_ocean_windstress_method)
      call mpas_pool_add_config(configPool, 'config_global_ocean_windstress_conversion_factor', &
config_global_ocean_windstress_conversion_factor)
      call mpas_pool_add_config(configPool, 'config_global_ocean_ecosys_file', config_global_ocean_ecosys_file)
      call mpas_pool_add_config(configPool, 'config_global_ocean_ecosys_forcing_file', config_global_ocean_ecosys_forcing_file)
      call mpas_pool_add_config(configPool, 'config_global_ocean_ecosys_nlat_dimname', config_global_ocean_ecosys_nlat_dimname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_ecosys_nlon_dimname', config_global_ocean_ecosys_nlon_dimname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_ecosys_ndepth_dimname', config_global_ocean_ecosys_ndepth_dimname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_ecosys_depth_conversion_factor', &
config_global_ocean_ecosys_depth_conversion_factor)
      call mpas_pool_add_config(configPool, 'config_global_ocean_ecosys_vert_levels', config_global_ocean_ecosys_vert_levels)
      call mpas_pool_add_config(configPool, 'config_global_ocean_ecosys_lat_varname', config_global_ocean_ecosys_lat_varname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_ecosys_lon_varname', config_global_ocean_ecosys_lon_varname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_ecosys_depth_varname', config_global_ocean_ecosys_depth_varname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_ecosys_latlon_degrees', config_global_ocean_ecosys_latlon_degrees)
      call mpas_pool_add_config(configPool, 'config_global_ocean_ecosys_method', config_global_ocean_ecosys_method)
      call mpas_pool_add_config(configPool, 'config_global_ocean_ecosys_forcing_time_dimname', &
config_global_ocean_ecosys_forcing_time_dimname)
      call mpas_pool_add_config(configPool, 'config_global_ocean_smooth_ecosys_iterations', &
config_global_ocean_smooth_ecosys_iterations)

   end subroutine ocn_setup_nmlrec_global_ocean


   subroutine ocn_setup_nmlrec_cvmix_WSwSBF(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      integer :: config_cvmix_WSwSBF_vert_levels = 100
      real (kind=RKIND) :: config_cvmix_WSwSBF_surface_temperature = 15.000000
      real (kind=RKIND) :: config_cvmix_WSwSBF_surface_salinity = 35.000000
      real (kind=RKIND) :: config_cvmix_WSwSBF_surface_restoring_temperature = 15.000000
      real (kind=RKIND) :: config_cvmix_WSwSBF_surface_restoring_salinity = 35.000000
      real (kind=RKIND) :: config_cvmix_WSwSBF_temperature_piston_velocity = 0.000004
      real (kind=RKIND) :: config_cvmix_WSwSBF_salinity_piston_velocity = 0.000004
      real (kind=RKIND) :: config_cvmix_WSwSBF_sensible_heat_flux = 0.000000
      real (kind=RKIND) :: config_cvmix_WSwSBF_latent_heat_flux = 0.000000
      real (kind=RKIND) :: config_cvmix_WSwSBF_shortwave_heat_flux = 0.000000
      real (kind=RKIND) :: config_cvmix_WSwSBF_rain_flux = 0.000000
      real (kind=RKIND) :: config_cvmix_WSwSBF_evaporation_flux = 0.000000
      real (kind=RKIND) :: config_cvmix_WSwSBF_interior_temperature_restoring_rate = 0.000001
      real (kind=RKIND) :: config_cvmix_WSwSBF_interior_salinity_restoring_rate = 0.000001
      real (kind=RKIND) :: config_cvmix_WSwSBF_temperature_gradient = 0.010000
      real (kind=RKIND) :: config_cvmix_WSwSBF_salinity_gradient = 0.000000
      real (kind=RKIND) :: config_cvmix_WSwSBF_temperature_gradient_mixed_layer = 0.000000
      real (kind=RKIND) :: config_cvmix_WSwSBF_salinity_gradient_mixed_layer = 0.000000
      real (kind=RKIND) :: config_cvmix_WSwSBF_mixed_layer_depth_temperature = 0.000000
      real (kind=RKIND) :: config_cvmix_WSwSBF_mixed_layer_depth_salinity = 0.000000
      real (kind=RKIND) :: config_cvmix_WSwSBF_mixed_layer_temperature_change = 0.000000
      real (kind=RKIND) :: config_cvmix_WSwSBF_mixed_layer_salinity_change = 0.000000
      character (len=StrKIND) :: config_cvmix_WSwSBF_vertical_grid = 'uniform'
      real (kind=RKIND) :: config_cvmix_WSwSBF_bottom_depth = 400.000000
      real (kind=RKIND) :: config_cvmix_WSwSBF_max_windstress = 0.100000
      real (kind=RKIND) :: config_cvmix_WSwSBF_coriolis_parameter = 0.000100

      namelist /cvmix_WSwSBF/ &
         config_cvmix_WSwSBF_vert_levels, &
         config_cvmix_WSwSBF_surface_temperature, &
         config_cvmix_WSwSBF_surface_salinity, &
         config_cvmix_WSwSBF_surface_restoring_temperature, &
         config_cvmix_WSwSBF_surface_restoring_salinity, &
         config_cvmix_WSwSBF_temperature_piston_velocity, &
         config_cvmix_WSwSBF_salinity_piston_velocity, &
         config_cvmix_WSwSBF_sensible_heat_flux, &
         config_cvmix_WSwSBF_latent_heat_flux, &
         config_cvmix_WSwSBF_shortwave_heat_flux, &
         config_cvmix_WSwSBF_rain_flux, &
         config_cvmix_WSwSBF_evaporation_flux, &
         config_cvmix_WSwSBF_interior_temperature_restoring_rate, &
         config_cvmix_WSwSBF_interior_salinity_restoring_rate, &
         config_cvmix_WSwSBF_temperature_gradient, &
         config_cvmix_WSwSBF_salinity_gradient, &
         config_cvmix_WSwSBF_temperature_gradient_mixed_layer, &
         config_cvmix_WSwSBF_salinity_gradient_mixed_layer, &
         config_cvmix_WSwSBF_mixed_layer_depth_temperature, &
         config_cvmix_WSwSBF_mixed_layer_depth_salinity, &
         config_cvmix_WSwSBF_mixed_layer_temperature_change, &
         config_cvmix_WSwSBF_mixed_layer_salinity_change, &
         config_cvmix_WSwSBF_vertical_grid, &
         config_cvmix_WSwSBF_bottom_depth, &
         config_cvmix_WSwSBF_max_windstress, &
         config_cvmix_WSwSBF_coriolis_parameter
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, cvmix_WSwSBF, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_int(dminfo, config_cvmix_WSwSBF_vert_levels)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_WSwSBF_surface_temperature)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_WSwSBF_surface_salinity)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_WSwSBF_surface_restoring_temperature)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_WSwSBF_surface_restoring_salinity)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_WSwSBF_temperature_piston_velocity)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_WSwSBF_salinity_piston_velocity)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_WSwSBF_sensible_heat_flux)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_WSwSBF_latent_heat_flux)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_WSwSBF_shortwave_heat_flux)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_WSwSBF_rain_flux)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_WSwSBF_evaporation_flux)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_WSwSBF_interior_temperature_restoring_rate)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_WSwSBF_interior_salinity_restoring_rate)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_WSwSBF_temperature_gradient)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_WSwSBF_salinity_gradient)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_WSwSBF_temperature_gradient_mixed_layer)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_WSwSBF_salinity_gradient_mixed_layer)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_WSwSBF_mixed_layer_depth_temperature)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_WSwSBF_mixed_layer_depth_salinity)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_WSwSBF_mixed_layer_temperature_change)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_WSwSBF_mixed_layer_salinity_change)
         call mpas_dmpar_bcast_char(dminfo, config_cvmix_WSwSBF_vertical_grid)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_WSwSBF_bottom_depth)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_WSwSBF_max_windstress)
         call mpas_dmpar_bcast_real(dminfo, config_cvmix_WSwSBF_coriolis_parameter)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record cvmix_WSwSBF')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_cvmix_WSwSBF_vert_levels = $i', intArgs=(/config_cvmix_WSwSBF_vert_levels/))
            call mpas_log_write('        config_cvmix_WSwSBF_surface_temperature = $r', &
realArgs=(/config_cvmix_WSwSBF_surface_temperature/))
            call mpas_log_write('        config_cvmix_WSwSBF_surface_salinity = $r', &
realArgs=(/config_cvmix_WSwSBF_surface_salinity/))
            call mpas_log_write('        config_cvmix_WSwSBF_surface_restoring_temperature = $r', &
realArgs=(/config_cvmix_WSwSBF_surface_restoring_temperature/))
            call mpas_log_write('        config_cvmix_WSwSBF_surface_restoring_salinity = $r', &
realArgs=(/config_cvmix_WSwSBF_surface_restoring_salinity/))
            call mpas_log_write('        config_cvmix_WSwSBF_temperature_piston_velocity = $r', &
realArgs=(/config_cvmix_WSwSBF_temperature_piston_velocity/))
            call mpas_log_write('        config_cvmix_WSwSBF_salinity_piston_velocity = $r', &
realArgs=(/config_cvmix_WSwSBF_salinity_piston_velocity/))
            call mpas_log_write('        config_cvmix_WSwSBF_sensible_heat_flux = $r', &
realArgs=(/config_cvmix_WSwSBF_sensible_heat_flux/))
            call mpas_log_write('        config_cvmix_WSwSBF_latent_heat_flux = $r', &
realArgs=(/config_cvmix_WSwSBF_latent_heat_flux/))
            call mpas_log_write('        config_cvmix_WSwSBF_shortwave_heat_flux = $r', &
realArgs=(/config_cvmix_WSwSBF_shortwave_heat_flux/))
            call mpas_log_write('        config_cvmix_WSwSBF_rain_flux = $r', realArgs=(/config_cvmix_WSwSBF_rain_flux/))
            call mpas_log_write('        config_cvmix_WSwSBF_evaporation_flux = $r', &
realArgs=(/config_cvmix_WSwSBF_evaporation_flux/))
            call mpas_log_write('        config_cvmix_WSwSBF_interior_temperature_restoring_rate = $r', &
realArgs=(/config_cvmix_WSwSBF_interior_temperature_restoring_rate/))
            call mpas_log_write('        config_cvmix_WSwSBF_interior_salinity_restoring_rate = $r', &
realArgs=(/config_cvmix_WSwSBF_interior_salinity_restoring_rate/))
            call mpas_log_write('        config_cvmix_WSwSBF_temperature_gradient = $r', &
realArgs=(/config_cvmix_WSwSBF_temperature_gradient/))
            call mpas_log_write('        config_cvmix_WSwSBF_salinity_gradient = $r', &
realArgs=(/config_cvmix_WSwSBF_salinity_gradient/))
            call mpas_log_write('        config_cvmix_WSwSBF_temperature_gradient_mixed_layer = $r', &
realArgs=(/config_cvmix_WSwSBF_temperature_gradient_mixed_layer/))
            call mpas_log_write('        config_cvmix_WSwSBF_salinity_gradient_mixed_layer = $r', &
realArgs=(/config_cvmix_WSwSBF_salinity_gradient_mixed_layer/))
            call mpas_log_write('        config_cvmix_WSwSBF_mixed_layer_depth_temperature = $r', &
realArgs=(/config_cvmix_WSwSBF_mixed_layer_depth_temperature/))
            call mpas_log_write('        config_cvmix_WSwSBF_mixed_layer_depth_salinity = $r', &
realArgs=(/config_cvmix_WSwSBF_mixed_layer_depth_salinity/))
            call mpas_log_write('        config_cvmix_WSwSBF_mixed_layer_temperature_change = $r', &
realArgs=(/config_cvmix_WSwSBF_mixed_layer_temperature_change/))
            call mpas_log_write('        config_cvmix_WSwSBF_mixed_layer_salinity_change = $r', &
realArgs=(/config_cvmix_WSwSBF_mixed_layer_salinity_change/))
            call mpas_log_write('        config_cvmix_WSwSBF_vertical_grid = '&
//''//mpas_log_escape_dollars(config_cvmix_WSwSBF_vertical_grid))
            call mpas_log_write('        config_cvmix_WSwSBF_bottom_depth = $r', realArgs=(/config_cvmix_WSwSBF_bottom_depth/))
            call mpas_log_write('        config_cvmix_WSwSBF_max_windstress = $r', realArgs=(/config_cvmix_WSwSBF_max_windstress/))
            call mpas_log_write('        config_cvmix_WSwSBF_coriolis_parameter = $r', &
realArgs=(/config_cvmix_WSwSBF_coriolis_parameter/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record cvmix_WSwSBF.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_cvmix_WSwSBF_vert_levels', config_cvmix_WSwSBF_vert_levels)
      call mpas_pool_add_config(configPool, 'config_cvmix_WSwSBF_surface_temperature', config_cvmix_WSwSBF_surface_temperature)
      call mpas_pool_add_config(configPool, 'config_cvmix_WSwSBF_surface_salinity', config_cvmix_WSwSBF_surface_salinity)
      call mpas_pool_add_config(configPool, 'config_cvmix_WSwSBF_surface_restoring_temperature', &
config_cvmix_WSwSBF_surface_restoring_temperature)
      call mpas_pool_add_config(configPool, 'config_cvmix_WSwSBF_surface_restoring_salinity', &
config_cvmix_WSwSBF_surface_restoring_salinity)
      call mpas_pool_add_config(configPool, 'config_cvmix_WSwSBF_temperature_piston_velocity', &
config_cvmix_WSwSBF_temperature_piston_velocity)
      call mpas_pool_add_config(configPool, 'config_cvmix_WSwSBF_salinity_piston_velocity', &
config_cvmix_WSwSBF_salinity_piston_velocity)
      call mpas_pool_add_config(configPool, 'config_cvmix_WSwSBF_sensible_heat_flux', config_cvmix_WSwSBF_sensible_heat_flux)
      call mpas_pool_add_config(configPool, 'config_cvmix_WSwSBF_latent_heat_flux', config_cvmix_WSwSBF_latent_heat_flux)
      call mpas_pool_add_config(configPool, 'config_cvmix_WSwSBF_shortwave_heat_flux', config_cvmix_WSwSBF_shortwave_heat_flux)
      call mpas_pool_add_config(configPool, 'config_cvmix_WSwSBF_rain_flux', config_cvmix_WSwSBF_rain_flux)
      call mpas_pool_add_config(configPool, 'config_cvmix_WSwSBF_evaporation_flux', config_cvmix_WSwSBF_evaporation_flux)
      call mpas_pool_add_config(configPool, 'config_cvmix_WSwSBF_interior_temperature_restoring_rate', &
config_cvmix_WSwSBF_interior_temperature_restoring_rate)
      call mpas_pool_add_config(configPool, 'config_cvmix_WSwSBF_interior_salinity_restoring_rate', &
config_cvmix_WSwSBF_interior_salinity_restoring_rate)
      call mpas_pool_add_config(configPool, 'config_cvmix_WSwSBF_temperature_gradient', config_cvmix_WSwSBF_temperature_gradient)
      call mpas_pool_add_config(configPool, 'config_cvmix_WSwSBF_salinity_gradient', config_cvmix_WSwSBF_salinity_gradient)
      call mpas_pool_add_config(configPool, 'config_cvmix_WSwSBF_temperature_gradient_mixed_layer', &
config_cvmix_WSwSBF_temperature_gradient_mixed_layer)
      call mpas_pool_add_config(configPool, 'config_cvmix_WSwSBF_salinity_gradient_mixed_layer', &
config_cvmix_WSwSBF_salinity_gradient_mixed_layer)
      call mpas_pool_add_config(configPool, 'config_cvmix_WSwSBF_mixed_layer_depth_temperature', &
config_cvmix_WSwSBF_mixed_layer_depth_temperature)
      call mpas_pool_add_config(configPool, 'config_cvmix_WSwSBF_mixed_layer_depth_salinity', &
config_cvmix_WSwSBF_mixed_layer_depth_salinity)
      call mpas_pool_add_config(configPool, 'config_cvmix_WSwSBF_mixed_layer_temperature_change', &
config_cvmix_WSwSBF_mixed_layer_temperature_change)
      call mpas_pool_add_config(configPool, 'config_cvmix_WSwSBF_mixed_layer_salinity_change', &
config_cvmix_WSwSBF_mixed_layer_salinity_change)
      call mpas_pool_add_config(configPool, 'config_cvmix_WSwSBF_vertical_grid', config_cvmix_WSwSBF_vertical_grid)
      call mpas_pool_add_config(configPool, 'config_cvmix_WSwSBF_bottom_depth', config_cvmix_WSwSBF_bottom_depth)
      call mpas_pool_add_config(configPool, 'config_cvmix_WSwSBF_max_windstress', config_cvmix_WSwSBF_max_windstress)
      call mpas_pool_add_config(configPool, 'config_cvmix_WSwSBF_coriolis_parameter', config_cvmix_WSwSBF_coriolis_parameter)

   end subroutine ocn_setup_nmlrec_cvmix_WSwSBF


   subroutine ocn_setup_nmlrec_iso(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      integer :: config_iso_vert_levels = 100
      real (kind=RKIND) :: config_iso_main_channel_depth = 4000.000000
      real (kind=RKIND) :: config_iso_north_wall_lat = -50.000000
      real (kind=RKIND) :: config_iso_south_wall_lat = -70.000000
      logical :: config_iso_ridge_flag = .true.
      real (kind=RKIND) :: config_iso_ridge_center_lon = 180.000000
      real (kind=RKIND) :: config_iso_ridge_height = 2000.000000
      real (kind=RKIND) :: config_iso_ridge_width = 2000000.000000
      logical :: config_iso_plateau_flag = .true.
      real (kind=RKIND) :: config_iso_plateau_center_lon = 300.000000
      real (kind=RKIND) :: config_iso_plateau_center_lat = -58.000000
      real (kind=RKIND) :: config_iso_plateau_height = 2000.000000
      real (kind=RKIND) :: config_iso_plateau_radius = 200000.000000
      real (kind=RKIND) :: config_iso_plateau_slope_width = 1000000.000000
      logical :: config_iso_shelf_flag = .true.
      real (kind=RKIND) :: config_iso_shelf_depth = 500.000000
      real (kind=RKIND) :: config_iso_shelf_width = 120000.000000
      logical :: config_iso_cont_slope_flag = .true.
      real (kind=RKIND) :: config_iso_max_cont_slope = 0.010000
      logical :: config_iso_embayment_flag = .true.
      real (kind=RKIND) :: config_iso_embayment_center_lon = 60.000000
      real (kind=RKIND) :: config_iso_embayment_center_lat = -71.000000
      real (kind=RKIND) :: config_iso_embayment_radius = 500000.000000
      real (kind=RKIND) :: config_iso_embayment_depth = 2000.000000
      logical :: config_iso_depression_flag = .true.
      real (kind=RKIND) :: config_iso_depression_center_lon = 60.000000
      real (kind=RKIND) :: config_iso_depression_south_lat = -72.000000
      real (kind=RKIND) :: config_iso_depression_north_lat = -65.000000
      real (kind=RKIND) :: config_iso_depression_width = 480000.000000
      real (kind=RKIND) :: config_iso_depression_depth = 800.000000
      real (kind=RKIND) :: config_iso_salinity = 35.000000
      real (kind=RKIND) :: config_iso_wind_stress_max = 0.010000
      real (kind=RKIND) :: config_iso_acc_wind = 0.200000
      real (kind=RKIND) :: config_iso_asf_wind = -0.050000
      real (kind=RKIND) :: config_iso_wind_trans = -65.000000
      real (kind=RKIND) :: config_iso_heat_flux_south = -5.000000
      real (kind=RKIND) :: config_iso_heat_flux_middle = 10.000000
      real (kind=RKIND) :: config_iso_heat_flux_north = -5.000000
      real (kind=RKIND) :: config_iso_heat_flux_lat_ss = -70.000000
      real (kind=RKIND) :: config_iso_heat_flux_lat_sm = -65.000000
      real (kind=RKIND) :: config_iso_heat_flux_lat_mn = -53.000000
      real (kind=RKIND) :: config_iso_region1_center_lon = 60.000000
      real (kind=RKIND) :: config_iso_region1_center_lat = -75.000000
      real (kind=RKIND) :: config_iso_region2_center_lon = 150.000000
      real (kind=RKIND) :: config_iso_region2_center_lat = -71.000000
      real (kind=RKIND) :: config_iso_region3_center_lon = 240.000000
      real (kind=RKIND) :: config_iso_region3_center_lat = -71.000000
      real (kind=RKIND) :: config_iso_region4_center_lon = 330.000000
      real (kind=RKIND) :: config_iso_region4_center_lat = -71.000000
      logical :: config_iso_heat_flux_region1_flag = .false.
      real (kind=RKIND) :: config_iso_heat_flux_region1 = -5.000000
      real (kind=RKIND) :: config_iso_heat_flux_region1_radius = 300000.000000
      logical :: config_iso_heat_flux_region2_flag = .false.
      real (kind=RKIND) :: config_iso_heat_flux_region2 = -5.000000
      real (kind=RKIND) :: config_iso_heat_flux_region2_radius = 240000.000000
      real (kind=RKIND) :: config_iso_surface_temperature_piston_velocity = 0.000058
      real (kind=RKIND) :: config_iso_initial_temp_t1 = 3.500000
      real (kind=RKIND) :: config_iso_initial_temp_t2 = 4.000000
      real (kind=RKIND) :: config_iso_initial_temp_h0 = 1200.000000
      real (kind=RKIND) :: config_iso_initial_temp_h1 = 500.000000
      real (kind=RKIND) :: config_iso_initial_temp_mt = 0.000075
      real (kind=RKIND) :: config_iso_initial_temp_latS = -75.000000
      real (kind=RKIND) :: config_iso_initial_temp_latN = -50.000000
      real (kind=RKIND) :: config_iso_temperature_sponge_t1 = 10.000000
      real (kind=RKIND) :: config_iso_temperature_sponge_h1 = 1000.000000
      real (kind=RKIND) :: config_iso_temperature_sponge_l1 = 120000.000000
      real (kind=RKIND) :: config_iso_temperature_sponge_tau1 = 10.000000
      logical :: config_iso_temperature_restore_region1_flag = .true.
      real (kind=RKIND) :: config_iso_temperature_restore_t1 = -1.000000
      real (kind=RKIND) :: config_iso_temperature_restore_lcx1 = 600000.000000
      real (kind=RKIND) :: config_iso_temperature_restore_lcy1 = 600000.000000
      logical :: config_iso_temperature_restore_region2_flag = .true.
      real (kind=RKIND) :: config_iso_temperature_restore_t2 = -1.000000
      real (kind=RKIND) :: config_iso_temperature_restore_lcx2 = 600000.000000
      real (kind=RKIND) :: config_iso_temperature_restore_lcy2 = 250000.000000
      logical :: config_iso_temperature_restore_region3_flag = .true.
      real (kind=RKIND) :: config_iso_temperature_restore_t3 = -1.000000
      real (kind=RKIND) :: config_iso_temperature_restore_lcx3 = 600000.000000
      real (kind=RKIND) :: config_iso_temperature_restore_lcy3 = 250000.000000
      logical :: config_iso_temperature_restore_region4_flag = .true.
      real (kind=RKIND) :: config_iso_temperature_restore_t4 = -1.000000
      real (kind=RKIND) :: config_iso_temperature_restore_lcx4 = 600000.000000
      real (kind=RKIND) :: config_iso_temperature_restore_lcy4 = 250000.000000

      namelist /iso/ &
         config_iso_vert_levels, &
         config_iso_main_channel_depth, &
         config_iso_north_wall_lat, &
         config_iso_south_wall_lat, &
         config_iso_ridge_flag, &
         config_iso_ridge_center_lon, &
         config_iso_ridge_height, &
         config_iso_ridge_width, &
         config_iso_plateau_flag, &
         config_iso_plateau_center_lon, &
         config_iso_plateau_center_lat, &
         config_iso_plateau_height, &
         config_iso_plateau_radius, &
         config_iso_plateau_slope_width, &
         config_iso_shelf_flag, &
         config_iso_shelf_depth, &
         config_iso_shelf_width, &
         config_iso_cont_slope_flag, &
         config_iso_max_cont_slope, &
         config_iso_embayment_flag, &
         config_iso_embayment_center_lon, &
         config_iso_embayment_center_lat, &
         config_iso_embayment_radius, &
         config_iso_embayment_depth, &
         config_iso_depression_flag, &
         config_iso_depression_center_lon, &
         config_iso_depression_south_lat, &
         config_iso_depression_north_lat, &
         config_iso_depression_width, &
         config_iso_depression_depth, &
         config_iso_salinity, &
         config_iso_wind_stress_max, &
         config_iso_acc_wind, &
         config_iso_asf_wind, &
         config_iso_wind_trans, &
         config_iso_heat_flux_south, &
         config_iso_heat_flux_middle, &
         config_iso_heat_flux_north, &
         config_iso_heat_flux_lat_ss, &
         config_iso_heat_flux_lat_sm, &
         config_iso_heat_flux_lat_mn, &
         config_iso_region1_center_lon, &
         config_iso_region1_center_lat, &
         config_iso_region2_center_lon, &
         config_iso_region2_center_lat, &
         config_iso_region3_center_lon, &
         config_iso_region3_center_lat, &
         config_iso_region4_center_lon, &
         config_iso_region4_center_lat, &
         config_iso_heat_flux_region1_flag, &
         config_iso_heat_flux_region1, &
         config_iso_heat_flux_region1_radius, &
         config_iso_heat_flux_region2_flag, &
         config_iso_heat_flux_region2, &
         config_iso_heat_flux_region2_radius, &
         config_iso_surface_temperature_piston_velocity, &
         config_iso_initial_temp_t1, &
         config_iso_initial_temp_t2, &
         config_iso_initial_temp_h0, &
         config_iso_initial_temp_h1, &
         config_iso_initial_temp_mt, &
         config_iso_initial_temp_latS, &
         config_iso_initial_temp_latN, &
         config_iso_temperature_sponge_t1, &
         config_iso_temperature_sponge_h1, &
         config_iso_temperature_sponge_l1, &
         config_iso_temperature_sponge_tau1, &
         config_iso_temperature_restore_region1_flag, &
         config_iso_temperature_restore_t1, &
         config_iso_temperature_restore_lcx1, &
         config_iso_temperature_restore_lcy1, &
         config_iso_temperature_restore_region2_flag, &
         config_iso_temperature_restore_t2, &
         config_iso_temperature_restore_lcx2, &
         config_iso_temperature_restore_lcy2, &
         config_iso_temperature_restore_region3_flag, &
         config_iso_temperature_restore_t3, &
         config_iso_temperature_restore_lcx3, &
         config_iso_temperature_restore_lcy3, &
         config_iso_temperature_restore_region4_flag, &
         config_iso_temperature_restore_t4, &
         config_iso_temperature_restore_lcx4, &
         config_iso_temperature_restore_lcy4
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, iso, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_int(dminfo, config_iso_vert_levels)
         call mpas_dmpar_bcast_real(dminfo, config_iso_main_channel_depth)
         call mpas_dmpar_bcast_real(dminfo, config_iso_north_wall_lat)
         call mpas_dmpar_bcast_real(dminfo, config_iso_south_wall_lat)
         call mpas_dmpar_bcast_logical(dminfo, config_iso_ridge_flag)
         call mpas_dmpar_bcast_real(dminfo, config_iso_ridge_center_lon)
         call mpas_dmpar_bcast_real(dminfo, config_iso_ridge_height)
         call mpas_dmpar_bcast_real(dminfo, config_iso_ridge_width)
         call mpas_dmpar_bcast_logical(dminfo, config_iso_plateau_flag)
         call mpas_dmpar_bcast_real(dminfo, config_iso_plateau_center_lon)
         call mpas_dmpar_bcast_real(dminfo, config_iso_plateau_center_lat)
         call mpas_dmpar_bcast_real(dminfo, config_iso_plateau_height)
         call mpas_dmpar_bcast_real(dminfo, config_iso_plateau_radius)
         call mpas_dmpar_bcast_real(dminfo, config_iso_plateau_slope_width)
         call mpas_dmpar_bcast_logical(dminfo, config_iso_shelf_flag)
         call mpas_dmpar_bcast_real(dminfo, config_iso_shelf_depth)
         call mpas_dmpar_bcast_real(dminfo, config_iso_shelf_width)
         call mpas_dmpar_bcast_logical(dminfo, config_iso_cont_slope_flag)
         call mpas_dmpar_bcast_real(dminfo, config_iso_max_cont_slope)
         call mpas_dmpar_bcast_logical(dminfo, config_iso_embayment_flag)
         call mpas_dmpar_bcast_real(dminfo, config_iso_embayment_center_lon)
         call mpas_dmpar_bcast_real(dminfo, config_iso_embayment_center_lat)
         call mpas_dmpar_bcast_real(dminfo, config_iso_embayment_radius)
         call mpas_dmpar_bcast_real(dminfo, config_iso_embayment_depth)
         call mpas_dmpar_bcast_logical(dminfo, config_iso_depression_flag)
         call mpas_dmpar_bcast_real(dminfo, config_iso_depression_center_lon)
         call mpas_dmpar_bcast_real(dminfo, config_iso_depression_south_lat)
         call mpas_dmpar_bcast_real(dminfo, config_iso_depression_north_lat)
         call mpas_dmpar_bcast_real(dminfo, config_iso_depression_width)
         call mpas_dmpar_bcast_real(dminfo, config_iso_depression_depth)
         call mpas_dmpar_bcast_real(dminfo, config_iso_salinity)
         call mpas_dmpar_bcast_real(dminfo, config_iso_wind_stress_max)
         call mpas_dmpar_bcast_real(dminfo, config_iso_acc_wind)
         call mpas_dmpar_bcast_real(dminfo, config_iso_asf_wind)
         call mpas_dmpar_bcast_real(dminfo, config_iso_wind_trans)
         call mpas_dmpar_bcast_real(dminfo, config_iso_heat_flux_south)
         call mpas_dmpar_bcast_real(dminfo, config_iso_heat_flux_middle)
         call mpas_dmpar_bcast_real(dminfo, config_iso_heat_flux_north)
         call mpas_dmpar_bcast_real(dminfo, config_iso_heat_flux_lat_ss)
         call mpas_dmpar_bcast_real(dminfo, config_iso_heat_flux_lat_sm)
         call mpas_dmpar_bcast_real(dminfo, config_iso_heat_flux_lat_mn)
         call mpas_dmpar_bcast_real(dminfo, config_iso_region1_center_lon)
         call mpas_dmpar_bcast_real(dminfo, config_iso_region1_center_lat)
         call mpas_dmpar_bcast_real(dminfo, config_iso_region2_center_lon)
         call mpas_dmpar_bcast_real(dminfo, config_iso_region2_center_lat)
         call mpas_dmpar_bcast_real(dminfo, config_iso_region3_center_lon)
         call mpas_dmpar_bcast_real(dminfo, config_iso_region3_center_lat)
         call mpas_dmpar_bcast_real(dminfo, config_iso_region4_center_lon)
         call mpas_dmpar_bcast_real(dminfo, config_iso_region4_center_lat)
         call mpas_dmpar_bcast_logical(dminfo, config_iso_heat_flux_region1_flag)
         call mpas_dmpar_bcast_real(dminfo, config_iso_heat_flux_region1)
         call mpas_dmpar_bcast_real(dminfo, config_iso_heat_flux_region1_radius)
         call mpas_dmpar_bcast_logical(dminfo, config_iso_heat_flux_region2_flag)
         call mpas_dmpar_bcast_real(dminfo, config_iso_heat_flux_region2)
         call mpas_dmpar_bcast_real(dminfo, config_iso_heat_flux_region2_radius)
         call mpas_dmpar_bcast_real(dminfo, config_iso_surface_temperature_piston_velocity)
         call mpas_dmpar_bcast_real(dminfo, config_iso_initial_temp_t1)
         call mpas_dmpar_bcast_real(dminfo, config_iso_initial_temp_t2)
         call mpas_dmpar_bcast_real(dminfo, config_iso_initial_temp_h0)
         call mpas_dmpar_bcast_real(dminfo, config_iso_initial_temp_h1)
         call mpas_dmpar_bcast_real(dminfo, config_iso_initial_temp_mt)
         call mpas_dmpar_bcast_real(dminfo, config_iso_initial_temp_latS)
         call mpas_dmpar_bcast_real(dminfo, config_iso_initial_temp_latN)
         call mpas_dmpar_bcast_real(dminfo, config_iso_temperature_sponge_t1)
         call mpas_dmpar_bcast_real(dminfo, config_iso_temperature_sponge_h1)
         call mpas_dmpar_bcast_real(dminfo, config_iso_temperature_sponge_l1)
         call mpas_dmpar_bcast_real(dminfo, config_iso_temperature_sponge_tau1)
         call mpas_dmpar_bcast_logical(dminfo, config_iso_temperature_restore_region1_flag)
         call mpas_dmpar_bcast_real(dminfo, config_iso_temperature_restore_t1)
         call mpas_dmpar_bcast_real(dminfo, config_iso_temperature_restore_lcx1)
         call mpas_dmpar_bcast_real(dminfo, config_iso_temperature_restore_lcy1)
         call mpas_dmpar_bcast_logical(dminfo, config_iso_temperature_restore_region2_flag)
         call mpas_dmpar_bcast_real(dminfo, config_iso_temperature_restore_t2)
         call mpas_dmpar_bcast_real(dminfo, config_iso_temperature_restore_lcx2)
         call mpas_dmpar_bcast_real(dminfo, config_iso_temperature_restore_lcy2)
         call mpas_dmpar_bcast_logical(dminfo, config_iso_temperature_restore_region3_flag)
         call mpas_dmpar_bcast_real(dminfo, config_iso_temperature_restore_t3)
         call mpas_dmpar_bcast_real(dminfo, config_iso_temperature_restore_lcx3)
         call mpas_dmpar_bcast_real(dminfo, config_iso_temperature_restore_lcy3)
         call mpas_dmpar_bcast_logical(dminfo, config_iso_temperature_restore_region4_flag)
         call mpas_dmpar_bcast_real(dminfo, config_iso_temperature_restore_t4)
         call mpas_dmpar_bcast_real(dminfo, config_iso_temperature_restore_lcx4)
         call mpas_dmpar_bcast_real(dminfo, config_iso_temperature_restore_lcy4)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record iso')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_iso_vert_levels = $i', intArgs=(/config_iso_vert_levels/))
            call mpas_log_write('        config_iso_main_channel_depth = $r', realArgs=(/config_iso_main_channel_depth/))
            call mpas_log_write('        config_iso_north_wall_lat = $r', realArgs=(/config_iso_north_wall_lat/))
            call mpas_log_write('        config_iso_south_wall_lat = $r', realArgs=(/config_iso_south_wall_lat/))
            call mpas_log_write('        config_iso_ridge_flag = $l', logicArgs=(/config_iso_ridge_flag/))
            call mpas_log_write('        config_iso_ridge_center_lon = $r', realArgs=(/config_iso_ridge_center_lon/))
            call mpas_log_write('        config_iso_ridge_height = $r', realArgs=(/config_iso_ridge_height/))
            call mpas_log_write('        config_iso_ridge_width = $r', realArgs=(/config_iso_ridge_width/))
            call mpas_log_write('        config_iso_plateau_flag = $l', logicArgs=(/config_iso_plateau_flag/))
            call mpas_log_write('        config_iso_plateau_center_lon = $r', realArgs=(/config_iso_plateau_center_lon/))
            call mpas_log_write('        config_iso_plateau_center_lat = $r', realArgs=(/config_iso_plateau_center_lat/))
            call mpas_log_write('        config_iso_plateau_height = $r', realArgs=(/config_iso_plateau_height/))
            call mpas_log_write('        config_iso_plateau_radius = $r', realArgs=(/config_iso_plateau_radius/))
            call mpas_log_write('        config_iso_plateau_slope_width = $r', realArgs=(/config_iso_plateau_slope_width/))
            call mpas_log_write('        config_iso_shelf_flag = $l', logicArgs=(/config_iso_shelf_flag/))
            call mpas_log_write('        config_iso_shelf_depth = $r', realArgs=(/config_iso_shelf_depth/))
            call mpas_log_write('        config_iso_shelf_width = $r', realArgs=(/config_iso_shelf_width/))
            call mpas_log_write('        config_iso_cont_slope_flag = $l', logicArgs=(/config_iso_cont_slope_flag/))
            call mpas_log_write('        config_iso_max_cont_slope = $r', realArgs=(/config_iso_max_cont_slope/))
            call mpas_log_write('        config_iso_embayment_flag = $l', logicArgs=(/config_iso_embayment_flag/))
            call mpas_log_write('        config_iso_embayment_center_lon = $r', realArgs=(/config_iso_embayment_center_lon/))
            call mpas_log_write('        config_iso_embayment_center_lat = $r', realArgs=(/config_iso_embayment_center_lat/))
            call mpas_log_write('        config_iso_embayment_radius = $r', realArgs=(/config_iso_embayment_radius/))
            call mpas_log_write('        config_iso_embayment_depth = $r', realArgs=(/config_iso_embayment_depth/))
            call mpas_log_write('        config_iso_depression_flag = $l', logicArgs=(/config_iso_depression_flag/))
            call mpas_log_write('        config_iso_depression_center_lon = $r', realArgs=(/config_iso_depression_center_lon/))
            call mpas_log_write('        config_iso_depression_south_lat = $r', realArgs=(/config_iso_depression_south_lat/))
            call mpas_log_write('        config_iso_depression_north_lat = $r', realArgs=(/config_iso_depression_north_lat/))
            call mpas_log_write('        config_iso_depression_width = $r', realArgs=(/config_iso_depression_width/))
            call mpas_log_write('        config_iso_depression_depth = $r', realArgs=(/config_iso_depression_depth/))
            call mpas_log_write('        config_iso_salinity = $r', realArgs=(/config_iso_salinity/))
            call mpas_log_write('        config_iso_wind_stress_max = $r', realArgs=(/config_iso_wind_stress_max/))
            call mpas_log_write('        config_iso_acc_wind = $r', realArgs=(/config_iso_acc_wind/))
            call mpas_log_write('        config_iso_asf_wind = $r', realArgs=(/config_iso_asf_wind/))
            call mpas_log_write('        config_iso_wind_trans = $r', realArgs=(/config_iso_wind_trans/))
            call mpas_log_write('        config_iso_heat_flux_south = $r', realArgs=(/config_iso_heat_flux_south/))
            call mpas_log_write('        config_iso_heat_flux_middle = $r', realArgs=(/config_iso_heat_flux_middle/))
            call mpas_log_write('        config_iso_heat_flux_north = $r', realArgs=(/config_iso_heat_flux_north/))
            call mpas_log_write('        config_iso_heat_flux_lat_ss = $r', realArgs=(/config_iso_heat_flux_lat_ss/))
            call mpas_log_write('        config_iso_heat_flux_lat_sm = $r', realArgs=(/config_iso_heat_flux_lat_sm/))
            call mpas_log_write('        config_iso_heat_flux_lat_mn = $r', realArgs=(/config_iso_heat_flux_lat_mn/))
            call mpas_log_write('        config_iso_region1_center_lon = $r', realArgs=(/config_iso_region1_center_lon/))
            call mpas_log_write('        config_iso_region1_center_lat = $r', realArgs=(/config_iso_region1_center_lat/))
            call mpas_log_write('        config_iso_region2_center_lon = $r', realArgs=(/config_iso_region2_center_lon/))
            call mpas_log_write('        config_iso_region2_center_lat = $r', realArgs=(/config_iso_region2_center_lat/))
            call mpas_log_write('        config_iso_region3_center_lon = $r', realArgs=(/config_iso_region3_center_lon/))
            call mpas_log_write('        config_iso_region3_center_lat = $r', realArgs=(/config_iso_region3_center_lat/))
            call mpas_log_write('        config_iso_region4_center_lon = $r', realArgs=(/config_iso_region4_center_lon/))
            call mpas_log_write('        config_iso_region4_center_lat = $r', realArgs=(/config_iso_region4_center_lat/))
            call mpas_log_write('        config_iso_heat_flux_region1_flag = $l', logicArgs=(/config_iso_heat_flux_region1_flag/))
            call mpas_log_write('        config_iso_heat_flux_region1 = $r', realArgs=(/config_iso_heat_flux_region1/))
            call mpas_log_write('        config_iso_heat_flux_region1_radius = $r', &
realArgs=(/config_iso_heat_flux_region1_radius/))
            call mpas_log_write('        config_iso_heat_flux_region2_flag = $l', logicArgs=(/config_iso_heat_flux_region2_flag/))
            call mpas_log_write('        config_iso_heat_flux_region2 = $r', realArgs=(/config_iso_heat_flux_region2/))
            call mpas_log_write('        config_iso_heat_flux_region2_radius = $r', &
realArgs=(/config_iso_heat_flux_region2_radius/))
            call mpas_log_write('        config_iso_surface_temperature_piston_velocity = $r', &
realArgs=(/config_iso_surface_temperature_piston_velocity/))
            call mpas_log_write('        config_iso_initial_temp_t1 = $r', realArgs=(/config_iso_initial_temp_t1/))
            call mpas_log_write('        config_iso_initial_temp_t2 = $r', realArgs=(/config_iso_initial_temp_t2/))
            call mpas_log_write('        config_iso_initial_temp_h0 = $r', realArgs=(/config_iso_initial_temp_h0/))
            call mpas_log_write('        config_iso_initial_temp_h1 = $r', realArgs=(/config_iso_initial_temp_h1/))
            call mpas_log_write('        config_iso_initial_temp_mt = $r', realArgs=(/config_iso_initial_temp_mt/))
            call mpas_log_write('        config_iso_initial_temp_latS = $r', realArgs=(/config_iso_initial_temp_latS/))
            call mpas_log_write('        config_iso_initial_temp_latN = $r', realArgs=(/config_iso_initial_temp_latN/))
            call mpas_log_write('        config_iso_temperature_sponge_t1 = $r', realArgs=(/config_iso_temperature_sponge_t1/))
            call mpas_log_write('        config_iso_temperature_sponge_h1 = $r', realArgs=(/config_iso_temperature_sponge_h1/))
            call mpas_log_write('        config_iso_temperature_sponge_l1 = $r', realArgs=(/config_iso_temperature_sponge_l1/))
            call mpas_log_write('        config_iso_temperature_sponge_tau1 = $r', realArgs=(/config_iso_temperature_sponge_tau1/))
            call mpas_log_write('        config_iso_temperature_restore_region1_flag = $l', &
logicArgs=(/config_iso_temperature_restore_region1_flag/))
            call mpas_log_write('        config_iso_temperature_restore_t1 = $r', realArgs=(/config_iso_temperature_restore_t1/))
            call mpas_log_write('        config_iso_temperature_restore_lcx1 = $r', &
realArgs=(/config_iso_temperature_restore_lcx1/))
            call mpas_log_write('        config_iso_temperature_restore_lcy1 = $r', &
realArgs=(/config_iso_temperature_restore_lcy1/))
            call mpas_log_write('        config_iso_temperature_restore_region2_flag = $l', &
logicArgs=(/config_iso_temperature_restore_region2_flag/))
            call mpas_log_write('        config_iso_temperature_restore_t2 = $r', realArgs=(/config_iso_temperature_restore_t2/))
            call mpas_log_write('        config_iso_temperature_restore_lcx2 = $r', &
realArgs=(/config_iso_temperature_restore_lcx2/))
            call mpas_log_write('        config_iso_temperature_restore_lcy2 = $r', &
realArgs=(/config_iso_temperature_restore_lcy2/))
            call mpas_log_write('        config_iso_temperature_restore_region3_flag = $l', &
logicArgs=(/config_iso_temperature_restore_region3_flag/))
            call mpas_log_write('        config_iso_temperature_restore_t3 = $r', realArgs=(/config_iso_temperature_restore_t3/))
            call mpas_log_write('        config_iso_temperature_restore_lcx3 = $r', &
realArgs=(/config_iso_temperature_restore_lcx3/))
            call mpas_log_write('        config_iso_temperature_restore_lcy3 = $r', &
realArgs=(/config_iso_temperature_restore_lcy3/))
            call mpas_log_write('        config_iso_temperature_restore_region4_flag = $l', &
logicArgs=(/config_iso_temperature_restore_region4_flag/))
            call mpas_log_write('        config_iso_temperature_restore_t4 = $r', realArgs=(/config_iso_temperature_restore_t4/))
            call mpas_log_write('        config_iso_temperature_restore_lcx4 = $r', &
realArgs=(/config_iso_temperature_restore_lcx4/))
            call mpas_log_write('        config_iso_temperature_restore_lcy4 = $r', &
realArgs=(/config_iso_temperature_restore_lcy4/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record iso.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_iso_vert_levels', config_iso_vert_levels)
      call mpas_pool_add_config(configPool, 'config_iso_main_channel_depth', config_iso_main_channel_depth)
      call mpas_pool_add_config(configPool, 'config_iso_north_wall_lat', config_iso_north_wall_lat)
      call mpas_pool_add_config(configPool, 'config_iso_south_wall_lat', config_iso_south_wall_lat)
      call mpas_pool_add_config(configPool, 'config_iso_ridge_flag', config_iso_ridge_flag)
      call mpas_pool_add_config(configPool, 'config_iso_ridge_center_lon', config_iso_ridge_center_lon)
      call mpas_pool_add_config(configPool, 'config_iso_ridge_height', config_iso_ridge_height)
      call mpas_pool_add_config(configPool, 'config_iso_ridge_width', config_iso_ridge_width)
      call mpas_pool_add_config(configPool, 'config_iso_plateau_flag', config_iso_plateau_flag)
      call mpas_pool_add_config(configPool, 'config_iso_plateau_center_lon', config_iso_plateau_center_lon)
      call mpas_pool_add_config(configPool, 'config_iso_plateau_center_lat', config_iso_plateau_center_lat)
      call mpas_pool_add_config(configPool, 'config_iso_plateau_height', config_iso_plateau_height)
      call mpas_pool_add_config(configPool, 'config_iso_plateau_radius', config_iso_plateau_radius)
      call mpas_pool_add_config(configPool, 'config_iso_plateau_slope_width', config_iso_plateau_slope_width)
      call mpas_pool_add_config(configPool, 'config_iso_shelf_flag', config_iso_shelf_flag)
      call mpas_pool_add_config(configPool, 'config_iso_shelf_depth', config_iso_shelf_depth)
      call mpas_pool_add_config(configPool, 'config_iso_shelf_width', config_iso_shelf_width)
      call mpas_pool_add_config(configPool, 'config_iso_cont_slope_flag', config_iso_cont_slope_flag)
      call mpas_pool_add_config(configPool, 'config_iso_max_cont_slope', config_iso_max_cont_slope)
      call mpas_pool_add_config(configPool, 'config_iso_embayment_flag', config_iso_embayment_flag)
      call mpas_pool_add_config(configPool, 'config_iso_embayment_center_lon', config_iso_embayment_center_lon)
      call mpas_pool_add_config(configPool, 'config_iso_embayment_center_lat', config_iso_embayment_center_lat)
      call mpas_pool_add_config(configPool, 'config_iso_embayment_radius', config_iso_embayment_radius)
      call mpas_pool_add_config(configPool, 'config_iso_embayment_depth', config_iso_embayment_depth)
      call mpas_pool_add_config(configPool, 'config_iso_depression_flag', config_iso_depression_flag)
      call mpas_pool_add_config(configPool, 'config_iso_depression_center_lon', config_iso_depression_center_lon)
      call mpas_pool_add_config(configPool, 'config_iso_depression_south_lat', config_iso_depression_south_lat)
      call mpas_pool_add_config(configPool, 'config_iso_depression_north_lat', config_iso_depression_north_lat)
      call mpas_pool_add_config(configPool, 'config_iso_depression_width', config_iso_depression_width)
      call mpas_pool_add_config(configPool, 'config_iso_depression_depth', config_iso_depression_depth)
      call mpas_pool_add_config(configPool, 'config_iso_salinity', config_iso_salinity)
      call mpas_pool_add_config(configPool, 'config_iso_wind_stress_max', config_iso_wind_stress_max)
      call mpas_pool_add_config(configPool, 'config_iso_acc_wind', config_iso_acc_wind)
      call mpas_pool_add_config(configPool, 'config_iso_asf_wind', config_iso_asf_wind)
      call mpas_pool_add_config(configPool, 'config_iso_wind_trans', config_iso_wind_trans)
      call mpas_pool_add_config(configPool, 'config_iso_heat_flux_south', config_iso_heat_flux_south)
      call mpas_pool_add_config(configPool, 'config_iso_heat_flux_middle', config_iso_heat_flux_middle)
      call mpas_pool_add_config(configPool, 'config_iso_heat_flux_north', config_iso_heat_flux_north)
      call mpas_pool_add_config(configPool, 'config_iso_heat_flux_lat_ss', config_iso_heat_flux_lat_ss)
      call mpas_pool_add_config(configPool, 'config_iso_heat_flux_lat_sm', config_iso_heat_flux_lat_sm)
      call mpas_pool_add_config(configPool, 'config_iso_heat_flux_lat_mn', config_iso_heat_flux_lat_mn)
      call mpas_pool_add_config(configPool, 'config_iso_region1_center_lon', config_iso_region1_center_lon)
      call mpas_pool_add_config(configPool, 'config_iso_region1_center_lat', config_iso_region1_center_lat)
      call mpas_pool_add_config(configPool, 'config_iso_region2_center_lon', config_iso_region2_center_lon)
      call mpas_pool_add_config(configPool, 'config_iso_region2_center_lat', config_iso_region2_center_lat)
      call mpas_pool_add_config(configPool, 'config_iso_region3_center_lon', config_iso_region3_center_lon)
      call mpas_pool_add_config(configPool, 'config_iso_region3_center_lat', config_iso_region3_center_lat)
      call mpas_pool_add_config(configPool, 'config_iso_region4_center_lon', config_iso_region4_center_lon)
      call mpas_pool_add_config(configPool, 'config_iso_region4_center_lat', config_iso_region4_center_lat)
      call mpas_pool_add_config(configPool, 'config_iso_heat_flux_region1_flag', config_iso_heat_flux_region1_flag)
      call mpas_pool_add_config(configPool, 'config_iso_heat_flux_region1', config_iso_heat_flux_region1)
      call mpas_pool_add_config(configPool, 'config_iso_heat_flux_region1_radius', config_iso_heat_flux_region1_radius)
      call mpas_pool_add_config(configPool, 'config_iso_heat_flux_region2_flag', config_iso_heat_flux_region2_flag)
      call mpas_pool_add_config(configPool, 'config_iso_heat_flux_region2', config_iso_heat_flux_region2)
      call mpas_pool_add_config(configPool, 'config_iso_heat_flux_region2_radius', config_iso_heat_flux_region2_radius)
      call mpas_pool_add_config(configPool, 'config_iso_surface_temperature_piston_velocity', &
config_iso_surface_temperature_piston_velocity)
      call mpas_pool_add_config(configPool, 'config_iso_initial_temp_t1', config_iso_initial_temp_t1)
      call mpas_pool_add_config(configPool, 'config_iso_initial_temp_t2', config_iso_initial_temp_t2)
      call mpas_pool_add_config(configPool, 'config_iso_initial_temp_h0', config_iso_initial_temp_h0)
      call mpas_pool_add_config(configPool, 'config_iso_initial_temp_h1', config_iso_initial_temp_h1)
      call mpas_pool_add_config(configPool, 'config_iso_initial_temp_mt', config_iso_initial_temp_mt)
      call mpas_pool_add_config(configPool, 'config_iso_initial_temp_latS', config_iso_initial_temp_latS)
      call mpas_pool_add_config(configPool, 'config_iso_initial_temp_latN', config_iso_initial_temp_latN)
      call mpas_pool_add_config(configPool, 'config_iso_temperature_sponge_t1', config_iso_temperature_sponge_t1)
      call mpas_pool_add_config(configPool, 'config_iso_temperature_sponge_h1', config_iso_temperature_sponge_h1)
      call mpas_pool_add_config(configPool, 'config_iso_temperature_sponge_l1', config_iso_temperature_sponge_l1)
      call mpas_pool_add_config(configPool, 'config_iso_temperature_sponge_tau1', config_iso_temperature_sponge_tau1)
      call mpas_pool_add_config(configPool, 'config_iso_temperature_restore_region1_flag', &
config_iso_temperature_restore_region1_flag)
      call mpas_pool_add_config(configPool, 'config_iso_temperature_restore_t1', config_iso_temperature_restore_t1)
      call mpas_pool_add_config(configPool, 'config_iso_temperature_restore_lcx1', config_iso_temperature_restore_lcx1)
      call mpas_pool_add_config(configPool, 'config_iso_temperature_restore_lcy1', config_iso_temperature_restore_lcy1)
      call mpas_pool_add_config(configPool, 'config_iso_temperature_restore_region2_flag', &
config_iso_temperature_restore_region2_flag)
      call mpas_pool_add_config(configPool, 'config_iso_temperature_restore_t2', config_iso_temperature_restore_t2)
      call mpas_pool_add_config(configPool, 'config_iso_temperature_restore_lcx2', config_iso_temperature_restore_lcx2)
      call mpas_pool_add_config(configPool, 'config_iso_temperature_restore_lcy2', config_iso_temperature_restore_lcy2)
      call mpas_pool_add_config(configPool, 'config_iso_temperature_restore_region3_flag', &
config_iso_temperature_restore_region3_flag)
      call mpas_pool_add_config(configPool, 'config_iso_temperature_restore_t3', config_iso_temperature_restore_t3)
      call mpas_pool_add_config(configPool, 'config_iso_temperature_restore_lcx3', config_iso_temperature_restore_lcx3)
      call mpas_pool_add_config(configPool, 'config_iso_temperature_restore_lcy3', config_iso_temperature_restore_lcy3)
      call mpas_pool_add_config(configPool, 'config_iso_temperature_restore_region4_flag', &
config_iso_temperature_restore_region4_flag)
      call mpas_pool_add_config(configPool, 'config_iso_temperature_restore_t4', config_iso_temperature_restore_t4)
      call mpas_pool_add_config(configPool, 'config_iso_temperature_restore_lcx4', config_iso_temperature_restore_lcx4)
      call mpas_pool_add_config(configPool, 'config_iso_temperature_restore_lcy4', config_iso_temperature_restore_lcy4)

   end subroutine ocn_setup_nmlrec_iso


   subroutine ocn_setup_nmlrec_soma(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      integer :: config_soma_vert_levels = 100
      real (kind=RKIND) :: config_soma_domain_width = 1250000.000000
      real (kind=RKIND) :: config_soma_center_latitude = 35.000000
      real (kind=RKIND) :: config_soma_center_longitude = 0.000000
      real (kind=RKIND) :: config_soma_phi = 0.100000
      real (kind=RKIND) :: config_soma_bottom_depth = 2500.000000
      real (kind=RKIND) :: config_soma_shelf_width = -0.400000
      real (kind=RKIND) :: config_soma_shelf_depth = 100.000000
      real (kind=RKIND) :: config_soma_ref_density = 1000.000000
      real (kind=RKIND) :: config_soma_density_difference = 4.000000
      real (kind=RKIND) :: config_soma_thermocline_depth = 300.000000
      real (kind=RKIND) :: config_soma_density_difference_linear = 0.050000
      real (kind=RKIND) :: config_soma_surface_temperature = 20.000000
      real (kind=RKIND) :: config_soma_surface_salinity = 33.000000
      logical :: config_soma_use_surface_temp_restoring = .false.
      real (kind=RKIND) :: config_soma_surface_temp_restoring_at_center_latitude = 7.500000
      real (kind=RKIND) :: config_soma_surface_temp_restoring_latitude_gradient = 0.500000
      real (kind=RKIND) :: config_soma_restoring_temp_piston_vel = 0.000010

      namelist /soma/ &
         config_soma_vert_levels, &
         config_soma_domain_width, &
         config_soma_center_latitude, &
         config_soma_center_longitude, &
         config_soma_phi, &
         config_soma_bottom_depth, &
         config_soma_shelf_width, &
         config_soma_shelf_depth, &
         config_soma_ref_density, &
         config_soma_density_difference, &
         config_soma_thermocline_depth, &
         config_soma_density_difference_linear, &
         config_soma_surface_temperature, &
         config_soma_surface_salinity, &
         config_soma_use_surface_temp_restoring, &
         config_soma_surface_temp_restoring_at_center_latitude, &
         config_soma_surface_temp_restoring_latitude_gradient, &
         config_soma_restoring_temp_piston_vel
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, soma, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_int(dminfo, config_soma_vert_levels)
         call mpas_dmpar_bcast_real(dminfo, config_soma_domain_width)
         call mpas_dmpar_bcast_real(dminfo, config_soma_center_latitude)
         call mpas_dmpar_bcast_real(dminfo, config_soma_center_longitude)
         call mpas_dmpar_bcast_real(dminfo, config_soma_phi)
         call mpas_dmpar_bcast_real(dminfo, config_soma_bottom_depth)
         call mpas_dmpar_bcast_real(dminfo, config_soma_shelf_width)
         call mpas_dmpar_bcast_real(dminfo, config_soma_shelf_depth)
         call mpas_dmpar_bcast_real(dminfo, config_soma_ref_density)
         call mpas_dmpar_bcast_real(dminfo, config_soma_density_difference)
         call mpas_dmpar_bcast_real(dminfo, config_soma_thermocline_depth)
         call mpas_dmpar_bcast_real(dminfo, config_soma_density_difference_linear)
         call mpas_dmpar_bcast_real(dminfo, config_soma_surface_temperature)
         call mpas_dmpar_bcast_real(dminfo, config_soma_surface_salinity)
         call mpas_dmpar_bcast_logical(dminfo, config_soma_use_surface_temp_restoring)
         call mpas_dmpar_bcast_real(dminfo, config_soma_surface_temp_restoring_at_center_latitude)
         call mpas_dmpar_bcast_real(dminfo, config_soma_surface_temp_restoring_latitude_gradient)
         call mpas_dmpar_bcast_real(dminfo, config_soma_restoring_temp_piston_vel)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record soma')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_soma_vert_levels = $i', intArgs=(/config_soma_vert_levels/))
            call mpas_log_write('        config_soma_domain_width = $r', realArgs=(/config_soma_domain_width/))
            call mpas_log_write('        config_soma_center_latitude = $r', realArgs=(/config_soma_center_latitude/))
            call mpas_log_write('        config_soma_center_longitude = $r', realArgs=(/config_soma_center_longitude/))
            call mpas_log_write('        config_soma_phi = $r', realArgs=(/config_soma_phi/))
            call mpas_log_write('        config_soma_bottom_depth = $r', realArgs=(/config_soma_bottom_depth/))
            call mpas_log_write('        config_soma_shelf_width = $r', realArgs=(/config_soma_shelf_width/))
            call mpas_log_write('        config_soma_shelf_depth = $r', realArgs=(/config_soma_shelf_depth/))
            call mpas_log_write('        config_soma_ref_density = $r', realArgs=(/config_soma_ref_density/))
            call mpas_log_write('        config_soma_density_difference = $r', realArgs=(/config_soma_density_difference/))
            call mpas_log_write('        config_soma_thermocline_depth = $r', realArgs=(/config_soma_thermocline_depth/))
            call mpas_log_write('        config_soma_density_difference_linear = $r', &
realArgs=(/config_soma_density_difference_linear/))
            call mpas_log_write('        config_soma_surface_temperature = $r', realArgs=(/config_soma_surface_temperature/))
            call mpas_log_write('        config_soma_surface_salinity = $r', realArgs=(/config_soma_surface_salinity/))
            call mpas_log_write('        config_soma_use_surface_temp_restoring = $l', &
logicArgs=(/config_soma_use_surface_temp_restoring/))
            call mpas_log_write('        config_soma_surface_temp_restoring_at_center_latitude = $r', &
realArgs=(/config_soma_surface_temp_restoring_at_center_latitude/))
            call mpas_log_write('        config_soma_surface_temp_restoring_latitude_gradient = $r', &
realArgs=(/config_soma_surface_temp_restoring_latitude_gradient/))
            call mpas_log_write('        config_soma_restoring_temp_piston_vel = $r', &
realArgs=(/config_soma_restoring_temp_piston_vel/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record soma.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_soma_vert_levels', config_soma_vert_levels)
      call mpas_pool_add_config(configPool, 'config_soma_domain_width', config_soma_domain_width)
      call mpas_pool_add_config(configPool, 'config_soma_center_latitude', config_soma_center_latitude)
      call mpas_pool_add_config(configPool, 'config_soma_center_longitude', config_soma_center_longitude)
      call mpas_pool_add_config(configPool, 'config_soma_phi', config_soma_phi)
      call mpas_pool_add_config(configPool, 'config_soma_bottom_depth', config_soma_bottom_depth)
      call mpas_pool_add_config(configPool, 'config_soma_shelf_width', config_soma_shelf_width)
      call mpas_pool_add_config(configPool, 'config_soma_shelf_depth', config_soma_shelf_depth)
      call mpas_pool_add_config(configPool, 'config_soma_ref_density', config_soma_ref_density)
      call mpas_pool_add_config(configPool, 'config_soma_density_difference', config_soma_density_difference)
      call mpas_pool_add_config(configPool, 'config_soma_thermocline_depth', config_soma_thermocline_depth)
      call mpas_pool_add_config(configPool, 'config_soma_density_difference_linear', config_soma_density_difference_linear)
      call mpas_pool_add_config(configPool, 'config_soma_surface_temperature', config_soma_surface_temperature)
      call mpas_pool_add_config(configPool, 'config_soma_surface_salinity', config_soma_surface_salinity)
      call mpas_pool_add_config(configPool, 'config_soma_use_surface_temp_restoring', config_soma_use_surface_temp_restoring)
      call mpas_pool_add_config(configPool, 'config_soma_surface_temp_restoring_at_center_latitude', &
config_soma_surface_temp_restoring_at_center_latitude)
      call mpas_pool_add_config(configPool, 'config_soma_surface_temp_restoring_latitude_gradient', &
config_soma_surface_temp_restoring_latitude_gradient)
      call mpas_pool_add_config(configPool, 'config_soma_restoring_temp_piston_vel', config_soma_restoring_temp_piston_vel)

   end subroutine ocn_setup_nmlrec_soma


   subroutine ocn_setup_nmlrec_ziso(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      integer :: config_ziso_vert_levels = 100
      logical :: config_ziso_add_easterly_wind_stress_ASF = .false.
      real (kind=RKIND) :: config_ziso_wind_transition_position = 800000.000000
      real (kind=RKIND) :: config_ziso_antarctic_shelf_front_width = 600000.000000
      real (kind=RKIND) :: config_ziso_wind_stress_shelf_front_max = -0.050000
      logical :: config_ziso_use_slopping_bathymetry = .false.
      real (kind=RKIND) :: config_ziso_meridional_extent = 2000000.000000
      real (kind=RKIND) :: config_ziso_zonal_extent = 1000000.000000
      real (kind=RKIND) :: config_ziso_bottom_depth = 2500.000000
      real (kind=RKIND) :: config_ziso_shelf_depth = 500.000000
      real (kind=RKIND) :: config_ziso_slope_half_width = 100000.000000
      real (kind=RKIND) :: config_ziso_slope_center_position = 500000.000000
      real (kind=RKIND) :: config_ziso_reference_coriolis = -0.000100
      real (kind=RKIND) :: config_ziso_coriolis_gradient = 0.000000
      real (kind=RKIND) :: config_ziso_wind_stress_max = 0.200000
      real (kind=RKIND) :: config_ziso_mean_restoring_temp = 3.000000
      real (kind=RKIND) :: config_ziso_restoring_temp_dev_ta = 2.000000
      real (kind=RKIND) :: config_ziso_restoring_temp_dev_tb = 2.000000
      real (kind=RKIND) :: config_ziso_restoring_temp_tau = 30.000000
      real (kind=RKIND) :: config_ziso_restoring_temp_piston_vel = 0.000019
      real (kind=RKIND) :: config_ziso_restoring_temp_ze = 1250.000000
      real (kind=RKIND) :: config_ziso_restoring_sponge_l = 80000.000000
      real (kind=RKIND) :: config_ziso_initial_temp_t1 = 6.000000
      real (kind=RKIND) :: config_ziso_initial_temp_t2 = 3.600000
      real (kind=RKIND) :: config_ziso_initial_temp_h1 = 300.000000
      real (kind=RKIND) :: config_ziso_initial_temp_mt = 0.000075
      logical :: config_ziso_frazil_enable = .false.
      real (kind=RKIND) :: config_ziso_frazil_temperature_anomaly = -3.000000

      namelist /ziso/ &
         config_ziso_vert_levels, &
         config_ziso_add_easterly_wind_stress_ASF, &
         config_ziso_wind_transition_position, &
         config_ziso_antarctic_shelf_front_width, &
         config_ziso_wind_stress_shelf_front_max, &
         config_ziso_use_slopping_bathymetry, &
         config_ziso_meridional_extent, &
         config_ziso_zonal_extent, &
         config_ziso_bottom_depth, &
         config_ziso_shelf_depth, &
         config_ziso_slope_half_width, &
         config_ziso_slope_center_position, &
         config_ziso_reference_coriolis, &
         config_ziso_coriolis_gradient, &
         config_ziso_wind_stress_max, &
         config_ziso_mean_restoring_temp, &
         config_ziso_restoring_temp_dev_ta, &
         config_ziso_restoring_temp_dev_tb, &
         config_ziso_restoring_temp_tau, &
         config_ziso_restoring_temp_piston_vel, &
         config_ziso_restoring_temp_ze, &
         config_ziso_restoring_sponge_l, &
         config_ziso_initial_temp_t1, &
         config_ziso_initial_temp_t2, &
         config_ziso_initial_temp_h1, &
         config_ziso_initial_temp_mt, &
         config_ziso_frazil_enable, &
         config_ziso_frazil_temperature_anomaly
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, ziso, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_int(dminfo, config_ziso_vert_levels)
         call mpas_dmpar_bcast_logical(dminfo, config_ziso_add_easterly_wind_stress_ASF)
         call mpas_dmpar_bcast_real(dminfo, config_ziso_wind_transition_position)
         call mpas_dmpar_bcast_real(dminfo, config_ziso_antarctic_shelf_front_width)
         call mpas_dmpar_bcast_real(dminfo, config_ziso_wind_stress_shelf_front_max)
         call mpas_dmpar_bcast_logical(dminfo, config_ziso_use_slopping_bathymetry)
         call mpas_dmpar_bcast_real(dminfo, config_ziso_meridional_extent)
         call mpas_dmpar_bcast_real(dminfo, config_ziso_zonal_extent)
         call mpas_dmpar_bcast_real(dminfo, config_ziso_bottom_depth)
         call mpas_dmpar_bcast_real(dminfo, config_ziso_shelf_depth)
         call mpas_dmpar_bcast_real(dminfo, config_ziso_slope_half_width)
         call mpas_dmpar_bcast_real(dminfo, config_ziso_slope_center_position)
         call mpas_dmpar_bcast_real(dminfo, config_ziso_reference_coriolis)
         call mpas_dmpar_bcast_real(dminfo, config_ziso_coriolis_gradient)
         call mpas_dmpar_bcast_real(dminfo, config_ziso_wind_stress_max)
         call mpas_dmpar_bcast_real(dminfo, config_ziso_mean_restoring_temp)
         call mpas_dmpar_bcast_real(dminfo, config_ziso_restoring_temp_dev_ta)
         call mpas_dmpar_bcast_real(dminfo, config_ziso_restoring_temp_dev_tb)
         call mpas_dmpar_bcast_real(dminfo, config_ziso_restoring_temp_tau)
         call mpas_dmpar_bcast_real(dminfo, config_ziso_restoring_temp_piston_vel)
         call mpas_dmpar_bcast_real(dminfo, config_ziso_restoring_temp_ze)
         call mpas_dmpar_bcast_real(dminfo, config_ziso_restoring_sponge_l)
         call mpas_dmpar_bcast_real(dminfo, config_ziso_initial_temp_t1)
         call mpas_dmpar_bcast_real(dminfo, config_ziso_initial_temp_t2)
         call mpas_dmpar_bcast_real(dminfo, config_ziso_initial_temp_h1)
         call mpas_dmpar_bcast_real(dminfo, config_ziso_initial_temp_mt)
         call mpas_dmpar_bcast_logical(dminfo, config_ziso_frazil_enable)
         call mpas_dmpar_bcast_real(dminfo, config_ziso_frazil_temperature_anomaly)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record ziso')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_ziso_vert_levels = $i', intArgs=(/config_ziso_vert_levels/))
            call mpas_log_write('        config_ziso_add_easterly_wind_stress_ASF = $l', &
logicArgs=(/config_ziso_add_easterly_wind_stress_ASF/))
            call mpas_log_write('        config_ziso_wind_transition_position = $r', &
realArgs=(/config_ziso_wind_transition_position/))
            call mpas_log_write('        config_ziso_antarctic_shelf_front_width = $r', &
realArgs=(/config_ziso_antarctic_shelf_front_width/))
            call mpas_log_write('        config_ziso_wind_stress_shelf_front_max = $r', &
realArgs=(/config_ziso_wind_stress_shelf_front_max/))
            call mpas_log_write('        config_ziso_use_slopping_bathymetry = $l', &
logicArgs=(/config_ziso_use_slopping_bathymetry/))
            call mpas_log_write('        config_ziso_meridional_extent = $r', realArgs=(/config_ziso_meridional_extent/))
            call mpas_log_write('        config_ziso_zonal_extent = $r', realArgs=(/config_ziso_zonal_extent/))
            call mpas_log_write('        config_ziso_bottom_depth = $r', realArgs=(/config_ziso_bottom_depth/))
            call mpas_log_write('        config_ziso_shelf_depth = $r', realArgs=(/config_ziso_shelf_depth/))
            call mpas_log_write('        config_ziso_slope_half_width = $r', realArgs=(/config_ziso_slope_half_width/))
            call mpas_log_write('        config_ziso_slope_center_position = $r', realArgs=(/config_ziso_slope_center_position/))
            call mpas_log_write('        config_ziso_reference_coriolis = $r', realArgs=(/config_ziso_reference_coriolis/))
            call mpas_log_write('        config_ziso_coriolis_gradient = $r', realArgs=(/config_ziso_coriolis_gradient/))
            call mpas_log_write('        config_ziso_wind_stress_max = $r', realArgs=(/config_ziso_wind_stress_max/))
            call mpas_log_write('        config_ziso_mean_restoring_temp = $r', realArgs=(/config_ziso_mean_restoring_temp/))
            call mpas_log_write('        config_ziso_restoring_temp_dev_ta = $r', realArgs=(/config_ziso_restoring_temp_dev_ta/))
            call mpas_log_write('        config_ziso_restoring_temp_dev_tb = $r', realArgs=(/config_ziso_restoring_temp_dev_tb/))
            call mpas_log_write('        config_ziso_restoring_temp_tau = $r', realArgs=(/config_ziso_restoring_temp_tau/))
            call mpas_log_write('        config_ziso_restoring_temp_piston_vel = $r', &
realArgs=(/config_ziso_restoring_temp_piston_vel/))
            call mpas_log_write('        config_ziso_restoring_temp_ze = $r', realArgs=(/config_ziso_restoring_temp_ze/))
            call mpas_log_write('        config_ziso_restoring_sponge_l = $r', realArgs=(/config_ziso_restoring_sponge_l/))
            call mpas_log_write('        config_ziso_initial_temp_t1 = $r', realArgs=(/config_ziso_initial_temp_t1/))
            call mpas_log_write('        config_ziso_initial_temp_t2 = $r', realArgs=(/config_ziso_initial_temp_t2/))
            call mpas_log_write('        config_ziso_initial_temp_h1 = $r', realArgs=(/config_ziso_initial_temp_h1/))
            call mpas_log_write('        config_ziso_initial_temp_mt = $r', realArgs=(/config_ziso_initial_temp_mt/))
            call mpas_log_write('        config_ziso_frazil_enable = $l', logicArgs=(/config_ziso_frazil_enable/))
            call mpas_log_write('        config_ziso_frazil_temperature_anomaly = $r', &
realArgs=(/config_ziso_frazil_temperature_anomaly/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record ziso.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_ziso_vert_levels', config_ziso_vert_levels)
      call mpas_pool_add_config(configPool, 'config_ziso_add_easterly_wind_stress_ASF', config_ziso_add_easterly_wind_stress_ASF)
      call mpas_pool_add_config(configPool, 'config_ziso_wind_transition_position', config_ziso_wind_transition_position)
      call mpas_pool_add_config(configPool, 'config_ziso_antarctic_shelf_front_width', config_ziso_antarctic_shelf_front_width)
      call mpas_pool_add_config(configPool, 'config_ziso_wind_stress_shelf_front_max', config_ziso_wind_stress_shelf_front_max)
      call mpas_pool_add_config(configPool, 'config_ziso_use_slopping_bathymetry', config_ziso_use_slopping_bathymetry)
      call mpas_pool_add_config(configPool, 'config_ziso_meridional_extent', config_ziso_meridional_extent)
      call mpas_pool_add_config(configPool, 'config_ziso_zonal_extent', config_ziso_zonal_extent)
      call mpas_pool_add_config(configPool, 'config_ziso_bottom_depth', config_ziso_bottom_depth)
      call mpas_pool_add_config(configPool, 'config_ziso_shelf_depth', config_ziso_shelf_depth)
      call mpas_pool_add_config(configPool, 'config_ziso_slope_half_width', config_ziso_slope_half_width)
      call mpas_pool_add_config(configPool, 'config_ziso_slope_center_position', config_ziso_slope_center_position)
      call mpas_pool_add_config(configPool, 'config_ziso_reference_coriolis', config_ziso_reference_coriolis)
      call mpas_pool_add_config(configPool, 'config_ziso_coriolis_gradient', config_ziso_coriolis_gradient)
      call mpas_pool_add_config(configPool, 'config_ziso_wind_stress_max', config_ziso_wind_stress_max)
      call mpas_pool_add_config(configPool, 'config_ziso_mean_restoring_temp', config_ziso_mean_restoring_temp)
      call mpas_pool_add_config(configPool, 'config_ziso_restoring_temp_dev_ta', config_ziso_restoring_temp_dev_ta)
      call mpas_pool_add_config(configPool, 'config_ziso_restoring_temp_dev_tb', config_ziso_restoring_temp_dev_tb)
      call mpas_pool_add_config(configPool, 'config_ziso_restoring_temp_tau', config_ziso_restoring_temp_tau)
      call mpas_pool_add_config(configPool, 'config_ziso_restoring_temp_piston_vel', config_ziso_restoring_temp_piston_vel)
      call mpas_pool_add_config(configPool, 'config_ziso_restoring_temp_ze', config_ziso_restoring_temp_ze)
      call mpas_pool_add_config(configPool, 'config_ziso_restoring_sponge_l', config_ziso_restoring_sponge_l)
      call mpas_pool_add_config(configPool, 'config_ziso_initial_temp_t1', config_ziso_initial_temp_t1)
      call mpas_pool_add_config(configPool, 'config_ziso_initial_temp_t2', config_ziso_initial_temp_t2)
      call mpas_pool_add_config(configPool, 'config_ziso_initial_temp_h1', config_ziso_initial_temp_h1)
      call mpas_pool_add_config(configPool, 'config_ziso_initial_temp_mt', config_ziso_initial_temp_mt)
      call mpas_pool_add_config(configPool, 'config_ziso_frazil_enable', config_ziso_frazil_enable)
      call mpas_pool_add_config(configPool, 'config_ziso_frazil_temperature_anomaly', config_ziso_frazil_temperature_anomaly)

   end subroutine ocn_setup_nmlrec_ziso


   subroutine ocn_setup_nmlrec_sub_ice_shelf_2D(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      integer :: config_sub_ice_shelf_2D_vert_levels = 20
      real (kind=RKIND) :: config_sub_ice_shelf_2D_bottom_depth = 2000.000000
      real (kind=RKIND) :: config_sub_ice_shelf_2D_cavity_thickness = 25.000000
      real (kind=RKIND) :: config_sub_ice_shelf_2D_slope_height = 500.000000
      real (kind=RKIND) :: config_sub_ice_shelf_2D_edge_width = 15000.000000
      real (kind=RKIND) :: config_sub_ice_shelf_2D_y1 = 30000.000000
      real (kind=RKIND) :: config_sub_ice_shelf_2D_y2 = 60000.000000
      real (kind=RKIND) :: config_sub_ice_shelf_2D_temperature = 1.000000
      real (kind=RKIND) :: config_sub_ice_shelf_2D_surface_salinity = 34.500000
      real (kind=RKIND) :: config_sub_ice_shelf_2D_bottom_salinity = 34.700000

      namelist /sub_ice_shelf_2D/ &
         config_sub_ice_shelf_2D_vert_levels, &
         config_sub_ice_shelf_2D_bottom_depth, &
         config_sub_ice_shelf_2D_cavity_thickness, &
         config_sub_ice_shelf_2D_slope_height, &
         config_sub_ice_shelf_2D_edge_width, &
         config_sub_ice_shelf_2D_y1, &
         config_sub_ice_shelf_2D_y2, &
         config_sub_ice_shelf_2D_temperature, &
         config_sub_ice_shelf_2D_surface_salinity, &
         config_sub_ice_shelf_2D_bottom_salinity
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, sub_ice_shelf_2D, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_int(dminfo, config_sub_ice_shelf_2D_vert_levels)
         call mpas_dmpar_bcast_real(dminfo, config_sub_ice_shelf_2D_bottom_depth)
         call mpas_dmpar_bcast_real(dminfo, config_sub_ice_shelf_2D_cavity_thickness)
         call mpas_dmpar_bcast_real(dminfo, config_sub_ice_shelf_2D_slope_height)
         call mpas_dmpar_bcast_real(dminfo, config_sub_ice_shelf_2D_edge_width)
         call mpas_dmpar_bcast_real(dminfo, config_sub_ice_shelf_2D_y1)
         call mpas_dmpar_bcast_real(dminfo, config_sub_ice_shelf_2D_y2)
         call mpas_dmpar_bcast_real(dminfo, config_sub_ice_shelf_2D_temperature)
         call mpas_dmpar_bcast_real(dminfo, config_sub_ice_shelf_2D_surface_salinity)
         call mpas_dmpar_bcast_real(dminfo, config_sub_ice_shelf_2D_bottom_salinity)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record sub_ice_shelf_2D')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_sub_ice_shelf_2D_vert_levels = $i', &
intArgs=(/config_sub_ice_shelf_2D_vert_levels/))
            call mpas_log_write('        config_sub_ice_shelf_2D_bottom_depth = $r', &
realArgs=(/config_sub_ice_shelf_2D_bottom_depth/))
            call mpas_log_write('        config_sub_ice_shelf_2D_cavity_thickness = $r', &
realArgs=(/config_sub_ice_shelf_2D_cavity_thickness/))
            call mpas_log_write('        config_sub_ice_shelf_2D_slope_height = $r', &
realArgs=(/config_sub_ice_shelf_2D_slope_height/))
            call mpas_log_write('        config_sub_ice_shelf_2D_edge_width = $r', realArgs=(/config_sub_ice_shelf_2D_edge_width/))
            call mpas_log_write('        config_sub_ice_shelf_2D_y1 = $r', realArgs=(/config_sub_ice_shelf_2D_y1/))
            call mpas_log_write('        config_sub_ice_shelf_2D_y2 = $r', realArgs=(/config_sub_ice_shelf_2D_y2/))
            call mpas_log_write('        config_sub_ice_shelf_2D_temperature = $r', &
realArgs=(/config_sub_ice_shelf_2D_temperature/))
            call mpas_log_write('        config_sub_ice_shelf_2D_surface_salinity = $r', &
realArgs=(/config_sub_ice_shelf_2D_surface_salinity/))
            call mpas_log_write('        config_sub_ice_shelf_2D_bottom_salinity = $r', &
realArgs=(/config_sub_ice_shelf_2D_bottom_salinity/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record sub_ice_shelf_2D.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_sub_ice_shelf_2D_vert_levels', config_sub_ice_shelf_2D_vert_levels)
      call mpas_pool_add_config(configPool, 'config_sub_ice_shelf_2D_bottom_depth', config_sub_ice_shelf_2D_bottom_depth)
      call mpas_pool_add_config(configPool, 'config_sub_ice_shelf_2D_cavity_thickness', config_sub_ice_shelf_2D_cavity_thickness)
      call mpas_pool_add_config(configPool, 'config_sub_ice_shelf_2D_slope_height', config_sub_ice_shelf_2D_slope_height)
      call mpas_pool_add_config(configPool, 'config_sub_ice_shelf_2D_edge_width', config_sub_ice_shelf_2D_edge_width)
      call mpas_pool_add_config(configPool, 'config_sub_ice_shelf_2D_y1', config_sub_ice_shelf_2D_y1)
      call mpas_pool_add_config(configPool, 'config_sub_ice_shelf_2D_y2', config_sub_ice_shelf_2D_y2)
      call mpas_pool_add_config(configPool, 'config_sub_ice_shelf_2D_temperature', config_sub_ice_shelf_2D_temperature)
      call mpas_pool_add_config(configPool, 'config_sub_ice_shelf_2D_surface_salinity', config_sub_ice_shelf_2D_surface_salinity)
      call mpas_pool_add_config(configPool, 'config_sub_ice_shelf_2D_bottom_salinity', config_sub_ice_shelf_2D_bottom_salinity)

   end subroutine ocn_setup_nmlrec_sub_ice_shelf_2D


   subroutine ocn_setup_nmlrec_periodic_planar(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      integer :: config_periodic_planar_vert_levels = 100
      real (kind=RKIND) :: config_periodic_planar_bottom_depth = 2500.000000
      real (kind=RKIND) :: config_periodic_planar_velocity_strength = 1.000000

      namelist /periodic_planar/ &
         config_periodic_planar_vert_levels, &
         config_periodic_planar_bottom_depth, &
         config_periodic_planar_velocity_strength
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, periodic_planar, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_int(dminfo, config_periodic_planar_vert_levels)
         call mpas_dmpar_bcast_real(dminfo, config_periodic_planar_bottom_depth)
         call mpas_dmpar_bcast_real(dminfo, config_periodic_planar_velocity_strength)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record periodic_planar')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_periodic_planar_vert_levels = $i', intArgs=(/config_periodic_planar_vert_levels/))
            call mpas_log_write('        config_periodic_planar_bottom_depth = $r', &
realArgs=(/config_periodic_planar_bottom_depth/))
            call mpas_log_write('        config_periodic_planar_velocity_strength = $r', &
realArgs=(/config_periodic_planar_velocity_strength/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record periodic_planar.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_periodic_planar_vert_levels', config_periodic_planar_vert_levels)
      call mpas_pool_add_config(configPool, 'config_periodic_planar_bottom_depth', config_periodic_planar_bottom_depth)
      call mpas_pool_add_config(configPool, 'config_periodic_planar_velocity_strength', config_periodic_planar_velocity_strength)

   end subroutine ocn_setup_nmlrec_periodic_planar


   subroutine ocn_setup_nmlrec_ecosys_column(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      integer :: config_ecosys_column_vert_levels = 100
      character (len=StrKIND) :: config_ecosys_column_vertical_grid = '100layerE3SMv1'
      character (len=StrKIND) :: config_ecosys_column_TS_filename = 'unknown'
      character (len=StrKIND) :: config_ecosys_column_ecosys_filename = 'unknown'
      real (kind=RKIND) :: config_ecosys_column_bottom_depth = 6000.000000

      namelist /ecosys_column/ &
         config_ecosys_column_vert_levels, &
         config_ecosys_column_vertical_grid, &
         config_ecosys_column_TS_filename, &
         config_ecosys_column_ecosys_filename, &
         config_ecosys_column_bottom_depth
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, ecosys_column, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_int(dminfo, config_ecosys_column_vert_levels)
         call mpas_dmpar_bcast_char(dminfo, config_ecosys_column_vertical_grid)
         call mpas_dmpar_bcast_char(dminfo, config_ecosys_column_TS_filename)
         call mpas_dmpar_bcast_char(dminfo, config_ecosys_column_ecosys_filename)
         call mpas_dmpar_bcast_real(dminfo, config_ecosys_column_bottom_depth)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record ecosys_column')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_ecosys_column_vert_levels = $i', intArgs=(/config_ecosys_column_vert_levels/))
            call mpas_log_write('        config_ecosys_column_vertical_grid = '&
//''//mpas_log_escape_dollars(config_ecosys_column_vertical_grid))
            call mpas_log_write('        config_ecosys_column_TS_filename = '&
//''//mpas_log_escape_dollars(config_ecosys_column_TS_filename))
            call mpas_log_write('        config_ecosys_column_ecosys_filename = '&
//''//mpas_log_escape_dollars(config_ecosys_column_ecosys_filename))
            call mpas_log_write('        config_ecosys_column_bottom_depth = $r', realArgs=(/config_ecosys_column_bottom_depth/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record ecosys_column.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_ecosys_column_vert_levels', config_ecosys_column_vert_levels)
      call mpas_pool_add_config(configPool, 'config_ecosys_column_vertical_grid', config_ecosys_column_vertical_grid)
      call mpas_pool_add_config(configPool, 'config_ecosys_column_TS_filename', config_ecosys_column_TS_filename)
      call mpas_pool_add_config(configPool, 'config_ecosys_column_ecosys_filename', config_ecosys_column_ecosys_filename)
      call mpas_pool_add_config(configPool, 'config_ecosys_column_bottom_depth', config_ecosys_column_bottom_depth)

   end subroutine ocn_setup_nmlrec_ecosys_column


   subroutine ocn_setup_nmlrec_sea_mount(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      integer :: config_sea_mount_vert_levels = 10
      character (len=StrKIND) :: config_sea_mount_layer_type = 'sigma'
      character (len=StrKIND) :: config_sea_mount_stratification_type = 'exponential'
      real (kind=RKIND) :: config_sea_mount_density_coef_linear = 1024.000000
      real (kind=RKIND) :: config_sea_mount_density_coef_exp = 1028.000000
      real (kind=RKIND) :: config_sea_mount_density_gradient_linear = 0.100000
      real (kind=RKIND) :: config_sea_mount_density_gradient_exp = 3.000000
      real (kind=RKIND) :: config_sea_mount_density_depth_linear = 4500.000000
      real (kind=RKIND) :: config_sea_mount_density_depth_exp = 500.000000
      real (kind=RKIND) :: config_sea_mount_density_ref = 1028.000000
      real (kind=RKIND) :: config_sea_mount_density_Tref = 5.000000
      real (kind=RKIND) :: config_sea_mount_density_alpha = 0.200000
      real (kind=RKIND) :: config_sea_mount_bottom_depth = 5000.000000
      real (kind=RKIND) :: config_sea_mount_height = 4500.000000
      real (kind=RKIND) :: config_sea_mount_radius = 10000.000000
      real (kind=RKIND) :: config_sea_mount_width = 40000.000000
      real (kind=RKIND) :: config_sea_mount_salinity = 35.000000
      real (kind=RKIND) :: config_sea_mount_coriolis_parameter = -0.000100

      namelist /sea_mount/ &
         config_sea_mount_vert_levels, &
         config_sea_mount_layer_type, &
         config_sea_mount_stratification_type, &
         config_sea_mount_density_coef_linear, &
         config_sea_mount_density_coef_exp, &
         config_sea_mount_density_gradient_linear, &
         config_sea_mount_density_gradient_exp, &
         config_sea_mount_density_depth_linear, &
         config_sea_mount_density_depth_exp, &
         config_sea_mount_density_ref, &
         config_sea_mount_density_Tref, &
         config_sea_mount_density_alpha, &
         config_sea_mount_bottom_depth, &
         config_sea_mount_height, &
         config_sea_mount_radius, &
         config_sea_mount_width, &
         config_sea_mount_salinity, &
         config_sea_mount_coriolis_parameter
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, sea_mount, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_int(dminfo, config_sea_mount_vert_levels)
         call mpas_dmpar_bcast_char(dminfo, config_sea_mount_layer_type)
         call mpas_dmpar_bcast_char(dminfo, config_sea_mount_stratification_type)
         call mpas_dmpar_bcast_real(dminfo, config_sea_mount_density_coef_linear)
         call mpas_dmpar_bcast_real(dminfo, config_sea_mount_density_coef_exp)
         call mpas_dmpar_bcast_real(dminfo, config_sea_mount_density_gradient_linear)
         call mpas_dmpar_bcast_real(dminfo, config_sea_mount_density_gradient_exp)
         call mpas_dmpar_bcast_real(dminfo, config_sea_mount_density_depth_linear)
         call mpas_dmpar_bcast_real(dminfo, config_sea_mount_density_depth_exp)
         call mpas_dmpar_bcast_real(dminfo, config_sea_mount_density_ref)
         call mpas_dmpar_bcast_real(dminfo, config_sea_mount_density_Tref)
         call mpas_dmpar_bcast_real(dminfo, config_sea_mount_density_alpha)
         call mpas_dmpar_bcast_real(dminfo, config_sea_mount_bottom_depth)
         call mpas_dmpar_bcast_real(dminfo, config_sea_mount_height)
         call mpas_dmpar_bcast_real(dminfo, config_sea_mount_radius)
         call mpas_dmpar_bcast_real(dminfo, config_sea_mount_width)
         call mpas_dmpar_bcast_real(dminfo, config_sea_mount_salinity)
         call mpas_dmpar_bcast_real(dminfo, config_sea_mount_coriolis_parameter)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record sea_mount')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_sea_mount_vert_levels = $i', intArgs=(/config_sea_mount_vert_levels/))
            call mpas_log_write('        config_sea_mount_layer_type = '//mpas_log_escape_dollars(config_sea_mount_layer_type))
            call mpas_log_write('        config_sea_mount_stratification_type = '&
//''//mpas_log_escape_dollars(config_sea_mount_stratification_type))
            call mpas_log_write('        config_sea_mount_density_coef_linear = $r', &
realArgs=(/config_sea_mount_density_coef_linear/))
            call mpas_log_write('        config_sea_mount_density_coef_exp = $r', realArgs=(/config_sea_mount_density_coef_exp/))
            call mpas_log_write('        config_sea_mount_density_gradient_linear = $r', &
realArgs=(/config_sea_mount_density_gradient_linear/))
            call mpas_log_write('        config_sea_mount_density_gradient_exp = $r', &
realArgs=(/config_sea_mount_density_gradient_exp/))
            call mpas_log_write('        config_sea_mount_density_depth_linear = $r', &
realArgs=(/config_sea_mount_density_depth_linear/))
            call mpas_log_write('        config_sea_mount_density_depth_exp = $r', realArgs=(/config_sea_mount_density_depth_exp/))
            call mpas_log_write('        config_sea_mount_density_ref = $r', realArgs=(/config_sea_mount_density_ref/))
            call mpas_log_write('        config_sea_mount_density_Tref = $r', realArgs=(/config_sea_mount_density_Tref/))
            call mpas_log_write('        config_sea_mount_density_alpha = $r', realArgs=(/config_sea_mount_density_alpha/))
            call mpas_log_write('        config_sea_mount_bottom_depth = $r', realArgs=(/config_sea_mount_bottom_depth/))
            call mpas_log_write('        config_sea_mount_height = $r', realArgs=(/config_sea_mount_height/))
            call mpas_log_write('        config_sea_mount_radius = $r', realArgs=(/config_sea_mount_radius/))
            call mpas_log_write('        config_sea_mount_width = $r', realArgs=(/config_sea_mount_width/))
            call mpas_log_write('        config_sea_mount_salinity = $r', realArgs=(/config_sea_mount_salinity/))
            call mpas_log_write('        config_sea_mount_coriolis_parameter = $r', &
realArgs=(/config_sea_mount_coriolis_parameter/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record sea_mount.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_sea_mount_vert_levels', config_sea_mount_vert_levels)
      call mpas_pool_add_config(configPool, 'config_sea_mount_layer_type', config_sea_mount_layer_type)
      call mpas_pool_add_config(configPool, 'config_sea_mount_stratification_type', config_sea_mount_stratification_type)
      call mpas_pool_add_config(configPool, 'config_sea_mount_density_coef_linear', config_sea_mount_density_coef_linear)
      call mpas_pool_add_config(configPool, 'config_sea_mount_density_coef_exp', config_sea_mount_density_coef_exp)
      call mpas_pool_add_config(configPool, 'config_sea_mount_density_gradient_linear', config_sea_mount_density_gradient_linear)
      call mpas_pool_add_config(configPool, 'config_sea_mount_density_gradient_exp', config_sea_mount_density_gradient_exp)
      call mpas_pool_add_config(configPool, 'config_sea_mount_density_depth_linear', config_sea_mount_density_depth_linear)
      call mpas_pool_add_config(configPool, 'config_sea_mount_density_depth_exp', config_sea_mount_density_depth_exp)
      call mpas_pool_add_config(configPool, 'config_sea_mount_density_ref', config_sea_mount_density_ref)
      call mpas_pool_add_config(configPool, 'config_sea_mount_density_Tref', config_sea_mount_density_Tref)
      call mpas_pool_add_config(configPool, 'config_sea_mount_density_alpha', config_sea_mount_density_alpha)
      call mpas_pool_add_config(configPool, 'config_sea_mount_bottom_depth', config_sea_mount_bottom_depth)
      call mpas_pool_add_config(configPool, 'config_sea_mount_height', config_sea_mount_height)
      call mpas_pool_add_config(configPool, 'config_sea_mount_radius', config_sea_mount_radius)
      call mpas_pool_add_config(configPool, 'config_sea_mount_width', config_sea_mount_width)
      call mpas_pool_add_config(configPool, 'config_sea_mount_salinity', config_sea_mount_salinity)
      call mpas_pool_add_config(configPool, 'config_sea_mount_coriolis_parameter', config_sea_mount_coriolis_parameter)

   end subroutine ocn_setup_nmlrec_sea_mount


   subroutine ocn_setup_nmlrec_isomip(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      integer :: config_isomip_vert_levels = 30
      character (len=StrKIND) :: config_isomip_vertical_level_distribution = 'constant'
      real (kind=RKIND) :: config_isomip_bottom_depth = -900.000000
      real (kind=RKIND) :: config_isomip_temperature = -1.900000
      real (kind=RKIND) :: config_isomip_salinity = 34.400000
      real (kind=RKIND) :: config_isomip_restoring_temperature = -1.900000
      real (kind=RKIND) :: config_isomip_temperature_piston_velocity = 0.000012
      real (kind=RKIND) :: config_isomip_restoring_salinity = 34.400000
      real (kind=RKIND) :: config_isomip_salinity_piston_velocity = 0.000012
      real (kind=RKIND) :: config_isomip_coriolis_parameter = -0.000140
      real (kind=RKIND) :: config_isomip_southern_boundary = 0.000000
      real (kind=RKIND) :: config_isomip_northern_boundary = 1000000.000000
      real (kind=RKIND) :: config_isomip_western_boundary = 0.000000
      real (kind=RKIND) :: config_isomip_eastern_boundary = 500000.000000
      real (kind=RKIND) :: config_isomip_y1 = 0.000000
      real (kind=RKIND) :: config_isomip_z1 = -700.000000
      real (kind=RKIND) :: config_isomip_ice_fraction1 = 1.000000
      real (kind=RKIND) :: config_isomip_y2 = 400000.000000
      real (kind=RKIND) :: config_isomip_z2 = -200.000000
      real (kind=RKIND) :: config_isomip_ice_fraction2 = 1.000000
      real (kind=RKIND) :: config_isomip_y3 = 1000000.000000
      real (kind=RKIND) :: config_isomip_z3 = -200.000000
      real (kind=RKIND) :: config_isomip_ice_fraction3 = 1.000000

      namelist /isomip/ &
         config_isomip_vert_levels, &
         config_isomip_vertical_level_distribution, &
         config_isomip_bottom_depth, &
         config_isomip_temperature, &
         config_isomip_salinity, &
         config_isomip_restoring_temperature, &
         config_isomip_temperature_piston_velocity, &
         config_isomip_restoring_salinity, &
         config_isomip_salinity_piston_velocity, &
         config_isomip_coriolis_parameter, &
         config_isomip_southern_boundary, &
         config_isomip_northern_boundary, &
         config_isomip_western_boundary, &
         config_isomip_eastern_boundary, &
         config_isomip_y1, &
         config_isomip_z1, &
         config_isomip_ice_fraction1, &
         config_isomip_y2, &
         config_isomip_z2, &
         config_isomip_ice_fraction2, &
         config_isomip_y3, &
         config_isomip_z3, &
         config_isomip_ice_fraction3
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, isomip, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_int(dminfo, config_isomip_vert_levels)
         call mpas_dmpar_bcast_char(dminfo, config_isomip_vertical_level_distribution)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_bottom_depth)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_temperature)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_salinity)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_restoring_temperature)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_temperature_piston_velocity)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_restoring_salinity)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_salinity_piston_velocity)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_coriolis_parameter)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_southern_boundary)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_northern_boundary)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_western_boundary)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_eastern_boundary)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_y1)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_z1)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_ice_fraction1)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_y2)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_z2)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_ice_fraction2)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_y3)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_z3)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_ice_fraction3)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record isomip')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_isomip_vert_levels = $i', intArgs=(/config_isomip_vert_levels/))
            call mpas_log_write('        config_isomip_vertical_level_distribution = '&
//''//mpas_log_escape_dollars(config_isomip_vertical_level_distribution))
            call mpas_log_write('        config_isomip_bottom_depth = $r', realArgs=(/config_isomip_bottom_depth/))
            call mpas_log_write('        config_isomip_temperature = $r', realArgs=(/config_isomip_temperature/))
            call mpas_log_write('        config_isomip_salinity = $r', realArgs=(/config_isomip_salinity/))
            call mpas_log_write('        config_isomip_restoring_temperature = $r', &
realArgs=(/config_isomip_restoring_temperature/))
            call mpas_log_write('        config_isomip_temperature_piston_velocity = $r', &
realArgs=(/config_isomip_temperature_piston_velocity/))
            call mpas_log_write('        config_isomip_restoring_salinity = $r', realArgs=(/config_isomip_restoring_salinity/))
            call mpas_log_write('        config_isomip_salinity_piston_velocity = $r', &
realArgs=(/config_isomip_salinity_piston_velocity/))
            call mpas_log_write('        config_isomip_coriolis_parameter = $r', realArgs=(/config_isomip_coriolis_parameter/))
            call mpas_log_write('        config_isomip_southern_boundary = $r', realArgs=(/config_isomip_southern_boundary/))
            call mpas_log_write('        config_isomip_northern_boundary = $r', realArgs=(/config_isomip_northern_boundary/))
            call mpas_log_write('        config_isomip_western_boundary = $r', realArgs=(/config_isomip_western_boundary/))
            call mpas_log_write('        config_isomip_eastern_boundary = $r', realArgs=(/config_isomip_eastern_boundary/))
            call mpas_log_write('        config_isomip_y1 = $r', realArgs=(/config_isomip_y1/))
            call mpas_log_write('        config_isomip_z1 = $r', realArgs=(/config_isomip_z1/))
            call mpas_log_write('        config_isomip_ice_fraction1 = $r', realArgs=(/config_isomip_ice_fraction1/))
            call mpas_log_write('        config_isomip_y2 = $r', realArgs=(/config_isomip_y2/))
            call mpas_log_write('        config_isomip_z2 = $r', realArgs=(/config_isomip_z2/))
            call mpas_log_write('        config_isomip_ice_fraction2 = $r', realArgs=(/config_isomip_ice_fraction2/))
            call mpas_log_write('        config_isomip_y3 = $r', realArgs=(/config_isomip_y3/))
            call mpas_log_write('        config_isomip_z3 = $r', realArgs=(/config_isomip_z3/))
            call mpas_log_write('        config_isomip_ice_fraction3 = $r', realArgs=(/config_isomip_ice_fraction3/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record isomip.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_isomip_vert_levels', config_isomip_vert_levels)
      call mpas_pool_add_config(configPool, 'config_isomip_vertical_level_distribution', config_isomip_vertical_level_distribution)
      call mpas_pool_add_config(configPool, 'config_isomip_bottom_depth', config_isomip_bottom_depth)
      call mpas_pool_add_config(configPool, 'config_isomip_temperature', config_isomip_temperature)
      call mpas_pool_add_config(configPool, 'config_isomip_salinity', config_isomip_salinity)
      call mpas_pool_add_config(configPool, 'config_isomip_restoring_temperature', config_isomip_restoring_temperature)
      call mpas_pool_add_config(configPool, 'config_isomip_temperature_piston_velocity', config_isomip_temperature_piston_velocity)
      call mpas_pool_add_config(configPool, 'config_isomip_restoring_salinity', config_isomip_restoring_salinity)
      call mpas_pool_add_config(configPool, 'config_isomip_salinity_piston_velocity', config_isomip_salinity_piston_velocity)
      call mpas_pool_add_config(configPool, 'config_isomip_coriolis_parameter', config_isomip_coriolis_parameter)
      call mpas_pool_add_config(configPool, 'config_isomip_southern_boundary', config_isomip_southern_boundary)
      call mpas_pool_add_config(configPool, 'config_isomip_northern_boundary', config_isomip_northern_boundary)
      call mpas_pool_add_config(configPool, 'config_isomip_western_boundary', config_isomip_western_boundary)
      call mpas_pool_add_config(configPool, 'config_isomip_eastern_boundary', config_isomip_eastern_boundary)
      call mpas_pool_add_config(configPool, 'config_isomip_y1', config_isomip_y1)
      call mpas_pool_add_config(configPool, 'config_isomip_z1', config_isomip_z1)
      call mpas_pool_add_config(configPool, 'config_isomip_ice_fraction1', config_isomip_ice_fraction1)
      call mpas_pool_add_config(configPool, 'config_isomip_y2', config_isomip_y2)
      call mpas_pool_add_config(configPool, 'config_isomip_z2', config_isomip_z2)
      call mpas_pool_add_config(configPool, 'config_isomip_ice_fraction2', config_isomip_ice_fraction2)
      call mpas_pool_add_config(configPool, 'config_isomip_y3', config_isomip_y3)
      call mpas_pool_add_config(configPool, 'config_isomip_z3', config_isomip_z3)
      call mpas_pool_add_config(configPool, 'config_isomip_ice_fraction3', config_isomip_ice_fraction3)

   end subroutine ocn_setup_nmlrec_isomip


   subroutine ocn_setup_nmlrec_isomip_plus(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      integer :: config_isomip_plus_vert_levels = 36
      character (len=StrKIND) :: config_isomip_plus_vertical_level_distribution = 'constant'
      real (kind=RKIND) :: config_isomip_plus_max_bottom_depth = -720.000000
      integer :: config_isomip_plus_minimum_levels = 3
      real (kind=RKIND) :: config_isomip_plus_min_column_thickness = 10.000000
      real (kind=RKIND) :: config_isomip_plus_min_ocean_fraction = 0.500000
      character (len=StrKIND) :: config_isomip_plus_topography_file = 'input_geometry_processed.nc'
      real (kind=RKIND) :: config_isomip_plus_init_top_temp = -1.900000
      real (kind=RKIND) :: config_isomip_plus_init_bot_temp = -1.900000
      real (kind=RKIND) :: config_isomip_plus_init_top_sal = 33.800000
      real (kind=RKIND) :: config_isomip_plus_init_bot_sal = 34.500000
      real (kind=RKIND) :: config_isomip_plus_restore_top_temp = -1.900000
      real (kind=RKIND) :: config_isomip_plus_restore_bot_temp = 1.000000
      real (kind=RKIND) :: config_isomip_plus_restore_top_sal = 33.800000
      real (kind=RKIND) :: config_isomip_plus_restore_bot_sal = 34.700000
      real (kind=RKIND) :: config_isomip_plus_restore_rate = 10.000000
      real (kind=RKIND) :: config_isomip_plus_restore_evap_rate = 200.000000
      real (kind=RKIND) :: config_isomip_plus_restore_xMin = 790000.000000
      real (kind=RKIND) :: config_isomip_plus_restore_xMax = 800000.000000
      real (kind=RKIND) :: config_isomip_plus_coriolis_parameter = -0.000141
      real (kind=RKIND) :: config_isomip_plus_effective_density = 1026.000000

      namelist /isomip_plus/ &
         config_isomip_plus_vert_levels, &
         config_isomip_plus_vertical_level_distribution, &
         config_isomip_plus_max_bottom_depth, &
         config_isomip_plus_minimum_levels, &
         config_isomip_plus_min_column_thickness, &
         config_isomip_plus_min_ocean_fraction, &
         config_isomip_plus_topography_file, &
         config_isomip_plus_init_top_temp, &
         config_isomip_plus_init_bot_temp, &
         config_isomip_plus_init_top_sal, &
         config_isomip_plus_init_bot_sal, &
         config_isomip_plus_restore_top_temp, &
         config_isomip_plus_restore_bot_temp, &
         config_isomip_plus_restore_top_sal, &
         config_isomip_plus_restore_bot_sal, &
         config_isomip_plus_restore_rate, &
         config_isomip_plus_restore_evap_rate, &
         config_isomip_plus_restore_xMin, &
         config_isomip_plus_restore_xMax, &
         config_isomip_plus_coriolis_parameter, &
         config_isomip_plus_effective_density
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, isomip_plus, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_int(dminfo, config_isomip_plus_vert_levels)
         call mpas_dmpar_bcast_char(dminfo, config_isomip_plus_vertical_level_distribution)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_plus_max_bottom_depth)
         call mpas_dmpar_bcast_int(dminfo, config_isomip_plus_minimum_levels)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_plus_min_column_thickness)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_plus_min_ocean_fraction)
         call mpas_dmpar_bcast_char(dminfo, config_isomip_plus_topography_file)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_plus_init_top_temp)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_plus_init_bot_temp)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_plus_init_top_sal)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_plus_init_bot_sal)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_plus_restore_top_temp)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_plus_restore_bot_temp)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_plus_restore_top_sal)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_plus_restore_bot_sal)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_plus_restore_rate)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_plus_restore_evap_rate)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_plus_restore_xMin)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_plus_restore_xMax)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_plus_coriolis_parameter)
         call mpas_dmpar_bcast_real(dminfo, config_isomip_plus_effective_density)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record isomip_plus')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_isomip_plus_vert_levels = $i', intArgs=(/config_isomip_plus_vert_levels/))
            call mpas_log_write('        config_isomip_plus_vertical_level_distribution = '&
//''//mpas_log_escape_dollars(config_isomip_plus_vertical_level_distribution))
            call mpas_log_write('        config_isomip_plus_max_bottom_depth = $r', &
realArgs=(/config_isomip_plus_max_bottom_depth/))
            call mpas_log_write('        config_isomip_plus_minimum_levels = $i', intArgs=(/config_isomip_plus_minimum_levels/))
            call mpas_log_write('        config_isomip_plus_min_column_thickness = $r', &
realArgs=(/config_isomip_plus_min_column_thickness/))
            call mpas_log_write('        config_isomip_plus_min_ocean_fraction = $r', &
realArgs=(/config_isomip_plus_min_ocean_fraction/))
            call mpas_log_write('        config_isomip_plus_topography_file = '&
//''//mpas_log_escape_dollars(config_isomip_plus_topography_file))
            call mpas_log_write('        config_isomip_plus_init_top_temp = $r', realArgs=(/config_isomip_plus_init_top_temp/))
            call mpas_log_write('        config_isomip_plus_init_bot_temp = $r', realArgs=(/config_isomip_plus_init_bot_temp/))
            call mpas_log_write('        config_isomip_plus_init_top_sal = $r', realArgs=(/config_isomip_plus_init_top_sal/))
            call mpas_log_write('        config_isomip_plus_init_bot_sal = $r', realArgs=(/config_isomip_plus_init_bot_sal/))
            call mpas_log_write('        config_isomip_plus_restore_top_temp = $r', &
realArgs=(/config_isomip_plus_restore_top_temp/))
            call mpas_log_write('        config_isomip_plus_restore_bot_temp = $r', &
realArgs=(/config_isomip_plus_restore_bot_temp/))
            call mpas_log_write('        config_isomip_plus_restore_top_sal = $r', realArgs=(/config_isomip_plus_restore_top_sal/))
            call mpas_log_write('        config_isomip_plus_restore_bot_sal = $r', realArgs=(/config_isomip_plus_restore_bot_sal/))
            call mpas_log_write('        config_isomip_plus_restore_rate = $r', realArgs=(/config_isomip_plus_restore_rate/))
            call mpas_log_write('        config_isomip_plus_restore_evap_rate = $r', &
realArgs=(/config_isomip_plus_restore_evap_rate/))
            call mpas_log_write('        config_isomip_plus_restore_xMin = $r', realArgs=(/config_isomip_plus_restore_xMin/))
            call mpas_log_write('        config_isomip_plus_restore_xMax = $r', realArgs=(/config_isomip_plus_restore_xMax/))
            call mpas_log_write('        config_isomip_plus_coriolis_parameter = $r', &
realArgs=(/config_isomip_plus_coriolis_parameter/))
            call mpas_log_write('        config_isomip_plus_effective_density = $r', &
realArgs=(/config_isomip_plus_effective_density/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record isomip_plus.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_isomip_plus_vert_levels', config_isomip_plus_vert_levels)
      call mpas_pool_add_config(configPool, 'config_isomip_plus_vertical_level_distribution', &
config_isomip_plus_vertical_level_distribution)
      call mpas_pool_add_config(configPool, 'config_isomip_plus_max_bottom_depth', config_isomip_plus_max_bottom_depth)
      call mpas_pool_add_config(configPool, 'config_isomip_plus_minimum_levels', config_isomip_plus_minimum_levels)
      call mpas_pool_add_config(configPool, 'config_isomip_plus_min_column_thickness', config_isomip_plus_min_column_thickness)
      call mpas_pool_add_config(configPool, 'config_isomip_plus_min_ocean_fraction', config_isomip_plus_min_ocean_fraction)
      call mpas_pool_add_config(configPool, 'config_isomip_plus_topography_file', config_isomip_plus_topography_file)
      call mpas_pool_add_config(configPool, 'config_isomip_plus_init_top_temp', config_isomip_plus_init_top_temp)
      call mpas_pool_add_config(configPool, 'config_isomip_plus_init_bot_temp', config_isomip_plus_init_bot_temp)
      call mpas_pool_add_config(configPool, 'config_isomip_plus_init_top_sal', config_isomip_plus_init_top_sal)
      call mpas_pool_add_config(configPool, 'config_isomip_plus_init_bot_sal', config_isomip_plus_init_bot_sal)
      call mpas_pool_add_config(configPool, 'config_isomip_plus_restore_top_temp', config_isomip_plus_restore_top_temp)
      call mpas_pool_add_config(configPool, 'config_isomip_plus_restore_bot_temp', config_isomip_plus_restore_bot_temp)
      call mpas_pool_add_config(configPool, 'config_isomip_plus_restore_top_sal', config_isomip_plus_restore_top_sal)
      call mpas_pool_add_config(configPool, 'config_isomip_plus_restore_bot_sal', config_isomip_plus_restore_bot_sal)
      call mpas_pool_add_config(configPool, 'config_isomip_plus_restore_rate', config_isomip_plus_restore_rate)
      call mpas_pool_add_config(configPool, 'config_isomip_plus_restore_evap_rate', config_isomip_plus_restore_evap_rate)
      call mpas_pool_add_config(configPool, 'config_isomip_plus_restore_xMin', config_isomip_plus_restore_xMin)
      call mpas_pool_add_config(configPool, 'config_isomip_plus_restore_xMax', config_isomip_plus_restore_xMax)
      call mpas_pool_add_config(configPool, 'config_isomip_plus_coriolis_parameter', config_isomip_plus_coriolis_parameter)
      call mpas_pool_add_config(configPool, 'config_isomip_plus_effective_density', config_isomip_plus_effective_density)

   end subroutine ocn_setup_nmlrec_isomip_plus


   subroutine ocn_setup_nmlrec_hurricane(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      integer :: config_hurricane_vert_levels = 3
      real (kind=RKIND) :: config_hurricane_min_depth = 10.000000
      real (kind=RKIND) :: config_hurricane_max_depth = 60.000000
      real (kind=RKIND) :: config_hurricane_gaussian_hump_amplitude = 1.000000
      logical :: config_hurricane_use_gaussian_hump = .false.
      real (kind=RKIND) :: config_hurricane_gaussian_lon_center = 286.000000
      real (kind=RKIND) :: config_hurricane_gaussian_lat_center = 38.000000
      real (kind=RKIND) :: config_hurricane_gaussian_width = 1.000000
      real (kind=RKIND) :: config_hurricane_gaussian_amplitude = 1.000000
      real (kind=RKIND) :: config_hurricane_gaussian_slr_amp = 0.000000
      real (kind=RKIND) :: config_hurricane_land_z_limit = -2.000000
      real (kind=RKIND) :: config_hurricane_marsh_z_limit = 2.000000
      real (kind=RKIND) :: config_hurricane_land_drag = 0.100000
      real (kind=RKIND) :: config_hurricane_marsh_drag = 0.050000
      real (kind=RKIND) :: config_hurricane_channel_drag = 0.020000
      real (kind=RKIND) :: config_hurricane_sea_level_rise_adjustment = 0.000000

      namelist /hurricane/ &
         config_hurricane_vert_levels, &
         config_hurricane_min_depth, &
         config_hurricane_max_depth, &
         config_hurricane_gaussian_hump_amplitude, &
         config_hurricane_use_gaussian_hump, &
         config_hurricane_gaussian_lon_center, &
         config_hurricane_gaussian_lat_center, &
         config_hurricane_gaussian_width, &
         config_hurricane_gaussian_amplitude, &
         config_hurricane_gaussian_slr_amp, &
         config_hurricane_land_z_limit, &
         config_hurricane_marsh_z_limit, &
         config_hurricane_land_drag, &
         config_hurricane_marsh_drag, &
         config_hurricane_channel_drag, &
         config_hurricane_sea_level_rise_adjustment
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, hurricane, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_int(dminfo, config_hurricane_vert_levels)
         call mpas_dmpar_bcast_real(dminfo, config_hurricane_min_depth)
         call mpas_dmpar_bcast_real(dminfo, config_hurricane_max_depth)
         call mpas_dmpar_bcast_real(dminfo, config_hurricane_gaussian_hump_amplitude)
         call mpas_dmpar_bcast_logical(dminfo, config_hurricane_use_gaussian_hump)
         call mpas_dmpar_bcast_real(dminfo, config_hurricane_gaussian_lon_center)
         call mpas_dmpar_bcast_real(dminfo, config_hurricane_gaussian_lat_center)
         call mpas_dmpar_bcast_real(dminfo, config_hurricane_gaussian_width)
         call mpas_dmpar_bcast_real(dminfo, config_hurricane_gaussian_amplitude)
         call mpas_dmpar_bcast_real(dminfo, config_hurricane_gaussian_slr_amp)
         call mpas_dmpar_bcast_real(dminfo, config_hurricane_land_z_limit)
         call mpas_dmpar_bcast_real(dminfo, config_hurricane_marsh_z_limit)
         call mpas_dmpar_bcast_real(dminfo, config_hurricane_land_drag)
         call mpas_dmpar_bcast_real(dminfo, config_hurricane_marsh_drag)
         call mpas_dmpar_bcast_real(dminfo, config_hurricane_channel_drag)
         call mpas_dmpar_bcast_real(dminfo, config_hurricane_sea_level_rise_adjustment)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record hurricane')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_hurricane_vert_levels = $i', intArgs=(/config_hurricane_vert_levels/))
            call mpas_log_write('        config_hurricane_min_depth = $r', realArgs=(/config_hurricane_min_depth/))
            call mpas_log_write('        config_hurricane_max_depth = $r', realArgs=(/config_hurricane_max_depth/))
            call mpas_log_write('        config_hurricane_gaussian_hump_amplitude = $r', &
realArgs=(/config_hurricane_gaussian_hump_amplitude/))
            call mpas_log_write('        config_hurricane_use_gaussian_hump = $l', &
logicArgs=(/config_hurricane_use_gaussian_hump/))
            call mpas_log_write('        config_hurricane_gaussian_lon_center = $r', &
realArgs=(/config_hurricane_gaussian_lon_center/))
            call mpas_log_write('        config_hurricane_gaussian_lat_center = $r', &
realArgs=(/config_hurricane_gaussian_lat_center/))
            call mpas_log_write('        config_hurricane_gaussian_width = $r', realArgs=(/config_hurricane_gaussian_width/))
            call mpas_log_write('        config_hurricane_gaussian_amplitude = $r', &
realArgs=(/config_hurricane_gaussian_amplitude/))
            call mpas_log_write('        config_hurricane_gaussian_slr_amp = $r', realArgs=(/config_hurricane_gaussian_slr_amp/))
            call mpas_log_write('        config_hurricane_land_z_limit = $r', realArgs=(/config_hurricane_land_z_limit/))
            call mpas_log_write('        config_hurricane_marsh_z_limit = $r', realArgs=(/config_hurricane_marsh_z_limit/))
            call mpas_log_write('        config_hurricane_land_drag = $r', realArgs=(/config_hurricane_land_drag/))
            call mpas_log_write('        config_hurricane_marsh_drag = $r', realArgs=(/config_hurricane_marsh_drag/))
            call mpas_log_write('        config_hurricane_channel_drag = $r', realArgs=(/config_hurricane_channel_drag/))
            call mpas_log_write('        config_hurricane_sea_level_rise_adjustment = $r', &
realArgs=(/config_hurricane_sea_level_rise_adjustment/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record hurricane.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_hurricane_vert_levels', config_hurricane_vert_levels)
      call mpas_pool_add_config(configPool, 'config_hurricane_min_depth', config_hurricane_min_depth)
      call mpas_pool_add_config(configPool, 'config_hurricane_max_depth', config_hurricane_max_depth)
      call mpas_pool_add_config(configPool, 'config_hurricane_gaussian_hump_amplitude', config_hurricane_gaussian_hump_amplitude)
      call mpas_pool_add_config(configPool, 'config_hurricane_use_gaussian_hump', config_hurricane_use_gaussian_hump)
      call mpas_pool_add_config(configPool, 'config_hurricane_gaussian_lon_center', config_hurricane_gaussian_lon_center)
      call mpas_pool_add_config(configPool, 'config_hurricane_gaussian_lat_center', config_hurricane_gaussian_lat_center)
      call mpas_pool_add_config(configPool, 'config_hurricane_gaussian_width', config_hurricane_gaussian_width)
      call mpas_pool_add_config(configPool, 'config_hurricane_gaussian_amplitude', config_hurricane_gaussian_amplitude)
      call mpas_pool_add_config(configPool, 'config_hurricane_gaussian_slr_amp', config_hurricane_gaussian_slr_amp)
      call mpas_pool_add_config(configPool, 'config_hurricane_land_z_limit', config_hurricane_land_z_limit)
      call mpas_pool_add_config(configPool, 'config_hurricane_marsh_z_limit', config_hurricane_marsh_z_limit)
      call mpas_pool_add_config(configPool, 'config_hurricane_land_drag', config_hurricane_land_drag)
      call mpas_pool_add_config(configPool, 'config_hurricane_marsh_drag', config_hurricane_marsh_drag)
      call mpas_pool_add_config(configPool, 'config_hurricane_channel_drag', config_hurricane_channel_drag)
      call mpas_pool_add_config(configPool, 'config_hurricane_sea_level_rise_adjustment', &
config_hurricane_sea_level_rise_adjustment)

   end subroutine ocn_setup_nmlrec_hurricane


   subroutine ocn_setup_nmlrec_tidal_boundary(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      integer :: config_tidal_boundary_vert_levels = 100
      integer :: config_tidal_boundary_min_vert_levels = 10
      character (len=StrKIND) :: config_tidal_boundary_layer_type = 'zstar'
      real (kind=RKIND) :: config_tidal_boundary_right_bottom_depth = 10.000000
      logical :: config_tidal_start_dry = .false.
      logical :: config_tidal_boundary_use_distances = .true.
      real (kind=RKIND) :: config_tidal_boundary_left_value = 0.000000
      real (kind=RKIND) :: config_tidal_boundary_right_value = 25000.000000
      real (kind=RKIND) :: config_tidal_boundary_left_bottom_depth = 10.000000
      real (kind=RKIND) :: config_tidal_boundary_salinity = 35.000000
      real (kind=RKIND) :: config_tidal_boundary_domain_temperature = 20.000000
      real (kind=RKIND) :: config_tidal_boundary_plug_temperature = 20.000000
      real (kind=RKIND) :: config_tidal_boundary_plug_width_frac = 0.100000
      real (kind=RKIND) :: config_tidal_forcing_left_Cd_or_n = 0.001000
      real (kind=RKIND) :: config_tidal_forcing_right_Cd_or_n = 0.001000
      logical :: config_use_idealized_transect = .false.
      real (kind=RKIND) :: config_idealized_transect_Lshore = 0.600000
      real (kind=RKIND) :: config_idealized_transect_Sshore = 0.001000
      real (kind=RKIND) :: config_idealized_transect_Lcoast = 0.300000
      real (kind=RKIND) :: config_idealized_transect_Scoast = 0.001000
      real (kind=RKIND) :: config_idealized_transect_Lmarsh = 0.100000
      real (kind=RKIND) :: config_idealized_transect_Smarsh = 0.000000
      real (kind=RKIND) :: config_idealized_transect_roughness = 0.025000
      real (kind=RKIND) :: config_idealized_transect_roughness_marsh = 0.075000
      real (kind=RKIND) :: config_idealized_vegetation_diameter = 0.050000
      real (kind=RKIND) :: config_idealized_vegetation_height = 0.200000
      real (kind=RKIND) :: config_idealized_vegetation_density = 1000.000000

      namelist /tidal_boundary/ &
         config_tidal_boundary_vert_levels, &
         config_tidal_boundary_min_vert_levels, &
         config_tidal_boundary_layer_type, &
         config_tidal_boundary_right_bottom_depth, &
         config_tidal_start_dry, &
         config_tidal_boundary_use_distances, &
         config_tidal_boundary_left_value, &
         config_tidal_boundary_right_value, &
         config_tidal_boundary_left_bottom_depth, &
         config_tidal_boundary_salinity, &
         config_tidal_boundary_domain_temperature, &
         config_tidal_boundary_plug_temperature, &
         config_tidal_boundary_plug_width_frac, &
         config_tidal_forcing_left_Cd_or_n, &
         config_tidal_forcing_right_Cd_or_n, &
         config_use_idealized_transect, &
         config_idealized_transect_Lshore, &
         config_idealized_transect_Sshore, &
         config_idealized_transect_Lcoast, &
         config_idealized_transect_Scoast, &
         config_idealized_transect_Lmarsh, &
         config_idealized_transect_Smarsh, &
         config_idealized_transect_roughness, &
         config_idealized_transect_roughness_marsh, &
         config_idealized_vegetation_diameter, &
         config_idealized_vegetation_height, &
         config_idealized_vegetation_density
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, tidal_boundary, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_int(dminfo, config_tidal_boundary_vert_levels)
         call mpas_dmpar_bcast_int(dminfo, config_tidal_boundary_min_vert_levels)
         call mpas_dmpar_bcast_char(dminfo, config_tidal_boundary_layer_type)
         call mpas_dmpar_bcast_real(dminfo, config_tidal_boundary_right_bottom_depth)
         call mpas_dmpar_bcast_logical(dminfo, config_tidal_start_dry)
         call mpas_dmpar_bcast_logical(dminfo, config_tidal_boundary_use_distances)
         call mpas_dmpar_bcast_real(dminfo, config_tidal_boundary_left_value)
         call mpas_dmpar_bcast_real(dminfo, config_tidal_boundary_right_value)
         call mpas_dmpar_bcast_real(dminfo, config_tidal_boundary_left_bottom_depth)
         call mpas_dmpar_bcast_real(dminfo, config_tidal_boundary_salinity)
         call mpas_dmpar_bcast_real(dminfo, config_tidal_boundary_domain_temperature)
         call mpas_dmpar_bcast_real(dminfo, config_tidal_boundary_plug_temperature)
         call mpas_dmpar_bcast_real(dminfo, config_tidal_boundary_plug_width_frac)
         call mpas_dmpar_bcast_real(dminfo, config_tidal_forcing_left_Cd_or_n)
         call mpas_dmpar_bcast_real(dminfo, config_tidal_forcing_right_Cd_or_n)
         call mpas_dmpar_bcast_logical(dminfo, config_use_idealized_transect)
         call mpas_dmpar_bcast_real(dminfo, config_idealized_transect_Lshore)
         call mpas_dmpar_bcast_real(dminfo, config_idealized_transect_Sshore)
         call mpas_dmpar_bcast_real(dminfo, config_idealized_transect_Lcoast)
         call mpas_dmpar_bcast_real(dminfo, config_idealized_transect_Scoast)
         call mpas_dmpar_bcast_real(dminfo, config_idealized_transect_Lmarsh)
         call mpas_dmpar_bcast_real(dminfo, config_idealized_transect_Smarsh)
         call mpas_dmpar_bcast_real(dminfo, config_idealized_transect_roughness)
         call mpas_dmpar_bcast_real(dminfo, config_idealized_transect_roughness_marsh)
         call mpas_dmpar_bcast_real(dminfo, config_idealized_vegetation_diameter)
         call mpas_dmpar_bcast_real(dminfo, config_idealized_vegetation_height)
         call mpas_dmpar_bcast_real(dminfo, config_idealized_vegetation_density)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record tidal_boundary')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_tidal_boundary_vert_levels = $i', intArgs=(/config_tidal_boundary_vert_levels/))
            call mpas_log_write('        config_tidal_boundary_min_vert_levels = $i', &
intArgs=(/config_tidal_boundary_min_vert_levels/))
            call mpas_log_write('        config_tidal_boundary_layer_type = '&
//''//mpas_log_escape_dollars(config_tidal_boundary_layer_type))
            call mpas_log_write('        config_tidal_boundary_right_bottom_depth = $r', &
realArgs=(/config_tidal_boundary_right_bottom_depth/))
            call mpas_log_write('        config_tidal_start_dry = $l', logicArgs=(/config_tidal_start_dry/))
            call mpas_log_write('        config_tidal_boundary_use_distances = $l', &
logicArgs=(/config_tidal_boundary_use_distances/))
            call mpas_log_write('        config_tidal_boundary_left_value = $r', realArgs=(/config_tidal_boundary_left_value/))
            call mpas_log_write('        config_tidal_boundary_right_value = $r', realArgs=(/config_tidal_boundary_right_value/))
            call mpas_log_write('        config_tidal_boundary_left_bottom_depth = $r', &
realArgs=(/config_tidal_boundary_left_bottom_depth/))
            call mpas_log_write('        config_tidal_boundary_salinity = $r', realArgs=(/config_tidal_boundary_salinity/))
            call mpas_log_write('        config_tidal_boundary_domain_temperature = $r', &
realArgs=(/config_tidal_boundary_domain_temperature/))
            call mpas_log_write('        config_tidal_boundary_plug_temperature = $r', &
realArgs=(/config_tidal_boundary_plug_temperature/))
            call mpas_log_write('        config_tidal_boundary_plug_width_frac = $r', &
realArgs=(/config_tidal_boundary_plug_width_frac/))
            call mpas_log_write('        config_tidal_forcing_left_Cd_or_n = $r', realArgs=(/config_tidal_forcing_left_Cd_or_n/))
            call mpas_log_write('        config_tidal_forcing_right_Cd_or_n = $r', realArgs=(/config_tidal_forcing_right_Cd_or_n/))
            call mpas_log_write('        config_use_idealized_transect = $l', logicArgs=(/config_use_idealized_transect/))
            call mpas_log_write('        config_idealized_transect_Lshore = $r', realArgs=(/config_idealized_transect_Lshore/))
            call mpas_log_write('        config_idealized_transect_Sshore = $r', realArgs=(/config_idealized_transect_Sshore/))
            call mpas_log_write('        config_idealized_transect_Lcoast = $r', realArgs=(/config_idealized_transect_Lcoast/))
            call mpas_log_write('        config_idealized_transect_Scoast = $r', realArgs=(/config_idealized_transect_Scoast/))
            call mpas_log_write('        config_idealized_transect_Lmarsh = $r', realArgs=(/config_idealized_transect_Lmarsh/))
            call mpas_log_write('        config_idealized_transect_Smarsh = $r', realArgs=(/config_idealized_transect_Smarsh/))
            call mpas_log_write('        config_idealized_transect_roughness = $r', &
realArgs=(/config_idealized_transect_roughness/))
            call mpas_log_write('        config_idealized_transect_roughness_marsh = $r', &
realArgs=(/config_idealized_transect_roughness_marsh/))
            call mpas_log_write('        config_idealized_vegetation_diameter = $r', &
realArgs=(/config_idealized_vegetation_diameter/))
            call mpas_log_write('        config_idealized_vegetation_height = $r', realArgs=(/config_idealized_vegetation_height/))
            call mpas_log_write('        config_idealized_vegetation_density = $r', &
realArgs=(/config_idealized_vegetation_density/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record tidal_boundary.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_tidal_boundary_vert_levels', config_tidal_boundary_vert_levels)
      call mpas_pool_add_config(configPool, 'config_tidal_boundary_min_vert_levels', config_tidal_boundary_min_vert_levels)
      call mpas_pool_add_config(configPool, 'config_tidal_boundary_layer_type', config_tidal_boundary_layer_type)
      call mpas_pool_add_config(configPool, 'config_tidal_boundary_right_bottom_depth', config_tidal_boundary_right_bottom_depth)
      call mpas_pool_add_config(configPool, 'config_tidal_start_dry', config_tidal_start_dry)
      call mpas_pool_add_config(configPool, 'config_tidal_boundary_use_distances', config_tidal_boundary_use_distances)
      call mpas_pool_add_config(configPool, 'config_tidal_boundary_left_value', config_tidal_boundary_left_value)
      call mpas_pool_add_config(configPool, 'config_tidal_boundary_right_value', config_tidal_boundary_right_value)
      call mpas_pool_add_config(configPool, 'config_tidal_boundary_left_bottom_depth', config_tidal_boundary_left_bottom_depth)
      call mpas_pool_add_config(configPool, 'config_tidal_boundary_salinity', config_tidal_boundary_salinity)
      call mpas_pool_add_config(configPool, 'config_tidal_boundary_domain_temperature', config_tidal_boundary_domain_temperature)
      call mpas_pool_add_config(configPool, 'config_tidal_boundary_plug_temperature', config_tidal_boundary_plug_temperature)
      call mpas_pool_add_config(configPool, 'config_tidal_boundary_plug_width_frac', config_tidal_boundary_plug_width_frac)
      call mpas_pool_add_config(configPool, 'config_tidal_forcing_left_Cd_or_n', config_tidal_forcing_left_Cd_or_n)
      call mpas_pool_add_config(configPool, 'config_tidal_forcing_right_Cd_or_n', config_tidal_forcing_right_Cd_or_n)
      call mpas_pool_add_config(configPool, 'config_use_idealized_transect', config_use_idealized_transect)
      call mpas_pool_add_config(configPool, 'config_idealized_transect_Lshore', config_idealized_transect_Lshore)
      call mpas_pool_add_config(configPool, 'config_idealized_transect_Sshore', config_idealized_transect_Sshore)
      call mpas_pool_add_config(configPool, 'config_idealized_transect_Lcoast', config_idealized_transect_Lcoast)
      call mpas_pool_add_config(configPool, 'config_idealized_transect_Scoast', config_idealized_transect_Scoast)
      call mpas_pool_add_config(configPool, 'config_idealized_transect_Lmarsh', config_idealized_transect_Lmarsh)
      call mpas_pool_add_config(configPool, 'config_idealized_transect_Smarsh', config_idealized_transect_Smarsh)
      call mpas_pool_add_config(configPool, 'config_idealized_transect_roughness', config_idealized_transect_roughness)
      call mpas_pool_add_config(configPool, 'config_idealized_transect_roughness_marsh', config_idealized_transect_roughness_marsh)
      call mpas_pool_add_config(configPool, 'config_idealized_vegetation_diameter', config_idealized_vegetation_diameter)
      call mpas_pool_add_config(configPool, 'config_idealized_vegetation_height', config_idealized_vegetation_height)
      call mpas_pool_add_config(configPool, 'config_idealized_vegetation_density', config_idealized_vegetation_density)

   end subroutine ocn_setup_nmlrec_tidal_boundary


   subroutine ocn_setup_nmlrec_cosine_bell(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      real (kind=RKIND) :: config_cosine_bell_temperature = 15.000000
      real (kind=RKIND) :: config_cosine_bell_salinity = 35.000000
      real (kind=RKIND) :: config_cosine_bell_lat_center = 0.000000
      real (kind=RKIND) :: config_cosine_bell_lon_center = 3.141592
      real (kind=RKIND) :: config_cosine_bell_psi0 = 1.000000
      real (kind=RKIND) :: config_cosine_bell_radius = 2123666.667000
      real (kind=RKIND) :: config_cosine_bell_vel_pd = 24.000000

      namelist /cosine_bell/ &
         config_cosine_bell_temperature, &
         config_cosine_bell_salinity, &
         config_cosine_bell_lat_center, &
         config_cosine_bell_lon_center, &
         config_cosine_bell_psi0, &
         config_cosine_bell_radius, &
         config_cosine_bell_vel_pd
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, cosine_bell, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_real(dminfo, config_cosine_bell_temperature)
         call mpas_dmpar_bcast_real(dminfo, config_cosine_bell_salinity)
         call mpas_dmpar_bcast_real(dminfo, config_cosine_bell_lat_center)
         call mpas_dmpar_bcast_real(dminfo, config_cosine_bell_lon_center)
         call mpas_dmpar_bcast_real(dminfo, config_cosine_bell_psi0)
         call mpas_dmpar_bcast_real(dminfo, config_cosine_bell_radius)
         call mpas_dmpar_bcast_real(dminfo, config_cosine_bell_vel_pd)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record cosine_bell')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_cosine_bell_temperature = $r', realArgs=(/config_cosine_bell_temperature/))
            call mpas_log_write('        config_cosine_bell_salinity = $r', realArgs=(/config_cosine_bell_salinity/))
            call mpas_log_write('        config_cosine_bell_lat_center = $r', realArgs=(/config_cosine_bell_lat_center/))
            call mpas_log_write('        config_cosine_bell_lon_center = $r', realArgs=(/config_cosine_bell_lon_center/))
            call mpas_log_write('        config_cosine_bell_psi0 = $r', realArgs=(/config_cosine_bell_psi0/))
            call mpas_log_write('        config_cosine_bell_radius = $r', realArgs=(/config_cosine_bell_radius/))
            call mpas_log_write('        config_cosine_bell_vel_pd = $r', realArgs=(/config_cosine_bell_vel_pd/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record cosine_bell.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_cosine_bell_temperature', config_cosine_bell_temperature)
      call mpas_pool_add_config(configPool, 'config_cosine_bell_salinity', config_cosine_bell_salinity)
      call mpas_pool_add_config(configPool, 'config_cosine_bell_lat_center', config_cosine_bell_lat_center)
      call mpas_pool_add_config(configPool, 'config_cosine_bell_lon_center', config_cosine_bell_lon_center)
      call mpas_pool_add_config(configPool, 'config_cosine_bell_psi0', config_cosine_bell_psi0)
      call mpas_pool_add_config(configPool, 'config_cosine_bell_radius', config_cosine_bell_radius)
      call mpas_pool_add_config(configPool, 'config_cosine_bell_vel_pd', config_cosine_bell_vel_pd)

   end subroutine ocn_setup_nmlrec_cosine_bell


   subroutine ocn_setup_nmlrec_mixed_layer_eddy(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      integer :: config_mixed_layer_eddy_vert_levels = 60
      real (kind=RKIND) :: config_mixed_layer_eddy_bottom_depth = 300.000000
      real (kind=RKIND) :: config_mixed_layer_eddy_mixed_layer_depth = 200.000000
      real (kind=RKIND) :: config_mixed_layer_eddy_base_temperature = 16.000000
      real (kind=RKIND) :: config_mixed_layer_eddy_temperature_stratification_mixed_layer = 0.000100
      real (kind=RKIND) :: config_mixed_layer_eddy_temperature_stratification_interior = 0.010000
      real (kind=RKIND) :: config_mixed_layer_eddy_temperature_horizontal_gradient = 0.000020
      real (kind=RKIND) :: config_mixed_layer_eddy_temperature_front_width = 10000.000000
      real (kind=RKIND) :: config_mixed_layer_eddy_temperature_perturbation_magnitude = 0.000010
      real (kind=RKIND) :: config_mixed_layer_eddy_salinity = 35.000000
      logical :: config_mixed_layer_eddy_two_fronts = .false.
      real (kind=RKIND) :: config_mixed_layer_eddy_restoring_width = 5000.000000
      real (kind=RKIND) :: config_mixed_layer_eddy_restoring_tau = 5.000000
      real (kind=RKIND) :: config_mixed_layer_eddy_heat_flux = 0.000000
      real (kind=RKIND) :: config_mixed_layer_eddy_evaporation_flux = 0.000000
      real (kind=RKIND) :: config_mixed_layer_eddy_wind_stress_zonal = 0.000000
      real (kind=RKIND) :: config_mixed_layer_eddy_wind_stress_meridional = 0.000000
      real (kind=RKIND) :: config_mixed_layer_eddy_coriolis_parameter = 0.000100

      namelist /mixed_layer_eddy/ &
         config_mixed_layer_eddy_vert_levels, &
         config_mixed_layer_eddy_bottom_depth, &
         config_mixed_layer_eddy_mixed_layer_depth, &
         config_mixed_layer_eddy_base_temperature, &
         config_mixed_layer_eddy_temperature_stratification_mixed_layer, &
         config_mixed_layer_eddy_temperature_stratification_interior, &
         config_mixed_layer_eddy_temperature_horizontal_gradient, &
         config_mixed_layer_eddy_temperature_front_width, &
         config_mixed_layer_eddy_temperature_perturbation_magnitude, &
         config_mixed_layer_eddy_salinity, &
         config_mixed_layer_eddy_two_fronts, &
         config_mixed_layer_eddy_restoring_width, &
         config_mixed_layer_eddy_restoring_tau, &
         config_mixed_layer_eddy_heat_flux, &
         config_mixed_layer_eddy_evaporation_flux, &
         config_mixed_layer_eddy_wind_stress_zonal, &
         config_mixed_layer_eddy_wind_stress_meridional, &
         config_mixed_layer_eddy_coriolis_parameter
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, mixed_layer_eddy, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_int(dminfo, config_mixed_layer_eddy_vert_levels)
         call mpas_dmpar_bcast_real(dminfo, config_mixed_layer_eddy_bottom_depth)
         call mpas_dmpar_bcast_real(dminfo, config_mixed_layer_eddy_mixed_layer_depth)
         call mpas_dmpar_bcast_real(dminfo, config_mixed_layer_eddy_base_temperature)
         call mpas_dmpar_bcast_real(dminfo, config_mixed_layer_eddy_temperature_stratification_mixed_layer)
         call mpas_dmpar_bcast_real(dminfo, config_mixed_layer_eddy_temperature_stratification_interior)
         call mpas_dmpar_bcast_real(dminfo, config_mixed_layer_eddy_temperature_horizontal_gradient)
         call mpas_dmpar_bcast_real(dminfo, config_mixed_layer_eddy_temperature_front_width)
         call mpas_dmpar_bcast_real(dminfo, config_mixed_layer_eddy_temperature_perturbation_magnitude)
         call mpas_dmpar_bcast_real(dminfo, config_mixed_layer_eddy_salinity)
         call mpas_dmpar_bcast_logical(dminfo, config_mixed_layer_eddy_two_fronts)
         call mpas_dmpar_bcast_real(dminfo, config_mixed_layer_eddy_restoring_width)
         call mpas_dmpar_bcast_real(dminfo, config_mixed_layer_eddy_restoring_tau)
         call mpas_dmpar_bcast_real(dminfo, config_mixed_layer_eddy_heat_flux)
         call mpas_dmpar_bcast_real(dminfo, config_mixed_layer_eddy_evaporation_flux)
         call mpas_dmpar_bcast_real(dminfo, config_mixed_layer_eddy_wind_stress_zonal)
         call mpas_dmpar_bcast_real(dminfo, config_mixed_layer_eddy_wind_stress_meridional)
         call mpas_dmpar_bcast_real(dminfo, config_mixed_layer_eddy_coriolis_parameter)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record mixed_layer_eddy')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_mixed_layer_eddy_vert_levels = $i', &
intArgs=(/config_mixed_layer_eddy_vert_levels/))
            call mpas_log_write('        config_mixed_layer_eddy_bottom_depth = $r', &
realArgs=(/config_mixed_layer_eddy_bottom_depth/))
            call mpas_log_write('        config_mixed_layer_eddy_mixed_layer_depth = $r', &
realArgs=(/config_mixed_layer_eddy_mixed_layer_depth/))
            call mpas_log_write('        config_mixed_layer_eddy_base_temperature = $r', &
realArgs=(/config_mixed_layer_eddy_base_temperature/))
            call mpas_log_write('        config_mixed_layer_eddy_temperature_stratification_mixed_layer = $r', &
realArgs=(/config_mixed_layer_eddy_temperature_stratification_mixed_layer/))
            call mpas_log_write('        config_mixed_layer_eddy_temperature_stratification_interior = $r', &
realArgs=(/config_mixed_layer_eddy_temperature_stratification_interior/))
            call mpas_log_write('        config_mixed_layer_eddy_temperature_horizontal_gradient = $r', &
realArgs=(/config_mixed_layer_eddy_temperature_horizontal_gradient/))
            call mpas_log_write('        config_mixed_layer_eddy_temperature_front_width = $r', &
realArgs=(/config_mixed_layer_eddy_temperature_front_width/))
            call mpas_log_write('        config_mixed_layer_eddy_temperature_perturbation_magnitude = $r', &
realArgs=(/config_mixed_layer_eddy_temperature_perturbation_magnitude/))
            call mpas_log_write('        config_mixed_layer_eddy_salinity = $r', realArgs=(/config_mixed_layer_eddy_salinity/))
            call mpas_log_write('        config_mixed_layer_eddy_two_fronts = $l', &
logicArgs=(/config_mixed_layer_eddy_two_fronts/))
            call mpas_log_write('        config_mixed_layer_eddy_restoring_width = $r', &
realArgs=(/config_mixed_layer_eddy_restoring_width/))
            call mpas_log_write('        config_mixed_layer_eddy_restoring_tau = $r', &
realArgs=(/config_mixed_layer_eddy_restoring_tau/))
            call mpas_log_write('        config_mixed_layer_eddy_heat_flux = $r', realArgs=(/config_mixed_layer_eddy_heat_flux/))
            call mpas_log_write('        config_mixed_layer_eddy_evaporation_flux = $r', &
realArgs=(/config_mixed_layer_eddy_evaporation_flux/))
            call mpas_log_write('        config_mixed_layer_eddy_wind_stress_zonal = $r', &
realArgs=(/config_mixed_layer_eddy_wind_stress_zonal/))
            call mpas_log_write('        config_mixed_layer_eddy_wind_stress_meridional = $r', &
realArgs=(/config_mixed_layer_eddy_wind_stress_meridional/))
            call mpas_log_write('        config_mixed_layer_eddy_coriolis_parameter = $r', &
realArgs=(/config_mixed_layer_eddy_coriolis_parameter/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record mixed_layer_eddy.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_mixed_layer_eddy_vert_levels', config_mixed_layer_eddy_vert_levels)
      call mpas_pool_add_config(configPool, 'config_mixed_layer_eddy_bottom_depth', config_mixed_layer_eddy_bottom_depth)
      call mpas_pool_add_config(configPool, 'config_mixed_layer_eddy_mixed_layer_depth', config_mixed_layer_eddy_mixed_layer_depth)
      call mpas_pool_add_config(configPool, 'config_mixed_layer_eddy_base_temperature', config_mixed_layer_eddy_base_temperature)
      call mpas_pool_add_config(configPool, 'config_mixed_layer_eddy_temperature_stratification_mixed_layer', &
config_mixed_layer_eddy_temperature_stratification_mixed_layer)
      call mpas_pool_add_config(configPool, 'config_mixed_layer_eddy_temperature_stratification_interior', &
config_mixed_layer_eddy_temperature_stratification_interior)
      call mpas_pool_add_config(configPool, 'config_mixed_layer_eddy_temperature_horizontal_gradient', &
config_mixed_layer_eddy_temperature_horizontal_gradient)
      call mpas_pool_add_config(configPool, 'config_mixed_layer_eddy_temperature_front_width', &
config_mixed_layer_eddy_temperature_front_width)
      call mpas_pool_add_config(configPool, 'config_mixed_layer_eddy_temperature_perturbation_magnitude', &
config_mixed_layer_eddy_temperature_perturbation_magnitude)
      call mpas_pool_add_config(configPool, 'config_mixed_layer_eddy_salinity', config_mixed_layer_eddy_salinity)
      call mpas_pool_add_config(configPool, 'config_mixed_layer_eddy_two_fronts', config_mixed_layer_eddy_two_fronts)
      call mpas_pool_add_config(configPool, 'config_mixed_layer_eddy_restoring_width', config_mixed_layer_eddy_restoring_width)
      call mpas_pool_add_config(configPool, 'config_mixed_layer_eddy_restoring_tau', config_mixed_layer_eddy_restoring_tau)
      call mpas_pool_add_config(configPool, 'config_mixed_layer_eddy_heat_flux', config_mixed_layer_eddy_heat_flux)
      call mpas_pool_add_config(configPool, 'config_mixed_layer_eddy_evaporation_flux', config_mixed_layer_eddy_evaporation_flux)
      call mpas_pool_add_config(configPool, 'config_mixed_layer_eddy_wind_stress_zonal', config_mixed_layer_eddy_wind_stress_zonal)
      call mpas_pool_add_config(configPool, 'config_mixed_layer_eddy_wind_stress_meridional', &
config_mixed_layer_eddy_wind_stress_meridional)
      call mpas_pool_add_config(configPool, 'config_mixed_layer_eddy_coriolis_parameter', &
config_mixed_layer_eddy_coriolis_parameter)

   end subroutine ocn_setup_nmlrec_mixed_layer_eddy


   subroutine ocn_setup_nmlrec_test_sht(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      integer :: config_test_sht_function_option = 1
      real (kind=RKIND) :: config_test_sht_cosine_bell_lat_center = 0.000000
      real (kind=RKIND) :: config_test_sht_cosine_bell_lon_center = 3.141592
      real (kind=RKIND) :: config_test_sht_cosine_bell_psi0 = 1.000000
      real (kind=RKIND) :: config_test_sht_cosine_bell_radius = 2123666.667000
      real (kind=RKIND) :: config_test_sht_function3_cell_width_equator = 30.000000
      real (kind=RKIND) :: config_test_sht_function3_cell_width_pole = 120.000000
      real (kind=RKIND) :: config_test_sht_function3_lat_transition_start = 1.000000
      real (kind=RKIND) :: config_test_sht_function3_lat_transition_width = 10.000000
      integer :: config_test_sht_n_iterations = 10

      namelist /test_sht/ &
         config_test_sht_function_option, &
         config_test_sht_cosine_bell_lat_center, &
         config_test_sht_cosine_bell_lon_center, &
         config_test_sht_cosine_bell_psi0, &
         config_test_sht_cosine_bell_radius, &
         config_test_sht_function3_cell_width_equator, &
         config_test_sht_function3_cell_width_pole, &
         config_test_sht_function3_lat_transition_start, &
         config_test_sht_function3_lat_transition_width, &
         config_test_sht_n_iterations
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, test_sht, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_int(dminfo, config_test_sht_function_option)
         call mpas_dmpar_bcast_real(dminfo, config_test_sht_cosine_bell_lat_center)
         call mpas_dmpar_bcast_real(dminfo, config_test_sht_cosine_bell_lon_center)
         call mpas_dmpar_bcast_real(dminfo, config_test_sht_cosine_bell_psi0)
         call mpas_dmpar_bcast_real(dminfo, config_test_sht_cosine_bell_radius)
         call mpas_dmpar_bcast_real(dminfo, config_test_sht_function3_cell_width_equator)
         call mpas_dmpar_bcast_real(dminfo, config_test_sht_function3_cell_width_pole)
         call mpas_dmpar_bcast_real(dminfo, config_test_sht_function3_lat_transition_start)
         call mpas_dmpar_bcast_real(dminfo, config_test_sht_function3_lat_transition_width)
         call mpas_dmpar_bcast_int(dminfo, config_test_sht_n_iterations)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record test_sht')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_test_sht_function_option = $i', intArgs=(/config_test_sht_function_option/))
            call mpas_log_write('        config_test_sht_cosine_bell_lat_center = $r', &
realArgs=(/config_test_sht_cosine_bell_lat_center/))
            call mpas_log_write('        config_test_sht_cosine_bell_lon_center = $r', &
realArgs=(/config_test_sht_cosine_bell_lon_center/))
            call mpas_log_write('        config_test_sht_cosine_bell_psi0 = $r', realArgs=(/config_test_sht_cosine_bell_psi0/))
            call mpas_log_write('        config_test_sht_cosine_bell_radius = $r', realArgs=(/config_test_sht_cosine_bell_radius/))
            call mpas_log_write('        config_test_sht_function3_cell_width_equator = $r', &
realArgs=(/config_test_sht_function3_cell_width_equator/))
            call mpas_log_write('        config_test_sht_function3_cell_width_pole = $r', &
realArgs=(/config_test_sht_function3_cell_width_pole/))
            call mpas_log_write('        config_test_sht_function3_lat_transition_start = $r', &
realArgs=(/config_test_sht_function3_lat_transition_start/))
            call mpas_log_write('        config_test_sht_function3_lat_transition_width = $r', &
realArgs=(/config_test_sht_function3_lat_transition_width/))
            call mpas_log_write('        config_test_sht_n_iterations = $i', intArgs=(/config_test_sht_n_iterations/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record test_sht.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_test_sht_function_option', config_test_sht_function_option)
      call mpas_pool_add_config(configPool, 'config_test_sht_cosine_bell_lat_center', config_test_sht_cosine_bell_lat_center)
      call mpas_pool_add_config(configPool, 'config_test_sht_cosine_bell_lon_center', config_test_sht_cosine_bell_lon_center)
      call mpas_pool_add_config(configPool, 'config_test_sht_cosine_bell_psi0', config_test_sht_cosine_bell_psi0)
      call mpas_pool_add_config(configPool, 'config_test_sht_cosine_bell_radius', config_test_sht_cosine_bell_radius)
      call mpas_pool_add_config(configPool, 'config_test_sht_function3_cell_width_equator', &
config_test_sht_function3_cell_width_equator)
      call mpas_pool_add_config(configPool, 'config_test_sht_function3_cell_width_pole', config_test_sht_function3_cell_width_pole)
      call mpas_pool_add_config(configPool, 'config_test_sht_function3_lat_transition_start', &
config_test_sht_function3_lat_transition_start)
      call mpas_pool_add_config(configPool, 'config_test_sht_function3_lat_transition_width', &
config_test_sht_function3_lat_transition_width)
      call mpas_pool_add_config(configPool, 'config_test_sht_n_iterations', config_test_sht_n_iterations)

   end subroutine ocn_setup_nmlrec_test_sht


   subroutine ocn_setup_nmlrec_parabolic_bowl(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      integer :: config_parabolic_bowl_vert_levels = 3
      real (kind=RKIND) :: config_parabolic_bowl_Coriolis_parameter = 0.000103
      real (kind=RKIND) :: config_parabolic_bowl_eta0 = 2.000000
      real (kind=RKIND) :: config_parabolic_bowl_b0 = 50.000000
      real (kind=RKIND) :: config_parabolic_bowl_omega = 0.000145
      real (kind=RKIND) :: config_parabolic_bowl_gravity = 9.810000
      logical :: config_parabolic_bowl_adjust_domain_center = .true.

      namelist /parabolic_bowl/ &
         config_parabolic_bowl_vert_levels, &
         config_parabolic_bowl_Coriolis_parameter, &
         config_parabolic_bowl_eta0, &
         config_parabolic_bowl_b0, &
         config_parabolic_bowl_omega, &
         config_parabolic_bowl_gravity, &
         config_parabolic_bowl_adjust_domain_center
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, parabolic_bowl, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_int(dminfo, config_parabolic_bowl_vert_levels)
         call mpas_dmpar_bcast_real(dminfo, config_parabolic_bowl_Coriolis_parameter)
         call mpas_dmpar_bcast_real(dminfo, config_parabolic_bowl_eta0)
         call mpas_dmpar_bcast_real(dminfo, config_parabolic_bowl_b0)
         call mpas_dmpar_bcast_real(dminfo, config_parabolic_bowl_omega)
         call mpas_dmpar_bcast_real(dminfo, config_parabolic_bowl_gravity)
         call mpas_dmpar_bcast_logical(dminfo, config_parabolic_bowl_adjust_domain_center)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record parabolic_bowl')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_parabolic_bowl_vert_levels = $i', intArgs=(/config_parabolic_bowl_vert_levels/))
            call mpas_log_write('        config_parabolic_bowl_Coriolis_parameter = $r', &
realArgs=(/config_parabolic_bowl_Coriolis_parameter/))
            call mpas_log_write('        config_parabolic_bowl_eta0 = $r', realArgs=(/config_parabolic_bowl_eta0/))
            call mpas_log_write('        config_parabolic_bowl_b0 = $r', realArgs=(/config_parabolic_bowl_b0/))
            call mpas_log_write('        config_parabolic_bowl_omega = $r', realArgs=(/config_parabolic_bowl_omega/))
            call mpas_log_write('        config_parabolic_bowl_gravity = $r', realArgs=(/config_parabolic_bowl_gravity/))
            call mpas_log_write('        config_parabolic_bowl_adjust_domain_center = $l', &
logicArgs=(/config_parabolic_bowl_adjust_domain_center/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record parabolic_bowl.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_parabolic_bowl_vert_levels', config_parabolic_bowl_vert_levels)
      call mpas_pool_add_config(configPool, 'config_parabolic_bowl_Coriolis_parameter', config_parabolic_bowl_Coriolis_parameter)
      call mpas_pool_add_config(configPool, 'config_parabolic_bowl_eta0', config_parabolic_bowl_eta0)
      call mpas_pool_add_config(configPool, 'config_parabolic_bowl_b0', config_parabolic_bowl_b0)
      call mpas_pool_add_config(configPool, 'config_parabolic_bowl_omega', config_parabolic_bowl_omega)
      call mpas_pool_add_config(configPool, 'config_parabolic_bowl_gravity', config_parabolic_bowl_gravity)
      call mpas_pool_add_config(configPool, 'config_parabolic_bowl_adjust_domain_center', &
config_parabolic_bowl_adjust_domain_center)

   end subroutine ocn_setup_nmlrec_parabolic_bowl


   subroutine ocn_setup_nmlrec_partial_cells(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_alter_ICs_for_pcs = .false.
      character (len=StrKIND) :: config_pc_alteration_type = 'full_cell'
      real (kind=RKIND) :: config_min_pc_fraction = 0.100000

      namelist /partial_cells/ &
         config_alter_ICs_for_pcs, &
         config_pc_alteration_type, &
         config_min_pc_fraction
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, partial_cells, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_alter_ICs_for_pcs)
         call mpas_dmpar_bcast_char(dminfo, config_pc_alteration_type)
         call mpas_dmpar_bcast_real(dminfo, config_min_pc_fraction)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record partial_cells')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_alter_ICs_for_pcs = $l', logicArgs=(/config_alter_ICs_for_pcs/))
            call mpas_log_write('        config_pc_alteration_type = '//mpas_log_escape_dollars(config_pc_alteration_type))
            call mpas_log_write('        config_min_pc_fraction = $r', realArgs=(/config_min_pc_fraction/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record partial_cells.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_alter_ICs_for_pcs', config_alter_ICs_for_pcs)
      call mpas_pool_add_config(configPool, 'config_pc_alteration_type', config_pc_alteration_type)
      call mpas_pool_add_config(configPool, 'config_min_pc_fraction', config_min_pc_fraction)

   end subroutine ocn_setup_nmlrec_partial_cells


   subroutine ocn_setup_nmlrec_init_setup(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      character (len=StrKIND) :: config_init_configuration = 'none'
      logical :: config_expand_sphere = .false.
      logical :: config_realistic_coriolis_parameter = .false.
      logical :: config_write_cull_cell_mask = .true.
      character (len=StrKIND) :: config_vertical_grid = 'uniform'

      namelist /init_setup/ &
         config_init_configuration, &
         config_expand_sphere, &
         config_realistic_coriolis_parameter, &
         config_write_cull_cell_mask, &
         config_vertical_grid
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, init_setup, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_char(dminfo, config_init_configuration)
         call mpas_dmpar_bcast_logical(dminfo, config_expand_sphere)
         call mpas_dmpar_bcast_logical(dminfo, config_realistic_coriolis_parameter)
         call mpas_dmpar_bcast_logical(dminfo, config_write_cull_cell_mask)
         call mpas_dmpar_bcast_char(dminfo, config_vertical_grid)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record init_setup')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_init_configuration = '//mpas_log_escape_dollars(config_init_configuration))
            call mpas_log_write('        config_expand_sphere = $l', logicArgs=(/config_expand_sphere/))
            call mpas_log_write('        config_realistic_coriolis_parameter = $l', &
logicArgs=(/config_realistic_coriolis_parameter/))
            call mpas_log_write('        config_write_cull_cell_mask = $l', logicArgs=(/config_write_cull_cell_mask/))
            call mpas_log_write('        config_vertical_grid = '//mpas_log_escape_dollars(config_vertical_grid))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record init_setup.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_init_configuration', config_init_configuration)
      call mpas_pool_add_config(configPool, 'config_expand_sphere', config_expand_sphere)
      call mpas_pool_add_config(configPool, 'config_realistic_coriolis_parameter', config_realistic_coriolis_parameter)
      call mpas_pool_add_config(configPool, 'config_write_cull_cell_mask', config_write_cull_cell_mask)
      call mpas_pool_add_config(configPool, 'config_vertical_grid', config_vertical_grid)

   end subroutine ocn_setup_nmlrec_init_setup


   subroutine ocn_setup_nmlrec_CVTgenerator(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      real (kind=RKIND) :: config_1dCVTgenerator_stretch1 = 1.077000
      real (kind=RKIND) :: config_1dCVTgenerator_stretch2 = 1.027500
      real (kind=RKIND) :: config_1dCVTgenerator_dzSeed = 1.200000

      namelist /CVTgenerator/ &
         config_1dCVTgenerator_stretch1, &
         config_1dCVTgenerator_stretch2, &
         config_1dCVTgenerator_dzSeed
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, CVTgenerator, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_real(dminfo, config_1dCVTgenerator_stretch1)
         call mpas_dmpar_bcast_real(dminfo, config_1dCVTgenerator_stretch2)
         call mpas_dmpar_bcast_real(dminfo, config_1dCVTgenerator_dzSeed)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record CVTgenerator')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_1dCVTgenerator_stretch1 = $r', realArgs=(/config_1dCVTgenerator_stretch1/))
            call mpas_log_write('        config_1dCVTgenerator_stretch2 = $r', realArgs=(/config_1dCVTgenerator_stretch2/))
            call mpas_log_write('        config_1dCVTgenerator_dzSeed = $r', realArgs=(/config_1dCVTgenerator_dzSeed/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record CVTgenerator.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_1dCVTgenerator_stretch1', config_1dCVTgenerator_stretch1)
      call mpas_pool_add_config(configPool, 'config_1dCVTgenerator_stretch2', config_1dCVTgenerator_stretch2)
      call mpas_pool_add_config(configPool, 'config_1dCVTgenerator_dzSeed', config_1dCVTgenerator_dzSeed)

   end subroutine ocn_setup_nmlrec_CVTgenerator


   subroutine ocn_setup_nmlrec_init_vertical_grid(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      character (len=StrKIND) :: config_init_vertical_grid_type = 'z-star'

      namelist /init_vertical_grid/ &
         config_init_vertical_grid_type
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, init_vertical_grid, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_char(dminfo, config_init_vertical_grid_type)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record init_vertical_grid')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_init_vertical_grid_type = '&
//''//mpas_log_escape_dollars(config_init_vertical_grid_type))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record init_vertical_grid.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_init_vertical_grid_type', config_init_vertical_grid_type)

   end subroutine ocn_setup_nmlrec_init_vertical_grid


   subroutine ocn_setup_nmlrec_constrain_Haney_number(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      integer :: config_rx1_outer_iter_count = 20
      integer :: config_rx1_inner_iter_count = 10
      real (kind=RKIND) :: config_rx1_init_inner_weight = 0.100000
      real (kind=RKIND) :: config_rx1_max = 5.000000
      real (kind=RKIND) :: config_rx1_horiz_smooth_weight = 1.000000
      real (kind=RKIND) :: config_rx1_vert_smooth_weight = 1.000000
      real (kind=RKIND) :: config_rx1_slope_weight = 0.100000
      real (kind=RKIND) :: config_rx1_zstar_weight = 1.000000
      integer :: config_rx1_horiz_smooth_open_ocean_cells = 20
      integer :: config_rx1_min_levels = 3
      real (kind=RKIND) :: config_rx1_min_layer_thickness = 1.000000

      namelist /constrain_Haney_number/ &
         config_rx1_outer_iter_count, &
         config_rx1_inner_iter_count, &
         config_rx1_init_inner_weight, &
         config_rx1_max, &
         config_rx1_horiz_smooth_weight, &
         config_rx1_vert_smooth_weight, &
         config_rx1_slope_weight, &
         config_rx1_zstar_weight, &
         config_rx1_horiz_smooth_open_ocean_cells, &
         config_rx1_min_levels, &
         config_rx1_min_layer_thickness
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, constrain_Haney_number, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_int(dminfo, config_rx1_outer_iter_count)
         call mpas_dmpar_bcast_int(dminfo, config_rx1_inner_iter_count)
         call mpas_dmpar_bcast_real(dminfo, config_rx1_init_inner_weight)
         call mpas_dmpar_bcast_real(dminfo, config_rx1_max)
         call mpas_dmpar_bcast_real(dminfo, config_rx1_horiz_smooth_weight)
         call mpas_dmpar_bcast_real(dminfo, config_rx1_vert_smooth_weight)
         call mpas_dmpar_bcast_real(dminfo, config_rx1_slope_weight)
         call mpas_dmpar_bcast_real(dminfo, config_rx1_zstar_weight)
         call mpas_dmpar_bcast_int(dminfo, config_rx1_horiz_smooth_open_ocean_cells)
         call mpas_dmpar_bcast_int(dminfo, config_rx1_min_levels)
         call mpas_dmpar_bcast_real(dminfo, config_rx1_min_layer_thickness)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record constrain_Haney_number')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_rx1_outer_iter_count = $i', intArgs=(/config_rx1_outer_iter_count/))
            call mpas_log_write('        config_rx1_inner_iter_count = $i', intArgs=(/config_rx1_inner_iter_count/))
            call mpas_log_write('        config_rx1_init_inner_weight = $r', realArgs=(/config_rx1_init_inner_weight/))
            call mpas_log_write('        config_rx1_max = $r', realArgs=(/config_rx1_max/))
            call mpas_log_write('        config_rx1_horiz_smooth_weight = $r', realArgs=(/config_rx1_horiz_smooth_weight/))
            call mpas_log_write('        config_rx1_vert_smooth_weight = $r', realArgs=(/config_rx1_vert_smooth_weight/))
            call mpas_log_write('        config_rx1_slope_weight = $r', realArgs=(/config_rx1_slope_weight/))
            call mpas_log_write('        config_rx1_zstar_weight = $r', realArgs=(/config_rx1_zstar_weight/))
            call mpas_log_write('        config_rx1_horiz_smooth_open_ocean_cells = $i', &
intArgs=(/config_rx1_horiz_smooth_open_ocean_cells/))
            call mpas_log_write('        config_rx1_min_levels = $i', intArgs=(/config_rx1_min_levels/))
            call mpas_log_write('        config_rx1_min_layer_thickness = $r', realArgs=(/config_rx1_min_layer_thickness/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record constrain_Haney_number.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_rx1_outer_iter_count', config_rx1_outer_iter_count)
      call mpas_pool_add_config(configPool, 'config_rx1_inner_iter_count', config_rx1_inner_iter_count)
      call mpas_pool_add_config(configPool, 'config_rx1_init_inner_weight', config_rx1_init_inner_weight)
      call mpas_pool_add_config(configPool, 'config_rx1_max', config_rx1_max)
      call mpas_pool_add_config(configPool, 'config_rx1_horiz_smooth_weight', config_rx1_horiz_smooth_weight)
      call mpas_pool_add_config(configPool, 'config_rx1_vert_smooth_weight', config_rx1_vert_smooth_weight)
      call mpas_pool_add_config(configPool, 'config_rx1_slope_weight', config_rx1_slope_weight)
      call mpas_pool_add_config(configPool, 'config_rx1_zstar_weight', config_rx1_zstar_weight)
      call mpas_pool_add_config(configPool, 'config_rx1_horiz_smooth_open_ocean_cells', config_rx1_horiz_smooth_open_ocean_cells)
      call mpas_pool_add_config(configPool, 'config_rx1_min_levels', config_rx1_min_levels)
      call mpas_pool_add_config(configPool, 'config_rx1_min_layer_thickness', config_rx1_min_layer_thickness)

   end subroutine ocn_setup_nmlrec_constrain_Haney_number


