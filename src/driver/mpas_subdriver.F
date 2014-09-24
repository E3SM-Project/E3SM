! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_subdriver

   use mpas_framework
   use mpas_core
   use mpas_kind_types

   ! TODO: This should be moved inside the domain type?
   !       Unless we re-arrange modules, doing so would create a circular dependence
   !       between the mpas_grid_types and mpas_stream_manager modules.
   use mpas_stream_manager, only : MPAS_streamManager_type

   type (dm_info), pointer :: dminfo
   type (domain_type), save, target :: domain
   type (domain_type), pointer :: domain_ptr
   type (io_output_object), save :: output_obj
   integer :: output_frame

   ! TODO: This should be moved inside the domain type?
   !       Unless we re-arrange modules, doing so would create a circular dependence
   !       between the mpas_grid_types and mpas_stream_manager modules.
   type (MPAS_streamManager_type), pointer :: stream_manager


   contains


   subroutine mpas_init()

      use mpas_stream_manager, only : MPAS_stream_mgr_init
      use iso_c_binding, only : c_char, c_loc, c_ptr, c_int
      use mpas_c_interfacing, only : mpas_f_to_c_string, mpas_c_to_f_string
      use mpas_timekeeping, only : MPAS_Time_type, mpas_get_clock_time, mpas_get_time
 
      implicit none

      character(len=StrKIND) :: timeStamp
      integer :: ierr
      integer, pointer :: config_frames_per_outfile

      character(kind=c_char), dimension(StrKIND+1) :: c_filename       ! StrKIND+1 for C null-termination character
      integer(kind=c_int) :: c_comm
      integer(kind=c_int) :: c_ierr
      type (c_ptr) :: mgr_p
      character(len=StrKIND) :: mesh_stream
      character(len=StrKIND) :: mesh_filename
      character(len=StrKIND) :: mesh_filename_temp
      character(kind=c_char), dimension(StrKIND+1) :: c_mesh_stream
      character(kind=c_char), dimension(StrKIND+1) :: c_mesh_filename_temp
      type (MPAS_Time_type) :: start_time
      character(len=StrKIND) :: start_timestamp


      interface
         subroutine xml_stream_parser(xmlname, mgr_p, comm, ierr) bind(c)
            use iso_c_binding, only : c_char, c_ptr, c_int
            character(kind=c_char), dimension(*), intent(in) :: xmlname
            type (c_ptr), intent(inout) :: mgr_p
            integer(kind=c_int), intent(inout) :: comm
            integer(kind=c_int), intent(out) :: ierr
         end subroutine xml_stream_parser

         subroutine xml_stream_get_filename(xmlname, streamname, comm, filename, ierr) bind(c)
            use iso_c_binding, only : c_char, c_int
            character(kind=c_char), dimension(*), intent(in) :: xmlname
            character(kind=c_char), dimension(*), intent(in) :: streamname
            integer(kind=c_int), intent(inout) :: comm
            character(kind=c_char), dimension(*), intent(in) :: filename
            integer(kind=c_int), intent(out) :: ierr
         end subroutine xml_stream_get_filename
      end interface


      domain_ptr => domain

      !
      ! Initialize infrastructure
      !
      call mpas_framework_init(dminfo, domain_ptr)
      call mpas_core_setup_packages(domain_ptr % configs, domain_ptr % packages, ierr)
      call mpas_core_setup_clock(domain_ptr % clock, domain_ptr % configs, ierr)

      call mpas_timer_start("total time")
      call mpas_timer_start("initialize")

      !
      ! Using information from the namelist, a graph.info file, and a file containing
      !    mesh fields, build halos and allocate blocks in the domain
      !
      call mpas_core_get_mesh_stream(domain_ptr % configs, mesh_stream, ierr)
      if (ierr /= 0) then
         call mpas_dmpar_abort(dminfo)
      end if
write(0,*) 'MGD DEBUG mesh stream is '//trim(mesh_stream)
      call mpas_f_to_c_string(domain % streams_filename, c_filename)
      call mpas_f_to_c_string(mesh_stream, c_mesh_stream)
      c_comm = dminfo % comm
      call xml_stream_get_filename(c_filename, c_mesh_stream, c_comm, c_mesh_filename_temp, c_ierr)
      if (c_ierr /= 0) then
         call mpas_dmpar_abort(dminfo)
      end if
      call mpas_c_to_f_string(c_mesh_filename_temp, mesh_filename_temp)
write(0,*) 'MGD DEBUG mesh filename template is '//trim(mesh_filename_temp)
      start_time = mpas_get_clock_time(domain_ptr % clock, MPAS_START_TIME, ierr)
      call mpas_get_time(start_time, dateTimeString=start_timestamp, ierr=ierr)
      call mpas_expand_string(start_timestamp, mesh_filename_temp, mesh_filename)
write(0,*) 'MGD DEBUG mesh filename is '//trim(mesh_filename)
      call mpas_bootstrap_framework(domain_ptr, mesh_filename)

      !
      ! Read initial fields into blocks
      !
      call mpas_input_state_for_domain(domain_ptr)


      !
      ! Set up run-time streams
      !
      call MPAS_stream_mgr_init(stream_manager, domain % clock, domain % blocklist % allFields, domain % packages)

      call mpas_setup_immutable_streams(stream_manager)

      mgr_p = c_loc(stream_manager)
      call xml_stream_parser(c_filename, mgr_p, c_comm, c_ierr)
      if (c_ierr /= 0) then
         call mpas_dmpar_abort(dminfo)
      end if


      !
      ! Initialize core
      !
      call mpas_core_init(domain_ptr, stream_manager, timeStamp)

      call mpas_timer_stop("initialize")


      !
      ! Set up output streams to be written to by the MPAS core
      !
      output_frame = 1

      call mpas_pool_get_config(domain_ptr % configs, 'config_frames_per_outfile', config_frames_per_outfile)

      if(config_frames_per_outfile > 0) then
         call mpas_output_state_init(output_obj, domain_ptr, "OUTPUT", trim(timeStamp))
      else
         call mpas_output_state_init(output_obj, domain_ptr, "OUTPUT")         
      end if

   end subroutine mpas_init


   subroutine mpas_run()

      implicit none

      call mpas_core_run(domain_ptr, stream_manager, output_obj, output_frame)

   end subroutine mpas_run


   subroutine mpas_finalize()

      use mpas_stream_manager, only : MPAS_stream_mgr_finalize
   
      implicit none

      !
      ! Finalize output streams
      !
      call mpas_output_state_finalize(output_obj, domain_ptr % dminfo)


      !
      ! Finalize core
      !
      call mpas_core_finalize(domain_ptr, stream_manager)

      call mpas_timer_stop("total time")
      call mpas_timer_write()


      !
      ! Finalize infrastructure
      !
      call MPAS_stream_mgr_finalize(stream_manager)

      call mpas_framework_finalize(dminfo, domain_ptr)

   end subroutine mpas_finalize

end module mpas_subdriver
