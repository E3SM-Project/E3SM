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
      use mpas_c_interfacing, only : mpas_f_to_c_string
 
      implicit none

      character(len=StrKIND) :: timeStamp
      integer :: ierr
      integer, pointer :: config_frames_per_outfile

      character(kind=c_char), dimension(StrKIND+1) :: c_filename       ! StrKIND+1 for C null-termination character
      integer(kind=c_int) :: c_comm
      integer(kind=c_int) :: c_ierr
      type (c_ptr) :: mgr_p


      interface
         subroutine xml_stream_parser(fname, mgr_p, comm, ierr) bind(c)
            use iso_c_binding, only : c_char, c_ptr, c_int
            character(kind=c_char), dimension(*), intent(in) :: fname
            type (c_ptr), intent(inout) :: mgr_p
            integer(kind=c_int), intent(inout) :: comm
            integer(kind=c_int), intent(out) :: ierr
         end subroutine xml_stream_parser
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

      call mpas_input_state_for_domain(domain_ptr)


      !
      ! Set up run-time streams
      !
      call MPAS_stream_mgr_init(stream_manager, domain % clock, domain % blocklist % allFields)

      call mpas_f_to_c_string(domain % streams_filename, c_filename)
      mgr_p = c_loc(stream_manager)
      c_comm = dminfo % comm
      call xml_stream_parser(c_filename, mgr_p, c_comm, c_ierr)
      if (c_ierr /= 0) then
         call mpas_dmpar_abort(dminfo)
      end if


      !
      ! Initialize core
      !
      call mpas_core_init(domain_ptr, timeStamp)

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

      call mpas_core_run(domain_ptr, output_obj, output_frame)

   end subroutine mpas_run


   subroutine mpas_finalize()

      use mpas_stream_manager, only : MPAS_stream_mgr_finalize
   
      implicit none

      !
      ! Finalize output streams
      !
      call mpas_output_state_finalize(output_obj, domain_ptr % dminfo)

      call MPAS_stream_mgr_finalize(stream_manager)


      !
      ! Finalize core
      !
      call mpas_core_finalize(domain_ptr)

      call mpas_timer_stop("total time")
      call mpas_timer_write()


      !
      ! Finalize infrastructure
      !
      call mpas_framework_finalize(dminfo, domain_ptr)

   end subroutine mpas_finalize

end module mpas_subdriver
