! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module init_atm_core


   contains


   function init_atm_core_init(domain, startTimeStamp) result(ierr)
   
      use mpas_derived_types
      use mpas_stream_manager
      use mpas_io_streams, only : MPAS_STREAM_NEAREST
      use init_atm_cases
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      character(len=*), intent(out) :: startTimeStamp

      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: state, mesh 
      character (len=StrKIND), pointer :: xtime
      character (len=StrKIND), pointer :: config_start_time
      real (kind=RKIND), pointer :: sphere_radius
      integer :: ierr


      ierr = 0

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', state)
         call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
         call mpas_pool_get_array(state, 'xtime', xtime)
         call mpas_pool_get_config(mesh, 'sphere_radius', sphere_radius)
         call mpas_pool_get_config(block % configs, 'config_start_time', config_start_time)

         startTimeStamp = config_start_time
         xtime = config_start_time
         domain % sphere_radius = a      ! Appears in output files
         sphere_radius = a               ! Used in setting up test cases

         block => block % next
      end do 

      call MPAS_stream_mgr_add_att(domain % streamManager, 'sphere_radius', domain % sphere_radius, streamID='output', ierr=ierr)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'sphere_radius', domain % sphere_radius, streamID='surface', ierr=ierr)

      !
      ! We don't actually expect the time in the (most likely 'static') file to
      !    match the time in the namelist, so just read whatever time we find in
      !    the input file.
      !
      call MPAS_stream_mgr_read(domain % streamManager, whence=MPAS_STREAM_NEAREST, ierr=ierr)
      call MPAS_stream_mgr_reset_alarms(domain % streamManager, direction=MPAS_STREAM_INPUT, ierr=ierr)

   end function init_atm_core_init
   
   
   function init_atm_core_run(domain) result(ierr)
   
      use mpas_derived_types
      use mpas_stream_manager
      use mpas_timer
      use init_atm_cases
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      integer :: ierr


      ierr = 0
      
      call init_atm_setup_case(domain, domain % streamManager)
   
   !
   ! Note: The following initialization calls have been moved to the mpas_init_atm_case_*() subroutines,
   !       since values computed by these routines are needed to produce initial fields
   !
   !   call atm_initialize_advection_rk(mesh)
   !   call atm_initialize_deformation_weights(mesh)
  
      call mpas_stream_mgr_write(domain % streamManager, ierr=ierr)
      call mpas_stream_mgr_reset_alarms(domain % streamManager, direction=MPAS_STREAM_OUTPUT, ierr=ierr)
   
   end function init_atm_core_run
   
   
   function init_atm_core_finalize(domain) result(ierr)
   
      use mpas_derived_types
      use mpas_decomp
      use mpas_stream_manager
      use mpas_io_units, only : stderrUnit
   
      implicit none
   
      type (domain_type), intent(inout) :: domain 
      integer :: ierr


      ierr = 0

      call mpas_decomp_destroy_decomp_list(domain % decompositions)

      write(stderrUnit,'(a)') ''
      write(stderrUnit,'(a)') '********************************************************'
      write(stderrUnit,'(a)') '   Finished running the init_atmosphere core'
      write(stderrUnit,'(a)') '********************************************************'
   
   end function init_atm_core_finalize
   
end module init_atm_core
