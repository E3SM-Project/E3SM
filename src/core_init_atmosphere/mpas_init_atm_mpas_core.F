! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_core


   contains


   subroutine mpas_core_init(domain, startTimeStamp)
   
      use mpas_grid_types
      use mpas_configure
      use init_atm_cases
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      character(len=*), intent(out) :: startTimeStamp

      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: state, mesh 
      character (len=StrKIND), pointer :: xtime
      character (len=StrKIND), pointer :: config_start_time
      real (kind=RKIND), pointer :: sphere_radius

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

   end subroutine mpas_core_init
   
   
   subroutine mpas_core_run(domain, output_obj, output_frame)
   
      use mpas_grid_types
      use mpas_io_output
      use mpas_timer
      use init_atm_cases
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      type (io_output_object), intent(inout) :: output_obj
      integer, intent(inout) :: output_frame
   
      
      call init_atm_setup_case(domain)
   
   !
   ! Note: The following initialization calls have been moved to the mpas_init_atm_case_*() subroutines,
   !       since values computed by these routines are needed to produce initial fields
   !
   !   call atm_initialize_advection_rk(mesh)
   !   call atm_initialize_deformation_weights(mesh)
  
      call mpas_output_state_for_domain(output_obj, domain, output_frame)
   
   end subroutine mpas_core_run
   
   
   subroutine mpas_core_finalize(domain)
   
      use mpas_grid_types
   
      implicit none
   
      type (domain_type), intent(inout) :: domain 
   
   end subroutine mpas_core_finalize


   subroutine mpas_core_setup_packages(configs, packages, ierr)

      use mpas_packages
      use mpas_grid_types

      implicit none

      type (mpas_pool_type), intent(inout) :: configs
      type (mpas_pool_type), intent(inout) :: packages
      integer, intent(out) :: ierr

      ierr = 0

   end subroutine mpas_core_setup_packages
   
end module mpas_core
