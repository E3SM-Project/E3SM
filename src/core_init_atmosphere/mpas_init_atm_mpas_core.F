! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
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
      use init_atm_test_cases
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      character(len=*), intent(out) :: startTimeStamp

      type (block_type), pointer :: block

      startTimeStamp = config_start_time

      block => domain % blocklist
      do while (associated(block))
         block % state % time_levs(1) % state % xtime % scalar = startTimeStamp
         block % mesh % sphere_radius = a
         block => block % next
      end do 

   end subroutine mpas_core_init
   
   
   subroutine mpas_core_run(domain, output_obj, output_frame)
   
      use mpas_grid_types
      use mpas_io_output
      use mpas_timer
      use init_atm_test_cases
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      type (io_output_object), intent(inout) :: output_obj
      integer, intent(inout) :: output_frame
   
      
      call init_atm_setup_test_case(domain)
   
   !
   ! Note: The following initialization calls have been moved to mpas_setup_test_case()
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
   
end module mpas_core
