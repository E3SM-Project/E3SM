! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  oac_TEMPLATE
!
!> \brief MPAS ocean analysis core member: TEMPLATE
!> \author FILL_IN_AUTHOR
!> \date   FILL_IN_DATE
!> \details
!>  MPAS ocean analysis core member: TEMPLATE
!>  In order to add a new analysis member, do the following:
!>  1. Copy these to your new analysis member name:
!>     mpas_oac_TEMPLATE.F
!>     Registry_oac_TEMPLATE.xml
!>  
!>  2. In those two new files, replace the following text:
!>     TEMPLATE, amTemplate, FILL_IN_AUTHOR, FILL_IN_DATE
!>  
!>  3. Add a line for your registry to
!>     Registry_shared_oac.xml
!>  
!>  4. In mpas_oac_driver.F, add calls for your analysis member
!>     by copying lines with TEMPLATE.
!>  
!>  5. In src/core_ocean_analysis/shared_oac/Makefile, add your
!>     new analysis member everywhere you see
!>     mpas_oac_global_stats.o
!>  
!
!-----------------------------------------------------------------------

module oac_TEMPLATE

   use mpas_configure
   use mpas_grid_types
   use mpas_timer
   use mpas_dmpar

   use ocn_diagnostics_routines

   implicit none
   private
   save

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

   public :: oac_init_TEMPLATE, &
             oac_compute_TEMPLATE, &
             oac_write_TEMPLATE, &
             oac_restart_TEMPLATE, &
             oac_finalize_TEMPLATE

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   type (timer_node), pointer :: amTemplateTimer

!***********************************************************************

contains

!***********************************************************************
!
!  routine oac_init_TEMPLATE
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  FILL_IN_AUTHOR
!> \date    FILL_IN_DATE
!> \details 
!>  This routine conducts all initializations required for the 
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine oac_init_TEMPLATE(err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

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

      integer :: err_tmp

      err = 0

   end subroutine oac_init_TEMPLATE!}}}

!***********************************************************************
!
!  routine oac_compute_TEMPLATE
!
!> \brief   Computeialize MPAS-Ocean analysis member
!> \author  FILL_IN_AUTHOR
!> \date    FILL_IN_DATE
!> \details 
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine oac_compute_TEMPLATE(amTemplate, state_ptr, dminfo, timeLevel, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (state_type), intent(in), pointer :: state_ptr
      type (dm_info), intent(in), pointer :: dminfo
      integer, intent(in) :: timeLevel

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      type (amTemplate_type), intent(inout) :: amTemplate
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: err_tmp
      type (block_type), pointer :: block

      type (state_type), pointer :: state
      type (mesh_type), pointer :: mesh
      type (scratch_type), pointer :: scratch

      err = 0

      call mpas_timer_start("compute_TEMPLATE", .false., amTemplateTimer)

      block => state_ptr % block 
      do while (associated(block))
         state => block % state % time_levs(timeLevel) % state
         mesh => block % mesh
         scratch => block % scratch

         ! block calculations may be placed here

         block => block % next
      end do

      ! mpi gather/scatter calls may be placed here

      call mpas_timer_stop("TEMPLATE", amTemplateTimer)

   end subroutine oac_compute_TEMPLATE!}}}

!***********************************************************************
!
!  routine oac_restart_TEMPLATE
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  FILL_IN_AUTHOR
!> \date    FILL_IN_DATE
!> \details 
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine oac_restart_TEMPLATE(err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

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

      integer :: err_tmp

      err = 0

   end subroutine oac_restart_TEMPLATE!}}}

!***********************************************************************
!
!  routine oac_write_TEMPLATE
!
!> \brief   Driver for MPAS-Ocean analysis output
!> \author  FILL_IN_AUTHOR
!> \date    FILL_IN_DATE
!> \details 
!>  This routine writes all output for this MPAS-Ocean analysis member.
!>  At this time this is just a stub, and all analysis output is written
!>  to the output file specified by config_output_name.
!
!-----------------------------------------------------------------------

   subroutine oac_write_TEMPLATE(err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

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

      integer :: err_tmp

      err = 0

   end subroutine oac_write_TEMPLATE!}}}

!***********************************************************************
!
!  routine oac_finalize_TEMPLATE
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  FILL_IN_AUTHOR
!> \date    FILL_IN_DATE
!> \details 
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine oac_finalize_TEMPLATE(err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

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

      integer :: err_tmp

      err = 0

print *, 'oac_finalize_TEMPLATE'

   end subroutine oac_finalize_TEMPLATE!}}}

end module oac_TEMPLATE

! vim: foldmethod=marker
