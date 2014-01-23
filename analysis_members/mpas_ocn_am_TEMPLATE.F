! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_am_TEMPLATE
!
!> \brief MPAS ocean analysis core member: TEMPLATE
!> \author FILL_IN_AUTHOR
!> \date   FILL_IN_DATE
!> \details
!>  MPAS ocean analysis core member: TEMPLATE
!>  In order to add a new analysis member, do the following:
!>  1. Copy these to your new analysis member name:
!>     mpas_ocn_am_TEMPLATE.F
!>     Registry_ocn_am_TEMPLATE.xml
!>  
!>  2. In those two new files, replace the following text:
!>     TEMPLATE, amTemplate, FILL_IN_AUTHOR, FILL_IN_DATE
!>     Note TEMPLATE uses underscores, like global_stats, while
!>     amTemplate uses caps, e.g. amGlobalStats.
!>  
!>  3. Add a #include line for your registry to
!>     Registry_shared_ocn_am.xml
!>  
!>  4. In mpas_ocn_am_driver.F, add calls for your analysis member
!>     by copying lines with TEMPLATE.
!>  
!>  5. In src/core_ocean_analysis/shared_ocn_am/Makefile, add your
!>     new analysis member everywhere you see
!>     mpas_ocn_am_global_stats.o
!>  
!
!-----------------------------------------------------------------------

module ocn_am_TEMPLATE

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

   public :: ocn_am_setup_packages_TEMPLATE, &
             ocn_am_init_TEMPLATE, &
             ocn_am_compute_TEMPLATE, &
             ocn_am_write_TEMPLATE, &
             ocn_am_restart_TEMPLATE, &
             ocn_am_finalize_TEMPLATE

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
!  routine ocn_am_setup_packages_TEMPLATE
!
!> \brief   Set up packages for MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    November 2013
!> \details 
!>  This routine is intended to configure the packages for this MPAS
!>   ocean analysis member
!
!-----------------------------------------------------------------------

   subroutine ocn_am_setup_packages_TEMPLATE(err)!{{{

      use mpas_packages

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

      err = 0

      ! turn on package for this analysis member
      amTemplateActive = .true.

   end subroutine ocn_am_setup_packages_TEMPLATE!}}}


!***********************************************************************
!
!  routine ocn_am_init_TEMPLATE
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  FILL_IN_AUTHOR
!> \date    FILL_IN_DATE
!> \details 
!>  This routine conducts all initializations required for the 
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_am_init_TEMPLATE(domain, err)!{{{

      use mpas_packages

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

      ! turn on package for this analysis member
      TEMPLATEActive = .true.

   end subroutine ocn_am_init_TEMPLATE!}}}

!***********************************************************************
!
!  routine ocn_am_compute_TEMPLATE
!
!> \brief   Computeialize MPAS-Ocean analysis member
!> \author  FILL_IN_AUTHOR
!> \date    FILL_IN_DATE
!> \details 
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_am_compute_TEMPLATE(domain, timeLevel, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

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

      type (amTemplate_type), pointer :: amTemplate
      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (state_type), pointer :: state
      type (mesh_type), pointer :: mesh
      type (scratch_type), pointer :: scratch
      type (diagnostics_type), pointer :: diagnostics

      err = 0

      dminfo = domain % dminfo

      call mpas_timer_start("compute_TEMPLATE", .false., amTemplateTimer)

      block => domain % block 
      do while (associated(block))
         state => block % state % time_levs(timeLevel) % state
         mesh => block % mesh
         scratch => block % scratch
         diagnostics => block % diagnostics
         amTemplate => block % amTemplate

         ! Computations which are functions of nCells, nEdges, or nVertices 
         ! must be placed within this block loop

         block => block % next
      end do

      ! mpi gather/scatter calls may be placed here.
      ! Here are some examples.  See mpas_oac_global_stats.F for further details.
      call mpas_dmpar_sum_real_array(dminfo, nVariables, sumSquares(1:nVariables), reductions(1:nVariables))
      call mpas_dmpar_min_real_array(dminfo, nMins, mins(1:nMins), reductions(1:nMins))
      call mpas_dmpar_max_real_array(dminfo, nMaxes, maxes(1:nMaxes), reductions(1:nMaxes))

      ! Computations which are not functions of nCells, nEdges, or nVertices 
      ! may be placed outside of block loops
      amTemplate => domain % blocklist % amTemplate

      call mpas_timer_stop("TEMPLATE", amTemplateTimer)

   end subroutine ocn_am_compute_TEMPLATE!}}}

!***********************************************************************
!
!  routine ocn_am_restart_TEMPLATE
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  FILL_IN_AUTHOR
!> \date    FILL_IN_DATE
!> \details 
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_am_restart_TEMPLATE(domain, err)!{{{

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

   end subroutine ocn_am_restart_TEMPLATE!}}}

!***********************************************************************
!
!  routine ocn_am_write_TEMPLATE
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

   subroutine ocn_am_write_TEMPLATE(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(in) :: domain

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

      err = 0

   end subroutine ocn_am_write_TEMPLATE!}}}

!***********************************************************************
!
!  routine ocn_am_finalize_TEMPLATE
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  FILL_IN_AUTHOR
!> \date    FILL_IN_DATE
!> \details 
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_am_finalize_TEMPLATE(domain, err)!{{{

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

   end subroutine ocn_am_finalize_TEMPLATE!}}}

end module ocn_am_TEMPLATE

! vim: foldmethod=marker
