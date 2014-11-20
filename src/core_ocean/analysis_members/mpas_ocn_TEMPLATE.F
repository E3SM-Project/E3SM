! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_TEMPLATE
!
!> \brief MPAS ocean analysis core member: TEMPLATE
!> \author FILL_IN_AUTHOR
!> \date   FILL_IN_DATE
!> \details
!>  MPAS ocean analysis core member: TEMPLATE
!>  In order to add a new analysis member, do the following:
!>  1. Copy these to your new analysis member name:
!>     cp mpas_ocn_TEMPLATE.F mpas_ocn_your_new_name.F
!>     cp Registry_ocn_TEMPLATE.xml Registry_ocn_your_new_name.xml
!>  
!>  2. In those two new files, replace the following text:
!>     TEMPLATE, amTemplate, FILL_IN_AUTHOR, FILL_IN_DATE
!>     Note TEMPLATE uses underscores, like global_stats, while
!>     amTemplate uses caps, e.g. amGlobalStats.
!>  
!>  3. Add a #include line for your registry to
!>     Registry_analysis_members.xml
!>  
!>  4. In mpas_ocn_analysis_driver.F, add calls for your analysis member
!>     by copying lines with TEMPLATE.
!>  
!>  5. In src/core_ocean/analysis_members/Makefile, add your
!>     new analysis member everywhere you see
!>     mpas_ocn_global_stats.o
!>  
!
!-----------------------------------------------------------------------

module ocn_TEMPLATE

   use mpas_grid_types
   use mpas_timer
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager

   use ocn_constants
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

   public :: ocn_setup_packages_TEMPLATE, &
             ocn_init_TEMPLATE, &
             ocn_init_alarms_TEMPLATE, &
             ocn_compute_TEMPLATE, &
             ocn_restart_TEMPLATE, &
             ocn_finalize_TEMPLATE

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
!  routine ocn_setup_packages_TEMPLATE
!
!> \brief   Set up packages for MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    November 2013
!> \details 
!>  This routine is intended to configure the packages for this MPAS
!>   ocean analysis member
!
!-----------------------------------------------------------------------

   subroutine ocn_setup_packages_TEMPLATE(configPool, packagePool, err)!{{{

      use mpas_packages

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: configPool
      type (mpas_pool_type), intent(in) :: packagePool

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
      logical, pointer :: amTemplateActive

      err = 0

      call mpas_pool_get_package(packagePool, 'amTemplateActive', amTemplateActive)

      ! turn on package for this analysis member
      amTemplateActive = .true.

   end subroutine ocn_setup_packages_TEMPLATE!}}}


!***********************************************************************
!
!  routine ocn_init_TEMPLATE
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  FILL_IN_AUTHOR
!> \date    FILL_IN_DATE
!> \details 
!>  This routine conducts all initializations required for the 
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_TEMPLATE(domain, err)!{{{

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

   end subroutine ocn_init_TEMPLATE!}}}

!***********************************************************************
!
!  routine ocn_init_alarms_TEMPLATE
!
!> \brief   Initialize alarms for MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    November 2013
!> \details 
!>  This routine conducts all alarm initializations required for the 
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_alarms_TEMPLATE(clock, startTime, computeAlarmID, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      integer, intent(in) :: computeAlarmID
      type (MPAS_Time_Type), intent(in) :: startTime

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (MPAS_Clock_type), intent(inout) :: clock

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

      integer :: ierr
      character(len=StrKIND) :: compute_interval
      type (MPAS_Time_Type) :: alarmStartTime
      type (MPAS_TimeInterval_type) :: alarmTimeStep

      character (len=StrKIND), pointer :: config_TEMPLATE_compute_interval, config_output_interval

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_TEMPLATE_compute_interval', config_TEMPLATE_compute_interval)
      call mpas_pool_get_config(ocnConfigs, 'config_output_interval', config_output_interval)

      ! set compute alarm for this analysis member
      if (config_TEMPLATE_compute_interval=='same_as_output') then
         compute_interval = config_output_interval
      else
         compute_interval = config_TEMPLATE_compute_interval
      endif

      call mpas_set_timeInterval(alarmTimeStep, timeString=compute_interval, ierr=ierr)
      alarmStartTime = startTime + alarmTimeStep
      call mpas_add_clock_alarm(clock, computeAlarmID, alarmStartTime, alarmTimeStep, ierr=ierr)

   end subroutine ocn_init_alarms_TEMPLATE!}}}

!***********************************************************************
!
!  routine ocn_compute_TEMPLATE
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  FILL_IN_AUTHOR
!> \date    FILL_IN_DATE
!> \details 
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_TEMPLATE(domain, timeLevel, err)!{{{

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

      type (mpas_pool_type), pointer :: amTemplatePool
      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: amTemplate

      ! Here are some example variables which may be needed for your analysis member
      integer, pointer :: nVertLevels, nCellsSolve, nEdgesSolve, nVerticesSolve, num_tracers
      integer :: iTracer, k, iCell
      integer, dimension(:), pointer :: maxLevelCell, maxLevelEdgeTop, maxLevelVertexBot

      real (kind=RKIND), dimension(:), pointer ::  areaCell, dcEdge, dvEdge

      err = 0

      dminfo = domain % dminfo

      call mpas_timer_start("compute_TEMPLATE", .false., amTemplateTimer)

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'amTemplate', amTemplatePool)

         ! Here are some example variables which may be needed for your analysis member
         call mpas_pool_get_dimension(statePool, 'num_tracers', num_tracers)

         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(block % dimensions, 'nEdgesSolve', nEdgesSolve)
         call mpas_pool_get_dimension(block % dimensions, 'nVerticesSolve', nVerticesSolve)

         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
         call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
         call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
         call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
         call mpas_pool_get_array(meshPool, 'maxLevelVertexBot', maxLevelVertexBot)

         ! Computations which are functions of nCells, nEdges, or nVertices 
         ! must be placed within this block loop
         ! Here are some example loops
         do iCell = 1,nCellsSolve
            do k = 1, maxLevelCell(iCell)
               do iTracer = 1, num_tracers
               ! computations on tracers(iTracer,k, iCell)
               end do
            end do
         end do

         block => block % next
      end do

      ! mpi gather/scatter calls may be placed here.
      ! Here are some examples.  See mpas_oac_global_stats.F for further details.
!      call mpas_dmpar_sum_real_array(dminfo, nVariables, sumSquares(1:nVariables), reductions(1:nVariables))
!      call mpas_dmpar_min_real_array(dminfo, nMins, mins(1:nMins), reductions(1:nMins))
!      call mpas_dmpar_max_real_array(dminfo, nMaxes, maxes(1:nMaxes), reductions(1:nMaxes))

      ! Even though some variables do not include an index that is decomposed amongst 
      ! domain partitions, we assign them within a block loop so that all blocks have the
      ! correct values for writing output.
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'amTemplate', amTemplatePool)

         ! assignment of final amTemplate variables could occur here.

         block => block % next
      end do

      call mpas_timer_stop("TEMPLATE", amTemplateTimer)

   end subroutine ocn_compute_TEMPLATE!}}}

!***********************************************************************
!
!  routine ocn_restart_TEMPLATE
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  FILL_IN_AUTHOR
!> \date    FILL_IN_DATE
!> \details 
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_TEMPLATE(domain, err)!{{{

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

   end subroutine ocn_restart_TEMPLATE!}}}

!***********************************************************************
!
!  routine ocn_finalize_TEMPLATE
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  FILL_IN_AUTHOR
!> \date    FILL_IN_DATE
!> \details 
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_TEMPLATE(domain, err)!{{{

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

   end subroutine ocn_finalize_TEMPLATE!}}}

end module ocn_TEMPLATE

! vim: foldmethod=marker
