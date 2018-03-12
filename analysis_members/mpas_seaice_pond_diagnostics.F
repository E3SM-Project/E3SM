! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_pond_diagnostics
!
!> \brief MPAS sea ice analysis mode member: pond_diagnostics
!> \author Adrian K. Turner
!> \date   6th December 2016
!> \details
!>  MPAS sea ice analysis mode member: pond_diagnostics
!>
!-----------------------------------------------------------------------

module cice_pond_diagnostics

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager

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

   public :: cice_bootstrap_pond_diagnostics, &
             cice_init_pond_diagnostics, &
             cice_precompute_pond_diagnostics, &
             cice_compute_pond_diagnostics, &
             cice_restart_pond_diagnostics, &
             cice_finalize_pond_diagnostics

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine cice_bootstrap_pond_diagnostics
!
!> \brief   Bootstrap MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    6th December 2016
!> \details
!>  This routine conducts all bootstraps required for the
!>  MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_bootstrap_pond_diagnostics(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

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

   end subroutine cice_bootstrap_pond_diagnostics!}}}

!***********************************************************************
!
!  routine cice_init_pond_diagnostics
!
!> \brief   Initialize MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    6th December 2016
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_init_pond_diagnostics(domain, instance, err)!{{{

      use cice_deactivate_unneeded_fields, only: &
          cice_deactivate_unneeded_stream_fields

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

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

      call cice_deactivate_unneeded_stream_fields(domain, "pondDiagnosticsAM")

   end subroutine cice_init_pond_diagnostics!}}}

!***********************************************************************
!
!  routine cice_precompute_pond_diagnostics
!
!> \brief   Precompute MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    6th December 2016
!> \details
!>  This routine conducts all pre-computation required for this
!>  MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_precompute_pond_diagnostics(domain, instance, timeLevel, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

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

      err = 0

   end subroutine cice_precompute_pond_diagnostics!}}}

!***********************************************************************
!
!  routine cice_compute_pond_diagnostics
!
!> \brief   Compute MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    6th December 2016
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_compute_pond_diagnostics(domain, instance, timeLevel, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

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

      type(block_type), pointer :: &
           block

      type(MPAS_pool_type), pointer :: &
           tracersAggregatePool, &
           tracersPool, &
           pondDiagnosticsAMPool

      real(kind=RKIND), dimension(:), pointer :: &
           pondAreaCell, &
           pondDepthCell, &
           pondLidThicknessCell, &
           iceAreaCell, &
           levelIceAreaCell

      real(kind=RKIND), dimension(:,:,:), pointer :: &
           pondArea, &
           pondDepth

      type(field1DReal), pointer :: &
           meltPondArea, &
           meltPondAreaFinalArea, &
           meltPondDepth, &
           meltPondDepthFinalArea, &
           meltPondLidThickness, &
           meltPondLidThicknessFinalArea

      type(field2DReal), pointer :: &
           meltPondDepthCategory

      logical, pointer :: &
           config_use_cesm_meltponds, &
           config_use_level_meltponds, &
           config_use_topo_meltponds

      integer, pointer :: &
           nCellsSolve

      integer :: &
           iCell

      err = 0

      block => domain % blocklist
      do while (associated(block))

         call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracersAggregatePool)
         call MPAS_pool_get_subpool(block % structs, "tracers", tracersPool)
         call MPAS_pool_get_subpool(block % structs, "pondDiagnosticsAM", pondDiagnosticsAMPool)

         call MPAS_pool_get_array(tracersAggregatePool, "pondAreaCell", pondAreaCell)
         call MPAS_pool_get_array(tracersAggregatePool, "pondDepthCell", pondDepthCell)
         call MPAS_pool_get_array(tracersAggregatePool, "pondLidThicknessCell", pondLidThicknessCell)
         call MPAS_pool_get_array(tracersAggregatePool, "iceAreaCell", iceAreaCell)
         call MPAS_pool_get_array(tracersAggregatePool, "levelIceAreaCell", levelIceAreaCell)

         call MPAS_pool_get_array(tracersPool, "pondArea", pondArea, 1)
         call MPAS_pool_get_array(tracersPool, "pondDepth", pondDepth, 1)

         call MPAS_pool_get_field(pondDiagnosticsAMPool, "meltPondArea", meltPondArea)
         call MPAS_pool_get_field(pondDiagnosticsAMPool, "meltPondAreaFinalArea", meltPondAreaFinalArea)
         call MPAS_pool_get_field(pondDiagnosticsAMPool, "meltPondDepth", meltPondDepth)
         call MPAS_pool_get_field(pondDiagnosticsAMPool, "meltPondDepthFinalArea", meltPondDepthFinalArea)
         call MPAS_pool_get_field(pondDiagnosticsAMPool, "meltPondLidThickness", meltPondLidThickness)
         call MPAS_pool_get_field(pondDiagnosticsAMPool, "meltPondLidThicknessFinalArea", meltPondLidThicknessFinalArea)

         call MPAS_pool_get_field(pondDiagnosticsAMPool, "meltPondDepthCategory", meltPondDepthCategory)

         call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)

         call MPAS_pool_get_config(block % configs, "config_use_cesm_meltponds", config_use_cesm_meltponds)
         call MPAS_pool_get_config(block % configs, "config_use_level_meltponds", config_use_level_meltponds)
         call MPAS_pool_get_config(block % configs, "config_use_topo_meltponds", config_use_topo_meltponds)

         if (config_use_cesm_meltponds) then
            ! cesm ponds

            ! meltPondArea - apond
            if (meltPondArea % isActive) then
               do iCell = 1, nCellsSolve
                  meltPondArea % array(iCell) = pondAreaCell(iCell)
               enddo ! iCell
            endif ! is active

            ! meltPondAreaFinalArea - apond_ai
            if (meltPondAreaFinalArea % isActive) then
               do iCell = 1, nCellsSolve
                  meltPondAreaFinalArea % array(iCell) = pondAreaCell(iCell) * iceAreaCell(iCell)
               enddo ! iCell
            endif ! is active

            ! meltPondDepth - hpond
            if (meltPondDepth % isActive) then
               do iCell = 1, nCellsSolve
                  meltPondDepth % array(iCell) = pondAreaCell(iCell) * pondDepthCell(iCell)
               enddo ! iCell
            endif ! is active

            ! meltPondDepthFinalArea - hpond_ai
            if (meltPondDepthFinalArea % isActive) then
               do iCell = 1, nCellsSolve
                  meltPondDepthFinalArea % array(iCell) = pondAreaCell(iCell) * pondDepthCell(iCell) * iceAreaCell(iCell)
               enddo ! iCell
            endif ! is active

         else if (config_use_level_meltponds) then
            ! level melt ponds

            ! meltPondArea - apond
            if (meltPondArea % isActive) then
               do iCell = 1, nCellsSolve
                  meltPondArea % array(iCell) = levelIceAreaCell(iCell) * pondAreaCell(iCell)
               enddo ! iCell
            endif ! is active

            ! meltPondAreaFinalArea - apond_ai
            if (meltPondAreaFinalArea % isActive) then
               do iCell = 1, nCellsSolve
                  meltPondAreaFinalArea % array(iCell) = levelIceAreaCell(iCell) * pondAreaCell(iCell) * iceAreaCell(iCell)
               enddo ! iCell
            endif ! is active

            ! meltPondDepth - hpond
            if (meltPondDepth % isActive) then
               do iCell = 1, nCellsSolve
                  meltPondDepth % array(iCell) = levelIceAreaCell(iCell) * pondAreaCell(iCell) * pondDepthCell(iCell)
               enddo ! iCell
            endif ! is active

            ! meltPondDepthFinalArea - hpond_ai
            if (meltPondDepthFinalArea % isActive) then
               do iCell = 1, nCellsSolve
                  meltPondDepthFinalArea % array(iCell) = &
                       levelIceAreaCell(iCell) * pondAreaCell(iCell) * pondDepthCell(iCell) * iceAreaCell(iCell)
               enddo ! iCell
            endif ! is active

            ! meltPondLidThickness - ipond
            if (meltPondLidThickness % isActive) then
               do iCell = 1, nCellsSolve
                  meltPondLidThickness % array(iCell) = levelIceAreaCell(iCell) * pondAreaCell(iCell) * pondLidThicknessCell(iCell)
               enddo ! iCell
            endif ! is active

            ! meltPondLidThicknessFinalArea - ipond_ai
            if (meltPondLidThicknessFinalArea % isActive) then
               do iCell = 1, nCellsSolve
                  meltPondLidThicknessFinalArea % array(iCell) = &
                       levelIceAreaCell(iCell) * pondAreaCell(iCell) * pondLidThicknessCell(iCell) * iceAreaCell(iCell)
               enddo ! iCell
            endif ! is active

         else if (config_use_topo_meltponds) then
            ! topo ponds

            ! meltPondArea - apond
            if (meltPondArea % isActive) then
               do iCell = 1, nCellsSolve
                  meltPondArea % array(iCell) = pondAreaCell(iCell)
               enddo ! iCell
            endif ! is active

            ! meltPondAreaFinalArea - apond_ai
            if (meltPondAreaFinalArea % isActive) then
               do iCell = 1, nCellsSolve
                  meltPondAreaFinalArea % array(iCell) = pondAreaCell(iCell) * iceAreaCell(iCell)
               enddo ! iCell
            endif ! is active

            ! meltPondDepth - hpond
            if (meltPondDepth % isActive) then
               do iCell = 1, nCellsSolve
                  meltPondDepth % array(iCell) = pondAreaCell(iCell) * pondDepthCell(iCell)
               enddo ! iCell
            endif ! is active

            ! meltPondDepthFinalArea - hpond_ai
            if (meltPondDepthFinalArea % isActive) then
               do iCell = 1, nCellsSolve
                  meltPondDepthFinalArea % array(iCell) = pondAreaCell(iCell) * pondDepthCell(iCell) * iceAreaCell(iCell)
               enddo ! iCell
            endif ! is active

            ! meltPondLidThickness - ipond
            if (meltPondLidThickness % isActive) then
               do iCell = 1, nCellsSolve
                  meltPondLidThickness % array(iCell) = pondAreaCell(iCell) * pondLidThicknessCell(iCell)
               enddo ! iCell
            endif ! is active

            ! meltPondLidThicknessFinalArea - ipond_ai
            if (meltPondLidThicknessFinalArea % isActive) then
               do iCell = 1, nCellsSolve
                  meltPondLidThicknessFinalArea % array(iCell) = &
                       pondAreaCell(iCell) * pondLidThicknessCell(iCell) * iceAreaCell(iCell)
               enddo ! iCell
            endif ! is active

         endif

         ! meltPondDepthCategory - hpondn
         if (meltPondDepthCategory % isActive) then
            do iCell = 1, nCellsSolve
               meltPondDepthCategory % array(:,iCell) = pondArea(1,:,iCell) * pondDepth(1,:,iCell)
            enddo ! iCell
         endif ! is active

         block => block % next
      enddo

   end subroutine cice_compute_pond_diagnostics!}}}

!***********************************************************************
!
!  routine cice_restart_pond_diagnostics
!
!> \brief   Save restart for MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    6th December 2016
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_restart_pond_diagnostics(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

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

   end subroutine cice_restart_pond_diagnostics!}}}

!***********************************************************************
!
!  routine cice_finalize_pond_diagnostics
!
!> \brief   Finalize MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    6th December 2016
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_finalize_pond_diagnostics(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

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

   end subroutine cice_finalize_pond_diagnostics!}}}

!-----------------------------------------------------------------------

end module cice_pond_diagnostics

! vim: foldmethod=marker
