! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_miscellaneous
!
!> \brief MPAS sea ice analysis mode member: miscellaneous
!> \author Adrian K. Turner
!> \date   22nd December 2015
!> \details
!>  MPAS sea ice analysis mode member: miscellaneous
!>
!-----------------------------------------------------------------------

module cice_miscellaneous

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

   public :: cice_bootstrap_miscellaneous, &
             cice_init_miscellaneous, &
             cice_precompute_miscellaneous, &
             cice_compute_miscellaneous, &
             cice_restart_miscellaneous, &
             cice_finalize_miscellaneous

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine cice_bootstrap_miscellaneous
!
!> \brief   Bootstrap MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    22nd December 2015
!> \details
!>  This routine conducts all bootstraps required for the
!>  MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_bootstrap_miscellaneous(domain, instance, err)!{{{

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

   end subroutine cice_bootstrap_miscellaneous!}}}

!***********************************************************************
!
!  routine cice_init_miscellaneous
!
!> \brief   Initialize MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    22nd December 2015
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_init_miscellaneous(domain, instance, err)!{{{

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

   end subroutine cice_init_miscellaneous!}}}

!***********************************************************************
!
!  routine cice_precompute_miscellaneous
!
!> \brief   Precompute MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    22nd December 2015
!> \details
!>  This routine conducts all pre-computation required for this
!>  MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_precompute_miscellaneous(domain, instance, timeLevel, err)!{{{

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

   end subroutine cice_precompute_miscellaneous!}}}

!***********************************************************************
!
!  routine cice_compute_miscellaneous
!
!> \brief   Compute MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    22nd December 2015
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_compute_miscellaneous(domain, instance, timeLevel, err)!{{{

     use ice_constants_colpkg, only: &
          awtvdr, &
          awtidr, &
          awtvdf, &
          awtidf

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
           miscellaneousAMPool, &
           tracersAggregatePool, &
           shortwavePool, &
           icestatePool, &
           velocitySolverPool

      real(kind=RKIND), dimension(:,:), pointer :: &
           iceSalinityCell

      real(kind=RKIND), dimension(:), pointer :: &
           bulkSalinity, &
           broadbandAlbedo, &
           iceAreaCellInitial, &
           albedoVisibleDirectCell, &
           albedoIRDirectCell, &
           albedoVisibleDiffuseCell, &
           albedoIRDiffuseCell, &
           coriolisStressVertexU, &
           coriolisStressVertexV, &
           uVelocity, &
           vVelocity, &
           totalMassVertexfVertex

      integer, pointer :: &
           nCellsSolve, &
           nVerticesSolve, &
           nIceLayers

      integer :: &
           iCell, &
           iVertex

      logical, pointer :: &
           config_use_column_package, &
           config_use_column_shortwave, &
           config_use_velocity_solver

      err = 0

      call MPAS_pool_get_config(domain % configs, "config_use_column_package", config_use_column_package)
      call MPAS_pool_get_config(domain % configs, "config_use_column_shortwave", config_use_column_shortwave)
      call MPAS_pool_get_config(domain % configs, "config_use_velocity_solver", config_use_velocity_solver)

      ! grid cell mean salinity
      if (config_use_column_package) then
         block => domain % blocklist
         do while(associated(block))

            call MPAS_pool_get_subpool(block % structs, "miscellaneousAM", miscellaneousAMPool)
            call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracersAggregatePool)

            call MPAS_pool_get_array(tracersAggregatePool, "iceSalinityCell", iceSalinityCell)
            call MPAS_pool_get_array(miscellaneousAMPool, "bulkSalinity", bulkSalinity)

            call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)
            call MPAS_pool_get_dimension(block % dimensions, "nIceLayers", nIceLayers)

            do iCell = 1, nCellsSolve

               bulkSalinity(iCell) = sum(iceSalinityCell(:,iCell)) / real(nIceLayers,RKIND)

            enddo ! iCell

            block => block % next
         enddo
      endif

      ! snow/ice broad band albedo
      if (config_use_column_package .and. config_use_column_shortwave) then
         block => domain % blocklist
         do while(associated(block))

            call MPAS_pool_get_subpool(block % structs, "miscellaneousAM", miscellaneousAMPool)
            call MPAS_pool_get_subpool(block % structs, "shortwave", shortwavePool)
            call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)

            call MPAS_pool_get_array(miscellaneousAMPool, "broadbandAlbedo", broadbandAlbedo)
            call MPAS_pool_get_array(shortwavePool, "albedoVisibleDirectCell", albedoVisibleDirectCell)
            call MPAS_pool_get_array(shortwavePool, "albedoIRDirectCell", albedoIRDirectCell)
            call MPAS_pool_get_array(shortwavePool, "albedoVisibleDiffuseCell", albedoVisibleDiffuseCell)
            call MPAS_pool_get_array(shortwavePool, "albedoIRDiffuseCell", albedoIRDiffuseCell)
            call MPAS_pool_get_array(icestatePool, "iceAreaCellInitial", iceAreaCellInitial)

            call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)

            do iCell = 1, nCellsSolve

               broadbandAlbedo(iCell) = &
                    (awtvdr * albedoVisibleDirectCell(iCell) + &
                    awtidr * albedoIRDirectCell(iCell) + &
                    awtvdf * albedoVisibleDiffuseCell(iCell) + &
                    awtidf * albedoIRDiffuseCell(iCell)) * &
                    iceAreaCellInitial(iCell)

            enddo ! iCell

            block => block % next
         enddo
      endif

      ! coriolis stress
      if (config_use_velocity_solver) then
         block => domain % blocklist
         do while(associated(block))

            call MPAS_pool_get_subpool(block % structs, "miscellaneousAM", miscellaneousAMPool)
            call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)

            call MPAS_pool_get_array(miscellaneousAMPool, "coriolisStressVertexU", coriolisStressVertexU)
            call MPAS_pool_get_array(miscellaneousAMPool, "coriolisStressVertexV", coriolisStressVertexV)
            call MPAS_pool_get_array(velocitySolverPool, "uVelocity", uVelocity)
            call MPAS_pool_get_array(velocitySolverPool, "vVelocity", vVelocity)
            call MPAS_pool_get_array(velocitySolverPool, "totalMassVertexfVertex", totalMassVertexfVertex)

            call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)

            do iVertex = 1, nVerticesSolve

               coriolisStressVertexU(iVertex) =  vVelocity(iVertex) * totalMassVertexfVertex(iVertex)
               coriolisStressVertexV(iVertex) = -uVelocity(iVertex) * totalMassVertexfVertex(iVertex)

            enddo ! iVertex

            block => block % next
         enddo
      endif

   end subroutine cice_compute_miscellaneous!}}}

!***********************************************************************
!
!  routine cice_restart_miscellaneous
!
!> \brief   Save restart for MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    22nd December 2015
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_restart_miscellaneous(domain, instance, err)!{{{

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

   end subroutine cice_restart_miscellaneous!}}}

!***********************************************************************
!
!  routine cice_finalize_miscellaneous
!
!> \brief   Finalize MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    22nd December 2015
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_finalize_miscellaneous(domain, instance, err)!{{{

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

   end subroutine cice_finalize_miscellaneous!}}}

!-----------------------------------------------------------------------

end module cice_miscellaneous

! vim: foldmethod=marker
