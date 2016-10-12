! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_area_variables
!
!> \brief MPAS sea ice analysis mode member: area_variables
!> \author Adrian K. Turner
!> \date   5th January 2016
!> \details
!>  MPAS sea ice analysis mode member: area_variables
!>
!-----------------------------------------------------------------------

module cice_area_variables

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

   public :: cice_bootstrap_area_variables, &
             cice_init_area_variables, &
             cice_precompute_area_variables, &
             cice_compute_area_variables, &
             cice_restart_area_variables, &
             cice_finalize_area_variables

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine cice_bootstrap_area_variables
!
!> \brief   Bootstrap MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    5th January 2016
!> \details
!>  This routine conducts all bootstraps required for the
!>  MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_bootstrap_area_variables(domain, instance, err)!{{{

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

   end subroutine cice_bootstrap_area_variables!}}}

!***********************************************************************
!
!  routine cice_init_area_variables
!
!> \brief   Initialize MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    5th January 2016
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_init_area_variables(domain, instance, err)!{{{

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

      call cice_deactivate_unneeded_stream_fields(domain, "areaVariablesAM")

   end subroutine cice_init_area_variables!}}}

!***********************************************************************
!
!  routine cice_precompute_area_variables
!
!> \brief   Precompute MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    5th January 2016
!> \details
!>  This routine conducts all pre-computation required for this
!>  MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_precompute_area_variables(domain, instance, timeLevel, err)!{{{

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

   end subroutine cice_precompute_area_variables!}}}

!***********************************************************************
!
!  routine cice_compute_area_variables
!
!> \brief   Compute MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    5th January 2016
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_compute_area_variables(domain, instance, timeLevel, err)!{{{

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
           areaVariablesAMPool, &
           icestatePool, &
           atmosCouplingPool, &
           shortwavePool, &
           atmosFluxesPool

      type(field1DReal), pointer :: &
           snowfallRateInitialArea, &
           rainfallRateInitialArea, &
           absorbedShortwaveFluxInitialArea, &
           latentHeatFluxInitialArea, &
           sensibleHeatFluxInitialArea, &
           longwaveUpInitialArea, &
           evaporativeWaterFluxInitialArea, &
           surfaceHeatFluxInitialArea, &
           surfaceConductiveFluxInitialArea

      type(field2DReal), pointer :: &
           surfaceHeatFluxCategoryInitialArea, &
           surfaceConductiveFluxCategoryInitialArea, &
           latentHeatFluxCategoryInitialArea, &
           sensibleHeatFluxCategoryInitialArea

      real(kind=RKIND), dimension(:), pointer :: &
           iceAreaCellInitial, &
           snowfallRate, &
           rainfallRate, &
           absorbedShortwaveFlux, &
           latentHeatFlux, &
           sensibleHeatFlux, &
           longwaveUp, &
           evaporativeWaterFlux, &
           surfaceHeatFlux, &
           surfaceConductiveFlux

      real(kind=RKIND), dimension(:,:), pointer :: &
           iceAreaCategoryInitial, &
           surfaceHeatFluxCategory, &
           surfaceConductiveFluxCategory, &
           latentHeatFluxCategory, &
           sensibleHeatFluxCategory

      integer, pointer :: &
           nCellsSolve

      integer :: &
           iCell

      err = 0

      block => domain % blocklist
      do while (associated(block))

         call MPAS_pool_get_subpool(block % structs, "areaVariablesAM", areaVariablesAMPool)
         call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)

         call MPAS_pool_get_array(icestatePool, "iceAreaCellInitial", iceAreaCellInitial)
         call MPAS_pool_get_array(icestatePool, "iceAreaCategoryInitial", iceAreaCategoryInitial)

         call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)

         call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmosCouplingPool)
         call MPAS_pool_get_subpool(block % structs, "shortwave", shortwavePool)
         call MPAS_pool_get_subpool(block % structs, "atmos_fluxes", atmosFluxesPool)

         ! snowfallRate - fsnow
         call MPAS_pool_get_field(areaVariablesAMPool, "snowfallRateInitialArea", snowfallRateInitialArea)
         if (snowfallRateInitialArea % isActive) then
            call MPAS_pool_get_array(atmosCouplingPool, "snowfallRate", snowfallRate)
            do iCell = 1, nCellsSolve
               snowfallRateInitialArea % array(iCell) = snowfallRate(iCell) * iceAreaCellInitial(iCell)
            enddo ! iCell
         endif

         ! rainfallRate - frain
         call MPAS_pool_get_field(areaVariablesAMPool, "rainfallRateInitialArea", rainfallRateInitialArea)
         if (rainfallRateInitialArea % isActive) then
            call MPAS_pool_get_array(atmosCouplingPool, "rainfallRate", rainfallRate)
            do iCell = 1, nCellsSolve
               rainfallRateInitialArea % array(iCell) = rainfallRate(iCell) * iceAreaCellInitial(iCell)
            enddo ! iCell
         endif

         ! absorbedShortwaveFlux - fswabs
         call MPAS_pool_get_field(areaVariablesAMPool, "absorbedShortwaveFluxInitialArea", absorbedShortwaveFluxInitialArea)
         if (absorbedShortwaveFluxInitialArea % isActive) then
            call MPAS_pool_get_array(shortwavePool, "absorbedShortwaveFlux", absorbedShortwaveFlux)
            do iCell = 1, nCellsSolve
               absorbedShortwaveFluxInitialArea % array(iCell) = absorbedShortwaveFlux(iCell) * iceAreaCellInitial(iCell)
            enddo ! iCell
         endif

         ! latentHeatFlux - flat
         call MPAS_pool_get_field(areaVariablesAMPool, "latentHeatFluxInitialArea", latentHeatFluxInitialArea)
         if (latentHeatFluxInitialArea % isActive) then
            call MPAS_pool_get_array(atmosFluxesPool, "latentHeatFlux", latentHeatFlux)
            do iCell = 1, nCellsSolve
               latentHeatFluxInitialArea % array(iCell) = latentHeatFlux(iCell) * iceAreaCellInitial(iCell)
            enddo ! iCell
         endif

         ! latentHeatFluxCategory - flatn
         call MPAS_pool_get_field(areaVariablesAMPool, "latentHeatFluxCategoryInitialArea", latentHeatFluxCategoryInitialArea)
         if (latentHeatFluxCategoryInitialArea % isActive) then
            call MPAS_pool_get_array(atmosFluxesPool, "latentHeatFluxCategory", latentHeatFluxCategory)
            do iCell = 1, nCellsSolve
               latentHeatFluxCategoryInitialArea % array(:,iCell) = latentHeatFluxCategory(:,iCell) * iceAreaCategoryInitial(:,iCell)
            enddo ! iCell
         endif

         ! sensibleHeatFlux - fsens
         call MPAS_pool_get_field(areaVariablesAMPool, "sensibleHeatFluxInitialArea", sensibleHeatFluxInitialArea)
         if (sensibleHeatFluxInitialArea % isActive) then
            call MPAS_pool_get_array(atmosFluxesPool, "sensibleHeatFlux", sensibleHeatFlux)
            do iCell = 1, nCellsSolve
               sensibleHeatFluxInitialArea % array(iCell) = sensibleHeatFlux(iCell) * iceAreaCellInitial(iCell)
            enddo ! iCell
         endif

         ! sensibleHeatFluxCategory - fsensn
         call MPAS_pool_get_field(areaVariablesAMPool, "sensibleHeatFluxCategoryInitialArea", sensibleHeatFluxCategoryInitialArea)
         if (sensibleHeatFluxCategoryInitialArea % isActive) then
            call MPAS_pool_get_array(atmosFluxesPool, "sensibleHeatFluxCategory", sensibleHeatFluxCategory)
            do iCell = 1, nCellsSolve
               sensibleHeatFluxCategoryInitialArea % array(:,iCell) = sensibleHeatFluxCategory(:,iCell) * iceAreaCategoryInitial(:,iCell)
            enddo ! iCell
         endif

         ! longwaveUp - flwout
         call MPAS_pool_get_field(areaVariablesAMPool, "longwaveUpInitialArea", longwaveUpInitialArea)
         if (longwaveUpInitialArea % isActive) then
            call MPAS_pool_get_array(atmosFluxesPool, "longwaveUp", longwaveUp)
            do iCell = 1, nCellsSolve
               longwaveUpInitialArea % array(iCell) = longwaveUp(iCell) * iceAreaCellInitial(iCell)
            enddo ! iCell
         endif

         ! evaporativeWaterFlux - evap
         call MPAS_pool_get_field(areaVariablesAMPool, "evaporativeWaterFluxInitialArea", evaporativeWaterFluxInitialArea)
         if (evaporativeWaterFluxInitialArea % isActive) then
            call MPAS_pool_get_array(atmosFluxesPool, "evaporativeWaterFlux", evaporativeWaterFlux)
            do iCell = 1, nCellsSolve
               evaporativeWaterFluxInitialArea % array(iCell) = evaporativeWaterFlux(iCell) * iceAreaCellInitial(iCell)
            enddo ! iCell
         endif

         ! surfaceHeatFlux - fsurf
         call MPAS_pool_get_field(areaVariablesAMPool, "surfaceHeatFluxInitialArea", surfaceHeatFluxInitialArea)
         if (surfaceHeatFluxInitialArea % isActive) then
            call MPAS_pool_get_array(atmosFluxesPool, "surfaceHeatFlux", surfaceHeatFlux)
            do iCell = 1, nCellsSolve
               surfaceHeatFluxInitialArea % array(iCell) = surfaceHeatFlux(iCell) * iceAreaCellInitial(iCell)
            enddo ! iCell
         endif

         ! surfaceHeatFluxCategory - fsurfn
         call MPAS_pool_get_field(areaVariablesAMPool, "surfaceHeatFluxCategoryInitialArea", surfaceHeatFluxCategoryInitialArea)
         if (surfaceHeatFluxCategoryInitialArea % isActive) then
            call MPAS_pool_get_array(atmosFluxesPool, "surfaceHeatFluxCategory", surfaceHeatFluxCategory)
            do iCell = 1, nCellsSolve
               surfaceHeatFluxCategoryInitialArea % array(:,iCell) = surfaceHeatFluxCategory(:,iCell) * iceAreaCategoryInitial(:,iCell)
            enddo ! iCell
         endif

         ! surfaceConductiveFlux - fcondtop
         call MPAS_pool_get_field(areaVariablesAMPool, "surfaceConductiveFluxInitialArea", surfaceConductiveFluxInitialArea)
         if (surfaceConductiveFluxInitialArea % isActive) then
            call MPAS_pool_get_array(atmosFluxesPool, "surfaceConductiveFlux", surfaceConductiveFlux)
            do iCell = 1, nCellsSolve
               surfaceConductiveFluxInitialArea % array(iCell) = surfaceConductiveFlux(iCell) * iceAreaCellInitial(iCell)
            enddo ! iCell
         endif

         ! surfaceConductiveFluxCategory - fcondtopn
         call MPAS_pool_get_field(areaVariablesAMPool, "surfaceConductiveFluxCategoryInitialArea", surfaceConductiveFluxCategoryInitialArea)
         if (surfaceConductiveFluxCategoryInitialArea % isActive) then
            call MPAS_pool_get_array(atmosFluxesPool, "surfaceConductiveFluxCategory", surfaceConductiveFluxCategory)
            do iCell = 1, nCellsSolve
               surfaceConductiveFluxCategoryInitialArea % array(:,iCell) = surfaceConductiveFluxCategory(:,iCell) * iceAreaCategoryInitial(:,iCell)
            enddo ! iCell
         endif

         block => block % next
      enddo

   end subroutine cice_compute_area_variables!}}}

!***********************************************************************
!
!  routine cice_restart_area_variables
!
!> \brief   Save restart for MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    5th January 2016
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_restart_area_variables(domain, instance, err)!{{{

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

   end subroutine cice_restart_area_variables!}}}

!***********************************************************************
!
!  routine cice_finalize_area_variables
!
!> \brief   Finalize MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    5th January 2016
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_finalize_area_variables(domain, instance, err)!{{{

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

   end subroutine cice_finalize_area_variables!}}}

!-----------------------------------------------------------------------

end module cice_area_variables

! vim: foldmethod=marker
