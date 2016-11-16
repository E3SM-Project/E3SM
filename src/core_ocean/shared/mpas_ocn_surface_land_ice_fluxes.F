! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_surface_land_ice_fluxes
!
!> \brief MPAS ocean surface land-ice fluxes
!> \author Xylar Asay-Davis
!> \date   10/02/2014
!> \details
!>  This module contains routines for computing surface flux related
!>  melting under land-ice.
!
!-----------------------------------------------------------------------

module ocn_surface_land_ice_fluxes

   use mpas_timer
   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines

   use ocn_constants
   use ocn_equation_of_state

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

   public :: ocn_surface_land_ice_fluxes_tracers, &
             ocn_surface_land_ice_fluxes_vel, &
             ocn_surface_land_ice_fluxes_thick, &
             ocn_surface_land_ice_fluxes_build_arrays, &
             ocn_surface_land_ice_fluxes_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: landIceFluxesOn, standaloneOn, isomipOn, jenkinsOn, hollandJenkinsOn

   real (kind=RKIND) :: cp_land_ice, rho_land_ice


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_surface_land_ice_fluxes_tracers
!
!> \brief   Determines the tracers melt fluxes under land ice
!> \author  Xylar Asay-Davis
!> \date    9 September 2015
!> \details
!>  This routine adds land-ice tracer fluxes to the surface flux array
!>  used to compute tracer tendencies later in MPAS.
!
!-----------------------------------------------------------------------

   subroutine ocn_surface_land_ice_fluxes_tracers(meshPool, groupName, forcingPool, tracersSurfaceFlux, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: meshPool !< Input: mesh information
      character (len=*) :: groupName !< Input: Name of tracer group

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(inout) :: forcingPool !< Input: Forcing information
      real (kind=RKIND), dimension(:,:), intent(inout) :: tracersSurfaceFlux !< Input/Output: Surface flux for tracer group

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: Error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

      if ( .not. landIceFluxesOn ) return

      call mpas_timer_start("land_ice_" // trim(groupName))

      if ( trim(groupName) == 'activeTracers' ) then
         call ocn_surface_land_ice_fluxes_active_tracers(meshPool, forcingPool, tracersSurfaceFlux, err)
      end if

      call mpas_timer_stop("land_ice_" // trim(groupName))

   end subroutine ocn_surface_land_ice_fluxes_tracers!}}}

!***********************************************************************
!
!  routine ocn_surface_land_ice_fluxes_vel
!
!> \brief   Computes tendency term for top drag
!> \author  Xylar Asay-Davis
!> \date    9 September 2015
!> \details
!>  This routine computes the top-drag tendency for momentum
!>  based on current state.
!
!-----------------------------------------------------------------------

   subroutine ocn_surface_land_ice_fluxes_vel(meshPool, diagnosticsPool, surfaceStress, surfaceStressMagnitude, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: meshPool !< Input: mesh information
      type (mpas_pool_type), intent(in) :: diagnosticsPool !< Input: Diagnostics information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      real (kind=RKIND), dimension(:), intent(inout) :: surfaceStress, & !< Input/Output: Array for total surface stress
                                                  surfaceStressMagnitude !< Input/Output: Array for magnitude of surface stress

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: Error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------
      integer :: iEdge, iCell
      integer, pointer :: nCells, nEdges

      real (kind=RKIND), dimension(:), pointer :: topDrag, topDragMagnitude

      err = 0

      if ( .not. landIceFluxesOn ) return

      call mpas_timer_start("top_drag", .false.)

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)

      call mpas_pool_get_array(diagnosticsPool, 'topDrag', topDrag)
      call mpas_pool_get_array(diagnosticsPool, 'topDragMagnitude', topDragMagnitude)

      !$omp do schedule(runtime)
      do iEdge = 1, nEdges
        surfaceStress(iEdge) = surfaceStress(iEdge) + topDrag(iEdge)
      end do
      !$omp end do

      ! Build surface stress magnitude at cell centers
      !$omp do schedule(runtime)
      do iCell = 1, nCells
        surfaceStressMagnitude(iCell) = surfaceStressMagnitude(iCell) + topDragMagnitude(iCell)
      end do
      !$omp end do

      call mpas_timer_stop("top_drag")

   !--------------------------------------------------------------------

   end subroutine ocn_surface_land_ice_fluxes_vel!}}}

!***********************************************************************
!
!  routine ocn_surface_land_ice_fluxes_thick
!
!> \brief   Add land-ice fluxes to surfaceThicknessFlux.
!> \author  Xylar Asay-Davis
!> \date    11 September 2015
!> \details
!>  This routine adds land-ice freshwater fluxes to the surface thickness flux
!>  to be converted into a thickness tendency later.
!
!-----------------------------------------------------------------------

   subroutine ocn_surface_land_ice_fluxes_thick(meshPool, forcingPool, surfaceThicknessFlux, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: meshPool !< Input: mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(inout) :: forcingPool !< Input: Forcing information
      real (kind=RKIND), dimension(:), intent(inout) :: surfaceThicknessFlux !< Input/Output: Array for surface thickness flux

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: Error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iCell
      integer, pointer :: nCells

      real (kind=RKIND), dimension(:), pointer :: landIceFreshwaterFlux

      err = 0

      if ( .not. landIceFluxesOn ) return

      call mpas_timer_start("land_ice_thick")

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

      call mpas_pool_get_array(forcingPool, 'landIceFreshwaterFlux', landIceFreshwaterFlux)

      ! Build surface fluxes at cell centers
      !$omp do schedule(runtime)
      do iCell = 1, nCells
        surfaceThicknessFlux(iCell) = surfaceThicknessFlux(iCell) + landIceFreshwaterFlux(iCell) / rho_sw
      end do
      !$omp end do

      call mpas_timer_stop("land_ice_thick")

   end subroutine ocn_surface_land_ice_fluxes_thick!}}}

!***********************************************************************
!
!  routine ocn_surface_land_ice_fluxes_active_tracers
!
!> \brief   Adds the active tracers fluxes from land-ice melting.
!> \author  Xylar Asay-Davis
!> \date    11 September 2015
!> \details
!>  This routine adds the active tracers fluxes to surface fluxes
!>  from which tracer tendencies are computed later.
!
!-----------------------------------------------------------------------

   subroutine ocn_surface_land_ice_fluxes_active_tracers(meshPool, forcingPool, tracersSurfaceFlux, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: meshPool !< Input: mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(inout) :: forcingPool !< Input: Forcing information
      real (kind=RKIND), dimension(:,:), intent(inout) :: tracersSurfaceFlux

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: Error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iCell
      integer, pointer :: nCells

      real (kind=RKIND), dimension(:), pointer :: landIceHeatFlux

      err = 0

      if ( .not. landIceFluxesOn ) return

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

      call mpas_pool_get_array(forcingPool, 'landIceHeatFlux', landIceHeatFlux)

      ! add to surface fluxes at cell centers
      !$omp do schedule(runtime)
      do iCell = 1, nCells
        tracersSurfaceFlux(1, iCell) = tracersSurfaceFlux(1, iCell) + landIceHeatFlux(iCell)/(rho_sw*cp_sw)
      end do
      !$omp end do

   end subroutine ocn_surface_land_ice_fluxes_active_tracers!}}}


!***********************************************************************
!
!  routine ocn_surface_land_ice_fluxes_build_arrays
!
!> \brief Builds the forcing array for land-ice forcing
!> \author Xylar Asay-Davis
!> \date   10/02/2014
!> \details
!>  This routine computes surface fluxes related to land-ice forcing based
!>  on diagnostics from the previous time step.
!
!-----------------------------------------------------------------------

   subroutine ocn_surface_land_ice_fluxes_build_arrays(meshPool, diagnosticsPool, &
      forcingPool, scratchPool, statePool, dt, err) !{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool, &     !< Input: mesh information
         diagnosticsPool !< Input: diagnostics information
      real(kind=RKIND), intent(in) :: dt ! the time step over which to accumulate fluxes

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(inout) :: &
         forcingPool, & !< Input: Forcing information
         scratchPool, & !< Input: scratch field information
         statePool      !< Input: state field information

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: Error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), pointer :: tracersPool

      integer :: iCell
      integer, pointer :: nCellsSolve

      real (kind=RKIND), pointer :: config_land_ice_flux_ISOMIP_gammaT

      logical, pointer :: config_land_ice_flux_useHollandJenkinsAdvDiff


      real (kind=RKIND) :: freshwaterFlux, heatFlux

      real (kind=RKIND), dimension(:), pointer :: landIcePressure, landIceFraction, &
                                                  landIceSurfaceTemperature, &
                                                  landIceFrictionVelocity, &
                                                  landIceFreshwaterFlux, &
                                                  landIceHeatFlux, heatFluxToLandIce, &
                                                  freezeInterfaceSalinity, freezeInterfaceTemperature, &
                                                  freezeFreshwaterFlux, freezeHeatFlux, &
                                                  freezeIceHeatFlux
      real (kind=RKIND), dimension(:), pointer :: accumulatedLandIceMassOld, &
                                                  accumulatedLandIceMassNew, &
                                                  accumulatedLandIceHeatOld, &
                                                  accumulatedLandIceHeatNew

      integer, dimension(:), pointer :: landIceMask

      real (kind=RKIND), dimension(:,:), pointer :: landIceBoundaryLayerTracers, &
                                                    landIceInterfaceTracers, &
                                                    landIceTracerTransferVelocities
      integer, pointer :: indexBLT, indexBLS, indexIT, indexIS, indexHeatTrans, indexSaltTrans

      type (field1DReal), pointer :: boundaryLayerTemperatureField, boundaryLayerSalinityField, &
                                     freezeInterfaceSalinityField, freezeInterfaceTemperatureField, &
                                     freezeFreshwaterFluxField, freezeHeatFluxField, &
                                     freezeIceHeatFluxField

      err = 0

      if ( .not. standaloneOn ) return

      call mpas_timer_start("land_ice_build_arrays")

      call mpas_pool_get_config(ocnConfigs, 'config_land_ice_flux_ISOMIP_gammaT', config_land_ice_flux_ISOMIP_gammaT)
      call mpas_pool_get_config(ocnConfigs, 'config_land_ice_flux_useHollandJenkinsAdvDiff', &
                                config_land_ice_flux_useHollandJenkinsAdvDiff)

      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

      call mpas_pool_get_array(diagnosticsPool, 'landIceFrictionVelocity', landIceFrictionVelocity)

      call mpas_pool_get_array(forcingPool, 'landIcePressure', landIcePressure)
      call mpas_pool_get_array(diagnosticsPool, 'landIceBoundaryLayerTracers', landIceBoundaryLayerTracers)
      call mpas_pool_get_dimension(diagnosticsPool, 'index_landIceBoundaryLayerTemperature', indexBLT)
      call mpas_pool_get_dimension(diagnosticsPool, 'index_landIceBoundaryLayerSalinity', indexBLS)

      if(jenkinsOn .or. hollandJenkinsOn) then
        call mpas_pool_get_array(diagnosticsPool, 'landIceTracerTransferVelocities', landIceTracerTransferVelocities)
        call mpas_pool_get_dimension(diagnosticsPool, 'index_landIceHeatTransferVelocity', indexHeatTrans)
        call mpas_pool_get_dimension(diagnosticsPool, 'index_landIceSaltTransferVelocity', indexSaltTrans)
      end if

      call mpas_pool_get_array(forcingPool, 'landIceFraction', landIceFraction)
      call mpas_pool_get_array(forcingPool, 'landIceMask', landIceMask)

      call mpas_pool_get_array(forcingPool, 'landIceFreshwaterFlux', landIceFreshwaterFlux)
      call mpas_pool_get_array(forcingPool, 'landIceHeatFlux', landIceHeatFlux)
      call mpas_pool_get_array(forcingPool, 'heatFluxToLandIce', heatFluxToLandIce)

      call mpas_pool_get_array(forcingPool, 'landIceInterfaceTracers', landIceInterfaceTracers)
      call mpas_pool_get_dimension(forcingPool, 'index_landIceInterfaceTemperature', indexIT)
      call mpas_pool_get_dimension(forcingPool, 'index_landIceInterfaceSalinity', indexIS)
      call mpas_pool_get_array(statePool, 'accumulatedLandIceMass', accumulatedLandIceMassNew, 2)
      call mpas_pool_get_array(statePool, 'accumulatedLandIceMass', accumulatedLandIceMassOld, 1)
      call mpas_pool_get_array(statePool, 'accumulatedLandIceHeat', accumulatedLandIceHeatNew, 2)
      call mpas_pool_get_array(statePool, 'accumulatedLandIceHeat', accumulatedLandIceHeatOld, 1)

      if(config_land_ice_flux_useHollandJenkinsAdvDiff) then
         call mpas_pool_get_array(forcingPool, 'landIceSurfaceTemperature', landIceSurfaceTemperature)

         call mpas_pool_get_field(scratchPool, 'freezeInterfaceSalinityScratch', freezeInterfaceSalinityField)
         call mpas_pool_get_field(scratchPool, 'freezeInterfaceTemperatureScratch', freezeInterfaceTemperatureField)
         call mpas_pool_get_field(scratchPool, 'freezeFreshwaterFluxScratch', freezeFreshwaterFluxField)
         call mpas_pool_get_field(scratchPool, 'freezeHeatFluxScratch', freezeHeatFluxField)
         call mpas_pool_get_field(scratchPool, 'freezeIceHeatFluxScratch', freezeIceHeatFluxField)
         call mpas_allocate_scratch_field(freezeInterfaceSalinityField, .true.)
         call mpas_allocate_scratch_field(freezeInterfaceTemperatureField, .true.)
         call mpas_allocate_scratch_field(freezeFreshwaterFluxField, .true.)
         call mpas_allocate_scratch_field(freezeHeatFluxField, .true.)
         call mpas_allocate_scratch_field(freezeIceHeatFluxField, .true.)
         freezeInterfaceSalinity => freezeInterfaceSalinityField % array
         freezeInterfaceTemperature => freezeInterfaceTemperatureField % array
         freezeFreshwaterFlux => freezeFreshwaterFluxField % array
         freezeHeatFlux => freezeHeatFluxField % array
         freezeIceHeatFlux => freezeIceHeatFluxField % array
      end if

      if(isomipOn) then
         !$omp do schedule(runtime) private(heatFlux)
         do iCell = 1, nCellsSolve
            if (landIceMask(iCell) == 0) cycle

            ! linearized equaiton for the S and p dependent potential freezing temperature
            landIceInterfaceTracers(indexIT,iCell) = ocn_freezing_temperature( &
               salinity=landIceBoundaryLayerTracers(indexBLT,iCell), &
               pressure=landIcePressure(iCell), &
               inLandIceCavity=.true.)

            ! using (3) and (4) from Hunter (2006)
            ! or (7) from Jenkins et al. (2001) if gamma constant
            ! and no heat flux into ice
            ! freshwater flux = density * melt rate is in kg/m^2/s
            freshwaterFlux = -rho_sw * config_land_ice_flux_ISOMIP_gammaT * (cp_sw/latent_heat_fusion_mks) &
                       * (landIceInterfaceTracers(indexIT,iCell)-landIceBoundaryLayerTracers(indexBLT,iCell))

            landIceFreshwaterFlux(iCell) = landIceFraction(iCell)*freshwaterFlux

            ! Using (13) from Jenkins et al. (2001)
            ! heat flux is in W/s
            heatFlux = cp_sw*(freshwaterFlux*landIceInterfaceTracers(indexIT,iCell) &
                              + rho_sw*config_land_ice_flux_ISOMIP_gammaT &
                                * (landIceInterfaceTracers(indexIT,iCell)-landIceBoundaryLayerTracers(indexBLT,iCell)))
            landIceHeatFlux(iCell) = landIceFraction(iCell)*heatFlux

            heatFluxToLandIce(iCell) = 0.0_RKIND

         end do
         !$omp end do
      end if

      if(jenkinsOn .or. hollandJenkinsOn) then
         if(config_land_ice_flux_useHollandJenkinsAdvDiff) then
            ! melting solution
            call compute_HJ99_melt_fluxes( &
               landIceMask, &
               landIceBoundaryLayerTracers(indexBLT,:), &
               landIceBoundaryLayerTracers(indexBLS,:), &
               landIceTracerTransferVelocities(indexHeatTrans,:), &
               landIceTracerTransferVelocities(indexSaltTrans,:), &
               landIceSurfaceTemperature, &
               landIcePressure, &
               landIceInterfaceTracers(indexIT,:), &
               landIceInterfaceTracers(indexIS,:), &
               landIceFreshwaterFlux, &
               landIceHeatFlux, &
               heatFluxToLandIce, &
               nCellsSolve, &
               err)
            if(err .ne. 0) then
               call mpas_dmpar_global_abort("MPAS-ocean: ERROR: compute_HJ99_melt_fluxes failed.")
            end if

            ! freezing solution
            call compute_melt_fluxes( &
               landIceMask, &
               landIceBoundaryLayerTracers(indexBLT,:), &
               landIceBoundaryLayerTracers(indexBLS,:), &
               landIceTracerTransferVelocities(indexHeatTrans,:), &
               landIceTracerTransferVelocities(indexSaltTrans,:), &
               landIcePressure, &
               freezeInterfaceTemperature, &
               freezeInterfaceSalinity, &
               freezeFreshwaterFlux, &
               freezeHeatFlux, &
               freezeIceHeatFlux, &
               nCellsSolve, &
               err)
            if(err .ne. 0) then
               call mpas_dmpar_global_abort("MPAS-ocean: ERROR: compute_melt_fluxes failed.")
            end if

            do iCell = 1, nCellsSolve
               if ((landIceMask(iCell) == 0) .or. (landIceFreshwaterFlux(iCell) >= 0.0_RKIND)) cycle

               landIceInterfaceTracers(indexIS,iCell) = freezeInterfaceSalinity(iCell)
               landIceInterfaceTracers(indexIT,iCell) = freezeInterfaceTemperature(iCell)
               landIceFreshwaterFlux(iCell) = freezeFreshwaterFlux(iCell)
               landIceHeatFlux(iCell) = freezeHeatFlux(iCell)
               heatFluxToLandIce(iCell) = freezeIceHeatFlux(iCell)
            end do
         else ! not using Holland and Jenkins advection/diffusion
            call compute_melt_fluxes( &
               landIceMask, &
               landIceBoundaryLayerTracers(indexBLT,:), &
               landIceBoundaryLayerTracers(indexBLS,:), &
               landIceTracerTransferVelocities(indexHeatTrans,:), &
               landIceTracerTransferVelocities(indexSaltTrans,:), &
               landIcePressure, &
               landIceInterfaceTracers(indexIT,:), &
               landIceInterfaceTracers(indexIS,:), &
               landIceFreshwaterFlux, &
               landIceHeatFlux, &
               heatFluxToLandIce, &
               nCellsSolve, &
               err)
            if(err .ne. 0) then
               call mpas_dmpar_global_abort("MPAS-ocean: ERROR: compute_melt_fluxes failed.")
            end if
         end if

         ! modulate the fluxes by the landIceFraction
         do iCell = 1, nCellsSolve
            if (landIceMask(iCell) == 0) cycle

            landIceFreshwaterFlux(iCell) = landIceFraction(iCell)*landIceFreshwaterFlux(iCell)
            landIceHeatFlux(iCell) = landIceFraction(iCell)*landIceHeatFlux(iCell)
            heatFluxToLandIce(iCell) = landIceFraction(iCell)*heatFluxToLandIce(iCell)
         end do

      end if

      if(config_land_ice_flux_useHollandJenkinsAdvDiff) then
         call mpas_deallocate_scratch_field(freezeInterfaceSalinityField, .true.)
         call mpas_deallocate_scratch_field(freezeInterfaceTemperatureField, .true.)
         call mpas_deallocate_scratch_field(freezeFreshwaterFluxField, .true.)
         call mpas_deallocate_scratch_field(freezeHeatFluxField, .true.)
         call mpas_deallocate_scratch_field(freezeIceHeatFluxField, .true.)
      end if

      ! accumulate land-ice mass and heat
      do iCell = 1, nCellsSolve
        accumulatedLandIceMassNew(iCell) = accumulatedLandIceMassOld(iCell) - dt*landIceFreshwaterFlux(iCell)
        accumulatedLandIceHeatNew(iCell) = accumulatedLandIceHeatOld(iCell) + dt*heatFluxToLandIce(iCell)
      end do

      call mpas_timer_stop("land_ice_build_arrays")

   !--------------------------------------------------------------------

   end subroutine ocn_surface_land_ice_fluxes_build_arrays!}}}

!***********************************************************************
!
!  routine ocn_surface_land_ice_fluxes_init
!
!> \brief Initializes land-ice forcing
!> \author Xylar Asay-Davis
!> \date   10/02/2014
!> \details
!>  This routine initializes a variety of quantities related to
!>  land-ice forcing.
!
!-----------------------------------------------------------------------

   subroutine ocn_surface_land_ice_fluxes_init(err)!{{{

      integer, intent(out) :: err !< Output: error flag

      character (len=StrKIND), pointer :: config_land_ice_flux_formulation, config_land_ice_flux_mode

      real (kind=RKIND), pointer :: config_land_ice_flux_cp_ice, &
                                    config_land_ice_flux_rho_ice


      err = 0
      isomipOn = .false.
      jenkinsOn = .false.
      hollandJenkinsOn = .false.

      call mpas_pool_get_config(ocnConfigs, 'config_land_ice_flux_mode', config_land_ice_flux_mode)
      landIceFluxesOn = (trim(config_land_ice_flux_mode) == 'standalone')  &
           .or. (trim(config_land_ice_flux_mode) == 'coupled')
      if(.not. landIceFluxesOn) return

      standaloneOn = trim(config_land_ice_flux_mode) == 'standalone'

      call mpas_pool_get_config(ocnConfigs, 'config_land_ice_flux_formulation', config_land_ice_flux_formulation)

      call mpas_pool_get_config(ocnConfigs, 'config_land_ice_flux_cp_ice', config_land_ice_flux_cp_ice)
      call mpas_pool_get_config(ocnConfigs, 'config_land_ice_flux_rho_ice', config_land_ice_flux_rho_ice)

      if ( trim(config_land_ice_flux_formulation) == 'ISOMIP' ) then
         isomipOn = .true.
      else if ( trim(config_land_ice_flux_formulation) == 'Jenkins' ) then
         jenkinsOn = .true.
      else if ( trim(config_land_ice_flux_formulation) == 'HollandJenkins' ) then
         hollandJenkinsOn = .true.
      else
         write(stderrUnit, *) "ERROR: config_land_ice_flux_formulation not one of 'ISOMIP', 'Jenkins', or 'HollandJenkins'."
         err = 1
         call mpas_dmpar_global_abort("MPAS-ocean: ERROR: config_land_ice_flux_formulation not one of 'ISOMIP', 'Jenkins', " &
                                      // "or 'HollandJenkins'.")
      end if

      cp_land_ice = config_land_ice_flux_cp_ice
      rho_land_ice = config_land_ice_flux_rho_ice

   !--------------------------------------------------------------------

   end subroutine ocn_surface_land_ice_fluxes_init!}}}

!***********************************************************************
!
!  routine ocn_forcing_compute_melt_fluxes
!
!> \brief   Computes ocean and ice melt fluxes, etc.
!> \author  Xylar Asay-Davis
!> \date    3/27/2015
!>  This routine computes melt fluxes (melt rate, temperature fluxes
!>  into the ice and the ocean, and salt flux) as well as the interface
!>  temperature and salinity.  This routine expects an ice temperature
!>  in the bottom layer of ice and ocean temperature and salinity in
!>  the top ocean layer as well as the pressure at the ice/ocean interface.
!>
!>  The ocean heat and salt transfer velocities are determined based on
!>  observations of turbulent mixing rates in the under-ice boundary layer.
!>  They should be the product of the friction velocity and a (possibly
!>  spatially variable) non-dimenional transfer coefficient.
!>
!>  The iceTemperatureDistance is the distance between the location
!>  where the iceTemperature is supplied and the ice-ocean interface,
!>  used to compute a temperature gradient.  The ice thermal conductivity,
!>  kappa_land_ice, is zero for the freezing solution from Holland and Jenkins
!>  (1999) in which the ice is purely insulating.
!
!-----------------------------------------------------------------------


  subroutine compute_melt_fluxes( &
    mask, &
    oceanTemperature, &
    oceanSalinity, &
    oceanHeatTransferVelocity, &
    oceanSaltTransferVelocity, &
    interfacePressure, &
    outInterfaceTemperature, &
    outInterfaceSalinity, &
    outFreshwaterFlux, &
    outOceanHeatFlux, &
    outIceHeatFlux, &
    nCells, &
    err, &
    iceTemperature, &
    iceTemperatureDistance, &
    kappa_land_ice) !{{{

    !-----------------------------------------------------------------
    !
    ! input variables
    !
    !-----------------------------------------------------------------

    integer, dimension(:), intent(in) :: &
      mask                         !< Input: mask for land-ice fluxes

    real (kind=RKIND), dimension(:), intent(in) :: &
      oceanTemperature, &          !< Input: ocean temperature in top layer
      oceanSalinity, &             !< Input: ocean salinity in top layer
      oceanHeatTransferVelocity, & !< Input: ocean heat transfer velocity
      oceanSaltTransferVelocity, & !< Input: ocean salt transfer velocity
      interfacePressure            !< Input: pressure at the ice-ocean interface

    integer, intent(in) :: nCells !< Input: number of cells in each array

    real (kind=RKIND), dimension(:), intent(in), optional:: &
      iceTemperature, &            !< Input: ice temperature in bottom layer
      iceTemperatureDistance       !< Input: distance to ice temperature from ice-ocean interface

    real (kind=RKIND), intent(in), optional:: &
      kappa_land_ice     !< Input: the diffusivity of heat in land ice

    !-----------------------------------------------------------------
    !
    ! output variables
    !
    !-----------------------------------------------------------------

    real (kind=RKIND), dimension(:), intent(out) :: &
      outInterfaceTemperature, & !< Output: ice/ocean temperature at the interface
      outInterfaceSalinity, &    !< Output: ocean salinity at the interface
      outFreshwaterFlux, &   !< Output: ocean thickness flux (melt rate)
      outOceanHeatFlux, & !< Output: the temperature flux into the ocean
      outIceHeatFlux      !< Output: the temperature flux into the ice

    integer, intent(out) :: err !< Output: Error flag

    !-----------------------------------------------------------------
    !
    ! local variables
    !
    !-----------------------------------------------------------------

    real (kind=RKIND) :: T0, transferVelocityRatio, Tlatent, nu, a, b, c, eta, &
                         iceHeatFluxCoeff, iceDeltaT, dTf_dS
    integer :: iCell

    logical :: coupled

    real (kind=RKIND), parameter :: minInterfaceSalinity = 0.001_RKIND

    err = 0
    coupled = present(iceTemperature) .and. present(iceTemperatureDistance) &
            .and. present(kappa_land_ice)
    Tlatent = latent_heat_fusion_mks/cp_sw

    !$omp do schedule(runtime) private(iceHeatFluxCoeff, nu, iceDeltaT, T0, transferVelocityRatio, a, b, c)
    do iCell = 1, nCells
      if (mask(iCell) == 0) cycle

      if(coupled) then
        iceHeatFluxCoeff = rho_land_ice*cp_land_ice*kappa_land_ice/iceTemperatureDistance(iCell)
        nu = iceHeatFluxCoeff/(rho_sw*cp_sw*oceanHeatTransferVelocity(iCell))
        iceDeltaT = T0 - iceTemperature(iCell)
      else
        nu = 0.0_RKIND
        iceDeltaT = 0.0_RKIND
      end if
      T0 = ocn_freezing_temperature(salinity=0.0_RKIND, pressure=interfacePressure(iCell), &
                                    inLandIceCavity=.true.)
      dTf_dS = ocn_freezing_temperature_salinity_deriv(salinity=0.0_RKIND, pressure=interfacePressure(iCell), &
                                                       inLandIceCavity=.true.)

      transferVelocityRatio = oceanSaltTransferVelocity(iCell)/oceanHeatTransferVelocity(iCell)

      a = -dTf_dS*(1.0_RKIND + nu)
      b = transferVelocityRatio*Tlatent - nu*iceDeltaT + oceanTemperature(iCell) - T0
      c = -transferVelocityRatio*Tlatent*max(oceanSalinity(iCell), 0.0_RKIND)

      ! a is non-negative; c is strictly non-positive so we never get imaginary roots.
      ! Since a can be zero, we need a solution of the quadratic equation for 1/Si instead of Si.
      ! Following: https://people.csail.mit.edu/bkph/articles/Quadratics.pdf
      ! Since a and -c are are non-negative, the term in the square root is also always >= |b|.
      ! In all reasonable cases, b will be strictly positive, since transferVelocityRatio*Tlatent ~ 2 C,
      ! T0 ~ -1.8 C and oceanTemperature should never be able to get below about -3 C
      ! As long as either b or both a and c are greater than zero, the strictly non-negative root is
      outInterfaceSalinity(iCell) = max(-(2.0_RKIND*c)/(b + sqrt(b**2 - 4.0_RKIND*a*c)), minInterfaceSalinity)

      outInterfaceTemperature(iCell) = dTf_dS*outInterfaceSalinity(iCell)+T0

      outFreshwaterFlux(iCell) = rho_sw*oceanSaltTransferVelocity(iCell) &
        * (oceanSalinity(iCell)/outInterfaceSalinity(iCell) - 1.0_RKIND)

      ! According to Jenkins et al. (2001), the temperature fluxes into the ocean are:
      !   1. the advection of meltwater into the top layer (or removal for freezing)
      !   2. the turbulent transfer of heat across the boundary layer, based on the termal driving
      outOceanHeatFlux(iCell) = cp_sw*(outFreshwaterFlux(iCell)*outInterfaceTemperature(iCell) &
        - rho_sw*oceanHeatTransferVelocity(iCell)*(oceanTemperature(iCell)-outInterfaceTemperature(iCell)))

      ! the temperature fluxes into the ice are:
      !   1. the advection of ice at the interface temperature out of the domain due to melting
      !      (or in due to freezing)
      !   2. the diffusion (if any) of heat into the ice, based on temperature difference between
      !      the reference point in the ice (either the surface or the middle of the bottom layer)
      !      and the interface
      outIceHeatFlux(iCell) = -cp_land_ice*outFreshwaterFlux(iCell)*outInterfaceTemperature(iCell)
      if(coupled) then
        outIceHeatFlux(iCell) = outIceHeatFlux(iCell) &
          - iceHeatFluxCoeff*(iceTemperature(iCell) - outInterfaceTemperature(iCell))
      end if
    end do
    !$omp end do

  !--------------------------------------------------------------------

  end subroutine compute_melt_fluxes !}}}


!***********************************************************************
!
!  routine compute_HJ99_melt_fluxes
!
!> \brief   Computes melt fluxes, etc. according to HJ99
!> \author  Xylar Asay-Davis
!> \date    3/28/2015
!> \details
!>  This routine computes melt fluxes (melt rate, temperature fluxes
!>  into the ice and the ocean, and salt flux) as well as the interface
!>  temperature and salinity.  Following Holland and Jenkins (1999),
!>  temperature is assumed to be vertically advected and diffused in
!>  the ice at a rate determined by the melt rate, so that no
!>  heat transfer velocity for the ice need be supplied.  Except for
!>  very small melt rates, the Holland and Jenkins advection/diffusion
!>  solution produces an ice temperature profile that is approximately
!>  constant with depth except near the ice-ocean interface.  The ice
!>  temperature supplied to this routine should be the far-field value,
!>  equal to the time-averaged surface temperature.
!>
!>  The solution is only appropriate for melting (positive ocean
!>  thickness flux). For freezing, the fluxes should be computed using
!>  ocn_forcing_compute_melt_fluxes with ``insulating'' ice where
!>  the iceHeatTransferVelocity is set to zero.
!
!-----------------------------------------------------------------------

  subroutine compute_HJ99_melt_fluxes( &
    mask, &
    oceanTemperature, &
    oceanSalinity, &
    oceanHeatTransferVelocity, &
    oceanSaltTransferVelocity, &
    iceTemperature, &
    interfacePressure, &
    outInterfaceTemperature, &
    outInterfaceSalinity, &
    outFreshwaterFlux, &
    outOceanHeatFlux, &
    outIceHeatFlux, &
    nCells, &
    err) !{{{

    !-----------------------------------------------------------------
    !
    ! input variables
    !
    !-----------------------------------------------------------------

    integer, dimension(:), intent(in) :: &
      mask                         !< Input: mask for land-ice fluxes

    real (kind=RKIND), dimension(:), intent(in) :: &
      oceanTemperature, &          !< Input: ocean temperature in top layer
      oceanSalinity, &             !< Input: ocean salinity in top layer
      oceanHeatTransferVelocity, & !< Input: ocean heat transfer velocity
      oceanSaltTransferVelocity, & !< Input: ocean salt transfer velocity
      iceTemperature, &            !< Input: ice temperature in bottom layer
      interfacePressure            !< Input: pressure at the ice-ocean interface

    integer, intent(in) :: nCells !< Input: number of cells in each array

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

    real (kind=RKIND), dimension(:), intent(out) :: &
      outInterfaceTemperature, & !< Output: ice/ocean temperature at the interface
      outInterfaceSalinity, &    !< Output: ocean salinity at the interface
      outFreshwaterFlux, &   !< Output: ocean thickness flux (melt rate)
      outOceanHeatFlux, & !< Output: the temperature flux into the ocean
      outIceHeatFlux      !< Output: the temperature flux into the ice

    integer, intent(out) :: err !< Output: Error flag

    !-----------------------------------------------------------------
    !
    ! local variables
    !
    !-----------------------------------------------------------------

    real (kind=RKIND) :: T0, cpRatio, transferVelocityRatio, Tlatent, a, b, c, eta, TlatentStar, dTf_dS

    integer :: iCell

    err = 0
    cpRatio = cp_land_ice/cp_sw
    !$omp do schedule(runtime) private(T0, transferVelocityRatio, Tlatent, eta, TlatentStar, a, b, c)
    do iCell = 1, nCells
      if (mask(iCell) == 0) cycle

      T0 = ocn_freezing_temperature(salinity=0.0_RKIND, pressure=interfacePressure(iCell), &
                                    inLandIceCavity=.true.)
      dTf_dS = ocn_freezing_temperature_salinity_deriv(salinity=0.0_RKIND, pressure=interfacePressure(iCell), &
                                                       inLandIceCavity=.true.)
      transferVelocityRatio = (rho_fw/rho_sw)*oceanSaltTransferVelocity(iCell)/oceanHeatTransferVelocity(iCell)
      Tlatent = latent_heat_fusion_mks/cp_sw

      eta = cpRatio * transferVelocityRatio
      TlatentStar = Tlatent + cpRatio*(T0-iceTemperature(iCell))
      a = -dTf_dS*(1.0_RKIND - eta)
      b = (transferVelocityRatio*TlatentStar - eta*dTf_dS*oceanSalinity(iCell) &
        + oceanTemperature(iCell) - T0)
      c = -transferVelocityRatio*TlatentStar

      ! a is strictly positive; c is strictly negative so we never get imaginary roots
      ! The positive root is the one we want (salinity is strictly positive)
      outInterfaceSalinity(iCell) = (-b + sqrt(b**2 - 4.0_RKIND*a*c*oceanSalinity(iCell)))/(2.0_RKIND*a)
      if (outInterfaceSalinity(iCell) .le. 0.0_RKIND) then
        err = 1
        call mpas_dmpar_global_abort('MPAS-ocean: ERROR: interfaceSalinity is negative...')
      end if
      outInterfaceTemperature(iCell) = dTf_dS*outInterfaceSalinity(iCell)+T0

      outFreshwaterFlux(iCell) = rho_sw*oceanSaltTransferVelocity(iCell) &
        * (oceanSalinity(iCell)/outInterfaceSalinity(iCell) - 1.0_RKIND)

      ! According to Jenkins et al. (2001), the temperature fluxes into the ocean are:
      !   1. the advection of meltwater into the top layer (or removal for freezing)
      !   2. the turbulent transfer of heat across the boundary layer, based on the termal driving
      outOceanHeatFlux(iCell) = cp_sw*(outFreshwaterFlux(iCell)*outInterfaceTemperature(iCell) &
        - rho_sw*oceanHeatTransferVelocity(iCell)*(oceanTemperature(iCell)-outInterfaceTemperature(iCell)))

      ! Since we're considering only melting and ignoring diffusion,
      ! the ice loses heat simply by the loss of ice mass at the prescribed
      ! (surface?) ice temperature
      outIceHeatFlux(iCell) = -cp_land_ice*outFreshwaterFlux(iCell)*iceTemperature(iCell)
    end do
    !$omp end do

  !--------------------------------------------------------------------

  end subroutine compute_HJ99_melt_fluxes !}}}


!***********************************************************************

end module ocn_surface_land_ice_fluxes

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
