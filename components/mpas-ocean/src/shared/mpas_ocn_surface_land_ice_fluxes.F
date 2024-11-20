! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_surface_land_ice_fluxes
!
!> \brief MPAS ocean surface land-ice fluxes
!> \author Xylar Asay-Davis
!> \date   29 Jan 2023
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
   use ocn_config
   use ocn_mesh
   use ocn_equation_of_state
   use ocn_diagnostics_variables

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   logical, public :: &
      landIceFluxesOn,   &! on/off switch for land ice fluxes
      landIcePressureOn, &! on/off switch for land ice pressure
      isomipOn,          &! use ISOMIP flux formulation
      jenkinsOn,         &! use Jenkins flux formulation
      hollandJenkinsOn    ! use Holland-Jenkins flux formulation

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_surface_land_ice_fluxes_tracers, &
             ocn_surface_land_ice_fluxes_vel, &
             ocn_surface_land_ice_fluxes_thick, &
             ocn_surface_land_ice_fluxes_build_arrays, &
             ocn_surface_land_ice_fluxes_accumulate_fluxes, &
             ocn_surface_land_ice_fluxes_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: &
      landIceDataOn ,         &! flag to determine if using data (prescribed) fluxes
      landIceStandaloneOn,    &! flag to determine if using standalone fluxes
      landIceCoupledOn,       &! flag to determine if using fluxes from the coupler
      useHollandJenkinsAdvDiff ! use Holland-Jenkins advection-diffusion

   real (kind=RKIND) :: &
      cpLandIce,        &! specific heat for land ice
      rhoLandIce,       &! density of land ice
      ISOMIPgammaT       ! parameter for ISOMIP formulation

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

      if (.not.landIceFluxesOn) return

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
!>  This routine computes the top-drag tendency under land ice for
!>  momentum based on current state.
!
!-----------------------------------------------------------------------

   subroutine ocn_surface_land_ice_fluxes_vel(&
                                      sfcStress, sfcStressMag, err)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(inout) :: &
         sfcStress,  &!< [inout] accumulated surface stress
         sfcStressMag !< [inout] accumulated surface stress magnitude

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< [out] Error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: &
         iEdge, iCell  ! loop indices for edge, cell loops

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !*** set error flag and return early if land ice fluxes not on
      !*** otherwise, start timer

      err = 0
      if (.not.landIceFluxesOn) return
      if (config_use_implicit_top_drag) return

      call mpas_timer_start("top_drag", .false.)

      !*** Transfer data to device
      !$acc update device (topDrag, topDragMag)

      !*** simply add input values to accumulated stresses

#ifdef MPAS_OPENACC
      !$acc parallel loop &
      !$acc    present(sfcStress, topDrag)
#else
      !$omp parallel
      !$omp do schedule(runtime)
#endif
      do iEdge = 1, nEdgesAll
        sfcStress(iEdge) = sfcStress(iEdge) + topDrag(iEdge)
      end do
#ifndef MPAS_OPENACC
      !$omp end do
#endif

#ifdef MPAS_OPENACC
      !$acc parallel loop &
      !$acc    present(sfcStressMag, topDragMag)
#else
      !$omp do schedule(runtime)
#endif
      do iCell = 1, nCellsAll
        sfcStressMag(iCell) = sfcStressMag(iCell) + topDragMag(iCell)
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      !*** Transfer any needed data back to host and delete inputs

      call mpas_timer_stop("top_drag")

   !--------------------------------------------------------------------

   end subroutine ocn_surface_land_ice_fluxes_vel!}}}

!***********************************************************************
!
!  routine ocn_surface_land_ice_fluxes_thick
!
!> \brief   Add land-ice fluxes to surfaceThicknessFlux.
!> \author  Xylar Asay-Davis
!> \date    29 Jan 2023
!> \details
!>  This routine adds land-ice freshwater fluxes to the surface thickness flux
!>  to be converted into a thickness tendency later.
!
!-----------------------------------------------------------------------

   subroutine ocn_surface_land_ice_fluxes_thick(forcingPool, surfaceThicknessFlux, err)!{{{

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

      real (kind=RKIND), dimension(:), pointer :: landIceFreshwaterFlux, &
                                                  dataLandIceFreshwaterFlux

      err = 0

      if (.not.landIceFluxesOn) return

      call mpas_timer_start("land_ice_thick")

      call mpas_pool_get_array(forcingPool, 'landIceFreshwaterFlux', landIceFreshwaterFlux)

      if ( landIceDataOn ) then
         call mpas_pool_get_array(forcingPool, 'dataLandIceFreshwaterFlux', dataLandIceFreshwaterFlux)
#ifdef MPAS_OPENACC
         !$acc enter data copyin(landIceFreshwaterFlux)

         !$acc parallel loop present(dataLandIceFreshwaterFlux)
#else
         !$omp parallel
         !$omp do schedule(runtime)
#endif
         do iCell = 1, nCellsAll
            landIceFreshwaterFlux(iCell) = dataLandIceFreshwaterFlux(iCell)
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif
      end if
      ! Build surface fluxes at cell centers
#ifdef MPAS_OPENACC
      !$acc enter data copyin(landIceFreshwaterFlux)

      !$acc parallel loop present(surfaceThicknessFlux, landIceFreshwaterFlux)
#else
      !$omp parallel
      !$omp do schedule(runtime)
#endif
      do iCell = 1, nCellsAll
         surfaceThicknessFlux(iCell) = surfaceThicknessFlux(iCell) + landIceFreshwaterFlux(iCell) / rho_sw
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

#ifdef MPAS_OPENACC
      !$acc exit data delete(landIceFreshwaterFlux)
#endif
      call mpas_timer_stop("land_ice_thick")

   end subroutine ocn_surface_land_ice_fluxes_thick!}}}

!***********************************************************************
!
!  routine ocn_surface_land_ice_fluxes_active_tracers
!
!> \brief   Adds the active tracers fluxes from land-ice melting.
!> \author  Xylar Asay-Davis
!> \date    29 Jan 2023
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

      real (kind=RKIND), dimension(:), pointer :: landIceHeatFlux, &
                                                  dataLandIceHeatFlux

      err = 0

      if (.not.landIceFluxesOn) return

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

      call mpas_pool_get_array(forcingPool, 'landIceHeatFlux', landIceHeatFlux)

      if ( landIceDataOn ) then
         call mpas_pool_get_array(forcingPool, 'dataLandIceHeatFlux', dataLandIceHeatFlux)
         !$omp parallel
         !$omp do schedule(runtime)
         do iCell = 1, nCells
            landIceHeatFlux(iCell) = dataLandIceHeatFlux(iCell)
         end do
         !$omp end do
         !$omp end parallel
      end if

      ! add to surface fluxes at cell centers
      !$omp parallel
      !$omp do schedule(runtime)
      do iCell = 1, nCells
         tracersSurfaceFlux(1, iCell) = tracersSurfaceFlux(1, iCell) + landIceHeatFlux(iCell)/(rho_sw*cp_sw)
      end do
      !$omp end do
      !$omp end parallel

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

   subroutine ocn_surface_land_ice_fluxes_build_arrays(meshPool, &
      forcingPool, scratchPool, statePool, err) !{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool        !< Input: mesh information

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

      integer :: iCell, nCells
      integer, dimension(:), pointer :: nCellsArray

      real (kind=RKIND) :: freshwaterFlux, heatFlux

      real (kind=RKIND), dimension(:), pointer :: landIcePressure, landIceFloatingFraction, &
                                                  landIceSurfaceTemperature, &
                                                  landIceFreshwaterFlux, &
                                                  landIceHeatFlux, heatFluxToLandIce

      real (kind=RKIND), dimension(:), pointer :: frazilIceFreshwaterFlux, landIceFreshwaterFluxTotal

      integer, dimension(:), pointer :: landIceFloatingMask

      real (kind=RKIND), dimension(:,:), pointer :: &
                                                    landIceInterfaceTracers

      integer, pointer :: indexITPtr, indexISPtr
      integer :: indexIT, indexIS

      ! Scratch Arrays
      !---
      ! freezeInterfaceSalinity: salinity at land ice-ocean interface where freezing is occurring
      !                   units: PSU
      real (kind=RKIND), dimension(:), allocatable :: freezeInterfaceSalinity
      ! freezeInterfaceTemperature: temperature at land ice-ocean interface where freezing is occurring
      !                      units: degrees C
      real (kind=RKIND), dimension(:), allocatable :: freezeInterfaceTemperature
      real (kind=RKIND), dimension(:), allocatable :: freezeFreshwaterFlux
      real (kind=RKIND), dimension(:), allocatable :: freezeHeatFlux
      real (kind=RKIND), dimension(:), allocatable :: freezeIceHeatFlux

      err = 0

      if (.not. (landIceStandaloneOn .or. landIceDataOn)) return

      call mpas_timer_start("land_ice_build_arrays")

      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)

      call mpas_pool_get_array(forcingPool, 'landIceFreshwaterFlux', landIceFreshwaterFlux)
      call mpas_pool_get_array(forcingPool, 'frazilIceFreshwaterFlux', frazilIceFreshwaterFlux)
      call mpas_pool_get_array(forcingPool, 'landIceFreshwaterFluxTotal', landIceFreshwaterFluxTotal)
      call mpas_pool_get_array(forcingPool, 'landIceFloatingMask', landIceFloatingMask)

      nCells = nCellsArray( size(nCellsArray) )

      if (landIceStandaloneOn) then

         call mpas_pool_get_array(forcingPool, 'landIcePressure', landIcePressure)
         call mpas_pool_get_array(forcingPool, 'landIceFloatingFraction', landIceFloatingFraction)
         call mpas_pool_get_array(forcingPool, 'landIceHeatFlux', landIceHeatFlux)
         call mpas_pool_get_array(forcingPool, 'heatFluxToLandIce', heatFluxToLandIce)


         call mpas_pool_get_array(forcingPool, 'landIceInterfaceTracers', landIceInterfaceTracers)
         call mpas_pool_get_dimension(forcingPool, &
                                 'index_landIceInterfaceTemperature', &
                                    indexITPtr)
         call mpas_pool_get_dimension(forcingPool, &
                                 'index_landIceInterfaceSalinity', &
                                    indexISPtr)
         indexIT = indexITPtr
         indexIS = indexISPtr

         if (useHollandJenkinsAdvDiff) then
            call mpas_pool_get_array(forcingPool, 'landIceSurfaceTemperature', landIceSurfaceTemperature)

            allocate(freezeInterfaceSalinity(nCells), &
                     freezeInterfaceTemperature(nCells), &
                     freezeFreshwaterFlux(nCells), &
                     freezeHeatFlux(nCells), &
                     freezeIceHeatFlux(nCells))
         end if

         if (isomipOn) then !*** ISOMIP formulation

            !$omp parallel
            !$omp do schedule(runtime) private(freshwaterFlux, heatFlux)
            do iCell = 1, nCells
               if (landIceFloatingMask(iCell) == 0) cycle

               ! linearized equaiton for the S and p dependent potential freezing temperature
               landIceInterfaceTracers(indexIT,iCell) = ocn_freezing_temperature( &
                  salinity=landIceBoundaryLayerTracers(indexBLT,iCell), &
                  pressure=landIcePressure(iCell), &
                  inLandIceCavity=.true.)

               ! using (3) and (4) from Hunter (2006)
               ! or (7) from Jenkins et al. (2001) if gamma constant
               ! and no heat flux into ice
               ! freshwater flux = density * melt rate is in kg/m^2/s
               freshwaterFlux = -rho_sw * ISOMIPgammaT * (cp_sw/latent_heat_fusion_mks) &
                        * (landIceInterfaceTracers(indexIT,iCell)-landIceBoundaryLayerTracers(indexBLT,iCell))

               landIceFreshwaterFlux(iCell) = landIceFloatingFraction(iCell)*freshwaterFlux

               ! Using (13) from Jenkins et al. (2001)
               ! heat flux is in W/s
               heatFlux = cp_sw*(freshwaterFlux*landIceInterfaceTracers(indexIT,iCell) &
                                 + rho_sw*ISOMIPgammaT &
                                 * (landIceInterfaceTracers(indexIT,iCell)-landIceBoundaryLayerTracers(indexBLT,iCell)))
               landIceHeatFlux(iCell) = landIceFloatingFraction(iCell)*heatFlux

               heatFluxToLandIce(iCell) = 0.0_RKIND

            end do
            !$omp end do
            !$omp end parallel
         end if ! isomipOn

         if (jenkinsOn .or. hollandJenkinsOn) then
            if(useHollandJenkinsAdvDiff) then
               ! melting solution
               call compute_HJ99_melt_fluxes( &
                  landIceFloatingMask, &
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
                  nCells, &
                  err)
               if(err .ne. 0) then
                  call mpas_log_write( &
                     'compute_HJ99_melt_fluxes failed.', &
                     MPAS_LOG_CRIT)
               end if

               ! freezing solution
               call compute_melt_fluxes( &
                  landIceFloatingMask, &
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
                  nCells, &
                  err)
               if(err .ne. 0) then
                  call mpas_log_write( &
                     'compute_melt_fluxes failed.', &
                     MPAS_LOG_CRIT)
               end if

               do iCell = 1, nCells
                  if ((landIceFloatingMask(iCell) == 0) .or. (landIceFreshwaterFlux(iCell) >= 0.0_RKIND)) cycle

                  landIceInterfaceTracers(indexIS,iCell) = freezeInterfaceSalinity(iCell)
                  landIceInterfaceTracers(indexIT,iCell) = freezeInterfaceTemperature(iCell)
                  landIceFreshwaterFlux(iCell) = freezeFreshwaterFlux(iCell)
                  landIceHeatFlux(iCell) = freezeHeatFlux(iCell)
                  heatFluxToLandIce(iCell) = freezeIceHeatFlux(iCell)
               end do
            else ! not using Holland and Jenkins advection/diffusion
               call compute_melt_fluxes( &
                  landIceFloatingMask, &
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
                  nCells, &
                  err)
               if(err .ne. 0) then
                  call mpas_log_write( &
                     'compute_melt_fluxes failed.', &
                     MPAS_LOG_CRIT)
               end if
            end if

            ! modulate the fluxes by the landIceFloatingFraction
            do iCell = 1, nCells
               if (landIceFloatingMask(iCell) == 0) cycle

               landIceFreshwaterFlux(iCell) = landIceFloatingFraction(iCell)*landIceFreshwaterFlux(iCell)
               landIceHeatFlux(iCell) = landIceFloatingFraction(iCell)*landIceHeatFlux(iCell)
               heatFluxToLandIce(iCell) = landIceFloatingFraction(iCell)*heatFluxToLandIce(iCell)
            end do

         end if ! jenkinsOn or hollandJenkinsOn

         if(useHollandJenkinsAdvDiff) then
            deallocate(freezeInterfaceSalinity, &
                     freezeInterfaceTemperature, &
                     freezeFreshwaterFlux, &
                     freezeHeatFlux, &
                     freezeIceHeatFlux)
         end if

      end if ! landIceStandaloneOn

      ! Add frazil and interface melt/freeze to get total fresh water flux
      if ( associated(frazilIceFreshwaterFlux) ) then
          do iCell = 1, nCells
                landIceFreshwaterFluxTotal(iCell) = landIceFreshwaterFlux(iCell) + landIceFloatingMask(iCell)*frazilIceFreshwaterFlux(iCell)
          end do
      else
          do iCell = 1, nCells
                landIceFreshwaterFluxTotal(iCell) = landIceFreshwaterFlux(iCell)
          end do
      end if

      call mpas_timer_stop("land_ice_build_arrays")

   !--------------------------------------------------------------------

   end subroutine ocn_surface_land_ice_fluxes_build_arrays!}}}


!***********************************************************************
!
!  routine ocn_surface_land_ice_fluxes_accumulate_fluxes
!
!> \brief Accumulate total land-ice mass and heat fluxes
!> \author Xylar Asay-Davis
!> \date   09/02/2023
!> \details
!>  This routine accumulates land-ice mass and heat fluxes into variables
!>  used to keep track of the total mass and energy budgets.
!-----------------------------------------------------------------------

   subroutine ocn_surface_land_ice_fluxes_accumulate_fluxes(meshPool, &
      forcingPool, statePool, dt, err) !{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool        !< Input: mesh information
      real(kind=RKIND), intent(in) :: dt ! the time step over which to accumulate fluxes

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(inout) :: &
         forcingPool, & !< Input: Forcing information
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

      integer :: iCell, nCells
      integer, dimension(:), pointer :: nCellsArray

      real (kind=RKIND), dimension(:), pointer :: landIceFreshwaterFlux, &
                                                  heatFluxToLandIce
      real (kind=RKIND), dimension(:), pointer :: accumulatedLandIceMassOld, &
                                                  accumulatedLandIceMassNew, &
                                                  accumulatedLandIceHeatOld, &
                                                  accumulatedLandIceHeatNew

      err = 0

      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      nCells = nCellsArray( size(nCellsArray) )

      if (landIceStandaloneOn .or. landIceDataOn) then
        call mpas_pool_get_array(statePool, 'accumulatedLandIceMass', accumulatedLandIceMassNew, 2)
        call mpas_pool_get_array(statePool, 'accumulatedLandIceMass', accumulatedLandIceMassOld, 1)
        call mpas_pool_get_array(forcingPool, 'landIceFreshwaterFlux', landIceFreshwaterFlux)
        ! accumulate land-ice mass
        do iCell = 1, nCells
          accumulatedLandIceMassNew(iCell) = accumulatedLandIceMassOld(iCell) - dt*landIceFreshwaterFlux(iCell)
        end do
      end if

      if (landIceStandaloneOn) then
        call mpas_pool_get_array(forcingPool, 'heatFluxToLandIce', heatFluxToLandIce)
        call mpas_pool_get_array(statePool, 'accumulatedLandIceHeat', accumulatedLandIceHeatNew, 2)
        call mpas_pool_get_array(statePool, 'accumulatedLandIceHeat', accumulatedLandIceHeatOld, 1)
        ! accumulate land-ice heat
        do iCell = 1, nCells
          accumulatedLandIceHeatNew(iCell) = accumulatedLandIceHeatOld(iCell) + dt*heatFluxToLandIce(iCell)
        end do
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_surface_land_ice_fluxes_accumulate_fluxes!}}}

!***********************************************************************
!
!  routine ocn_surface_land_ice_fluxes_init
!
!> \brief Initializes land-ice forcing
!> \author Xylar Asay-Davis
!> \date   29 Jan 2023
!> \details
!>  This routine initializes a variety of quantities related to
!>  land-ice forcing.
!
!-----------------------------------------------------------------------

   subroutine ocn_surface_land_ice_fluxes_init(err)!{{{

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< [out] error flag

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !*** Set error code and default values for module variables

      err = 0

      isomipOn         = .false.
      jenkinsOn        = .false.
      hollandJenkinsOn = .false.
      cpLandIce        = 0.0_RKIND
      rhoLandIce       = 0.0_RKIND
      ISOMIPgammaT     = 0.0_RKIND

      !*** Determine whether land ice fluxes are on

      select case (trim(config_land_ice_flux_mode))
      case ('data','Data','DATA')
         landIcePressureOn = .true.
         landIceFluxesOn = .true.
         landIceDataOn = .true.
         landIceStandaloneOn = .false.
         landIceCoupledOn = .false.
      case ('standalone','Standalone','STANDALONE')
         landIcePressureOn = .true.
         landIceFluxesOn = .true.
         landIceDataOn = .false.
         landIceStandaloneOn = .true.
         landIceCoupledOn = .false.
      case ('coupled','Coupled','COUPLED')
         landIcePressureOn = .true.
         landIceFluxesOn = .true.
         landIceDataOn = .false.
         landIceStandaloneOn = .false.
         landIceCoupledOn = .true.
      case ('pressure_only','Pressure_Only','PRESSURE_ONLY')
         landIcePressureOn = .true.
         landIceFluxesOn = .false.
         landIceDataOn = .false.
         landIceStandaloneOn = .false.
         landIceCoupledOn = .false.
      case ('off','Off','OFF')
         landIcePressureOn = .false.
         landIceFluxesOn = .false.
         landIceDataOn = .false.
         landIceStandaloneOn = .false.
         landIceCoupledOn = .false.
      case default
         call mpas_log_write( ' Error: Incorrect choice of config_land_ice_flux_mode: '// config_land_ice_flux_mode)
         err = 1
         return
      end select

      !*** Determine which flux formulation to use

      select case (trim(config_land_ice_flux_formulation))

      case ('isomip','Isomip','ISOMIP')
         isomipOn = .true.
      case ('jenkins', 'Jenkins','JENKINS')
         jenkinsOn = .true.
      case ('hollandjenkins','HollandJenkins','HOLLANDJENKINS')
         hollandJenkinsOn = .true.
      case default
         call mpas_log_write( &
            "config_land_ice_flux_formulation not one of 'ISOMIP', 'Jenkins', " &
               // "or 'HollandJenkins'.", &
               MPAS_LOG_CRIT)
         err = 1
      end select

      useHollandJenkinsAdvDiff = &
                       config_land_ice_flux_useHollandJenkinsAdvDiff

      !*** Set physical values

      cpLandIce    = config_land_ice_flux_cp_ice
      rhoLandIce   = config_land_ice_flux_rho_ice
      ISOMIPgammaT = config_land_ice_flux_ISOMIP_gammaT

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
      interfacePressure            !< Input: pressure at the ice-ocean interface

    integer, intent(in) :: nCells !< Input: number of cells in each array

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

    real (kind=RKIND) :: T0, transferVelocityRatio, Tlatent, a, b, c, &
                         dTf_dS
    integer :: iCell

    real (kind=RKIND), parameter :: minInterfaceSalinity = 0.001_RKIND

    err = 0
    Tlatent = latent_heat_fusion_mks/cp_sw

    !$omp parallel
    !$omp do schedule(runtime) private(T0, dTf_dS, transferVelocityRatio, a, b, c)
    do iCell = 1, nCells
      if (mask(iCell) == 0) cycle

      T0 = ocn_freezing_temperature(salinity=0.0_RKIND, pressure=interfacePressure(iCell), &
                                    inLandIceCavity=.true.)
      dTf_dS = ocn_freezing_temperature_salinity_deriv(salinity=0.0_RKIND, pressure=interfacePressure(iCell), &
                                                       inLandIceCavity=.true.)

      transferVelocityRatio = oceanSaltTransferVelocity(iCell)/oceanHeatTransferVelocity(iCell)

      a = -dTf_dS
      b = transferVelocityRatio*Tlatent + oceanTemperature(iCell) - T0
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
      outIceHeatFlux(iCell) = -cpLandIce*outFreshwaterFlux(iCell)*outInterfaceTemperature(iCell)
    end do
    !$omp end do
    !$omp end parallel

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
    cpRatio = cpLandIce/cp_sw
    !$omp parallel
    !$omp do schedule(runtime) &
    !$omp private(T0, dTf_dS, transferVelocityRatio, Tlatent, eta, TlatentStar, a, b, c, err)
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
         call mpas_log_write('interfaceSalinity is negative.', MPAS_LOG_CRIT)
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
      ! the ice loses heat simply by the loss of ice mass at the data
      ! (surface?) ice temperature
      outIceHeatFlux(iCell) = -cpLandIce*outFreshwaterFlux(iCell)*iceTemperature(iCell)
    end do
    !$omp end do
    !$omp end parallel

  !--------------------------------------------------------------------

  end subroutine compute_HJ99_melt_fluxes !}}}


!***********************************************************************

end module ocn_surface_land_ice_fluxes

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
