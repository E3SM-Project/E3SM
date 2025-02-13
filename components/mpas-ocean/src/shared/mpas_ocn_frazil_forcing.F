! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_frazil_forcing
!
!> \brief MPAS ocean frazil formation module
!> \author Todd Ringler
!> \date   10/19/2015
!> \details
!>  This module contains routines for the formation of frazil ice.
!
!-----------------------------------------------------------------------

module ocn_frazil_forcing

   use mpas_kind_types
   use mpas_constants
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timekeeping
   use mpas_timer
   use ocn_constants
   use ocn_config
   use ocn_mesh
   use ocn_diagnostics_variables
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

   public :: ocn_frazil_forcing_build_arrays, &
             ocn_frazil_forcing_tracers, &
             ocn_frazil_forcing_layer_thickness, &
             ocn_frazil_forcing_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: frazilFormationOn
   integer :: usePorosity   ! 1 when computing internal frazil salinity using porosity
                            ! 0 when using config_frazil_sea_ice_reference_salinity

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_frazil_forcing_tracers
!
!> \brief   Determines the tracer tendency due to frazil
!> \author  Todd Ringler
!> \date    18 October 2015
!> \details
!>  This routine adds to the tracer tendency arrays
!>  used to compute tracer at n+1.
!
!-----------------------------------------------------------------------

   subroutine ocn_frazil_forcing_tracers(meshPool, tracersPool, groupName, forcingPool, tracersTend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      character (len=*) :: groupName !< Input: Name of tracer group

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: meshPool !< Input/Output: mesh information
      type (mpas_pool_type), intent(inout) :: tracersPool !< Input/Output: tracer tendency pool
      type (mpas_pool_type), intent(inout) :: forcingPool  !< Input/Output: forcing pool holding frazil-induced tendencies
      real (kind=RKIND), dimension(:,:,:), intent(inout) ::  tracersTend
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

      if ( .not. frazilFormationOn ) return

      call mpas_timer_start("tracer frazil")

      if ( trim(groupName) == 'activeTracers' ) then
         call ocn_frazil_forcing_active_tracers(meshPool, tracersPool, forcingPool, tracersTend, err)
      end if

      call mpas_timer_stop("tracer frazil")

   end subroutine ocn_frazil_forcing_tracers!}}}

!***********************************************************************
!
!  routine ocn_frazil_forcing_layer_thickness
!
!> \brief   Add tendency due to frazil processes
!> \author  Todd Ringler
!> \date    18 October 2015
!> \details
!>  This routine adds a tendency to layer thickness due to frazil formation
!
!-----------------------------------------------------------------------

   subroutine ocn_frazil_forcing_layer_thickness(forcingPool, layerThicknessTend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: forcingPool !< Input: Forcing information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      real (kind=RKIND), intent(inout), dimension(:,:) :: layerThicknessTend

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

      integer :: iCell, k, nCells
      real (kind=RKIND), dimension(:,:), pointer :: frazilLayerThicknessTendency

      err = 0

      if ( .not. frazilFormationOn ) return

      call mpas_timer_start("frazil thickness tendency")

      call mpas_pool_get_array(forcingPool, 'frazilLayerThicknessTendency', frazilLayerThicknessTendency)

      ! frazil fields are needed only over 0 and 1 halos
      nCells = nCellsHalo( 1 )

      ! Build surface fluxes at cell centers
#ifdef MPAS_OPENACC
      !$acc enter data copyin(frazilLayerThicknessTendency)

      !$acc parallel loop &
      !$acc     present(layerThicknessTend, minLevelCell, maxLevelCell, frazilLayerThicknessTendency) &
      !$acc     private(k)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(k)
#endif
      do iCell = 1, nCells
        do k = minLevelCell(iCell), maxLevelCell(iCell)
          layerThicknessTend(k,iCell) = layerThicknessTend(k,iCell) + frazilLayerThicknessTendency(k,iCell)
        end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

#ifdef MPAS_OPENACC
      !$acc exit data delete(frazilLayerThicknessTendency)
#endif

      call mpas_timer_stop("frazil thickness tendency")

   end subroutine ocn_frazil_forcing_layer_thickness!}}}


!***********************************************************************
!
!  routine ocn_frazil_forcing_active_tracers
!
!> \brief   Adds the active tracers forcing due to frazil
!> \author  Todd Ringler
!> \date    18 October 2015
!> \details
!>  This routine adds the active tracers forcing due to frazil
!>  from which tracer tendencies are computed later.
!
!-----------------------------------------------------------------------

   subroutine ocn_frazil_forcing_active_tracers(meshPool, tracersPool, forcingPool, activeTracersTend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: meshPool !< Input: mesh information
      type (mpas_pool_type), intent(inout) :: tracersPool !< Input: tracer tendency pool
      type (mpas_pool_type), intent(inout) :: forcingPool !< Input: Forcing information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      real (kind=RKIND), dimension(:,:,:), intent(inout) :: activeTracersTend

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

      integer :: iCell, k, nCells
      integer, pointer :: indexTemperaturePtr
      integer, pointer :: indexSalinityPtr
      integer :: indexTemperature
      integer :: indexSalinity
      integer, pointer, dimension(:) :: nCellsArray
      integer, pointer, dimension(:) :: minLevelCell, maxLevelCell

      real (kind=RKIND), dimension(:,:), pointer :: frazilTemperatureTendency
      real (kind=RKIND), dimension(:,:), pointer :: frazilSalinityTendency

      err = 0

      if ( .not. frazilFormationOn ) return

      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_dimension(tracersPool, 'index_temperature', &
                                                 indexTemperaturePtr)
      call mpas_pool_get_dimension(tracersPool, 'index_salinity', &
                                                 indexSalinityPtr)
      indexTemperature = indexTemperaturePtr
      indexSalinity    = indexSalinityPtr

      call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(forcingPool, 'frazilTemperatureTendency', frazilTemperatureTendency)
      call mpas_pool_get_array(forcingPool, 'frazilSalinityTendency', frazilSalinityTendency)

      ! frazil fields are needed only over 0 and 1 halos
      nCells = nCellsArray( 2 )

      ! add to surface fluxes at cell centers
      !$omp parallel
      !$omp do schedule(runtime) private(k)
      do iCell = 1, nCells
        do k = minLevelCell(iCell), maxLevelCell(iCell)
          activeTracersTend(indexTemperature,k,iCell) = activeTracersTend(indexTemperature,k,iCell) &
                                                      + frazilTemperatureTendency(k,iCell)
          activeTracersTend(indexSalinity,k,iCell) = activeTracersTend(indexSalinity,k,iCell) + frazilSalinityTendency(k,iCell)
        end do
      end do
      !$omp end do
      !$omp end parallel

   end subroutine ocn_frazil_forcing_active_tracers!}}}

!***********************************************************************
!
!  routine ocn_frazil_forcing_build_arrays
!
!> \brief   Performs the formation of frazil within the ocean.
!> \author  Todd Ringler
!> \date    10/19/2015
!> \details
!>   ocn_frazil_forcing_build_arrays computes the tendencies to layer thickness, temperature and salinity
!>   due to the creation and possible melting of frazil ice
!>
!>   these tendencies can be retrieved at any point by calling into ocn_frazil_forcing_{tracers, thickness} routines
!>
!>   the pressure exerted by the frazil on the ocean "surface" is added to the pressure computation in diagnostics
!>
!>   this routine should be call at the beginning of whatever time stepping method is utilized
!>      and the tendencies should be retieved when building up the RHS of the thickess, temperature
!>      and salinity equations.
!>
!>   this routine is only applicable to the surface pressure, thickness and active tracer fields
!
!-----------------------------------------------------------------------

   subroutine ocn_frazil_forcing_build_arrays(domain, meshPool, forcingPool, statePool, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), pointer, intent(in) :: meshPool !< Input: Mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: domain
      type (mpas_pool_type), pointer, intent(inout) :: statePool !< Input: State information
      type (mpas_pool_type), pointer, intent(inout) :: forcingPool !< Input: Forcing information
      integer, intent(inout) :: err !< Error flag

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), pointer :: tracersPool

      real (kind=RKIND), dimension(:,:), pointer :: frazilLayerThicknessTendency
      real (kind=RKIND), dimension(:,:), pointer :: frazilTemperatureTendency
      real (kind=RKIND), dimension(:,:), pointer :: frazilSalinityTendency
      real (kind=RKIND), dimension(:), pointer :: frazilSurfacePressure
      real (kind=RKIND), dimension(:), pointer :: lowSalinityFrazilIce
      integer, dimension(:), pointer :: landIceMask, landIceFloatingMask

      integer :: iCell, k, nCells
      integer, dimension(:), pointer :: nCellsArray
      integer, pointer :: nVertLevels
      real (kind=RKIND) :: dt, columnTemperatureMin

      type (MPAS_timeInterval_type) :: timeStep

      real (kind=RKIND) :: newFrazilIceThickness, newThicknessWeightedSaltContent
      real (kind=RKIND) :: sumNewFrazilIceThickness, sumNewThicknessWeightedSaltContent
      real (kind=RKIND) :: meltedFrazilIceThickness, meltedThicknessWeightedSaltContent
      real (kind=RKIND) :: oceanFreezingTemperature, frazilSalinity

      real (kind=RKIND), pointer, dimension(:)     :: accumulatedFrazilIceMassNew
      real (kind=RKIND), pointer, dimension(:)     :: accumulatedFrazilIceMassOld
      real (kind=RKIND), pointer, dimension(:)     :: accumulatedFrazilIceSalinityNew
      real (kind=RKIND), pointer, dimension(:)     :: accumulatedFrazilIceSalinityOld
      real (kind=RKIND), pointer, dimension(:)     :: accumulatedLandIceFrazilMassNew
      real (kind=RKIND), pointer, dimension(:)     :: accumulatedLandIceFrazilMassOld
      real (kind=RKIND), pointer, dimension(:)     :: frazilIceFreshwaterFlux
      real (kind=RKIND), pointer, dimension(:)     :: ssh
      real (kind=RKIND), pointer, dimension(:,:)   :: layerThickness
      real (kind=RKIND), pointer, dimension(:,:,:) :: activeTracers

      integer, dimension(:), pointer :: minLevelCell, maxLevelCell
      integer, pointer :: indexTemperaturePtr !< index in tracers array for temperature
      integer, pointer :: indexSalinityPtr !< index in tracers array for salinity
      integer :: indexTemperature !< index in tracers array for temperature
      integer :: indexSalinity !< index in tracers array for salinity
      integer :: kBottomFrazil ! k index where testing for frazil begins

      logical :: underFloatingLandIce, underLandIce ! indicates if we are under land ice

      real (kind=RKIND) :: potential  ! scalar holding freezing/melt potential
      real (kind=RKIND) :: freezingEnergy  ! energy available for freezing, positive definite
      real (kind=RKIND) :: meltingEnergy ! energy available for melting, positive definite

      ! if frazil is not enabled, return
      if(.not. frazilFormationOn) return

      call mpas_timer_start("frazil")

      ! get pool pointers
      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      ! get dimensions
      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(tracersPool, 'index_temperature', &
                                                 indexTemperaturePtr)
      call mpas_pool_get_dimension(tracersPool, 'index_salinity', &
                                                 indexSalinityPtr)
      indexTemperature = indexTemperaturePtr
      indexSalinity    = indexSalinityPtr

      ! get mesh information
      call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

      ! get arrays
      ! note: state information is used to produce tendencies, so always grab "new" time level
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
      call mpas_pool_get_array(statePool, 'accumulatedFrazilIceMass', accumulatedFrazilIceMassNew, 2)
      call mpas_pool_get_array(statePool, 'accumulatedFrazilIceMass', accumulatedFrazilIceMassOld, 1)
      call mpas_pool_get_array(statePool, 'accumulatedFrazilIceSalinity', accumulatedFrazilIceSalinityNew, 2)
      call mpas_pool_get_array(statePool, 'accumulatedFrazilIceSalinity', accumulatedFrazilIceSalinityOld, 1)
      call mpas_pool_get_array(statePool, 'accumulatedLandIceFrazilMass', accumulatedLandIceFrazilMassNew, 2)
      call mpas_pool_get_array(statePool, 'accumulatedLandIceFrazilMass', accumulatedLandIceFrazilMassOld, 1)
      call mpas_pool_get_array(statePool, 'ssh', ssh, 1)
      call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
      call mpas_pool_get_array(forcingPool, 'frazilTemperatureTendency', frazilTemperatureTendency)
      call mpas_pool_get_array(forcingPool, 'frazilSalinityTendency', frazilSalinityTendency)
      call mpas_pool_get_array(forcingPool, 'frazilLayerThicknessTendency', frazilLayerThicknessTendency)
      call mpas_pool_get_array(forcingPool, 'frazilSurfacePressure', frazilSurfacePressure)
      call mpas_pool_get_array(forcingPool, 'lowSalinityFrazilIce', lowSalinityFrazilIce)
      call mpas_pool_get_array(forcingPool, 'landIceMask', landIceMask)
      call mpas_pool_get_array(forcingPool, 'landIceFloatingMask', landIceFloatingMask)
      call mpas_pool_get_array(forcingPool, 'frazilIceFreshwaterFlux', frazilIceFreshwaterFlux)

      ! get time step in units of seconds
      timeStep = mpas_get_clock_timestep(domain % clock, ierr=err)
      call mpas_get_timeInterval(timeStep, dt=dt)

      ! initialize frazil tendency fields
      frazilTemperatureTendency = 0.0_RKIND
      frazilSalinityTendency = 0.0_RKIND
      frazilLayerThicknessTendency = 0.0_RKIND

      ! initialize frazil salinity conservation field
      lowSalinityFrazilIce = 0.0_RKIND

      ! frazil fields are needed only over 0 and 1 halos
      nCells = nCellsArray( 2 )

      ! loop over all columns
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp private(underFloatingLandIce, underLandIce, kBottomFrazil, columnTemperatureMin, k, sumNewFrazilIceThickness, &
      !$omp         sumNewThicknessWeightedSaltContent, oceanFreezingTemperature, potential, &
      !$omp         freezingEnergy, meltingEnergy, frazilSalinity, newFrazilIceThickness, &
      !$omp         newThicknessWeightedSaltContent, meltedFrazilIceThickness, &
      !$omp         meltedThicknessWeightedSaltContent)
      do iCell=1,nCells

        underLandIce = .false.
        if ( associated(landIceMask) ) then
           underLandIce = (landIceMask(iCell) == 1)
        end if

        underFloatingLandIce = .false.
        if ( associated(landIceFloatingMask) ) then
           underFloatingLandIce = (landIceFloatingMask(iCell) == 1)
        end if

        if (underFloatingLandIce .and. .not. config_frazil_under_land_ice) then
          ! skip this cell because we're not doing frazil under floating land ice
          cycle
        end if

        if (.not. underLandIce .and. .not. config_frazil_in_open_ocean) then
          ! skip this cell because we're not doing frazil in the open ocean
          cycle
        end if

        if (underLandIce .and. .not. underFloatingLandIce) then
          ! skip this cell because we're not doing frazil under grounded land ice
          cycle
        end if

        ! find deepest level where frazil can be created
        kBottomFrazil=maxLevelCell(iCell)
        do k=maxLevelCell(iCell), minLevelCell(iCell), -1
         ! add the ssh so frazil can form below land ice, where the ssh is depressed
          if (-zMid(k,iCell) < -ssh(iCell) + config_frazil_maximum_depth) then
            kBottomFrazil=k
            exit
          endif
        enddo

        ! find minimum temperature between 1:kBottomFrazil
        columnTemperatureMin = 1.0e30_RKIND
        do k=minLevelCell(iCell),kBottomFrazil
          if ( activeTracers(indexTemperature,k,iCell) < columnTemperatureMin ) then
             columnTemperatureMin = activeTracers(indexTemperature,k,iCell)
          end if
        enddo

        ! test min temperature agains max freezing temperature to see if we should even consider creating frazil
        if (columnTemperatureMin > config_frazil_maximum_freezing_temperature) cycle

        ! initialize the sum of new frazil ice created
        sumNewFrazilIceThickness = 0.0_RKIND
        sumNewThicknessWeightedSaltContent = 0.0_RKIND

        ! loop from maximum depth of frazil creation to surface
        do k = kBottomFrazil, minLevelCell(iCell), -1

          ! get freezing temperature
          oceanFreezingTemperature = ocn_freezing_temperature(salinity=activeTracers(indexSalinity,k,iCell), &
                                                              pressure=pressure(k,iCell), &
                                                              inLandIceCavity=underLandIce)

          potential = layerThickness(k,iCell) * config_specific_heat_sea_water &
                    * rho_sw * (activeTracers(indexTemperature,k,iCell) - oceanFreezingTemperature)

          freezingEnergy = max(0.0_RKIND, -potential)
          meltingEnergy = max(0.0_RKIND, potential)

          if (freezingEnergy > 0) then

            ! get frazil salinity
            if (underLandIce) then
              frazilSalinity = usePorosity * config_frazil_ice_porosity * activeTracers(indexSalinity, k, iCell)
            else
              frazilSalinity = usePorosity * config_frazil_ice_porosity * activeTracers(indexSalinity, k, iCell) &
                 + config_frazil_sea_ice_reference_salinity * (1 - usePorosity)
            end if
            frazilSalinity = min( frazilSalinity, activeTracers(indexSalinity, k, iCell) )

            ! new frazil ice formation measured in meters
            newFrazilIceThickness = freezingEnergy / (config_frazil_heat_of_fusion * config_frazil_ice_density)

            ! limit the frazil formed appropriately
            newFrazilIceThickness = min(newFrazilIceThickness, layerThickness(k,iCell) * config_frazil_fractional_thickness_limit)

            ! Determine new salt in frazil
            newThicknessWeightedSaltContent = newFrazilIceThickness * frazilSalinity

            ! compute tendency to thickness, temperature and salinity
            ! layerTendency is scaled so that mass of ice created == mass of ocean water removed

            ! layer thickness decreased due to creation of frazil
            ! note: -- this has to be density (not rho_sw) to keep buoyancy equal
            frazilLayerThicknessTendency(k,iCell) = - newFrazilIceThickness * config_frazil_ice_density / density(k,iCell) / dt

            ! salt is extracted with the frazil
            frazilSalinityTendency(k,iCell) = - newThicknessWeightedSaltContent / dt

            ! ocean fluid temperature is warmed due to creation of frazil
            frazilTemperatureTendency(k,iCell) =  + ( newFrazilIceThickness * config_frazil_heat_of_fusion &
                                               * config_frazil_ice_density ) / (config_specific_heat_sea_water * rho_sw) / dt

            ! keep track of sum of frazil ice
            sumNewFrazilIceThickness = sumNewFrazilIceThickness + newFrazilIceThickness
            sumNewThicknessWeightedSaltContent = sumNewThicknessWeightedSaltContent + newThicknessWeightedSaltContent

          else
            ! ocean water is warm enough to melt frazil

            ! test to see if there is frazil to be melted
            if (sumNewFrazilIceThickness > 0.0_RKIND) then

              ! Frazil melting
              meltedFrazilIceThickness = meltingEnergy / (config_frazil_heat_of_fusion * config_frazil_ice_density)

              ! limit melting by what there is to melt
              meltedFrazilIceThickness = min(meltedFrazilIceThickness, sumNewFrazilIceThickness)

              ! limit melting by fraction of layer thickness
              meltedFrazilIceThickness = min(meltedFrazilIceThickness, layerThickness(k,iCell) &
                                       * config_frazil_fractional_thickness_limit)

              ! assign some salt to the melted ice
              meltedThicknessWeightedSaltContent = meltedFrazilIceThickness &
                                                 * ( sumNewThicknessWeightedSaltContent / sumNewFrazilIceThickness )

              ! compute tendency to thickness, temperature and salinity

              ! layer thickness increases due to melting of frazil
              ! note -- scaling by local ocean density to mimimize surface pressure forcing errors
              frazilLayerThicknessTendency(k,iCell) = + meltedFrazilIceThickness * config_frazil_ice_density &
                                                    / density(k,iCell) / dt

              ! salt is released into ocean with the melting frazil
              frazilSalinityTendency(k,iCell) = + meltedThicknessWeightedSaltContent/ dt

              ! ocean fluid temperature is cooled due to melting of frazil
              frazilTemperatureTendency(k,iCell) = - ( meltedFrazilIceThickness * config_frazil_heat_of_fusion &
                                                 * config_frazil_ice_density ) / (config_specific_heat_sea_water * rho_sw) &
                                                 / dt

              ! ocean fluid temperature is modified due to mixing of existing water with meltwater at the freezing point
              frazilTemperatureTendency(k,iCell) = frazilTemperatureTendency(k,iCell) + &
                   meltedFrazilIceThickness * oceanFreezingTemperature / dt

              ! keep track of new frazil ice
              sumNewThicknessWeightedSaltContent = max(0.0_RKIND,sumNewThicknessWeightedSaltContent - meltedThicknessWeightedSaltContent)
              sumNewFrazilIceThickness = max(0.0_RKIND,sumNewFrazilIceThickness - meltedFrazilIceThickness)

            endif  ! if (sumNewFrazilIceThickness > 0.0_RKIND)

          endif  ! if (freezingEnergy < 0)

        enddo   ! do k=kBottom,minLevelCell,-1

        ! accumulate frazil mass to column total
        ! note: the accumulatedFrazilIceMass (at both time levels) is reset to zero after being sent to the coupler
        accumulatedFrazilIceMassNew(iCell) = accumulatedFrazilIceMassOld(iCell) + sumNewFrazilIceThickness &
           * config_frazil_ice_density

        accumulatedFrazilIceSalinityNew(iCell) = accumulatedFrazilIceSalinityOld(iCell) + &
           (1-usePorosity) * sumNewThicknessWeightedSaltContent

        !NJ: Assume frazil salinity is the coupling salinity
        if ( underLandIce ) then
           ! frazilIceFreshwaterFlux forms part of the landIceFreshwaterFluxTotal computed in surface land ice fluxes
           frazilIceFreshwaterFlux(iCell) = - (sumNewFrazilIceThickness * config_frazil_ice_density) / dt
           ! accumulate frazil formed under land ice in case we're not coupling and we need to keep track of it
           ! for freshwater budgets
           accumulatedLandIceFrazilMassNew(iCell) = accumulatedLandIceFrazilMassOld(iCell) &
                                                  + sumNewFrazilIceThickness * config_frazil_ice_density
           ! There is no accumulatedLandIceFrazilSalinity because it is always 0
           frazilSalinity = 0.0_RKIND
        else
           frazilSalinity = config_frazil_sea_ice_reference_salinity
        end if

        !Determine salt expelled in the upper layer to arrive at the assumed frazil salt for coupling
        if (sumNewFrazilIceThickness > 0.0_RKIND .and. usePorosity .ge. 1) then
           newThicknessWeightedSaltContent = sumNewFrazilIceThickness * frazilSalinity
           if (newThicknessWeightedSaltContent > sumNewThicknessWeightedSaltContent) then
            !Warning: salt conservation error. Track with lowSalinityFrazilIce
            lowSalinityFrazilIce(iCell) = newThicknessWeightedSaltContent - sumNewThicknessWeightedSaltContent
            accumulatedFrazilIceSalinityNew(iCell) = accumulatedFrazilIceSalinityOld(iCell) + max(0.0_RKIND, sumNewThicknessWeightedSaltContent)
           else
            !salt expelled at the surface
            frazilSalinityTendency(minLevelCell(iCell),iCell) = frazilSalinityTendency(minLevelCell(iCell),iCell) + &
                max(0.0_RKIND,(sumNewThicknessWeightedSaltContent -  newThicknessWeightedSaltContent) ) / dt
            accumulatedFrazilIceSalinityNew(iCell) = accumulatedFrazilIceSalinityOld(iCell) + newThicknessWeightedSaltContent
           endif
         endif ! if frazil ice is formed
      enddo   ! do iCell = 1, nCells
      !$omp end do
      !$omp end parallel

      if ( config_frazil_use_surface_pressure ) then
        !$omp parallel
        !$omp do schedule(runtime)
        do iCell = 1, nCells
           frazilSurfacePressure(iCell) = accumulatedFrazilIceMassNew(iCell) * gravity
        end do
        !$omp end do
        !$omp end parallel
      end if

      call mpas_timer_stop("frazil")

   end subroutine ocn_frazil_forcing_build_arrays!}}}

!***********************************************************************
!
!  routine ocn_frazil_forcing_init
!
!> \brief   Initializes ocean frazil ice module.
!> \author  Todd Ringler
!> \date    10/19/2015
!> \details
!>  This routine initializes the ocean frazil ice module and variables..
!
!-----------------------------------------------------------------------

   subroutine ocn_frazil_forcing_init(err)!{{{

      integer, intent(out) :: err !< Output: error flag

      err = 0

      frazilFormationOn = .false.

      if(config_use_frazil_ice_formation) then
        frazilFormationOn = .true.
      end if

      usePorosity = 0
      if(config_use_frazil_ice_porosity) usePorosity = 1

   end subroutine ocn_frazil_forcing_init!}}}

!***********************************************************************

end module ocn_frazil_forcing

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
