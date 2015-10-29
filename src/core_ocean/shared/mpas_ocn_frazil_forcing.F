! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
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
   type (timer_node), pointer :: timer_frazil

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
      type (mpas_pool_type), intent(inout) :: meshPool !< Input/Output: mesh information
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

      if ( trim(groupName) == 'activeTracers' ) then
         call ocn_frazil_forcing_active_tracers(meshPool, tracersPool, forcingPool, tracersTend, err)
      end if

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

   subroutine ocn_frazil_forcing_layer_thickness(meshPool, forcingPool, layerThicknessTend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: meshPool !< Input: mesh information
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

      integer :: iCell, k
      integer, pointer :: nCells
      integer, dimension(:), pointer :: maxLevelCell
      real (kind=RKIND), dimension(:,:), pointer :: frazilLayerThicknessTendency

      err = 0

      if ( .not. frazilFormationOn ) return

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(forcingPool, 'frazilLayerThicknessTendency', frazilLayerThicknessTendency)

      ! Build surface fluxes at cell centers
      !$omp do schedule(runtime) private(k)
      do iCell = 1, nCells
        do k = 1, maxLevelCell(iCell)
          layerThicknessTend(k,iCell) = layerThicknessTend(k,iCell) + frazilLayerThicknessTendency(k,iCell)
        end do
      end do
      !$omp end do

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

      integer :: iCell, k
      integer, pointer :: nCells
      integer, pointer :: indexTemperature
      integer, pointer :: indexSalinity
      integer, pointer, dimension(:) :: maxLevelCell

      real (kind=RKIND), dimension(:,:), pointer :: frazilTemperatureTendency
      real (kind=RKIND), dimension(:,:), pointer :: frazilSalinityTendency

      err = 0

      if ( .not. frazilFormationOn ) return

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(tracersPool, 'index_temperature', indexTemperature)
      call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexSalinity)

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(forcingPool, 'frazilTemperatureTendency', frazilTemperatureTendency)
      call mpas_pool_get_array(forcingPool, 'frazilSalinityTendency', frazilSalinityTendency)

      ! add to surface fluxes at cell centers
      !$omp do schedule(runtime) private(k)
      do iCell = 1, nCells
        do k = 1, maxLevelCell(iCell)
          activeTracersTend(indexTemperature,k,iCell) = activeTracersTend(indexTemperature,k,iCell) + frazilTemperatureTendency(k,iCell)
          activeTracersTend(indexSalinity,k,iCell) = activeTracersTend(indexSalinity,k,iCell) + frazilSalinityTendency(k,iCell)
        end do
      end do
      !$omp end do

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
!>	due to the creation and possible melting of frazil ice
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

   subroutine ocn_frazil_forcing_build_arrays(domain, meshPool, forcingPool, diagnosticsPool, statePool, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), pointer, intent(in) :: meshPool !< Input: Mesh information
      type (mpas_pool_type), pointer, intent(in) :: diagnosticsPool !< Input: Diagnostic information

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
      real (kind=RKIND),   dimension(:), pointer :: frazilSurfacePressure

      integer :: iCell, k
      integer, pointer :: nCells, nVertLevels
      real (kind=RKIND) :: dt, columnTemperatureMin

      type (MPAS_timeInterval_type) :: timeStep
      real (kind=RKIND), pointer :: config_frazil_heat_of_fusion
      real (kind=RKIND), pointer :: config_frazil_sea_ice_density
      real (kind=RKIND), pointer :: config_frazil_fractional_thickness_limit
      real (kind=RKIND), pointer :: config_frazil_maximum_depth
      real (kind=RKIND), pointer :: config_specific_heat_sea_water 
      real (kind=RKIND), pointer :: config_frazil_ice_reference_salinity
      real (kind=RKIND), pointer :: config_frazil_maximum_freezing_temperature

      real (kind=RKIND) :: newFrazilIceThickness
      real (kind=RKIND) :: sumNewFrazilIceThickness
      real (kind=RKIND) :: meltedFrazilIceThickness
      real (kind=RKIND) :: oceanFreezingTemperature

      real (kind=RKIND), pointer, dimension(:)     :: accumulatedFrazilIceMassNew
      real (kind=RKIND), pointer, dimension(:)     :: accumulatedFrazilIceMassOld
      real (kind=RKIND), pointer, dimension(:,:)   :: zMid
      real (kind=RKIND), pointer, dimension(:,:)   :: layerThickness
      real (kind=RKIND), pointer, dimension(:,:)   :: density
      real (kind=RKIND), pointer, dimension(:,:,:) :: activeTracers

      integer, dimension(:), pointer :: maxLevelCell
      integer, pointer :: indexTemperature !< index in tracers array for temperature
      integer, pointer :: indexSalinity !< index in tracers array for salinity
      integer :: kBottomFrazil  ! k index where testing for frazil begins

      real (kind=RKIND) :: potential  ! scalar holding freezing/melt potential
      real (kind=RKIND) :: freezingEnergy  ! energy available for freezing, positive definite
      real (kind=RKIND) :: meltingEnergy ! energy available for melting, positive definite

      ! if frazil is not enabled, return
      if(.not. frazilFormationOn) return

      call mpas_timer_start("frazil", .false., timer_frazil)

      ! get pool pointers
      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      ! get dimensions
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(tracersPool, 'index_temperature', indexTemperature)
      call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexSalinity)

      ! get mesh information
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

      ! get arrays
      ! note: state information is used to produce tendencies, so always grab "new" time level
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
      call mpas_pool_get_array(statePool, 'accumulatedFrazilIceMass', accumulatedFrazilIceMassNew, 2)
      call mpas_pool_get_array(statePool, 'accumulatedFrazilIceMass', accumulatedFrazilIceMassOld, 1)
      call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
      call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)
      call mpas_pool_get_array(diagnosticsPool, 'density', density)
      call mpas_pool_get_array(forcingPool, 'frazilTemperatureTendency', frazilTemperatureTendency)
      call mpas_pool_get_array(forcingPool,'frazilSalinityTendency', frazilSalinityTendency)
      call mpas_pool_get_array(forcingPool,'frazilLayerThicknessTendency', frazilLayerThicknessTendency)
      call mpas_pool_get_array(forcingPool,'frazilSurfacePressure', frazilSurfacePressure)

      ! get configure parameters
      call mpas_pool_get_config(ocnConfigs, 'config_frazil_maximum_depth', config_frazil_maximum_depth)
      call mpas_pool_get_config(ocnConfigs, 'config_frazil_fractional_thickness_limit', config_frazil_fractional_thickness_limit)
      call mpas_pool_get_config(ocnConfigs, 'config_specific_heat_sea_water', config_specific_heat_sea_water)
      call mpas_pool_get_config(ocnConfigs, 'config_frazil_heat_of_fusion', config_frazil_heat_of_fusion)
      call mpas_pool_get_config(ocnConfigs, 'config_frazil_sea_ice_density', config_frazil_sea_ice_density)
      call mpas_pool_get_config(ocnConfigs, 'config_frazil_ice_reference_salinity', config_frazil_ice_reference_salinity)
      call mpas_pool_get_config(ocnConfigs, 'config_frazil_maximum_freezing_temperature', config_frazil_maximum_freezing_temperature)

      ! get time step in units of seconds
      timeStep = mpas_get_clock_timestep(domain % clock, ierr=err)
      call mpas_get_timeInterval(timeStep, dt=dt)

      ! initialize frazil tendency fields
      frazilTemperatureTendency = 0.0_RKIND
      frazilSalinityTendency = 0.0_RKIND
      frazilLayerThicknessTendency = 0.0_RKIND

      ! loop over all columns
      !$omp do schedule(runtime) private(kBottomFrazil, k, columnTemperatureMin, sumNewFrazilIceThickness, oceanFreezingTemperature, &
      !$omp                              potential, freezingEnergy, meltingEnergy, newFrazilIceThickness, meltedFrazilIceThickness)
      do iCell=1,nCells
      
        ! find deepest level where frazil can be created
        kBottomFrazil=maxLevelCell(iCell)
        do k=maxLevelCell(iCell), 1, -1
          if(-zMid(k,iCell).lt.config_frazil_maximum_depth) then
            kBottomFrazil=k
            exit
          endif
        enddo
        
        ! find minimum temperature between 1:kBottomFrazil
        columnTemperatureMin = 1.0e30_RKIND
        do k=1,kBottomFrazil
          if(activeTracers(indexTemperature,k,iCell).lt.columnTemperatureMin) columnTemperatureMin=activeTracers(indexTemperature,k,iCell)
        enddo
        
        ! test min temperature agains max freezing temperature to see if we should even consider creating frazil
        if(columnTemperatureMin.gt.config_frazil_maximum_freezing_temperature) cycle
        
        ! initialize the sum of new frazil ice created
        sumNewFrazilIceThickness = 0.0_RKIND
        
        ! loop from maximum depth of frazil creation to surface
        do k = kBottomFrazil, 1, -1
        
          ! get freezing temperature
          oceanFreezingTemperature = ocn_freezing_temperature(activeTracers(indexSalinity,k,iCell))
          
          potential = layerThickness(k,iCell) * config_specific_heat_sea_water &
                    * rho_sw * (activeTracers(indexTemperature,k,iCell) - oceanFreezingTemperature)

          freezingEnergy = max(0.0_RKIND, -potential)
          meltingEnergy = max(0.0_RKIND, potential)
        
          if (freezingEnergy > 0) then
          
            ! new frazil ice formation measured in meters
            newFrazilIceThickness = freezingEnergy / (config_frazil_heat_of_fusion * config_frazil_sea_ice_density)
            
            ! limit the frazil formed appropriately
            newFrazilIceThickness = min(newFrazilIceThickness, layerThickness(k,iCell) * config_frazil_fractional_thickness_limit)
            
            ! compute tendency to thickness, temperature and salinity
            ! layerTendency is scaled so that mass of ice created == mass of ocean water removed
            
            ! layer thickness decreased due to creation of frazil
            ! note: -- this has to be density (not rho_sw) to keep buoyancy equal
            frazilLayerThicknessTendency(k,iCell) = - newFrazilIceThickness * config_frazil_sea_ice_density / density(k,iCell) / dt
            
            ! salt is extracted with the frazil
            frazilSalinityTendency(k,iCell) = - newFrazilIceThickness * config_frazil_ice_reference_salinity / dt
            
            ! ocean fluid temperature is warmed due to creation of frazil
            frazilTemperatureTendency(k,iCell) =  + ( newFrazilIceThickness * config_frazil_heat_of_fusion * config_frazil_sea_ice_density )  &
                                               / (config_specific_heat_sea_water * rho_sw) / dt
            
            ! keep track of sum of frazil ice
            sumNewFrazilIceThickness = sumNewFrazilIceThickness + newFrazilIceThickness
          
          else
          
            ! ocean water is warm enough to melt frazil
          
            ! test to see if there is frazil to be melted
            if (sumNewFrazilIceThickness > 0.0_RKIND) then
            
              ! Frazil melting
              meltedFrazilIceThickness = meltingEnergy / (config_frazil_heat_of_fusion * config_frazil_sea_ice_density)
              
              ! limit melting by what there is to melt
              meltedFrazilIceThickness = min(meltedFrazilIceThickness, sumNewFrazilIceThickness)
              
              ! limit melting by fraction of layer thickness
              meltedFrazilIceThickness = min(meltedFrazilIceThickness, layerThickness(k,iCell)*config_frazil_fractional_thickness_limit)
              
              ! compute tendency to thickness, temperature and salinity
              
              ! layer thickness increases due to melting of frazil
              ! note -- scaling by local ocean density to mimimize surface pressure forcing errors
              frazilLayerThicknessTendency(k,iCell) = + meltedFrazilIceThickness * config_frazil_sea_ice_density / density(k,iCell) / dt
              
              ! salt is released into ocean with the melting frazil
              frazilSalinityTendency(k,iCell) = + meltedFrazilIceThickness * config_frazil_ice_reference_salinity / dt
              
              ! ocean fluid temperature is cooled due to melting of frazil
              frazilTemperatureTendency(k,iCell) = - ( meltedFrazilIceThickness * config_frazil_heat_of_fusion * config_frazil_sea_ice_density )  &
                                                 / (config_specific_heat_sea_water * rho_sw) / dt
              
              ! keep track of new frazil ice
              sumNewFrazilIceThickness = sumNewFrazilIceThickness - meltedFrazilIceThickness
            
            endif  ! if (sumNewFrazilIceThickness > 0.0_RKIND) 
          
          endif  ! if (freezingEnergy < 0)
        
        enddo   ! do k=kBottom,1-1
        
        ! accumulate frazil mass to column total
        ! note: the accumulatedFrazilIceMass (at both time levels) is reset to zero after being sent to the coupler
        accumulatedFrazilIceMassNew(iCell) = accumulatedFrazilIceMassOld(iCell) + sumNewFrazilIceThickness*config_frazil_sea_ice_density
        
        ! sea surface pressure due to the net production of frazil ice
        frazilSurfacePressure(iCell) = accumulatedFrazilIceMassNew(iCell) * gravity
      
      enddo   ! do iCell = 1, nCells
      
      call mpas_timer_stop("frazil", timer_frazil)

   end subroutine ocn_frazil_forcing_build_arrays!}}}

!***********************************************************************
!
!  function ocn_freezing_temperature
!
!> \brief   Computes the freezing temperature of the ocean.
!> \author  Todd Ringler
!> \date    10/29/2015
!> \details 
!>  This routine computes the freezing temperature of the ocean at a given
!>  salinity value.
!
!-----------------------------------------------------------------------
    real (kind=RKIND) function ocn_freezing_temperature(salinity)!{{{
      real (kind=RKIND) :: salinity !< Input: Salinity value of water for freezing temperature
      ocn_freezing_temperature = -1.8
    end function ocn_freezing_temperature!}}}


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
      logical, pointer :: config_use_frazil_ice_formation

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_use_frazil_ice_formation', config_use_frazil_ice_formation)

      frazilFormationOn = .false.

      if(config_use_frazil_ice_formation) then
        frazilFormationOn = .true.
      end if

   end subroutine ocn_frazil_forcing_init!}}}

!***********************************************************************

end module ocn_frazil_forcing

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
