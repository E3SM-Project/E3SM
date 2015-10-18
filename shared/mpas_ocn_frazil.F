! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_frazil
!
!> \brief MPAS ocean frazil formation module
!> \author Todd Ringler
!> \date   10/19/2015
!> \details
!>  This module contains routines for the formation of frazil ice.
!
!-----------------------------------------------------------------------

module ocn_frazil

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timekeeping
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

   public :: ocn_frazil_formation, &
             ocn_frazil_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: frazilFormationOn

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_frazil_formation
!
!> \brief   Performs the formation of frazil within the ocean.
!> \author  Todd Ringler
!> \date    10/19/2015
!> \details 
!>   ocn_frazil_formation computes the tendencies to layer thickness, temperature and salinity
!>	due to the creation and possible melting of frazil ice
!>
!>   these tendencies can be retrieved at any point by calling into ocn_frazil_*_tendency routines
!>	where * is layerThickness, temperature or salinity
!>
!>   the pressure exerted by the frazil on the ocean "surface" can be retrieved by calling into
!> 	ocn_frazil_surface_pressure
!>
!>   this routine should be call at the beginning of whatever time stepping method is utilized
!>      and the tendencies should be retieved when building up the RHS of the thickess, temperature
!>      and salinity equations.
!>
!>   this routine is only applicable to the thickness and active tracer fields
!
!-----------------------------------------------------------------------

   subroutine ocn_frazil_formation(meshPool, statePool, tendPool, tracers, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information
      type (mpas_pool_type), intent(in) :: statePool !< Input: State information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(out) :: tendPool !< Output: Tendency information
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
      integer :: iCell, k, kBottomFrazil
      integer, pointer :: nCells, nVertLevels

config_frazil_heat_of_fusion
config_frazil_sea_ice_density
config_frazil_fractional_thickness_limit


real (kind=RKIND) :: newFrazilIceThickness
real (kind=RKIND) :: meltedFrazilIceThickness

      type (mpas_pool_type) :: tracerPool

      real (kind=RKIND), pointer, dimension(:) :: accumulatedFrazilIceThickness
      real (kind=RKIND), pointer, dimension(:) :: accumulatedFrazilIceMass
      real (kind=RKIND), pointer, dimension(:,:) :: zMid
      real (kind=RKIND), pointer, dimension(:,:) :: density
      real (kind=RKIND), pointer, dimension(:,:) :: layerThickness


layerTendencyFrazil
temperatureTendencyFrazil
salinityTendencyFrazil
surfacePressureTendencyFrazil

      integer, dimension(:), pointer :: maxLevelCell
      integer :: indexTemperature !< Input: Index in tracers array for temperature
      integer :: indexSalinity !< Input: Index in tracers array for salinity

      real (kind=RKIND) :: kBottomFrazil  ! k index where testing for frazil begins
      real (kind=RKIND) :: potential  ! scalar holding freezing/melt potential
      real (kind=RKIND) :: freezingEnergy  ! energy available for freezing, positive definite
      real (kind=RKIND) :: meltingEnergy ! energy available for melting, positive definite


      ! if frazil is not enabled, return
      if(.not. frazilFormationOn) return


      do block -----

         ! get dimensions
         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
         call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

         ! get configure parameters
		config_frazil_maximum_depth

         ! get mesh fields
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

         ! get state fields
            layerThickness
            activeTracers
            temperature index
            salinity index

         ! get diagnostic fields
		zMid
        density
        accumulatedFrazilIceThickness

        ! get tendency fields from pool


         ! loop over all columns
         do iCell=1,nCells

            ! reset frazil thickness and mass to be zero for each column
            accumulatedFrazilIceThickness(iCell) = 0.0_RKIND
            accumulatedFrazilIceMass(iCell) = 0.0_RKIND

            ! find deepest level where frazil can be created
            do k=maxLevelCell(iCell), 1, -1
               if(-zMid(k,iCell).lt.config_frazil_maximum_depth) then
                  kBottomFrazil=k
                  exit
               endif
            enddo

            ! loop from maximum depth of frazil creation to surface
            do k = kBottomFrazil, 1, -1

               potential =   layerThickness(k,iCell) * config_specific_heat_sea_water &
                           * density(k,iCell) * (temperature(k,iCell) - oceanFreezingTemperature)
               freezingEnergy = max(0.0_RKIND, -potential)
               meltingEnergy = max(0.0_RKIND, potential)

              if (freezingEnergy < 0) then

                 ! new frazil ice formation measured in meters
                 newFrazilIceThickness = freezingEnergy / (config_frazil_heat_of_fusion * config_frazil_sea_ice_density)

                 ! limit the frazil formed appropriately
                 newFrazilIceThickness = min(newFrazilIceThickness, layerThickness(k,iCell) * config_frazil_fractional_thickness_limit)

                 ! compute increments to thickness, temperature and salinity
                 ! layerTendencyFrazil is scaled so that mass of ice created == mass of ocean water removed
                 layerTendencyFrazil(k,iCell) = -newFrazilIceThickness * config_frazil_sea_ice_density / density(k,iCell)
                 saltTendencyFrazil(k,iCell) = -newFrazilIceThickness * config_frazil_iceReferenceSalinity
                 temperatureTendencyFrazil(k,iCell) =  (newFrazilIceThickness * config_frazil_heat_of_fusion * config_frazil_sea_ice_density)  &
                                                    / (config_specific_heat_sea_water * density(k,iCell))

                 ! accumulate frazil
                 accumulatedFrazilIceThickness(iCell) = accumulatedFrazilIceThickness(iCell) + newFrazilIceThickness
                 accumulatedFrazilIceMass(iCell) = accumulatedFrazilIceMass(iCell) + newFrazilIceThickness*config_frazil_sea_ice_density

              else

              ! ocean water is warm enough to melt frazil

              ! test to see if there is frazil to be melted
              if (accumulatedFrazilIceThickness(iCell) > 0) then

                 ! Frazil melting
                 meltedFrazilIceThickness = meltingEnergy / (config_frazil_heat_of_fusion * config_frazil_sea_ice_density)

                 ! limit melting by what there is to melt
                 meltedFrazilIceThickness = min(meltedFrazilIceThickness, accumulatedFrazilIceThickness(iCell))
 
                 ! limit  melting by fraction of layer thickness
                 meltedFrazilIceThickness = min(meltedFrazilIceThickness, layerThickness(k,iCell)*config_frazil_fractional_thickness_limit)

                 ! compute increments to thickness, temperature and salinity
                 layerTendencyFrazil(k,iCell) = meltedFrazilIceThickness
                 saltTendencyFrazil(k,iCell) = meltedFrazilIceThickness * config_frazil_sea_ice_reference_salinity
                 temperatureTendencyFrazil(k,iCell) = -(meltedFrazilIceThickness * config_frazil_heat_of_fusion * config_frazil_sea_ice_density)  &
                                                    / (config_specific_heat_sea_water * density(k,iCell))

                 ! deaccumulate frazil
                 accumulatedFrazilIceThickness(iCell) = accumulatedFrazilIceThickness(iCell) - meltedFrazilThickness
                 accumulatedFrazilIceMass(iCell) = accumulatedFrazilIceMass(iCell) - meltedFrazilIceThickness*config_frazil_sea_ice_density
 
              endif  ! if (freezingEnergy < 0)

            ! convert tendencies to rates
            ! each of these tendencies can be access through public subroutines below
            layerTendencyFrazil(k,iCell) = layerTendencyFrazil(k,iCell) / dt
            saltTendencyFrazil(k,iCell) = saltTendencyFrazil(k,iCell) / dt
            temperatureTendencyFrazil(k,iCell) = temperatureTendencyFrazil(k,iCell) / dt

          enddo   ! do k=kBottom,1-1

          ! sea surface pressure tendency from frazil ice
          ! note: surfacePressureFrazil should incrememented by surfacePressureTendencyFrazil * dt
          ! note: surfacePressureFrazil should be reset to zero after sending to coupler
          surfacePressureTendencyFrazil(iCell) = accumulatedFrazilIceMass(iCell) * gravity / dt

        enddo   ! do iCell = 1, nCells

     enddo  ! iBlock

   end subroutine ocn_frazil_formation!}}}

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
!  routine ocn_frazil_init
!
!> \brief   Initializes ocean frazil ice module.
!> \author  Todd Ringler
!> \date    10/19/2015
!> \details 
!>  This routine initializes the ocean frazil ice module and variables..
!
!-----------------------------------------------------------------------

   subroutine ocn_frazil_init(err)!{{{

      integer, intent(out) :: err !< Output: error flag
      logical, pointer :: config_use_frazil_ice_formation

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_use_frazil_ice_formation', config_use_frazil_ice_formation)

      frazilFormationOn = .false.

      if(config_use_frazil_ice_formation) then
        frazilFormationOn = .true.
      end if

   end subroutine ocn_frazil_init!}}}

!***********************************************************************

end module ocn_sea_ice

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
