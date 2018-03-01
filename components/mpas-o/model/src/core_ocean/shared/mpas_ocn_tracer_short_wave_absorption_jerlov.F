!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_short_wave_absorption_jerlov
!
!> \brief MPAS ocean tracer short wave
!> \author Doug Jacobsen
!> \date   12/17/12
!> \details
!>  This module contains the routine for computing
!>  short wave tendencies using Jerlov
!
!-----------------------------------------------------------------------

module ocn_tracer_short_wave_absorption_jerlov

   use mpas_derived_types
   use mpas_pool_routines
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

   public :: ocn_tracer_short_wave_absorption_jerlov_tend, &
             ocn_get_jerlov_fraction

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------


   integer, parameter :: num_water_types = 5

   !-----------------------------------------------------------------------
   !
   !   define Jerlov water properties with rfac, depth1, depth2
   !     Jerlov water type :  I       IA      IB      II      III
   !     jerlov_water_type :  1       2       3       4       5
   !
   !-----------------------------------------------------------------------

   real (kind=RKIND), dimension(num_water_types) ::                       &
      rfac   = (/ 0.58_RKIND, 0.62_RKIND, 0.67_RKIND, 0.77_RKIND, 0.78_RKIND /), &
      depth1 = (/ 0.35_RKIND, 0.60_RKIND, 1.00_RKIND, 1.50_RKIND, 1.40_RKIND /), &
      depth2 = (/ 23.0_RKIND, 20.0_RKIND, 17.0_RKIND, 14.0_RKIND, 7.90_RKIND /)

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_short_wave_absorption_jerlov_tend
!
!> \brief   Computes tendency term for surface fluxes
!> \author  Doug Jacobsen
!> \date    12/17/12
!> \details
!>  This routine computes the tendency for tracers based on surface fluxes.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_short_wave_absorption_jerlov_tend(meshPool, forcingPool, index_temperature, layerThickness, &
                                        penetrativeTemperatureFlux, penetrativeTemperatureFluxOBL, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool, forcingPool          !< Input: mesh information

      real (kind=RKIND), dimension(:), intent(in) :: &
        penetrativeTemperatureFlux !< Input: penetrative temperature flux through the surface

      real (kind=RKIND), dimension(:), intent(out) :: &
        penetrativeTemperatureFluxOBL

      real (kind=RKIND), dimension(:,:), intent(in) :: layerThickness !< Input: Layer thicknesses

      integer, intent(in) :: index_temperature

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
         tend          !< Input/Output: velocity tendency

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

      integer :: iCell, k, depLev, nCells
      integer, pointer :: nVertLevels
      integer, dimension(:), pointer :: nCellsArray

      integer, dimension(:), pointer :: maxLevelCell

      real (kind=RKIND) :: depth
      real (kind=RKIND), pointer :: config_surface_buoyancy_depth
      real (kind=RKIND), dimension(:), pointer :: refBottomDepth
      real (kind=RKIND), dimension(:), allocatable :: weights

      err = 0

      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
      call mpas_pool_get_config(ocnConfigs, 'config_surface_buoyancy_depth', config_surface_buoyancy_depth)

      allocate(weights(nVertLevels+1))
      weights = 0.0_RKIND
      weights(1) = 1.0_RKIND

      nCells = nCellsArray( 3 )

      !$omp do schedule(runtime) private(depth, k, depLev)
      do iCell = 1, nCells
         depth = 0.0_RKIND
         do k = 1, maxLevelCell(iCell)
            depth = depth + layerThickness(k, iCell)

            call ocn_get_jerlov_fraction(depth, weights(k+1))
            tend(index_temperature, k, iCell) = tend(index_temperature, k, iCell) + penetrativeTemperatureFlux(iCell) &
                                              * (weights(k) - weights(k+1))
         end do

         depth = 0.0_RKIND
         do k=1,maxLevelCell(iCell)
           depth = depth + layerThickness(k,iCell)
           if(depth > abs(config_surface_buoyancy_depth)) exit
         enddo

         if(k == maxLevelCell(iCell) .or. k == 1) then
           depLev=2
         else
           depLev=k
         endif
         penetrativeTemperatureFluxOBL(iCell)=penetrativeTemperatureFlux(iCell)*weights(depLev)

      end do
      !$omp end do

      deallocate(weights)

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_short_wave_absorption_jerlov_tend!}}}

!***********************************************************************
!
!  routine ocn_init_jerlov_fractions
!
!> \brief   Initializes short wave absorption fractions
!> \author  Doug Jacobsen
!> \date    12/17/12
!> \details
!>  Computes fraction of solar short-wave flux penetrating to
!>  specified depth due to exponential decay in Jerlov water type.
!>  Reference : two band solar absorption model of Simpson and
!>     Paulson (1977)
!
!-----------------------------------------------------------------------
   subroutine ocn_get_jerlov_fraction(depth, weight)!{{{
!  Note: below 200m the solar penetration gets set to zero,
!     otherwise the limit for the exponent ($+/- 5678$) needs to be
!     taken care of.

      real (kind=RKIND), intent(in) :: depth !< Input: Depth of bottom of cell
      real (kind=RKIND), intent(out) :: weight !< Output: Weight for Jerlov absorption

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

!
      integer :: k, nVertLevels
      integer,  parameter :: num_water_types = 5  ! max number of different water types

      ! I don't understand what the previous two lines are for.  They appear unnecessary

      real (kind=RKIND), parameter :: depth_cutoff = -200.0_RKIND

      integer, pointer :: config_jerlov_water_type

!-----------------------------------------------------------------------
!
!  compute absorption fraction
!
!-----------------------------------------------------------------------

      call mpas_pool_get_config(ocnConfigs, 'config_jerlov_water_type', config_jerlov_water_type)

      if (-depth < depth_cutoff) then
         weight = 0.0_RKIND
      else
         weight = rfac(config_jerlov_water_type) * exp(-depth/depth1(config_jerlov_water_type)) &
                  + (1.0_RKIND - rfac(config_jerlov_water_type)) * exp(-depth/depth2(config_jerlov_water_type))
      endif
   end subroutine ocn_get_jerlov_fraction!}}}

end module ocn_tracer_short_wave_absorption_jerlov

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
