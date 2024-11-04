! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tidal_forcing
!
!> \brief MPAS ocean tidal forcing module
!> \author Phillip J. Wolfram
!> \date   04/05/2019
!> \details
!>  This module contains routines for the tidal forcing.
!
!-----------------------------------------------------------------------

module ocn_tidal_forcing

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
   use ocn_subgrid

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

   public :: ocn_tidal_forcing_build_array, &
             ocn_tidal_forcing_layer_thickness, &
             ocn_tidal_forcing_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: tidalFluxOn

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tidal_forcing_layer_thickness
!
!> \brief   Add tendency due to tidal forcing
!> \author  Phillip J. Wolfram
!> \date    04/05/2019
!> \details
!>  This routine adds a tendency to layer thickness due to tidal forcing
!
!-----------------------------------------------------------------------

   subroutine ocn_tidal_forcing_layer_thickness(forcingPool, layerThicknessTend, err)!{{{

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
      real (kind=RKIND), dimension(:,:), pointer :: tidalLayerThicknessTendency

      err = 0

      if ( .not. tidalFluxOn ) return

      call mpas_timer_start("tidal thickness tendency")

      call mpas_pool_get_array(forcingPool, 'tidalLayerThicknessTendency', &
                                             tidalLayerThicknessTendency)

      ! frazil fields are needed only over 0 and 1 halos
      nCells = nCellsHalo( 1 )

      ! Build surface fluxes at cell centers
#ifdef MPAS_OPENACC
      !$acc enter data copyin(tidalLayerThicknessTendency)

      !$acc parallel loop present(layerThicknessTend, maxLevelCell, minLevelCell, &
      !$acc     tidalLayerThicknessTendency) &
      !$acc     private(k)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(k)
#endif
      do iCell = 1, nCells
        do k = minLevelCell(iCell), maxLevelCell(iCell)
          layerThicknessTend(k,iCell) = layerThicknessTend(k,iCell) + &
                                        tidalLayerThicknessTendency(k,iCell)

        end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

#ifdef MPAS_OPENACC
      !$acc exit data delete(tidalLayerThicknessTendency)
#endif

      call mpas_timer_stop("tidal thickness tendency")

   end subroutine ocn_tidal_forcing_layer_thickness !}}}


!***********************************************************************
!
!  routine ocn_tidal_forcing_build_array
!
!> \brief   Computes tidal forcing functions.
!> \author  Phillip J. Wolfram
!> \date    04/05/2019
!> \details
!>   ocn_tidal_forcing_build_array computes the tendencies to layer thickness
!>   due to tidal forcing of layer thickness
!>
!>   this routine should be call at the beginning of whatever time stepping method is utilized
!>      and the tendencies should be retieved when building up the RHS of the thickess, temperature
!>      and salinity equations.
!>
!>   this routine is only applicable to the thickness fields
!
!-----------------------------------------------------------------------

   subroutine ocn_tidal_forcing_build_array(domain, meshPool, forcingPool, statePool, err)!{{{

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

      real (kind=RKIND), dimension(:,:), pointer :: tidalLayerThicknessTendency
      real (kind=RKIND), dimension(:), pointer :: tidalInputMask
      real (kind=RKIND), dimension(:), pointer :: tidalBCValue
      real (kind=RKIND), dimension(:), pointer :: bottomDepth

      integer :: iCell, k, nCells
      integer, dimension(:), pointer :: nCellsArray
      integer, pointer, dimension(:) :: minLevelCell, maxLevelCell

      real (kind=RKIND), pointer, dimension(:)     :: ssh
      real (kind=RKIND), pointer, dimension(:,:)   :: layerThickness

      ! local variables
      real (kind=RKIND) :: totalDepth, tidalHeight

      character (len=StrKIND), pointer :: xtime

      if ( .not. tidalFluxOn ) return

      call mpas_pool_get_array(statePool, 'ssh', ssh, 1)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
      call mpas_pool_get_array(forcingPool, 'tidalLayerThicknessTendency', &
                                             tidalLayerThicknessTendency)
      call mpas_pool_get_array(forcingPool, 'tidalInputMask', tidalInputMask)
      call mpas_pool_get_array(forcingPool, 'tidalBCValue', tidalBCValue)

      call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
      call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)

      ! tidal fields are needed only over 0 and 1 halos
      nCells = nCellsArray( 2 )

      ! loop over all columns
      do iCell = 1, nCells

        ! compute the tidalHeight
        if (trim(config_tidal_forcing_model) == 'monochromatic') then
          tidalHeight = config_tidal_forcing_monochromatic_amp * &
            SIN(2.0_RKIND*pi/config_tidal_forcing_monochromatic_period * daysSinceStartOfSim - &
            pi*config_tidal_forcing_monochromatic_phaseLag/180.0_RKIND) - &
            config_tidal_forcing_monochromatic_baseline
        elseif (trim(config_tidal_forcing_model) == 'linear') then
           tidalHeight = max(config_tidal_forcing_linear_min, &
                             config_tidal_forcing_linear_baseline + &
                             config_tidal_forcing_linear_rate * daysSinceStartOfSim)
        !else if (trim(config_tidal_forcing_type) == 'data') then
        !  ! data option
        !  ! pass
        else
          tidalHeight = 0.0_RKIND
        end if

        ! ensure tidalHeight can't push flow below anticipated minimum
        if (config_use_wetting_drying .and. tidalInputMask(iCell) == 1.0_RKIND) then
          ! ensure that tidal height can't force below total minimum thickness
          ! condition wrong to ensure that there isn't any drying according to criteria
          if ( config_use_subgrid_wetting_drying ) then
            tidalHeight = max(tidalHeight, subgridSShCellTableRange(3,iCell) ) ;      
          else
           tidalHeight = max(-bottomDepth(iCell) + (float(maxLevelCell(iCell))+1.0_RKIND)*config_drying_min_cell_height, tidalHeight)
          end if
        end if

        ! compute total depth for relative thickness contribution
        totalDepth = 0.0_RKIND
        do k = minLevelCell(iCell), maxLevelCell(iCell)
          totalDepth = totalDepth + layerThickness(k,iCell)
        end do

        tidalLayerThicknessTendency(:,iCell) = 0.0_RKIND
        if (trim(config_tidal_forcing_type) == 'thickness_source') then
          ! distribute tidal forcing tendency fractionally over water column
          do k = minLevelCell(iCell), maxLevelCell(iCell)
            tidalLayerThicknessTendency(k,iCell) = tidalInputMask(iCell) / config_use_tidal_forcing_tau &
              * (layerThickness(k,iCell)/totalDepth) * (tidalHeight - ssh(iCell))
          end do
          ! else if
          ! using the 'direct' option to immediately force the free surface as a prescribe boundary
          ! in mpas_ocn_time_integration_rk4.F
        end if

        ! store tidal boundary cell values (e.g., for config_tidal_forcing_type == 'direct')
        tidalBCValue(iCell) = tidalHeight

      end do   ! do iCell = 1, nCells

   end subroutine ocn_tidal_forcing_build_array !}}}


!***********************************************************************
!
!  routine ocn_tidal_forcing_init
!
!> \brief   Initializes ocean tidal forcing module.
!> \author  Phillip J. Wolfram
!> \date    04/05/2019
!> \details
!>  This routine initializes the ocean tidal forcing module and variables.
!
!-----------------------------------------------------------------------

   subroutine ocn_tidal_forcing_init(err)!{{{

      integer, intent(out) :: err !< Output: error flag

      err = 0

      tidalFluxOn = .false.

      if(config_use_tidal_forcing) then
        tidalFluxOn = .true.
      end if

   end subroutine ocn_tidal_forcing_init!}}}


!***********************************************************************

end module ocn_tidal_forcing

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
