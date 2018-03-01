! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_equation_of_state_linear
!
!> \brief MPAS ocean equation of state driver
!> \author Mark Petersen, Todd Ringler
!> \date   September 2011
!> \details
!>  This module contains the main driver routine for calling
!>  the equation of state.
!
!-----------------------------------------------------------------------

module ocn_equation_of_state_linear

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use ocn_constants
   use mpas_dmpar
   use mpas_log

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

   public :: ocn_equation_of_state_linear_density, &
             ocn_equation_of_state_linear_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------
   real (kind=RKIND), pointer :: config_eos_linear_densityref
   real (kind=RKIND), pointer :: config_eos_linear_alpha
   real (kind=RKIND), pointer :: config_eos_linear_beta
   real (kind=RKIND), pointer :: config_eos_linear_Tref
   real (kind=RKIND), pointer :: config_eos_linear_Sref



!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_equation_of_state_linear_density
!
!> \brief   Calls equation of state
!> \author  Mark Petersen, Todd Ringler
!> \date    September 2011
!> \details 
!>  This routine uses a linear equation of state to update the density
!
!-----------------------------------------------------------------------

   subroutine ocn_equation_of_state_linear_density(meshPool, k_displaced, displacement_type, indexT, indexS, tracers, density, err, &
      tracersSurfaceLayerValue, thermalExpansionCoeff, salineContractionCoeff)!{{{
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !  This module contains routines necessary for computing the density
   !  from model temperature and salinity using an equation of state.
   !
   ! Input: mesh - mesh metadata
   !        s - state: tracers
   !
   ! Output: s - state: computed density
   !
   !
   !>  While somewhat unnecessary, we make the interface and capability
   !>  of linear eos to be identical to nonlinear eos
   !>
   !>  Density can be computed in-situ using k_displaced=0 and
   !>      displacement_type = 'relative'.
   !>
   !>  Potential density (referenced to top layer) can be computed
   !>      using k_displaced=1 and displacement_type = 'absolute'.
   !>
   !>  The density of SST/SSS after adiabatic displacement to each layer
   !>      can be computed using displacement_type = 'surfaceDisplaced'.
   !>
   !>  When using displacement_type = 'surfaceDisplaced', k_displaced is
   !>      ignored and tracersSurfaceLayerValue must be present.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit none

      type (mpas_pool_type), intent(in) :: meshPool
      character(len=*), intent(in) :: displacement_type
      integer, intent(in) :: k_displaced, indexT, indexS
      real (kind=RKIND), dimension(:,:,:), intent(in) :: tracers
      real (kind=RKIND), dimension(:,:), intent(inout) :: density
      integer, intent(out) :: err
      real (kind=RKIND), dimension(:,:), intent(in), optional :: tracersSurfaceLayerValue
      real (kind=RKIND), dimension(:,:), intent(out), optional :: &
         thermalExpansionCoeff,  &! Thermal expansion coefficient (alpha), defined as $-1/\rho d\rho/dT$ (note negative sign)
         salineContractionCoeff   ! Saline contraction coefficient (beta), defined as $1/\rho d\rho/dS$

      integer, dimension(:), pointer :: maxLevelCell
      integer :: iCell, k, k_displaced_local, k_ref
      integer, pointer :: nCells, nVertLevels
      character(len=60) :: displacement_type_local
      type (dm_info) :: dminfo

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      err = 0

      ! copy some intent(in) into local work space
      displacement_type_local = trim(displacement_type)
      k_displaced_local = k_displaced

      ! test of request to address out of bounds
      if (displacement_type_local == 'absolute' .and.   &
         (k_displaced_local <= 0 .or. k_displaced_local > nVertLevels) ) then
         call mpas_log_write('Abort: In equation_of_state_linear' // &
             ' k_displaced must be between 1 and nVertLevels for ' // &
             'displacement_type = absolute', MPAS_LOG_CRIT)
      endif

      ! if surfaceDisplaced, then compute density at all levels based on surface values
      if (displacement_type_local == 'surfaceDisplaced') then
         do iCell = 1, nCells
            do k = 1, maxLevelCell(iCell)
               ! Linear equation of state
               density(k,iCell) =  config_eos_linear_densityref &
                     - config_eos_linear_alpha * (tracersSurfaceLayerValue(indexT,iCell) - config_eos_linear_Tref) &
                     + config_eos_linear_beta  * (tracersSurfaceLayerValue(indexS,iCell) - config_eos_linear_Sref)
            end do
         end do
      endif


      ! if absolute, then compute density at all levels based on pressure of k_displaced value
      ! but since linear EOS does not (at present) have a pressure dependency, this just returns density
      if (displacement_type_local == 'absolute') then
         do iCell = 1, nCells
            do k = 1, maxLevelCell(iCell)
               ! Linear equation of state
               density(k,iCell) =  config_eos_linear_densityref &
                     - config_eos_linear_alpha * (tracers(indexT,k,iCell) - config_eos_linear_Tref) &
                     + config_eos_linear_beta  * (tracers(indexS,k,iCell) - config_eos_linear_Sref)
            end do
         end do
      endif

      ! if relative, then compute density at all levels based on k+k_displaced pressure value
      ! but since (at present) linear EOS has not dependence on pressure, it returns density
      if (displacement_type_local == 'relative') then
         do iCell = 1, nCells
            do k = 1, maxLevelCell(iCell)
               ! Linear equation of state
               density(k,iCell) =  config_eos_linear_densityref &
                     - config_eos_linear_alpha * (tracers(indexT,k,iCell) - config_eos_linear_Tref) &
                     + config_eos_linear_beta  * (tracers(indexS,k,iCell) - config_eos_linear_Sref)
            end do
         end do
      endif

      if (present(thermalExpansionCoeff)) then
         do iCell = 1, nCells
            do k = 1, maxLevelCell(iCell)
               thermalExpansionCoeff(k,iCell) = config_eos_linear_alpha / density(k,iCell)
            end do
         end do
      endif

      if (present(salineContractionCoeff)) then
         do iCell = 1, nCells
            do k = 1, maxLevelCell(iCell)
               salineContractionCoeff(k,iCell) = config_eos_linear_beta / density(k,iCell)
            end do
         end do
      endif

   end subroutine ocn_equation_of_state_linear_density!}}}

!***********************************************************************
!
!  routine ocn_equation_of_state_linear_init
!
!> \brief   Initializes ocean momentum horizontal mixing quantities
!> \author  Mark Petersen, Todd Ringler
!> \date    September 2011
!> \details 
!>  This routine initializes a variety of quantities related to 
!>  horizontal velocity mixing in the ocean. Since a variety of 
!>  parameterizations are available, this routine primarily calls the
!>  individual init routines for each parameterization. 
!
!-----------------------------------------------------------------------

   subroutine ocn_equation_of_state_linear_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err

      integer :: err1, err2

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_eos_linear_densityref', config_eos_linear_densityref)
      call mpas_pool_get_config(ocnConfigs, 'config_eos_linear_alpha', config_eos_linear_alpha)
      call mpas_pool_get_config(ocnConfigs, 'config_eos_linear_beta', config_eos_linear_beta)
      call mpas_pool_get_config(ocnConfigs, 'config_eos_linear_Tref', config_eos_linear_Tref)
      call mpas_pool_get_config(ocnConfigs, 'config_eos_linear_Sref', config_eos_linear_Sref)

   !--------------------------------------------------------------------

   end subroutine ocn_equation_of_state_linear_init!}}}

!***********************************************************************

end module ocn_equation_of_state_linear

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
