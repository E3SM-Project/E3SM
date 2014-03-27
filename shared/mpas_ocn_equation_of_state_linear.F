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

   use mpas_grid_types
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

   subroutine ocn_equation_of_state_linear_density(meshPool, indexT, indexS, tracers, density, err, &
      thermalExpansionCoeff, salineContractionCoeff)!{{{
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !  This module contains routines necessary for computing the density
   !  from model temperature and salinity using an equation of state.
   !
   ! Input: mesh - mesh metadata
   !        s - state: tracers
   !        k_displaced 
   !  If k_displaced<=0, state % density is returned with no displaced
   !  If k_displaced>0,the state % densityDisplaced is returned, and is for
   !  a parcel adiabatically displaced from its original level to level 
   !  k_displaced.  This does not effect the linear EOS.
   !
   ! Output: s - state: computed density
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit none

      type (mpas_pool_type), intent(in) :: meshPool
      integer, intent(in) :: indexT, indexS
      real (kind=RKIND), dimension(:,:,:), intent(in) :: tracers
      real (kind=RKIND), dimension(:,:), intent(inout) :: density
      integer, intent(out) :: err
      real (kind=RKIND), dimension(:,:), intent(out), optional :: &
         thermalExpansionCoeff,  &! Thermal expansion coefficient (alpha), defined as $-1/\rho d\rho/dT$ (note negative sign)
         salineContractionCoeff   ! Saline contraction coefficient (beta), defined as $1/\rho d\rho/dS$

      integer, dimension(:), pointer :: maxLevelCell
      integer :: iCell, k
      integer, pointer :: nCells
      type (dm_info) :: dminfo

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

      err = 0

      do iCell = 1, nCells
         do k = 1, maxLevelCell(iCell)
            ! Linear equation of state
            density(k,iCell) =  config_eos_linear_densityref &
                  - config_eos_linear_alpha * (tracers(indexT,k,iCell) - config_eos_linear_Tref) &
                  + config_eos_linear_beta  * (tracers(indexS,k,iCell) - config_eos_linear_Sref)
         end do
      end do

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
