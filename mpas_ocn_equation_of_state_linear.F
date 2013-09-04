! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
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
   use mpas_configure

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

   subroutine ocn_equation_of_state_linear_density(grid, indexT, indexS, tracers, density, err, &
       dDensityDT, dDensityDS)!{{{
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !  This module contains routines necessary for computing the density
   !  from model temperature and salinity using an equation of state.
   !
   ! Input: grid - grid metadata
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

      type (mesh_type), intent(in) :: grid
      integer, intent(in) :: indexT, indexS
      real (kind=RKIND), dimension(:,:,:), intent(in) :: tracers
      real (kind=RKIND), dimension(:,:), intent(inout) :: density
      integer, intent(out) :: err
      real (kind=RKIND), dimension(:,:), intent(out), optional :: &
         dDensityDT,  &! derivative of density with respect to temperature
         dDensityDS    ! derivative of density with respect to salinity

      integer, dimension(:), pointer :: maxLevelCell
      integer :: nCells, iCell, k
      type (dm_info) :: dminfo

      maxLevelCell  => grid % maxLevelCell % array
      nCells        =  grid % nCells

      err = 0

      do iCell=1,nCells
         do k=1,maxLevelCell(iCell)
            ! Linear equation of state
            density(k,iCell) =  config_eos_linear_densityref &
                  - config_eos_linear_alpha * (tracers(indexT,k,iCell)-config_eos_linear_Tref) &
                  + config_eos_linear_beta  * (tracers(indexS,k,iCell)-config_eos_linear_Sref)
         end do
      end do

      if (present(dDensityDT)) then
         dDensityDT = - config_eos_linear_alpha
      endif

      if (present(dDensityDS)) then
         dDensityDS = config_eos_linear_beta
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

   !--------------------------------------------------------------------

   end subroutine ocn_equation_of_state_linear_init!}}}

!***********************************************************************

end module ocn_equation_of_state_linear

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
