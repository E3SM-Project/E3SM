! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_tidal_potential
!
!> \brief MPAS ocean pressure gradient module
!> \author Steven Brus
!> \date April 2020
!> \details
!>  This module contains the routine for computing
!>  tendency from the tidal potential.
!>
!
!-----------------------------------------------------------------------

module ocn_vel_tidal_potential

   use mpas_timer
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_log

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

   public :: ocn_vel_tidal_potential_tend, &
             ocn_vel_tidal_potential_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: tidalPotentialOn

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vel_tidal_potential_tend
!
!> \brief   Computes tendency term for tidal potential
!> \author  Steven Brus
!> \date    April 2020
!> \details
!>  This routine computes the tidal potential tendency for momentum
!>  based on current state.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_tidal_potential_tend(meshPool, forcingPool, ssh, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(in) :: &
         ssh !< Input: Sea surface height

      type (mpas_pool_type), intent(in) :: meshPool          !< Input: mesh information
      type (mpas_pool_type), intent(in) :: forcingPool       !< Input: forcinginformation

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
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

      integer :: iEdge, k, cell1, cell2, nEdges
      integer, dimension(:), pointer :: nEdgesArray
      integer, dimension(:), pointer :: maxLevelEdgeTop, maxLevelCell
      integer, dimension(:,:), pointer :: cellsOnEdge, edgeMask

      real (kind=RKIND), dimension(:), pointer :: dcEdge
      real (kind=RKIND), dimension(:), pointer :: tidalPotentialEta
      real (kind=RKIND), pointer :: config_self_attraction_and_loading_beta
      real (kind=RKIND) :: invdcEdge
      real (kind=RKIND) :: potentialGrad
      logical, pointer :: config_use_tidal_potential_forcing

      err = 0

      if (.not. tidalPotentialOn) return

      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)
      call mpas_pool_get_array(forcingPool, 'tidalPotentialEta', tidalPotentialEta)
      call mpas_pool_get_config(ocnConfigs, 'config_use_tidal_potential_forcing', config_use_tidal_potential_forcing)
      call mpas_pool_get_config(ocnConfigs, 'config_self_attraction_and_loading_beta', config_self_attraction_and_loading_beta)

      nEdges = nEdgesArray( 1 )

      !$omp do schedule(runtime) private(cell1, cell2, invdcEdge, potentialGrad, k)
      do iEdge=1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         invdcEdge = 1.0_RKIND / dcEdge(iEdge)

         potentialGrad = - gravity * invdcEdge * ( tidalPotentialEta(cell2) - tidalPotentialEta(cell1) &
                                   + config_self_attraction_and_loading_beta*(ssh(cell2) - ssh(cell1)))

         do k=1,maxLevelEdgeTop(iEdge)
            tend(k,iEdge) = tend(k,iEdge) -  edgeMask(k,iEdge) * potentialGrad
         end do
      end do
      !$omp end do

   !--------------------------------------------------------------------

   end subroutine ocn_vel_tidal_potential_tend!}}}

!***********************************************************************
!
!  routine ocn_vel_tidal_potential_init
!
!> \brief   Initializes ocean momentum tidal potential gradient
!> \author  Steven Brus
!> \date    April 2020
!> \details
!>  This routine initializes parameters required for the computation of the
!>  tidal potential gradient.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_tidal_potential_init(err)!{{{

   !--------------------------------------------------------------------


      !-----------------------------------------------------------------
      !
      ! Output Variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag


      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------
      logical, pointer :: config_use_tidal_potential_forcing

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_use_tidal_potential_forcing', config_use_tidal_potential_forcing)

      tidalPotentialOn= .false.

      if (config_use_tidal_potential_forcing) tidalPotentialOn= .true.

   !--------------------------------------------------------------------

   end subroutine ocn_vel_tidal_potential_init!}}}

!***********************************************************************

end module ocn_vel_tidal_potential

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
