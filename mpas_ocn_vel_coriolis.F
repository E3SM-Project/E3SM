! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_coriolis
!
!> \brief MPAS ocean horizontal momentum mixing driver
!> \author Doug Jacobsen
!> \date   15 September 2011
!> \details
!>  This module contains the routine for computing 
!>  tendencies from the coriolis force.  
!>
!
!-----------------------------------------------------------------------

module ocn_vel_coriolis

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

   public :: ocn_vel_coriolis_tend, &
             ocn_vel_coriolis_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: coriolisOn

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vel_coriolis_tend
!
!> \brief   Computes tendency term for coriolis force
!> \author  Doug Jacobsen
!> \date    15 September 2011
!> \details 
!>  This routine computes the coriolis tendency for momentum
!>  based on current state.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_coriolis_tend(grid, Vor_edge, layerThicknessEdge, normalVelocity, kineticEnergy, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         Vor_edge  !< Input: Potential vorticity on edge
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThicknessEdge  !< Input: Thickness on edge
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity  !< Input: Horizontal velocity
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         kineticEnergy  !< Input: Kinetic Energy

      type (mesh_type), intent(in) :: &
         grid          !< Input: grid information

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

      integer, dimension(:), pointer :: maxLevelEdgeTop, nEdgesOnEdge
      integer, dimension(:,:), pointer :: cellsOnEdge, edgesOnEdge, edgeMask
      real (kind=RKIND), dimension(:,:), pointer :: weightsOnEdge
      real (kind=RKIND), dimension(:), pointer :: dcEdge

      integer :: j, k
      integer :: cell1, cell2, nEdgesSolve, iEdge, eoe
      real (kind=RKIND) :: workpv, q, invLength

      err = 0

      if(.not.coriolisOn) return

      maxLevelEdgeTop => grid % maxLevelEdgeTop % array
      nEdgesOnEdge => grid % nEdgesOnEdge % array
      cellsOnEdge => grid % cellsOnEdge % array
      edgesOnEdge => grid % edgesOnEdge % array
      weightsOnEdge => grid % weightsOnEdge % array
      dcEdge => grid % dcEdge % array

      edgeMask => grid % edgeMask % array

      nEdgesSolve = grid % nEdgesSolve

      do iEdge=1,grid % nEdgesSolve
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         invLength = 1.0 / dcEdge(iEdgE)

         do k=1,maxLevelEdgeTop(iEdge)

            q = 0.0
            do j = 1,nEdgesOnEdge(iEdge)
               eoe = edgesOnEdge(j,iEdge)
               workpv = 0.5 * (Vor_edge(k,iEdge) + Vor_edge(k,eoe))
               q = q + weightsOnEdge(j,iEdge) * normalVelocity(k,eoe) * workpv * layerThicknessEdge(k,eoe) 
            end do

           tend(k,iEdge) = tend(k,iEdge) + edgeMask(k, iEdge) * (q - (   kineticEnergy(k,cell2) - kineticEnergy(k,cell1) ) * invLength )

         end do
      end do

   !--------------------------------------------------------------------

   end subroutine ocn_vel_coriolis_tend!}}}

!***********************************************************************
!
!  routine ocn_vel_coriolis_init
!
!> \brief   Initializes ocean momentum horizontal mixing quantities
!> \author  Doug Jacobsen
!> \date    15 September 2011
!> \details 
!>  This routine initializes a variety of quantities related to 
!>  horizontal velocity mixing in the ocean. Since a variety of 
!>  parameterizations are available, this routine primarily calls the
!>  individual init routines for each parameterization. 
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_coriolis_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! Output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      err = 0

      coriolisOn = .true.

      if(config_disable_vel_coriolis) coriolisOn = .false.

   !--------------------------------------------------------------------

   end subroutine ocn_vel_coriolis_init!}}}

!***********************************************************************

end module ocn_vel_coriolis

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
