! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_forcing_windstress
!
!> \brief MPAS ocean wind stress
!> \author Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date   September 2011
!> \details
!>  This module contains the routine for computing 
!>  tendencies from wind stress.  
!
!-----------------------------------------------------------------------

module ocn_vel_forcing_windstress

   use mpas_grid_types
   use mpas_configure

   use ocn_forcing

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

   public :: ocn_vel_forcing_windstress_tend, &
             ocn_vel_forcing_windstress_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: windStressOn

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vel_forcing_windstress_tend
!
!> \brief   Computes tendency term from wind stress
!> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date    15 September 2011
!> \details 
!>  This routine computes the wind stress tendency for momentum
!>  based on current state.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_forcing_windstress_tend(mesh, normalVelocityForcing, surfaceWindStress, layerThicknessEdge, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocityForcing    !< Input: Normal velocity forcing at edges

      real (kind=RKIND), dimension(:), intent(in) :: &
         surfaceWindStress    !< Input: Wind stress at surface

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThicknessEdge     !< Input: thickness at edge

      type (mesh_type), intent(in) :: &
         mesh          !< Input: mesh information

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

      integer :: iEdge, nEdgesSolve, k
      integer, dimension(:), pointer :: maxLevelEdgeTop
      integer, dimension(:,:), pointer :: edgeMask

      real (kind=RKIND) :: transmissionCoeffTop, transmissionCoeffBot, zTop, zBot, remainingStress

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing tendencies
      ! note that the user can choose multiple options and the 
      !   tendencies will be added together
      !
      !-----------------------------------------------------------------

      err = 0

      if(.not.windStressOn) return

      nEdgesSolve = mesh % nEdgesSolve
      maxLevelEdgeTop => mesh % maxLevelEdgeTop % array
      edgeMask => mesh % edgeMask % array

      do iEdge=1,nEdgesSolve
        ! k = min(maxLevelEdgeTop(iEdge), 1)
        ! tend(k,iEdge) =  tend(k,iEdge) + edgeMask(k, iEdge) * surfaceWindStress(iEdge) / config_density0 / layerThicknessEdge(k, iEdge)

        zTop = 0.0_RKIND
        transmissionCoeffTop = ocn_forcing_transmission(zTop)
        remainingStress = 1.0_RKIND
        do k = 1,maxLevelEdgeTop(iEdge)
           zBot = zTop - layerThicknessEdge(k, iEdge)

           transmissionCoeffBot = ocn_forcing_transmission(zBot)

           remainingStress = remainingStress - (transmissionCoeffTop - transmissionCoeffBot)

           tend(k,iEdge) =  tend(k,iEdge) + edgeMask(k, iEdge) * (normalVelocityForcing(k,iEdge) + surfaceWindStress(iEdge) * (transmissionCoeffTop - transmissionCoeffBot)) / config_density0 / layerThicknessEdge(k,iEdge)

           zTop = zBot
           transmissionCoeffTop = transmissionCoeffBot
        enddo

        if ( maxLevelEdgeTop(iEdge) > 0 .and. remainingStress > 0.0_RKIND) then
           tend(maxLevelEdgeTop(iEdge), iEdge) = tend(maxLevelEdgeTop(iEdge), iEdge) + edgeMask(maxLevelEdgeTop(iEdge), iEdge) * surfaceWindStress(iEdge) * remainingStress &
                                                                                     / config_density0 / layerThicknessEdge(maxLevelEdgeTop(iEdge), iEdge)
        end if
      enddo


   !--------------------------------------------------------------------

   end subroutine ocn_vel_forcing_windstress_tend!}}}

!***********************************************************************
!
!  routine ocn_vel_forcing_windstress_init
!
!> \brief   Initializes ocean wind stress forcing
!> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date    September 2011
!> \details 
!>  This routine initializes quantities related to wind stress 
!>  in the ocean. 
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_forcing_windstress_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      windStressOn = .true.

      if(config_disable_vel_windstress) windStressOn = .false.

      err = 0

   !--------------------------------------------------------------------

   end subroutine ocn_vel_forcing_windstress_init!}}}

!***********************************************************************

end module ocn_vel_forcing_windstress

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
