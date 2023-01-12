! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vvel_vert_advection
!
!> \brief MPAS vertical advection for (non-hydrostatic)
!> \author Luke Van Roekel, Darren Engwirda, Sara Calandrini
!> \date   January 2021
!> \details
!>  This module contains routines for advection of (non-hydrostatic) vertical 
!>  velocity.
!
!-------------------------------------------------------------------------------

module ocn_vvel_vert_advection

   ! module includes
   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_io_units
   use mpas_threading
   use ocn_mesh
   use ocn_tracer_advection_shared

   implicit none
   private
   save

   ! public method interfaces
   public :: ocn_vvel_vert_advection_tend, &
             ocn_vvel_vert_advection_init

   contains

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_vvel_vert_advection_tend
!
!> \brief MPAS vertical advection tendency for (non-hydrostatic) vertical velocity
!> \author Luke Van Roekel, Darren Engwirda, Sara Calandrini
!> \date   12/22/2020, updated 11/20/2021
!> \details
!>  NOTE this is only active for non-hydrostatic mode
!
!-----------------------------------------------------------------------
   
   subroutine ocn_vvel_vert_advection_tend(vertVelocityNonhydro, layerThickness, dt, &
      vertAleTransportTop, tend)!{{{

      implicit none 

      !*** input parameters 

      real (kind=RKIND), dimension(:,:), intent(in) :: & 
         vertVelocityNonhydro !< Input: non-hydrostatic vertical velocity
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness !< Input: layer thickness at cells
      real (kind=RKIND), intent(in) :: dt !< Input: time-step
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         vertAleTransportTop !< Input: velocity through layer top levels
      
      !*** output parameters

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend !< Output: vvel tendencies

      integer :: iCell, k 
      real (kind=RKIND) :: levelThickness, vertVelocityUp, vertVelocityDn

      ! -w_ale * dw/dz
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp private(k, levelThickness, vertVelocityUp, vertVelocityDn)
      do iCell = 1, nCellsOwned
         ! w_ale = 0.0 at k==1
         do k = 2, maxLevelCell(iCell)
            levelThickness = &
               0.5_RKIND * (layerThickness(k-1,iCell) + layerThickness(k,iCell))

            vertVelocityUp = &
               0.5_RKIND * (vertVelocityNonhydro(k-1,iCell) + vertVelocityNonhydro(k,iCell))
            vertVelocityDn = &
               0.5_RKIND * (vertVelocityNonhydro(k+1,iCell) + vertVelocityNonhydro(k,iCell))

            tend(k,iCell) = tend(k,iCell) - &
               vertAleTransportTop(k,iCell) * &
                  (vertVelocityUp - vertVelocityDn) / levelThickness
         end do
      end do
      !$omp end do
      !$omp end parallel

   end subroutine ocn_vvel_vert_advection_tend!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_vvel_vert_advection_init
!
!> \brief MPAS initialize vertical velocity advection tendency.
!> \author Luke Van Roekel, Darren Engwirda, Sara Calandrini
!> \date   January 2021, updated November 2021
!> \details
!>  This routine initializes constants and choices for the (non-hydrostatic) 
!>  vertical velocity advection tendency
!
!-------------------------------------------------------------------------------

   subroutine ocn_vvel_vert_advection_init(err) !{{{

      implicit none

      !*** input parameters

      !*** output parameters

      integer, intent(out) :: err !< [out] Error Flag  

      err = 0 ! set error code to success

      ! don't have actual choices for vvel advection here, unlike tracers

   end subroutine ocn_vvel_vert_advection_init!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

end module ocn_vvel_vert_advection

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
