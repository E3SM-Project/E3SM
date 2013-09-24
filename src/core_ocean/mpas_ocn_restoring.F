! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_restoring
!
!> \brief MPAS ocean restoring
!> \author Doug Jacobsen
!> \date   19 September 2011
!> \details
!>  This module contains the main driver routine for computing 
!>  tendencies for restoring.
!
!-----------------------------------------------------------------------

module ocn_restoring

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

   public :: ocn_restoring_tend, &
             ocn_restoring_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: restoringOn !< Flag to turn on/off resotring

   real (kind=RKIND) :: temperatureTimeScale, salinityTimeScale !< restoring timescales


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_restoring_tend
!
!> \brief   Computes tendency term for restoring
!> \author  Doug Jacobsen
!> \date    19 September 2011
!> \details 
!>  This routine computes the restoring tendency for tracers
!>  based on current state.
!
!-----------------------------------------------------------------------

   subroutine ocn_restoring_tend(mesh, layerThickness, indexT, indexS, tracers, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness    !< Input: thickness

      type (mesh_type), intent(in) :: &
         mesh          !< Input: mesh information

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
        tracers !< Input: tracer quantities

      integer, intent(in) :: indexT !< Input: index for temperature
      integer, intent(in) :: indexS !< Input: index for salinity

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

      integer, intent(out) :: err !< Output: Error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iCell, nCellsSolve, k

      real (kind=RKIND), dimension(:), pointer :: temperatureRestore, salinityRestore
      real (kind=RKIND) :: invTemp, invSalinity

      err = 0

      if(.not.restoringOn) return

      nCellsSolve = mesh % nCellsSolve

      temperatureRestore => mesh % temperatureRestore % array
      salinityRestore => mesh % salinityRestore % array

      invTemp = 1.0 / (temperatureTimeScale * 86400.0)
      invSalinity = 1.0 / (salinityTimeScale * 86400.0)

      k = 1  ! restoring only in top layer
      do iCell=1,nCellsSolve
        tend(indexT, k, iCell) = tend(indexT, k, iCell) - layerThickness(k,iCell)*(tracers(indexT, k, iCell) - temperatureRestore(iCell)) * invTemp
        tend(indexS, k, iCell) = tend(indexS, k, iCell) - layerThickness(k,iCell)*(tracers(indexS, k, iCell) - salinityRestore(iCell)) * invSalinity

!       write(6,10) iCell, tracers(indexT, k, iCell), &
!              temperatureRestore(iCell), tracers(indexT, k, iCell), &
!             (tracers(indexT, k, iCell) - temperatureRestore(iCell)) &
!             / (config_restoreT_timescale * 86400.0)

      enddo

   !--------------------------------------------------------------------

   end subroutine ocn_restoring_tend!}}}

!***********************************************************************
!
!  routine ocn_restoring_init
!
!> \brief   Initializes ocean tracer restoring
!> \author  Doug Jacobsen
!> \date    19 September 2011
!> \details 
!>  This routine initializes a variety of quantities related to 
!>  restoring in the ocean. 
!
!-----------------------------------------------------------------------

   subroutine ocn_restoring_init(err)!{{{

      integer, intent(out) :: err !< Output: error flag

      err = 0
      restoringOn = .false.

      if(config_restoreTS) then
          restoringOn = .true.
          temperatureTimeScale = config_restoreT_timescale
          salinityTimeScale = config_restoreS_timescale
      endif

   !--------------------------------------------------------------------

   end subroutine ocn_restoring_init!}}}

!***********************************************************************

end module ocn_restoring

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
