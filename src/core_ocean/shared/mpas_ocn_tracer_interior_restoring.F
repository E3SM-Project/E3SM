! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_interior_restoring
!
!> \brief MPAS ocean restoring
!> \author Todd Ringler
!> \date   06/08/2015
!> \details
!>  This module contains routines for computing the tracer tendency due to restoring
!
!-----------------------------------------------------------------------

module ocn_tracer_interior_restoring

   use mpas_kind_types
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

   public :: ocn_tracer_interior_restoring_compute, &
             ocn_tracer_interior_restoring_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_interior_restoring_compute
!
!> \brief   computes a tracer tendency due to interior restoring
!> \author  Todd Ringler
!> \date    06/09/2015
!> \details
!>  This routine computes a tracer tendency due to interior restoring
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_interior_restoring_compute(nTracers, nCellsSolve, maxLevelCell, layerThickness, &
                  tracers, tracersInteriorRestoringRate, tracersInteriorRestoringValue, tracer_tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      ! one dimensional arrays
      integer, dimension(:), intent(in) :: &
         maxLevelCell

      ! two dimensional arrays
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness

      ! three dimensional arrays
      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         tracers, &
         tracersInteriorRestoringRate, &
         tracersInteriorRestoringValue

      ! scalars
      integer, intent(in) :: nTracers, nCellsSolve

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
        tracer_tend

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

      integer :: iCell, iLevel, iTracer

      err = 0

      !$omp do schedule(runtime) private(iLevel, iTracer)
      do iCell=1,nCellsSolve
        do iLevel=1,maxLevelCell(iCell)
          do iTracer=1,nTracers
             tracer_tend(iTracer, iLevel, iCell) = tracer_tend(iTracer, iLevel, iCell)  - layerThickness(iLevel,iCell) &
                                                 * ( tracers(iTracer, iLevel, iCell) &
                                                 - tracersInteriorRestoringValue(iTracer, iLevel, iCell) )  &
                                                 * tracersInteriorRestoringRate(iTracer, iLevel, iCell)
          enddo
        enddo
      enddo
      !$omp end do

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_interior_restoring_compute!}}}

!***********************************************************************
!
!  routine ocn_tracer_interior_restoring_init
!
!> \brief   Initializes ocean interior restoring
!> \author  Todd Ringler
!> \date    06/09/2015
!> \details
!>  This routine initializes fields required for tracer interior restoring
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_interior_restoring_init(err)!{{{

      integer, intent(out) :: err !< Output: error flag

      err = 0

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_interior_restoring_init!}}}

!***********************************************************************

end module ocn_tracer_interior_restoring

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
