! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_TTD
!
!> \brief MPAS ocean restoring
!> \author Todd Ringler
!> \date   06/08/2015
!> \details
!>  This module contains routines for computing the tracer tendency due to
!>    to transit time distribution
!
!-----------------------------------------------------------------------

module ocn_tracer_TTD

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

   public :: ocn_tracer_TTD_compute, &
             ocn_tracer_TTD_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_TTD_compute
!
!> \brief   computes a tracer tendency to approximate transit time distribution
!> \author  Todd Ringler
!> \date    06/09/2015
!> \details
!>  This routine computes a tracer tendency to approximate transit time distribution
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_TTD_compute(nTracers, nCellsSolve, maxLevelCell, layerThickness, &
                  TTDMask, tracers, err)!{{{

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
         layerThickness, &
         TTDMask

      integer, intent(in) :: nTracers, nCellsSolve

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      ! three dimensional arrays
      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
        tracers

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

      !move to ocean constants
      real (kind=RKIND), parameter :: c0 = 0.0_RKIND
      real (kind=RKIND), parameter :: c1 = 1.0_RKIND

      err = 0

      ! zero tracers at surface to TTDMask at top-most layer
      ! TTDMask should be 1 within region of interest and zero elsewhere
      !$omp do schedule(runtime)
      do iCell = 1, nCellsSolve
         tracers(:, 1, iCell) = TTDMask(:, iCell)
      end do
      !$omp end do

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_TTD_compute!}}}

!***********************************************************************
!
!  routine ocn_tracer_TTD_init
!
!> \brief   Initializes ocean ideal age
!> \author  Todd Ringler
!> \date    06/09/2015
!> \details
!>  This routine initializes fields required for tracer ideal age
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_TTD_init(err)!{{{

      integer, intent(out) :: err !< Output: error flag

      err = 0

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_TTD_init!}}}

!***********************************************************************

end module ocn_tracer_TTD

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
