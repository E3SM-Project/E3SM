! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_ideal_age
!
!> \brief MPAS ocean restoring
!> \author Todd Ringler
!> \date   06/08/2015
!> \details
!>  This module contains routines for computing the tracer tendency due to restoring
!
!-----------------------------------------------------------------------

module ocn_tracer_ideal_age

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

   public :: ocn_tracer_ideal_age_compute, &
             ocn_tracer_ideal_age_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_ideal_age_compute
!
!> \brief   computes a tracer tendency to approximate ideal age
!> \author  Todd Ringler
!> \date    06/09/2015
!> \details
!>  This routine computes a tracer tendency to approximate ideal age
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_ideal_age_compute(nTracers, nCellsSolve, maxLevelCell, layerThickness, &
                  idealAgeMask, tracers, tracer_tend, err)!{{{

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
         idealAgeMask

      integer, intent(in) :: nTracers, nCellsSolve

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      ! three dimensional arrays
      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
        tracers, &
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

      !move to ocean constants
      real (kind=RKIND), parameter :: c0 = 0.0_RKIND
      real (kind=RKIND), parameter :: c1 = 1.0_RKIND

      err = 0

      !$omp do schedule(runtime) private(iLevel, iTracer)
      do iCell=1,nCellsSolve
        do iLevel=1,maxLevelCell(iCell)
          do iTracer=1,nTracers
             ! zero tracers at surface to zero where idealAgeMask == zero
             ! idealAgeMask should be equal to 1.0 elsewhere
             tracers(iTracer, iLevel, iCell) = idealAgeMask(iTracer, iCell) * tracers(iTracer, iLevel, iCell)

             ! add a tendency increment equivalent to "dt" to entire domain
             tracer_tend(iTracer, iLevel, iCell) =   tracer_tend(iTracer, iLevel, iCell) +  &
                                                     layerThickness(iLevel,iCell) * c1
          enddo
        enddo
      enddo
      !$omp end do

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_ideal_age_compute!}}}

!***********************************************************************
!
!  routine ocn_tracer_ideal_age_init
!
!> \brief   Initializes ocean ideal age
!> \author  Todd Ringler
!> \date    06/09/2015
!> \details
!>  This routine initializes fields required for tracer ideal age
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_ideal_age_init(err)!{{{

      integer, intent(out) :: err !< Output: error flag

      err = 0

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_ideal_age_init!}}}

!***********************************************************************

end module ocn_tracer_ideal_age

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
