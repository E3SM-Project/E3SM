! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_exponential_decay
!
!> \brief MPAS ocean exponential decay
!> \author Todd Ringler
!> \date   06/08/2015
!> \details
!>  This module contains routines for computing tracer forcing due to exponential decay
!
!-----------------------------------------------------------------------

module ocn_tracer_exponential_decay

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

   public :: ocn_tracer_exponential_decay_compute, &
             ocn_tracer_exponential_decay_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_exponential_decay_compute
!
!> \brief   computes a tracer tendency due to exponential decay
!> \author  Todd Ringler
!> \date    06/09/2015
!> \details
!>  This routine computes a tracer tendency due to exponential decay
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_exponential_decay_compute(nTracers, nCellsSolve, maxLevelCell, layerThickness, tracers, & !{{{
                                                   tracersExponentialDecayRate, tracer_tend, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      ! one dimensional arrays
      integer, dimension(:), intent(in) :: &
         maxLevelCell

      real (kind=RKIND), dimension(:), intent(in) :: &
         tracersExponentialDecayRate

      ! two dimensional arrays
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness

      ! three dimensional arrays
      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         tracers

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
               tracer_tend(iTracer,iLevel,iCell) =   tracer_tend(iTracer,iLevel,iCell)   &
                                                    - (   layerThickness(iLevel,iCell)    &
                                                        * tracers(iTracer,iLevel,iCell)   &
                                                        * exp(-tracersExponentialDecayRate(iTracer)) )
            enddo
         enddo
      enddo
      !$omp end do

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_exponential_decay_compute!}}}

!***********************************************************************
!
!  routine ocn_tracer_exponential_decay_init
!
!> \brief   Initializes ocean surface restoring
!> \author  Todd Ringler
!> \date    06/09/2015
!> \details
!>  This routine initializes fields required for tracer surface flux restoring
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_exponential_decay_init(err)!{{{

      integer, intent(out) :: err !< Output: error flag

      err = 0

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_exponential_decay_init!}}}

!***********************************************************************

end module ocn_tracer_exponential_decay

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
