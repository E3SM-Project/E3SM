! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_surface_restoring
!
!> \brief MPAS ocean restoring
!> \author Todd Ringler
!> \date   06/08/2015
!> \details
!>  This module contains routines for computing the surface tracer flux due to restoring
!
!-----------------------------------------------------------------------

module ocn_tracer_surface_restoring

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

   public :: ocn_tracer_surface_restoring_compute, &
             ocn_tracer_surface_restoring_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_surface_restoring_compute
!
!> \brief   computes a surface tracer flux due to surface restoring
!> \author  Todd Ringler
!> \date    06/09/2015
!> \details 
!>  This routine computes a surface tracer flux due to surface restoring
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_surface_restoring_compute(nTracers, nCellsSolve, tracers, pistonVelocity, tracersSurfaceRestoringValue, tracersSurfaceFlux, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      ! scalars
      integer, intent(in) ::  &
          nTracers,           &
          nCellsSolve

      ! three dimensional arrays
      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         tracers

      ! two dimensional ararys
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         pistonVelocity,     &
         tracersSurfaceRestoringValue

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
        tracersSurfaceFlux 

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

      iLevel = 1  ! base surface flux restoring on tracer fields in the top layer
      do iCell=1,nCellsSolve
        do iTracer=1,nTracers
          tracersSurfaceFlux(iTracer, iCell) =   tracersSurfaceFlux(iTracer, iCell) - &
                                                 pistonVelocity(iTracer,iCell) *      &
                                                 (tracers(iTracer, iLevel, iCell) - tracersSurfaceRestoringValue(iTracer,iCell))
        enddo
      enddo

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_surface_restoring_compute!}}}

!***********************************************************************
!
!  routine ocn_tracer_surface_restoring_init
!
!> \brief   Initializes ocean surface restoring
!> \author  Todd Ringler
!> \date    06/09/2015
!> \details 
!>  This routine initializes fields required for tracer surface flux restoring
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_surface_restoring_init(err)!{{{

      integer, intent(out) :: err !< Output: error flag

      err = 0

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_surface_restoring_init!}}}

!***********************************************************************

end module ocn_tracer_surface_restoring

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
