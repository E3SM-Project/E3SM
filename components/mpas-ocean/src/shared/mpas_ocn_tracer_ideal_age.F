! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_ideal_age
!
!> \brief MPAS ideal age tracer routines
!> \author Mathew Maltrud
!> \date   08/19/2021
!> \details
!>  This module contains routines for computing the tracer tendency due to restoring
!
!-----------------------------------------------------------------------

module ocn_tracer_ideal_age

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use ocn_constants
   use ocn_config

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
!> \author  Mathew Maltrud
!> \date    08/19/2021
!> \details
!>  This routine computes a tracer tendency to approximate ideal age
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_ideal_age_compute(nTracers, nCellsSolve, minLevelCell, maxLevelCell, layerThickness, &
                  idealAgeMask, tracers, tracer_tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      ! one dimensional arrays
      integer, dimension(:), intent(in) :: &
         minLevelCell, maxLevelCell

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

      err = 0

      ! add a tendency increment equivalent to "dt" to entire domain below top level
      !$omp parallel
      !$omp do schedule(runtime) private(iLevel, iTracer)
      do iCell=1,nCellsSolve
        do iLevel=minLevelCell(iCell) + 1, maxLevelCell(iCell)
          do iTracer=1,nTracers
             tracer_tend(iTracer, iLevel, iCell) = tracer_tend(iTracer, iLevel, iCell) +  &
                layerThickness(iLevel,iCell)
          enddo
        enddo
      enddo
      !$omp end do
      !$omp end parallel

      ! zero surface tendency where mask == 0
      !$omp parallel
      !$omp do schedule(runtime) private(iLevel, iTracer)
      do iCell=1,nCellsSolve
        iLevel = minLevelCell(iCell)
        do iTracer=1,nTracers
           tracer_tend(iTracer, iLevel, iCell) = tracer_tend(iTracer, iLevel, iCell) +  &
              idealAgeMask(iTracer, iCell)*layerThickness(iLevel,iCell)
        enddo
      enddo
      !$omp end do
      !$omp end parallel

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_ideal_age_compute!}}}

!***********************************************************************
!
!  routine ocn_tracer_ideal_age_init
!
!> \brief   Initializes ocean ideal age
!> \author  Mathew Maltrud
!> \date    08/19/2021
!> \details
!>  This routine initializes fields required for tracer ideal age
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_ideal_age_init(domain,err)!{{{

      type (domain_type), intent(inout) :: domain !< Input/Output: domain information

      integer, intent(out) :: err !< Output: error flag

      type (mpas_pool_type), pointer :: forcingPool
      type (mpas_pool_type), pointer :: tracersIdealAgeFieldsPool

      integer, dimension(:), pointer :: &
         landIceMask

      real (kind=RKIND), dimension(:,:), pointer :: &
         idealAgeMask

      err = 0

      if (.not. config_use_idealAgeTracers) return

      call mpas_pool_get_subpool(domain % blocklist % structs, 'forcing', forcingPool)
      call mpas_pool_get_subpool(forcingPool, &
                                 'tracersIdealAgeFields', &
                                  tracersIdealAgeFieldsPool)
      call mpas_pool_get_array(tracersIdealAgeFieldsPool, &
                               'idealAgeTracersIdealAgeMask', &
                                idealAgeMask)

      call mpas_pool_get_array(forcingPool, 'landIceMask', landIceMask)

      if (config_use_idealAgeTracers_idealAge_forcing) then
         ! set mask = 0 for open ocean
         ! set mask = 1 under land ice
         idealAgeMask = 0.0_RKIND
         if (associated(landIceMask)) then
            where (landIceMask /= 0) idealAgeMask(1,:) = 1.0_RKIND
         end if
      end if

   end subroutine ocn_tracer_ideal_age_init!}}}

!***********************************************************************

end module ocn_tracer_ideal_age

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
