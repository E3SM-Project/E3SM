! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_nonlocalflux
!
!> \brief MPAS ocean tracer non-local flux
!> \author Todd Ringler
!> \date   11/25/13
!> \version SVN:$Id:$
!> \details
!>  This module contains the routine for computing 
!>  tracer tendencies due to non-local vertical fluxes computed in CVMix KPP
!
!-----------------------------------------------------------------------

module ocn_tracer_nonlocalflux

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

   public :: ocn_tracer_nonlocalflux_tend, &
             ocn_tracer_nonlocalflux_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: nonLocalFluxOn

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_nonlocalflux_tend
!
!> \brief   Computes tendency term due to non-local flux transport
!> \author  Todd Ringler
!> \date    11/25/13
!> \details 
!>  This routine computes the tendency for tracers based the vertical divergence of non-local fluxes.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_nonlocalflux_tend(meshPool, vertNonLocalFlux, surfaceTracerFlux, tend, err)!{{{
      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
        surfaceTracerFlux !< Input: surface tracer fluxes

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
        vertNonLocalFlux !< Input: non-local flux of tracers defined at layer interfaces

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

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iCell, k, iTracer, nTracers
      integer, pointer :: nCells, nVertLevels
      integer, dimension(:), pointer :: maxLevelCell
      integer, dimension(:,:), pointer :: cellMask
      real (kind=RKIND) :: fluxTopOfCell, fluxBottomOfCell

      err = 0

      if (.not. nonLocalFluxOn) return

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      nTracers = size(tend, dim=1)

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'cellMask', cellMask)

      do iCell = 1, nCells
        do k = 2, maxLevelCell(iCell)-1

          ! NOTE: at the moment, all tracers are based on the flux-profile used for temperature, i.e. vertNonLocalFlux(1,:,:)
          do iTracer = 1, nTracers
            fluxTopOfCell = surfaceTracerFlux(iTracer, iCell) * vertNonLocalFlux(1, k, iCell)
            fluxBottomOfCell = surfaceTracerFlux(iTracer, iCell) * vertNonLocalFlux(1, k+1, iCell)
            tend(iTracer, k, iCell) = tend(iTracer, k, iCell) + cellMask(k, icell) * (fluxTopOfCell-fluxBottomOfCell)
          end do
        end do

        ! enforce boundary conditions at bottom of column
        k = maxLevelCell(iCell)
        do iTracer = 1, nTracers
          fluxTopOfCell = surfaceTracerFlux(iTracer, iCell) * vertNonLocalFlux(1, k, iCell)
          fluxBottomOfCell = 0.0
          tend(iTracer, k, iCell) = tend(iTracer, k, iCell) + cellMask(k, icell) * (fluxTopOfCell-fluxBottomOfCell)
        end do

        ! enforce boundary conditions at top of column
        k = 1
        do iTracer = 1, nTracers
          fluxTopOfCell = 0.0
          fluxBottomOfCell = surfaceTracerFlux(iTracer, iCell) * vertNonLocalFlux(1, k+1, iCell)
          tend(iTracer, k, iCell) = tend(iTracer, k, iCell) + cellMask(k, icell) * (fluxTopOfCell-fluxBottomOfCell)
        end do

      end do

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_nonlocalflux_tend!}}}

!***********************************************************************
!
!  routine ocn_tracer_nonlocalflux_init
!
!> \brief   Initializes ocean tracer nonlocal flux computation
!> \author  Todd Ringler
!> \date    11/25/13
!> \version SVN:$Id$
!> \details 
!>  This routine initializes quantities related to nonlocal flux computation
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_nonlocalflux_init(err)!{{{

   !--------------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag
      logical, pointer :: config_disable_tr_nonlocalflux, config_use_cvmix_kpp

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_disable_tr_nonlocalflux', config_disable_tr_nonlocalflux)
      call mpas_pool_get_config(ocnConfigs, 'config_use_cvmix_kpp', config_use_cvmix_kpp)

      nonLocalFluxOn = .true.

      if (config_disable_tr_nonlocalflux) then
         nonLocalFluxOn = .false.
      end if

      if (.not.config_use_cvmix_kpp) then
         nonLocalFluxOn = .false.
      end if

   end subroutine ocn_tracer_nonlocalflux_init!}}}

!***********************************************************************

end module ocn_tracer_nonlocalflux

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
