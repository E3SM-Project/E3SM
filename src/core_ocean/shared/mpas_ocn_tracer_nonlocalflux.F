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
!> \version SVN:$Id$
!> \details 
!>  This routine computes the tendency for tracers based the vertical divergence of non-local fluxes.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_nonlocalflux_tend(mesh, layerThickness, vertNonLocalFlux, surfaceTracerFlux, tend, err)!{{{
      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mesh_type), intent(in) :: &
         mesh          !< Input: mesh information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
        layerThickness !< Input: Layer thickness

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

      integer :: iCell, nCells, k, iTracer, nTracers, nVertLevels
      integer, dimension(:), pointer :: maxLevelCell
      integer, dimension(:,:), pointer :: cellMask

      err = 0

      if (.not. nonLocalFluxOn) return

      nCells = mesh % nCells
      nVertLevels = mesh % nVertLevels
      nTracers = size(tend, dim=1)

      maxLevelCell => mesh % maxLevelCell % array
      cellMask => mesh % cellMask % array

      do iCell = 1, nCells
        do k = 1, maxLevelCell(iCell)

          ! NOTE: at the moment, all tracers are based on the flux-profile used for temperature, i.e. vertNonLocalFlux(1,:,:)
          do iTracer = 1, nTracers
            tend(iTracer, k, iCell) = tend(iTracer, k, iCell) + cellMask(k, icell) * surfaceTracerFlux(iTracer, iCell) * vertNonLocalFlux(1, k, iCell)
          end do
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

      err = 0

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
