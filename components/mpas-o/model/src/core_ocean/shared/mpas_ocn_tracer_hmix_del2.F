! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_hmix_del2
!
!> \brief MPAS ocean horizontal tracer mixing driver
!> \author Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date   September 2011
!> \details
!>  This module contains the main driver routine for computing
!>  horizontal mixing tendencies.
!>
!>  It provides an init and a tend function. Each are described below.
!
!-----------------------------------------------------------------------

module ocn_tracer_hmix_del2

   use mpas_timer
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_threading

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

   public :: ocn_tracer_hmix_del2_tend, &
             ocn_tracer_hmix_del2_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: del2On
   real (kind=RKIND) :: eddyDiff2


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_hmix_del2_tend
!
!> \brief   Computes Laplacian tendency term for horizontal tracer mixing
!> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date    September 2011
!> \details
!>  This routine computes the horizontal mixing tendency for tracers
!>  based on current state using a Laplacian parameterization.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_hmix_del2_tend(meshPool, layerThicknessEdge, tracers, tend, err)!{{{
      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThicknessEdge !< Input: thickness at edges

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
        tracers !< Input: tracer quantities

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

      integer :: iCell, iEdge, cell1, cell2
      integer :: i, k, iTracer, num_tracers, nCells
      integer, pointer :: nVertLevels
      integer, dimension(:), pointer :: nCellsArray

      integer, dimension(:), pointer :: maxLevelEdgeTop, nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge, edgesOnCell, edgeSignOnCell

      real (kind=RKIND) :: invAreaCell
      real (kind=RKIND) :: tracer_turb_flux, flux, r_tmp

      real (kind=RKIND), dimension(:), pointer :: areaCell, dvEdge, dcEdge
      real (kind=RKIND), dimension(:), pointer :: meshScalingDel2

      err = 0

      if (.not.del2On) return

      call mpas_timer_start("tracer del2")

      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      num_tracers = size(tracers, dim=1)

      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'meshScalingDel2', meshScalingDel2)

      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)

      nCells = nCellsArray( 1 )

      !
      ! compute a boundary mask to enforce insulating boundary conditions in the horizontal
      !
      !$omp do schedule(runtime) private(invAreaCell, i, iEdge, cell1, cell2, r_tmp, k, iTracer, tracer_turb_flux, flux)
      do iCell = 1, nCells
        invAreaCell = 1.0_RKIND / areaCell(iCell)
        do i = 1, nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i, iCell)
          cell1 = cellsOnEdge(1,iEdge)
          cell2 = cellsOnEdge(2,iEdge)

          r_tmp = meshScalingDel2(iEdge) * eddyDiff2 * dvEdge(iEdge) / dcEdge(iEdge)

          do k = 1, maxLevelEdgeTop(iEdge)
            do iTracer = 1, num_tracers
              ! \kappa_2 \nabla \phi on edge
              tracer_turb_flux = tracers(iTracer, k, cell2) - tracers(iTracer, k, cell1)

              ! div(h \kappa_2 \nabla \phi) at cell center
              flux = layerThicknessEdge(k, iEdge) * tracer_turb_flux * r_tmp

              tend(iTracer, k, iCell) = tend(iTracer, k, iCell) - edgeSignOnCell(i, iCell) * flux * invAreaCell
            end do
          end do

        end do
      end do
      !$omp end do

      call mpas_timer_stop("tracer del2")

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_hmix_del2_tend!}}}

!***********************************************************************
!
!  routine ocn_tracer_hmix_del2_init
!
!> \brief   Initializes ocean tracer horizontal mixing quantities
!> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date    September 2011
!> \details
!>  This routine initializes a variety of quantities related to
!>  Laplacian horizontal velocity mixing in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_hmix_del2_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      logical, pointer :: config_use_tracer_del2
      real (kind=RKIND), pointer :: config_tracer_del2

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_use_tracer_del2', config_use_tracer_del2)
      call mpas_pool_get_config(ocnConfigs, 'config_tracer_del2', config_tracer_del2)

      del2on = .false.

      if ( config_use_tracer_del2 ) then
         if ( config_tracer_del2 > 0.0_RKIND ) then
            del2On = .true.
            eddyDiff2 = config_tracer_del2
         endif
      endif

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_hmix_del2_init!}}}

!***********************************************************************

end module ocn_tracer_hmix_del2

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
