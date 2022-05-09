! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_hmix_del4
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

module ocn_tracer_hmix_del4

   use mpas_timer
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_threading
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

   public :: ocn_tracer_hmix_del4_tend, &
             ocn_tracer_hmix_del4_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: del4On

   real (kind=RKIND) :: eddyDiff4

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_hmix_del4_tend
!
!> \brief   Computes biharmonic tendency term for horizontal tracer mixing
!> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date    September 2011
!> \details
!>  This routine computes the horizontal mixing tendency for tracers
!>  based on current state using a biharmonic parameterization.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_hmix_del4_tend(meshPool, layerThickEdgeMean, tracers, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickEdgeMean    !< Input: mean thickness at edge

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

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

      integer :: iEdge, num_tracers, nCells
      integer :: iTracer, k, iCell, cell1, cell2, i
      integer, pointer :: nVertLevels
      integer, dimension(:), pointer :: nCellsArray, nEdgesArray

      integer, dimension(:), pointer :: minLevelEdgeBot, maxLevelEdgeTop, minLevelCell, maxLevelCell, nEdgesOnCell
      integer, dimension(:,:), pointer :: edgeMask, cellsOnEdge, edgesOnCell, edgeSignOnCell

      real (kind=RKIND) :: invAreaCell1, tracer_turb_flux, flux, invdcEdge, r_tmp1, r_tmp2

      real (kind=RKIND), dimension(:), pointer :: dcEdge, dvEdge, areaCell, meshScalingDel4

      ! Scratch Arrays
      real (kind=RKIND), dimension(:,:,:), allocatable :: delsq_tracer

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing tendencies
      ! note that the user can choose multiple options and the
      !   tendencies will be added together
      !
      !-----------------------------------------------------------------

      err = 0

      if ( .not. del4On ) return

      call mpas_timer_start("tracer del4")

      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)
      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      num_tracers = size(tracers, dim=1)

      call mpas_pool_get_array(meshPool, 'minLevelEdgeBot', minLevelEdgeBot)
      call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)

      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'meshScalingDel4', meshScalingDel4)

      call mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)

      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)

      ! Need 1 halo around owned cells
      nCells = nCellsArray( 2 )

      allocate(delsq_tracer(num_tracers, nVertLevels, nCells))

      ! first del2: div(h \nabla \phi) at cell center
      !$omp parallel
      !$omp do schedule(runtime) private(invAreaCell1, i, iEdge, invdcEdge, cell1, cell2, k, iTracer, r_tmp1, r_tmp2)
      do iCell = 1, nCells
        delsq_tracer(:, :, iCell) = 0.0_RKIND
        invAreaCell1 = 1.0_RKIND / areaCell(iCell)
        do i = 1, nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i, iCell)
          invdcEdge = dvEdge(iEdge) / dcEdge(iEdge)
          cell1 = cellsOnEdge(1,iEdge)
          cell2 = cellsOnEdge(2,iEdge)

          do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
            do iTracer = 1, num_tracers * edgeMask(k, iEdge)

              r_tmp1 = invdcEdge * layerThickEdgeMean(k, iEdge) * tracers(iTracer, k, cell1)
              r_tmp2 = invdcEdge * layerThickEdgeMean(k, iEdge) * tracers(iTracer, k, cell2)

              delsq_tracer(iTracer, k, iCell) = delsq_tracer(iTracer, k, iCell) - edgeSignOnCell(i, iCell) &
                                              * (r_tmp2 - r_tmp1) * invAreaCell1
            end do
          end do
        end do
      end do
      !$omp end do
      !$omp end parallel

      ! Only need tendency on owned cells
      nCells = nCellsArray( 1 )

      ! second del2: div(h \nabla [delsq_tracer]) at cell center
      !$omp parallel
      !$omp do schedule(runtime) private(invAreaCell1, i, iEdge, cell1, cell2, invdcEdge, k, iTracer, tracer_turb_flux, flux)
      do iCell = 1, nCells
        invAreaCell1 = 1.0_RKIND / areaCell(iCell)
        do i = 1, nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i, iCell)
          cell1 = cellsOnEdge(1, iEdge)
          cell2 = cellsOnedge(2, iEdge)

          invdcEdge = meshScalingDel4(iEdge) * dvEdge(iEdge) * eddyDiff4 / dcEdge(iEdge)

          do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
            do iTracer = 1, num_tracers * edgeMask(k, iEdge)
              tracer_turb_flux = (delsq_tracer(iTracer, k, cell2) - delsq_tracer(iTracer, k, cell1))

              flux = tracer_turb_flux * invdcEdge

              tend(iTracer, k, iCell) = tend(iTracer, k, iCell) + edgeSignOnCell(i, iCell) * flux * invAreaCell1
            end do
          end do
        end do
      end do
      !$omp end do
      !$omp end parallel

      deallocate(delsq_tracer)

      call mpas_timer_stop("tracer del4")

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_hmix_del4_tend!}}}

!***********************************************************************
!
!  routine ocn_tracer_hmix_del4_init
!
!> \brief   Initializes ocean tracer horizontal mixing quantities
!> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date    September 2011
!> \details
!>  This routine initializes a variety of quantities related to
!>  biharmonic horizontal velocity mixing in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_hmix_del4_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      err = 0

      del4on = .false.

      if ( config_tracer_del4 > 0.0_RKIND ) then
          del4On = .true.
          eddyDiff4 = config_tracer_del4
      endif

      if ( .not. config_use_tracer_del4 ) del4on = .false.

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_hmix_del4_init!}}}

!***********************************************************************

end module ocn_tracer_hmix_del4

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
