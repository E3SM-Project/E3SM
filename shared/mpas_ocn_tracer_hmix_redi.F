! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_hmix_redi
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

module ocn_tracer_hmix_redi

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

   public :: ocn_tracer_hmix_redi_tend, &
             ocn_tracer_hmix_redi_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: rediOn
   logical, pointer :: config_disable_redi_horizontal_term1
   logical, pointer :: config_disable_redi_horizontal_term2
   logical, pointer :: config_disable_redi_horizontal_term3
   real (kind=RKIND), pointer :: config_Redi_kappa


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_hmix_redi_tend
!
!> \brief   Computes Laplacian tendency term for horizontal tracer mixing
!> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date    September 2011
!> \details
!>  This routine computes the horizontal mixing tendency for tracers
!>  based on current state using a Laplacian parameterization.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_hmix_redi_tend(meshPool, scratchPool, layerThicknessEdge, zMid, tracers, &
                                        relativeSlopeTopOfEdge, relativeSlopeTapering, relativeSlopeTaperingCell, tend, err)!{{{


      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information
      type (mpas_pool_type), intent(in) :: scratchPool !< Input: Scratch information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThicknessEdge, &!< Input: thickness at edge
         zMid,               &!< Input: Z coordinate at the center of a cell
         relativeSlopeTopOfEdge,    &!< Input: slope of coordinate relative to neutral surface at edges
         relativeSlopeTapering,     &!< Input: tapering of slope of coordinate relative to neutral surface at edges
         relativeSlopeTaperingCell   !< Input: tapering of slope of coordinate relative to neutral surface at cells

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
      integer :: i, k, iTracer, num_tracers, nCells, nEdges
      integer, pointer :: nVertLevels
      integer, dimension(:), pointer :: nCellsArray, nEdgesArray

      integer, dimension(:,:), allocatable :: boundaryMask

      integer, dimension(:), pointer :: maxLevelEdgeTop, nEdgesOnCell, maxLevelCell
      integer, dimension(:,:), pointer :: cellsOnEdge, edgesOnCell, edgeSignOnCell

      real (kind=RKIND) :: invAreaCell1, invAreaCell2, invAreaCell, areaEdge
      real (kind=RKIND) :: tracer_turb_flux, flux, s_tmp, r_tmp, h1, h2, s_tmpU, s_tmpD

      real (kind=RKIND), dimension(:), pointer :: areaCell, dvEdge, dcEdge

      real (kind=RKIND), dimension(:,:), pointer :: gradTracerEdge, gradTracerTopOfEdge, gradHTracerSlopedTopOfCell, &
         dTracerdZTopOfCell, dTracerdZTopOfEdge, areaCellSum

      type (field2DReal), pointer :: gradTracerEdgeField, gradTracerTopOfEdgeField, gradHTracerSlopedTopOfCellField, &
                                     dTracerdZTopOfCellField, dTracerdZTopOfEdgeField, areaCellSumField

      err = 0

      if (.not.rediOn) return

      call mpas_timer_start("tracer redi")

      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      num_tracers = size(tracers, dim=1)

      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)

      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)

      call mpas_pool_get_config(ocnConfigs, 'config_Redi_kappa',config_Redi_kappa)
      call mpas_pool_get_config(ocnConfigs, 'config_disable_redi_horizontal_term1',config_disable_redi_horizontal_term1)
      call mpas_pool_get_config(ocnConfigs, 'config_disable_redi_horizontal_term2',config_disable_redi_horizontal_term2)
      call mpas_pool_get_config(ocnConfigs, 'config_disable_redi_horizontal_term3',config_disable_redi_horizontal_term3)

      !
      ! COMPUTE the extra terms arising due to mismatch between the constant coordinate surfaces and the
      ! isopycnal surfaces.
      !

      call mpas_pool_get_field(scratchPool, 'gradTracerEdge', gradTracerEdgeField)
      call mpas_pool_get_field(scratchPool, 'gradTracerTopOfEdge', gradTracerTopOfEdgeField)
      call mpas_pool_get_field(scratchPool, 'gradHTracerSlopedTopOfCell', gradHTracerSlopedTopOfCellField)
      call mpas_pool_get_field(scratchPool, 'dTracerdZTopOfCell', dTracerdZTopOfCellField)
      call mpas_pool_get_field(scratchPool, 'dTracerdZTopOfEdge', dTracerdZTopOfEdgeField)
      call mpas_pool_get_field(scratchPool, 'areaCellSum', areaCellSumField)

      call mpas_allocate_scratch_field(gradTracerEdgeField, .true., .false.)
      call mpas_allocate_scratch_field(gradTracerTopOfEdgeField, .true., .false.)
      call mpas_allocate_scratch_field(gradHTracerSlopedTopOfCellField, .true., .false.)
      call mpas_allocate_scratch_field(dTracerdZTopOfCellField, .true., .false.)
      call mpas_allocate_scratch_field(dTracerdZTopOfEdgeField, .true., .false.)
      call mpas_allocate_scratch_field(areaCellSumField, .true., .false.)
      call mpas_threading_barrier()

      gradTracerEdge => gradTracerEdgeField % array
      gradTracerTopOfEdge => gradTracerTopOfEdgeField % array
      gradHTracerSlopedTopOfCell => gradHTracerSlopedTopOfCellField % array
      dTracerdZTopOfCell => dTracerdZTopOfCellField % array
      dTracerdZTopOfEdge => dTracerdZTopOfEdgeField % array
      areaCellSum => areaCellSumField % array

      nCells = nCellsArray( size(nCellsArray) )
      nEdges = nEdgesArray( size(nEdgesArray) )

      !$omp do schedule(runtime)
      do iCell = 1, nCells
        gradHTracerSlopedTopOfCell(:, iCell) = 0.0_RKIND
        dTracerdZTopOfCell(:, iCell) = 0.0_RKIND
      end do
      !$omp end do

      !$omp do schedule(runtime)
      do iEdge = 1, nEdges
         gradTracerEdge(:, iEdge) = 0.0_RKIND
         gradTracerTopOfEdge(:, iEdge) = 0.0_RKIND
         dTracerdZTopOfEdge(:, iEdge) = 0.0_RKIND
      end do
      !$omp end do

      ! this is the "standard" del2 term, but forced to use config_redi_kappa
      if(.not.config_disable_redi_horizontal_term1) then

         nCells = nCellsArray( 1 )
         !$omp do schedule(runtime) private(invAreaCell, i, iEdge, cell1, cell2, r_tmp, k, s_tmp, iTracer, tracer_turb_flux, flux)
         do iCell = 1, nCells
            invAreaCell = 1.0_RKIND / areaCell(iCell)
            do i = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i, iCell)
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)

               r_tmp = config_redi_kappa * dvEdge(iEdge) / dcEdge(iEdge)

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

      endif

      ! Compute vertical derivative of tracers at cell center and top of layer
      do iTracer = 1, num_tracers
         ! Sync threads before starting on tracers
         call mpas_threading_barrier()

         nCells = nCellsArray( 2 )
         !$omp do schedule(runtime) private(k)
         do iCell = 1, nCells
            do k = 2, maxLevelCell(iCell)
               dTracerdZTopOfCell(k,iCell) = (tracers(iTracer,k-1,iCell) - tracers(iTracer,k,iCell)) &
                                           / (zMid(k-1,iCell) - zMid(k,iCell))
            end do

            ! Approximation of dTracerdZTopOfCell on the top and bottom interfaces through the idea of having
            ! ghost cells above the top and below the bottom layers of the same depths and tracer density.
            ! Essentially, this enforces the boundary condition (d tracer)/dz = 0 at the top and bottom.
            dTracerdZTopOfCell(1,iCell) = 0.0_RKIND
            dTracerdZTopOfCell(maxLevelCell(iCell)+1,iCell) = 0.0_RKIND
         end do
         !$omp end do

         nEdges = nEdgesArray( 2 )
         ! Compute tracer gradient (gradTracerEdge) along the constant coordinate surface.
         ! The computed variables lives at edge and mid-layer depth
         !$omp do schedule(runtime) private(cell1, cell2, k)
         do iEdge = 1, nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            do k=1,maxLevelEdgeTop(iEdge)
               gradTracerEdge(k,iEdge) = (tracers(iTracer,k,cell2) - tracers(iTracer,k,cell1)) / dcEdge(iEdge)
            end do
         end do
         !$omp end do

         nEdges = nEdgesArray( 2 )
         ! Interpolate dTracerdZTopOfCell to edge and top of layer
         !$omp do schedule(runtime) private(cell1, cell2, k)
         do iEdge = 1, nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            do k = 1, maxLevelEdgeTop(iEdge)
               dTracerdZTopOfEdge(k,iEdge) = 0.5_RKIND * (dTracerdZTopOfCell(k,cell1) + dTracerdZTopOfCell(k,cell2))
            end do
            dTracerdZTopOfEdge(maxLevelEdgeTop(iEdge)+1,iEdge) = 0.0_RKIND
         end do
         !$omp end do

         nEdges = nEdgesArray( 2 )
         ! Interpolate gradTracerEdge to edge and top of layer
         !$omp do schedule(runtime) private(k, h1, h2)
         do iEdge = 1, nEdges
            do k = 2, maxLevelEdgeTop(iEdge)
               h1 = layerThicknessEdge(k-1,iEdge)
               h2 = layerThicknessEdge(k,iEdge)

               ! Using second-order interpolation below
               gradTracerTopOfEdge(k,iEdge) = (h2 * gradTracerEdge(k-1,iEdge) + h1 * gradTracerEdge(k,iEdge)) / (h1 + h2)
            end do

            ! Approximation of values on the top and bottom interfaces through the idea of having ghost cells above
            ! the top and below the bottom layers of the same depths and tracer concentration.
            gradTracerTopOfEdge(1,iEdge) = gradTracerEdge(2,iEdge)
            gradTracerTopOfEdge(maxLevelEdgeTop(iEdge)+1,iEdge) = gradTracerEdge(max(maxLevelEdgeTop(iEdge),1),iEdge)
         end do
         !$omp end do

         ! Compute \nabla\cdot(relativeSlope d\phi/dz)
         if(.not.config_disable_redi_horizontal_term2) then
            nCells = nCellsArray( 1 )
            !$omp do schedule(runtime) private(invAreaCell, i, iEdge, k, s_tmpU, s_tmpD, flux)
            do iCell = 1, nCells
               invAreaCell = 1.0_RKIND / areaCell(iCell)
               do i = 1, nEdgesOnCell(iCell)
                  iEdge = edgesOnCell(i, iCell)
                  do k = 1, maxLevelEdgeTop(iEdge)
                     s_tmpU = relativeSlopeTapering(k, iEdge) * relativeSlopeTopOfEdge(k, iEdge) * dTracerdZTopOfEdge(k, iEdge)
                     s_tmpD = relativeSlopeTapering(k+1, iEdge) * relativeSlopeTopOfEdge(k+1, iEdge) &
                            * dTracerdZTopOfEdge(k+1, iEdge)

                     flux = 0.5 * dvEdge(iEdge) * ( s_tmpU + s_tmpD )
                     flux = flux * layerThicknessEdge(k, iEdge)
                     tend(iTracer, k, iCell) = tend(iTracer, k, iCell) + edgeSignOnCell(i, iCell) * config_Redi_kappa * flux &
                                             * invAreaCell
                  end do
               end do
            end do
            !$omp end do
         endif

         ! Compute dz * d(relativeSlope\cdot\nabla\phi)/dz  (so the dz cancel out)

         ! Compute relativeSlope\cdot\nabla\phi (variable gradHTracerSlopedTopOfCell) at non-boundary edges

         nCells = nCellsArray( 1 )
         !$omp do schedule(runtime) private(i, iedge, areaEdge, k, r_tmp)
         do iCell = 1, nCells
            areaCellSum(:, iCell) = 1.0e-34_RKIND
            gradHTracerSlopedTopOfCell(:, iCell) = 0.0_RKIND
            do i = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i, iCell)
               areaEdge = 0.5_RKIND * dcEdge(iEdge) * dvEdge(iEdge)
               do k = 1, maxLevelEdgeTop(iEdge)
                  r_tmp = areaEdge * relativeSlopeTapering(k, iEdge) * relativeSlopeTopOfEdge(k,iEdge) &
                        * gradTracerTopOfEdge(k,iEdge)
                  gradHTracerSlopedTopOfCell(k, iCell) = gradHTracerSlopedTopOfCell(k, iCell) + r_tmp
                  areaCellSum(k, iCell) = areaCellSum(k, iCell) + areaEdge
               end do
            end do
         end do
         !$omp end do

         nCells = nCellsArray( 1 )
         !$omp do schedule(runtime) private(k)
         do iCell=1,nCells
            do k = 1, maxLevelCell(iCell)
               gradHTracerSlopedTopOfCell(k,iCell) = gradHTracerSlopedTopOfCell(k,iCell)/areaCellSum(k,iCell)
            end do
         end do
         !$omp end do

         if(.not.config_disable_redi_horizontal_term3) then
            nCells = nCellsArray( 1 )
            !$omp do schedule(runtime) private(k, s_tmp)
            do iCell = 1, nCells
               ! impose no-flux boundary conditions at top and bottom of column
               gradHTracerSlopedTopOfCell(1,iCell) = 0.0_RKIND
               gradHTracerSlopedTopOfCell(maxLevelCell(iCell)+1,iCell) = 0.0_RKIND
               do k = 1, maxLevelCell(iCell)
                  tend(iTracer,k,iCell) = tend(iTracer,k,iCell) + config_Redi_kappa * &
                      (gradHTracerSlopedTopOfCell(k,iCell) - gradHTracerSlopedTopOfCell(k+1,iCell))
               end do
            end do
            !$omp end do
         endif

      end do  ! iTracer

      call mpas_threading_barrier()
      call mpas_deallocate_scratch_field(gradTracerEdgeField, .true.)
      call mpas_deallocate_scratch_field(gradTracerTopOfEdgeField, .true.)
      call mpas_deallocate_scratch_field(gradHTracerSlopedTopOfCellField, .true.)
      call mpas_deallocate_scratch_field(dTracerdZTopOfCellField, .true.)
      call mpas_deallocate_scratch_field(dTracerdZTopOfEdgeField, .true.)

      call mpas_timer_stop("tracer redi")

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_hmix_redi_tend!}}}

!***********************************************************************
!
!  routine ocn_tracer_hmix_redi_init
!
!> \brief   Initializes ocean tracer horizontal mixing quantities
!> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date    September 2011
!> \details
!>  This routine initializes a variety of quantities related to
!>  Laplacian horizontal velocity mixing in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_hmix_redi_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      logical, pointer :: config_use_standardGM

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_use_standardGM', config_use_standardGM)

      rediOn = .false.

      if ( config_use_standardGM ) then
          rediOn = .true.
      endif


   !--------------------------------------------------------------------

   end subroutine ocn_tracer_hmix_redi_init!}}}

!***********************************************************************

end module ocn_tracer_hmix_redi

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
