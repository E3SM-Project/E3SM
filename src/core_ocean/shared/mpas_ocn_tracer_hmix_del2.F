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

   use mpas_grid_types

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
   logical, pointer :: config_use_standardGM
   logical, pointer :: config_disable_redi_horizontal_term1
   logical, pointer :: config_disable_redi_horizontal_term2
   logical, pointer :: config_disable_redi_horizontal_term3
   real (kind=RKIND) :: eddyDiff2
   real (kind=RKIND), pointer :: config_Redi_kappa


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

   subroutine ocn_tracer_hmix_del2_tend(meshPool, scratchPool, layerThickness, layerThicknessEdge, zMid, tracers, &
                                        relativeSlopeTopOfEdge, relativeSlopeTapering, relativeSlopeTaperingCell, tend, err)!{{{


      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information
      type (mpas_pool_type), intent(in) :: scratchPool !< Input: Scratch information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness,     &!< Input: thickness at cell centers
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
      integer :: i, k, iTracer, num_tracers
      integer, pointer :: nCells, nVertLevels, nEdges

      integer, dimension(:,:), allocatable :: boundaryMask

      integer, dimension(:), pointer :: maxLevelEdgeTop, nEdgesOnCell, maxLevelCell
      integer, dimension(:,:), pointer :: cellsOnEdge, edgesOnCell, edgeSignOnCell

      real (kind=RKIND) :: invAreaCell1, invAreaCell2, invAreaCell, areaEdge
      real (kind=RKIND) :: tracer_turb_flux, flux, s_tmp, r_tmp, h1, h2, s_tmpU, s_tmpD

      real (kind=RKIND), dimension(:), pointer :: areaCell, dvEdge, dcEdge
      real (kind=RKIND), dimension(:), pointer :: meshScalingDel2

      real (kind=RKIND), dimension(:,:), pointer :: gradTracerEdge, gradTracerTopOfEdge, gradHTracerSlopedTopOfCell, &
         dTracerdZTopOfCell, dTracerdZTopOfEdge, areaCellSum

      type (field2DReal), pointer :: gradTracerEdgeField, gradTracerTopOfEdgeField, gradHTracerSlopedTopOfCellField, dTracerdZTopOfCellField, dTracerdZTopOfEdgeField, &
         areaCellSumField

      err = 0

      if (.not.del2On) return

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      num_tracers = size(tracers, dim=1)

      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'meshScalingDel2', meshScalingDel2)

      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)

      !
      ! compute a boundary mask to enforce insulating boundary conditions in the horizontal
      !
      do iCell = 1, nCells
        invAreaCell = 1.0 / areaCell(iCell)
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

     !
     ! COMPUTE the extra terms arising due to mismatch between the constant coordinate surfaces and the
     ! isopycnal surfaces.
     !
     ! mrp note: Redi diffusion should be put in a separate subroutine
      if (config_use_standardGM) then

         call mpas_pool_get_field(scratchPool, 'gradTracerEdge', gradTracerEdgeField)
         call mpas_pool_get_field(scratchPool, 'gradTracerTopOfEdge', gradTracerTopOfEdgeField)
         call mpas_pool_get_field(scratchPool, 'gradHTracerSlopedTopOfCell', gradHTracerSlopedTopOfCellField)
         call mpas_pool_get_field(scratchPool, 'dTracerdZTopOfCell', dTracerdZTopOfCellField)
         call mpas_pool_get_field(scratchPool, 'dTracerdZTopOfEdge', dTracerdZTopOfEdgeField)
         call mpas_pool_get_field(scratchPool, 'areaCellSum', areaCellSumField)

         call mpas_allocate_scratch_field(gradTracerEdgeField, .true.)
         call mpas_allocate_scratch_field(gradTracerTopOfEdgeField, .true.)
         call mpas_allocate_scratch_field(gradHTracerSlopedTopOfCellField, .true.)
         call mpas_allocate_scratch_field(dTracerdZTopOfCellField, .true.)
         call mpas_allocate_scratch_field(dTracerdZTopOfEdgeField, .true.)
         call mpas_allocate_scratch_field(areaCellSumField, .True.)

         gradTracerEdge => gradTracerEdgeField % array
         gradTracerTopOfEdge => gradTracerTopOfEdgeField % array
         gradHTracerSlopedTopOfCell => gradHTracerSlopedTopOfCellField % array
         dTracerdZTopOfCell => dTracerdZTopOfCellField % array
         dTracerdZTopOfEdge => dTracerdZTopOfEdgeField % array
         areaCellSum => areaCellSumField % array

         gradTracerEdge = 0.0
         gradTracerTopOfEdge = 0.0
         gradHTracerSlopedTopOfCell = 0.0
         dTracerdZTopOfCell = 0.0
         dTracerdZTopOfEdge = 0.0

         ! this is the "standard" del2 term, but forced to use config_redi_kappa
         if(.not.config_disable_redi_horizontal_term1) then
         do iCell = 1, nCells
            invAreaCell = 1.0 / areaCell(iCell)
            do i = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i, iCell)
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)

               r_tmp = config_redi_kappa * dvEdge(iEdge) / dcEdge(iEdge)

               do k = 1, maxLevelEdgeTop(iEdge)

                  ! this is the tapering of config_redi_kappa where abs(slope) > config_max_relative_slope
                  s_tmp = relativeSlopeTapering(k,iEdge)

                  do iTracer = 1, num_tracers
                     ! \kappa_2 \nabla \phi on edge
                     tracer_turb_flux = tracers(iTracer, k, cell2) - tracers(iTracer, k, cell1)

                     ! div(h \kappa_2 \nabla \phi) at cell center
                     flux = layerThicknessEdge(k, iEdge) * tracer_turb_flux * r_tmp * s_tmp

                     tend(iTracer, k, iCell) = tend(iTracer, k, iCell) - edgeSignOnCell(i, iCell) * flux * invAreaCell
                  end do
               end do

            end do
         end do
         endif

         ! Compute vertical derivative of tracers at cell center and top of layer
         do iTracer = 1, num_tracers

            do iCell = 1, nCells
               do k = 2, maxLevelCell(iCell)
                  dTracerdZTopOfCell(k,iCell) = (tracers(iTracer,k-1,iCell) - tracers(iTracer,k,iCell)) / (zMid(k-1,iCell) - zMid(k,iCell))
               end do

               ! Approximation of dTracerdZTopOfCell on the top and bottom interfaces through the idea of having
               ! ghost cells above the top and below the bottom layers of the same depths and tracer density.
               ! Essentially, this enforces the boundary condition (d tracer)/dz = 0 at the top and bottom.
               dTracerdZTopOfCell(1,iCell) = 0.0
               dTracerdZTopOfCell(maxLevelCell(iCell)+1,iCell) = 0.0
            end do

            ! Compute tracer gradient (gradTracerEdge) along the constant coordinate surface.
            ! The computed variables lives at edge and mid-layer depth
            do iEdge = 1, nEdges
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)

               do k=1,maxLevelEdgeTop(iEdge)
                  gradTracerEdge(k,iEdge) = (tracers(iTracer,k,cell2) - tracers(iTracer,k,cell1)) / dcEdge(iEdge)
               end do
            end do

            ! Interpolate dTracerdZTopOfCell to edge and top of layer
            do iEdge = 1, nEdges
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               do k = 1, maxLevelEdgeTop(iEdge)
                  dTracerdZTopOfEdge(k,iEdge) = 0.5 * (dTracerdZTopOfCell(k,cell1) + dTracerdZTopOfCell(k,cell2)) 
               end do
               dTracerdZTopOfEdge(maxLevelEdgeTop(iEdge)+1,iEdge) = 0.0
            end do

            ! Interpolate gradTracerEdge to edge and top of layer
            do iEdge = 1, nEdges
               do k = 2, maxLevelEdgeTop(iEdge)
                  h1 = layerThicknessEdge(k-1,iEdge)
                  h2 = layerThicknessEdge(k,iEdge)

                  ! Using second-order interpolation below
                  gradTracerTopOfEdge(k,iEdge) = (h2 * gradTracerEdge(k-1,iEdge) + h1 * gradTracerEdge(k,iEdge)) / (h1 + h2)
               end do

               ! Approximation of values on the top and bottom interfaces through the idea of having ghost cells above
               ! the top and below the bottom layers of the same depths and tracer concentration.
               gradTracerTopOfEdge(1,iEdge) = gradTracerEdge(1,iEdge)
               gradTracerTopOfEdge(maxLevelEdgeTop(iEdge)+1,iEdge) = gradTracerEdge(max(maxLevelEdgeTop(iEdge),1),iEdge)
            end do

            ! Compute \nabla\cdot(relativeSlope d\phi/dz)
            if(.not.config_disable_redi_horizontal_term2) then
            do iEdge = 1, nEdges
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               invAreaCell1 = 1./areaCell(cell1)
               invAreaCell2 = 1./areaCell(cell2)

               do k = 1, maxLevelEdgeTop(iEdge)
                  s_tmpU = relativeSlopeTapering(k  , iEdge) * relativeSlopeTopOfEdge(k,iEdge)*dTracerdZTopOfEdge(k,iEdge)
                  s_tmpD = relativeSlopeTapering(k+1, iEdge) * relativeSlopeTopOfEdge(k+1,iEdge)*dTracerdZTopOfEdge(k+1,iEdge)
                  flux = 0.5*dvEdge(iEdge)*(s_tmpU + s_tmpD)
                  flux = flux * layerThicknessEdge(k, iEdge)
                  tend(iTracer,k,cell1) = tend(iTracer,k,cell1) + config_Redi_kappa * flux * invAreaCell1
                  tend(iTracer,k,cell2) = tend(iTracer,k,cell2) - config_Redi_kappa * flux * invAreaCell2
               end do

            end do
            endif

            ! Compute dz * d(relativeSlope\cdot\nabla\phi)/dz  (so the dz cancel out)
            gradHTracerSlopedTopOfCell = 0.0

            ! Compute relativeSlope\cdot\nabla\phi (variable gradHTracerSlopedTopOfCell) at non-boundary edges
            areaCellSum = 1.0e-34
            do iEdge = 1, nEdges
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               ! contribution of cell area from this edge:
               areaEdge = 0.25 * dcEdge(iEdge) * dvEdge(iEdge)

               do k = 1, maxLevelEdgeTop(iEdge) 
                  r_tmp = 2.0 * areaEdge * relativeSlopeTopOfEdge(k,iEdge) * gradTracerTopOfEdge(k,iEdge)
                  gradHTracerSlopedTopOfCell(k,cell1) = gradHTracerSlopedTopOfCell(k,cell1) + r_tmp
                  gradHTracerSlopedTopOfCell(k,cell2) = gradHTracerSlopedTopOfCell(k,cell2) + r_tmp

                  areaCellSum(k,cell1) = areaCellSum(k,cell1) + areaEdge
                  areaCellSum(k,cell2) = areaCellSum(k,cell2) + areaEdge

               end do
            end do
            do iCell=1,nCells
               do k = 1, maxLevelCell(iCell)
                  gradHTracerSlopedTopOfCell(k,iCell) = gradHTracerSlopedTopOfCell(k,iCell)/areaCellSum(k,iCell)
               end do
            end do

            if(.not.config_disable_redi_horizontal_term3) then
            do iCell = 1, nCells
               ! impose no-flux boundary conditions at top and bottom of column
               gradHTracerSlopedTopOfCell(1,iCell) = 0.0
               gradHTracerSlopedTopOfCell(maxLevelCell(iCell)+1,iCell) = 0.0
               do k = 1, maxLevelCell(iCell)
                  s_tmp = relativeSlopeTaperingCell(k,iCell)
                  tend(iTracer,k,iCell) = tend(iTracer,k,iCell) + s_tmp * config_Redi_kappa * (gradHTracerSlopedTopOfCell(k,iCell) - gradHTracerSlopedTopOfCell(k+1,iCell))
               end do
            end do
            endif

         end do  ! iTracer

         call mpas_deallocate_scratch_field(gradTracerEdgeField, .true.)
         call mpas_deallocate_scratch_field(gradTracerTopOfEdgeField, .true.)
         call mpas_deallocate_scratch_field(gradHTracerSlopedTopOfCellField, .true.)
         call mpas_deallocate_scratch_field(dTracerdZTopOfCellField, .true.)
         call mpas_deallocate_scratch_field(dTracerdZTopOfEdgeField, .true.)

      end if ! config_use_standardGM

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
      call mpas_pool_get_config(ocnConfigs, 'config_use_standardGM',config_use_standardGM)
      call mpas_pool_get_config(ocnConfigs, 'config_Redi_kappa',config_Redi_kappa)
      call mpas_pool_get_config(ocnConfigs, 'config_disable_redi_horizontal_term1',config_disable_redi_horizontal_term1)
      call mpas_pool_get_config(ocnConfigs, 'config_disable_redi_horizontal_term2',config_disable_redi_horizontal_term2)
      call mpas_pool_get_config(ocnConfigs, 'config_disable_redi_horizontal_term3',config_disable_redi_horizontal_term3)

      del2on = .false.

      if ( config_use_tracer_del2 ) then
      if ( config_tracer_del2 > 0.0 ) then
          del2On = .true.
          eddyDiff2 = config_tracer_del2
      endif
      endif

     if ( config_use_standardGM ) then
     if ( config_Redi_kappa > 0.0 ) then
          del2On = .true.
     endif
     endif


   !--------------------------------------------------------------------

   end subroutine ocn_tracer_hmix_del2_init!}}}

!***********************************************************************

end module ocn_tracer_hmix_del2

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
