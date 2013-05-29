! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
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

   subroutine ocn_tracer_hmix_del4_tend(grid, layerThicknessEdge, tracers, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThicknessEdge    !< Input: thickness at edge

      type (mesh_type), intent(in) :: &
         grid          !< Input: grid information

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

      integer :: iEdge, nEdges, num_tracers, nVertLevels, nCells
      integer :: iTracer, k, iCell, cell1, cell2, i

      integer, dimension(:), pointer :: maxLevelEdgeTop, maxLevelCell, nEdgesOnCell
      integer, dimension(:,:), pointer :: edgeMask, cellsOnEdge, edgesOnCell, edgeSignOnCell

      real (kind=RKIND) :: invAreaCell1, invAreaCell2, tracer_turb_flux, flux, invdcEdge, r_tmp1, r_tmp2

      real (kind=RKIND), dimension(:,:,:), allocatable :: delsq_tracer

      real (kind=RKIND), dimension(:), pointer :: dcEdge, dvEdge, areaCell, meshScalingDel4


      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing tendencies
      ! note that the user can choose multiple options and the 
      !   tendencies will be added together
      !
      !-----------------------------------------------------------------

      err = 0

      if (.not.del4On) return

      nEdges = grid % nEdges
      nCells = grid % nCells
      num_tracers = size(tracers, dim=1)
      nVertLevels = grid % nVertLevels

      maxLevelEdgeTop => grid % maxLevelEdgeTop % array
      maxLevelCell => grid % maxLevelCell % array
      cellsOnEdge => grid % cellsOnEdge % array

      dcEdge => grid % dcEdge % array
      dvEdge => grid % dvEdge % array
      areaCell => grid % areaCell % array
      meshScalingDel4 => grid % meshScalingDel4 % array

      edgeMask => grid % edgeMask % array

      nEdgesOnCell => grid % nEdgesOnCell % array
      edgesOnCell => grid % edgesOnCell % array
      edgeSignOnCell => grid % edgeSignOnCell % array

      allocate(delsq_tracer(num_tracers,nVertLevels, nCells+1))

      delsq_tracer(:,:,:) = 0.0

      ! first del2: div(h \nabla \phi) at cell center
      do iCell = 1, nCells
        invAreaCell1 = 1.0 / areaCell(iCell)
        do i = 1, nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i, iCell)
          invdcEdge = dvEdge(iEdge) / dcEdge(iEdge)
          cell1 = cellsOnEdge(1,iEdge)
          cell2 = cellsOnEdge(2,iEdge)

          do k = 1, maxLevelEdgeTop(iEdge)
            do iTracer = 1, num_tracers * edgeMask(k, iEdge)

              r_tmp1 = invdcEdge * layerThicknessEdge(k, iEdge) * tracers(iTracer, k, cell1)
              r_tmp2 = invdcEdge * layerThicknessEdge(k, iEdge) * tracers(iTracer, k, cell2)

              delsq_tracer(iTracer, k, iCell) = delsq_tracer(iTracer, k, iCell) - edgeSignOnCell(i, iCell) * (r_tmp2 - r_tmp1) * invAreaCell1
            end do
          end do
        end do
      end do

      ! second del2: div(h \nabla [delsq_tracer]) at cell center
      do iCell = 1, nCells
        invAreaCell1 = 1.0 / areaCell(iCell)
        do i = 1, nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i, iCell)
          cell1 = cellsOnEdge(1, iEdge)
          cell2 = cellsOnedge(2, iEdge)

          invdcEdge = meshScalingDel4(iEdge) * dvEdge(iEdge) * eddyDiff4 / dcEdge(iEdge)

          do k = 1, maxLevelEdgeTop(iEdge)
            do iTracer = 1, num_tracers * edgeMask(k, iEdge)
              tracer_turb_flux = (delsq_tracer(iTracer, k, cell2) - delsq_tracer(iTracer, k, cell1))
                
              flux = tracer_turb_flux * invdcEdge

              tend(iTracer, k, iCell) = tend(iTracer, k, iCell) + edgeSignOnCell(i, iCell) * flux * invAreaCell1
            end do
          end do
        end do
      end do

      deallocate(delsq_tracer)
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

      if ( config_tracer_del4 > 0.0 ) then
          del4On = .true.
          eddyDiff4 = config_tracer_del4
      endif

      if(.not.config_use_tracer_del4) del4on = .false.

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_hmix_del4_init!}}}

!***********************************************************************

end module ocn_tracer_hmix_del4

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
