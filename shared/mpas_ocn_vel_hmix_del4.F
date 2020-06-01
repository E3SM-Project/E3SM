! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_hmix_del4
!
!> \brief Ocean horizontal mixing - biharmonic parameterization
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011
!> \details
!>  This module contains routines and variables for computing
!>  horizontal mixing tendencies using a biharmonic formulation.
!
!-----------------------------------------------------------------------

module ocn_vel_hmix_del4

   use mpas_timer
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_threading
   use mpas_vector_operations
   use mpas_matrix_operations
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

   public :: ocn_vel_hmix_del4_tend, &
             ocn_vel_hmix_del4_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: hmixDel4On       !< local flag to determine whether del4 chosen

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vel_hmix_del4_tend
!
!> \brief   Computes tendency term for biharmonic horizontal momentum mixing
!> \author  Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date    September 2011
!> \details
!>  This routine computes the horizontal mixing tendency for momentum
!>  based on a biharmonic form for the mixing.  This mixing tendency
!>  takes the form  \f$-\nu_4 \nabla^4 u\f$
!>  but is computed as
!>  \f$\nabla^2 u = \nabla divergence + k \times \nabla relativeVorticity\f$
!>  applied recursively.
!>  This formulation is only valid for constant \f$\nu_4\f$ .
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_hmix_del4_tend(meshPool, scratchPool, divergence, relativeVorticity, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         divergence      !< Input: velocity divergence

      type (mpas_pool_type), intent(in) :: scratchPool !< Input: scratch variables

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         relativeVorticity       !< Input: relative vorticity

      type (mpas_pool_type), intent(in) :: &
         meshPool           !< Input: mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend       !< Input/Output: velocity tendency

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

      integer :: iEdge, cell1, cell2, vertex1, vertex2, k, i
      integer :: iCell, iVertex, nEdges, nCells, nVertices
      integer, pointer :: nVertLevels, vertexDegree
      integer, dimension(:), pointer :: nEdgesArray, nCellsArray, nVerticesArray

      integer, dimension(:), pointer :: maxLevelEdgeTop, maxLevelVertexTop, &
            maxLevelCell, nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge, verticesOnEdge, edgeMask, edgesOnVertex, edgesOnCell, edgeSignOnVertex, &
                                          edgeSignOnCell


      real (kind=RKIND) :: u_diffusion, invAreaCell1, invAreaCell2, invAreaTri1, &
            invAreaTri2, invDcEdge, invDvEdge, r_tmp
      real (kind=RKIND), dimension(:), pointer :: dcEdge, dvEdge, areaTriangle, &
            meshScalingDel4, areaCell

      real (kind=RKIND), dimension(:,:), pointer :: delsq_divergence, delsq_relativeVorticity, delsq_u
      type (field2DReal), pointer :: delsq_uField, delsq_divergenceField, delsq_relativeVorticityField

      real (kind=RKIND), pointer :: config_mom_del4, config_mom_del4_div_factor

      err = 0

      if(.not.hmixDel4On) return

      call mpas_timer_start("vel del4")

      call mpas_pool_get_config(ocnConfigs, 'config_mom_del4', config_mom_del4)
      call mpas_pool_get_config(ocnConfigs, 'config_mom_del4_div_factor', config_mom_del4_div_factor)

      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)
      call mpas_pool_get_dimension(meshPool, 'nVerticesArray', nVerticesArray)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree)

      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'maxLevelVertexTop', maxLevelVertexTop)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'areaTriangle', areaTriangle)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'meshScalingDel4', meshScalingDel4)
      call mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnVertex', edgesOnVertex)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgeSignOnVertex', edgeSignOnVertex)
      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)

      call mpas_pool_get_field(scratchPool, 'delsq_u', delsq_uField)
      call mpas_pool_get_field(scratchPool, 'delsq_divergence', delsq_divergenceField)
      call mpas_pool_get_field(scratchPool, 'delsq_relativeVorticity', delsq_relativeVorticityField)
      call mpas_allocate_scratch_field(delsq_uField, .true., .false.)
      call mpas_allocate_scratch_field(delsq_divergenceField, .true., .false.)
      call mpas_allocate_scratch_field(delsq_relativeVorticityField, .true., .false.)

      delsq_u => delsq_uField % array
      delsq_divergence => delsq_divergenceField % array
      delsq_relativeVorticity => delsq_relativeVorticityField % array

      nEdges = nEdgesArray( 3 )

      !Compute delsq_u
      !$omp do schedule(runtime)
      do iEdge = 1, nEdges
         delsq_u(:, iEdge) = 0.0_RKIND
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         vertex1 = verticesOnEdge(1,iEdge)
         vertex2 = verticesOnEdge(2,iEdge)

         invDcEdge = 1.0_RKIND / dcEdge(iEdge)
         invDvEdge = 1.0_RKIND / max(dvEdge(iEdge), 0.25_RKIND*dcEdge(iEdge))

         do k=1,maxLevelEdgeTop(iEdge)
            ! Compute \nabla^2 u = \nabla divergence + k \times \nabla relativeVorticity
            delsq_u(k, iEdge) = ( divergence(k,cell2)  - divergence(k,cell1) ) * invDcEdge  &
                               -( relativeVorticity(k,vertex2) - relativeVorticity(k,vertex1)) * invDvEdge
         end do
      end do
      !$omp end do

      nVertices = nVerticesArray( 2 )

      ! Compute delsq_relativeVorticity
      !$omp do schedule(runtime)
      do iVertex = 1, nVertices
         delsq_relativeVorticity(:, iVertex) = 0.0_RKIND
         invAreaTri1 = 1.0_RKIND / areaTriangle(iVertex)
         do i = 1, vertexDegree
            iEdge = edgesOnVertex(i, iVertex)
            do k = 1, maxLevelVertexTop(iVertex)
               delsq_relativeVorticity(k, iVertex) = delsq_relativeVorticity(k, iVertex) + edgeSignOnVertex(i, iVertex) &
                                                   * dcEdge(iEdge) * delsq_u(k, iEdge) * invAreaTri1
            end do
         end do
      end do
      !$omp end do

      nCells = nCellsArray( 2 )

      ! Compute delsq_divergence
      !$omp do schedule(runtime)
      do iCell = 1, nCells
         delsq_divergence(:, iCell) = 0.0_RKIND
         invAreaCell1 = 1.0_RKIND / areaCell(iCell)
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            do k = 1, maxLevelCell(iCell)
               delsq_divergence(k, iCell) = delsq_divergence(k, iCell) - edgeSignOnCell(i, iCell) * dvEdge(iEdge) &
                                          * delsq_u(k, iEdge) * invAreaCell1
            end do
         end do
      end do
      !$omp end do

      nEdges = nEdgesArray( 1 )

      ! Compute - \kappa \nabla^4 u
      ! as  \nabla div(\nabla^2 u) + k \times \nabla ( k \cross curl(\nabla^2 u) )
      !$omp do schedule(runtime)
      do iEdge = 1, nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         vertex1 = verticesOnEdge(1,iEdge)
         vertex2 = verticesOnEdge(2,iEdge)

         invDcEdge = 1.0_RKIND / dcEdge(iEdge)
         invDvEdge = 1.0_RKIND / dvEdge(iEdge)
         r_tmp = config_mom_del4 * meshScalingDel4(iEdge)

         do k=1,maxLevelEdgeTop(iEdge)
            u_diffusion = config_mom_del4_div_factor*(delsq_divergence(k,cell2) - delsq_divergence(k,cell1)) * invDcEdge  &
                        - (delsq_relativeVorticity(k,vertex2) - delsq_relativeVorticity(k,vertex1) ) * invDvEdge

            tend(k,iEdge) = tend(k,iEdge) - edgeMask(k, iEdge) * u_diffusion * r_tmp
         end do
      end do
      !$omp end do

      call mpas_deallocate_scratch_field(delsq_uField, .true.)
      call mpas_deallocate_scratch_field(delsq_divergenceField, .true.)
      call mpas_deallocate_scratch_field(delsq_relativeVorticityField, .true.)

      call mpas_timer_stop("vel del4")

   !--------------------------------------------------------------------

   end subroutine ocn_vel_hmix_del4_tend!}}}

!***********************************************************************
!
!  routine ocn_vel_hmix_del4_init
!
!> \brief   Initializes ocean momentum biharmonic horizontal mixing
!> \author  Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date    September 2011
!> \details
!>  This routine initializes a variety of quantities related to
!>  biharmonic horizontal tracer mixing in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_hmix_del4_init(err)!{{{

   integer, intent(out) :: err !< Output: error flag

   !--------------------------------------------------------------------
   !
   ! set some local module variables based on input config choices
   !
   !--------------------------------------------------------------------

   real (kind=RKIND), pointer :: config_mom_del4
   logical, pointer :: config_use_mom_del4

   err = 0

   call mpas_pool_get_config(ocnConfigs, 'config_mom_del4', config_mom_del4)
   call mpas_pool_get_config(ocnConfigs, 'config_use_mom_del4', config_use_mom_del4)

   hmixDel4On = .false.

   if ( config_mom_del4 > 0.0_RKIND ) then
      hmixDel4On = .true.
   endif

   if(.not.config_use_mom_del4) hmixDel4On = .false.

   !--------------------------------------------------------------------

   end subroutine ocn_vel_hmix_del4_init!}}}

!***********************************************************************

end module ocn_vel_hmix_del4

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
