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

   use mpas_grid_types
   use mpas_vector_operations
   use mpas_matrix_operations
   use mpas_tensor_operations
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
             ocn_vel_hmix_del4_tensor_tend, &
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

   subroutine ocn_vel_hmix_del4_tend(meshPool, divergence, relativeVorticity, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         divergence      !< Input: velocity divergence

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
      integer :: iCell, iVertex
      integer, pointer :: nVertices, nVertLevels, nCells, nEdges, nEdgesSolve, vertexDegree

      integer, dimension(:), pointer :: maxLevelEdgeTop, maxLevelVertexTop, &
            maxLevelCell, nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge, verticesOnEdge, edgeMask, edgesOnVertex, edgesOnCell, edgeSignOnVertex, edgeSignOnCell


      real (kind=RKIND) :: u_diffusion, invAreaCell1, invAreaCell2, invAreaTri1, &
            invAreaTri2, invDcEdge, invDvEdge, r_tmp
      real (kind=RKIND), dimension(:), pointer :: dcEdge, dvEdge, areaTriangle, &
            meshScalingDel4, areaCell

      real (kind=RKIND), dimension(:,:), allocatable :: delsq_divergence, &
            delsq_circulation, delsq_relativeVorticity, delsq_u

      real (kind=RKIND), pointer :: config_mom_del4

      err = 0

      if(.not.hmixDel4On) return

      call mpas_pool_get_config(ocnConfigs, 'config_mom_del4', config_mom_del4)

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
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

      allocate(delsq_u(nVertLEvels, nEdges+1))
      allocate(delsq_divergence(nVertLevels, nCells+1))
      allocate(delsq_relativeVorticity(nVertLevels, nVertices+1))

      delsq_u(:,:) = 0.0
      delsq_relativeVorticity(:,:) = 0.0
      delsq_divergence(:,:) = 0.0

      !Compute delsq_u
      do iEdge = 1, nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         vertex1 = verticesOnEdge(1,iEdge)
         vertex2 = verticesOnEdge(2,iEdge)

         invDcEdge = 1.0 / dcEdge(iEdge)
         invDvEdge = 1.0 / dvEdge(iEdge)

         do k=1,maxLevelEdgeTop(iEdge)
            ! Compute \nabla^2 u = \nabla divergence + k \times \nabla relativeVorticity
            delsq_u(k, iEdge) = ( divergence(k,cell2)  - divergence(k,cell1) ) * invDcEdge  &
                               -( relativeVorticity(k,vertex2) - relativeVorticity(k,vertex1)) * invDcEdge * sqrt(3.0)   
         end do
      end do

      ! Compute delsq_relativeVorticity
      do iVertex = 1, nVertices
         invAreaTri1 = 1.0 / areaTriangle(iVertex)
         do i = 1, vertexDegree
            iEdge = edgesOnVertex(i, iVertex)
            do k = 1, maxLevelVertexTop(iVertex)
               delsq_relativeVorticity(k, iVertex) = delsq_relativeVorticity(k, iVertex) + edgeSignOnVertex(i, iVertex) * dcEdge(iEdge) * delsq_u(k, iEdge) * invAreaTri1
            end do
         end do
      end do

      ! Compute delsq_divergence
      do iCell = 1, nCells
         invAreaCell1 = 1.0 / areaCell(iCell)
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            do k = 1, maxLevelCell(iCell)
               delsq_divergence(k, iCell) = delsq_divergence(k, iCell) - edgeSignOnCell(i, iCell) * dvEdge(iEdge) * delsq_u(k, iEdge) * invAreaCell1
            end do
         end do
      end do

      ! Compute - \kappa \nabla^4 u 
      ! as  \nabla div(\nabla^2 u) + k \times \nabla ( k \cross curl(\nabla^2 u) )
      do iEdge=1,nEdgesSolve
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         vertex1 = verticesOnEdge(1,iEdge)
         vertex2 = verticesOnEdge(2,iEdge)

         invDcEdge = 1.0 / dcEdge(iEdge)
         invDvEdge = 1.0 / dvEdge(iEdge)
         r_tmp = config_mom_del4 * meshScalingDel4(iEdge)

         do k=1,maxLevelEdgeTop(iEdge)
            u_diffusion = (delsq_divergence(k,cell2) - delsq_divergence(k,cell1)) * invDcEdge  &
                        - (delsq_relativeVorticity(k,vertex2) - delsq_relativeVorticity(k,vertex1) ) * invDcEdge * sqrt(3.0) 

            tend(k,iEdge) = tend(k,iEdge) - edgeMask(k, iEdge) * u_diffusion * r_tmp
         end do
      end do

      deallocate(delsq_u)
      deallocate(delsq_divergence)
      deallocate(delsq_relativeVorticity)

   !--------------------------------------------------------------------

   end subroutine ocn_vel_hmix_del4_tend!}}}

!***********************************************************************
!
!  routine ocn_vel_hmix_del4_tensor_tend
!
!> \brief   Computes tendency term for Laplacian horizontal momentum mixing
!> \author  Mark Petersen
!> \date    July 2013
!> \details 
!>  This routine computes the horizontal mixing tendency for momentum
!>  using tensor operations, 
!>  based on a Laplacian form for the mixing, 
!>  \f$-\nabla\cdot( \sqrt{\nu_4} \nabla(\nabla\cdot( \sqrt{\nu_4} \nabla(u))))\f$
!>  where \f$\nu_4\f$ is the del4 viscosity.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_hmix_del4_tensor_tend(meshPool, normalVelocity, tangentialVelocity, viscosity, scratchPool, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity     !< Input: velocity normal to an edge

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         tangentialVelocity     !< Input: velocity, tangent to an edge

      type (mpas_pool_type), intent(in) :: &
         meshPool            !< Input: mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         viscosity       !< Input/Output: viscosity

      type (mpas_pool_type), intent(inout) :: &
         scratchPool !< Input/Output: Scratch structure

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend             !< Input/Output: velocity tendency

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

      integer :: iEdge, k
      integer, pointer :: nEdgesSolve, nEdges, nVertLevels
      integer, dimension(:), pointer :: maxLevelEdgeTop
      integer, dimension(:,:), pointer :: edgeMask, edgeSignOnCell

      real (kind=RKIND) :: visc4_sqrt
      real (kind=RKIND), dimension(:), pointer :: meshScalingDel4
      real (kind=RKIND), dimension(:,:), pointer :: normalVectorEdge, tangentialVectorEdge, edgeTangentVectors
      real (kind=RKIND), dimension(:,:,:), pointer :: &
         strainRateR3Cell, strainRateR3Edge, divTensorR3Cell, outerProductEdge

      type (field2DReal), pointer :: normalVectorEdgeField, tangentialVectorEdgeField
      type (field3DReal), pointer :: strainRateR3CellField, strainRateR3EdgeField, divTensorR3CellField, outerProductEdgeField

      logical, pointer :: config_use_mom_del4_tensor
      real (kind=RKIND), pointer :: config_mom_del4_tensor


      !-----------------------------------------------------------------
      !
      ! exit if this mixing is not selected
      !
      !-----------------------------------------------------------------

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_use_mom_del4_tensor', config_use_mom_del4_tensor)

      if(.not.config_use_mom_del4_tensor) return

      call mpas_pool_get_config(ocnConfigs, 'config_mom_del4_tensor', config_mom_del4_tensor)

      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'meshScalingDel4', meshScalingDel4)
      call mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)
      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
      call mpas_pool_get_array(meshPool, 'edgeTangentVectors', edgeTangentVectors)

      call mpas_pool_get_field(scratchPool, 'strainRateR3Cell', strainRateR3CellField)
      call mpas_pool_get_field(scratchPool, 'strainRateR3Edge', strainRateR3EdgeField)
      call mpas_pool_get_field(scratchPool, 'divTensorR3Cell', divTensorR3CellField)
      call mpas_pool_get_field(scratchPool, 'outerProductEdge', outerProductEdgeField)
      call mpas_pool_get_field(scratchPool, 'normalVectorEdge', normalVectorEdgeField)
      call mpas_pool_get_field(scratchPool, 'tangentialVectorEdge', tangentialVectorEdgeField)

      call mpas_allocate_scratch_field(strainRateR3CellField, .true.)
      call mpas_allocate_scratch_field(strainRateR3EdgeField, .true.)
      call mpas_allocate_scratch_field(divTensorR3CellField, .true.)
      call mpas_allocate_scratch_field(outerProductEdgeField, .true.)
      call mpas_allocate_scratch_field(normalVectorEdgeField, .true.)
      call mpas_allocate_scratch_field(tangentialVectorEdgeField, .true.)

      strainRateR3Cell => strainRateR3CellField % array
      strainRateR3Edge => strainRateR3EdgeField % array
      divTensorR3Cell  => divTensorR3CellField % array
      outerProductEdge => outerProductEdgeField % array
      normalVectorEdge => normalVectorEdgeField % array
      tangentialVectorEdge => tangentialVectorEdgeField % array

      !!!!!!! first div(grad())

      call mpas_strain_rate_R3Cell(normalVelocity, tangentialVelocity, &
         meshPool, edgeSignOnCell, edgeTangentVectors, .true., &
         outerProductEdge, strainRateR3Cell)

      call mpas_matrix_cell_to_edge(strainRateR3Cell, meshPool, .true., strainRateR3Edge)

      ! The following loop could possibly be reduced to nEdgesSolve
      do iEdge = 1, nEdges 
         visc4_sqrt = sqrt(config_mom_del4_tensor * meshScalingDel4(iEdge))
         do k = 1, maxLevelEdgeTop(iEdge)
            strainRateR3Edge(:,k,iEdge) = visc4_sqrt * strainRateR3Edge(:,k,iEdge) 
         end do
         ! Impose zero strain rate at land boundaries
         do k = maxLevelEdgeTop(iEdge)+1, nVertLevels
            strainRateR3Edge(:,k,iEdge) = 0.0
         end do
      end do

      ! may change boundaries to false later
      call mpas_divergence_of_tensor_R3Cell(strainRateR3Edge, meshPool, edgeSignOnCell, .true., divTensorR3Cell)

      call mpas_vector_R3Cell_to_2DEdge(divTensorR3Cell, meshPool, edgeTangentVectors, .true., normalVectorEdge, tangentialVectorEdge)

      !!!!!!! second div(grad())

      call mpas_strain_rate_R3Cell(normalVectorEdge, tangentialVectorEdge, &
         meshPool, edgeSignOnCell, edgeTangentVectors, .true., &
         outerProductEdge, strainRateR3Cell)

      call mpas_matrix_cell_to_edge(strainRateR3Cell, meshPool, .true., strainRateR3Edge)

      ! The following loop could possibly be reduced to nEdgesSolve
      do iEdge = 1, nEdges  
         visc4_sqrt = sqrt(config_mom_del4_tensor * meshScalingDel4(iEdge))
         viscosity(:,iEdge) = viscosity(:,iEdge) + config_mom_del4_tensor * meshScalingDel4(iEdge)
         do k = 1, maxLevelEdgeTop(iEdge)
            strainRateR3Edge(:,k,iEdge) = visc4_sqrt * strainRateR3Edge(:,k,iEdge) 
         end do
         ! Impose zero strain rate at land boundaries
         do k = maxLevelEdgeTop(iEdge)+1, nVertLevels
            strainRateR3Edge(:,k,iEdge) = 0.0
         end do
      end do

      ! may change boundaries to false later
      call mpas_divergence_of_tensor_R3Cell(strainRateR3Edge, meshPool, edgeSignOnCell, .true., divTensorR3Cell)

      call mpas_vector_R3Cell_to_normalVectorEdge(divTensorR3Cell, meshPool, .true., normalVectorEdge)

      ! The following loop could possibly be reduced to nEdgesSolve
      do iEdge = 1,nEdges
         do k = 1,maxLevelEdgeTop(iEdge)
            tend(k,iEdge) = tend(k,iEdge) - edgeMask(k, iEdge) * normalVectorEdge(k,iEdge)
         end do
      end do

      call mpas_deallocate_scratch_field(strainRateR3CellField, .true.)
      call mpas_deallocate_scratch_field(strainRateR3EdgeField, .true.)
      call mpas_deallocate_scratch_field(divTensorR3CellField, .true.)
      call mpas_deallocate_scratch_field(outerProductEdgeField, .true.)
      call mpas_deallocate_scratch_field(normalVectorEdgeField, .true.)
      call mpas_deallocate_scratch_field(tangentialVectorEdgeField, .true.)

   !--------------------------------------------------------------------

   end subroutine ocn_vel_hmix_del4_tensor_tend!}}}

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

   if ( config_mom_del4 > 0.0 ) then
      hmixDel4On = .true.
   endif

   if(.not.config_use_mom_del4) hmixDel4On = .false.

   !--------------------------------------------------------------------

   end subroutine ocn_vel_hmix_del4_init!}}}

!***********************************************************************

end module ocn_vel_hmix_del4

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
