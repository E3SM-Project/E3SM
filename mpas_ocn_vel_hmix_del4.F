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
   use mpas_configure
   use mpas_vector_operations
   use mpas_matrix_operations
   use mpas_tensor_operations

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

   subroutine ocn_vel_hmix_del4_tend(mesh, divergence, relativeVorticity, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         divergence      !< Input: velocity divergence

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         relativeVorticity       !< Input: relative vorticity

      type (mesh_type), intent(in) :: &
         mesh           !< Input: mesh information

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
      integer :: nVertices, nVertLevels, nCells, nEdges, nEdgesSolve, vertexDegree

      integer, dimension(:), pointer :: maxLevelEdgeTop, maxLevelVertexTop, &
            maxLevelCell, nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge, verticesOnEdge, edgeMask, edgesOnVertex, edgesOnCell, edgeSignOnVertex, edgeSignOnCell


      real (kind=RKIND) :: u_diffusion, invAreaCell1, invAreaCell2, invAreaTri1, &
            invAreaTri2, invDcEdge, invDvEdge, r_tmp
      real (kind=RKIND), dimension(:), pointer :: dcEdge, dvEdge, areaTriangle, &
            meshScalingDel4, areaCell

      real (kind=RKIND), dimension(:,:), allocatable :: delsq_divergence, &
            delsq_circulation, delsq_relativeVorticity, delsq_u

      err = 0

      if(.not.hmixDel4On) return

      nCells = mesh % nCells
      nEdges = mesh % nEdges
      nEdgesSolve = mesh % nEdgessolve
      nVertices = mesh % nVertices
      nVertLevels = mesh % nVertLevels
      vertexDegree = mesh % vertexDegree

      maxLevelEdgeTop => mesh % maxLevelEdgeTop % array
      maxLevelVertexTop => mesh % maxLevelVertexTop % array
      maxLevelCell => mesh % maxLevelCell % array
      cellsOnEdge => mesh % cellsOnEdge % array
      verticesOnEdge => mesh % verticesOnEdge % array
      dcEdge => mesh % dcEdge % array
      dvEdge => mesh % dvEdge % array
      areaTriangle => mesh % areaTriangle % array
      areaCell => mesh % areaCell % array
      meshScalingDel4 => mesh % meshScalingDel4 % array
      edgeMask => mesh % edgeMask % array
      nEdgesOnCell => mesh % nEdgesOnCell % array
      edgesOnVertex => mesh % edgesOnVertex % array
      edgesOnCell => mesh % edgesOnCell % array
      edgeSignOnVertex => mesh % edgeSignOnVertex % array
      edgeSignOnCell => mesh % edgeSignOnCell % array

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

   subroutine ocn_vel_hmix_del4_tensor_tend(mesh, normalVelocity, tangentialVelocity, viscosity, scratch, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity     !< Input: velocity normal to an edge

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         tangentialVelocity     !< Input: velocity, tangent to an edge

      type (mesh_type), intent(in) :: &
         mesh            !< Input: mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         viscosity       !< Input/Output: viscosity

      type (scratch_type), intent(inout) :: &
         scratch !< Input/Output: Scratch structure

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

      integer :: iEdge, nEdgesSolve, nEdges, k, nVertLevels
      integer, dimension(:), pointer :: maxLevelEdgeTop
      integer, dimension(:,:), pointer :: edgeMask, edgeSignOnCell

      real (kind=RKIND) :: visc4_sqrt
      real (kind=RKIND), dimension(:), pointer :: meshScalingDel4
      real (kind=RKIND), dimension(:,:), pointer :: normalVectorEdge, tangentialVectorEdge, edgeTangentVectors
      real (kind=RKIND), dimension(:,:,:), pointer :: &
         strainRateR3Cell, strainRateR3Edge, divTensorR3Cell, outerProductEdge

      !-----------------------------------------------------------------
      !
      ! exit if this mixing is not selected
      !
      !-----------------------------------------------------------------

      err = 0

      if(.not.config_use_mom_del4_tensor) return

      nEdges = mesh % nEdges
      nVertLevels = mesh % nVertLevels
      nEdgesSolve = mesh % nEdgesSolve
      maxLevelEdgeTop => mesh % maxLevelEdgeTop % array
      meshScalingDel4 => mesh % meshScalingDel4 % array
      edgeMask => mesh % edgeMask % array
      edgeSignOnCell => mesh % edgeSignOnCell % array
      edgeTangentVectors => mesh % edgeTangentVectors % array

      call mpas_allocate_scratch_field(scratch % strainRateR3Cell, .true.)
      call mpas_allocate_scratch_field(scratch % strainRateR3Edge, .true.)
      call mpas_allocate_scratch_field(scratch % divTensorR3Cell, .true.)
      call mpas_allocate_scratch_field(scratch % outerProductEdge, .true.)
      call mpas_allocate_scratch_field(scratch % normalVectorEdge, .true.)
      call mpas_allocate_scratch_field(scratch % tangentialVectorEdge, .true.)

      strainRateR3Cell => scratch % strainRateR3Cell % array
      strainRateR3Edge => scratch % strainRateR3Edge % array
      divTensorR3Cell  => scratch % divTensorR3Cell % array
      outerProductEdge => scratch % outerProductEdge % array
      normalVectorEdge => scratch % normalVectorEdge % array
      tangentialVectorEdge => scratch % tangentialVectorEdge % array

      !!!!!!! first div(grad())

      call mpas_strain_rate_R3Cell(normalVelocity, tangentialVelocity, &
         mesh, edgeSignOnCell, edgeTangentVectors, .true., &
         outerProductEdge, strainRateR3Cell)

      call mpas_matrix_cell_to_edge(strainRateR3Cell, mesh, .true., strainRateR3Edge)

      ! The following loop could possibly be reduced to nEdgesSolve
      do iEdge=1,nEdges 
         visc4_sqrt = sqrt(config_mom_del4_tensor * meshScalingDel4(iEdge))
         do k=1,maxLevelEdgeTop(iEdge)
            strainRateR3Edge(:,k,iEdge) = visc4_sqrt * strainRateR3Edge(:,k,iEdge) 
         end do
         ! Impose zero strain rate at land boundaries
         do k=maxLevelEdgeTop(iEdge)+1,nVertLevels
            strainRateR3Edge(:,k,iEdge) = 0.0
         end do
      end do

      ! may change boundaries to false later
      call mpas_divergence_of_tensor_R3Cell(strainRateR3Edge, mesh, edgeSignOnCell, .true., divTensorR3Cell)

      call mpas_vector_R3Cell_to_2DEdge(divTensorR3Cell, mesh, edgeTangentVectors, .true., normalVectorEdge, tangentialVectorEdge)

      !!!!!!! second div(grad())

      call mpas_strain_rate_R3Cell(normalVectorEdge, tangentialVectorEdge, &
         mesh, edgeSignOnCell, edgeTangentVectors, .true., &
         outerProductEdge, strainRateR3Cell)

      call mpas_matrix_cell_to_edge(strainRateR3Cell, mesh, .true., strainRateR3Edge)

      ! The following loop could possibly be reduced to nEdgesSolve
      do iEdge=1,nEdges  
         visc4_sqrt = sqrt(config_mom_del4_tensor * meshScalingDel4(iEdge))
         viscosity(:,iEdge) = viscosity(:,iEdge) + config_mom_del4_tensor * meshScalingDel4(iEdge)
         do k=1,maxLevelEdgeTop(iEdge)
            strainRateR3Edge(:,k,iEdge) = visc4_sqrt * strainRateR3Edge(:,k,iEdge) 
         end do
         ! Impose zero strain rate at land boundaries
         do k=maxLevelEdgeTop(iEdge)+1,nVertLevels
            strainRateR3Edge(:,k,iEdge) = 0.0
         end do
      end do

      ! may change boundaries to false later
      call mpas_divergence_of_tensor_R3Cell(strainRateR3Edge, mesh, edgeSignOnCell, .true., divTensorR3Cell)

      call mpas_vector_R3Cell_to_normalVectorEdge(divTensorR3Cell, mesh, .true., normalVectorEdge)

      ! The following loop could possibly be reduced to nEdgesSolve
      do iEdge=1,nEdges
         do k=1,maxLevelEdgeTop(iEdge)
            tend(k,iEdge) = tend(k,iEdge) - edgeMask(k, iEdge) * normalVectorEdge(k,iEdge)
         end do
      end do

      call mpas_deallocate_scratch_field(scratch % strainRateR3Cell, .true.)
      call mpas_deallocate_scratch_field(scratch % strainRateR3Edge, .true.)
      call mpas_deallocate_scratch_field(scratch % divTensorR3Cell, .true.)
      call mpas_deallocate_scratch_field(scratch % outerProductEdge, .true.)
      call mpas_deallocate_scratch_field(scratch % normalVectorEdge, .true.)
      call mpas_deallocate_scratch_field(scratch % tangentialVectorEdge, .true.)

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

   err = 0

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
