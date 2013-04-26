! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_tensor_operations
!
!> \brief MPAS tensor operations
!> \author Mark Petersen
!> \date    17 April 2013
!> \details
!>  This module contains the routines for computing
!>  the strain rate tensor, the divergence of a tensor,
!>  and a testing routine to verify these work properly.
!
!-----------------------------------------------------------------------

module mpas_tensor_operations

   use mpas_grid_types
   use mpas_constants
   use mpas_vector_operations
   use mpas_matrix_operations

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

   public :: mpas_strain_rate, &
             mpas_divergence_of_tensor, &
             mpas_tensor_edge_3D_to_2D, &
             mpas_tensor_edge_2D_to_3D, &
             mpas_test_tensor

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine mpas_strain_rate
!
!> \brief   Computes 2D strain rate at edges
!> \author  Mark Petersen
!> \date    17 April 2013
!> \details 
!>  This routine computes the strain rate at cell centers using the weak 
!>  derivative.  Output is a 3D strain rate tensor in 6-index format.
!
!-----------------------------------------------------------------------

   subroutine mpas_strain_rate(grid, normalVelocity, tangentialVelocity, &
      strainRate3DCell, outerProductEdge)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity, &!< Input: Horizontal velocity normal to edge
         tangentialVelocity   !< Input: Horizontal velocity tangent to edge

      type (mesh_type), intent(in) :: &
         grid          !< Input: grid information

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(out) :: &
         strainRate3DCell   !< Output: strain rate tensor at cell center, 3D, in symmetric 6-index form

      real (kind=RKIND), dimension(:,:,:), intent(out) :: outerProductEdge


      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iEdge, nEdges, iCell, nCells, i,j,k

      integer, dimension(:), pointer :: maxLevelEdgeBot, nEdgesOnCell, maxLevelCell
      integer, dimension(:,:), pointer :: edgesOnCell, edgeSignOnCell

      real (kind=RKIND) :: invAreaCell
      real (kind=RKIND), dimension(3,3) :: outerProductEdge3x3
      real (kind=RKIND), dimension(:), pointer :: dvEdge, areaCell, angleEdge
      real (kind=RKIND), dimension(:), pointer :: xCell,yCell,xEdge,yEdge ! temp
      real (kind=RKIND), dimension(:,:), pointer :: edgeNormalVectors, edgeTangentVectors

      nEdges = grid % nEdges
      nCells = grid % nCells

      maxLevelEdgeBot   => grid % maxLevelEdgeBot % array
      maxLevelCell      => grid % maxLevelCell % array
      nEdgesOnCell      => grid % nEdgesOnCell % array
      edgesOnCell       => grid % edgesOnCell % array
      edgeSignOnCell    => grid % edgeSignOnCell % array
      dvEdge            => grid % dvEdge % array
      angleEdge         => grid % angleEdge % array
      areaCell          => grid % areaCell % array
      edgeNormalVectors  => grid % edgeNormalVectors % array
      edgeTangentVectors => grid % edgeTangentVectors % array

! temp:
      xCell  => grid % xCell % array
      yCell  => grid % yCell % array
      xEdge  => grid % xEdge % array
      yEdge  => grid % yEdge % array

if (1==2) then
      print '(a,3f20.10)', 'min max mean edgeNormalVectors1  ',minval(edgeNormalVectors(1,1:nEdges)),maxval(edgeNormalVectors(1,1:nEdges)),sum(edgeNormalVectors(1,1:nEdges))/nEdges/grid%nVertLevels
      print '(a,3f20.10)', 'min max mean edgeNormalVectors2  ',minval(edgeNormalVectors(2,1:nEdges)),maxval(edgeNormalVectors(2,1:nEdges)),sum(edgeNormalVectors(2,1:nEdges))/nEdges/grid%nVertLevels
      print '(a,3f20.10)', 'min max mean edgeNormalVectors3  ',minval(edgeNormalVectors(3,1:nEdges)),maxval(edgeNormalVectors(3,1:nEdges)),sum(edgeNormalVectors(3,1:nEdges))/nEdges/grid%nVertLevels
      print '(a,3f20.10)', 'min max mean edgeTangentVectors1 ',minval(edgeTangentVectors(1,1:nEdges)),maxval(edgeTangentVectors(1,1:nEdges)),sum(edgeTangentVectors(1,1:nEdges))/nEdges/grid%nVertLevels
      print '(a,3f20.10)', 'min max mean edgeTangentVectors2 ',minval(edgeTangentVectors(2,1:nEdges)),maxval(edgeTangentVectors(2,1:nEdges)),sum(edgeTangentVectors(2,1:nEdges))/nEdges/grid%nVertLevels
      print '(a,3f20.10)', 'min max mean edgeTangentVectors3 ',minval(edgeTangentVectors(3,1:nEdges)),maxval(edgeTangentVectors(3,1:nEdges)),sum(edgeTangentVectors(3,1:nEdges))/nEdges/grid%nVertLevels
endif
      do iEdge=1,nEdges
         ! mrp question: for this to be general across cores, we need to use nVertLevels
         do k=1,maxLevelEdgeBot(iEdge)
           do i=1,3
             do j=1,3
               ! outer produce at each edge:
               ! u_e n_e n_e* + v_e n_e \tilde{n}_e* 
               outerProductEdge3x3(i,j) = edgeNormalVectors(i,iEdge) &
                       *(  normalVelocity(k,iEdge) *edgeNormalVectors(j,iEdge) &
                         + tangentialVelocity(k,iEdge)*edgeTangentVectors(j,iEdge) &
                           )
             enddo
           enddo
           call mpas_matrix_3x3_to_sym6index(outerProductEdge3x3,outerProductEdge(:,k,iEdge))
if (1==2) then
  print '(a,i5,30f8.2)', '*********** iEdge, angleEdge',iEdge, angleEdge(iEdge)
  print '(a,30es10.2)', 'nor',edgeNormalVectors(:,iEdge)
  print '(a,30es10.2)', 'tan',edgeTangentVectors(:,iEdge)
  print '(a,30es10.2)', 'vel',normalVelocity(k,iEdge),tangentialVelocity(k,iEdge)
  print '(a,30es10.2)', '3x3',outerProductEdge3x3
  print '(a,30es10.2)', 'sm6', outerProductEdge(:,k,iEdge)
endif
         enddo
      enddo

if (1==2) then
      print '(a,3f20.10)', 'min max mean outerProductEdge 1  ',minval(outerProductEdge(1,1,1:nEdges)),maxval(outerProductEdge(1,1,1:nEdges)),sum(outerProductEdge(1,1,1:nEdges))/nEdges/grid%nVertLevels
      print '(a,3f20.10)', 'min max mean outerProductEdge 2  ',minval(outerProductEdge(2,1,1:nEdges)),maxval(outerProductEdge(2,1,1:nEdges)),sum(outerProductEdge(2,1,1:nEdges))/nEdges/grid%nVertLevels
      print '(a,3f20.10)', 'min max mean outerProductEdge 4  ',minval(outerProductEdge(4,1,1:nEdges)),maxval(outerProductEdge(4,1,1:nEdges)),sum(outerProductEdge(4,1,1:nEdges))/nEdges/grid%nVertLevels
print *, 'iCell, i, k, edgeSignOnCell(i,iCell),invAreaCell,outerProductEdge(:,k,iEdge)'
endif
      strainRate3DCell = 0.0
      do iCell = 1, nCells
!  print '(a,i5,30f8.2)', '***************************************************** iCell,xCell,yCell',iCell,xCell(iCell),yCell(iCell)
!  print '(a,30i5,30f8.2)', '** nCells,cellsOnCell',nCells,grid % cellsOnCell % array(:,iCell)
         invAreaCell = 1.0 / areaCell(iCell)
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            do k = 1, maxLevelCell(iCell)
               ! mrp edgeSignOnCell is to get outward unit normal on edgeNormalVectors
               ! minus sign in front is to match form on divergence, below
               strainRate3DCell(:,k,iCell) = strainRate3DCell(:,k,iCell) &
                 - edgeSignOnCell(i,iCell)*outerProductEdge(:,k,iEdge)*invAreaCell*dvEdge(iEdge) 
if (1==2) then
!if (iCell>=15.and.iCell<=17) then
  print '(a,2i5,30f8.2)', '*********** iEdge, edgeSignOnCell, angleEdge,xEdge,yEdge',iEdge, edgeSignOnCell(i,iCell), angleEdge(iEdge),xEdge(iEdge),yEdge(iEdge)
  print '(a,2i5,30f8.2)', 'cell1, cell2',grid%cellsOnEdge%array(:,iEdge)
  print '(a,30es10.2)', 'nor',edgeNormalVectors(:,iEdge)
  print '(a,30es10.2)', 'tan',edgeTangentVectors(:,iEdge)
  print '(a,30es10.2)', 'vel',normalVelocity(k,iEdge),tangentialVelocity(k,iEdge)
  print '(a,30es10.2)', 'sm6', outerProductEdge(:,k,iEdge)
  print '(a,30es10.2)', 'src', strainRate3DCell(:,k,iCell)
endif
            end do
         end do
if (1==2) then
  print '(a,30es10.2)', 'nor',edgeNormalVectors(:,iEdge)
  print '(a,30es10.2)', 'tan',edgeTangentVectors(:,iEdge)
  print '(a,30es10.2)', 'vel',normalVelocity(k,iEdge),tangentialVelocity(k,iEdge)
  print '(a,30es10.2)', '3x3',outerProductEdge3x3
  print '(a,30es10.2)', 'sm6', outerProductEdge(:,k,iEdge)
endif
      end do

   end subroutine mpas_strain_rate!}}}

!***********************************************************************
!
!  routine mpas_divergence_of_tensor
!
!> \brief   Computes divergence of the stress tensor
!> \author  Mark Petersen
!> \date    17 April 2013
!> \details 
!>  This routine computes the divergence of the stress tensor
!
!-----------------------------------------------------------------------

   subroutine mpas_divergence_of_tensor(grid, strainRate3DEdge, divTensor3DCell)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         strainRate3DEdge  !< Input: tensor at edge, 3D, in symmetric 6-index form

      type (mesh_type), intent(in) :: &
         grid          !< Input: grid information

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(out) :: &
         divTensor3DCell  !< Output: divergence of the tensor at cell center, 
                          !< as a 3-vector in x,y,z space

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iEdge, nEdges, iCell, nCells, i,j,k,p,q

      integer, dimension(:), pointer :: nEdgesOnCell, maxLevelCell
      integer, dimension(:,:), pointer :: edgesOnCell, edgeSignOnCell

      real (kind=RKIND) :: invAreaCell
      real (kind=RKIND), dimension(3) :: edgeNormalDotTensor
      real (kind=RKIND), dimension(3,3) :: strainRate3DEdge3x3
      real (kind=RKIND), dimension(:), pointer :: dvEdge, areaCell
      real (kind=RKIND), dimension(:,:), pointer :: edgeNormalVectors

      nEdges = grid % nEdges
      nCells = grid % nCells

      maxLevelCell      => grid % maxLevelCell % array
      edgesOnCell       => grid % edgesOnCell % array
      nEdgesOnCell      => grid % nEdgesOnCell % array
      dvEdge            => grid % dvEdge % array
      areaCell          => grid % areaCell % array
      edgeNormalVectors  => grid % edgeNormalVectors % array
      edgeSignOnCell    => grid % edgeSignOnCell % array

      divTensor3DCell(:,:,:) = 0.0
      do iCell = 1, nCells
         invAreaCell = 1.0 / areaCell(iCell)
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            do k = 1, maxLevelCell(iCell)
               call mpas_matrix_sym6index_to_3x3(strainRate3DEdge(:,k,iEdge),strainRate3DEdge3x3)
               edgeNormalDotTensor(:) = 0.0
               do q=1,3
                 do p=1,3
                   edgeNormalDotTensor(q) = edgeNormalDotTensor(q) + edgeNormalVectors(p,iEdge)*strainRate3DEdge3x3(p,q)
                 enddo
               enddo
               divTensor3DCell(:,k,iCell) = divTensor3DCell(:,k,iCell) &
                 - edgeSignOnCell(i,iCell) * dvEdge(iEdge) * edgeNormalDotTensor(:) * invAreaCell
            end do
         end do
      end do

   end subroutine mpas_divergence_of_tensor!}}}

!***********************************************************************
!
!  routine mpas_tensor_edge_3D_to_2D
!
!> \brief   Convert a 3D tensor to a 2D tensor, at an edge
!> \author  Mark Petersen
!> \date    17 April 2013
!> \details 
!>  Given a 3D tensor in symetric 6-index form, this routine rotates
!>  the tensor so that the 1-direction is towards the edge normal, and 
!>  the 2-direction is towards the edge tangent, and returns a 2D
!>  tensor in symmetric 3-index form. 
!
!-----------------------------------------------------------------------

   subroutine mpas_tensor_edge_3D_to_2D(strainRate3DEdge, strainRate2DEdge, grid)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         strainRate3DEdge  !< Input: strain rate tensor at edge, 3D, in symmetric 6-index form

      type (mesh_type), intent(in) :: &
         grid          !< Input: grid information

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(out) :: &
         strainRate2DEdge   !< Output: strain rate tensor at edge, 2D, in symmetric 3-index form

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iEdge, nEdges, i,j,k,p,q, nVertLevels
      integer, dimension(:), pointer :: maxLevelEdgeBot

      real (kind=RKIND), dimension(3) :: edgeVerticalVector
      real (kind=RKIND), dimension(3,3) :: rotationMatrix, strainRate3DEdge3x3, strainRate3DEdge3x3Rotated
      real (kind=RKIND), dimension(:,:), pointer :: edgeNormalVectors, edgeTangentVectors

      nEdges = grid % nEdges
      nVertLevels = grid % nVertLevels

      maxLevelEdgeBot   => grid % maxLevelEdgeBot % array
      edgeNormalVectors  => grid % edgeNormalVectors % array
      edgeTangentVectors => grid % edgeTangentVectors % array

      do iEdge=1,nEdges

         ! compute vertical vector at edge
         call mpas_cross_product_in_r3(edgeNormalVectors(:,iEdge),edgeTangentVectors(:,iEdge),edgeVerticalvector)

         rotationMatrix(:,1) = edgeNormalVectors(:,iEdge)
         rotationMatrix(:,2) = edgeTangentVectors(:,iEdge)
         rotationMatrix(:,3) = edgeVerticalvector

         ! mrp question: for this to be general across cores, we need to use nVertLevels
         do k=1,maxLevelEdgeBot(iEdge)
 
           call mpas_matrix_sym6index_to_3x3(strainRate3DEdge(:,k,iEdge),strainRate3DEdge3x3)

           ! If I were convinced everything is correct, could just compute i=1,2 and j=1,i
           do i=1,3
             do j=1,3
               do p=1,3
                 do q=1,3
                    ! I'm not sure if the star should be on the first or third matrix.
                    ! Here it is on the first
                    strainRate3DEdge3x3Rotated(i,j) = rotationMatrix(p,i)*strainRate3DEdge3x3(p,q)*rotationMatrix(q,j)
                 enddo
               enddo
             enddo
           enddo

           strainRate2DEdge(1,k,iEdge) = strainRate3DEdge3x3Rotated(1,1)
           strainRate2DEdge(2,k,iEdge) = strainRate3DEdge3x3Rotated(2,2)
           strainRate2DEdge(3,k,iEdge) = 0.5*(strainRate3DEdge3x3Rotated(1,2) + strainRate3DEdge3x3Rotated(2,1))
           ! mrp set land edges to zero?
         enddo

      enddo

   end subroutine mpas_tensor_edge_3D_to_2D!}}}

!***********************************************************************
!
!  routine mpas_tensor_edge_2D_to_3D
!
!> \brief   Convert a 3D tensor to a 2D tensor, at an edge
!> \author  Mark Petersen
!> \date    17 April 2013
!> \details 
!>  Given a 2D tensor in symetric 3-index form that is rotated such that
!>  the 1-direction is towards the edge normal, and 
!>  the 2-direction is towards the edge tangent, this routine rotates
!>  the tensor to x-y-z, and returns a 3D tensor in symetric 6-index form.
!
!-----------------------------------------------------------------------

   subroutine mpas_tensor_edge_2D_to_3D(strainRate2DEdge, strainRate3DEdge, grid)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         strainRate2DEdge   !< Input: strain rate tensor at edge, 2D, in symmetric 3-index form

      type (mesh_type), intent(in) :: &
         grid          !< Input: grid information

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(out) :: &
         strainRate3DEdge  !< Output: strain rate tensor at edge, 3D, in symmetric 6-index form

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iEdge, nEdges, i,j,k,p,q, nVertLevels
      integer, dimension(:), pointer :: maxLevelEdgeBot

      real (kind=RKIND), dimension(3) :: edgeVerticalVector
      real (kind=RKIND), dimension(3,3) :: rotationMatrix, strainRate3DEdge3x3, strainRate3DEdge3x3Rotated
      real (kind=RKIND), dimension(:,:), pointer :: edgeNormalVectors, edgeTangentVectors

      nEdges = grid % nEdges
      nVertLevels = grid % nVertLevels
      maxLevelEdgeBot   => grid % maxLevelEdgeBot % array
      edgeNormalVectors  => grid % edgeNormalVectors % array
      edgeTangentVectors => grid % edgeTangentVectors % array

      do iEdge=1,nEdges

         ! compute vertical vector at edge
         call mpas_cross_product_in_r3(edgeNormalVectors(:,iEdge),edgeTangentVectors(:,iEdge),edgeVerticalVector)

         rotationMatrix(:,1) = edgeNormalVectors(:,iEdge)
         rotationMatrix(:,2) = edgeTangentVectors(:,iEdge)
         rotationMatrix(:,3) = edgeVerticalVector

         ! mrp question: for this to be general across cores, we need to use nVertLevels
         do k=1,maxLevelEdgeBot(iEdge)

           strainRate3DEdge3x3Rotated = 0.0
           strainRate3DEdge3x3Rotated(1,1) = strainRate2DEdge(1,k,iEdge) 
           strainRate3DEdge3x3Rotated(2,2) = strainRate2DEdge(2,k,iEdge) 
           strainRate3DEdge3x3Rotated(1,2) = strainRate2DEdge(3,k,iEdge) 
           strainRate3DEdge3x3Rotated(2,1) = strainRate2DEdge(3,k,iEdge) 

           ! If I were convinced everything is correct, could just compute i=1,2 and j=1,i
           do i=1,3
             do j=1,3
               do p=1,3
                 do q=1,3
                    ! I'm not sure if the star should be on the first or third matrix.
                    ! Here it is on the third
                    strainRate3DEdge3x3(i,j) = rotationMatrix(i,p)*strainRate3DEdge3x3Rotated(p,q)*rotationMatrix(j,q)
                 enddo
               enddo
             enddo
           enddo
 
           call mpas_matrix_3x3_to_sym6index(strainRate3DEdge3x3,strainRate3DEdge(:,k,iEdge))

           ! mrp set land edges to zero?
         enddo

      enddo

   end subroutine mpas_tensor_edge_2D_to_3D!}}}

!***********************************************************************
!
!  routine mpas_test_tensor
!
!> \brief   Tests strain rate and tensor divergence operators
!> \author  Mark Petersen
!> \date    17 April 2013
!> \details 
!>  This routine tests strain rate and tensor divergence operators.
!
!-----------------------------------------------------------------------

   subroutine mpas_test_tensor(domain) !{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      type (block_type), pointer :: block

      integer :: nCells, nEdges, nVertices, nVertLevels, iCell, iEdge, k, i, p, minCell, maxCell

      real (kind=RKIND) :: xVelocity, yVelocity, cn, cs, r, theta, rot, f, g1, g2, fcos, pi2l, ld, &
        divSigma_x, divSigma_y
      real (kind=RKIND), dimension(:), pointer :: angleEdge, xCell, yCell, xEdge, yEdge
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity, tangentialVelocity
      real (kind=RKIND), dimension(:,:), pointer :: edgeNormalVectors, edgeTangentVectors
      real (kind=RKIND), dimension(:,:,:), pointer :: strainRate3DCell, strainRate3DEdge, divTensor3DCell

      real (kind=RKIND), dimension(:,:,:), pointer :: strainRate3DCellSolution, divTensor3DCellSolution
      real (kind=RKIND), dimension(:,:,:), pointer :: strainRate3DCellDiff, divTensor3DCellDiff, outerProductEdge
      character(len=16) :: test_function_type

    ! Initialize z-level grid variables from h, read in from input file.
    block => domain % blocklist
    do while (associated(block))

      nCells      = block % mesh % nCells
      nEdges      = block % mesh % nEdges
      nVertices   = block % mesh % nVertices
      nVertLevels = block % mesh % nVertLevels

      xCell  => block % mesh % xCell % array
      yCell  => block % mesh % yCell % array
      xEdge  => block % mesh % xEdge % array
      yEdge  => block % mesh % yEdge % array
      angleEdge => block % mesh % angleEdge % array

      edgeNormalVectors  => block % mesh % edgeNormalVectors % array
      edgeTangentVectors => block % mesh % edgeTangentVectors % array


      normalVelocity      => block % state % time_levs(1) % state % normalVelocity % array
      tangentialVelocity  => block % state % time_levs(1) % state % tangentialVelocity % array

      ! mrp convert these to scratch so they are not required for other cores
      strainRate3DCell    => block % state % time_levs(1) % state % strainRate3DCell % array
      strainRate3DEdge    => block % state % time_levs(1) % state % strainRate3DEdge % array
      divTensor3DCell     => block % state % time_levs(1) % state % divTensor3DCell % array
      outerProductEdge    => block % state % time_levs(1) % state % outerProductEdge % array

      strainRate3DCellSolution => block % state % time_levs(1) % state % strainRate3DCellSolution % array
      divTensor3DCellSolution => block % state % time_levs(1) % state % divTensor3DCellSolution % array
      strainRate3DCellDiff => block % state % time_levs(1) % state % strainRate3DCellDiff % array
      divTensor3DCellDiff => block % state % time_levs(1) % state % divTensor3DCellDiff % array

      ! create test functions for normalVelocity and tangentialVelocity
      normalVelocity = 0.0
      tangentialVelocity = 0.0
      strainRate3DCellSolution = 0.0
      divTensor3DCellSolution = 0.0

      cn = 15.0e4
      cs = 20.0e4
      rot = 1.0 !pii/4.0 !0.0*pii/2
      p = 2 ! power for polynomial function
!      test_function_type = 'constant'
!      test_function_type = 'linear_x'
!      test_function_type = 'linear_y'
!      test_function_type = 'linear_arb_rot'
      test_function_type = 'sin_arb_rot'
!      test_function_type = 'power_x'
!      test_function_type = 'power_y'
!      test_function_type = 'power_arb_rot'

      ld = 100.0e3  ! wavelength in meters
      pi2l = pii/ld
      g1 = cn*cos(rot) - cs*sin(rot)
      g2 = cn*sin(rot) + cs*cos(rot)

     if (test_function_type.eq.'constant') then
       ! linear function in x
        do iEdge = 1,nEdges
           xVelocity = cn
           yVelocity = cs
           normalVelocity(:,iEdge) = xVelocity*edgeNormalVectors(1,iEdge) + yVelocity*edgeNormalVectors(2,iEdge)
           tangentialVelocity(:,iEdge) = xVelocity*edgeTangentVectors(1,iEdge) + yVelocity*edgeTangentVectors(2,iEdge)
        enddo

        strainRate3DCellSolution = 0.0
        divTensor3DCellSolution = 0.0

     elseif (test_function_type.eq.'linear_x') then
       ! linear function in x
        do iEdge = 1,nEdges
           xVelocity = cn*xEdge(iEdge)
           yVelocity = cs*xEdge(iEdge)
           normalVelocity(:,iEdge) = xVelocity*edgeNormalVectors(1,iEdge) + yVelocity*edgeNormalVectors(2,iEdge)
           tangentialVelocity(:,iEdge) = xVelocity*edgeTangentVectors(1,iEdge) + yVelocity*edgeTangentVectors(2,iEdge)
        enddo

        strainRate3DCellSolution(1,:,:) = cn
        strainRate3DCellSolution(2,:,:) = 0.0
        strainRate3DCellSolution(4,:,:) = 0.5*cs

        divTensor3DCellSolution = 0.0

     elseif (test_function_type.eq.'linear_y') then

        ! linear function in y
        do iEdge = 1,nEdges
          xVelocity = -cs*yEdge(iEdge)
          yVelocity =  cn*yEdge(iEdge)
          normalVelocity(:,iEdge) = xVelocity*edgeNormalVectors(1,iEdge) + yVelocity*edgeNormalVectors(2,iEdge)
          tangentialVelocity(:,iEdge) = xVelocity*edgeTangentVectors(1,iEdge) + yVelocity*edgeTangentVectors(2,iEdge)
        enddo

        strainRate3DCellSolution(1,:,:) = 0.0
        strainRate3DCellSolution(2,:,:) = cn
        strainRate3DCellSolution(4,:,:) = -0.5*cs

        divTensor3DCellSolution = 0.0

     elseif (test_function_type.eq.'linear_arb_rot') then

        ! linear function, arbitrary rotation

        do iEdge = 1,nEdges
           ! arbitrary rotation
           r = sqrt(xEdge(iEdge)**2 + yEdge(iEdge)**2)
           theta = atan(yEdge(iEdge)/xEdge(iEdge))

           f = r*cos(theta-rot)
           xVelocity = f*g1
           yVelocity = f*g2

           normalVelocity(:,iEdge) = xVelocity*edgeNormalVectors(1,iEdge) + yVelocity*edgeNormalVectors(2,iEdge)
           tangentialVelocity(:,iEdge) = xVelocity*edgeTangentVectors(1,iEdge) + yVelocity*edgeTangentVectors(2,iEdge)

        enddo

        strainRate3DCellSolution(1,:,:) = cos(rot)*g1
        strainRate3DCellSolution(2,:,:) = sin(rot)*g2
        strainRate3DCellSolution(4,:,:) = 0.5*(cos(rot)*g2 + sin(rot)*g1)

        divTensor3DCellSolution = 0.0

     elseif (test_function_type.eq.'power_x') then

       ! power function in x: x^p
        do iEdge = 1,nEdges
           xVelocity = cn*xEdge(iEdge)**p
           yVelocity = cs*xEdge(iEdge)**p
           normalVelocity(:,iEdge) = xVelocity*edgeNormalVectors(1,iEdge) + yVelocity*edgeNormalVectors(2,iEdge)
           tangentialVelocity(:,iEdge) = xVelocity*edgeTangentVectors(1,iEdge) + yVelocity*edgeTangentVectors(2,iEdge)
        enddo

        do iCell = 1,nCells
           strainRate3DCellSolution(1,:,iCell) = cn    *p*xCell(iCell)**(p-1)
           strainRate3DCellSolution(2,:,iCell) = 0.0
           strainRate3DCellSolution(4,:,iCell) = 0.5*cs*p*xCell(iCell)**(p-1)

           divTensor3DCellSolution(1,:,iCell) = cn    *p*(p-1)*xCell(iCell)**(p-2)
           divTensor3DCellSolution(2,:,iCell) = 0.5*cs*p*(p-1)*xCell(iCell)**(p-2)  
        end do

     elseif (test_function_type.eq.'power_y') then

        ! power function in y: y^n
        do iEdge = 1,nEdges
           xVelocity = -cs*yEdge(iEdge)**p
           yVelocity =  cn*yEdge(iEdge)**p
           normalVelocity(:,iEdge) = xVelocity*edgeNormalVectors(1,iEdge) + yVelocity*edgeNormalVectors(2,iEdge)
           tangentialVelocity(:,iEdge) = xVelocity*edgeTangentVectors(1,iEdge) + yVelocity*edgeTangentVectors(2,iEdge)
        enddo

        do iCell = 1,nCells
           strainRate3DCellSolution(1,:,iCell) = 0.0
           strainRate3DCellSolution(2,:,iCell) =  cn    *p*yCell(iCell)**(p-1)
           strainRate3DCellSolution(4,:,iCell) = -0.5*cs*p*yCell(iCell)**(p-1)

           divTensor3DCellSolution(1,:,iCell) = -0.5*cs*p*(p-1)*yCell(iCell)**(p-2)  
           divTensor3DCellSolution(2,:,iCell) =  cn    *p*(p-1)*yCell(iCell)**(p-2)
        end do

     elseif (test_function_type.eq.'power_arb_rot') then

        ! power function, arbitrary rotation

        do iEdge = 1,nEdges
           ! arbitrary rotation
           r = sqrt(xEdge(iEdge)**2 + yEdge(iEdge)**2)
           theta = atan(yEdge(iEdge)/xEdge(iEdge))
           f = r*cos(theta-rot)

           xVelocity = g1*f**p
           yVelocity = g2*f**p
           normalVelocity(:,iEdge) = xVelocity*edgeNormalVectors(1,iEdge) + yVelocity*edgeNormalVectors(2,iEdge)
           tangentialVelocity(:,iEdge) = xVelocity*edgeTangentVectors(1,iEdge) + yVelocity*edgeTangentVectors(2,iEdge)
        enddo

        do iCell = 1,nCells
           r = sqrt(xCell(iCell)**2 + yCell(iCell)**2)
           theta = atan(yCell(iCell)/xCell(iCell))
           f = r*cos(theta-rot)

           strainRate3DCellSolution(1,:,iCell) = p *f**(p-1) *cos(rot)*g1
           strainRate3DCellSolution(2,:,iCell) = p *f**(p-1) *sin(rot)*g2
           strainRate3DCellSolution(4,:,iCell) = p *f**(p-1) *(cos(rot)*g2+sin(rot)*g1)/2.0

           divTensor3DCellSolution(1,:,iCell) = p*(p-1)*f**(p-2) *(cos(rot)**2*g1 + 0.5*(sin(rot)**2*g1 + sin(rot)*cos(rot)*g2) )
           divTensor3DCellSolution(2,:,iCell) = p*(p-1)*f**(p-2) *(sin(rot)**2*g2 + 0.5*(cos(rot)**2*g2 + sin(rot)*cos(rot)*g1) )
        end do

     elseif (test_function_type.eq.'sin_arb_rot') then

        ! sin function, arbitrary rotation

        do iEdge = 1,nEdges
           ! arbitrary rotation
           r = sqrt(xEdge(iEdge)**2 + yEdge(iEdge)**2)
           theta = atan(yEdge(iEdge)/xEdge(iEdge))
           f = sin(pi2l*r*cos(theta-rot))

           xVelocity = f*g1
           yVelocity = f*g2
           normalVelocity(:,iEdge) = xVelocity*edgeNormalVectors(1,iEdge) + yVelocity*edgeNormalVectors(2,iEdge)
           tangentialVelocity(:,iEdge) = xVelocity*edgeTangentVectors(1,iEdge) + yVelocity*edgeTangentVectors(2,iEdge)
        enddo

        do iCell = 1,nCells
           r = sqrt(xCell(iCell)**2 + yCell(iCell)**2)
           theta = atan(yCell(iCell)/xCell(iCell))
           f = sin(pi2l*r*cos(theta-rot))
           fcos = cos(pi2l*r*cos(theta-rot))

           strainRate3DCellSolution(1,:,iCell) = pi2l*fcos*cos(rot)*g1
           strainRate3DCellSolution(2,:,iCell) = pi2l*fcos*sin(rot)*g2
           strainRate3DCellSolution(4,:,iCell) = pi2l*fcos*(cos(rot)*g2+sin(rot)*g1)/2.0

           divTensor3DCellSolution(1,:,iCell) = -pi2l**2*f*(cos(rot)**2*g1 + 0.5*(sin(rot)**2*g1 + sin(rot)*cos(rot)*g2) )
           divTensor3DCellSolution(2,:,iCell) = -pi2l**2*f*(sin(rot)**2*g2 + 0.5*(cos(rot)**2*g2 + sin(rot)*cos(rot)*g1) )
        end do

     else
       print *, 'bad choice of test_function_type:',test_function_type
       stop
     endif

!     call mpas_tangential_velocity(block % mesh, normalVelocity, tangentialVelocity)

     call mpas_strain_rate(block % mesh, normalVelocity, tangentialVelocity, strainRate3DCell, outerProductEdge)

     call mpas_matrix_cell_to_edge(block % mesh, strainRate3DCell, strainRate3DEdge)

     call mpas_divergence_of_tensor(block % mesh, strainRate3DEdge, divTensor3DCell)


minCell=17
maxCell=17 
if (1==1) then
!      print '(a,3es20.10)', 'min max mean normalVelocity      ',minval(normalVelocity(:,1:nEdges)),maxval(normalVelocity(:,1:nEdges)),sum(normalVelocity(:,1:nEdges))/nEdges/nVertLevels
!      print '(a,3es20.10)', 'min max mean tangentialVelocity  ',minval(tangentialVelocity(:,1:nEdges)),maxval(tangentialVelocity(:,1:nEdges)),sum(tangentialVelocity(:,1:nEdges))/nEdges/nVertLevels
      print '(a,3es20.10)', 'min max mean strainRate3DCell1   ',minval(strainRate3DCell(1,:,minCell:maxCell)),maxval(strainRate3DCell(1,:,minCell:maxCell)),sum(strainRate3DCell(1,:,minCell:maxCell))/(maxCell-minCell+1)/nVertLevels
      print '(a,3es20.10)', 'min max mean strainRate3DCellSol1',minval(strainRate3DCellSolution(1,:,minCell:maxCell)),maxval(strainRate3DCellSolution(1,:,minCell:maxCell)),sum(strainRate3DCellSolution(1,:,minCell:maxCell))/(maxCell-minCell+1)/nVertLevels
      print '(/,a,3es20.10)', 'min max mean strainRate3DCell2   ',minval(strainRate3DCell(2,:,minCell:maxCell)),maxval(strainRate3DCell(2,:,minCell:maxCell)),sum(strainRate3DCell(2,:,minCell:maxCell))/(maxCell-minCell+1)/nVertLevels
      print '(a,3es20.10)', 'min max mean strainRate3DCell2Sol',minval(strainRate3DCellSolution(2,:,minCell:maxCell)),maxval(strainRate3DCellSolution(2,:,minCell:maxCell)),sum(strainRate3DCellSolution(2,:,minCell:maxCell))/(maxCell-minCell+1)/nVertLevels
      print '(/,a,3es20.10)', 'min max mean strainRate3DCell4   ',minval(strainRate3DCell(4,:,minCell:maxCell)),maxval(strainRate3DCell(4,:,minCell:maxCell)),sum(strainRate3DCell(4,:,minCell:maxCell))/(maxCell-minCell+1)/nVertLevels
      print '(a,3es20.10)', 'min max mean strainRate3DCell4Sol',minval(strainRate3DCellSolution(4,:,minCell:maxCell)),maxval(strainRate3DCellSolution(4,:,minCell:maxCell)),sum(strainRate3DCellSolution(4,:,minCell:maxCell))/(maxCell-minCell+1)/nVertLevels
      print '(/,a,3es20.10)', 'min max mean divTensor3DCell1    ',minval(divTensor3DCell(1,:,minCell:maxCell)),maxval(divTensor3DCell(1,:,minCell:maxCell)),sum(divTensor3DCell(1,:,minCell:maxCell))/(maxCell-minCell+1)/nVertLevels
      print '(a,3es20.10)',   'min max mean divTensor3DCellSol1 ',minval(divTensor3DCellSolution(1,:,minCell:maxCell)),maxval(divTensor3DCellSolution(1,:,minCell:maxCell)),sum(divTensor3DCellSolution(1,:,minCell:maxCell))/(maxCell-minCell+1)/nVertLevels
!      print '(/,a,3es20.10)', 'min max mean strainRate3DEdge1   ',minval(strainRate3DEdge(1,:,minCell:maxCell)),maxval(strainRate3DEdge(1,:,minCell:maxCell)),sum(strainRate3DEdge(1,:,minCell:maxCell))/(maxCell-minCell+1)/nVertLevels
!      print '(a,3es20.10)', 'min max mean strainRate3DEdge2   ',minval(strainRate3DEdge(2,:,minCell:maxCell)),maxval(strainRate3DEdge(2,:,minCell:maxCell)),sum(strainRate3DEdge(2,:,minCell:maxCell))/(maxCell-minCell+1)/nVertLevels
!      print '(a,3es20.10)', 'min max mean strainRate3DEdge3   ',minval(strainRate3DEdge(3,:,minCell:maxCell)),maxval(strainRate3DEdge(3,:,minCell:maxCell)),sum(strainRate3DEdge(3,:,minCell:maxCell))/(maxCell-minCell+1)/nVertLevels

      print '(/,a,3es20.10)', 'min max mean divTensor3DCell2    ',minval(divTensor3DCell(2,:,minCell:maxCell)),maxval(divTensor3DCell(2,:,minCell:maxCell)),sum(divTensor3DCell(2,:,minCell:maxCell))/(maxCell-minCell+1)/nVertLevels
      print '(a,3es20.10)',   'min max mean divTensor3DCellSol2 ',minval(divTensor3DCellSolution(2,:,minCell:maxCell)),maxval(divTensor3DCellSolution(2,:,minCell:maxCell)),sum(divTensor3DCellSolution(2,:,minCell:maxCell))/(maxCell-minCell+1)/nVertLevels
      print '(/,a,3es20.10)', 'min max mean divTensor3DCell3    ',minval(divTensor3DCell(3,:,minCell:maxCell)),maxval(divTensor3DCell(3,:,minCell:maxCell)),sum(divTensor3DCell(3,:,minCell:maxCell))/(maxCell-minCell+1)/nVertLevels
      print '(a,3es20.10)',   'min max mean divTensor3DCellSol3 ',minval(divTensor3DCellSolution(3,:,minCell:maxCell)),maxval(divTensor3DCellSolution(3,:,minCell:maxCell)),sum(divTensor3DCellSolution(3,:,minCell:maxCell))/(maxCell-minCell+1)/nVertLevels
endif
 
      block % state % time_levs(1) % state % strainRate3DCellDiff % array = abs(strainRate3DCell- strainRate3DCellSolution)
      block % state % time_levs(1) % state % divTensor3DCellDiff % array = abs(divTensor3DCell- divTensor3DCellSolution)

      print '(/,a,3es20.10)', 'max abs diff strainRate3DCell1',maxval(abs(strainRate3DCell(1,:,minCell:maxCell) - strainRate3DCellSolution(1,:,minCell:maxCell)))
      print '(a,3es20.10)', 'max abs diff strainRate3DCell2',maxval(abs(strainRate3DCell(2,:,minCell:maxCell) - strainRate3DCellSolution(2,:,minCell:maxCell)))
      print '(a,3es20.10)', 'max abs diff strainRate3DCell4',maxval(abs(strainRate3DCell(4,:,minCell:maxCell) - strainRate3DCellSolution(4,:,minCell:maxCell)))
      print '(a,3es20.10)', 'max abs diff divTensor3DCell  ',maxval(abs(divTensor3DCell(:,:,minCell:maxCell) - divTensor3DCellSolution(:,:,minCell:maxCell)))

        block => block % next
      end do

   end subroutine mpas_test_tensor!}}}

!***********************************************************************
!
!  routine mpas_tensor_operations_init
!
!> \brief   Initializes flags used within tendency routines.
!> \author  Mark Petersen
!> \date    17 April 2013
!> \details 
!>  This routine initializes flags related to quantities computed within
!>  other tendency routines.
!
!-----------------------------------------------------------------------
!    subroutine mpas_tensor_operations_init(err)!{{{
!        integer, intent(out) :: err !< Output: Error flag

!    end subroutine mpas_tensor_operations_init!}}}

!***********************************************************************

end module mpas_tensor_operations


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker

