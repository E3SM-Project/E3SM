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
             mpas_tensor_edge_R3_to_2D, &
             mpas_tensor_edge_2D_to_R3, &
             mpas_tensor_R3_to_LatLon, &
             mpas_tensor_LatLon_to_R3, &
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
!>  derivative.  Output is an R3 strain rate tensor in 6-index format.
!
!-----------------------------------------------------------------------

   subroutine mpas_strain_rate(grid, normalVelocity, tangentialVelocity, &
      strainRateR3Cell, outerProductEdge)!{{{

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
         strainRateR3Cell   !< Output: strain rate tensor at cell center, R3, in symmetric 6-index form

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
                       *(  normalVelocity(k,iEdge)    *edgeNormalVectors(j,iEdge) &
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
      strainRateR3Cell = 0.0
      do iCell = 1, nCells
!  print '(a,i5,30f8.2)', '***************************************************** iCell,xCell,yCell',iCell,xCell(iCell),yCell(iCell)
!  print '(a,30i5,30f8.2)', '** nCells,cellsOnCell',nCells,grid % cellsOnCell % array(:,iCell)
         invAreaCell = 1.0 / areaCell(iCell)
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            do k = 1, maxLevelCell(iCell)
               ! mrp edgeSignOnCell is to get outward unit normal on edgeNormalVectors
               ! minus sign in front is to match form on divergence, below
               strainRateR3Cell(:,k,iCell) = strainRateR3Cell(:,k,iCell) &
                 - edgeSignOnCell(i,iCell)*outerProductEdge(:,k,iEdge)*invAreaCell*dvEdge(iEdge) 
if (1==2) then
!if (iCell>=15.and.iCell<=17) then
  print '(a,2i5,30f8.2)', '*********** iEdge, edgeSignOnCell, angleEdge,xEdge,yEdge',iEdge, edgeSignOnCell(i,iCell), angleEdge(iEdge),xEdge(iEdge),yEdge(iEdge)
  print '(a,2i5,30f8.2)', 'cell1, cell2',grid%cellsOnEdge%array(:,iEdge)
  print '(a,30es10.2)', 'nor',edgeNormalVectors(:,iEdge)
  print '(a,30es10.2)', 'tan',edgeTangentVectors(:,iEdge)
  print '(a,30es10.2)', 'vel',normalVelocity(k,iEdge),tangentialVelocity(k,iEdge)
  print '(a,30es10.2)', 'sm6', outerProductEdge(:,k,iEdge)
  print '(a,30es10.2)', 'src', strainRateR3Cell(:,k,iCell)
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

   subroutine mpas_divergence_of_tensor(grid, strainRateR3Edge, divTensorR3Cell)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         strainRateR3Edge  !< Input: tensor at edge, R3, in symmetric 6-index form

      type (mesh_type), intent(in) :: &
         grid          !< Input: grid information

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(out) :: &
         divTensorR3Cell  !< Output: divergence of the tensor at cell center, 
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
      real (kind=RKIND), dimension(3,3) :: strainRateR3Edge3x3
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

      divTensorR3Cell(:,:,:) = 0.0
      do iCell = 1, nCells
         invAreaCell = 1.0 / areaCell(iCell)
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            do k = 1, maxLevelCell(iCell)
               call mpas_matrix_sym6index_to_3x3(strainRateR3Edge(:,k,iEdge),strainRateR3Edge3x3)
               edgeNormalDotTensor(:) = 0.0
               do q=1,3
                 do p=1,3
                   edgeNormalDotTensor(q) = edgeNormalDotTensor(q) + edgeNormalVectors(p,iEdge)*strainRateR3Edge3x3(p,q)
                 enddo
               enddo
               divTensorR3Cell(:,k,iCell) = divTensorR3Cell(:,k,iCell) &
                 - edgeSignOnCell(i,iCell) * dvEdge(iEdge) * edgeNormalDotTensor(:) * invAreaCell
            end do
         end do
      end do

   end subroutine mpas_divergence_of_tensor!}}}

!***********************************************************************
!
!  routine mpas_tensor_edge_R3_to_2D
!
!> \brief   Convert an R3 tensor to a 2D tensor, at an edge
!> \author  Mark Petersen
!> \date    17 April 2013
!> \details 
!>  Given an R3 tensor in symetric 6-index form, this routine rotates
!>  the tensor so that the 1-direction is towards the edge normal, and 
!>  the 2-direction is towards the edge tangent, and returns a 2D
!>  tensor in symmetric 3-index form. 
!
!-----------------------------------------------------------------------

   subroutine mpas_tensor_edge_R3_to_2D(strainRateR3Edge, strainRate2DEdge, grid)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         strainRateR3Edge  !< Input: strain rate tensor at edge, R3, in symmetric 6-index form

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
      real (kind=RKIND), dimension(3,3) :: rotationMatrix, strainRateR3Edge3x3, strainRateR3Edge3x3Rotated
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
 
           call mpas_matrix_sym6index_to_3x3(strainRateR3Edge(:,k,iEdge),strainRateR3Edge3x3)

           ! If I were convinced everything is correct, could just compute i=1,2 and j=1,i
           strainRateR3Edge3x3Rotated = 0.0
           do i=1,3
             do j=1,3
               do p=1,3
                 do q=1,3
                    ! I'm not sure if the star should be on the first or third matrix.
                    ! Here it is on the first
                    strainRateR3Edge3x3Rotated(i,j) = strainRateR3Edge3x3Rotated(i,j) + rotationMatrix(p,i)*strainRateR3Edge3x3(p,q)*rotationMatrix(q,j)
                 enddo
               enddo
             enddo
           enddo

           strainRate2DEdge(1,k,iEdge) = strainRateR3Edge3x3Rotated(1,1)
           strainRate2DEdge(2,k,iEdge) = strainRateR3Edge3x3Rotated(2,2)
           strainRate2DEdge(3,k,iEdge) = 0.5*(strainRateR3Edge3x3Rotated(1,2) + strainRateR3Edge3x3Rotated(2,1))
           ! mrp set land edges to zero?
         enddo

      enddo

   end subroutine mpas_tensor_edge_R3_to_2D!}}}

!***********************************************************************
!
!  routine mpas_tensor_edge_2D_to_R3
!
!> \brief   Convert a 2D tensor to a tensor in R3, at an edge
!> \author  Mark Petersen
!> \date    17 April 2013
!> \details 
!>  Given a 2D tensor in symetric 3-index form that is rotated such that
!>  the 1-direction is towards the edge normal, and 
!>  the 2-direction is towards the edge tangent, this routine rotates
!>  the tensor to R3, and returns an R3 tensor in symetric 6-index form.
!
!-----------------------------------------------------------------------

   subroutine mpas_tensor_edge_2D_to_R3(strainRate2DEdge, strainRateR3Edge, grid)!{{{

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
         strainRateR3Edge  !< Output: strain rate tensor at edge, R3, in symmetric 6-index form

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iEdge, nEdges, i,j,k,p,q, nVertLevels
      integer, dimension(:), pointer :: maxLevelEdgeBot

      real (kind=RKIND), dimension(3) :: edgeVerticalVector
      real (kind=RKIND), dimension(3,3) :: rotationMatrix, strainRateR3Edge3x3, strainRateR3Edge3x3Rotated
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

           strainRateR3Edge3x3Rotated = 0.0
           strainRateR3Edge3x3Rotated(1,1) = strainRate2DEdge(1,k,iEdge) 
           strainRateR3Edge3x3Rotated(2,2) = strainRate2DEdge(2,k,iEdge) 
           strainRateR3Edge3x3Rotated(1,2) = strainRate2DEdge(3,k,iEdge) 
           strainRateR3Edge3x3Rotated(2,1) = strainRate2DEdge(3,k,iEdge) 

           strainRateR3Edge3x3 = 0.0
           do i=1,3
             do j=1,3
               do p=1,3
                 do q=1,3
                    ! I'm not sure if the star should be on the first or third matrix.
                    ! Here it is on the third
                    strainRateR3Edge3x3(i,j) = strainRateR3Edge3x3(i,j) + rotationMatrix(i,p)*strainRateR3Edge3x3Rotated(p,q)*rotationMatrix(j,q)
                 enddo
               enddo
             enddo
           enddo
 
           call mpas_matrix_3x3_to_sym6index(strainRateR3Edge3x3,strainRateR3Edge(:,k,iEdge))

           ! mrp set land edges to zero?
         enddo

      enddo

   end subroutine mpas_tensor_edge_2D_to_R3!}}}



!***********************************************************************
!
!  routine mpas_tensor_LatLon_to_R3
!
!> \brief   Convert an R3 tensor to a 2D tensor
!> \author  Mark Petersen
!> \date    17 April 2013
!> \details 
!>  Given a 2D tensor in symetric 3-index form that is rotated such that
!>  the 1-direction is zonal
!>  the 2-direction is meridional, this routine rotates
!>  the tensor to R3, and returns an R3 tensor in symetric 6-index form.
!
!-----------------------------------------------------------------------

   subroutine mpas_tensor_LatLon_to_R3(strainRateLatLon, strainRateR3, lat, lon)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(3), intent(in) :: &
         strainRateLatLon   !< Input: strain rate tensor, 2D, in symmetric 3-index form

      real (kind=RKIND), intent(in) :: &
         lat, &!< Input: latitude,  in radians, ranging [-pi,pi]
         lon   !< Input: longitude, in radians, ranging [0,2*pi]

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(6), intent(out) :: &
         strainRateR3  !< Output: strain rate tensor, R3, in symmetric 6-index form

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: i,j,k,p,q, nVertLevels

      real (kind=RKIND) :: determinant
      real (kind=RKIND), dimension(3) :: zonalUnitVector, meridionalUnitVector, verticalUnitVector
      real (kind=RKIND), dimension(3,3) :: rotationMatrix, strainRateR3_3x3, strainRateLatLon3x3, A

      call mpas_zonal_meridional_vectors(lat, lon, zonalUnitVector, meridionalUnitVector, verticalUnitVector)

      rotationMatrix(:,1) = zonalUnitVector
      rotationMatrix(:,2) = meridionalUnitVector
      rotationMatrix(:,3) = verticalUnitVector

      ! just for checking
      A(1,1) = -2
      A(1,2) = 2
      A(1,3) =3
      A(2,1) =-1
      A(2,2) =1
      A(2,3) =3
      A(3,1) = 2
      A(3,2) = 0
      A(3,3) = -1
      A = rotationMatrix
!      determinant = A(1,1)*A(2,2)*A(3,3) + A(1,2)*A(2,3)*A(3,1) + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1) - A(1,2)*A(2,1)*A(3,3) - A(1,1)*A(2,3)*A(3,2)
!print '(/,a,100es10.2)', 'determinant', determinant

      strainRateLatLon3x3 = 0.0
      strainRateLatLon3x3(1,1) = strainRateLatLon(1)
      strainRateLatLon3x3(2,2) = strainRateLatLon(2)
      strainRateLatLon3x3(1,2) = strainRateLatLon(3)
      strainRateLatLon3x3(2,1) = strainRateLatLon(3)

!print '(/,a,100es10.2)', 'lat,lon, zonalUnitVector, meridionalUnitVector, verticalUnitVector',lat,lon, zonalUnitVector, meridionalUnitVector, verticalUnitVector
!print '(a,100es10.2)', 'strainRateLatLon3x3',strainRateLatLon3x3
!print '(a,100es10.2)', 'rotationMatrix',rotationMatrix
      strainRateR3_3x3 = 0.0
      do i=1,3
        do j=1,3
          do p=1,3
            do q=1,3
               strainRateR3_3x3(i,j) = strainRateR3_3x3(i,j) + rotationMatrix(i,p)*strainRateLatLon3x3(p,q)*rotationMatrix(j,q)
            enddo
          enddo
        enddo
      enddo
!print '(a,100es10.2)', 'strainRateR3_3x3',strainRateR3_3x3
 
      call mpas_matrix_3x3_to_sym6index(strainRateR3_3x3,strainRateR3)
!print '(a,100es10.2)', 'strainRateR3',strainRateR3

   end subroutine mpas_tensor_LatLon_to_R3!}}}


!***********************************************************************
!
!  routine mpas_tensor_3DLatLon_to_R3
!
!> \brief   Convert an R3 tensor to a 2D tensor
!> \author  Mark Petersen
!> \date    17 April 2013
!> \details 
!>  Given a 2D tensor in symetric 3-index form that is rotated such that
!>  the 1-direction is zonal
!>  the 2-direction is meridional, this routine rotates
!>  the tensor to R3, and returns an R3 tensor in symetric 6-index form.
!
!-----------------------------------------------------------------------

   subroutine mpas_tensor_3DLatLon_to_R3(strainRateLatLon, strainRateR3, lat, lon)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(6), intent(in) :: &
         strainRateLatLon   !< Input: latlon strain rate tensor, 3D, in symmetric 6-index form

      real (kind=RKIND), intent(in) :: &
         lat, &!< Input: latitude,  in radians, ranging [-pi,pi]
         lon   !< Input: longitude, in radians, ranging [0,2*pi]

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(6), intent(out) :: &
         strainRateR3  !< Output: strain rate tensor, R3, in symmetric 6-index form

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: i,j,k,p,q, nVertLevels

      real (kind=RKIND) :: determinant
      real (kind=RKIND), dimension(3) :: zonalUnitVector, meridionalUnitVector, verticalUnitVector
      real (kind=RKIND), dimension(3,3) :: rotationMatrix, strainRateR3_3x3, strainRateLatLon3x3, A

      call mpas_zonal_meridional_vectors(lat, lon, zonalUnitVector, meridionalUnitVector, verticalUnitVector)

      rotationMatrix(:,1) = zonalUnitVector
      rotationMatrix(:,2) = meridionalUnitVector
      rotationMatrix(:,3) = verticalUnitVector

      call mpas_matrix_sym6index_to_3x3(strainRateLatLon,strainRateLatLon3x3)

!print '(/,a,100es10.2)', 'lat,lon, zonalUnitVector, meridionalUnitVector, verticalUnitVector',lat,lon, zonalUnitVector, meridionalUnitVector, verticalUnitVector
!print '(a,100es10.2)', 'strainRateLatLon3x3',strainRateLatLon3x3
!print '(a,100es10.2)', 'rotationMatrix',rotationMatrix
      strainRateR3_3x3 = 0.0
      do i=1,3
        do j=1,3
          do p=1,3
            do q=1,3
               strainRateR3_3x3(i,j) = strainRateR3_3x3(i,j) + rotationMatrix(i,p)*strainRateLatLon3x3(p,q)*rotationMatrix(j,q)
            enddo
          enddo
        enddo
      enddo
!print '(a,100es10.2)', 'strainRateR3_3x3',strainRateR3_3x3
 
      call mpas_matrix_3x3_to_sym6index(strainRateR3_3x3,strainRateR3)
!print '(a,100es10.2)', 'strainRateR3',strainRateR3

   end subroutine mpas_tensor_3DLatLon_to_R3!}}}

!***********************************************************************
!
!  routine mpas_tensor_R3_to_LatLon
!
!> \brief   Convert an R3 tensor to a 2D tensor
!> \author  Mark Petersen
!> \date    17 April 2013
!> \details 
!>  Given a 2D tensor in symetric 3-index form that is rotated such that
!>  the 1-direction is zonal
!>  the 2-direction is meridional, this routine rotates
!>  the tensor to R3, and returns an R3 tensor in symetric 6-index form.
!
!-----------------------------------------------------------------------

   subroutine mpas_tensor_R3_to_LatLon(strainRateR3, strainRateLatLon, lat, lon)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(6), intent(in) :: &
         strainRateR3  !< Input: strain rate tensor at, R3, in symmetric 6-index form

      real (kind=RKIND), intent(in) :: &
         lat, &!< Input: latitude,  in radians, ranging [-pi,pi]
         lon   !< Input: longitude, in radians, ranging [0,2*pi]

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(3), intent(out) :: &
         strainRateLatLon   !< Output: strain rate tensor, 2D, in symmetric 3-index form

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: i,j,k,p,q

      real (kind=RKIND), dimension(3) :: zonalUnitVector, meridionalUnitVector, verticalUnitVector
      real (kind=RKIND), dimension(3,3) :: rotationMatrix, strainRateR3_3x3, strainRateLatLon3x3

      call mpas_zonal_meridional_vectors(lat, lon, zonalUnitVector, meridionalUnitVector, verticalUnitVector)

      rotationMatrix(:,1) = zonalUnitVector
      rotationMatrix(:,2) = meridionalUnitVector
      rotationMatrix(:,3) = verticalUnitVector

      call mpas_matrix_sym6index_to_3x3(strainRateR3,strainRateR3_3x3)

      strainRateLatLon3x3 = 0
      do i=1,3
        do j=1,3
          do p=1,3
            do q=1,3
               strainRateLatLon3x3(i,j) = strainRateLatLon3x3(i,j) + rotationMatrix(p,i)*strainRateR3_3x3(p,q)*rotationMatrix(q,j)
            enddo
          enddo
        enddo
      enddo
 
      strainRateLatLon(1) = strainRateLatLon3x3(1,1)
      strainRateLatLon(2) = strainRateLatLon3x3(2,2)
      strainRateLatLon(3) = 0.5*(strainRateLatLon3x3(1,2)+strainRateLatLon3x3(2,1))

   end subroutine mpas_tensor_R3_to_LatLon!}}}


!***********************************************************************
!
!  routine mpas_tensor_R3_to_3DLatLon
!
!> \brief   Convert an R3 tensor to a 2D tensor
!> \author  Mark Petersen
!> \date    17 April 2013
!> \details 
!>  Given a 2D tensor in symetric 3-index form that is rotated such that
!>  the 1-direction is zonal
!>  the 2-direction is meridional, this routine rotates
!>  the tensor to R3, and returns an R3 tensor in symetric 6-index form.
!
!-----------------------------------------------------------------------

   subroutine mpas_tensor_R3_to_3DLatLon(strainRateR3, strainRate3DLatLon, lat, lon)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(6), intent(in) :: &
         strainRateR3  !< Input: strain rate tensor at, R3, in symmetric 6-index form

      real (kind=RKIND), intent(in) :: &
         lat, &!< Input: latitude,  in radians, ranging [-pi,pi]
         lon   !< Input: longitude, in radians, ranging [0,2*pi]

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(6), intent(out) :: &
         strainRate3DLatLon   !< Output: strain rate tensor, 3D lat-lon coord, in symmetric 6-index form

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: i,j,k,p,q

      real (kind=RKIND), dimension(3) :: zonalUnitVector, meridionalUnitVector, verticalUnitVector
      real (kind=RKIND), dimension(3,3) :: rotationMatrix, strainRateR3_3x3, strainRateLatLon3x3

      call mpas_zonal_meridional_vectors(lat, lon, zonalUnitVector, meridionalUnitVector, verticalUnitVector)

      rotationMatrix(:,1) = zonalUnitVector
      rotationMatrix(:,2) = meridionalUnitVector
      rotationMatrix(:,3) = verticalUnitVector

      call mpas_matrix_sym6index_to_3x3(strainRateR3,strainRateR3_3x3)

      strainRateLatLon3x3 = 0
      do i=1,3
        do j=1,3
          do p=1,3
            do q=1,3
               strainRateLatLon3x3(i,j) = strainRateLatLon3x3(i,j) + rotationMatrix(p,i)*strainRateR3_3x3(p,q)*rotationMatrix(q,j)
            enddo
          enddo
        enddo
      enddo

      call mpas_matrix_3x3_to_sym6index(strainRateLatLon3x3,strainRate3DLatLon)

   end subroutine mpas_tensor_R3_to_3DLatLon!}}}


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
        divSigma_x, divSigma_y, kLon, kLat, &
        lon, lat, sin_kLon_lon, cos_kLon_lon, sin_kLat_lat, cos_kLat_lat, sin_lat, sec_lat, tan_lat, &
        velocityZonal, velocityMeridional
      real (kind=RKIND), dimension(6) :: strainRateLatLonCellSolution
      real (kind=RKIND), dimension(:), pointer :: angleEdge, xCell, yCell, xEdge, yEdge, latCell, lonCell, latEdge, lonEdge
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity, tangentialVelocity
      real (kind=RKIND), dimension(:,:), pointer :: edgeNormalVectors, edgeTangentVectors
      real (kind=RKIND), dimension(:,:,:), pointer :: strainRateR3Cell, strainRateR3Edge, divTensorR3Cell

      real (kind=RKIND), dimension(:,:,:), pointer :: strainRateR3CellSolution, strainRate2DCellSolution, divTensorR3CellSolution, strainRate3DLatLonCellSolution
      real (kind=RKIND), dimension(:,:,:), pointer :: strainRateR3CellDiff, divTensorR3CellDiff, outerProductEdge
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
      latCell=> block % mesh % latCell % array
      lonCell=> block % mesh % lonCell % array
      latEdge=> block % mesh % latEdge % array
      lonEdge=> block % mesh % lonEdge % array
      xEdge  => block % mesh % xEdge % array
      yEdge  => block % mesh % yEdge % array
      angleEdge => block % mesh % angleEdge % array

      edgeNormalVectors  => block % mesh % edgeNormalVectors % array
      edgeTangentVectors => block % mesh % edgeTangentVectors % array


      normalVelocity      => block % state % time_levs(1) % state % normalVelocity % array
      tangentialVelocity  => block % state % time_levs(1) % state % tangentialVelocity % array

      ! mrp convert these to scratch so they are not required for other cores
      strainRateR3Cell    => block % state % time_levs(1) % state % strainRateR3Cell % array
      strainRateR3Edge    => block % state % time_levs(1) % state % strainRateR3Edge % array
      divTensorR3Cell     => block % state % time_levs(1) % state % divTensorR3Cell % array
      outerProductEdge    => block % state % time_levs(1) % state % outerProductEdge % array

      strainRate2DCellSolution => block % state % time_levs(1) % state % strainRate2DCellSolution % array
      strainRate3DLatLonCellSolution => block % state % time_levs(1) % state % strainRate3DLatLonCellSolution % array
      strainRateR3CellSolution => block % state % time_levs(1) % state % strainRateR3CellSolution % array
      divTensorR3CellSolution => block % state % time_levs(1) % state % divTensorR3CellSolution % array
      strainRateR3CellDiff => block % state % time_levs(1) % state % strainRateR3CellDiff % array
      divTensorR3CellDiff => block % state % time_levs(1) % state % divTensorR3CellDiff % array

      ! create test functions for normalVelocity and tangentialVelocity
      normalVelocity = 0.0
      tangentialVelocity = 0.0
      strainRateR3CellSolution = 0.0
      divTensorR3CellSolution = 0.0

      cn = 15.0e4
      cs = 20.0e4
      rot = 1.0 !pii/4.0 !0.0*pii/2
      p = 2 ! power for polynomial function
!      test_function_type = 'constant'
!      test_function_type = 'linear_x'
!      test_function_type = 'linear_y'
!      test_function_type = 'linear_arb_rot'
!      test_function_type = 'sin_arb_rot'
!      test_function_type = 'power_x'
!      test_function_type = 'power_y'
!      test_function_type = 'power_arb_rot'
!      test_function_type = 'spherical_sinx'
!      test_function_type = 'spherical_u1'
      test_function_type = 'spherical_v1'
!     test_function_type = 'spherical_solid'
!      test_function_type = 'spherical'

      ld = 100.0e3  ! wavelength in meters
      pi2l = pii/ld
      g1 = cn*cos(rot) - cs*sin(rot)
      g2 = cn*sin(rot) + cs*cos(rot)

      R = block % mesh % sphere_radius

     if (test_function_type.eq.'constant') then
       ! linear function in x
        do iEdge = 1,nEdges
           xVelocity = cn
           yVelocity = cs
           normalVelocity(:,iEdge) = xVelocity*edgeNormalVectors(1,iEdge) + yVelocity*edgeNormalVectors(2,iEdge)
           tangentialVelocity(:,iEdge) = xVelocity*edgeTangentVectors(1,iEdge) + yVelocity*edgeTangentVectors(2,iEdge)
        enddo

        strainRateR3CellSolution = 0.0
        divTensorR3CellSolution = 0.0

     elseif (test_function_type.eq.'linear_x') then
       ! linear function in x
        do iEdge = 1,nEdges
           xVelocity = cn*xEdge(iEdge)
           yVelocity = cs*xEdge(iEdge)
           normalVelocity(:,iEdge) = xVelocity*edgeNormalVectors(1,iEdge) + yVelocity*edgeNormalVectors(2,iEdge)
           tangentialVelocity(:,iEdge) = xVelocity*edgeTangentVectors(1,iEdge) + yVelocity*edgeTangentVectors(2,iEdge)
        enddo

        strainRateR3CellSolution(1,:,:) = cn
        strainRateR3CellSolution(2,:,:) = 0.0
        strainRateR3CellSolution(4,:,:) = 0.5*cs

        divTensorR3CellSolution = 0.0

     elseif (test_function_type.eq.'linear_y') then

        ! linear function in y
        do iEdge = 1,nEdges
          xVelocity = -cs*yEdge(iEdge)
          yVelocity =  cn*yEdge(iEdge)
          normalVelocity(:,iEdge) = xVelocity*edgeNormalVectors(1,iEdge) + yVelocity*edgeNormalVectors(2,iEdge)
          tangentialVelocity(:,iEdge) = xVelocity*edgeTangentVectors(1,iEdge) + yVelocity*edgeTangentVectors(2,iEdge)
        enddo

        strainRateR3CellSolution(1,:,:) = 0.0
        strainRateR3CellSolution(2,:,:) = cn
        strainRateR3CellSolution(4,:,:) = -0.5*cs

        divTensorR3CellSolution = 0.0

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

        strainRateR3CellSolution(1,:,:) = cos(rot)*g1
        strainRateR3CellSolution(2,:,:) = sin(rot)*g2
        strainRateR3CellSolution(4,:,:) = 0.5*(cos(rot)*g2 + sin(rot)*g1)

        divTensorR3CellSolution = 0.0

     elseif (test_function_type.eq.'power_x') then

       ! power function in x: x^p
        do iEdge = 1,nEdges
           xVelocity = cn*xEdge(iEdge)**p
           yVelocity = cs*xEdge(iEdge)**p
           normalVelocity(:,iEdge) = xVelocity*edgeNormalVectors(1,iEdge) + yVelocity*edgeNormalVectors(2,iEdge)
           tangentialVelocity(:,iEdge) = xVelocity*edgeTangentVectors(1,iEdge) + yVelocity*edgeTangentVectors(2,iEdge)
        enddo

        do iCell = 1,nCells
           strainRateR3CellSolution(1,:,iCell) = cn    *p*xCell(iCell)**(p-1)
           strainRateR3CellSolution(2,:,iCell) = 0.0
           strainRateR3CellSolution(4,:,iCell) = 0.5*cs*p*xCell(iCell)**(p-1)

           divTensorR3CellSolution(1,:,iCell) = cn    *p*(p-1)*xCell(iCell)**(p-2)
           divTensorR3CellSolution(2,:,iCell) = 0.5*cs*p*(p-1)*xCell(iCell)**(p-2)  
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
           strainRateR3CellSolution(1,:,iCell) = 0.0
           strainRateR3CellSolution(2,:,iCell) =  cn    *p*yCell(iCell)**(p-1)
           strainRateR3CellSolution(4,:,iCell) = -0.5*cs*p*yCell(iCell)**(p-1)

           divTensorR3CellSolution(1,:,iCell) = -0.5*cs*p*(p-1)*yCell(iCell)**(p-2)  
           divTensorR3CellSolution(2,:,iCell) =  cn    *p*(p-1)*yCell(iCell)**(p-2)
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

           strainRateR3CellSolution(1,:,iCell) = p *f**(p-1) *cos(rot)*g1
           strainRateR3CellSolution(2,:,iCell) = p *f**(p-1) *sin(rot)*g2
           strainRateR3CellSolution(4,:,iCell) = p *f**(p-1) *(cos(rot)*g2+sin(rot)*g1)/2.0

           divTensorR3CellSolution(1,:,iCell) = p*(p-1)*f**(p-2) *(cos(rot)**2*g1 + 0.5*(sin(rot)**2*g1 + sin(rot)*cos(rot)*g2) )
           divTensorR3CellSolution(2,:,iCell) = p*(p-1)*f**(p-2) *(sin(rot)**2*g2 + 0.5*(cos(rot)**2*g2 + sin(rot)*cos(rot)*g1) )
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

           strainRateR3CellSolution(1,:,iCell) = pi2l*fcos*cos(rot)*g1
           strainRateR3CellSolution(2,:,iCell) = pi2l*fcos*sin(rot)*g2
           strainRateR3CellSolution(4,:,iCell) = pi2l*fcos*(cos(rot)*g2+sin(rot)*g1)/2.0

           divTensorR3CellSolution(1,:,iCell) = -pi2l**2*f*(cos(rot)**2*g1 + 0.5*(sin(rot)**2*g1 + sin(rot)*cos(rot)*g2) )
           divTensorR3CellSolution(2,:,iCell) = -pi2l**2*f*(sin(rot)**2*g2 + 0.5*(cos(rot)**2*g2 + sin(rot)*cos(rot)*g1) )
        end do

     elseif (test_function_type.eq.'spherical_solid') then

        do iEdge = 1,nEdges
           lon = lonEdge(iEdge)
           lat = latEdge(iEdge)

           velocityZonal      = cos(lat)
           velocityMeridional = 0.0
           normalVelocity(:,iEdge) = velocityZonal*cos(angleEdge(iEdge)) + velocityMeridional*sin(angleEdge(iEdge))
        enddo

        call mpas_tangential_velocity(block % mesh, normalVelocity, tangentialVelocity)

        do iCell = 1,nCells

           lon = lonCell(iCell)
           lat = latCell(iCell)

           block % state % time_levs(1) % state % streamFunctionCell     % array(:,iCell) = 0.0 ! actually -lat/2 - sin(lat)**2/4
           block % state % time_levs(1) % state % zonalVelocityCell      % array(:,iCell) = cos(lat)
           block % state % time_levs(1) % state % meridionalVelocityCell % array(:,iCell) = 0.0

           strainRate3DLatLonCellSolution(1,1,iCell) = & ! Elonlon
                 0.0
           strainRate3DLatLonCellSolution(2,1,iCell) = & ! Elatlat
                 0.0
           strainRate3DLatLonCellSolution(3,1,iCell) = & ! Err
                 0.0
           strainRate3DLatLonCellSolution(4,1,iCell) = & ! Elonlat
                 0.0
           strainRate3DLatLonCellSolution(5,1,iCell) = & ! Elatr
                 0.0
           strainRate3DLatLonCellSolution(6,1,iCell) = & ! Elonr
                 - 1.0/R/2.0*cos(lat)

           call mpas_tensor_3DLatLon_to_R3(strainRate3DLatLonCellSolution(:,1,iCell), strainRateR3CellSolution(:,1,iCell), lat, lon)
        enddo

        divTensorR3CellSolution = 0.0

     elseif (test_function_type.eq.'spherical_u1') then

        kLon = 1.0
        kLat = 0.0

        do iEdge = 1,nEdges
           lon = lonEdge(iEdge)
           lat = latEdge(iEdge)

           if (abs(lat)>0.8*pii/2) then
             velocityZonal      = 0.0
             velocityMeridional = 0.0
           else
             velocityZonal      = 1.0
             velocityMeridional = 0.0
           endif
           normalVelocity(:,iEdge) = velocityZonal*cos(angleEdge(iEdge)) + velocityMeridional*sin(angleEdge(iEdge))
        enddo

        call mpas_tangential_velocity(block % mesh, normalVelocity, tangentialVelocity)

        do iCell = 1,nCells

           lon = lonCell(iCell)
           lat = latCell(iCell)

           sin_kLon_lon = sin(kLon*lon)
           cos_kLon_lon = cos(kLon*lon)
           sin_kLat_lat = sin(kLat*lat)
           cos_kLat_lat = cos(kLat*lat)
           sin_lat = sin(lat)
           sec_lat = 1/cos(lat)
           tan_lat = sin_lat*sec_lat

           if (abs(lat)>0.8*pii/2) then
             block % state % time_levs(1) % state % streamFunctionCell     % array(:,iCell) = 0.0
             block % state % time_levs(1) % state % zonalVelocityCell      % array(:,iCell) = 0.0
             block % state % time_levs(1) % state % meridionalVelocityCell % array(:,iCell) = 0.0
           else
             block % state % time_levs(1) % state % streamFunctionCell     % array(:,iCell) = 0.0
             block % state % time_levs(1) % state % zonalVelocityCell      % array(:,iCell) = 1.0 
             block % state % time_levs(1) % state % meridionalVelocityCell % array(:,iCell) = 0.0
           endif

           strainRate3DLatLonCellSolution(1,1,iCell) = & ! Elonlon
                 0.0
           strainRate3DLatLonCellSolution(2,1,iCell) = & ! Elatlat
                 0.0
           strainRate3DLatLonCellSolution(3,1,iCell) = & ! Err
                 0.0
           strainRate3DLatLonCellSolution(4,1,iCell) = & ! Elonlat
                 -0.5/R*tan_lat
           strainRate3DLatLonCellSolution(5,1,iCell) = & ! Elatr
                 0.0
           strainRate3DLatLonCellSolution(6,1,iCell) = & ! Elonr
                 - 1.0/R/2.0

           call mpas_tensor_3DLatLon_to_R3(strainRate3DLatLonCellSolution(:,1,iCell), strainRateR3CellSolution(:,1,iCell), lat, lon)

!           divTensorR3CellSolution(1,:,iCell) = 
!           divTensorR3CellSolution(2,:,iCell) = 

        end do


     elseif (test_function_type.eq.'spherical_v1') then

        kLon = 1.0
        kLat = 0.0

        do iEdge = 1,nEdges
           lon = lonEdge(iEdge)
           lat = latEdge(iEdge)

           if (abs(lat)>0.8*pii/2) then
             velocityZonal      = 0.0
             velocityMeridional = 0.0
           else
             velocityZonal      = 0.0
             velocityMeridional = 1.0
           endif
           normalVelocity(:,iEdge) = velocityZonal*cos(angleEdge(iEdge)) + velocityMeridional*sin(angleEdge(iEdge))
        enddo

        call mpas_tangential_velocity(block % mesh, normalVelocity, tangentialVelocity)

        do iCell = 1,nCells

           lon = lonCell(iCell)
           lat = latCell(iCell)

           sin_kLon_lon = sin(kLon*lon)
           cos_kLon_lon = cos(kLon*lon)
           sin_kLat_lat = sin(kLat*lat)
           cos_kLat_lat = cos(kLat*lat)
           sin_lat = sin(lat)
           sec_lat = 1/cos(lat)
           tan_lat = sin_lat*sec_lat

           if (abs(lat)>0.8*pii/2) then
             block % state % time_levs(1) % state % streamFunctionCell     % array(:,iCell) = 0.0
             block % state % time_levs(1) % state % zonalVelocityCell      % array(:,iCell) = 0.0
             block % state % time_levs(1) % state % meridionalVelocityCell % array(:,iCell) = 0.0
           else
             block % state % time_levs(1) % state % streamFunctionCell     % array(:,iCell) = 0.0
             block % state % time_levs(1) % state % zonalVelocityCell      % array(:,iCell) = 0.0 
             block % state % time_levs(1) % state % meridionalVelocityCell % array(:,iCell) = 1.0
           endif

           strainRate3DLatLonCellSolution(1,1,iCell) = & ! Elonlon
                 1.0/R*tan_lat
           strainRate3DLatLonCellSolution(2,1,iCell) = & ! Elatlat
                 0.0
           strainRate3DLatLonCellSolution(3,1,iCell) = & ! Err
                 0.0
           strainRate3DLatLonCellSolution(4,1,iCell) = & ! Elonlat
                 0.0
           strainRate3DLatLonCellSolution(5,1,iCell) = & ! Elatr
                 - 1.0/R/2.0
           strainRate3DLatLonCellSolution(6,1,iCell) = & ! Elonr
                 0.0

           call mpas_tensor_3DLatLon_to_R3(strainRate3DLatLonCellSolution(:,1,iCell), strainRateR3CellSolution(:,1,iCell), lat, lon)

!           divTensorR3CellSolution(1,:,iCell) = 
!           divTensorR3CellSolution(2,:,iCell) = 

        end do


     elseif (test_function_type.eq.'spherical_sinx') then

        kLon = 1.0
        kLat = 0.0

        do iEdge = 1,nEdges
           lon = lonEdge(iEdge)
           lat = latEdge(iEdge)

           if (abs(lat)>0.8*pii/2) then
             velocityZonal      = 0.0
             velocityMeridional = 0.0
           else
             velocityZonal      = sin(kLon *lon) 
             velocityMeridional = 0.0
           endif
           normalVelocity(:,iEdge) = velocityZonal*cos(angleEdge(iEdge)) + velocityMeridional*sin(angleEdge(iEdge))
        enddo

        call mpas_tangential_velocity(block % mesh, normalVelocity, tangentialVelocity)

        do iCell = 1,nCells

           lon = lonCell(iCell)
           lat = latCell(iCell)

           sin_kLon_lon = sin(kLon*lon)
           cos_kLon_lon = cos(kLon*lon)
           sin_kLat_lat = sin(kLat*lat)
           cos_kLat_lat = cos(kLat*lat)
           sin_lat = sin(lat)
           sec_lat = 1/cos(lat)
           tan_lat = sin_lat*sec_lat

           if (abs(lat)>0.8*pii/2) then
             block % state % time_levs(1) % state % streamFunctionCell     % array(:,iCell) = 0.0
             block % state % time_levs(1) % state % zonalVelocityCell      % array(:,iCell) = 0.0
             block % state % time_levs(1) % state % meridionalVelocityCell % array(:,iCell) = 0.0
           else
             block % state % time_levs(1) % state % streamFunctionCell     % array(:,iCell) = 0.0
             block % state % time_levs(1) % state % zonalVelocityCell      % array(:,iCell) = sin_kLon_lon 
             block % state % time_levs(1) % state % meridionalVelocityCell % array(:,iCell) = 0.0
           endif

           strainRate2DCellSolution(1,:,iCell) = & ! Elonlon
                 1/R *sec_lat * kLon *cos_kLon_lon
           strainRate2DCellSolution(2,:,iCell) = & ! Elatlat
                 0.0
           strainRate2DCellSolution(3,:,iCell) = & ! Elonlat
                -1/R/2 *tan_lat *sin_kLon_lon 

           call mpas_tensor_LatLon_to_R3(strainRate2DCellSolution(:,1,iCell), strainRateR3CellSolution(:,1,iCell), lat, lon)

!           divTensorR3CellSolution(1,:,iCell) = 
!           divTensorR3CellSolution(2,:,iCell) = 

        end do

     elseif (test_function_type.eq.'spherical') then

        kLon = 0.0
        kLat = 1.0

        do iEdge = 1,nEdges
           lon = lonEdge(iEdge)
           lat = latEdge(iEdge)

           if (abs(lat)>0.8*pii/2) then
             velocityZonal      = 0.0
             velocityMeridional = 0.0
           else
             velocityZonal      = kLat *cos(kLon *lon) *sin(kLat *lat)/cos(lat) 
             velocityMeridional = kLon *sin(kLon *lon) *cos(kLat *lat) 
           endif
           normalVelocity(:,iEdge) = velocityZonal*cos(angleEdge(iEdge)) + velocityMeridional*sin(angleEdge(iEdge))
        enddo

        call mpas_tangential_velocity(block % mesh, normalVelocity, tangentialVelocity)

        do iCell = 1,nCells

           lon = lonCell(iCell)
           lat = latCell(iCell)

           sin_kLon_lon = sin(kLon*lon)
           cos_kLon_lon = cos(kLon*lon)
           sin_kLat_lat = sin(kLat*lat)
           cos_kLat_lat = cos(kLat*lat)
           sin_lat = sin(lat)
           sec_lat = 1/cos(lat)
           tan_lat = sin_lat*sec_lat

           if (abs(lat)>0.8*pii/2) then
             block % state % time_levs(1) % state % streamFunctionCell     % array(:,iCell) = 0.0
             block % state % time_levs(1) % state % zonalVelocityCell      % array(:,iCell) = 0.0
             block % state % time_levs(1) % state % meridionalVelocityCell % array(:,iCell) = 0.0
           else
             block % state % time_levs(1) % state % streamFunctionCell     % array(:,iCell) =       cos_kLon_lon *cos_kLat_lat
             block % state % time_levs(1) % state % zonalVelocityCell      % array(:,iCell) = kLat *cos_kLon_lon *sin_kLat_lat *sec_lat
             block % state % time_levs(1) % state % meridionalVelocityCell % array(:,iCell) = kLon *cos_kLat_lat *sin_kLon_lon
           endif

           strainRate2DCellSolution(1,:,iCell) = & ! Elonlon
                 1/R *sec_lat * (  kLon      *sin_lat *sin_kLon_lon *cos_kLat_lat  &
                                 - kLon*kLat *sec_lat *sin_kLon_lon *sin_kLat_lat)
           strainRate2DCellSolution(2,:,iCell) = & ! Elatlat
                 1/R *kLon *kLat *sin_kLon_lon *sin_kLat_lat
           strainRate2DCellSolution(3,:,iCell) = & ! Elonlat
                 1/R/2 * (  kLon**2          *sec_lat *cos_kLon_lon *cos_kLat_lat  &
                          - kLat**2          *sec_lat *cos_kLon_lon *cos_kLat_lat  &
                          - 2 *kLat *tan_lat *sec_lat *cos_kLon_lon *sin_kLat_lat)

           call mpas_tensor_LatLon_to_R3(strainRate2DCellSolution(:,1,iCell), strainRateR3CellSolution(:,1,iCell), lat, lon)

!           divElon = 1/R *sec_lat ((*sin_lat (kLon**2 *cos_kLon_lon *cos_kLat_lat *sec_lat - *kLat**2 *cos_kLon_lon *cos_kLat_lat *sec_lat &
!   + 2 *kLat *cos_kLon_lon *tan_lat *sec_lat *sin_kLat_lat))/(2 R) + 1/(2 R)
!   *sin(lat) (-kLon**2 *cos_kLon_lon *cos_kLat_lat *tan_lat *sec_lat + 3 *kLat**2 *cos_kLon_lon *cos_kLat_lat *tan_lat *sec_lat - kLon**2 *kLat *cos_kLon_lon *sec_lat *sin_kLat_lat + 
!      *kLat**3 *cos_kLon_lon *sec_lat *sin_kLat_lat -  2 *kLat *cos_kLon_lon *tan_lat**2 *sec_lat *sin_kLat_lat - 2 *kLat *cos_kLon_lon *sec_lat**3 *sin_kLat_lat) &
!    - ( *tan_lat *sec_lat (kLon *sin_lat *cos_kLat_lat *sin_kLon_lon + kLon *kLat *sec_lat *sin_kLon_lon *sin_kLat_lat))/R &
!     + 1/R *sec_lat (kLon *kLat**2 *cos_kLat_lat *sec_lat *sin_kLon_lon - 
!      kLon *cos_kLat_lat *sin_kLon_lon *sin(lat) - 
!      kLon *kLat *sin_lat *sin_kLon_lon *sin_kLat_lat - 
!      kLon *kLat *tan_lat *sec_lat *sin_kLon_lon *sin_kLat_lat))

!           divElat = 1/R**2 *sec_lat (-(kLon *kLat**2 *cos_kLat_lat *sin_kLon_lon *sin(lat))/R &
! - (  kLon *kLat *sin_lat *sin_kLon_lon *sin_kLat_lat)/R + (-kLon**3 *cos_kLat_lat *sec_lat *sin_kLon_lon + 
!    kLon *kLat**2 *cos_kLat_lat *sec_lat *sin_kLon_lon - 
!    2 kLon *kLat *tan_lat *sec_lat *sin_kLon_lon *sin_kLat_lat)/(2 R))

!           divTensorR3CellSolution(1,:,iCell) = 
!           divTensorR3CellSolution(2,:,iCell) = 

        end do

     else
       print *, 'bad choice of test_function_type:',test_function_type
       stop
     endif


     call mpas_strain_rate(block % mesh, normalVelocity, tangentialVelocity, strainRateR3Cell, outerProductEdge)

     call mpas_matrix_cell_to_edge(block % mesh, strainRateR3Cell, strainRateR3Edge)

     call mpas_divergence_of_tensor(block % mesh, strainRateR3Edge, divTensorR3Cell)

     do iCell = 1,nCells
        call mpas_tensor_R3_to_3DLatLon(strainRateR3Cell(:,1,iCell), &
            block % state % time_levs(1) % state % strainRate3DLatLonCell % array(:,1,iCell), latCell(iCell), lonCell(iCell))
     enddo

minCell=17
maxCell=17 
if (1==1) then
!      print '(a,3es20.10)', 'min max mean normalVelocity      ',minval(normalVelocity(:,1:nEdges)),maxval(normalVelocity(:,1:nEdges)),sum(normalVelocity(:,1:nEdges))/nEdges/nVertLevels
!      print '(a,3es20.10)', 'min max mean tangentialVelocity  ',minval(tangentialVelocity(:,1:nEdges)),maxval(tangentialVelocity(:,1:nEdges)),sum(tangentialVelocity(:,1:nEdges))/nEdges/nVertLevels
      print '(a,3es20.10)', 'min max mean strainRateR3Cell1   ',minval(strainRateR3Cell(1,:,minCell:maxCell)),maxval(strainRateR3Cell(1,:,minCell:maxCell)),sum(strainRateR3Cell(1,:,minCell:maxCell))/(maxCell-minCell+1)/nVertLevels
      print '(a,3es20.10)', 'min max mean strainRateR3CellSol1',minval(strainRateR3CellSolution(1,:,minCell:maxCell)),maxval(strainRateR3CellSolution(1,:,minCell:maxCell)),sum(strainRateR3CellSolution(1,:,minCell:maxCell))/(maxCell-minCell+1)/nVertLevels
      print '(/,a,3es20.10)', 'min max mean strainRateR3Cell2   ',minval(strainRateR3Cell(2,:,minCell:maxCell)),maxval(strainRateR3Cell(2,:,minCell:maxCell)),sum(strainRateR3Cell(2,:,minCell:maxCell))/(maxCell-minCell+1)/nVertLevels
      print '(a,3es20.10)', 'min max mean strainRateR3Cell2Sol',minval(strainRateR3CellSolution(2,:,minCell:maxCell)),maxval(strainRateR3CellSolution(2,:,minCell:maxCell)),sum(strainRateR3CellSolution(2,:,minCell:maxCell))/(maxCell-minCell+1)/nVertLevels
      print '(/,a,3es20.10)', 'min max mean strainRateR3Cell4   ',minval(strainRateR3Cell(4,:,minCell:maxCell)),maxval(strainRateR3Cell(4,:,minCell:maxCell)),sum(strainRateR3Cell(4,:,minCell:maxCell))/(maxCell-minCell+1)/nVertLevels
      print '(a,3es20.10)', 'min max mean strainRateR3Cell4Sol',minval(strainRateR3CellSolution(4,:,minCell:maxCell)),maxval(strainRateR3CellSolution(4,:,minCell:maxCell)),sum(strainRateR3CellSolution(4,:,minCell:maxCell))/(maxCell-minCell+1)/nVertLevels
      print '(/,a,3es20.10)', 'min max mean divTensorR3Cell1    ',minval(divTensorR3Cell(1,:,minCell:maxCell)),maxval(divTensorR3Cell(1,:,minCell:maxCell)),sum(divTensorR3Cell(1,:,minCell:maxCell))/(maxCell-minCell+1)/nVertLevels
      print '(a,3es20.10)',   'min max mean divTensorR3CellSol1 ',minval(divTensorR3CellSolution(1,:,minCell:maxCell)),maxval(divTensorR3CellSolution(1,:,minCell:maxCell)),sum(divTensorR3CellSolution(1,:,minCell:maxCell))/(maxCell-minCell+1)/nVertLevels
!      print '(/,a,3es20.10)', 'min max mean strainRateR3Edge1   ',minval(strainRateR3Edge(1,:,minCell:maxCell)),maxval(strainRateR3Edge(1,:,minCell:maxCell)),sum(strainRateR3Edge(1,:,minCell:maxCell))/(maxCell-minCell+1)/nVertLevels
!      print '(a,3es20.10)', 'min max mean strainRateR3Edge2   ',minval(strainRateR3Edge(2,:,minCell:maxCell)),maxval(strainRateR3Edge(2,:,minCell:maxCell)),sum(strainRateR3Edge(2,:,minCell:maxCell))/(maxCell-minCell+1)/nVertLevels
!      print '(a,3es20.10)', 'min max mean strainRateR3Edge3   ',minval(strainRateR3Edge(3,:,minCell:maxCell)),maxval(strainRateR3Edge(3,:,minCell:maxCell)),sum(strainRateR3Edge(3,:,minCell:maxCell))/(maxCell-minCell+1)/nVertLevels

      print '(/,a,3es20.10)', 'min max mean divTensorR3Cell2    ',minval(divTensorR3Cell(2,:,minCell:maxCell)),maxval(divTensorR3Cell(2,:,minCell:maxCell)),sum(divTensorR3Cell(2,:,minCell:maxCell))/(maxCell-minCell+1)/nVertLevels
      print '(a,3es20.10)',   'min max mean divTensorR3CellSol2 ',minval(divTensorR3CellSolution(2,:,minCell:maxCell)),maxval(divTensorR3CellSolution(2,:,minCell:maxCell)),sum(divTensorR3CellSolution(2,:,minCell:maxCell))/(maxCell-minCell+1)/nVertLevels
      print '(/,a,3es20.10)', 'min max mean divTensorR3Cell3    ',minval(divTensorR3Cell(3,:,minCell:maxCell)),maxval(divTensorR3Cell(3,:,minCell:maxCell)),sum(divTensorR3Cell(3,:,minCell:maxCell))/(maxCell-minCell+1)/nVertLevels
      print '(a,3es20.10)',   'min max mean divTensorR3CellSol3 ',minval(divTensorR3CellSolution(3,:,minCell:maxCell)),maxval(divTensorR3CellSolution(3,:,minCell:maxCell)),sum(divTensorR3CellSolution(3,:,minCell:maxCell))/(maxCell-minCell+1)/nVertLevels
endif
 
      block % state % time_levs(1) % state % strainRateR3CellDiff % array = abs(strainRateR3Cell- strainRateR3CellSolution)
      block % state % time_levs(1) % state % divTensorR3CellDiff % array = abs(divTensorR3Cell- divTensorR3CellSolution)

      print '(/,a,3es20.10)', 'max abs diff strainRateR3Cell1',maxval(abs(strainRateR3Cell(1,:,minCell:maxCell) - strainRateR3CellSolution(1,:,minCell:maxCell)))
      print '(a,3es20.10)', 'max abs diff strainRateR3Cell2',maxval(abs(strainRateR3Cell(2,:,minCell:maxCell) - strainRateR3CellSolution(2,:,minCell:maxCell)))
      print '(a,3es20.10)', 'max abs diff strainRateR3Cell4',maxval(abs(strainRateR3Cell(4,:,minCell:maxCell) - strainRateR3CellSolution(4,:,minCell:maxCell)))
      print '(a,3es20.10)', 'max abs diff divTensorR3Cell  ',maxval(abs(divTensorR3Cell(:,:,minCell:maxCell) - divTensorR3CellSolution(:,:,minCell:maxCell)))

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

