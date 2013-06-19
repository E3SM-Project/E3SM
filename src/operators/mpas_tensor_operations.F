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
!> \date    April 2013
!> \details
!>  This module contains the routines for computing
!>  the strain rate tensor, the divergence of a tensor,
!>  and a testing routine to verify these work properly.
!
!-----------------------------------------------------------------------

module mpas_tensor_operations

   use mpas_grid_types
   use mpas_constants
   use mpas_configure
   use mpas_vector_operations
   use mpas_matrix_operations
   use mpas_dmpar

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

   public :: mpas_strain_rate_R3Cell, &
             mpas_divergence_of_tensor_R3Cell, &
             mpas_tensor_edge_R3_to_2D, &
             mpas_tensor_edge_2D_to_R3, &
             mpas_tensor_R3_to_LonLat, &
             mpas_tensor_LonLat_to_R3, &
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
!  routine mpas_strain_rate_R3Cell
!
!> \brief   Computes strain rate at cell centers, in R3
!> \author  Mark Petersen
!> \date    April 2013
!> \details 
!>  This routine computes the strain rate at cell centers using the weak 
!>  derivative.  Output is an R3 strain rate tensor in 6-index format.
!
!-----------------------------------------------------------------------

   subroutine mpas_strain_rate_R3Cell(normalVelocity, tangentialVelocity, &
      strainRateR3Cell, outerProductEdge, &
      grid, edgeSignOnCell, edgeTangentVectors)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         edgeTangentVectors,   &!< Input: unit vector tangent to an edge
         normalVelocity,      &!< Input: Horizontal velocity normal to edge
         tangentialVelocity    !< Input: Horizontal velocity tangent to edge

      integer, dimension(:,:), intent(in) :: &
         edgeSignOnCell        !< Input: Direction of vector connecting cells

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

      integer :: iEdge, nEdges, iCell, nCells, i,j,k, nVertLevels

      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnCell

      real (kind=RKIND) :: invAreaCell
      real (kind=RKIND), dimension(3,3) :: outerProductEdge3x3
      real (kind=RKIND), dimension(:), pointer :: dvEdge, areaCell, angleEdge
      real (kind=RKIND), dimension(:,:), pointer :: edgeNormalVectors

      nEdges = grid % nEdges
      nCells = grid % nCells
      nVertLevels = grid % nVertLevels

      nEdgesOnCell      => grid % nEdgesOnCell % array
      edgesOnCell       => grid % edgesOnCell % array
      dvEdge            => grid % dvEdge % array
      angleEdge         => grid % angleEdge % array
      areaCell          => grid % areaCell % array
      edgeNormalVectors  => grid % edgeNormalVectors % array

      do iEdge=1,nEdges
         do k=1,nVertLevels
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
         enddo
      enddo

      strainRateR3Cell = 0.0
      do iCell = 1, nCells
         invAreaCell = 1.0 / areaCell(iCell)
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            do k = 1, nVertLevels
               ! edgeSignOnCell is to get outward unit normal on edgeNormalVectors
               ! minus sign in front is to match form on divergence operator
               strainRateR3Cell(:,k,iCell) = strainRateR3Cell(:,k,iCell) &
                 - edgeSignOnCell(i,iCell)*outerProductEdge(:,k,iEdge)*invAreaCell*dvEdge(iEdge) 
            end do
         end do
      end do

   end subroutine mpas_strain_rate_R3Cell!}}}

!***********************************************************************
!
!  routine mpas_divergence_of_tensor_R3Cell
!
!> \brief   Computes divergence of the stress tensor
!> \author  Mark Petersen
!> \date    April 2013
!> \details 
!>  This routine computes the divergence of the stress tensor
!
!-----------------------------------------------------------------------

   subroutine mpas_divergence_of_tensor_R3Cell(grid, edgeSignOnCell, strainRateR3Edge, divTensorR3Cell)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         strainRateR3Edge  !< Input: tensor at edge, R3, in symmetric 6-index form

      type (mesh_type), intent(in) :: &
         grid          !< Input: grid information

      integer, dimension(:,:), intent(in) :: &
         edgeSignOnCell        !< Input: Direction of vector connecting cells

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

      integer :: iEdge, iCell, nCells, i,j,k,p,q, nVertLevels

      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnCell

      real (kind=RKIND) :: invAreaCell
      real (kind=RKIND), dimension(3) :: edgeNormalDotTensor
      real (kind=RKIND), dimension(3,3) :: strainRateR3Edge3x3
      real (kind=RKIND), dimension(:), pointer :: dvEdge, areaCell
      real (kind=RKIND), dimension(:,:), pointer :: edgeNormalVectors

      nCells = grid % nCells
      nVertLevels = grid % nVertLevels

      edgesOnCell       => grid % edgesOnCell % array
      nEdgesOnCell      => grid % nEdgesOnCell % array
      dvEdge            => grid % dvEdge % array
      areaCell          => grid % areaCell % array
      edgeNormalVectors  => grid % edgeNormalVectors % array

      divTensorR3Cell(:,:,:) = 0.0
      do iCell = 1, nCells
         invAreaCell = 1.0 / areaCell(iCell)
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            do k = 1, nVertLevels
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

   end subroutine mpas_divergence_of_tensor_R3Cell!}}}

!***********************************************************************
!
!  routine mpas_tensor_edge_R3_to_2D
!
!> \brief   Convert an R3 tensor to a 2D tensor, at an edge
!> \author  Mark Petersen
!> \date    April 2013
!> \details 
!>  Given an R3 tensor in symetric 6-index form, this routine rotates
!>  the tensor so that the 1-direction is towards the edge normal, and 
!>  the 2-direction is towards the edge tangent, and returns a 2D
!>  tensor in symmetric 3-index form. 
!
!-----------------------------------------------------------------------

   subroutine mpas_tensor_edge_R3_to_2D(strainRateR3Edge, strainRate2DEdge, grid, edgeTangentVectors)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         strainRateR3Edge  !< Input: strain rate tensor at edge, R3, in symmetric 6-index form

      type (mesh_type), intent(in) :: &
         grid          !< Input: grid information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         edgeTangentVectors   !< Input: unit vector tangent to an edge

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

      real (kind=RKIND), dimension(3) :: edgeVerticalVector
      real (kind=RKIND), dimension(3,3) :: rotationMatrix, strainRateR3Edge3x3, strainRateR3Edge3x3Rotated
      real (kind=RKIND), dimension(:,:), pointer :: edgeNormalVectors

      nEdges = grid % nEdges
      nVertLevels = grid % nVertLevels

      edgeNormalVectors  => grid % edgeNormalVectors % array

      do iEdge=1,nEdges

         ! compute vertical vector at edge
         call mpas_cross_product_in_r3(edgeNormalVectors(:,iEdge),edgeTangentVectors(:,iEdge),edgeVerticalvector)

         rotationMatrix(:,1) = edgeNormalVectors(:,iEdge)
         rotationMatrix(:,2) = edgeTangentVectors(:,iEdge)
         rotationMatrix(:,3) = edgeVerticalvector

         do k=1,nVertLevels
 
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
!> \date    April 2013
!> \details 
!>  Given a 2D tensor in symetric 3-index form that is rotated such that
!>  the 1-direction is towards the edge normal, and 
!>  the 2-direction is towards the edge tangent, this routine rotates
!>  the tensor to R3, and returns an R3 tensor in symetric 6-index form.
!
!-----------------------------------------------------------------------

   subroutine mpas_tensor_edge_2D_to_R3(strainRate2DEdge, strainRateR3Edge, grid, edgeTangentVectors)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         strainRate2DEdge   !< Input: strain rate tensor at edge, 2D, in symmetric 3-index form

      type (mesh_type), intent(in) :: &
         grid          !< Input: grid information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         edgeTangentVectors   !< Input: unit vector tangent to an edge

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

      real (kind=RKIND), dimension(3) :: edgeVerticalVector
      real (kind=RKIND), dimension(3,3) :: rotationMatrix, strainRateR3Edge3x3, strainRateR3Edge3x3Rotated
      real (kind=RKIND), dimension(:,:), pointer :: edgeNormalVectors

      nEdges = grid % nEdges
      nVertLevels = grid % nVertLevels

      edgeNormalVectors  => grid % edgeNormalVectors % array

      do iEdge=1,nEdges

         ! compute vertical vector at edge
         call mpas_cross_product_in_r3(edgeNormalVectors(:,iEdge),edgeTangentVectors(:,iEdge),edgeVerticalVector)

         rotationMatrix(:,1) = edgeNormalVectors(:,iEdge)
         rotationMatrix(:,2) = edgeTangentVectors(:,iEdge)
         rotationMatrix(:,3) = edgeVerticalVector

         do k=1,nVertLevels

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
!  routine mpas_tensor_LonLat_to_R3
!
!> \brief   Convert an R3 tensor to a 2D tensor
!> \author  Mark Petersen
!> \date    April 2013
!> \details 
!>  Given a 2D tensor in symetric 3-index form that is rotated such that
!>  the 1-direction is zonal
!>  the 2-direction is meridional, this routine rotates
!>  the tensor to R3, and returns an R3 tensor in symetric 6-index form.
!
!-----------------------------------------------------------------------

   subroutine mpas_tensor_LonLat_to_R3(strainRateLonLat, strainRateR3, lon, lat)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(3), intent(in) :: &
         strainRateLonLat   !< Input: strain rate tensor, 2D, in symmetric 3-index form

      real (kind=RKIND), intent(in) :: &
         lon, &!< Input: longitude, in radians, ranging [0,2*pi]
         lat   !< Input: latitude,  in radians, ranging [-pi,pi]
 
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

      integer :: i,j,k,p,q

      real (kind=RKIND) :: determinant
      real (kind=RKIND), dimension(3) :: zonalUnitVector, meridionalUnitVector, verticalUnitVector
      real (kind=RKIND), dimension(3,3) :: rotationMatrix, strainRateR3_3x3, strainRateLonLat3x3

      call mpas_zonal_meridional_vectors(lon, lat, zonalUnitVector, meridionalUnitVector, verticalUnitVector)

      rotationMatrix(:,1) = zonalUnitVector
      rotationMatrix(:,2) = meridionalUnitVector
      rotationMatrix(:,3) = verticalUnitVector

      strainRateLonLat3x3 = 0.0
      strainRateLonLat3x3(1,1) = strainRateLonLat(1)
      strainRateLonLat3x3(2,2) = strainRateLonLat(2)
      strainRateLonLat3x3(1,2) = strainRateLonLat(3)
      strainRateLonLat3x3(2,1) = strainRateLonLat(3)

      strainRateR3_3x3 = 0.0
      do i=1,3
        do j=1,3
          do p=1,3
            do q=1,3
               strainRateR3_3x3(i,j) = strainRateR3_3x3(i,j) + rotationMatrix(i,p)*strainRateLonLat3x3(p,q)*rotationMatrix(j,q)
            enddo
          enddo
        enddo
      enddo
 
      call mpas_matrix_3x3_to_sym6index(strainRateR3_3x3,strainRateR3)

   end subroutine mpas_tensor_LonLat_to_R3!}}}


!***********************************************************************
!
!  routine mpas_tensor_LonLatR_to_R3
!
!> \brief   Convert an R3 tensor to a 2D tensor
!> \author  Mark Petersen
!> \date    April 2013
!> \details 
!>  Given a 3D tensor in symetric 6-index form that is rotated such that
!>  the 1-direction is zonal
!>  the 2-direction is meridional, 
!>  the 3-direction is radial, this routine rotates
!>  the tensor to R3, and returns an R3 tensor in symetric 6-index form.
!
!-----------------------------------------------------------------------

   subroutine mpas_tensor_LonLatR_to_R3(tensorLonLatR, strainRateR3, lon, lat)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(6), intent(in) :: &
         tensorLonLatR   !< Input: latlon strain rate tensor, 3D, in symmetric 6-index form

      real (kind=RKIND), intent(in) :: &
         lon, &!< Input: longitude, in radians, ranging [0,2*pi]
         lat   !< Input: latitude,  in radians, ranging [-pi,pi]
 
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

      integer :: i,j,k,p,q

      real (kind=RKIND) :: determinant
      real (kind=RKIND), dimension(3) :: zonalUnitVector, meridionalUnitVector, verticalUnitVector
      real (kind=RKIND), dimension(3,3) :: rotationMatrix, strainRateR3_3x3, tensorLonLatR3x3, A

      call mpas_zonal_meridional_vectors(lon, lat, zonalUnitVector, meridionalUnitVector, verticalUnitVector)

      rotationMatrix(:,1) = zonalUnitVector
      rotationMatrix(:,2) = meridionalUnitVector
      rotationMatrix(:,3) = verticalUnitVector

      call mpas_matrix_sym6index_to_3x3(tensorLonLatR,tensorLonLatR3x3)

      strainRateR3_3x3 = 0.0
      do i=1,3
        do j=1,3
          do p=1,3
            do q=1,3
               strainRateR3_3x3(i,j) = strainRateR3_3x3(i,j) + rotationMatrix(i,p)*tensorLonLatR3x3(p,q)*rotationMatrix(j,q)
            enddo
          enddo
        enddo
      enddo
 
      call mpas_matrix_3x3_to_sym6index(strainRateR3_3x3,strainRateR3)

   end subroutine mpas_tensor_LonLatR_to_R3!}}}

!***********************************************************************
!
!  routine mpas_tensor_R3_to_LonLat
!
!> \brief   Convert an R3 tensor to a 2D tensor
!> \author  Mark Petersen
!> \date    April 2013
!> \details 
!>  Given a 2D tensor in symetric 3-index form that is rotated such that
!>  the 1-direction is zonal
!>  the 2-direction is meridional, this routine rotates
!>  the tensor to R3, and returns an R3 tensor in symetric 6-index form.
!
!-----------------------------------------------------------------------

   subroutine mpas_tensor_R3_to_LonLat(strainRateR3, strainRateLonLat, lon, lat)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(6), intent(in) :: &
         strainRateR3  !< Input: strain rate tensor at, R3, in symmetric 6-index form

      real (kind=RKIND), intent(in) :: &
         lon, &!< Input: longitude, in radians, ranging [0,2*pi]
         lat   !< Input: latitude,  in radians, ranging [-pi,pi]
 
      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(3), intent(out) :: &
         strainRateLonLat   !< Output: strain rate tensor, 2D, in symmetric 3-index form

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: i,j,k,p,q

      real (kind=RKIND), dimension(3) :: zonalUnitVector, meridionalUnitVector, verticalUnitVector
      real (kind=RKIND), dimension(3,3) :: rotationMatrix, strainRateR3_3x3, strainRateLonLat3x3

      call mpas_zonal_meridional_vectors(lon, lat, zonalUnitVector, meridionalUnitVector, verticalUnitVector)

      rotationMatrix(:,1) = zonalUnitVector
      rotationMatrix(:,2) = meridionalUnitVector
      rotationMatrix(:,3) = verticalUnitVector

      call mpas_matrix_sym6index_to_3x3(strainRateR3,strainRateR3_3x3)

      strainRateLonLat3x3 = 0
      do i=1,3
        do j=1,3
          do p=1,3
            do q=1,3
               strainRateLonLat3x3(i,j) = strainRateLonLat3x3(i,j) + rotationMatrix(p,i)*strainRateR3_3x3(p,q)*rotationMatrix(q,j)
            enddo
          enddo
        enddo
      enddo
 
      strainRateLonLat(1) = strainRateLonLat3x3(1,1)
      strainRateLonLat(2) = strainRateLonLat3x3(2,2)
      strainRateLonLat(3) = 0.5*(strainRateLonLat3x3(1,2)+strainRateLonLat3x3(2,1))

   end subroutine mpas_tensor_R3_to_LonLat!}}}


!***********************************************************************
!
!  routine mpas_tensor_R3_to_LonLatR
!
!> \brief   Convert an R3 tensor to a 2D tensor
!> \author  Mark Petersen
!> \date    April 2013
!> \details 
!>  Given a 2D tensor in symetric 3-index form that is rotated such that
!>  the 1-direction is zonal
!>  the 2-direction is meridional, this routine rotates
!>  the tensor to R3, and returns an R3 tensor in symetric 6-index form.
!
!-----------------------------------------------------------------------

   subroutine mpas_tensor_R3_to_LonLatR(strainRateR3, tensorLonLatR, lon, lat)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(6), intent(in) :: &
         strainRateR3  !< Input: strain rate tensor at, R3, in symmetric 6-index form

      real (kind=RKIND), intent(in) :: &
         lon, &!< Input: longitude, in radians, ranging [0,2*pi]
         lat   !< Input: latitude,  in radians, ranging [-pi,pi]
 
      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(6), intent(out) :: &
         tensorLonLatR   !< Output: strain rate tensor, 3D lat-lon coord, in symmetric 6-index form

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: i,j,k,p,q

      real (kind=RKIND), dimension(3) :: zonalUnitVector, meridionalUnitVector, verticalUnitVector
      real (kind=RKIND), dimension(3,3) :: rotationMatrix, strainRateR3_3x3, strainRateLonLat3x3

      call mpas_zonal_meridional_vectors(lon, lat, zonalUnitVector, meridionalUnitVector, verticalUnitVector)

      rotationMatrix(:,1) = zonalUnitVector
      rotationMatrix(:,2) = meridionalUnitVector
      rotationMatrix(:,3) = verticalUnitVector

      call mpas_matrix_sym6index_to_3x3(strainRateR3,strainRateR3_3x3)

      strainRateLonLat3x3 = 0
      do i=1,3
        do j=1,3
          do p=1,3
            do q=1,3
               strainRateLonLat3x3(i,j) = strainRateLonLat3x3(i,j) + rotationMatrix(p,i)*strainRateR3_3x3(p,q)*rotationMatrix(q,j)
            enddo
          enddo
        enddo
      enddo

      call mpas_matrix_3x3_to_sym6index(strainRateLonLat3x3,tensorLonLatR)

   end subroutine mpas_tensor_R3_to_LonLatR!}}}

!***********************************************************************
!
!  routine mpas_test_tensor
!
!> \brief   Tests strain rate and tensor divergence operators
!> \author  Mark Petersen
!> \date    April 2013
!> \details 
!>  This routine tests strain rate and tensor divergence operators.
!
!-----------------------------------------------------------------------

   subroutine mpas_test_tensor(domain, tensor_test_function, &
         edgeSignOnCell_field, edgeTangentVectors_field, normalVelocity_field, tangentialVelocity_field, &
         strainRateR3Cell_field, &
         strainRateR3CellSolution_field, &
         strainRateR3Edge_field, &
         strainRateLonLatRCell_field, &
         strainRateLonLatRCellSolution_field, &
         strainRateLonLatREdge_field, &
         divTensorR3Cell_field, &
         divTensorR3CellSolution_field, &
         divTensorLonLatRCell_field, &
         divTensorLonLatRCellSolution_field, &
         outerProductEdge_field  ) !{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (field2dInteger), pointer :: &
         edgeSignOnCell_field

      type (field2dReal), pointer :: &
         normalVelocity_field, tangentialVelocity_field, edgeTangentVectors_field

      type (field3dReal), pointer :: &
         strainRateLonLatRCell_field, strainRateLonLatRCellSolution_field, &
         divTensorLonLatRCell_field,  divTensorLonLatRCellSolution_field, &
         strainRateR3Cell_field, strainRateR3CellSolution_field, strainRateR3Edge_field, &
         divTensorR3Cell_field, divTensorR3CellSolution_field, &
         outerProductEdge_field,  strainRateLonLatREdge_field

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      character(len=*) :: tensor_test_function
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
      type (dm_info), pointer :: dminfo

      integer :: nCells, nCellsSolve, nCellsSolveSum, nCellsGlobal, nEdges, nVertices, nVertLevels, iCell, iEdge, k, i, p, strainRateLonLatRIndex
      integer, dimension(:,:), pointer :: edgeSignOnCell

      real (kind=RKIND) :: xVelocity, yVelocity, cn, cs, r, theta, rot, f, g1, g2, fcos, pi2l, ld, &
        lon, lat, velocityZonal, velocityMeridional, rotAngle
      real (kind=RKIND), dimension(6) ::  strainRateLonLatRDiffSum, strainRateR3DiffSum, strainRateLonLatRDiffSumGlobal, strainRateR3DiffSumGlobal
      real (kind=RKIND), dimension(3) ::  divTensorLonLatRDiffSum, divTensorR3DiffSum, divTensorLonLatRDiffSumGlobal, divTensorR3DiffSumGlobal

      real (kind=RKIND), dimension(:), pointer :: angleEdge, xCell, yCell, xEdge, yEdge, latCell, lonCell, latEdge, lonEdge
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity, tangentialVelocity, edgeNormalVectors, edgeTangentVectors
      real (kind=RKIND), dimension(:,:,:), pointer :: &
         strainRateLonLatRCell, strainRateLonLatRCellSolution, &
         divTensorLonLatRCell,  divTensorLonLatRCellSolution, &
         strainRateR3Cell, strainRateR3CellSolution, strainRateR3Edge, divTensorR3Cell, divTensorR3CellSolution, &
         outerProductEdge,  strainRateLonLatREdge
      type (field2dInteger), pointer :: &
         edgeSignOnCell_field_ptr
      type (field2dReal), pointer :: &
         normalVelocity_field_ptr, tangentialVelocity_field_ptr, edgeTangentVectors_field_ptr
      type (field3dReal), pointer :: &
         strainRateLonLatRCell_field_ptr, strainRateLonLatRCellSolution_field_ptr, &
         divTensorLonLatRCell_field_ptr,  divTensorLonLatRCellSolution_field_ptr, &
         strainRateR3Cell_field_ptr, strainRateR3CellSolution_field_ptr, strainRateR3Edge_field_ptr, &
         divTensorR3Cell_field_ptr, divTensorR3CellSolution_field_ptr, &
         outerProductEdge_field_ptr,  strainRateLonLatREdge_field_ptr
      logical :: computeStrainRate

    nCellsSolveSum = 0
    strainRateLonLatRDiffSum = 0.0
    divTensorLonLatRDiffSum = 0.0
    strainRateR3DiffSum = 0.0
    divTensorR3DiffSum = 0.0

    edgeSignOnCell_field_ptr                => edgeSignOnCell_field
    edgeTangentVectors_field_ptr            => edgeTangentVectors_field
    normalVelocity_field_ptr                => normalVelocity_field
    tangentialVelocity_field_ptr            => tangentialVelocity_field
    strainRateR3Cell_field_ptr              => strainRateR3Cell_field
    strainRateR3CellSolution_field_ptr      => strainRateR3CellSolution_field
    strainRateR3Edge_field_ptr              => strainRateR3Edge_field
    strainRateLonLatRCell_field_ptr         => strainRateLonLatRCell_field
    strainRateLonLatRCellSolution_field_ptr => strainRateLonLatRCellSolution_field
    strainRateLonLatREdge_field_ptr         => strainRateLonLatREdge_field
    divTensorR3Cell_field_ptr               => divTensorR3Cell_field
    divTensorR3CellSolution_field_ptr       => divTensorR3CellSolution_field
    divTensorLonLatRCell_field_ptr          => divTensorLonLatRCell_field
    divTensorLonLatRCellSolution_field_ptr  => divTensorLonLatRCellSolution_field
    outerProductEdge_field_ptr              => outerProductEdge_field

    block => domain % blocklist
    dminfo => domain % dminfo
    do while (associated(block))

      nCells      = block % mesh % nCells
      nCellsSolve = block % mesh % nCellsSolve
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

      edgeSignOnCell                => edgeSignOnCell_field_ptr % array
      edgeTangentVectors            => edgeTangentVectors_field_ptr % array
      normalVelocity                => normalVelocity_field_ptr % array
      tangentialVelocity            => tangentialVelocity_field_ptr % array
      strainRateR3Cell              => strainRateR3Cell_field_ptr % array
      strainRateR3CellSolution      => strainRateR3CellSolution_field_ptr % array
      strainRateR3Edge              => strainRateR3Edge_field_ptr % array
      strainRateLonLatRCell         => strainRateLonLatRCell_field_ptr % array
      strainRateLonLatRCellSolution => strainRateLonLatRCellSolution_field_ptr % array
      strainRateLonLatREdge         => strainRateLonLatREdge_field_ptr % array
      divTensorR3Cell               => divTensorR3Cell_field_ptr % array
      divTensorR3CellSolution       => divTensorR3CellSolution_field_ptr % array
      divTensorLonLatRCell          => divTensorLonLatRCell_field_ptr % array
      divTensorLonLatRCellSolution  => divTensorLonLatRCellSolution_field_ptr % array
      outerProductEdge              => outerProductEdge_field_ptr % array

      strainRateR3Cell              = 0.0
      strainRateR3CellSolution      = 0.0
      strainRateR3Edge              = 0.0
      strainRateLonLatRCell         = 0.0
      strainRateLonLatRCellSolution = 0.0
      strainRateLonLatREdge         = 0.0
      divTensorR3Cell               = 0.0
      divTensorR3CellSolution       = 0.0
      divTensorLonLatRCell          = 0.0
      divTensorLonLatRCellSolution  = 0.0
      outerProductEdge              = 0.0

      ! create test functions for normalVelocity and tangentialVelocity
      normalVelocity = 0.0
      tangentialVelocity = 0.0
      strainRateR3CellSolution = 0.0
      divTensorR3CellSolution = 0.0

      cn = 15.0e4
      cs = 20.0e4
      rot = 1.0 !pii/4.0 !0.0*pii/2
      p = 2 ! power for polynomial function

      ld = 100.0e3  ! wavelength in meters
      pi2l = pii/ld
      g1 = cn*cos(rot) - cs*sin(rot)
      g2 = cn*sin(rot) + cs*cos(rot)

      R = block % mesh % sphere_radius
      computeStrainRate = .true.

     if (tensor_test_function.eq.'constant') then
       ! linear function in x
        do iEdge = 1,nEdges
           xVelocity = cn
           yVelocity = cs
           normalVelocity(:,iEdge) = xVelocity*edgeNormalVectors(1,iEdge) + yVelocity*edgeNormalVectors(2,iEdge)
           tangentialVelocity(:,iEdge) = xVelocity*edgeTangentVectors(1,iEdge) + yVelocity*edgeTangentVectors(2,iEdge)
        enddo

        strainRateR3CellSolution = 0.0
        divTensorR3CellSolution = 0.0

     elseif (tensor_test_function.eq.'linear_x') then
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

     elseif (tensor_test_function.eq.'linear_y') then

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

     elseif (tensor_test_function.eq.'linear_arb_rot') then

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

     elseif (tensor_test_function.eq.'power_x') then

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

     elseif (tensor_test_function.eq.'power_y') then

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

     elseif (tensor_test_function.eq.'power_arb_rot') then

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

     elseif (tensor_test_function.eq.'sin_arb_rot') then

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

     elseif (tensor_test_function.eq.'spherical_solid') then

        do iEdge = 1,nEdges
           lon = lonEdge(iEdge)
           lat = latEdge(iEdge)

           velocityZonal      = cos(lat)
           velocityMeridional = 0.0
           normalVelocity(:,iEdge) = velocityZonal*cos(angleEdge(iEdge)) + velocityMeridional*sin(angleEdge(iEdge))
           tangentialVelocity(:,iEdge) = -velocityZonal*sin(angleEdge(iEdge)) + velocityMeridional*cos(angleEdge(iEdge))
        enddo

        do iCell = 1,nCells

           lon = lonCell(iCell)
           lat = latCell(iCell)

           strainRateLonLatRCellSolution(1,1,iCell) = & ! Elonlon
                 0.0
           strainRateLonLatRCellSolution(2,1,iCell) = & ! Elatlat
                 0.0
           strainRateLonLatRCellSolution(3,1,iCell) = & ! Err
                 0.0
           strainRateLonLatRCellSolution(4,1,iCell) = & ! Elonlat
                 0.0
           strainRateLonLatRCellSolution(5,1,iCell) = & ! Elatr
                 0.0
           strainRateLonLatRCellSolution(6,1,iCell) = & ! Elonr
                 - 3.0/R/2.0*cos(lat) 

           call mpas_tensor_LonLatR_to_R3(strainRateLonLatRCellSolution(:,1,iCell), strainRateR3CellSolution(:,1,iCell), lon, lat)
        enddo

        divTensorR3CellSolution = 0.0

     elseif (tensor_test_function.eq.'sph_Williamson') then

        rotAngle=pii/4

        do iEdge = 1,nEdges
           lon = lonEdge(iEdge)
           lat = latEdge(iEdge)

           velocityZonal      =  cos(lat)*cos(rotAngle) + sin(lat)*cos(lon)*sin(rotAngle)
           velocityMeridional = -sin(lat)*sin(rotAngle)

           normalVelocity(:,iEdge) = velocityZonal*cos(angleEdge(iEdge)) + velocityMeridional*sin(angleEdge(iEdge))
           tangentialVelocity(:,iEdge) = -velocityZonal*sin(angleEdge(iEdge)) + velocityMeridional*cos(angleEdge(iEdge))
        enddo

        do iCell = 1,nCells

           lon = lonCell(iCell)
           lat = latCell(iCell)

           strainRateLonLatRCellSolution(1,1,iCell) = & ! Elonlon
                 1/R/cos(lat)*( sin(rotAngle)*sin(lat)**2 - sin(rotAngle)*sin(lat)*sin(lon))
           strainRateLonLatRCellSolution(2,1,iCell) = & ! Elatlat
                 -1/R*cos(lat)*sin(rotAngle)
           strainRateLonLatRCellSolution(3,1,iCell) = & ! Err
                 0.0
           strainRateLonLatRCellSolution(4,1,iCell) = & ! Elonlat
                 1/2.0/R*(  cos(lat)*cos(lon)*sin(rotAngle) - cos(rotAngle)*sin(lat) + tan(lat)*cos(rotAngle)*cos(lat) + tan(lat)*cos(lon)*sin(rotAngle)*sin(lat))
           strainRateLonLatRCellSolution(5,1,iCell) = & ! Elatr
                 3.0/R/2.0*sin(lat)*sin(rotAngle)
           strainRateLonLatRCellSolution(6,1,iCell) = & ! Elonr
                 - 3.0/R/2.0*( cos(lat)*cos(rotAngle) + cos(lon)*sin(rotAngle)*sin(lat))

           call mpas_tensor_LonLatR_to_R3(strainRateLonLatRCellSolution(:,1,iCell), strainRateR3CellSolution(:,1,iCell), lon, lat)

        end do

     elseif (tensor_test_function.eq.'sph_uCosCos') then

        do iEdge = 1,nEdges
           lon = lonEdge(iEdge)
           lat = latEdge(iEdge)

           velocityZonal      = cos(lon)*(1+cos(2*lat))
           velocityMeridional = 0.0

           normalVelocity(:,iEdge) = velocityZonal*cos(angleEdge(iEdge)) + velocityMeridional*sin(angleEdge(iEdge))
           tangentialVelocity(:,iEdge) = -velocityZonal*sin(angleEdge(iEdge)) + velocityMeridional*cos(angleEdge(iEdge))
        enddo

        do iCell = 1,nCells

           lon = lonCell(iCell)
           lat = latCell(iCell)

           strainRateLonLatRCellSolution(1,1,iCell) = & ! Elonlon
                 -1/R*(1+cos(2*lat))*sin(lon)/cos(lat)
           strainRateLonLatRCellSolution(2,1,iCell) = & ! Elatlat
                 0.0
           strainRateLonLatRCellSolution(3,1,iCell) = & ! Err
                 0.0
           strainRateLonLatRCellSolution(4,1,iCell) = & ! Elonlat
                 1/2.0/R*( -2*cos(lon)*sin(2*lat) + (1+cos(2*lat))*cos(lon)*tan(lat))
           strainRateLonLatRCellSolution(5,1,iCell) = & ! Elatr
                 0.0
           strainRateLonLatRCellSolution(6,1,iCell) = & ! Elonr
                 - 3.0/R/2.0*(1+cos(2*lat))*cos(lon)

           call mpas_tensor_LonLatR_to_R3(strainRateLonLatRCellSolution(:,1,iCell), strainRateR3CellSolution(:,1,iCell), lon, lat)

        end do

     elseif (tensor_test_function.eq.'sph_vCosCos') then

        do iEdge = 1,nEdges
           lon = lonEdge(iEdge)
           lat = latEdge(iEdge)

           velocityZonal      = 0.0
           velocityMeridional = cos(lon)*(1+cos(2*lat))

           normalVelocity(:,iEdge) = velocityZonal*cos(angleEdge(iEdge)) + velocityMeridional*sin(angleEdge(iEdge))
           tangentialVelocity(:,iEdge) = -velocityZonal*sin(angleEdge(iEdge)) + velocityMeridional*cos(angleEdge(iEdge))
        enddo

        do iCell = 1,nCells

           lon = lonCell(iCell)
           lat = latCell(iCell)

           strainRateLonLatRCellSolution(1,1,iCell) = & ! Elonlon
                 -1/R*(1+cos(2*lat))*cos(lon)*tan(lat)
           strainRateLonLatRCellSolution(2,1,iCell) = & ! Elatlat
                 -2/R*cos(lon)*sin(2*lat)
           strainRateLonLatRCellSolution(3,1,iCell) = & ! Err
                 0.0
           strainRateLonLatRCellSolution(4,1,iCell) = & ! Elonlat
                 -1/2.0/R*(1+cos(2*lat))*sin(lon)/cos(lat)
           strainRateLonLatRCellSolution(5,1,iCell) = & ! Elatr
                 - 3.0/R/2.0*(1+cos(2*lat))*cos(lon)
           strainRateLonLatRCellSolution(6,1,iCell) = & ! Elonr
                 0.0

           call mpas_tensor_LonLatR_to_R3(strainRateLonLatRCellSolution(:,1,iCell), strainRateR3CellSolution(:,1,iCell), lon, lat)

        end do

     elseif (tensor_test_function.eq.'sph_ELonLon_CosCos') then
        computeStrainRate = .false.
        strainRateLonLatRIndex = 1

        strainRateLonLatREdge = 0.0
        do iEdge = 1,nEdges
           lon = lonEdge(iEdge)
           lat = latEdge(iEdge)
           strainRateLonLatREdge(strainRateLonLatRIndex,1,iEdge) = cos(lon)*(1+cos(2*lat))
           call mpas_tensor_LonLatR_to_R3(strainRateLonLatREdge(:,1,iEdge), strainRateR3Edge(:,1,iEdge), lon, lat)
        enddo

        strainRateLonLatRCell = 0.0
        do iCell = 1,nCells
           lon = lonCell(iCell)
           lat = latCell(iCell)

           divTensorLonLatRCellSolution(1,:,iCell) = -1/R*(1+cos(2*lat))*sin(lon)/cos(lat)
           divTensorLonLatRCellSolution(2,:,iCell) =  1/R*(1+cos(2*lat))*cos(lon)*tan(lat)
           divTensorLonLatRCellSolution(3,:,iCell) = -1/R*(1+cos(2*lat))*cos(lon)

        end do

     elseif (tensor_test_function.eq.'sph_ELatLat_CosCos') then
        computeStrainRate = .false.
        strainRateLonLatRIndex = 2

        strainRateLonLatREdge = 0.0
        do iEdge = 1,nEdges
           lon = lonEdge(iEdge)
           lat = latEdge(iEdge)
           strainRateLonLatREdge(strainRateLonLatRIndex,1,iEdge) = cos(lon)*(1+cos(2*lat))
           call mpas_tensor_LonLatR_to_R3(strainRateLonLatREdge(:,1,iEdge), strainRateR3Edge(:,1,iEdge), lon, lat)
        enddo

        strainRateLonLatRCell = 0.0
        do iCell = 1,nCells
           lon = lonCell(iCell)
           lat = latCell(iCell)

           divTensorLonLatRCellSolution(1,:,iCell) =  0.0
           divTensorLonLatRCellSolution(2,:,iCell) = -2/R*cos(lon)*sin(2*lat) -1/R*(1+cos(2*lat))*cos(lon)*tan(lat)
           divTensorLonLatRCellSolution(3,:,iCell) = -1/R*(1+cos(2*lat))*cos(lon)

        end do

     elseif (tensor_test_function.eq.'sph_ELonLat_CosCos') then
        computeStrainRate = .false.
        strainRateLonLatRIndex = 4

        strainRateLonLatREdge = 0.0
        do iEdge = 1,nEdges
           lon = lonEdge(iEdge)
           lat = latEdge(iEdge)
           strainRateLonLatREdge(strainRateLonLatRIndex,1,iEdge) = cos(lon)*(1+cos(2*lat))
           call mpas_tensor_LonLatR_to_R3(strainRateLonLatREdge(:,1,iEdge), strainRateR3Edge(:,1,iEdge), lon, lat)
        enddo

        strainRateLonLatRCell = 0.0
        do iCell = 1,nCells
           lon = lonCell(iCell)
           lat = latCell(iCell)

           divTensorLonLatRCellSolution(1,:,iCell) = -2/R*cos(lon)*sin(2*lat) -2/R*(1+cos(2*lat))*cos(lon)*tan(lat)
           divTensorLonLatRCellSolution(2,:,iCell) = -1/R*(1+cos(2*lat))*sin(lon)/cos(lat)
           divTensorLonLatRCellSolution(3,:,iCell) =  0.0

        end do

     else
       print *, 'bad choice of tensor_test_function: ',tensor_test_function
       stop
     endif

     if (computeStrainRate) then
        call mpas_strain_rate_R3Cell(normalVelocity, tangentialVelocity, strainRateR3Cell, outerProductEdge, &
          block % mesh, edgeSignOnCell, edgeTangentVectors)

        call mpas_matrix_cell_to_edge(block % mesh, strainRateR3Cell, strainRateR3Edge)
     endif

     call mpas_divergence_of_tensor_R3Cell(block % mesh, edgeSignOnCell, strainRateR3Edge, divTensorR3Cell)

     do iCell = 1,nCells

        call mpas_tensor_R3_to_LonLatR(strainRateR3Cell(:,1,iCell), strainRateLonLatRCell(:,1,iCell), lonCell(iCell), latCell(iCell))

        call mpas_vector_R3_to_LonLatR(divTensorR3Cell(:,1,iCell), divTensorLonLatRCell(:,1,iCell), lonCell(iCell), latCell(iCell))

        call mpas_vector_LonLatR_to_R3(divTensorLonLatRCellSolution(:,1,iCell),  divTensorR3CellSolution(:,1,iCell), lonCell(iCell), latCell(iCell))

     enddo

     nCellsSolveSum = nCellsSolveSum + nCellsSolve
     do iCell = 1,nCellsSolve
        strainRateLonLatRDiffSum(:) = strainRateLonLatRDiffSum(:) + (strainRateLonLatRCell(:,1,iCell) - strainRateLonLatRCellSolution(:,1,iCell))**2
        divTensorLonLatRDiffSum(:) = divTensorLonLatRDiffSum(:) + (divTensorLonLatRCell(:,1,iCell) - divTensorLonLatRCellSolution(:,1,iCell))**2
        strainRateR3DiffSum(:) = strainRateR3DiffSum(:) + (strainRateR3Cell(:,1,iCell) - strainRateR3CellSolution(:,1,iCell))**2
        divTensorR3DiffSum(:) = divTensorR3DiffSum(:) + (divTensorR3Cell(:,1,iCell) - divTensorR3CellSolution(:,1,iCell))**2
     enddo

     block => block % next
     edgeSignOnCell_field_ptr                => edgeSignOnCell_field_ptr % next
     normalVelocity_field_ptr                => normalVelocity_field_ptr % next
     tangentialVelocity_field_ptr            => tangentialVelocity_field_ptr % next
     edgeTangentVectors_field_ptr            => edgeTangentVectors_field_ptr % next
     strainRateR3Cell_field_ptr              => strainRateR3Cell_field_ptr % next
     strainRateR3CellSolution_field_ptr      => strainRateR3CellSolution_field_ptr % next
     strainRateR3Edge_field_ptr              => strainRateR3Edge_field_ptr % next
     strainRateLonLatRCell_field_ptr         => strainRateLonLatRCell_field_ptr % next
     strainRateLonLatRCellSolution_field_ptr => strainRateLonLatRCellSolution_field_ptr % next
     strainRateLonLatREdge_field_ptr         => strainRateLonLatREdge_field_ptr % next
     divTensorR3Cell_field_ptr               => divTensorR3Cell_field_ptr % next
     divTensorR3CellSolution_field_ptr       => divTensorR3CellSolution_field_ptr % next
     divTensorLonLatRCell_field_ptr          => divTensorLonLatRCell_field_ptr % next
     divTensorLonLatRCellSolution_field_ptr  => divTensorLonLatRCellSolution_field_ptr % next
     outerProductEdge_field_ptr              => outerProductEdge_field_ptr % next
   end do

   call mpas_dmpar_sum_int(dminfo, nCellsSolveSum, nCellsGlobal)
   call mpas_dmpar_sum_real_array(dminfo, 6, strainRateLonLatRDiffSum, strainRateLonLatRDiffSumGlobal)
   call mpas_dmpar_sum_real_array(dminfo, 3, divTensorLonLatRDiffSum, divTensorLonLatRDiffSumGlobal)
   call mpas_dmpar_sum_real_array(dminfo, 6, strainRateR3DiffSum, strainRateR3DiffSumGlobal)
   call mpas_dmpar_sum_real_array(dminfo, 3, divTensorR3DiffSum, divTensorR3DiffSumGlobal)

   if (dminfo % my_proc_id == IO_NODE) then
     if (computeStrainRate) then
       print '(a,10es14.6)', 'rms error, strainRateLonLatRCell:',sqrt(strainRateLonLatRDiffSumGlobal/nCellsGlobal)
       print '(a,10es14.6)', 'rms error, strainRateR3Cell:     ',sqrt(strainRateR3DiffSumGlobal/nCellsGlobal)
     else
       print '(a,10es14.6)', 'rms error, divTensorLonLatRCell: ',sqrt(divTensorLonLatRDiffSumGlobal/nCellsGlobal)
       print '(a,10es14.6)', 'rms error, divTensorR3Cell:      ',sqrt(divTensorR3DiffSumGlobal/nCellsGlobal)
     endif
   endif


   end subroutine mpas_test_tensor!}}}

end module mpas_tensor_operations


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker

