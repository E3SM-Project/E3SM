! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!***********************************************************************
!
!  mpas_rbf_interpolation
!
!> \brief   MPAS Radial basis function interpolation module
!> \author  Xylar Asay-Davis
!> \date    03/28/13
!> \details 
!> This module provides routines for performing interpolation with radial basis functions.
!> It performs interpolation of scalar and vector functions in 2 and 3 dimensions.
!
!-----------------------------------------------------------------------
module mpas_rbf_interpolation
   use mpas_dmpar
   use mpas_grid_types

   implicit none
   private
   save

   ! Initialize the geometry that will be useful from interpolation
  public :: mpas_rbf_interp_initialize

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! Routines for perofrming interpolation in 2D (including Jacobian and Hessian)
  !  at locations that vary using a function that is fixed.  This is useful
  !  for finding the location on the the RBF reconstruction of a function
  !  (e.g., a height field) that minimizes the distance to a point in 3D space
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  public :: mpas_rbf_interp_loc_2D_sca_const_comp_coeffs, &
    mpas_rbf_interp_loc_2D_sca_lin_comp_coeffs, &
    mpas_rbf_interp_loc_2D_sca_const_eval_with_derivs, &
    mpas_rbf_interp_loc_2D_sca_lin_eval_with_derivs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! Routines for computing scalar interpolaiton coefficients in 3D (coplanar points
  !  in 3D) with support for either constant or constant and linear basis
  !  functions in addition to RBFs.  In constrast to the two subroutines
  !  above, these coefficients are valid for computing the value of the
  !  interpolant at a fixe point in space using function values that may
  !  vary (e.g., in time) at each of the interpolation "source" points.
  !  The last 3 routines can be used to compute coefficients for imposing both Neumann
  !  and Dirichlet boundary conditions.
  ! Pseudocode for function reconstruction at the destinationPoint is as follows
  !  Dirichlet:  functionAtDestination = sum(functionAtSources*dirichletCoefficients)
  !  Neumann:    functionAtDestination = sum(functionOrDerivAtSources*neumannCoefficients)
  !    where functionOrDerivAtSource = functionAtSources where .not.(isInterface)
  !                                  = functionNormalDerivAtSources where isInterface
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  public :: mpas_rbf_interp_func_3D_sca_const_dir_comp_coeffs, &
    mpas_rbf_interp_func_3D_plane_sca_lin_dir_comp_coeffs, &
    mpas_rbf_interp_func_3D_sca_lin_dir_comp_coeffs, &
    mpas_rbf_interp_func_3D_sca_const_dir_neu_comp_coeffs, &
    mpas_rbf_interp_func_3D_plane_sca_lin_dir_neu_comp_coeffs, &
    mpas_rbf_interp_func_3D_sca_lin_dir_neu_comp_coeffs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! Routines for computing vector interpolaiton coefficients in 3D (coplanar points
  !  in 3D) with support for only constant basis functions in addition to RBFs. 
  !  (Linear basis functions can cause problems: a linear vortex flow u = y xHat - x yHat
  !  cannot be reconstructed from vectors that are only normal to the cell edges in 2D.
  !  Therefore, we don't support them). As with the scalar 3D subroutines
  !  above, these coefficients are valid for computing the value of the
  !  interpolant at a fixe point in space using function values that may
  !  vary (e.g., in time) at each of the interpolation "source" points.
  !  The user supplies to these routines a set of sourcePoints and unitVectors
  !  as well as a destinationPoint and, for the last 2 routines, flags for
  !  which points are tangent to the interface and which of the supplied unitVectors
  !  is the normal at the corresponding point.
  !
  ! The first two routines compute coefficients that can be used to interpolate
  !  a vector function to the destination point given "function dot unitVector" values
  !  at each source point.  These routines are useful, for example, for reconstructing
  !  the full vector velocity at cell centers from the normal components of the velocity
  !  at cell faces (or cell edges in 2D), or for computing the full velocity at an
  !  immersed boundary image point based on the normal velocity at several faces and
  !  the full velocity at boundary points (e.g., a no-slip boundary condition).
  !
  ! The last two routines compute coefficients that can be used to interpolate
  !  a vector function to the destination point given "function dot unitVector" values
  !  at non-tangent source point and "dFunction/dn dot unitVector" values at
  !  tangent source point.  These routines are useful, for example, for computing the 
  !  full velocity at an immersed boundary image point based on the normal velocity at
  !  several faces, the normal velocity at boundary points (e.g., u dot n = 0 for a
  !  no-penetration boundary condition on a fixed boundary), and the normal derivative
  !  of the tangential components of velocity at the boundary points (e.g., a free-slip
  !  boundary condition).
  ! Pseudocode for function reconstruction at the destinationPoint is as follows
  !  dirichlet:  functionAtDestination_i = sum_j(functionDotUnitVectorAtSources_j*coefficients_j,i)
  !    for i = x,y,z
  !  tangentNeumann:    functionAtDestination_i &
  !    = sum_(j where .not isTangent) (functionDotUnitVectorAtSources_j*coefficients_j,i) &
  !    + sum_(j where isTangent) ((dFunction/dn dot UnitVector)_j*coefficients_j,i)
  !    for i = x,y,z
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  public :: mpas_rbf_interp_func_3D_vec_const_dir_comp_coeffs, &
    mpas_rbf_interp_func_3D_plane_vec_const_dir_comp_coeffs!, &
    !mpas_rbf_interp_func_3D_vec_const_tan_neu_comp_coeffs, &
    !mpas_rbf_interp_func_3D_plane_vec_const_tan_neu_comp_coeffs

  contains


!***********************************************************************
!
!  routine mpas_rbf_interp_initialize
!
!> \brief   MPAS RBF interpolation initialization routine
!> \author  Xylar Asay-Davis
!> \date    03/28/13
!> \details 
!> This routine computes geometric fields that will be potentially useful for calling
!> the interpolation routines.
!> Input: the grid
!> Output:
!>      edgeNormalVectors - the unit vector at the center of each edge tangent to the sphere
!>      cellTangentPlane - 2 orthogonal unit vectors in the tangent plane of each cell
!>                         The first unit vector is chosen to point toward the center of the first
!>                         edge on the cell.
!>      localVerticalUnitVectors - the unit normal vector of the tangent plane at the center 
!>                         of each cell
!-----------------------------------------------------------------------
  subroutine mpas_rbf_interp_initialize(grid)!{{{

    implicit none

    type (mesh_type), intent(inout) :: grid  !< Input/Output: Grid information

    integer :: nCells, nEdges
    integer, dimension(:,:), pointer :: cellsOnEdge, edgesOnCell
    integer :: iEdge, iCell, cell1, cell2
    real(kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell, xEdge, yEdge, zEdge
    real(kind=RKIND), dimension(:,:), pointer :: localVerticalUnitVectors, edgeNormalVectors
    real(kind=RKIND), dimension(:,:,:), pointer :: cellTangentPlane
    real(kind=RKIND), dimension(3) :: xHatPlane, yHatPlane, rHat
    real(kind=RKIND) :: normalDotRHat
    logical :: on_a_sphere

    xCell       => grid % xCell % array
    yCell       => grid % yCell % array
    zCell       => grid % zCell % array
    xEdge       => grid % xEdge % array
    yEdge       => grid % yEdge % array
    zEdge       => grid % zEdge % array
    cellsOnEdge => grid % cellsOnEdge % array
    edgesOnCell => grid % edgesOnCell % array
    nCells      = grid % nCells
    nEdges      = grid % nEdges
    on_a_sphere = grid % on_a_sphere

    localVerticalUnitVectors => grid % localVerticalUnitVectors % array
    edgeNormalVectors => grid % edgeNormalVectors % array
    cellTangentPlane => grid % cellTangentPlane % array

    ! init arrays
    edgeNormalVectors = 0
    localVerticalUnitVectors = 0

    ! loop over all cells to be solved on this block
    do iCell=1,nCells
      if(on_a_sphere) then
        localVerticalUnitVectors(1,iCell) = xCell(iCell)
        localVerticalUnitVectors(2,iCell) = yCell(iCell)
        localVerticalUnitVectors(3,iCell) = zCell(iCell)
        call mpas_unit_vec_in_r3(localVerticalUnitVectors(:,iCell))
      else ! on a plane
        localVerticalUnitVectors(:,iCell) = (/ 0., 0., 1. /)
      end if
    end do

    do iEdge = 1,nEdges
      iCell = cellsOnEdge(1,iEdge) ! the normal vector points from the first cell toward the edge
      if(iCell == nCells+1) then ! this is a boundary edge
        ! the first cell bordering this edge is not real, use the second cell
        !  The normal should always point outward at boundaries, away from the valid cell center
        iCell = cellsOnEdge(2,iEdge)
      end if
      ! the normal points from the cell location to the edge location
      edgeNormalVectors(1,iEdge) = xEdge(iEdge) - xCell(iCell)
      edgeNormalVectors(2,iEdge) = yEdge(iEdge) - yCell(iCell)
      edgeNormalVectors(3,iEdge) = zEdge(iEdge) - zCell(iCell)
      call mpas_unit_vec_in_r3(edgeNormalVectors(:,iEdge))
    end do

    do iCell=1,nCells
      iEdge = edgesOnCell(1,iCell)
      ! xHat and yHat are a local basis in the plane of the horizontal cell
      ! we arbitrarily choose xHat to point toward the first edge
      rHat = localVerticalUnitVectors(:,iCell)
      normalDotRHat = sum(edgeNormalVectors(:,iEdge)*rHat)
      xHatPlane = edgeNormalVectors(:,iEdge) - normalDotRHat*rHat
      call mpas_unit_vec_in_r3(xHatPlane)

      call mpas_cross_product_in_r3(rHat, xHatPlane, yHatPlane)
      call mpas_unit_vec_in_r3(yHatPlane) ! just to be sure...
      cellTangentPlane(:,1,iCell) = xHatPlane
      cellTangentPlane(:,2,iCell) = yHatPlane
    end do

  end subroutine mpas_rbf_interp_initialize!}}}

!***********************************************************************
!
!  routine mpas_rbf_interp_loc_2D_sca_const_comp_coeffs
!
!> \brief   MPAS 2D scalar constant interpolation coefficient routine
!> \author  Xylar Asay-Davis
!> \date    03/28/13
!> \details 
!>  This routine computes interpolation coefficients in 2D that can be used to reconstruct a given scalar function at varying locations.
!>  This is useful for finding the location on the RBF reconstruction of a function (e.g. a heigh field) that minimizes the distantce
!>  to a point in 3D space. The reconstruction is performed with basis functions that are RBFs and constant.
!>  Input:
!>   pointCount - the number of "source" points and functionValues supplied
!>   coeffCount - the size of coefficients, must be at least pointCount + 1
!>   points - the location of the "source" points in the 2D space where the values of
!>     the function are known
!>   fieldValues - the values of the function of interest at the points
!>   alpha - a constant that give the characteristic length scale of the RBFs,
!>     should be on the order of the distance between points
!>  Output:
!>   coefficients - the coefficients needed to perform interpolation of the funciton
!>     at destination points yet to be specified       
!-----------------------------------------------------------------------
  subroutine mpas_rbf_interp_loc_2D_sca_const_comp_coeffs(pointCount, coeffCount, &!{{{
    points, fieldValues, alpha, coefficients)
 
    integer, intent(in) :: pointCount !< Input: Number of points
    integer, intent(in) :: coeffCount !< Input: Number of coefficients
    real(kind=RKIND), dimension(pointCount,2), intent(in) :: points !< Input: List of points
    real(kind=RKIND), dimension(pointCount), intent(in) :: fieldValues !< Input: Value at points
    real(kind=RKIND), intent(in) :: alpha !< Input: Charachteristic length scale of RBFs
    real(kind=RKIND), dimension(coeffCount), intent(out) :: coefficients !< Output: List of coefficients

    integer :: i, j, matrixSize
    real(kind=RKIND), dimension(pointCount+1,pointCount+1) :: matrix
    real(kind=RKIND), dimension(pointCount+1) :: rhs
    integer, dimension(pointCount+1) :: pivotIndices
    real(kind=RKIND) :: rSquared

    matrixSize = pointCount+1
    coefficients = 0.0
    matrix = 0.0
    rhs = 0.0

    rhs(1:pointCount) = fieldValues

    do j=1,pointCount
      do i=j,pointCount
        rSquared = sum((points(j,:) - points(i,:))**2)/alpha**2
        matrix(i,j) = evaluate_rbf(rSquared)
        matrix(j,i) = matrix(i,j)
      end do
    end do
    do j=1,pointCount
      matrix(pointCount+1,j) = 1.0
      matrix(j,pointCount+1) = 1.0
    end do

    call mpas_legs(matrix(1:matrixSize,1:matrixSize), matrixSize, rhs(1:matrixSize), &
      coefficients(1:matrixSize), pivotIndices(1:matrixSize))

  end subroutine mpas_rbf_interp_loc_2D_sca_const_comp_coeffs!}}}

!***********************************************************************
!
!  routine mpas_rbf_interp_loc_2D_sca_lin_comp_coeffs
!
!> \brief   MPAS 2D scalar linear interpolation coefficient routine
!> \author  Xylar Asay-Davis
!> \date    03/28/13
!> \details 
!> This routine computes interpolation coefficients in 2D that can be used to
!>  reconstruct a given scalar function at varying locations. This is useful
!>  for finding the location on the the RBF reconstruction of a function
!>  (e.g., a height field) that minimizes the distance to a point in 3D space.
!>  The reconstruction is performed with basis functions that are RBFs plus constant
!>  and linear 
!> Input:
!>  pointCount - the number of "source" points and functionValues supplied
!>  coeffCount - the size of coefficients, must be at least pointCount + 3
!>  points - the location of the "source" points in the 2D space where the values of
!>    the function are known
!>  fieldValues - the values of the function of interest at the points
!>  alpha - a constant that give the characteristic length scale of the RBFs,
!>    should be on the order of the distance between points
!> Output:
!>  coefficients - the coefficients needed to perform interpolation of the funciton
!>    at destination points yet to be specified
!-----------------------------------------------------------------------
  subroutine mpas_rbf_interp_loc_2D_sca_lin_comp_coeffs(pointCount, coeffCount, &!{{{
    points, fieldValues, alpha, coefficients)
 
    integer, intent(in) :: pointCount !< Input: Number of points
    integer, intent(in) :: coeffCount !< Input: Number of coefficients
    real(kind=RKIND), dimension(pointCount,2), intent(in) :: points !< Input: List of points
    real(kind=RKIND), dimension(pointCount), intent(in) :: fieldValues !< Input: List of values at points
    real(kind=RKIND), intent(in) :: alpha !< Input: Characteristic length scale for RBFs 
    real(kind=RKIND), dimension(coeffCount), intent(out) :: coefficients !< Output: List of coefficients

    integer :: i, j, matrixSize
    real(kind=RKIND), dimension(pointCount+3,pointCount+3) :: matrix
    real(kind=RKIND), dimension(pointCount+3) :: rhs
    integer, dimension(pointCount+3) :: pivotIndices
    real(kind=RKIND) :: rSquared

    coefficients = 0.0
    matrix = 0.0
    rhs = 0.0

    rhs(1:pointCount) = fieldValues

    do j=1,pointCount
      do i=j,pointCount
        rSquared = sum((points(j,:) - points(i,:))**2)/alpha**2
        matrix(i,j) = evaluate_rbf(rSquared)
        matrix(j,i) = matrix(i,j)
      end do
    end do
    matrixSize = pointCount+3
    do j=1,pointCount
      matrix(pointCount+1,j) = 1.0
      matrix(pointCount+2,j) = points(j,1)
      matrix(pointCount+3,j) = points(j,2)
      matrix(j,pointCount+1:pointCount+3) = matrix(pointCount+1:pointCount+3, j)
    end do
    call mpas_legs(matrix(1:matrixSize,1:matrixSize), matrixSize, rhs(1:matrixSize), &
      coefficients(1:matrixSize), pivotIndices(1:matrixSize))

  end subroutine mpas_rbf_interp_loc_2D_sca_lin_comp_coeffs!}}}

!***********************************************************************
!
!  routine mpas_rbf_interp_loc_2D_sca_const_eval_with_derivs
!
!> \brief   MPAS 2D scalar constant evaulation routine
!> \author  Xylar Asay-Davis
!> \date    03/28/13
!> \details 
!>  This routine evalutes a scalar function in 2D using coefficients computed in
!>   rbfInterp_loc_2D_sca_const_compCoeffs.  This 
!>   function can be called repeatedly with different destination points
!>   to quickly evaluate the interpolating function using the same
!>   coefficients.  This is useful for finding the location on the the 
!>   RBF reconstruction of a function (e.g., a height field) that minimizes
!>   the distance to a point in 3D space. The reconstruction is performed
!>   with basis functions that are RBFs and constant 
!>  Input:
!>   fieldCount - the number fields to be evaluated.  This is useful for reconstructing,
!>     for example, the x-, y- and z-components of a vector field at the same
!>     point in 2D
!>   coeffCount - the size of coefficients, must be at least pointCount + 1
!>   pointCount - the number of "source" points and functionValues supplied
!>   coefficients - the coefficients needed to perform interpolation of the funciton
!>     at the evaluationPoint
!>   evaluationPoint - the point in 2D where the function is to be reconstructed
!>   points - the location of the "source" points in the 2D space where the values of
!>     the function are known
!>   alpha - a constant that give the characteristic length scale of the RBFs,
!>     should be on the order of the distance between points
!>  Output:
!>   derivs - the value of the function, the 2 components of its Jacobian and
!>     the 3 unique components of its Hessian at the evaluationPoint
!-----------------------------------------------------------------------
  subroutine mpas_rbf_interp_loc_2D_sca_const_eval_with_derivs(fieldCount, coeffCount, &!{{{
    pointCount, coefficients, evaluationPoint, points, alpha, derivs)
    integer, intent(in) :: fieldCount !< Input: Number of fields
    integer, intent(in) :: coeffCount !< Input: Number of coefficients
    integer, intent(in) :: pointCount !< Input: Number of points
    real(kind=RKIND), dimension(coeffCount, fieldCount), intent(in) :: coefficients !< Input: List of coefficients
    real(kind=RKIND), dimension(2), intent(in) :: evaluationPoint !< Input: Location for evaluation
    real(kind=RKIND), dimension(pointCount,2), intent(in) :: points !< Input: List of points
    real(kind=RKIND), intent(in) :: alpha !< Input: Characteristic length scale for RBFs

    real(kind=RKIND), dimension(6,fieldCount), intent(out) :: derivs !< Output: List of derivatives

    integer :: pointIndex
    real(kind=RKIND) :: x, y, rSquared, rbfValue, rbfDerivOverR, rbfSecondDeriv

    derivs = 0.0
    do pointIndex = 1, pointCount
      x = (evaluationPoint(1) - points(pointIndex,1))
      y = (evaluationPoint(2) - points(pointIndex,2))
      rSquared = x**2 + y**2
      call mpas_evaluate_rbf_and_derivs(rSquared/alpha**2,  rbfValue, rbfDerivOverR, rbfSecondDeriv)
      rbfDerivOverR = rbfDerivOverR/alpha**2
      rbfSecondDeriv = rbfSecondDeriv/alpha**2
      if(rSquared/alpha**2 < 1e-7) then
        ! terms 2,3 and 5 are zero by radial symmetry of the RBF
        derivs(1,:) = derivs(1,:) + coefficients(pointIndex,:)*rbfValue
        derivs(4,:) = derivs(4,:) + coefficients(pointIndex,:)*rbfSecondDeriv
        derivs(6,:) = derivs(6,:) + coefficients(pointIndex,:)*rbfSecondDeriv
      else
        call mpas_evaluate_rbf_and_derivs(rSquared,  rbfValue, rbfDerivOverR, rbfSecondDeriv)
        derivs(1,:) = derivs(1,:) + coefficients(pointIndex,:)*rbfValue
        derivs(2,:) = derivs(2,:) + coefficients(pointIndex,:)*rbfDerivOverR*x
        derivs(3,:) = derivs(3,:) + coefficients(pointIndex,:)*rbfDerivOverR*y
        derivs(4,:) = derivs(4,:) + coefficients(pointIndex,:) &
          * (rbfSecondDeriv*x**2 + rbfDerivOverR*y**2)/rSquared
        derivs(5,:) = derivs(5,:) + coefficients(pointIndex,:) &
          * (rbfSecondDeriv - rbfDerivOverR)*x*y/rSquared
        derivs(6,:) = derivs(6,:) + coefficients(pointIndex,:) &
          * (rbfSecondDeriv*y**2 + rbfDerivOverR*x**2)/rSquared
      end if
    end do
    derivs(1,:) = derivs(1,:) + coefficients(pointCount+1,:)
  end subroutine mpas_rbf_interp_loc_2D_sca_const_eval_with_derivs!}}}

!***********************************************************************
!
!  routine mpas_rbf_interp_loc_2D_sca_lin_eval_with_derivs
!
!> \brief   MPAS 2D scalar linear evaluation routine
!> \author  Xylar Asay-Davis
!> \date    03/28/13
!> \details 
!>  This routine evalutes a scalar function in 2D using coefficients computed in
!>   rbfInterp_loc_2D_sca_const_compCoeffs.  This 
!>   function can be called repeatedly with different destination points
!>   to quickly evaluate the interpolating function using the same
!>   coefficients.  This is useful for finding the location on the the 
!>   RBF reconstruction of a function (e.g., a height field) that minimizes
!>   the distance to a point in 3D space. The reconstruction is performed
!>   with basis functions that are RBFs, constant and linear
!>  Input:
!>   fieldCount - the number fields to be evaluated.  This is useful for reconstructing,
!>     for example, the x-, y- and z-components of a vector field at the same
!>     point in 2D
!>   coeffCount - the size of coefficients, must be at least pointCount + 1
!>   pointCount - the number of "source" points and functionValues supplied
!>   coefficients - the coefficients needed to perform interpolation of the funciton
!>     at the evaluationPoint
!>   evaluationPoint - the point in 2D where the function is to be reconstructed
!>   points - the location of the "source" points in the 2D space where the values of
!>     the function are known
!>   alpha - a constant that give the characteristic length scale of the RBFs,
!>     should be on the order of the distance between points
!>  Output:
!>   derivs - the value of the function, the 2 components of its Jacobian and
!>     the 3 unique components of its Hessian at the evaluationPoint
!-----------------------------------------------------------------------
  subroutine mpas_rbf_interp_loc_2D_sca_lin_eval_with_derivs(fieldCount, coeffCount, &!{{{
    pointCount, coefficients, evaluationPoint, points, alpha, derivs)
    integer, intent(in) :: fieldCount !< Input: Number of fields
    integer, intent(in) :: coeffCount !< Input: Number of coefficients
    integer, intent(in) :: pointCount !< Input: Number of points
    real(kind=RKIND), dimension(coeffCount, fieldCount), intent(in) :: coefficients !< Input: List of coefficients
    real(kind=RKIND), dimension(2), intent(in) :: evaluationPoint !< Input: Point for evaluation
    real(kind=RKIND), dimension(pointCount,2), intent(in) :: points !< Input: List of points
    real(kind=RKIND), intent(in) :: alpha !< Input: Characteristic length scale of RBFs

    real(kind=RKIND), dimension(6,fieldCount), intent(out) :: derivs !< Output: Derivatives

    integer :: pointIndex
    real(kind=RKIND) :: x, y, rSquared, rbfValue, rbfDerivOverR, rbfSecondDeriv

    derivs = 0.0
    do pointIndex = 1, pointCount
      x = (evaluationPoint(1) - points(pointIndex,1))
      y = (evaluationPoint(2) - points(pointIndex,2))
      rSquared = x**2 + y**2
      call mpas_evaluate_rbf_and_derivs(rSquared/alpha**2,  rbfValue, rbfDerivOverR, rbfSecondDeriv)
      rbfDerivOverR = rbfDerivOverR/alpha**2
      rbfSecondDeriv = rbfSecondDeriv/alpha**2
      if(rSquared/alpha**2 < 1e-7) then
        ! terms 2,3 and 5 are zero by radial symmetry of the RBF
        derivs(1,:) = derivs(1,:) + coefficients(pointIndex,:)*rbfValue
        derivs(4,:) = derivs(4,:) + coefficients(pointIndex,:)*rbfSecondDeriv
        derivs(6,:) = derivs(6,:) + coefficients(pointIndex,:)*rbfSecondDeriv
      else
        call mpas_evaluate_rbf_and_derivs(rSquared,  rbfValue, rbfDerivOverR, rbfSecondDeriv)
        derivs(1,:) = derivs(1,:) + coefficients(pointIndex,:)*rbfValue
        derivs(2,:) = derivs(2,:) + coefficients(pointIndex,:)*rbfDerivOverR*x
        derivs(3,:) = derivs(3,:) + coefficients(pointIndex,:)*rbfDerivOverR*y
        derivs(4,:) = derivs(4,:) + coefficients(pointIndex,:) &
          * (rbfSecondDeriv*x**2 + rbfDerivOverR*y**2)/rSquared
        derivs(5,:) = derivs(5,:) + coefficients(pointIndex,:) &
          * (rbfSecondDeriv - rbfDerivOverR)*x*y/rSquared
        derivs(6,:) = derivs(6,:) + coefficients(pointIndex,:) &
          * (rbfSecondDeriv*y**2 + rbfDerivOverR*x**2)/rSquared
      end if
    end do
    derivs(1,:) = derivs(1,:) + coefficients(pointCount+1,:) &
      + coefficients(pointCount+2,:)*evaluationPoint(1) &
      + coefficients(pointCount+3,:)*evaluationPoint(2)
    derivs(2,:) = derivs(2,:) + coefficients(pointCount+2,:)
    derivs(3,:) = derivs(3,:) + coefficients(pointCount+3,:)

  end subroutine mpas_rbf_interp_loc_2D_sca_lin_eval_with_derivs!}}}

!***********************************************************************
!
!  routine mpas_rbf_interp_func_3D_sca_const_dir_comp_coeffs
!
!> \brief   MPAS 3D scalar constant coefficients routine
!> \author  Xylar Asay-Davis
!> \date    03/28/13
!> \details 
!>  This routine computes interpolation coefficients in 3D that can be used to
!>   interpolate a number of scalar functions at a given locations. This is useful
!>   if the interpolation location does not change with time, or if several
!>   fields are to be interpolated at a given time step.  (If both the fields
!>   and the interpolation locations vary with time, there is no clear advantage in
!>   using either this method or the method for 2D interpoaltion above; for simplicity
!>   and because we foresee more uses for the method of this subroutine, we have not
!>   implemented a 3D version of the fixed field, variable interpolation location method
!>   as we have in 2D.) Coefficients are produced for handling Dirichlet boundary
!>   conditions (or no boundaries).  The interpolation is performed with basis functions
!>   that are RBFs plus a constant.
!>  Input:
!>   pointCount - the number of "source" points and functionValues supplied
!>   sourcePoints - the location of the "source" points in the 3D space where the values of
!>     the function are known
!>   destinationPoint - the point where the interpolation will be performed
!>   alpha - a constant that give the characteristic length scale of the RBFs,
!>     should be on the order of the distance between points
!>  Output:
!>   coefficients - the coefficients used to interpolate a function with Dirichlet
!>     boundary conditions to the specified destinationPoint
!-----------------------------------------------------------------------
  subroutine mpas_rbf_interp_func_3D_sca_const_dir_comp_coeffs( &!{{{
    pointCount, sourcePoints, destinationPoint, alpha, coefficients)

    integer, intent(in) :: pointCount !< Input: Number of points
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: sourcePoints !< Input: List of source points
    real(kind=RKIND), dimension(3), intent(in) :: destinationPoint !< Input: List of destination points
    real(kind=RKIND), intent(in) :: alpha !< Input: Characteristic length scale of RBFs
    real(kind=RKIND), dimension(pointCount), intent(out) :: coefficients !< Output: List of coefficients

    integer :: i, j
    integer :: matrixSize

    real(kind=RKIND), dimension(:,:), pointer :: dirichletMatrix
    real(kind=RKIND), dimension(:), pointer :: rhs, coeffs
    integer, dimension(:), pointer :: pivotIndices

    matrixSize = pointCount+1 !! 1 extra space for constant 

    allocate(dirichletMatrix(matrixSize,matrixSize))  
    allocate(rhs(matrixSize))  
    allocate(coeffs(matrixSize))  
    allocate(pivotIndices(matrixSize))  

    dirichletMatrix = 0.0
    rhs = 0.0
    coeffs = 0.0

    call mpas_set_up_scalar_rbf_dirichlet_matrix_and_rhs(pointCount, sourcePoints, destinationPoint, &
      alpha, dirichletMatrix(1:pointCount,1:pointCount), rhs(1:pointCount))

    do i = 1, pointCount
      dirichletMatrix(i,pointCount+1) = 1.0
      dirichletMatrix(pointCount+1,i) = dirichletMatrix(i,pointCount+1)
    end do

    rhs(pointCount+1) = 1.0

    ! solve each linear system
    call mpas_legs(dirichletMatrix, matrixSize, rhs, coeffs, pivotIndices)
    coefficients = coeffs(1:pointCount)

    deallocate(dirichletMatrix)  
    deallocate(rhs)  
    deallocate(coeffs)  
    deallocate(pivotIndices) 

  end subroutine mpas_rbf_interp_func_3D_sca_const_dir_comp_coeffs!}}}

!***********************************************************************
!
!  routine mpas_rbf_interp_func_3D_plane_sca_lin_dir_comp_coeffs
!
!> \brief   MPAS 3D planar scalar linear coefficients routine
!> \author  Xylar Asay-Davis
!> \date    03/28/13
!> \details 
!>  This routine computes interpolation coefficients in a plane in 3D that can be used to
!>   interpolate a number of scalar functions at a given locations. This is useful
!>   if the interpolation location does not change with time, or if several
!>   fields are to be interpolated at a given time step.  (If both the fields
!>   and the interpolation locations vary with time, there is no clear advantage in
!>   using either this method or the method for 2D interpoaltion above; for simplicity
!>   and because we foresee more uses for the method of this subroutine, we have not
!>   implemented a 3D version of the fixed field, variable interpolation location method
!>   as we have in 2D.) Coefficients are produced for handling both Dirichlet (or no)
!>   boundary conditions.  The interpolation is performed with basis functions that are
!>   RBFs plus constant and linear.  All points are projected into the plane given by the
!>   planeBasisVectors.
!>  Input:
!>   pointCount - the number of "source" points and functionValues supplied
!>   sourcePoints - the location of the "source" points in the 3D space where the values of
!>     the function are known.  The points will be projected into the plane given by 
!>     planeBasisVectors
!>   destinationPoint - the point in 3D where the interpolation will be performed.  The
!>     destinationPoint will be projected into the plane given by planeBasisVectors.
!>   alpha - a constant that give the characteristic length scale of the RBFs,
!>     should be on the order of the distance between points
!>   planeBasisVectors - the basis fectors for the plane where interpolation is performed.
!>     All points are projected into this plane. 
!>  Output:
!>   coefficients - the coefficients used to interpolate a function with Dirichlet
!>     boundary conditions to the specified destinationPoint
!-----------------------------------------------------------------------
  subroutine mpas_rbf_interp_func_3D_plane_sca_lin_dir_comp_coeffs( &!{{{
    pointCount, sourcePoints, destinationPoint, &
    alpha, planeBasisVectors, coefficients)

    integer, intent(in) :: pointCount !< Input: Number of points
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: sourcePoints !< Input: List of source points
    real(kind=RKIND), dimension(3), intent(in) :: destinationPoint !< Input: Destination point
    real(kind=RKIND), intent(in) :: alpha !< Input: Characteristic length scale of RBFs
    real(kind=RKIND), dimension(2,3) :: planeBasisVectors !< Input: Basis vectors for the interpolation plane
    real(kind=RKIND), dimension(pointCount), intent(out) :: coefficients !< Output: List of coefficients

    integer :: i, j
    integer :: matrixSize

    real(kind=RKIND), dimension(:,:), pointer :: dirichletMatrix
    real(kind=RKIND), dimension(:), pointer :: rhs, coeffs
    integer, dimension(:), pointer :: pivotIndices

    matrixSize = pointCount+3 !! 3 extra space for constant and 2 planar dimensions

    allocate(dirichletMatrix(matrixSize,matrixSize))  
    allocate(rhs(matrixSize))  
    allocate(coeffs(matrixSize))  
    allocate(pivotIndices(matrixSize))  

    dirichletMatrix = 0.0
    rhs = 0.0
    coeffs = 0.0

    call mpas_set_up_scalar_rbf_dirichlet_matrix_and_rhs(pointCount, sourcePoints, destinationPoint, &
      alpha, dirichletMatrix(1:pointCount,1:pointCount), rhs(1:pointCount))

    do i = 1, pointCount
      dirichletMatrix(i,pointCount+1) = 1.0
      dirichletMatrix(i,pointCount+2) = sum(sourcePoints(i,1:3)*planeBasisVectors(1,:))
      dirichletMatrix(i,pointCount+3) = sum(sourcePoints(i,1:3)*planeBasisVectors(2,:))
      dirichletMatrix(pointCount+1:pointCount+3,i) &
        = dirichletMatrix(i,pointCount+1:pointCount+3)
    end do

    rhs(pointCount+1) = 1.0
    rhs(pointCount+2) = sum(destinationPoint(1:3)*planeBasisVectors(1,:))
    rhs(pointCount+3) = sum(destinationPoint(1:3)*planeBasisVectors(2,:))

    ! solve each linear system
    call mpas_legs(dirichletMatrix, matrixSize, rhs, coeffs, pivotIndices)
    coefficients = coeffs(1:pointCount)

    deallocate(dirichletMatrix)  
    deallocate(rhs)  
    deallocate(coeffs)  
    deallocate(pivotIndices) 

  end subroutine mpas_rbf_interp_func_3D_plane_sca_lin_dir_comp_coeffs!}}}

!***********************************************************************
!
!  routine mpas_rbf_interp_func_3D_sca_lin_dir_comp_coeffs
!
!> \brief   MPAS 3D scalar linear coefficients routine
!> \author  Xylar Asay-Davis
!> \date    03/28/13
!> \details 
!>  This routine computes interpolation coefficients in 3D that can be used to
!>   interpolate a number of scalar functions at a given locations. This is useful
!>   if the interpolation location does not change with time, or if several
!>   fields are to be interpolated at a given time step.  (If both the fields
!>   and the interpolation locations vary with time, there is no clear advantage in
!>   using either this method or the method for 2D interpoaltion above; for simplicity
!>   and because we foresee more uses for the method of this subroutine, we have not
!>   implemented a 3D version of the fixed field, variable interpolation location method
!>   as we have in 2D.) Coefficients are produced for handling both Dirichlet (or no)
!>   boundary conditions.  The interpolation is performed with basis functions that are
!>   RBFs plus constant and linear.
!>  Input:
!>   pointCount - the number of "source" points and functionValues supplied
!>   sourcePoints - the location of the "source" points in the 3D space where the values of
!>     the function are known
!>   destinationPoint - the point where the interpolation will be performed
!>   alpha - a constant that give the characteristic length scale of the RBFs,
!>     should be on the order of the distance between points
!>  Output:
!>   coefficients - the coefficients used to interpolate a function with Dirichlet
!>     boundary conditions to the specified destinationPoint
!-----------------------------------------------------------------------
  subroutine mpas_rbf_interp_func_3D_sca_lin_dir_comp_coeffs(pointCount, &!{{{
    sourcePoints, destinationPoint, alpha, coefficients)

    integer, intent(in) :: pointCount !< Input: Number of points
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: sourcePoints !< Input: List of source points
    real(kind=RKIND), dimension(3), intent(in) :: destinationPoint !< Input: Destination point
    real(kind=RKIND), intent(in) :: alpha !< Input: Characteristic length scale for RBFs
    real(kind=RKIND), dimension(pointCount), intent(out) :: coefficients !< Output: List of coefficients

    integer :: i, j
    integer :: matrixSize

    real(kind=RKIND), dimension(:,:), pointer :: dirichletMatrix
    real(kind=RKIND), dimension(:), pointer :: rhs, coeffs
    integer, dimension(:), pointer :: pivotIndices

    matrixSize = pointCount+4 !! 4 extra space for constant and linear in 3D

    allocate(dirichletMatrix(matrixSize,matrixSize))  
    allocate(rhs(matrixSize))  
    allocate(coeffs(matrixSize))  
    allocate(pivotIndices(matrixSize))  

    dirichletMatrix = 0.0
    rhs = 0.0
    coeffs = 0.0

    call mpas_set_up_scalar_rbf_dirichlet_matrix_and_rhs(pointCount, sourcePoints, destinationPoint, &
      alpha, dirichletMatrix(1:pointCount,1:pointCount), rhs(1:pointCount))

    do i = 1, pointCount
      dirichletMatrix(i,pointCount+1) = 1.0
      dirichletMatrix(i,pointCount+2:pointCount+4) = sourcePoints(i,1:3)
      dirichletMatrix(pointCount+1:pointCount+4,i) &
        = dirichletMatrix(i,pointCount+1:pointCount+4)
    end do

    rhs(pointCount+1) = 1.0
    rhs(pointCount+2:pointCount+4) = destinationPoint(1:3)

    ! solve each linear system
    call mpas_legs(dirichletMatrix, matrixSize, rhs, coeffs, pivotIndices)
    coefficients = coeffs(1:pointCount)

    deallocate(dirichletMatrix)  
    deallocate(rhs)  
    deallocate(coeffs)  
    deallocate(pivotIndices) 

  end subroutine mpas_rbf_interp_func_3D_sca_lin_dir_comp_coeffs!}}}

!***********************************************************************
!
!  routine mpas_rbf_interp_func_3D_sca_const_dir_neu_comp_coeffs
!
!> \brief   MPAS 3D scalar constant Dirichlet and Neumann coefficients routine
!> \author  Xylar Asay-Davis
!> \date    03/28/13
!> \details 
!>  This routine computes interpolation coefficients in 3D that can be used to
!>   interpolate a number of scalar functions at a given locations. This is useful
!>   if the interpolation location does not change with time, or if several
!>   fields are to be interpolated at a given time step.  (If both the fields
!>   and the interpolation locations vary with time, there is no clear advantage in
!>   using either this method or the method for 2D interpoaltion above; for simplicity
!>   and because we foresee more uses for the method of this subroutine, we have not
!>   implemented a 3D version of the fixed field, variable interpolation location method
!>   as we have in 2D.) Coefficients are produced for handling both Dirichlet and Neumann
!>   boundary conditions.  The interpolation is performed with basis functions that are
!>   RBFs plus a constant.
!>  Input:
!>   pointCount - the number of "source" points and functionValues supplied
!>   sourcePoints - the location of the "source" points in the 3D space where the values of
!>     the function are known
!>   isInterface - a logical array indicating which of the source points (if any) are at
!>     at the domain interface.  These points and their normals will be used to compute the
!>     neumannCoefficients below
!>   interfaceNormals - a 3D normal vector for each sourcePoint.  These vectors are only used
!>     at points where isInterface == .true., and can take arbitrary values elsewehere.  The
!>     normal vector is used to compute coefficients for the normal derivative of the
!>     interpolating function in order to impose the Neumann Boundary condition
!>   destinationPoint - the point where the interpolation will be performed
!>   alpha - a constant that give the characteristic length scale of the RBFs,
!>     should be on the order of the distance between points
!>  Output:
!>   dirichletCoefficients - the coefficients used to interpolate a function with Dirichlet
!>     boundary conditions to the specified destinationPoint
!>   neumannCoefficients - the coefficients used to interpolate a function with Neumann
!>     boundary conditions to the specified destinationPoint
!-----------------------------------------------------------------------
  subroutine mpas_rbf_interp_func_3D_sca_const_dir_neu_comp_coeffs( &!{{{
    pointCount, sourcePoints, isInterface, interfaceNormals, destinationPoint, &
    alpha, dirichletCoefficients, neumannCoefficients)

    integer, intent(in) :: pointCount !< Input: Number of points
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: sourcePoints !< Input: List of source points
    logical, dimension(pointCount), intent(in) :: isInterface !< Input: Logicals determining if a source point is at an interface
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: interfaceNormals !< Input: Normal vector at interface for each source point
    real(kind=RKIND), dimension(3), intent(in) :: destinationPoint !< Input: Destination point
    real(kind=RKIND), intent(in) :: alpha !< Input: Characteristic length scale of RBFs
    real(kind=RKIND), dimension(pointCount), intent(out) :: dirichletCoefficients !< Output: Coefficients with Dirichlet BCs
    real(kind=RKIND), dimension(pointCount), intent(out) :: neumannCoefficients !< Output: Coefficients with Neumann BCs

    integer :: i, j
    integer :: matrixSize

    real(kind=RKIND), dimension(:,:), pointer :: dirichletMatrix, neumannMatrix
    real(kind=RKIND), dimension(:), pointer :: rhs, rhsCopy, coeffs
    integer, dimension(:), pointer :: pivotIndices

    matrixSize = pointCount+1 !! 1 extra space for constant 

    allocate(dirichletMatrix(matrixSize,matrixSize))  
    allocate(neumannMatrix(matrixSize,matrixSize))  
    allocate(rhs(matrixSize))  
    allocate(rhsCopy(matrixSize))  
    allocate(coeffs(matrixSize))  
    allocate(pivotIndices(matrixSize))  

    dirichletMatrix = 0.0
    neumannMatrix = 0.0
    rhs = 0.0
    rhsCopy = 0.0
    coeffs = 0.0

    call mpas_set_up_scalar_rbf_matrix_and_rhs(pointCount, &
      sourcePoints, isInterface, interfaceNormals, destinationPoint, &
      alpha, dirichletMatrix(1:pointCount,1:pointCount), &
      neumannMatrix(1:pointCount,1:pointCount), rhs(1:pointCount))

    do i = 1, pointCount
      dirichletMatrix(i,pointCount+1) = 1.0
      if(isInterface(i)) then
        neumannMatrix(i,pointCount+1) = 0.0
      else
        neumannMatrix(i,pointCount+1) = dirichletMatrix(i,pointCount+1)
      end if
      dirichletMatrix(pointCount+1,i) = dirichletMatrix(i,pointCount+1)
      neumannMatrix(pointCount+1,i) = neumannMatrix(i,pointCount+1)
    end do

    rhs(pointCount+1) = 1.0

    ! solve each linear system
    rhsCopy = rhs
    call mpas_legs(dirichletMatrix, matrixSize, rhs, coeffs, pivotIndices)
    dirichletCoefficients = coeffs(1:pointCount)

    call mpas_legs(neumannMatrix, matrixSize, rhsCopy, coeffs, pivotIndices)
    neumannCoefficients = coeffs(1:pointCount)

    deallocate(dirichletMatrix)  
    deallocate(neumannMatrix)  
    deallocate(rhs)  
    deallocate(rhsCopy)  
    deallocate(coeffs)  
    deallocate(pivotIndices) 

  end subroutine mpas_rbf_interp_func_3D_sca_const_dir_neu_comp_coeffs!}}}

!***********************************************************************
!
!  routine mpas_rbf_interp_func_3D_plane_sca_lin_dir_neu_comp_coeffs
!
!> \brief   MPAS 3D scalar planar linear Dirichlet and Neumann coefficients routine
!> \author  Xylar Asay-Davis
!> \date    03/28/13
!> \details 
!>  This routine computes interpolation coefficients in a plane in 3D that can be used to
!>   interpolate a number of scalar functions at a given locations. This is useful
!>   if the interpolation location does not change with time, or if several
!>   fields are to be interpolated at a given time step.  (If both the fields
!>   and the interpolation locations vary with time, there is no clear advantage in
!>   using either this method or the method for 2D interpoaltion above; for simplicity
!>   and because we foresee more uses for the method of this subroutine, we have not
!>   implemented a 3D version of the fixed field, variable interpolation location method
!>   as we have in 2D.) Coefficients are produced for handling both Dirichlet and Neumann
!>   boundary conditions.  The interpolation is performed with basis functions that are
!>   RBFs plus constant and linear.  All points are projected into the plane given by the
!>   planeBasisVectors.
!>  Input:
!>   pointCount - the number of "source" points and functionValues supplied
!>   sourcePoints - the location of the "source" points in the 3D space where the values of
!>     the function are known.  The sourcePoints will be projected into the plane given by
!>     planeBasisVectors
!>   isInterface - a logical array indicating which of the source points (if any) are at
!>     at the domain interface.  These points and their normals will be used to compute the
!>     neumannCoefficients below
!>   interfaceNormals - a 3D normal vector for each sourcePoint.  These vectors are only used
!>     at points where isInterface == .true., and can take arbitrary values elsewehere.  The
!>     normal vector is used to compute coefficients for the normal derivative of the
!>     interpolating function in order to impose the Neumann Boundary condition
!>   destinationPoint - the point in 3D where the interpolation will be performed.  The
!>     destinationPoint will be projected into the plane given by planeBasisVectors.
!>   alpha - a constant that give the characteristic length scale of the RBFs,
!>     should be on the order of the distance between points
!>   planeBasisVectors - the basis fectors for the plane where interpolation is performed.
!>     All points are projected into this plane. 
!>  Output:
!>   dirichletCoefficients - the coefficients used to interpolate a function with Dirichlet
!>     boundary conditions to the specified destinationPoint
!>   neumannCoefficients - the coefficients used to interpolate a function with Neumann
!>     boundary conditions to the specified destinationPoint
!-----------------------------------------------------------------------
  subroutine mpas_rbf_interp_func_3D_plane_sca_lin_dir_neu_comp_coeffs( &!{{{
    pointCount, sourcePoints, isInterface, interfaceNormals, destinationPoint, &
    alpha, planeBasisVectors, dirichletCoefficients, neumannCoefficients)

    integer, intent(in) :: pointCount !< Input: Number of points
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: sourcePoints !< Input: List of points
    logical, dimension(pointCount), intent(in) :: isInterface !< Input: List of logicals determining if point is at an interface
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: interfaceNormals !< Input: List of interface normals
    real(kind=RKIND), dimension(3), intent(in) :: destinationPoint !< Input: Destination point
    real(kind=RKIND), intent(in) :: alpha !< Input: Characteristic length scale of RBFs
    real(kind=RKIND), dimension(2,3) :: planeBasisVectors !< Input: Basis vectors for interpolation plane
    real(kind=RKIND), dimension(pointCount), intent(out) :: dirichletCoefficients !< Output: List of Dirichlet coefficients
    real(kind=RKIND), dimension(pointCount), intent(out) :: neumannCoefficients !< Output: List of Neumann coefficients

    integer :: i, j
    integer :: matrixSize

    real(kind=RKIND), dimension(:,:), pointer :: dirichletMatrix, neumannMatrix
    real(kind=RKIND), dimension(:), pointer :: rhs, rhsCopy, coeffs
    integer, dimension(:), pointer :: pivotIndices

    matrixSize = pointCount+3 !! 3 extra space for constant and 2 planar dimensions

    allocate(dirichletMatrix(matrixSize,matrixSize))  
    allocate(neumannMatrix(matrixSize,matrixSize))  
    allocate(rhs(matrixSize))  
    allocate(rhsCopy(matrixSize))  
    allocate(coeffs(matrixSize))  
    allocate(pivotIndices(matrixSize))  

    dirichletMatrix = 0.0
    neumannMatrix = 0.0
    rhs = 0.0
    rhsCopy = 0.0
    coeffs = 0.0

    call mpas_set_up_scalar_rbf_matrix_and_rhs(pointCount, &
      sourcePoints, isInterface, interfaceNormals, destinationPoint, &
      alpha, dirichletMatrix(1:pointCount,1:pointCount), &
      neumannMatrix(1:pointCount,1:pointCount), rhs(1:pointCount))

    do i = 1, pointCount
      dirichletMatrix(i,pointCount+1) = 1.0
      dirichletMatrix(i,pointCount+2) = sum(sourcePoints(i,1:3)*planeBasisVectors(1,:))
      dirichletMatrix(i,pointCount+3) = sum(sourcePoints(i,1:3)*planeBasisVectors(2,:))
      if(isInterface(i)) then
        neumannMatrix(i,pointCount+1) = 0.0
        neumannMatrix(i,pointCount+2) = sum(interfaceNormals(i,1:3)*planeBasisVectors(1,:))
        neumannMatrix(i,pointCount+3) = sum(interfaceNormals(i,1:3)*planeBasisVectors(2,:))
      else
        neumannMatrix(i,pointCount+1:pointCount+3) &
          = dirichletMatrix(i,pointCount+1:pointCount+3)
      end if
      dirichletMatrix(pointCount+1:pointCount+3,i) &
        = dirichletMatrix(i,pointCount+1:pointCount+3)
      neumannMatrix(pointCount+1:pointCount+3,i) &
        = neumannMatrix(i,pointCount+1:pointCount+3)
    end do

    rhs(pointCount+1) = 1.0
    rhs(pointCount+2) = sum(destinationPoint(1:3)*planeBasisVectors(1,:))
    rhs(pointCount+3) = sum(destinationPoint(1:3)*planeBasisVectors(2,:))

    ! solve each linear system
    rhsCopy = rhs
    call mpas_legs(dirichletMatrix, matrixSize, rhs, coeffs, pivotIndices)
    dirichletCoefficients = coeffs(1:pointCount)

    call mpas_legs(neumannMatrix, matrixSize, rhsCopy, coeffs, pivotIndices)
    neumannCoefficients = coeffs(1:pointCount)

    deallocate(dirichletMatrix)  
    deallocate(neumannMatrix)  
    deallocate(rhs)  
    deallocate(rhsCopy)  
    deallocate(coeffs)  
    deallocate(pivotIndices) 

  end subroutine mpas_rbf_interp_func_3D_plane_sca_lin_dir_neu_comp_coeffs!}}}

!***********************************************************************
!
!  routine mpas_rbf_interp_func_3D_sca_lin_dir_neu_comp_coeffs
!
!> \brief   MPAS 3D scalar linear Dirichlet and Neumann coefficients routine
!> \author  Xylar Asay-Davis
!> \date    03/28/13
!> \details 
!>  This routine computes interpolation coefficients in 3D that can be used to
!>   interpolate a number of scalar functions at a given locations. This is useful
!>   if the interpolation location does not change with time, or if several
!>   fields are to be interpolated at a given time step.  (If both the fields
!>   and the interpolation locations vary with time, there is no clear advantage in
!>   using either this method or the method for 2D interpoaltion above; for simplicity
!>   and because we foresee more uses for the method of this subroutine, we have not
!>   implemented a 3D version of the fixed field, variable interpolation location method
!>   as we have in 2D.) Coefficients are produced for handling both Dirichlet and Neumann
!>   boundary conditions.  The interpolation is performed with basis functions that are
!>   RBFs plus constant and linear.
!>  Input:
!>   pointCount - the number of "source" points and functionValues supplied
!>   sourcePoints - the location of the "source" points in the 3D space where the values of
!>     the function are known
!>   isInterface - a logical array indicating which of the source points (if any) are at
!>     at the domain interface.  These points and their normals will be used to compute the
!>     neumannCoefficients below
!>   interfaceNormals - a 3D normal vector for each sourcePoint.  These vectors are only used
!>     at points where isInterface == .true., and can take arbitrary values elsewehere.  The
!>     normal vector is used to compute coefficients for the normal derivative of the
!>     interpolating function in order to impose the Neumann Boundary condition
!>   destinationPoint - the point where the interpolation will be performed
!>   alpha - a constant that give the characteristic length scale of the RBFs,
!>     should be on the order of the distance between points
!>  Output:
!>   dirichletCoefficients - the coefficients used to interpolate a function with Dirichlet
!>     boundary conditions to the specified destinationPoint
!>   neumannCoefficients - the coefficients used to interpolate a function with Neumann
!>     boundary conditions to the specified destinationPoint
!-----------------------------------------------------------------------
  subroutine mpas_rbf_interp_func_3D_sca_lin_dir_neu_comp_coeffs(pointCount, &!{{{
    sourcePoints, isInterface, interfaceNormals, destinationPoint, &
    alpha, dirichletCoefficients, neumannCoefficients)

    integer, intent(in) :: pointCount !< Input: Number of points
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: sourcePoints !< Input: List of points
    logical, dimension(pointCount), intent(in) :: isInterface !< Input: List of logicals determining if point as at an interface
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: interfaceNormals !< Input: List of interface normals
    real(kind=RKIND), dimension(3), intent(in) :: destinationPoint !< Input: Destination point
    real(kind=RKIND), intent(in) :: alpha !< Input: Characteristic length scale of RBFs
    real(kind=RKIND), dimension(pointCount), intent(out) :: dirichletCoefficients !< Output: List of Dirichlet coefficients
    real(kind=RKIND), dimension(pointCount), intent(out) :: neumannCoefficients !< Outut: List of Neumann coefficients

    integer :: i, j
    integer :: matrixSize

    real(kind=RKIND), dimension(:,:), pointer :: dirichletMatrix, neumannMatrix
    real(kind=RKIND), dimension(:), pointer :: rhs, rhsCopy, coeffs
    integer, dimension(:), pointer :: pivotIndices

    matrixSize = pointCount+4 !! 4 extra space for constant and linear in 3D

    allocate(dirichletMatrix(matrixSize,matrixSize))  
    allocate(neumannMatrix(matrixSize,matrixSize))  
    allocate(rhs(matrixSize))  
    allocate(rhsCopy(matrixSize))  
    allocate(coeffs(matrixSize))  
    allocate(pivotIndices(matrixSize))  

    dirichletMatrix = 0.0
    neumannMatrix = 0.0
    rhs = 0.0
    rhsCopy = 0.0
    coeffs = 0.0

    call mpas_set_up_scalar_rbf_matrix_and_rhs(pointCount, &
      sourcePoints, isInterface, interfaceNormals, destinationPoint, &
      alpha, dirichletMatrix(1:pointCount,1:pointCount), &
      neumannMatrix(1:pointCount,1:pointCount), rhs(1:pointCount))

    do i = 1, pointCount
      dirichletMatrix(i,pointCount+1) = 1.0
      dirichletMatrix(i,pointCount+2:pointCount+4) = sourcePoints(i,1:3)
      if(isInterface(i)) then
        neumannMatrix(i,pointCount+1) = 0.0
        neumannMatrix(i,pointCount+2:pointCount+4) = interfaceNormals(i,1:3)
      else
        neumannMatrix(i,pointCount+1:pointCount+4) &
          = dirichletMatrix(i,pointCount+1:pointCount+4)
      end if
      dirichletMatrix(pointCount+1:pointCount+4,i) &
        = dirichletMatrix(i,pointCount+1:pointCount+4)
      neumannMatrix(pointCount+1:pointCount+4,i) &
        = neumannMatrix(i,pointCount+1:pointCount+4)
    end do

    rhs(pointCount+1) = 1.0
    rhs(pointCount+2:pointCount+4) = destinationPoint(1:3)

    ! solve each linear system
    rhsCopy = rhs
    call mpas_legs(dirichletMatrix, matrixSize, rhs, coeffs, pivotIndices)
    dirichletCoefficients = coeffs(1:pointCount)

    call mpas_legs(neumannMatrix, matrixSize, rhsCopy, coeffs, pivotIndices)
    neumannCoefficients = coeffs(1:pointCount)

    deallocate(dirichletMatrix)  
    deallocate(neumannMatrix)  
    deallocate(rhs)  
    deallocate(rhsCopy)  
    deallocate(coeffs)  
    deallocate(pivotIndices) 

  end subroutine mpas_rbf_interp_func_3D_sca_lin_dir_neu_comp_coeffs!}}}

!***********************************************************************
!
!  routine mpas_rbf_interp_func_3D_vec_const_dir_comp_coeffs
!
!> \brief   MPAS 3D vector constant Dirichlet coefficients routine
!> \author  Xylar Asay-Davis
!> \date    03/28/13
!> \details 
!>  This routine computes interpolation coefficients in 3D that can be used to
!>   interpolate a number of vector functions at a given locations. This is useful
!>   if the interpolation location does not change with time, or if several
!>   fields are to be interpolated at a given time step.  (If both the vector fields
!>   and the interpolation locations vary with time, there is no clear advantage in
!>   using either this method or the method for 2D interpoaltion above; for simplicity
!>   and because we foresee more uses for the method of this subroutine, we have not
!>   implemented a 3D version of the fixed field, variable interpolation location method
!>   as we have in 2D.) Coefficients are produced for handling Dirichlet boundary
!>   conditions (or no boundaries).  The interpolation is performed with basis functions
!>   that are RBFs plus a constant.
!>  Input:
!>   pointCount - the number of "source" points and functionValues supplied
!>   sourcePoints - the location of the "source" points in the 3D space where the values of
!>     the function are known
!>   unitVectors - the unit vectors associated with each of the sourcePoints.  Interpolation
!>     is performed by supplying the value of the vector function dotted into each of these unit
!>     vectors.  If multiple unit vectors are supplied at the same sourcePoint, they *must* be
!>     orthogonal for the interpolation to succeed.
!>   destinationPoint - the point where the interpolation will be performed
!>   alpha - a constant that give the characteristic length scale of the RBFs,
!>     should be on the order of the distance between points
!>  Output:
!>   coefficients - the coefficients used to interpolate a function with Dirichlet
!>     boundary conditions to the specified destinationPoint
!-----------------------------------------------------------------------
  subroutine mpas_rbf_interp_func_3D_vec_const_dir_comp_coeffs(pointCount, &!{{{
    sourcePoints, unitVectors, destinationPoint, &
    alpha, coefficients)

    integer, intent(in) :: pointCount !< Input: Number of points
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: sourcePoints !< Input: List of points
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: unitVectors !< Input: List of unit vectors
    real(kind=RKIND), dimension(3), intent(in) :: destinationPoint !< Input: Destination point
    real(kind=RKIND), intent(in) :: alpha !< Input: Characteristic length scale of RBFs
    real(kind=RKIND), dimension(pointCount, 3), intent(out) :: coefficients !< Output: List of coefficients

    integer :: i, j
    integer :: matrixSize

    real(kind=RKIND), dimension(:,:), pointer :: matrix, matrixCopy
    real(kind=RKIND), dimension(:,:), pointer :: rhs, coeffs
    integer, dimension(:), pointer :: pivotIndices

    matrixSize = pointCount+3 ! extra space for constant vector 

    allocate(matrix(matrixSize,matrixSize))  
    allocate(matrixCopy(matrixSize,matrixSize))  
    allocate(rhs(matrixSize,3))  
    allocate(coeffs(matrixSize,3))  
    allocate(pivotIndices(matrixSize))  

    matrix = 0.0
    rhs = 0.0
    coeffs = 0.0

    call mpas_set_up_vector_dirichlet_rbf_matrix_and_rhs(pointCount, 3, &
      sourcePoints, unitVectors, destinationPoint, &
      alpha, matrix(1:pointCount,1:pointCount), rhs(1:pointCount,:))

    do i = 1, pointCount
      matrix(i,pointCount+1:pointCount+3) = unitVectors(i,:)
      matrix(pointCount+1:pointCount+3,i) &
        = matrix(i,pointCount+1:pointCount+3)
    end do
    do i = 1, 3
      rhs(pointCount+i,i) = 1.0 ! the unit vector in the ith direction
    end do

    ! solve each linear system
    do i = 1, 3
      matrixCopy = matrix
      call mpas_legs(matrixCopy, matrixSize, rhs(:,i), coeffs(:,i), pivotIndices)
    end do
    coefficients = coeffs(1:pointCount,:)

    deallocate(matrix)  
    deallocate(matrixCopy)  
    deallocate(rhs)  
    deallocate(coeffs)  
    deallocate(pivotIndices) 

  end subroutine mpas_rbf_interp_func_3D_vec_const_dir_comp_coeffs !}}}

!***********************************************************************
!
!  routine mpas_rbf_interp_func_3D_plane_vec_const_dir_comp_coeffs
!
!> \brief   MPAS 3D vector planar constant Dirichlet coefficients routine
!> \author  Xylar Asay-Davis
!> \date    03/28/13
!> \details 
!>   This routine computes interpolation coefficients in 3D that can be used to
!>    interpolate a number of vector functions at a given locations. This is useful
!>    if the interpolation location does not change with time, or if several
!>    fields are to be interpolated at a given time step.  (If both the vector fields
!>    and the interpolation locations vary with time, there is no clear advantage in
!>    using either this method or the method for 2D interpoaltion above; for simplicity
!>    and because we foresee more uses for the method of this subroutine, we have not
!>    implemented a 3D version of the fixed field, variable interpolation location method
!>    as we have in 2D.) Coefficients are produced for handling Dirichlet boundary
!>    conditions (or no boundaries).  The interpolation is performed with basis functions
!>    that are RBFs plus a constant.
!>   Input:
!>    pointCount - the number of "source" points and functionValues supplied
!>    sourcePoints - the location of the "source" points in the 3D space where the values of
!>      the function are known.  The sourcePoints are projected into the plane given by
!>      planeBasisVectors
!>    unitVectors - the unit vectors associated with each of the sourcePoints.  Interpolation
!>      is performed by supplying the value of the vector function dotted into each of these unit
!>      vectors.  If multiple unit vectors are supplied at the same sourcePoint, they *must* be
!>      orthogonal for the interpolation to succeed.  The unitVectors are projected into the
!>      plane given by planeBasisVectors
!>    destinationPoint - the point where the interpolation will be performed.  The destinationPoint
!>      is projected into the plane given by planeBasisVectors
!>    alpha - a constant that give the characteristic length scale of the RBFs,
!>      should be on the order of the distance between points
!>    planeBasisVectors - the basis fectors for the plane where interpolation is performed.
!>      All points are projected into this plane. 
!>   Output:
!>    coefficients - the coefficients used to interpolate a function with Dirichlet
!>      boundary conditions to the specified destinationPoint
!-----------------------------------------------------------------------
  subroutine mpas_rbf_interp_func_3D_plane_vec_const_dir_comp_coeffs(pointCount, &!{{{
    sourcePoints, unitVectors, destinationPoint, &
    alpha, planeBasisVectors, coefficients)

    integer, intent(in) :: pointCount !< Input: Number of points
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: sourcePoints !< Input: List of points
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: unitVectors !< Input: List of unit vectors
    real(kind=RKIND), dimension(3), intent(in) :: destinationPoint !< Input: Destination point
    real(kind=RKIND), intent(in) :: alpha !< Input: Characteristic length scale of RBFs
    real(kind=RKIND), dimension(2,3) :: planeBasisVectors !< Input: Basis vectors for interpolation plane
    real(kind=RKIND), dimension(pointCount, 3), intent(out) :: coefficients !< Output: List of coefficients

    integer :: i, j
    integer :: matrixSize

    real(kind=RKIND) :: rSquared, rbfValue, unitVectorDotProduct

    real(kind=RKIND), dimension(pointCount,2) :: planarSourcePoints
    real(kind=RKIND), dimension(pointCount,2) :: planarUnitVectors
    real(kind=RKIND), dimension(2) :: planarDestinationPoint

    real(kind=RKIND), dimension(:,:), pointer :: matrix, matrixCopy
    real(kind=RKIND), dimension(:,:), pointer :: rhs, coeffs
    integer, dimension(:), pointer :: pivotIndices

    matrixSize = pointCount+2 ! space for constant vector in plane

    allocate(matrix(matrixSize,matrixSize))  
    allocate(matrixCopy(matrixSize,matrixSize))  
    allocate(rhs(matrixSize,2))  
    allocate(coeffs(matrixSize,2))  
    allocate(pivotIndices(matrixSize)) 

    matrix = 0.0
    rhs = 0.0
    coeffs = 0.0

    do i = 1, pointCount
      planarSourcePoints(i,1) = sum(sourcePoints(i,:)*planeBasisVectors(1,:)) 
      planarSourcePoints(i,2) = sum(sourcePoints(i,:)*planeBasisVectors(2,:)) 
      planarUnitVectors(i,1) = sum(unitVectors(i,:)*planeBasisVectors(1,:)) 
      planarUnitVectors(i,2) = sum(unitVectors(i,:)*planeBasisVectors(2,:)) 
    end do
    planarDestinationPoint(1) = sum(destinationPoint*planeBasisVectors(1,:)) 
    planarDestinationPoint(2) = sum(destinationPoint*planeBasisVectors(2,:)) 

    call mpas_set_up_vector_dirichlet_rbf_matrix_and_rhs(pointCount, 2, &
      planarSourcePoints, planarUnitVectors, planarDestinationPoint, &
      alpha, matrix(1:pointCount,1:pointCount), rhs(1:pointCount,:))

    do i = 1, pointCount
      matrix(i,pointCount+1:pointCount+2) = planarUnitVectors(i,:) 
      matrix(pointCount+1:pointCount+2,i) = matrix(i,pointCount+1:pointCount+2)
    end do
    do i = 1,2 
      rhs(pointCount+i,i) = 1.0 ! the unit vector in the ith direction
    end do

    ! solve each linear system
    matrixCopy = matrix
    call mpas_legs(matrix, matrixSize, rhs(:,1), coeffs(:,1), pivotIndices)
    call mpas_legs(matrixCopy, matrixSize, rhs(:,2), coeffs(:,2), pivotIndices)


    do i = 1,3 
      coefficients(:,i) = planeBasisVectors(1,i)*coeffs(1:pointCount,1) &
        + planeBasisVectors(2,i)*coeffs(1:pointCount,2) 
    end do

    deallocate(matrix)  
    deallocate(matrixCopy)  
    deallocate(rhs)  
    deallocate(coeffs)  
    deallocate(pivotIndices)

  end subroutine mpas_rbf_interp_func_3D_plane_vec_const_dir_comp_coeffs !}}}

!***********************************************************************
!
!  routine mpas_rbf_interp_func_3D_vec_const_tan_neu_comp_coeffs
!
!> \brief   MPAS 3D vector constant tangent Neumann coefficients routine
!> \author  Xylar Asay-Davis
!> \date    03/28/13
!> \details 
!>  This routine computes interpolation coefficients in 3D that can be used to
!>   interpolate a number of vector functions at a given locations. This is useful
!>   if the interpolation location does not change with time, or if several
!>   fields are to be interpolated at a given time step.  (If both the vector fields
!>   and the interpolation locations vary with time, there is no clear advantage in
!>   using either this method or the method for 2D interpoaltion above; for simplicity
!>   and because we foresee more uses for the method of this subroutine, we have not
!>   implemented a 3D version of the fixed field, variable interpolation location method
!>   as we have in 2D.) Coefficients are produced for handling Dirichlet normal /
!>   Neumann tangential boundary conditions (such as free slip).  The interpolation is 
!>   performed with basis functions that are RBFs plus a constant.
!>  Input:
!>   pointCount - the number of "source" points and functionValues supplied
!>   sourcePoints - the location of the "source" points in the 3D space where the values of
!>     the function are known
!>   isTangentToInterface - a logical array indicating which sourcePoints/unitVectors are
!>     tangent to the interface where the boundary condition will be applied.  A Neumann
!>     boundary condition will be applied at these points in these directions.
!>   normalVectorIndex - where isTangentToInterface == .true., the index into unitVectors that
!>     gives the normal vector at the same sourcePoint.  This information is needed to compute
!>     the Neumann boundary condition at this point.
!>   unitVectors - the unit vectors associated with each of the sourcePoints.  Interpolation
!>     is performed by supplying the value of the vector function dotted into each of these unit
!>     vectors.  If multiple unit vectors are supplied at the same sourcePoint, they *must* be
!>     orthogonal for the interpolation to succeed.  A normal vector and two tangential vectors
!>     are needed at each interface point in order to satisfy the Dirichlet normal boundary
!>     condition and the Neumann tangential boundary conditions at these points.
!>   destinationPoint - the point where the interpolation will be performed
!>   alpha - a constant that give the characteristic length scale of the RBFs,
!>     should be on the order of the distance between points
!>  Output:
!>   coefficients - the coefficients used to interpolate a function with Dirichlet
!>     boundary conditions to the specified destinationPoint
!-----------------------------------------------------------------------
  subroutine mpas_rbf_interp_func_3D_vec_const_tan_neu_comp_coeffs(pointCount, &!{{{
    sourcePoints, isTangentToInterface, normalVectorIndex, unitVectors, destinationPoint, &
    alpha, coefficients)

    integer, intent(in) :: pointCount !< Input: Number of points
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: sourcePoints !< Input: List of points
    logical, dimension(pointCount), intent(in) :: isTangentToInterface !< Input: List of logicals determining if point is tangent to interface
    integer, dimension(pointCount), intent(in) :: normalVectorIndex !< Input: Index of for normal vectors
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: unitVectors !< Input: List of unit vectors
    real(kind=RKIND), dimension(3), intent(in) :: destinationPoint !< Input: Destination point
    real(kind=RKIND), intent(in) :: alpha !< Input: Characteristic length scale of RBFs
    real(kind=RKIND), dimension(pointCount, 3), intent(out) :: coefficients !< Output: List of coefficients

    integer :: i, j
    integer :: matrixSize

    real(kind=RKIND), dimension(:,:), pointer :: matrix, matrixCopy
    real(kind=RKIND), dimension(:,:), pointer :: rhs, coeffs
    integer, dimension(:), pointer :: pivotIndices

    matrixSize = pointCount+3 ! extra space for constant vector 

    allocate(matrix(matrixSize,matrixSize))  
    allocate(matrixCopy(matrixSize,matrixSize))  
    allocate(rhs(matrixSize,3))  
    allocate(coeffs(matrixSize,3))  
    allocate(pivotIndices(matrixSize))  

    matrix = 0.0
    rhs = 0.0
    coeffs = 0.0

    call mpas_set_up_vector_free_slip_rbf_matrix_and_rhs(pointCount, 3, &
      sourcePoints, isTangentToInterface, normalVectorIndex, unitVectors, destinationPoint, &
      alpha, matrix(1:pointCount,1:pointCount), rhs(1:pointCount,:))

    do i = 1, pointCount
      matrix(pointCount+1:pointCount+3,i) = unitVectors(i,:)
      if(.not. isTangentToInterface(i)) then
        matrix(i,pointCount+1:pointCount+3) = matrix(pointCount+1:pointCount+3,i)
      end if
    end do
    do i = 1, 3
      rhs(pointCount+i,i) = 1.0 ! the unit vector in the ith direction
    end do

    ! solve each linear system
    do i = 1, 3
      matrixCopy = matrix
      call mpas_legs(matrixCopy, matrixSize, rhs(:,i), coeffs(:,i), pivotIndices)
    end do
    coefficients = coeffs(1:pointCount,:)

    deallocate(matrix)  
    deallocate(matrixCopy)  
    deallocate(rhs)  
    deallocate(coeffs)  
    deallocate(pivotIndices)

  end subroutine mpas_rbf_interp_func_3D_vec_const_tan_neu_comp_coeffs !}}}

!***********************************************************************
!
!  routine mpas_rbf_interp_func_3D_plane_vec_const_tan_neu_comp_coeffs
!
!> \brief   MPAS 3D vector planar constant tangent Neumann coefficients routine
!> \author  Xylar Asay-Davis
!> \date    03/28/13
!> \details 
!>  This routine computes interpolation coefficients in 3D that can be used to
!>   interpolate a number of vector functions at a given locations. This is useful
!>   if the interpolation location does not change with time, or if several
!>   fields are to be interpolated at a given time step.  (If both the vector fields
!>   and the interpolation locations vary with time, there is no clear advantage in
!>   using either this method or the method for 2D interpoaltion above; for simplicity
!>   and because we foresee more uses for the method of this subroutine, we have not
!>   implemented a 3D version of the fixed field, variable interpolation location method
!>   as we have in 2D.) Coefficients are produced for handling Dirichlet normal /
!>   Neumann tangential boundary conditions (such as free slip).  The interpolation is 
!>   performed with basis functions that are RBFs plus a constant.
!>  Input:
!>   pointCount - the number of "source" points and functionValues supplied
!>   sourcePoints - the location of the "source" points in the 3D space where the values of
!>     the function are known.  The sourcePoints are projected into the plane given by
!>     planeBasisVectors
!>   isTangentToInterface - a logical array indicating which sourcePoints/unitVectors are
!>     tangent to the interface where the boundary condition will be applied.  A Neumann
!>     boundary condition will be applied at these points in these directions.
!>   normalVectorIndex - where isTangentToInterface == .true., the index into unitVectors that
!>     gives the normal vector at the same sourcePoint.  This information is needed to compute
!>     the Neumann boundary condition at this point.
!>   unitVectors - the unit vectors associated with each of the sourcePoints.  Interpolation
!>     is performed by supplying the value of the vector function dotted into each of these unit
!>     vectors.  If multiple unit vectors are supplied at the same sourcePoint, they *must* be
!>     orthogonal for the interpolation to succeed.  A normal vector and two tangential vectors
!>     are needed at each interface point in order to satisfy the Dirichlet normal boundary
!>     condition and the Neumann tangential boundary conditions at these points. The unitVectors
!>     are projected into the plane given by planeBasisVectors
!>   destinationPoint - the point where the interpolation will be performed.  The destinationPoint
!>     is projected into the plane given by planeBasisVectors
!>   alpha - a constant that give the characteristic length scale of the RBFs,
!>     should be on the order of the distance between points
!>   planeBasisVectors - the basis fectors for the plane where interpolation is performed.
!>     All points are projected into this plane. 
!>  Output:
!>   coefficients - the coefficients used to interpolate a function with Dirichlet
!>     boundary conditions to the specified destinationPoint
!-----------------------------------------------------------------------

  subroutine mpas_rbf_interp_func_3D_plane_vec_const_tan_neu_comp_coeffs(&!{{{
    pointCount, sourcePoints, isTangentToInterface, normalVectorIndex, unitVectors, &
    destinationPoint, alpha, planeBasisVectors, coefficients)

    integer, intent(in) :: pointCount !< input: Number of points
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: sourcePoints !< Input: List of points
    logical, dimension(pointCount), intent(in) :: isTangentToInterface !< Input: List of logicals determining if point is tangent to interface
    integer, dimension(pointCount), intent(in) :: normalVectorIndex !< Input: Index for normal vectors
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: unitVectors !< Input: List of unit vectors
    real(kind=RKIND), dimension(3), intent(in) :: destinationPoint !< Input: Destination point
    real(kind=RKIND), intent(in) :: alpha !< Input: Characteristic length scale of RBFs
    real(kind=RKIND), dimension(2,3), intent(in) :: planeBasisVectors !< Input: Basis vectors for interpolation plane
    real(kind=RKIND), dimension(pointCount, 3), intent(out) :: coefficients !< Output: List of coefficients

    integer :: i, j
    integer :: matrixSize

    real(kind=RKIND), dimension(pointCount,2) :: planarSourcePoints
    real(kind=RKIND), dimension(pointCount,2) :: planarUnitVectors
    real(kind=RKIND), dimension(2) :: planarDestinationPoint

    real(kind=RKIND), dimension(:,:), pointer :: matrix, matrixCopy
    real(kind=RKIND), dimension(:,:), pointer :: rhs, coeffs
    integer, dimension(:), pointer :: pivotIndices

    matrixSize = pointCount+2 ! space for constant vector in plane

    allocate(matrix(matrixSize,matrixSize))  
    allocate(matrixCopy(matrixSize,matrixSize))  
    allocate(rhs(matrixSize,2))  
    allocate(coeffs(matrixSize,2))  
    allocate(pivotIndices(matrixSize)) 

    matrix = 0.0
    rhs = 0.0
    coeffs = 0.0

    do i = 1, pointCount
      planarSourcePoints(i,1) = sum(sourcePoints(i,:)*planeBasisVectors(1,:)) 
      planarSourcePoints(i,2) = sum(sourcePoints(i,:)*planeBasisVectors(2,:)) 
      planarUnitVectors(i,1) = sum(unitVectors(i,:)*planeBasisVectors(1,:)) 
      planarUnitVectors(i,2) = sum(unitVectors(i,:)*planeBasisVectors(2,:)) 
    end do
    planarDestinationPoint(1) = sum(destinationPoint*planeBasisVectors(1,:)) 
    planarDestinationPoint(2) = sum(destinationPoint*planeBasisVectors(2,:)) 
    call mpas_set_up_vector_free_slip_rbf_matrix_and_rhs(pointCount, 2, &
      planarSourcePoints, isTangentToInterface, normalVectorIndex, planarUnitVectors, &
      planarDestinationPoint, alpha, matrix(1:pointCount,1:pointCount), rhs(1:pointCount,:))

    do i = 1, pointCount
      matrix(pointCount+1,i) = sum(unitVectors(i,:)*planeBasisVectors(1,:))
      matrix(pointCount+2,i) = sum(unitVectors(i,:)*planeBasisVectors(2,:))
      if(.not. isTangentToInterface(i)) then
        matrix(i,pointCount+1:pointCount+2) = matrix(pointCount+1:pointCount+2,i)
      end if
    end do
    do i = 1,2 
      rhs(pointCount+i,i) = 1.0 ! the unit vector in the ith direction
    end do

    ! solve each linear system
    matrixCopy = matrix
    call mpas_legs(matrix, matrixSize, rhs(:,1), coeffs(:,1), pivotIndices)
    call mpas_legs(matrixCopy, matrixSize, rhs(:,2), coeffs(:,2), pivotIndices)

    coefficients(:,1) = planeBasisVectors(1,1)*coeffs(1:pointCount,1) &
      + planeBasisVectors(2,1)*coeffs(1:pointCount,2) 
    coefficients(:,2) = planeBasisVectors(1,2)*coeffs(1:pointCount,1) &
      + planeBasisVectors(2,2)*coeffs(1:pointCount,2) 
    coefficients(:,3) = planeBasisVectors(1,3)*coeffs(1:pointCount,1) &
      + planeBasisVectors(2,3)*coeffs(1:pointCount,2) 

    deallocate(matrix)  
    deallocate(matrixCopy)  
    deallocate(rhs)  
    deallocate(coeffs)  
    deallocate(pivotIndices)

   end subroutine mpas_rbf_interp_func_3D_plane_vec_const_tan_neu_comp_coeffs !}}}


!!!!!!!!!!!!!!!!!!!!!
! private subroutines
!!!!!!!!!!!!!!!!!!!!!

!***********************************************************************
!
!  function evaluate_rbf
!
!> \brief   MPAS RBF Evaluation function
!> \author  Xylar Asay-Davis
!> \date    03/28/13
!> \details 
!> This function evaluates an RBF and returns the value.
!-----------------------------------------------------------------------
  function evaluate_rbf(rSquared) result(rbfValue)!{{{
    real(kind=RKIND), intent(in) :: rSquared !< Input: Squared value of r
    real(kind=RKIND) :: rbfValue

    ! inverse multiquadratic
    rbfValue = 1/sqrt(1 + rSquared)

  end function evaluate_rbf!}}}

!***********************************************************************
!
!  routine mpas_evaluate_rbf_and_deriv
!
!> \brief   MPAS RBF Evaluation and derivative routine
!> \author  Xylar Asay-Davis
!> \date    03/28/13
!> \details 
!> This routine computes the value and derivative of a RBF.
!-----------------------------------------------------------------------
  subroutine mpas_evaluate_rbf_and_deriv(rSquared, rbfValue, rbfDerivOverR)!{{{
    real(kind=RKIND), intent(in) :: rSquared !< Input: Squared value of R
    real(kind=RKIND), intent(out) :: rbfValue !< Output: Value of RBF
    real(kind=RKIND), intent(out) :: rbfDerivOverR  !< Outut: Derivative of RBF over R

    ! inverse multiquadratic
    rbfValue = 1/sqrt(1 + rSquared)
    rbfDerivOverR = -rbfValue**3

  end subroutine mpas_evaluate_rbf_and_deriv!}}}

!***********************************************************************
!
!  routine mpas_evaluate_rbf_and_derivs
!
!> \brief   MPAS RBF Evaluation and first and second derivative routine
!> \author  Xylar Asay-Davis
!> \date    03/28/13
!> \details 
!> This routine computes the value and the first two derivatives of a RBF.
!-----------------------------------------------------------------------
  subroutine mpas_evaluate_rbf_and_derivs(rSquared, rbfValue, rbfDerivOverR, rbfSecondDeriv)!{{{
    real(kind=RKIND), intent(in) :: rSquared !< Input: Squared value of R
    real(kind=RKIND), intent(out) :: rbfValue !< Output: Value of RBF
    real(kind=RKIND), intent(out) :: rbfDerivOverR !< Output: Value of first derivative of RBF
    real(kind=RKIND), intent(out) :: rbfSecondDeriv !< Output: Value of second derivative of RBF

    ! inverse multiquadratic
    rbfValue = 1/sqrt(1 + rSquared)
    rbfDerivOverR = -rbfValue**3
    rbfSecondDeriv = (2*rSquared-1)*rbfValue**5

  end subroutine mpas_evaluate_rbf_and_derivs!}}}

!***********************************************************************
!
!  routine mpas_set_up_scalar_rbf_dirichlet_matrix_and_rhs
!
!> \brief   MPAS RBF Scalar Matrix and RHS setup routine
!> \author  Xylar Asay-Davis
!> \date    03/28/13
!> \details 
!> This routine sets up the matrix and RHS for scalar Dirichlet RBF interpolation.
!-----------------------------------------------------------------------
  subroutine mpas_set_up_scalar_rbf_dirichlet_matrix_and_rhs(pointCount, sourcePoints, destinationPoint, &!{{{
    alpha, dirichletMatrix, rhs)

    integer, intent(in) :: pointCount !< Input: Number of points
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: sourcePoints !< Input: List of points
    real(kind=RKIND), dimension(3), intent(in) :: destinationPoint !< Input: Destination point
    real(kind=RKIND), intent(in) :: alpha !< Input: Characteristic length scale of RBF
    real(kind=RKIND), dimension(pointCount,pointCount), intent(out) :: dirichletMatrix !< Output: Matrix
    real(kind=RKIND), dimension(pointCount), intent(out) :: rhs !< Output: Right hand side

    integer :: i, j

    real(kind=RKIND) :: rSquared, rbfValue

    do j = 1, pointCount
      do i = 1, pointCount
        rSquared = sum((sourcePoints(i,:)-sourcePoints(j,:))**2)/alpha**2
        rbfValue = evaluate_rbf(rSquared)
        dirichletMatrix(i,j) = rbfValue
      end do
    end do

    do j = 1, pointCount
      rSquared = sum((destinationPoint-sourcePoints(j,:))**2)/alpha**2
      rhs(j) = evaluate_rbf(rSquared)
    end do

  end subroutine mpas_set_up_scalar_rbf_dirichlet_matrix_and_rhs!}}}

!***********************************************************************
!
!  routine mpas_set_up_scalar_rbf_matrix_and_rhs
!
!> \brief   MPAS RBF Scalar Matrix and RHS setup routine
!> \author  Xylar Asay-Davis
!> \date    03/28/13
!> \details 
!> This routine sets up the matrix and RHS for scalar Dirichlet and Neumann RBF interpolation.
!-----------------------------------------------------------------------
  subroutine mpas_set_up_scalar_rbf_matrix_and_rhs(pointCount, &!{{{
    sourcePoints, isInterface, interfaceNormals, destinationPoint, &
    alpha, dirichletMatrix, neumannMatrix, rhs)

    integer, intent(in) :: pointCount !< Input: Number of points
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: sourcePoints !< Input: List of points
    logical, dimension(pointCount), intent(in) :: isInterface !< Input: Logicals determining if point is an interface
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: interfaceNormals !< Input: Normals at interfaces
    real(kind=RKIND), dimension(3), intent(in) :: destinationPoint !< Input: Destination point
    real(kind=RKIND), intent(in) :: alpha !< Input: Characteristic length scale of RBF
    real(kind=RKIND), dimension(pointCount,pointCount), intent(out) :: dirichletMatrix !< Output: Dirichlet Matrix
    real(kind=RKIND), dimension(pointCount,pointCount), intent(out) :: neumannMatrix !< Output: Neumann Matrix
    real(kind=RKIND), dimension(pointCount), intent(out) :: rhs !< Output: Right hand side

    integer :: i, j

    real(kind=RKIND) :: rSquared, rbfValue, rbfDerivOverR, normalDotX

    do j = 1, pointCount
      if(isInterface(j)) then
        do i = 1, pointCount 
          rSquared = sum((sourcePoints(i,:)-sourcePoints(j,:))**2)/alpha**2
          normalDotX = sum(interfaceNormals(j,:) &
            * (sourcePoints(j,:)-sourcePoints(i,:)))
          call mpas_evaluate_rbf_and_deriv(rSquared, rbfValue, rbfDerivOverR)
          rbfDerivOverR = rbfDerivOverR/alpha**2
          dirichletMatrix(i,j) = rbfValue
          neumannMatrix(i,j) = rbfDerivOverR*normalDotX
        end do
      else
        do i = 1, pointCount
          rSquared = sum((sourcePoints(i,:)-sourcePoints(j,:))**2)/alpha**2
          rbfValue = evaluate_rbf(rSquared)
          dirichletMatrix(i,j) = rbfValue
          neumannMatrix(i,j) = rbfValue
        end do
      end if
    end do

    do j = 1, pointCount
      rSquared = sum((destinationPoint-sourcePoints(j,:))**2)/alpha**2
      rhs(j) = evaluate_rbf(rSquared)
    end do

  end subroutine mpas_set_up_scalar_rbf_matrix_and_rhs!}}}

!***********************************************************************
!
!  routine mpas_set_up_vector_dirichlet_rbf_matrix_and_rhs
!
!> \brief   MPAS RBF Vector Matrix and RHS setup routine
!> \author  Xylar Asay-Davis
!> \date    03/28/13
!> \details 
!> This routine sets up the matrix and RHS for vector Dirichlet RBF interpolation.
!-----------------------------------------------------------------------
  subroutine mpas_set_up_vector_dirichlet_rbf_matrix_and_rhs(pointCount, dimensions, &!{{{
    sourcePoints, unitVectors, destinationPoint, &
    alpha, matrix, rhs)

    integer, intent(in) :: pointCount !< Input: Number of points
    integer, intent(in) :: dimensions !< Input: Number of dimensions
    real(kind=RKIND), dimension(pointCount,dimensions), intent(in) :: sourcePoints !< Input: List of points
    real(kind=RKIND), dimension(pointCount,dimensions), intent(in) :: unitVectors !< Input: List of unit vectors
    real(kind=RKIND), dimension(dimensions), intent(in) :: destinationPoint !< Input: Destination point
    real(kind=RKIND), intent(in) :: alpha !< Input: Characteristic length scale of RBFs
    real(kind=RKIND), dimension(pointCount,pointCount), intent(out) :: matrix !< Output: Matrix
    real(kind=RKIND), dimension(pointCount,dimensions), intent(out) :: rhs !< Output: Right hand side

    integer :: i, j

    real(kind=RKIND) :: rSquared, rbfValue, unitVectorDotProduct

    do j = 1, pointCount
      do i = j, pointCount
        rSquared = sum((sourcePoints(i,:)-sourcePoints(j,:))**2)/alpha**2
        rbfValue = evaluate_rbf(rSquared)
        unitVectorDotProduct = sum(unitVectors(i,:)*unitVectors(j,:))
        matrix(i,j) = rbfValue*unitVectorDotProduct
        matrix(j,i) = matrix(i,j)
      end do
    end do

    do j = 1, pointCount
      rSquared = sum((destinationPoint-sourcePoints(j,:))**2)/alpha**2
      rhs(j,:) = evaluate_rbf(rSquared)*unitVectors(j,:)
    end do

  end subroutine mpas_set_up_vector_dirichlet_rbf_matrix_and_rhs!}}}

!***********************************************************************
!
!  routine mpas_set_up_vector_free_slip_rbf_matrix_and_rhs
!
!> \brief   MPAS RBF Vector Matrix and RHS setup routine
!> \author  Xylar Asay-Davis
!> \date    03/28/13
!> \details 
!> This routine sets up the matrix and RHS for vector Free Slip RBF interpolation.
!-----------------------------------------------------------------------
  subroutine mpas_set_up_vector_free_slip_rbf_matrix_and_rhs(pointCount, dimensions, &!{{{
    sourcePoints, isTangentToInterface, normalVectorIndex, unitVectors, destinationPoint, &
    alpha, matrix, rhs)

    integer, intent(in) :: pointCount !< Input: Number of points
    integer, intent(in) :: dimensions !< Input: Number of dimensions
    real(kind=RKIND), dimension(pointCount,dimensions), intent(in) :: sourcePoints !< Input: List of points
    logical, dimension(pointCount), intent(in) :: isTangentToInterface !< Input: Logical to determine if point is tangent to interface
    integer, dimension(pointCount), intent(in) :: normalVectorIndex !< Input: Index to normal vector
    real(kind=RKIND), dimension(pointCount,dimensions), intent(in) :: unitVectors !< Input: List of unit vectors
    real(kind=RKIND), dimension(dimensions), intent(in) :: destinationPoint !< Input: Destination point
    real(kind=RKIND), intent(in) :: alpha !< Input: Characteristic length scale of RBF
    real(kind=RKIND), dimension(pointCount,pointCount), intent(out) :: matrix !< Output: Matrix
    real(kind=RKIND), dimension(pointCount,dimensions), intent(out) :: rhs !< Output: Right hand side

    integer :: i, j

    real(kind=RKIND) :: rSquared, rbfValue, rbfDerivOverR, normalVector(dimensions), &
      normalDotX, unitVectorDotProduct

    do j = 1, pointCount
      if(isTangentToInterface(j)) then
         do i = 1, pointCount
          rSquared = sum((sourcePoints(i,:)-sourcePoints(j,:))**2)/alpha**2
          normalVector = unitVectors(normalVectorIndex(j),:) 
          normalDotX = sum(normalVector * (sourcePoints(j,:)-sourcePoints(i,:)))
          call mpas_evaluate_rbf_and_deriv(rSquared, rbfValue, rbfDerivOverR)
          rbfDerivOverR = rbfDerivOverR/alpha**2
          unitVectorDotProduct = sum(unitVectors(i,:)*unitVectors(j,:))
          matrix(i,j) = rbfDerivOverR*normalDotX*unitVectorDotProduct
        end do
      else
        do i = 1, pointCount
          rSquared = sum((sourcePoints(i,:)-sourcePoints(j,:))**2)/alpha**2
          rbfValue = evaluate_rbf(rSquared)
          unitVectorDotProduct = sum(unitVectors(i,:)*unitVectors(j,:))
          matrix(i,j) = rbfValue*unitVectorDotProduct
        end do
      end if
    end do

    do j = 1, pointCount
      rSquared = sum((destinationPoint-sourcePoints(j,:))**2)/alpha**2
      rhs(j,:) = evaluate_rbf(rSquared)*unitVectors(j,:)
    end do

  end subroutine mpas_set_up_vector_free_slip_rbf_matrix_and_rhs!}}}

!***********************************************************************
!
!  routine mpas_unit_vec_in_r3
!
!> \brief   MPAS 3D unit vector routine
!> \author  Xylar Asay-Davis
!> \date    03/28/13
!> \details 
!> This routine creates a unit vector out of an input point.
!-----------------------------------------------------------------------
  subroutine mpas_unit_vec_in_r3(xin)!{{{
    implicit none
    real (kind=RKIND), intent(inout) :: xin(3) !< Input/Output: Vector and unit vector
    real (kind=RKIND) :: mag
    mag = sqrt(xin(1)**2+xin(2)**2+xin(3)**2)
    xin(:) = xin(:) / mag
  end subroutine mpas_unit_vec_in_r3!}}}

!***********************************************************************
!
!  routine mpas_cross_product_in_r3
!
!> \brief   MPAS 3D cross product routine
!> \author  Xylar Asay-Davis
!> \date    03/28/13
!> \details 
!> This routine computes the cross product of two input vectors.
!-----------------------------------------------------------------------
  subroutine mpas_cross_product_in_r3(p_1,p_2,p_out)!{{{
    real (kind=RKIND), intent(in) :: p_1 (3) !< Input: Vector 1
    real (kind=RKIND), intent(in) :: p_2 (3) !< Input: Vector 2
    real (kind=RKIND), intent(out) :: p_out (3) !< Output: Cross product of vector 1 and vector 2

    p_out(1) = p_1(2)*p_2(3)-p_1(3)*p_2(2)
    p_out(2) = p_1(3)*p_2(1)-p_1(1)*p_2(3)
    p_out(3) = p_1(1)*p_2(2)-p_1(2)*p_2(1)
  end subroutine mpas_cross_product_in_r3!}}}

! Updated 10/24/2001.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!   Program 4.3   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
! Please Note:                                                          !
!                                                                       !
! (1) This computer program is written by Tao Pang in conjunction with  !
!     his book, "An Introduction to Computational Physics," published   !
!     by Cambridge University Press in 1997.                            !
!                                                                       !
! (2) No warranties, express or implied, are made for this program.     !
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!PROGRAM EX43
!
!
! An example of solving linear equation set A(N,N)*X(N) = B(N)
! with the partial-pivoting Gaussian elimination scheme.  The
! numerical values are for the Wheatstone bridge example discussed
! in Section 4.1 in the book with all resistances being 100 ohms
! and the voltage 200 volts.
!
!  IMPLICIT NONE
!  INTEGER, PARAMETER :: N=3
!  INTEGER :: I,J
!  INTEGER, DIMENSION (N) :: INDX
!  REAL, DIMENSION (N) :: X,B
!  REAL, DIMENSION (N,N) :: A
!  DATA B /200.0,0.0,0.0/, &
!       ((A(I,J), J=1,N),I=1,N) /100.0,100.0,100.0,-100.0, &
!                         300.0,-100.0,-100.0,-100.0, 300.0/
!
!  call mpas_legs (A,N,B,X,INDX)
!
!  WRITE (6, "(F16.8)") (X(I), I=1,N)
!END PROGRAM EX43

!***********************************************************************
!
!  routine mpas_legs
!
!> \brief   MPAS Gaussian elimination solver routine
!> \author  Xylar Asay-Davis
!> \date    03/28/13
!> \details 
!> This routine solves the equation A(N,N)*X(N) = B(N) with the partial-pivoting
!> Gaussian Elimination scheme. Copyright (c) Tao Pang 2001.
!-----------------------------------------------------------------------
subroutine mpas_legs (A,N,B,X,INDX)!{{{

  IMPLICIT NONE
  integer, INTENT (IN) :: N !< Input: Size of matrix and vectors
  integer, INTENT (OUT), DIMENSION (N) :: INDX !< Output: Pivot vector
  real(kind=RKIND), INTENT (INOUT), DIMENSION (N,N) :: A !< Input/Output: Matrix
  real(kind=RKIND), INTENT (INOUT), DIMENSION (N) :: B !< Input/Output: Right hand side vector
  real(kind=RKIND), INTENT (OUT), DIMENSION (N) :: X !< Output: Solution

  integer :: I,J 
!
  CALL elgs (A,N,INDX)
!
  DO I = 1, N-1
    DO J = I+1, N
      B(INDX(J)) = B(INDX(J))-A(INDX(J),I)*B(INDX(I))
    END DO
  END DO
!
  X(N) = B(INDX(N))/A(INDX(N),N)
  DO I = N-1, 1, -1
    X(I) = B(INDX(I))
    DO J = I+1, N
      X(I) = X(I)-A(INDX(I),J)*X(J)
    END DO
    X(I) =  X(I)/A(INDX(I),I)
  END DO
!
END subroutine mpas_legs!}}}
!



! Updated 10/24/2001.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!   Program 4.4   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
! Please Note:                                                          !
!                                                                       !
! (1) This computer program is written by Tao Pang in conjunction with  !
!     his book, "An Introduction to Computational Physics," published   !
!     by Cambridge University Press in 1997.                            !
!                                                                       !
! (2) No warranties, express or implied, are made for this program.     !
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!***********************************************************************
!
!  routine migs
!
!> \brief   Matrix inversion routine
!> \author  Xylar Asay-Davis
!> \date    03/28/13
!> \details 
!> This routine inverts the matrix A(N,N) and stores it in X(N,B)
!> Copyright (c) Tao Pang 2001.
!-----------------------------------------------------------------------
subroutine migs (A,N,X,INDX)!{{{
  IMPLICIT NONE
  integer, INTENT (IN) :: N !< Input: Size of matrix and inverse
  integer, INTENT (OUT), DIMENSION (N) :: INDX !< Output: Pivot vector
  real(kind=RKIND), INTENT (INOUT), DIMENSION (N,N):: A !< Input/Output: Matrix to invert
  real(kind=RKIND), INTENT (OUT), DIMENSION (N,N):: X !< Output: Inverse of Matrix
  real(kind=RKIND), DIMENSION (N,N) :: B
  integer :: I,J,K
!
  DO I = 1, N
    DO J = 1, N
      B(I,J) = 0.0
    END DO
  END DO
  DO I = 1, N
    B(I,I) = 1.0
  END DO
!
  CALL elgs (A,N,INDX)
!
  DO I = 1, N-1
    DO J = I+1, N
      DO K = 1, N
        B(INDX(J),K) = B(INDX(J),K)-A(INDX(J),I)*B(INDX(I),K)
      END DO
    END DO
  END DO
!
  DO I = 1, N
    X(N,I) = B(INDX(N),I)/A(INDX(N),N)
    DO J = N-1, 1, -1
      X(J,I) = B(INDX(J),I)
      DO K = J+1, N
        X(J,I) = X(J,I)-A(INDX(J),K)*X(K,I)
      END DO
      X(J,I) =  X(J,I)/A(INDX(J),J)
    END DO
  END DO
END subroutine migs!}}}

!***********************************************************************
!
!  routine elgs
!
!> \brief   Partial-pivoting Gaussian elimination routine
!> \author  Xylar Asay-Davis
!> \date    03/28/13
!> \details 
!> This routine performs the partial-pivoting Gaussian elimination.
!> Copyright (c) Tao Pang 2001.
!-----------------------------------------------------------------------

subroutine elgs (A,N,INDX)!{{{
!
! subroutine to perform the partial-pivoting Gaussian elimination.
! A(N,N) is the original matrix in the input and transformed matrix
! plus the pivoting element ratios below the diagonal in the output.
! INDX(N) records the pivoting order.  Copyright (c) Tao Pang 2001.
!
  IMPLICIT NONE
  integer, INTENT (IN) :: N !< Input: Size of matrix
  integer, INTENT (OUT), DIMENSION (N) :: INDX !< Output: Pivot vector
  real(kind=RKIND), INTENT (INOUT), DIMENSION (N,N) :: A !< Input/Output: Matrix and solution
  integer :: I,J,K,ITMP
  real(kind=RKIND) :: C1,PI,PI1,PJ
  real(kind=RKIND), DIMENSION (N) :: C
!
! Initialize the index
!
  DO I = 1, N
    INDX(I) = I
  END DO
!
! Find the rescaling factors, one from each row
!
  DO I = 1, N
    C1= 0.0
    DO J = 1, N
      !C1 = AMAX1(C1,ABS(A(I,J)))
      C1 = MAX(C1,ABS(A(I,J)))
    END DO
    C(I) = C1
  END DO
!
! Search the pivoting (largest) element from each column
!
  DO J = 1, N-1
    PI1 = 0.0
    DO I = J, N
      PI = ABS(A(INDX(I),J))/C(INDX(I))
      IF (PI.GT.PI1) THEN
        PI1 = PI
        K   = I
      ENDIF
    END DO
!
! Interchange the rows via INDX(N) to record pivoting order
!
    ITMP    = INDX(J)
    INDX(J) = INDX(K)
    INDX(K) = ITMP
    DO I = J+1, N
      PJ  = A(INDX(I),J)/A(INDX(J),J)
!
! Record pivoting ratios below the diagonal
!
      A(INDX(I),J) = PJ
!
! Modify other elements accordingly
!
      DO K = J+1, N
        A(INDX(I),K) = A(INDX(I),K)-PJ*A(INDX(J),K)
      END DO
    END DO
  END DO
!
END subroutine elgs!}}}

end module mpas_rbf_interpolation

