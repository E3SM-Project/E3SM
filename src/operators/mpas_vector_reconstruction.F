! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!***********************************************************************
!
!  mpas_vector_reconstruction
!
!> \brief   MPAS Vector reconstruction module
!> \author  Xylar Asay-Davis, Todd Ringler
!> \date    03/28/13
!> \details 
!> This module provides routines for performing vector reconstruction from edges to cell centers.
!
!-----------------------------------------------------------------------
module mpas_vector_reconstruction

  use mpas_grid_types
  use mpas_configure
  use mpas_constants
  use mpas_rbf_interpolation

  implicit none

  public :: mpas_init_reconstruct, mpas_reconstruct

  contains

!***********************************************************************
!
!  routine mpas_init_reconstruct
!
!> \brief   MPAS Vector reconstruction initialization routine
!> \author  Xylar Asay-Davis, Todd Ringler
!> \date    03/28/13
!> \details 
!>  Purpose: pre-compute coefficients used by the reconstruct() routine
!>  Input: grid meta data
!>  Output: grid % coeffs_reconstruct - coefficients used to reconstruct 
!>                                      velocity vectors at cell centers 
!-----------------------------------------------------------------------
  subroutine mpas_init_reconstruct(grid)!{{{

    implicit none

    type (mesh_type), intent(inout) :: grid  !< Input/Output: Grid information

    ! temporary arrays needed in the (to be constructed) init procedure
    integer :: nCellsSolve
    integer, dimension(:,:), pointer :: edgesOnCell
    integer, dimension(:), pointer :: nEdgesOnCell
    integer :: i, iCell, iEdge, pointCount, maxEdgeCount
    real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell, xEdge, yEdge, zEdge
    real (kind=RKIND) :: r, cellCenter(3), alpha, tangentPlane(2,3)
    real (kind=RKIND), allocatable, dimension(:,:) :: edgeOnCellLocations, edgeOnCellNormals, coeffs, &
       edgeOnCellLocationsWork, edgeOnCellNormalsWork, coeffsWork
    real(kind=RKIND), dimension(:,:), pointer :: edgeNormalVectors
    real(kind=RKIND), dimension(:,:,:), pointer :: cellTangentPlane

    real (kind=RKIND), dimension(:,:,:), pointer :: coeffs_reconstruct

    !========================================================
    ! arrays filled and saved during init procedure
    !========================================================
    coeffs_reconstruct => grid % coeffs_reconstruct % array

    !========================================================
    ! temporary variables needed for init procedure
    !========================================================
    xCell       => grid % xCell % array
    yCell       => grid % yCell % array
    zCell       => grid % zCell % array
    xEdge       => grid % xEdge % array
    yEdge       => grid % yEdge % array
    zEdge       => grid % zEdge % array
    edgesOnCell => grid % edgesOnCell % array
    nEdgesOnCell=> grid % nEdgesOnCell % array
    nCellsSolve = grid % nCellsSolve
    edgeNormalVectors => grid % edgeNormalVectors % array
    cellTangentPlane => grid % cellTangentPlane % array


    ! init arrays
    coeffs_reconstruct = 0.0

    maxEdgeCount = maxval(nEdgesOnCell)

    allocate(edgeOnCellLocations(maxEdgeCount,3))
    allocate(edgeOnCellNormals(maxEdgeCount,3))
    allocate(coeffs(maxEdgeCount,3))

    ! loop over all cells to be solved on this block
    do iCell=1,nCellsSolve
      pointCount = nEdgesOnCell(iCell)
      cellCenter(1) = xCell(iCell)
      cellCenter(2) = yCell(iCell)
      cellCenter(3) = zCell(iCell)

      do i=1,pointCount
        iEdge = edgesOnCell(i,iCell)
        edgeOnCellLocations(i,1)  = xEdge(iEdge)
        edgeOnCellLocations(i,2)  = yEdge(iEdge)
        edgeOnCellLocations(i,3)  = zEdge(iEdge)
        edgeOnCellNormals(i,:)  = edgeNormalVectors(:, iEdge)
      end do

      alpha = 0.0
      do i=1,pointCount
        r = sqrt(sum((cellCenter - edgeOnCellLocations(i,:))**2))
        alpha = alpha + r
      enddo
      alpha = alpha/pointCount

      tangentPlane(1,:) = cellTangentPlane(:,1,iCell)
      tangentPlane(2,:) = cellTangentPlane(:,2,iCell)

      allocate(edgeOnCellLocationsWork(pointCount,3))
      allocate(edgeOnCellNormalsWork(pointCount,3))
      allocate(coeffsWork(pointCount,3))

      edgeOnCellLocationsWork = edgeOnCellLocations(1:pointCount,:)
      edgeOnCellNormalsWork = edgeOnCellNormals(1:pointCount,:)

      call mpas_rbf_interp_func_3D_plane_vec_const_dir_comp_coeffs(pointCount, &
        edgeOnCellLocationsWork, edgeOnCellNormalsWork, &
        cellCenter, alpha, tangentPlane, coeffsWork)

      coeffs(1:pointCount,:) = coeffsWork

      deallocate(edgeOnCellLocationsWork)
      deallocate(edgeOnCellNormalsWork)
      deallocate(coeffsWork)

      
      do i=1,pointCount
        coeffs_reconstruct(:,i,iCell) = coeffs(i,:)
      end do

    enddo   ! iCell

    deallocate(edgeOnCellLocations)
    deallocate(edgeOnCellNormals)
    deallocate(coeffs)

  end subroutine mpas_init_reconstruct!}}}

!***********************************************************************
!
!  routine mpas_reconstruct
!
!> \brief   MPAS Vector reconstruction routine
!> \author  Xylar Asay-Davis, Todd Ringler
!> \date    03/28/13
!> \details 
!>  Purpose: reconstruct vector field at cell centers based on radial basis functions
!>  Input: grid meta data and vector component data residing at cell edges
!>  Output: reconstructed vector field (measured in X,Y,Z) located at cell centers
!-----------------------------------------------------------------------
  subroutine mpas_reconstruct(grid, u, uReconstructX, uReconstructY, uReconstructZ, uReconstructZonal, uReconstructMeridional)!{{{

    implicit none

    type (mesh_type), intent(in) :: grid !< Input: Grid information
    real (kind=RKIND), dimension(:,:), intent(in) :: u !< Input: Velocity field on edges
    real (kind=RKIND), dimension(:,:), intent(out) :: uReconstructX !< Output: X Component of velocity reconstructed to cell centers
    real (kind=RKIND), dimension(:,:), intent(out) :: uReconstructY !< Output: Y Component of velocity reconstructed to cell centers
    real (kind=RKIND), dimension(:,:), intent(out) :: uReconstructZ !< Output: Z Component of velocity reconstructed to cell centers
    real (kind=RKIND), dimension(:,:), intent(out) :: uReconstructZonal !< Output: Zonal Component of velocity reconstructed to cell centers
    real (kind=RKIND), dimension(:,:), intent(out) :: uReconstructMeridional !< Output: Meridional Component of velocity reconstructed to cell centers

    !   temporary arrays needed in the compute procedure
    integer :: nCellsSolve
    integer, dimension(:,:), pointer :: edgesOnCell
    integer, dimension(:), pointer :: nEdgesOnCell
    integer :: iCell,iEdge, i
    real(kind=RKIND), dimension(:), pointer :: latCell, lonCell

    real (kind=RKIND), dimension(:,:,:), pointer :: coeffs_reconstruct

    logical :: on_a_sphere

    real (kind=RKIND) :: clat, slat, clon, slon


    ! stored arrays used during compute procedure
    coeffs_reconstruct => grid % coeffs_reconstruct % array

    ! temporary variables
    edgesOnCell => grid % edgesOnCell % array
    nEdgesOnCell=> grid % nEdgesOnCell % array
    nCellsSolve = grid % nCellsSolve

    latCell       => grid % latCell % array
    lonCell       => grid % lonCell % array
    on_a_sphere = grid % on_a_sphere

    ! init the intent(out)
    uReconstructX = 0.0
    uReconstructY = 0.0
    uReconstructZ = 0.0

    ! loop over cell centers
    do iCell=1,nCellsSolve
      ! a more efficient reconstruction where rbf_values*matrix_reconstruct has been precomputed
      ! in coeffs_reconstruct
      do i=1,nEdgesOnCell(iCell)
        iEdge = edgesOnCell(i,iCell)
        uReconstructX(:,iCell) = uReconstructX(:,iCell) &
          + coeffs_reconstruct(1,i,iCell) * u(:,iEdge)
        uReconstructY(:,iCell) = uReconstructY(:,iCell) &
          + coeffs_reconstruct(2,i,iCell) * u(:,iEdge)
        uReconstructZ(:,iCell) = uReconstructZ(:,iCell) &
          + coeffs_reconstruct(3,i,iCell) * u(:,iEdge)

      enddo
    enddo   ! iCell

    if(on_a_sphere) then
      do iCell=1,nCellsSolve
        clat = cos(latCell(iCell))
        slat = sin(latCell(iCell))
        clon = cos(lonCell(iCell))
        slon = sin(lonCell(iCell))
        uReconstructZonal(:,iCell) = -uReconstructX(:,iCell)*slon + uReconstructY(:,iCell)*clon
        uReconstructMeridional(:,iCell) = -(uReconstructX(:,iCell)*clon &
          + uReconstructY(:,iCell)*slon)*slat &
          + uReconstructZ(:,iCell)*clat
      end do
    else
      uReconstructZonal = uReconstructX
      uReconstructMeridional = uReconstructY
    end if

  end subroutine mpas_reconstruct!}}}

end module mpas_vector_reconstruction
