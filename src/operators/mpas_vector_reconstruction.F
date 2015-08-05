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

  use mpas_derived_types
  use mpas_pool_routines
  use mpas_constants
  use mpas_rbf_interpolation
  use mpas_vector_operations

  implicit none

  public :: mpas_init_reconstruct, mpas_reconstruct

  interface mpas_reconstruct
     module procedure mpas_reconstruct_1d
     module procedure mpas_reconstruct_2d
  end interface

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
  subroutine mpas_init_reconstruct(meshPool, includeHalos)!{{{

    implicit none

    type (mpas_pool_type), intent(in) :: &
         meshPool         !< Input: Mesh information

    logical, optional, intent(in) :: includeHalos

    ! temporary arrays needed in the (to be constructed) init procedure
    integer, pointer :: nCells
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
    logical, pointer :: is_periodic
    real(kind=RKIND), pointer :: x_period, y_period

    logical :: includeHalosLocal

    call mpas_pool_get_config(meshPool, 'is_periodic', is_periodic)
    call mpas_pool_get_config(meshPool, 'x_period', x_period)
    call mpas_pool_get_config(meshPool, 'y_period', y_period)

    if ( present(includeHalos) ) then
       includeHalosLocal = includeHalos
    else
       includeHalosLocal = .false.
    end if

    !========================================================
    ! arrays filled and saved during init procedure
    !========================================================
    call mpas_pool_get_array(meshPool, 'coeffs_reconstruct', coeffs_reconstruct)

    !========================================================
    ! temporary variables needed for init procedure
    !========================================================
    call mpas_pool_get_array(meshPool, 'xCell', xCell)
    call mpas_pool_get_array(meshPool, 'yCell', yCell)
    call mpas_pool_get_array(meshPool, 'zCell', zCell)
    call mpas_pool_get_array(meshPool, 'xEdge', xEdge)
    call mpas_pool_get_array(meshPool, 'yEdge', yEdge)
    call mpas_pool_get_array(meshPool, 'zEdge', zEdge)
    call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
    call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
    call mpas_pool_get_array(meshPool, 'edgeNormalVectors', edgeNormalVectors)
    call mpas_pool_get_array(meshPool, 'cellTangentPlane', cellTangentPlane)

    if ( includeHalosLocal ) then
       call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
    else
       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCells)
    end if

    ! init arrays
    coeffs_reconstruct = 0.0

    maxEdgeCount = maxval(nEdgesOnCell)

    allocate(edgeOnCellLocations(maxEdgeCount,3))
    allocate(edgeOnCellNormals(maxEdgeCount,3))
    allocate(coeffs(maxEdgeCount,3))

    ! loop over all cells to be solved on this block
    do iCell=1,nCells
      pointCount = nEdgesOnCell(iCell)
      cellCenter(1) = xCell(iCell)
      cellCenter(2) = yCell(iCell)
      cellCenter(3) = zCell(iCell)

      do i=1,pointCount
        iEdge = edgesOnCell(i,iCell)
        if (is_periodic) then
          edgeOnCellLocations(i,1)  = mpas_fix_periodicity(xEdge(iEdge), cellCenter(1), x_period)
          edgeOnCellLocations(i,2)  = mpas_fix_periodicity(yEdge(iEdge), cellCenter(2), y_period)
          edgeOnCellLocations(i,3)  = zEdge(iEdge)
        else
          edgeOnCellLocations(i,1)  = xEdge(iEdge)
          edgeOnCellLocations(i,2)  = yEdge(iEdge)
          edgeOnCellLocations(i,3)  = zEdge(iEdge)
        end if
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
!  routine mpas_reconstruct_2d
!
!> \brief   2d MPAS Vector reconstruction routine
!> \author  Xylar Asay-Davis, Todd Ringler
!> \date    03/28/13
!> \details 
!>  Purpose: reconstruct vector field at cell centers based on radial basis functions
!>  Input: grid meta data and vector component data residing at cell edges
!>  Output: reconstructed vector field (measured in X,Y,Z) located at cell centers
!-----------------------------------------------------------------------
  subroutine mpas_reconstruct_2d(meshPool, u, uReconstructX, uReconstructY, uReconstructZ, uReconstructZonal, uReconstructMeridional, includeHalos)!{{{

    implicit none

    type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information
    real (kind=RKIND), dimension(:,:), intent(in) :: u !< Input: Velocity field on edges
    real (kind=RKIND), dimension(:,:), intent(out) :: uReconstructX !< Output: X Component of velocity reconstructed to cell centers
    real (kind=RKIND), dimension(:,:), intent(out) :: uReconstructY !< Output: Y Component of velocity reconstructed to cell centers
    real (kind=RKIND), dimension(:,:), intent(out) :: uReconstructZ !< Output: Z Component of velocity reconstructed to cell centers
    real (kind=RKIND), dimension(:,:), intent(out) :: uReconstructZonal !< Output: Zonal Component of velocity reconstructed to cell centers
    real (kind=RKIND), dimension(:,:), intent(out) :: uReconstructMeridional !< Output: Meridional Component of velocity reconstructed to cell centers
    logical, optional, intent(in) :: includeHalos !< Input: Optional logical that allows reconstruction over halo regions

    !   temporary arrays needed in the compute procedure
    logical :: includeHalosLocal
    integer, pointer :: nCells
    integer, dimension(:,:), pointer :: edgesOnCell
    integer, dimension(:), pointer :: nEdgesOnCell
    integer :: iCell,iEdge, i
    real(kind=RKIND), dimension(:), pointer :: latCell, lonCell

    real (kind=RKIND), dimension(:,:,:), pointer :: coeffs_reconstruct

    logical, pointer :: on_a_sphere

    real (kind=RKIND) :: clat, slat, clon, slon

    if ( present(includeHalos) ) then
       includeHalosLocal = includeHalos
    else
       includeHalosLocal = .false.
    end if

    ! stored arrays used during compute procedure
    call mpas_pool_get_array(meshPool, 'coeffs_reconstruct', coeffs_reconstruct)

    ! temporary variables
    call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
    call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)

    if ( includeHalosLocal ) then
       call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
    else
       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCells)
    end if

    call mpas_pool_get_array(meshPool, 'latCell', latCell)
    call mpas_pool_get_array(meshPool, 'lonCell', lonCell)

    call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

    ! init the intent(out)
    uReconstructX = 0.0
    uReconstructY = 0.0
    uReconstructZ = 0.0

    ! loop over cell centers
    do iCell = 1, nCells
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

    if (on_a_sphere) then
      do iCell = 1, nCells
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

  end subroutine mpas_reconstruct_2d!}}}


!***********************************************************************
!
!  routine mpas_reconstruct_1d
!
!> \brief   1d MPAS Vector reconstruction routine
!> \author  Xylar Asay-Davis, Todd Ringler, Matt Hoffman
!> \date    03/28/13
!> \details 
!>  Purpose: reconstruct vector field at cell centers based on radial basis functions
!>  Input: grid meta data and vector component data residing at cell edges
!>  Output: reconstructed vector field (measured in X,Y,Z) located at cell centers
!-----------------------------------------------------------------------
  subroutine mpas_reconstruct_1d(meshPool, u, uReconstructX, uReconstructY, uReconstructZ, uReconstructZonal, uReconstructMeridional, includeHalos)!{{{

    implicit none

    type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information
    real (kind=RKIND), dimension(:), intent(in) :: u !< Input: Velocity field on edges
    real (kind=RKIND), dimension(:), intent(out) :: uReconstructX !< Output: X Component of velocity reconstructed to cell centers
    real (kind=RKIND), dimension(:), intent(out) :: uReconstructY !< Output: Y Component of velocity reconstructed to cell centers
    real (kind=RKIND), dimension(:), intent(out) :: uReconstructZ !< Output: Z Component of velocity reconstructed to cell centers
    real (kind=RKIND), dimension(:), intent(out) :: uReconstructZonal !< Output: Zonal Component of velocity reconstructed to cell centers
    real (kind=RKIND), dimension(:), intent(out) :: uReconstructMeridional !< Output: Meridional Component of velocity reconstructed to cell centers
    logical, optional, intent(in) :: includeHalos !< Input: Logical flag that allows reconstructing over halo regions

    !   temporary arrays needed in the compute procedure
    integer, pointer :: nCells
    integer, dimension(:,:), pointer :: edgesOnCell
    integer, dimension(:), pointer :: nEdgesOnCell
    integer :: iCell,iEdge, i
    real(kind=RKIND), dimension(:), pointer :: latCell, lonCell

    real (kind=RKIND), dimension(:,:,:), pointer :: coeffs_reconstruct

    logical, pointer :: on_a_sphere
    logical :: includeHalosLocal

    real (kind=RKIND) :: clat, slat, clon, slon

    if ( present(includeHalos) ) then
       includeHalosLocal = includeHalos
    else
       includeHalosLocal = .false.
    end if

    ! stored arrays used during compute procedure
    call mpas_pool_get_array(meshPool, 'coeffs_reconstruct', coeffs_reconstruct)

    ! temporary variables
    call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
    call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)

    if ( includeHalosLocal ) then
       call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
    else
       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCells)
    end if

    call mpas_pool_get_array(meshPool, 'latCell', latCell)
    call mpas_pool_get_array(meshPool, 'lonCell', lonCell)

    call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

    ! init the intent(out)
    uReconstructX = 0.0
    uReconstructY = 0.0
    uReconstructZ = 0.0

    ! loop over cell centers
    do iCell = 1, nCells
      ! a more efficient reconstruction where rbf_values*matrix_reconstruct has been precomputed
      ! in coeffs_reconstruct
      do i=1,nEdgesOnCell(iCell)
        iEdge = edgesOnCell(i,iCell)
        uReconstructX(iCell) = uReconstructX(iCell) &
          + coeffs_reconstruct(1,i,iCell) * u(iEdge)
        uReconstructY(iCell) = uReconstructY(iCell) &
          + coeffs_reconstruct(2,i,iCell) * u(iEdge)
        uReconstructZ(iCell) = uReconstructZ(iCell) &
          + coeffs_reconstruct(3,i,iCell) * u(iEdge)

      enddo
    enddo   ! iCell

    if (on_a_sphere) then
      do iCell = 1, nCells
        clat = cos(latCell(iCell))
        slat = sin(latCell(iCell))
        clon = cos(lonCell(iCell))
        slon = sin(lonCell(iCell))
        uReconstructZonal(iCell) = -uReconstructX(iCell)*slon + uReconstructY(iCell)*clon
        uReconstructMeridional(iCell) = -(uReconstructX(iCell)*clon &
          + uReconstructY(iCell)*slon)*slat &
          + uReconstructZ(iCell)*clat
      end do
    else
      uReconstructZonal = uReconstructX
      uReconstructMeridional = uReconstructY
    end if

  end subroutine mpas_reconstruct_1d!}}}

end module mpas_vector_reconstruction

