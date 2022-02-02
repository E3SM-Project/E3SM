! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_spherical_utils
!
!> \brief MPAS ocean spherical utilities
!> \author Doug Jacobsen
!> \date   03/20/2015
!> \details
!>  This module contains the routines for updating mesh quantities based on a spherical radius
!
!-----------------------------------------------------------------------

module ocn_init_spherical_utils

   use mpas_kind_types
   use mpas_io_units
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_stream_manager

   use ocn_config

   implicit none
   private

   public :: ocn_init_expand_sphere, ocn_init_realistic_coriolis_parameter
   public :: ocn_transform_from_lonlat_to_xyz
   public :: transform_from_xyz_to_lonlat, ocn_unit_vector_in_3space
   public :: ocn_vector_on_tangent_plane, ocn_cross_product_in_3space
   public :: ocn_init_set_pools_sphere_radius

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

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_expand_sphere
!
!> \brief   MPAS-Ocean Spherical Expansion Routine
!> \author  Doug Jacobsen
!> \date    03/20/2015
!> \details
!>  This routine expands mesh quantities to sphere of radius newRadius.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_expand_sphere(domain, stream_manager, newRadius, err)!{{{

   !--------------------------------------------------------------------

      type (domain_type), intent(inout) :: domain
      type (mpas_streamManager_type), intent(inout) :: stream_manager
      real (kind=RKIND), intent(in) :: newRadius
      integer, intent(out) :: err

      type (block_type), pointer :: block_ptr

      type (mpas_pool_type), pointer :: meshPool

      character (len=StrKIND) :: streamID
      integer :: directionProperty

      logical, pointer :: on_a_sphere
      real (kind=RKIND), pointer :: sphere_radius

      integer, pointer :: nCells, nCellsSolve, nEdgesSolve, nVerticesSolve, vertexDegree

      integer, dimension(:, :), pointer :: cellsOnVertex

      real (kind=RKIND), dimension(:), pointer :: areaCell, areaTriangle
      real (kind=RKIND), dimension(:), pointer :: dvEdge, dcEdge
      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell, latCell, lonCell
      real (kind=RKIND), dimension(:), pointer :: xEdge, yEdge, zEdge, latEdge, lonEdge
      real (kind=RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex, latVertex, lonVertex
      real (kind=RKIND), dimension(:, :), pointer :: kiteAreasOnVertex

      real (kind=RKIND) :: oldRadius, ratio
      real (kind=RKIND) :: norm
      real (kind=RKIND) :: oldX, oldY, oldZ
      integer :: iCell, iEdge, iVertex, i

      err = 0

      if ( .not. config_expand_sphere ) return

      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)

      call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)
      call mpas_pool_get_config(meshPool, 'sphere_radius', sphere_radius)

      if ( .not. on_a_sphere ) then
         call mpas_log_write( 'Warning: Only spherical meshes can been expanded.')
         call mpas_log_write( 'Skipping expansion')
         return
      end if

      if ( sphere_radius == 0.0_RKIND ) then
         call mpas_log_write( 'Sphere radius is 0.0', MPAS_LOG_CRIT)
         err = 1
         return
      end if

      call mpas_log_write( 'Expanding mesh to a radius of size: $f m.', realArgs=(/ newRadius /) )

      call mpas_stream_mgr_begin_iteration(stream_manager)
      do while (mpas_stream_mgr_get_next_stream(stream_manager, streamID, directionProperty))
         if ( directionProperty == MPAS_STREAM_OUTPUT .or. directionProperty == MPAS_STREAM_INPUT_OUTPUT ) then
            call mpas_stream_mgr_add_att(stream_manager, 'sphere_radius', newRadius, streamID)
         end if
      end do

      oldRadius = sphere_radius
      ratio = newRadius / oldRadius

      block_ptr => domain % blocklist
      do while(associated(block_ptr))
        call ocn_init_set_pools_sphere_radius(block_ptr % structs, newRadius)

        ! Expand cell quantities
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)

        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
        call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
        call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
        call mpas_pool_get_dimension(meshPool, 'nVerticesSolve', nVerticesSolve)
        call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree)

        call mpas_pool_get_array(meshPool, 'xCell', xCell)
        call mpas_pool_get_array(meshPool, 'yCell', yCell)
        call mpas_pool_get_array(meshPool, 'zCell', zCell)
        call mpas_pool_get_array(meshPool, 'latCell', latCell)
        call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
        call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

        call mpas_pool_get_array(meshPool, 'xEdge', xEdge)
        call mpas_pool_get_array(meshPool, 'yEdge', yEdge)
        call mpas_pool_get_array(meshPool, 'zEdge', zEdge)
        call mpas_pool_get_array(meshPool, 'latEdge', latEdge)
        call mpas_pool_get_array(meshPool, 'lonEdge', lonEdge)
        call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
        call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)

        call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
        call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
        call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
        call mpas_pool_get_array(meshPool, 'latVertex', latVertex)
        call mpas_pool_get_array(meshPool, 'lonVertex', lonVertex)
        call mpas_pool_get_array(meshPool, 'areaTriangle', areaTriangle)
        call mpas_pool_get_array(meshPool, 'kiteAreasOnVertex', kiteAreasOnVertex)
        call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)

        do iCell = 1, nCellsSolve
           oldX = xCell(iCell)
           oldY = yCell(iCell)
           oldZ = zCell(iCell)

           norm = sqrt(oldX**2 + oldY**2 + oldZ**2)

           xCell(iCell) = (oldX / norm) * newRadius
           yCell(iCell) = (oldY / norm) * newRadius
           zCell(iCell) = (oldZ / norm) * newRadius
           areaCell(iCell) = (areaCell(iCell) / (oldRadius**2) )* newRadius**2
        end do

        ! Expand vertex quantities
        do iVertex = 1, nVerticesSolve
           oldX = xVertex(iVertex)
           oldY = yVertex(iVertex)
           oldZ = zVertex(iVertex)

           norm = sqrt(oldX**2 + oldY**2 + oldZ**2)

           xVertex(iVertex) = (oldX / norm) * newRadius
           yVertex(iVertex) = (oldY / norm) * newRadius
           zVertex(iVertex) = (oldZ / norm) * newRadius
           areaTriangle(iVertex) = 0.0_RKIND

           do i = 1, vertexDegree
              if (cellsOnVertex(i, iVertex) < nCells+1) then
                 kiteAreasOnVertex(i, iVertex) = ( kiteAreasOnVertex(i, iVertex) / oldRadius**2) * newRadius**2
              else
                 kiteAreasOnVertex(i, iVertex) = 0.0_RKIND
              end if
              areaTriangle(iVertex) = areaTriangle(iVertex) + kiteAreasOnVertex(i, iVertex)
           end do
        end do

        ! Expand edge quantities
        do iEdge = 1, nEdgesSolve
           oldX = xEdge(iEdge)
           oldY = yEdge(iEdge)
           oldZ = zEdge(iEdge)

           norm = sqrt(oldX**2 + oldY**2 + oldZ**2)

           xEdge(iEdge) = (oldX / norm) * newRadius
           yEdge(iEdge) = (oldY / norm) * newRadius
           zEdge(iEdge) = (oldZ / norm) * newRadius
           dvEdge(iEdge) = (dvEdge(iEdge) / oldRadius) * newRadius
           dcEdge(iEdge) = (dcEdge(iEdge) / oldRadius) * newRadius
        end do

        block_ptr % domain % sphere_radius = newRadius
        sphere_radius = newRadius
        block_ptr => block_ptr % next

      end do

   !--------------------------------------------------------------------

   end subroutine ocn_init_expand_sphere!}}}


!***********************************************************************
!
!  routine ocn_init_realistic_coriolis_parameter
!
!> \brief   MPAS-Ocean Coriolis Initialization Routine
!> \author  Xylar Asay-Davis
!> \date    01/23/2022
!> \details
!>  This routine computes the Coriolis parameter at cells, edges and
!>  vertices if the configuration requests a realistic Coriolis
!>  parameter.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_realistic_coriolis_parameter(domain, err)!{{{

   !--------------------------------------------------------------------

      type (domain_type), intent(inout) :: domain
      integer, intent(out) :: err

      type (block_type), pointer :: block_ptr

      type (mpas_pool_type), pointer :: meshPool

      integer, pointer :: nCellsSolve, nEdgesSolve, nVerticesSolve

      real (kind=RKIND), dimension(:), pointer :: fCell, fEdge, fVertex
      real (kind=RKIND), dimension(:), pointer :: latCell, latEdge, latVertex
      integer :: iCell, iEdge, iVertex

      err = 0

      if ( .not. config_realistic_coriolis_parameter ) return

      block_ptr => domain % blocklist
      do while(associated(block_ptr))

        ! Expand cell quantities
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)

        call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
        call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
        call mpas_pool_get_dimension(meshPool, 'nVerticesSolve', nVerticesSolve)

        call mpas_pool_get_array(meshPool, 'latCell', latCell)
        call mpas_pool_get_array(meshPool, 'fCell', fCell)

        call mpas_pool_get_array(meshPool, 'latEdge', latEdge)
        call mpas_pool_get_array(meshPool, 'fEdge', fEdge)

        call mpas_pool_get_array(meshPool, 'latVertex', latVertex)
        call mpas_pool_get_array(meshPool, 'fVertex', fVertex)

        do iCell = 1, nCellsSolve
           fCell(iCell) = 2.0_RKIND * omega * sin(latCell(iCell))
        end do

        do iVertex = 1, nVerticesSolve
           fVertex(iVertex) = 2.0_RKIND * omega * sin( latVertex(iVertex) )
        end do

        do iEdge = 1, nEdgesSolve
           fEdge(iEdge) = 2.0_RKIND * omega * sin( latEdge(iEdge) )
        end do

        block_ptr => block_ptr % next

      end do

   !--------------------------------------------------------------------

   end subroutine ocn_init_realistic_coriolis_parameter!}}}

!***********************************************************************
!
!  recursive routine ocn_init_set_pools_sphere_radius
!
!> \brief   MPAS-Ocean Sphere radius update routine
!> \author  Doug Jacobsen
!> \date    09/09/2015
!> \details
!>  This routine updates the value of sphere_radius in all pools that contain
!>  it.
!
!-----------------------------------------------------------------------
   recursive subroutine ocn_init_set_pools_sphere_radius(inPool, newRadius)!{{{
      type (mpas_pool_type), intent(inout) :: inPool
      real (kind=RKIND), intent(in) :: newRadius

      type (mpas_pool_type), pointer :: subPool
      type (mpas_pool_iterator_type) :: poolItr
      real (kind=RKIND), pointer :: sphere_radius

      call mpas_pool_begin_iteration(inPool)

      do while ( mpas_pool_get_next_member(inPool, poolItr) )
         if ( poolItr % memberType == MPAS_POOL_SUBPOOL ) then
            call mpas_pool_get_subpool(inPool, poolItr % memberName, subPool)
            call ocn_init_set_pools_sphere_radius(subPool, newRadius)
         else if ( poolItr % memberType == MPAS_POOL_CONFIG ) then

            if ( poolItr % memberName == 'sphere_radius' ) then
               call mpas_pool_get_config(inPool, poolItr % memberName, sphere_radius)
               sphere_radius = newRadius
            end if

         end if
      end do

   end subroutine ocn_init_set_pools_sphere_radius!}}}

!***********************************************************************
!
!  routine ocn_transform_from_lonlat_to_xyz
!
!> \brief   MPAS-Ocean Tranform LatLon to XYZ
!> \author  Todd Ringler
!> \date    02/19/2014
!> \details
!>  This routine converts a (lat, lon) coordinate into an (x, y, z) coordinate
!>    INTENT(IN)
!>    xin = x position
!>    yin = y position
!>    zin = z position
!>    ulon = east component of vector
!>    ulat = north component of vector
!>
!>    INTENT(OUT)
!>    ux = x component of vector
!>    uy = y component of vector
!>    uz = z component of vector
!
!-----------------------------------------------------------------------
    subroutine ocn_transform_from_lonlat_to_xyz(xin, yin, zin, ulon, ulat, ux, uy, uz)!{{{
        implicit none
        real, intent(in) :: xin, yin, zin, ulon, ulat
        real, intent(out) :: ux, uy, uz
        real :: h(3,3), p(3), q(3), g(3), X1(3,3), X2(3,3), trans_X2_to_X1(3,3), r
        integer :: i,j,k
        logical :: l_Pole
        real, parameter :: epsvt = 1.0e-10_RKIND

        !-----------------------------------------------------------------------
        ! define the e1, e2, and e3 directions
        !-----------------------------------------------------------------------
        X1(1,1) = 1.0_RKIND; X1(1,2) = 0.0_RKIND; X1(1,3) = 0.0_RKIND
        X1(2,1) = 0.0_RKIND; X1(2,2) = 1.0_RKIND; X1(2,3) = 0.0_RKIND
        X1(3,1) = 0.0_RKIND; X1(3,2) = 0.0_RKIND; X1(3,3) = 1.0_RKIND

        !-----------------------------------------------------------------------
        ! find the vectors (measured in X1) that point in the local
        !   east (h(1,:)), north (h(2,:)), and vertical (h(3,:)) direction
        !-----------------------------------------------------------------------
        h(3,1) = xin; h(3,2) = yin; h(3,3) = zin
        call ocn_unit_vector_in_3space(h(3,:))

        !-----------------------------------------------------------------------
        ! g(:) is a work array and holds the vector pointing to the North Pole.
        ! measured in X1
        !-----------------------------------------------------------------------
        g(:) = X1(3,:)

        !-----------------------------------------------------------------------
        ! determine if the local vertical hits a pole
        !-----------------------------------------------------------------------
        l_Pole = .false.
        r = g(1)*h(3,1) + g(2)*h(3,2) + g(3)*h(3,3)
        r = abs(r) + epsvt
        if(r.gt.1.0_RKIND) then
           l_Pole = .true.
           h(3,:) = h(3,:) + epsvt
           call ocn_unit_vector_in_3space(h(3,:))
        endif

        !-----------------------------------------------------------------------
        ! find the vector that is perpendicular to the local vertical vector
        ! and points in the direction of of the North pole, this defines the local
        ! north direction. measured in X1
        !-----------------------------------------------------------------------
        call ocn_vector_on_tangent_plane ( h(3,:), g(:), h(2,:) )

        !-----------------------------------------------------------------------
        ! take the cross product of the local North direction and the local vertical
        ! to find the local east vector. still in X1
        !-----------------------------------------------------------------------
        call ocn_cross_product_in_3space ( h(2,:), h(3,:), h(1,:) )

        !-----------------------------------------------------------------------
        ! put these 3 vectors into a matrix X2
        !-----------------------------------------------------------------------
        X2(1,:) = h(1,:)              ! local east     (measured in X1)
        X2(2,:) = h(2,:)              ! local north    (measured in X1)
        X2(3,:) = h(3,:)              ! local vertical (measured in X1)

        !-----------------------------------------------------------------------
        ! compute the transformation matrix
        !-----------------------------------------------------------------------
        trans_X2_to_X1(:,:) = matmul(X1,transpose(X2))

        !-----------------------------------------------------------------------
        ! transform (ulon, ulat) into (x,y,z)
        !-----------------------------------------------------------------------
        p(1) = ulon; p(2) = ulat; p(3) = 0
        g(:) = matmul(trans_X2_to_X1(:, :), p(:))
        ux = g(1); uy = g(2); uz = g(3)

    end subroutine ocn_transform_from_lonlat_to_xyz!}}}

!***********************************************************************
!
!  routine ocn_transform_from_xyz_to_lonlat
!
!> \brief   MPAS-Ocean transform XYZ to LatLon
!> \author  Todd Ringler
!> \date    02/19/2014
!> \details
!>  This routine converts an (x, y, z) coordinate into a (lat, lon) coordinate
!>    INTENT(IN)
!>    xin = x position
!>    yin = y position
!>    zin = z position
!>    ux = x component of vector
!>    uy = y component of vector
!>    uz = z component of vector
!>
!>    INTENT(OUT)
!>    ulon = east component of vector
!>    ulat = north component of vector
!
!-----------------------------------------------------------------------
    subroutine transform_from_xyz_to_lonlat(xin, yin, zin, ux, uy, uz, ulon, ulat)!{{{
        implicit none
        real, intent(in) :: xin, yin, zin, ux, uy, uz
        real, intent(out) :: ulon, ulat
        real :: h(3,3), p(3), q(3), g(3), X1(3,3), X2(3,3), trans_X1_to_X2(3,3), r
        integer :: i,j,k
        logical :: l_Pole
        real, parameter :: epsvt = 1.0e-10_RKIND

        !-----------------------------------------------------------------------
        ! define the e1, e2, and e3 directions
        !-----------------------------------------------------------------------
        X1(1,1) = 1.0_RKIND; X1(1,2) = 0.0_RKIND; X1(1,3) = 0.0_RKIND
        X1(2,1) = 0.0_RKIND; X1(2,2) = 1.0_RKIND; X1(2,3) = 0.0_RKIND
        X1(3,1) = 0.0_RKIND; X1(3,2) = 0.0_RKIND; X1(3,3) = 1.0_RKIND

        !-----------------------------------------------------------------------
        ! find the vectors (measured in X1) that point in the local
        !   east (h(1,:)), north (h(2,:)), and vertical (h(3,:)) direction
        !-----------------------------------------------------------------------
        h(3,1) = xin; h(3,2) = yin; h(3,3) = zin
        call ocn_unit_vector_in_3space(h(3,:))

        !-----------------------------------------------------------------------
        ! g(:) is a work array and holds the vector pointing to the North Pole.
        ! measured in X1
        !-----------------------------------------------------------------------
        g(:) = X1(3,:)

        !-----------------------------------------------------------------------
        ! determine if the local vertical hits a pole
        !-----------------------------------------------------------------------
        l_Pole = .false.
        r = g(1)*h(3,1) + g(2)*h(3,2) + g(3)*h(3,3)
        r = abs(r) + epsvt
        if(r.gt.1.0_RKIND) then
           l_Pole = .true.
           h(3,:) = h(3,:) + epsvt
           call ocn_unit_vector_in_3space(h(3,:))
        endif

        !-----------------------------------------------------------------------
        ! find the vector that is perpendicular to the local vertical vector
        ! and points in the direction of of the North pole, this defines the local
        ! north direction. measured in X1
        !-----------------------------------------------------------------------
        call ocn_vector_on_tangent_plane ( h(3,:), g(:), h(2,:) )

        !-----------------------------------------------------------------------
        ! take the cross product of the local North direction and the local vertical
        ! to find the local east vector. still in X1
        !-----------------------------------------------------------------------
        call ocn_cross_product_in_3space ( h(2,:), h(3,:), h(1,:) )

        !-----------------------------------------------------------------------
        ! put these 3 vectors into a matrix X2
        !-----------------------------------------------------------------------
        X2(1,:) = h(1,:)              ! local east     (measured in X1)
        X2(2,:) = h(2,:)              ! local north    (measured in X1)
        X2(3,:) = h(3,:)              ! local vertical (measured in X1)

        !-----------------------------------------------------------------------
        ! compute the transformation matrix
        !-----------------------------------------------------------------------
        trans_X1_to_X2(:,:) = matmul(X2,transpose(X1))

        !-----------------------------------------------------------------------
        ! transform (ulon, ulat) into (x,y,z)
        !-----------------------------------------------------------------------
        p(1) = ux; p(2) = uy; p(3) = uz
        g(:) = matmul(trans_X1_to_X2(:, :), p(:))
        ulon = g(1); ulat= g(2);

    end subroutine transform_from_xyz_to_lonlat!}}}

!***********************************************************************
!
!  routine ocn_unit_vector_in_3space
!
!> \brief   MPAS-Ocean 3D unit vector
!> \author  Todd Ringler
!> \date    02/19/2014
!> \details
!>  This routine normalizes a vector in 3space.
!
!-----------------------------------------------------------------------
    subroutine ocn_unit_vector_in_3space (p_1)!{{{

        !-----------------------------------------------------------------------
        ! PURPOSE : normalize p_1 to unit length and overwrite p_1
        !-----------------------------------------------------------------------

        !-----------------------------------------------------------------------
        ! intent(inout)
        !-----------------------------------------------------------------------
        real , intent(inout) ::                         &
                        p_1 (:)

        !-----------------------------------------------------------------------
        ! local
        !-----------------------------------------------------------------------
        real  :: length

        length = SQRT (p_1(1)**2 + p_1(2)**2 + p_1(3)**2 )
        length = 1.0_RKIND/length
        p_1(1) = p_1(1)*length
        p_1(2) = p_1(2)*length
        p_1(3) = p_1(3)*length

    end subroutine ocn_unit_vector_in_3space!}}}

!***********************************************************************
!
!  routine ocn_vector_on_tangent_plane
!
!> \brief   MPAS-Ocean Vector on a tangent plane
!> \author  Todd Ringler
!> \date    02/19/2014
!> \details
!> Given two points measured in (x,y,z) and lying on
!> the unit sphere, find the vector (p_out) that lies on the plane
!> perpendicular to the p_1 vector and points in the direction of
!> the projection of p_2 onto the tangent plane.
!>
!>  NOTE : p_1 and p_2 are assumed to be of unit length
!>  NOTE : p_out is normalized to unit length
!
!-----------------------------------------------------------------------
    subroutine ocn_vector_on_tangent_plane(p_1, p_2, p_out)!{{{
!-----------------------------------------------------------------------
! intent(in)
!-----------------------------------------------------------------------
        real , intent(in) ::                            &
                        p_1 (:),                                      &
                        p_2 (:)

!-----------------------------------------------------------------------
! intent(out)
!-----------------------------------------------------------------------
        real , intent(out) ::                           &
                        p_out (:)

!-----------------------------------------------------------------------
! local
!-----------------------------------------------------------------------
        real  ::                                        &
                        work (3), t1(3), t2(3)

!       work (1) = - p_1(2) * ( -p_1(2) * p_2(1) + p_1(1) * p_2(2) )   &
!                  + p_1(3) * (  p_1(3) * p_2(1) - p_1(1) * p_2(3) )

!       work (2) = + p_1(1) * ( -p_1(2) * p_2(1) + p_1(1) * p_2(2) )   &
!                  - p_1(3) * ( -p_1(3) * p_2(2) + p_1(2) * p_2(3) )

!       work (3) = - p_1(1) * (  p_1(3) * p_2(1) - p_1(1) * p_2(3) )   &
!                  + p_1(2) * ( -p_1(3) * p_2(2) + p_1(2) * p_2(3) )


        t1(:) = p_2(:) - p_1(:)
        t2(:) = p_1

        call ocn_unit_vector_in_3space (t1)
        call ocn_unit_vector_in_3space (t2)

        call ocn_cross_product_in_3space(t1(:), t2(:), work(:))
        call ocn_unit_vector_in_3space (work)
        call ocn_cross_product_in_3space(t2(:),work(:),p_out(:))
        call ocn_unit_vector_in_3space (p_out)

    end subroutine ocn_vector_on_tangent_plane!}}}

!***********************************************************************
!
!  routine ocn_cross_product_in_3space
!
!> \brief   MPAS-Ocean Cross product in 3D
!> \author  Todd Ringler
!> \date    02/19/2014
!> \details
!>  compute p_1 cross p_2 and place in p_out
!
!-----------------------------------------------------------------------
        subroutine ocn_cross_product_in_3space(p_1,p_2,p_out)!{{{
!-----------------------------------------------------------------------
! intent(in)
!-----------------------------------------------------------------------
        real , intent(in) ::                            &
                        p_1 (:),                                      &
                        p_2 (:)

!-----------------------------------------------------------------------
! intent(out)
!-----------------------------------------------------------------------
        real , intent(out) ::                           &
                        p_out (:)

        p_out(1) = p_1(2)*p_2(3)-p_1(3)*p_2(2)
        p_out(2) = p_1(3)*p_2(1)-p_1(1)*p_2(3)
        p_out(3) = p_1(1)*p_2(2)-p_1(2)*p_2(1)

        end subroutine ocn_cross_product_in_3space!}}}

!***********************************************************************

end module ocn_init_spherical_utils

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
