#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module planar_mesh_mod

 use kinds, only : real_kind, long_kind
 use control_mod, only : MAX_FILE_LEN
 use parallel_mod, only : abortmp

#ifndef HOMME_WITHOUT_PIOLIBRARY
  use netcdf ! _EXTERNAL
#endif

  implicit none

  ! ===============================
  ! Public methods for mesh_mod
  ! ===============================

  public  :: MeshPlaneTopology   ! called afer MeshOpen
  public  :: MeshPlaneTopologyCoords   ! called afer MeshOpen

#ifndef HOMME_WITHOUT_PIOLIBRARY
  public  :: PlaneMeshSetCoordinates ! called after MeshPlaneTopology
#endif

contains

  subroutine MeshPlaneTopologyCoords(GridEdge, GridVertex, coord_dim1, coord_dim2, coord_dim3, coord_dimension)

  use parallel_mod,           only : abortmp
  use dimensions_mod,         only : np,  max_elements_attached_to_node
  use coordinate_systems_mod, only : cartesian2D_t, cartesian3D_t, cube_face_number_from_cart, &
                                     cube_face_number_from_sphere, cart2cubedsphere
  use gridgraph_mod,          only : GridVertex_t
  use gridgraph_mod,          only : GridEdge_t
  use cube_mod,               only : CubeSetupEdgeIndex
  use gridgraph_mod,          only : initgridedge, num_neighbors
  use control_mod,            only : north, south, east, west, neast, seast, swest, nwest, partmethod, &
                                     z2_map_method, coord_transform_method
  use params_mod,             only : SFCURVE, SPHERE_COORDS, CUBE_COORDS, FACE_2D_LB_COORDS

  use kinds, only : iulog, real_kind
  use zoltan_mod, only : is_zoltan_partition, is_zoltan_task_mapping

  implicit none
  type (GridEdge_t),   intent(inout) :: GridEdge(:)
  type (GridVertex_t), intent(inout) :: GridVertex(:)
  real(kind=real_kind),allocatable, intent(inout) :: coord_dim1(:)
  real(kind=real_kind),allocatable, intent(inout) :: coord_dim2(:)
  real(kind=real_kind),allocatable, intent(inout) :: coord_dim3(:)
  integer, intent(inout) :: coord_dimension

  call abortmp('MeshPlaneTopologyCoords not yet implemented')

  end subroutine MeshPlaneTopologyCoords


  !======================================================================
  ! subroutine MeshCubeTopology
  !======================================================================
     subroutine MeshPlaneTopology(GridEdge, GridVertex)
      use parallel_mod,           only : abortmp
      use dimensions_mod,         only : np,  max_elements_attached_to_node
      use coordinate_systems_mod, only : cartesian3D_t, cube_face_number_from_cart, cube_face_number_from_sphere
      use gridgraph_mod,          only : GridVertex_t
      use gridgraph_mod,          only : GridEdge_t
      use cube_mod,               only : CubeSetupEdgeIndex
      use gridgraph_mod,          only : initgridedge, num_neighbors
      use control_mod,            only : north, south, east, west, neast, seast, swest, nwest

      implicit none
      type (GridEdge_t),   intent(inout) :: GridEdge(:)
      type (GridVertex_t), intent(inout) :: GridVertex(:)

      call abortmp('MeshPlaneTopology not yet implemented')

      end subroutine MeshPlaneTopology


      subroutine PlaneMeshSetCoordinates(elem)
        use element_mod, only : element_t
        use parallel_mod, only : abortmp
        use coordinate_systems_mod, only   : cartesian3D_t, cartesian2d_t, spherical_polar_t, &
                                             change_coordinates, sphere2cubedsphere
        implicit none

        type (element_t), intent(inout)  :: elem(:)

        call abortmp('PlaneMeshSetCoordinates not yet implemented')

        end subroutine PlaneMeshSetCoordinates

end module planar_mesh_mod
