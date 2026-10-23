!--------------------------------------------------------------------------------
!
! 08/2016: O. Guba Modifying metric_atomic routine to add support for 'epsilon
! bubble' reference element map with GLL area = geometric area
!

! This produces periodic indexing in the range 1,ne_x or 1,ne_y
#define PERX(ii) (MODULO(ii-1,ne_x) + 1)
#define PERY(ii) (MODULO(ii-1,ne_y) + 1)

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module planar_mod
  use kinds, only : real_kind, long_kind, longdouble_kind
  use parallel_mod, only : abortmp
  use dimensions_mod, only : np,ne_x,ne_y
  use coordinate_systems_mod, only : cartesian3D_t, cartesian2d_t
  use metric_mod, only : metric_atomic
  use spacecurve_mod, only : GilbertCurve

  use physical_constants, only: Lx, Ly, Sx, Sy, dx, dy, dx_ref, dy_ref

  implicit none
  private

! We work with reference domain [0,1]^2 and physical domain [Sx,Lx+Sx] x [Sy,Ly+Sy]

  ! ==========================================
  ! Public Interfaces
  ! ==========================================

  public :: PlaneTopology

  ! ===============================
  ! Public methods for cube
  ! ===============================

  public  :: plane_init_atomic
  public  :: plane_set_corner_coordinates

  public  :: PlaneEdgeCount
  public  :: PlaneElemCount
  public  :: convert_gbl_index_plane
  public  :: plane_Dmap

  ! ===============================
  ! Private methods
  ! ===============================
  private :: coordinates_atomic
  private :: coriolis_init_atomic

contains



! =======================================
! TOPOLOGY RELATED ROUTINES
! =======================================



  subroutine PlaneTopology(GridEdge, GridVertex)
    use gridgraph_mod, only : GridEdge_t, GridVertex_t, initgridedge, PrintGridEdge, &
         allocate_gridvertex_nbrs, deallocate_gridvertex_nbrs
    use spacecurve_mod, only :  GilbertCurve
    use control_mod, only : north, south, east, west, neast, seast, swest, nwest
    !-----------------------
    implicit none

    type (GridEdge_t),   intent(inout),target     :: GridEdge(:)
    type (GridVertex_t), intent(inout),target     :: GridVertex(:)

    type (GridVertex_t),allocatable        :: GridElem(:,:)
    integer,allocatable       :: Mesh(:,:)
    integer                   :: i,j,k,number,EdgeWgtP,CornerWgt,ierr,ielem,ll,loc

    if (ne_x==0 .or. ne_y ==0) call abortmp('Error in PlaneTopology: ne_x or ne_y is zero')

    allocate(GridElem(ne_x,ne_y),stat=ierr)
    do j = 1, ne_y
      do i = 1, ne_x
         call allocate_gridvertex_nbrs(GridElem(i,j))
      end do
    end do

    if(ierr/=0) then
       call abortmp('error in allocation of GridElem structure')
    end if


    number=1
    EdgeWgtP   = np
    CornerWgt = 1
     do j=1,ne_y
        do i=1,ne_x
           ! ====================================
           ! Number elements
           ! ====================================
           GridElem(i,j)%nbrs(:)=0
           GridElem(i,j)%nbrs_wgt(:)=0
           GridElem(i,j)%nbrs_ptr(:)=0
           GridElem(i,j)%nbrs_wgt_ghost(:)=1  ! always this value
           GridElem(i,j)%SpaceCurve=0
           GridElem(i,j)%number=number
           number=number+1
      end do
    end do

! Use a generalized Hilbert space filling curve (Gilbert curve)

    allocate(Mesh(ne_x,ne_y))
    call GilbertCurve(Mesh)

   ! -------------------------------------------
   !  Setup the space-filling curve
   ! -------------------------------------------
   do j=1,ne_y
      do i=1,ne_x
         GridElem(i,j)%SpaceCurve = Mesh(i,j)
      enddo
   enddo


   ! ==================
   ! Neighbors
   ! ==================
    do j=1,ne_y
       do i=1,ne_x
          GridElem(i,j)%nbrs(west)  = GridElem(PERX(i-1),j)%number
          GridElem(i,j)%nbrs_face(west)  = 1
          GridElem(i,j)%nbrs_wgt(west)       = EdgeWgtP

          GridElem(i,j)%nbrs(south) = GridElem(i,PERY(j-1))%number
          GridElem(i,j)%nbrs_face(south) = 1
          GridElem(i,j)%nbrs_wgt(south)      = EdgeWgtP

          GridElem(i,j)%nbrs(east)   = GridElem(PERX(i+1),j)%number
          GridElem(i,j)%nbrs_face(east)   = 1
          GridElem(i,j)%nbrs_wgt(east)        = EdgeWgtP

          GridElem(i,j)%nbrs(north)  = GridElem(i,PERY(j+1))%number
          GridElem(i,j)%nbrs_face(north)  = 1
          GridElem(i,j)%nbrs_wgt(north)       = EdgeWgtP

          GridElem(i,j)%nbrs(swest) = GridElem(PERX(i-1),PERY(j-1))%number
          GridElem(i,j)%nbrs_face(swest) = 1
          GridElem(i,j)%nbrs_wgt(swest)      = CornerWgt

          GridElem(i,j)%nbrs(neast) = GridElem(PERX(i+1),PERY(j+1))%number
          GridElem(i,j)%nbrs_face(neast)  = 1
          GridElem(i,j)%nbrs_wgt(neast)       = CornerWgt

          GridElem(i,j)%nbrs(seast)  = GridElem(PERX(i+1),PERY(j-1))%number
          GridElem(i,j)%nbrs_face(seast)  = 1
          GridElem(i,j)%nbrs_wgt(seast)       = CornerWgt

          GridElem(i,j)%nbrs(nwest)  = GridElem(PERX(i-1),PERY(j+1))%number
          GridElem(i,j)%nbrs_face(nwest)  = 1
          GridElem(i,j)%nbrs_wgt(nwest)       = CornerWgt

       end do
    end do


ielem = 1                       ! Element counter
 do j=1,ne_y
    do i=1,ne_x
       GridVertex(ielem)%nbrs_ptr(1) = 1
       do ll=1,8
          loc =  GridVertex(ielem)%nbrs_ptr(ll)
          GridVertex(ielem)%nbrs(loc)       = GridElem(i,j)%nbrs(ll)
          GridVertex(ielem)%nbrs_face(loc)  = GridElem(i,j)%nbrs_face(ll)
          GridVertex(ielem)%nbrs_wgt(loc)       = GridElem(i,j)%nbrs_wgt(ll)
          GridVertex(ielem)%nbrs_wgt_ghost(loc) = GridElem(i,j)%nbrs_wgt_ghost(ll)
          GridVertex(ielem)%nbrs_ptr(ll+1) = GridVertex(ielem)%nbrs_ptr(ll)+1
       end do
       GridVertex(ielem)%number     = GridElem(i,j)%number
       GridVertex(ielem)%processor_number  = 0
       GridVertex(ielem)%SpaceCurve = GridElem(i,j)%SpaceCurve
       ielem=ielem+1
    end do
 end do

DEALLOCATE(Mesh)
   do j = 1, ne_y
      do i = 1, ne_x
         call deallocate_gridvertex_nbrs(GridElem(i,j))
      end do
  end do
DEALLOCATE(GridElem)

! =======================================
! Generate cube graph...
! =======================================
! ============================================
!  Setup the Grid edges (topology independent)
! ============================================
call initgridedge(GridEdge,GridVertex)

  end subroutine PlaneTopology





  function PlaneEdgeCount()  result(nedge)
    implicit none
    integer                     :: nedge

    if (ne_x == 0 .or. ne_y == 0) call abortmp('Error in PlaneEdgeCount: ne_x or ne_y is zero')
    nedge = ne_x*ne_y*8

  end function PlaneEdgeCount


  function PlaneElemCount()  result(nelem)
    implicit none
    integer                     :: nelem
    if (ne_x == 0 .or. ne_y == 0) call abortmp('Error in PlaneEdgeCount: ne_x or ne_y is zero')

    nelem = ne_x*ne_y
  end function PlaneElemCount

! =======================================
! GEOMETRY RELATED ROUTINES
! =======================================

  ! =======================================
  !  plane_init_atomic:
  !
  ! Initialize element descriptors for
  ! planar case for each element ...
  ! =======================================
  subroutine plane_init_atomic(elem,gll_points,alpha_in)
    use element_mod, only : element_t
    type (element_t),intent(inout) :: elem
    real (kind=real_kind),optional :: alpha_in
    real (kind=real_kind)          :: alpha=1
    real (kind=longdouble_kind)      :: gll_points(np)

    if(present(alpha_in)) alpha=alpha_in

    elem%FaceNum=0
    call coordinates_atomic(elem,gll_points)

    call metric_atomic(elem,gll_points,alpha,plane_Dmap)

    call coriolis_init_atomic(elem)


  end subroutine plane_init_atomic


  ! =======================================
  ! coordinates_atomic:
  !
  ! Initialize element coordinates for
  ! planar case ... (atomic)
  !
  ! =======================================

  subroutine coordinates_atomic(elem,gll_points)
    use element_mod, only : element_t

    type (element_t) :: elem
    real (kind=longdouble_kind)      :: gll_points(np)
    integer i,j


    ! compute the corners in Cartesian coordinates
    ! This maps [0,1] to [Sx, Lx + Sx] or [Sy, Ly + Sy]
    do i=1,4
       elem%corners3D(i)%x=Sx + elem%corners(i)%x * Lx
       elem%corners3D(i)%y=Sy + elem%corners(i)%y * Ly
       elem%corners3D(i)%z=0.0D0
    enddo

    ! =========================================
    ! compute x/y coordinates of each GLL point
    ! lat=y, lon=x, r=z
    ! =========================================

    do i=1,np
    do j=1,np
      ! this converts [-1,1] GLL points to [0, dx_ref] or [0, dy_ref], and then adds corner coords
       elem%cartp(i,j)%x= gll_points(i) * dx_ref/2.0D0 + dx_ref/2.0D0 + elem%corners(1)%x
       elem%cartp(i,j)%y= gll_points(j) * dy_ref/2.0D0 + dy_ref/2.0D0 + elem%corners(1)%y

    ! This maps [0,1] to [Sx, Lx + Sx] or [Sy, Ly + Sy]
       elem%spherep(i,j)%lon = Sx + elem%cartp(i,j)%x * Lx
       elem%spherep(i,j)%lat = Sy + elem%cartp(i,j)%y * Ly
       elem%spherep(i,j)%r= 0.0D0
    enddo
    enddo

    ! Matrix describing vector conversion to cartesian
    ! Basically just "identity" matrix
    ! x direction = zonal direction
    elem%vec_sphere2cart(:,:,1,1) = 1.0_real_kind
    elem%vec_sphere2cart(:,:,2,1) = 0.0_real_kind
    elem%vec_sphere2cart(:,:,3,1) = 0.0_real_kind
    ! y direction = meridional direction
    elem%vec_sphere2cart(:,:,1,2) = 0.0_real_kind
    elem%vec_sphere2cart(:,:,2,2) = 1.0_real_kind
    elem%vec_sphere2cart(:,:,3,2) = 0.0_real_kind
    ! z direction = vertical direction
    elem%vec_sphere2cart(:,:,1,3) = 0.0_real_kind
    elem%vec_sphere2cart(:,:,2,3) = 0.0_real_kind
    elem%vec_sphere2cart(:,:,3,3) = 1.0_real_kind

  end subroutine coordinates_atomic


  subroutine plane_Dmap(D, a,b, corners3D, ref_map, cartp, facenum)
    real (kind=real_kind), intent(out)  :: D(2,2)
    real (kind=real_kind), intent(in)     :: a,b
    type (cartesian3D_t)   :: corners3D(4)  !x,y,z coords of element corners
    integer :: ref_map
    ! only needed for ref_map=0,1
    type (cartesian2D_t),optional   :: cartp(np,np)    ! gnomonic coords of element corners
    integer,optional  :: facenum

! factor 1/2 required since HOMME reference element is [-1,1]^2 instead of [0,1]
! this is composition of 2 maps:
!  a map from [-1,1]^2 reference ELEMENT to [0,1]^2 reference DOMAIN composed of [ne_x, ne_y] elements
!  a map from [0,1]^2 reference DOMAIN to physical [Sx,Lx+Sx] x [Sy,Ly+Sy] DOMAIN
    D(1,1) = dx/2.0d0
    D(1,2) = 0.0D0
    D(2,1) = 0.0D0
    D(2,2) = dy/2.0d0

  end subroutine plane_Dmap

  subroutine coriolis_init_atomic(elem)
    use element_mod, only : element_t

    type (element_t) :: elem

    ! Local variables
    integer                  :: i,j

! default to no rotation at all
    do j=1,np
       do i=1,np
          elem%fcor(i,j)= 0.0D0
       end do
    end do

  end subroutine coriolis_init_atomic

  subroutine plane_set_corner_coordinates(elem)
    use element_mod,    only : element_t

    type (element_t) :: elem
    ! Local variables
    integer  i,ie,je,nn
    real (kind=real_kind)  :: startx, starty

    if (ne_x ==0 .or. ne_y == 0) call abortmp('Error in set_corner_coordinates: ne_x or ne_y is zero')

    ! ========================================
    ! compute planar coordinates of element
    !=========================================

    call convert_gbl_index_plane(elem%vertex%number,ie,je)

    elem%vertex%face_number = 1

    startx = ie*dx_ref
    starty = je*dy_ref

    elem%corners(1)%x = startx
    elem%corners(1)%y = starty
    elem%corners(2)%x = startx+dx_ref
    elem%corners(2)%y = starty
    elem%corners(3)%x = startx+dx_ref
    elem%corners(3)%y = starty+dy_ref
    elem%corners(4)%x = startx
    elem%corners(4)%y = starty+dy_ref

  end subroutine plane_set_corner_coordinates


  ! ================================================
  ! convert_gbl_index_plane:
  !
  ! Convert global element index to plane index
  ! ================================================

  subroutine convert_gbl_index_plane(number,ie,je)
    integer, intent(in)  :: number
    integer, intent(out) :: ie,je

    if (ne_x ==0 .or. ne_y == 0) call abortmp('Error in planar_mod:convert_gbl_index_plane: ne_x or ne_y is zero')

    !  inverse of the function:      number = 1 + ie + ne_x*je
    ie=MODULO(number-1,ne_x)
    je=(number-1)/ne_x

  end subroutine convert_gbl_index_plane

end module planar_mod
