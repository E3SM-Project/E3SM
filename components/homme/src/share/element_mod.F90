#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module element_mod

  use kinds,                  only: real_kind, long_kind, int_kind
  use coordinate_systems_mod, only: spherical_polar_t, cartesian2D_t, cartesian3D_t, distance
  use dimensions_mod,         only: np, npsq, nlev, nlevp, qsize_d, max_neigh_edges
  use edgetype_mod,           only: edgedescriptor_t
  use gridgraph_mod,          only: gridvertex_t
  use element_state,          only: elem_state_t, derived_state_t, elem_accum_t
  implicit none
  private

  !___________________________________________________________________
  public :: element_coordinates
  public :: element_var_coordinates
  public :: element_var_coordinates3D
  public :: GetColumnIdP
  public :: allocate_element_desc
  public :: setup_element_pointers

  type, public :: index_t
     integer(kind=int_kind) :: ia(npsq),ja(npsq)
     integer(kind=int_kind) :: is,ie
     integer(kind=int_kind) :: NumUniquePts
     integer(kind=int_kind) :: UniquePtOffset
  end type index_t

  !___________________________________________________________________
  type, public :: element_t
     integer(kind=int_kind) :: LocalId
     integer(kind=int_kind) :: GlobalId

     ! Coordinate values of element points
     type (spherical_polar_t) :: spherep(np,np)                       ! Spherical coords of GLL points

     ! Equ-angular gnomonic projection coordinates
     type (cartesian2D_t)     :: cartp(np,np)                         ! gnomonic coords of GLL points
     type (cartesian2D_t)     :: corners(4)                           ! gnomonic coords of element corners

     ! 3D cartesian coordinates
     type (cartesian3D_t)     :: corners3D(4)

     ! Element diagnostics
     real (kind=real_kind)    :: area                                 ! Area of element
     real (kind=real_kind)    :: normDinv                             ! some type of norm of Dinv used for CFL
     real (kind=real_kind)    :: dx_short                             ! short length scale in km
     real (kind=real_kind)    :: dx_long                              ! long length scale in km

     real (kind=real_kind)    :: variable_hyperviscosity(np,np)       ! hyperviscosity based on above
     real (kind=real_kind)    :: hv_courant                           ! hyperviscosity courant number
     real (kind=real_kind)    :: tensorVisc(np,np,2,2)                !og, matrix V for tensor viscosity

     type (GridVertex_t)      :: vertex                               ! element grid vertex information
     type (EdgeDescriptor_t)  :: desc

     type (elem_state_t)      :: state
     type (derived_state_t)   :: derived
     type (elem_accum_t)       :: accum

     ! Metric terms
     real (kind=real_kind)    :: met(np,np,2,2)                       ! metric tensor on velocity and pressure grid
     real (kind=real_kind)    :: metinv(np,np,2,2)                    ! metric tensor on velocity and pressure grid
     real (kind=real_kind)    :: metdet(np,np)                        ! g = SQRT(det(g_ij)) on velocity and pressure grid
     real (kind=real_kind)    :: rmetdet(np,np)                       ! 1/metdet on velocity pressure grid
     real (kind=real_kind)    :: D(np,np,2,2)                         ! Map covariant field on cube to vector field on the sphere
     real (kind=real_kind)    :: Dinv(np,np,2,2)                      ! Map vector field on the sphere to covariant v on cube


     ! Mass flux across the sides of each sub-element.
     ! The storage is redundent since the mass across shared sides
     ! must be equal in magnitude and opposite in sign.
     ! The layout is like:
     !   --------------------------------------------------------------
     ! ^|    (1,4,3)     |                |              |    (4,4,3) |
     ! ||                |                |              |            |
     ! ||(1,4,4)         |                |              |(4,4,4)     |
     ! ||         (1,4,2)|                |              |     (4,4,2)|
     ! ||                |                |              |            |
     ! ||   (1,4,1)      |                |              |  (4,4,1)   |
     ! |---------------------------------------------------------------
     ! S|                |                |              |            |
     ! e|                |                |              |            |
     ! c|                |                |              |            |
     ! o|                |                |              |            |
     ! n|                |                |              |            |
     ! d|                |                |              |            |
     !  ---------------------------------------------------------------
     ! C|                |                |              |            |
     ! o|                |                |              |            |
     ! o|                |                |              |            |
     ! r|                |                |              |            |
     ! d|                |                |              |            |
     ! i|                |                |              |            |
     ! n---------------------------------------------------------------
     ! a|    (1,1,3)     |                |              |    (4,1,3) |
     ! t|                |                |              |(4,1,4)     |
     ! e|(1,1,4)         |                |              |            |
     !  |         (1,1,2)|                |              |     (4,1,2)|
     !  |                |                |              |            |
     !  |    (1,1,1)     |                |              |  (4,1,1)   |
     !  ---------------------------------------------------------------
     !          First Coordinate ------->
     ! Convert vector fields from spherical to rectangular components
     ! The transpose of this operation is its pseudoinverse.
     real (kind=real_kind)    :: vec_sphere2cart(np,np,3,2)

     ! Mass matrix terms for an element on a cube face
     real (kind=real_kind)    :: mp(np,np)                            ! mass matrix on v and p grid
     real (kind=real_kind)    :: rmp(np,np)                           ! inverse mass matrix on v and p grid

     ! Mass matrix terms for an element on the sphere
     ! This mass matrix is used when solving the equations in weak form
     ! with the natural (surface area of the sphere) inner product
     real (kind=real_kind)    :: spheremp(np,np)                      ! mass matrix on v and p grid
     real (kind=real_kind)    :: rspheremp(np,np)                     ! inverse mass matrix on v and p grid

     integer(kind=long_kind)  :: gdofP(np,np)                         ! global degree of freedom (P-grid)

     real (kind=real_kind)    :: fcor(np,np)                          ! Coreolis term

     type (index_t) :: idxP
     integer :: FaceNum

     ! force element_t to be a multiple of 8 bytes.
     ! on BGP, code will crash (signal 7, or signal 15) if 8 byte alignment is off
     ! check core file for:
     ! core.63:Generated by interrupt..(Alignment Exception DEAR=0xa1ef671c ESR=0x01800000 CCR0=0x4800a002)
     !integer :: dummy
  end type element_t


contains

  function GetColumnIdP(elem,i,j) result(col_id)

    ! Get unique identifier for a Physics column on the P-grid

    type(element_t), intent(in) :: elem
    integer, intent(in) :: i,j
    integer :: col_id
    col_id = elem%gdofP(i,j)
  end function GetColumnIdP

  !___________________________________________________________________
  function element_coordinates(start,end,points) result(cart)

    ! Initialize 2D rectilinear element colocation points

    use kinds, only : longdouble_kind
    type (cartesian2D_t), intent(in) :: start
    type (cartesian2D_t), intent(in) :: end
    real (kind=longdouble_kind), intent(in) :: points(:)

    type (cartesian2D_t) :: cart(SIZE(points),SIZE(points))
    type (cartesian2D_t) :: length, centroid
    real (kind=longdouble_kind) :: y
    integer i,j

    length%x   = 0.50D0*(end%x-start%x)
    length%y   = 0.50D0*(end%y-start%y)
    centroid%x = 0.50D0*(end%x+start%x)
    centroid%y = 0.50D0*(end%y+start%y)
    do j=1,SIZE(points)
       y = centroid%y + length%y*points(j)
       do i=1,SIZE(points)
          cart(i,j)%x = centroid%x + length%x*points(i)
          cart(i,j)%y = y
       end do
    end do
  end function element_coordinates

  !___________________________________________________________________
  function element_var_coordinates(c,points) result(cart)

    use kinds, only : longdouble_kind
    type (cartesian2D_t), intent(in) :: c(4)
    real (kind=longdouble_kind), intent(in) :: points(:)
    type (cartesian2D_t) :: cart(SIZE(points),SIZE(points))

    real (kind=longdouble_kind) :: p(size(points))
    real (kind=longdouble_kind) :: q(size(points))
    integer i,j

    p(:) = (1.0D0-points(:))/2.0D0
    q(:) = (1.0D0+points(:))/2.0D0

    do j=1,SIZE(points)
       do i=1,SIZE(points)
          cart(i,j)%x = p(i)*p(j)*c(1)%x &
                      + q(i)*p(j)*c(2)%x &
                      + q(i)*q(j)*c(3)%x &
                      + p(i)*q(j)*c(4)%x
          cart(i,j)%y = p(i)*p(j)*c(1)%y &
                      + q(i)*p(j)*c(2)%y &
                      + q(i)*q(j)*c(3)%y &
                      + p(i)*q(j)*c(4)%y
       end do
    end do
  end function element_var_coordinates

  !___________________________________________________________________
  function element_var_coordinates3d(c,points) result(cart)

    use kinds, only : longdouble_kind
    type (cartesian3D_t), intent(in) :: c(4)
    real (kind=longdouble_kind), intent(in) :: points(:)
    type (cartesian3D_t) :: cart(SIZE(points),SIZE(points))

    real (kind=longdouble_kind) :: p(size(points))
    real (kind=longdouble_kind) :: q(size(points)),r
    integer i,j

    p(:) = (1.0D0-points(:))/2.0D0
    q(:) = (1.0D0+points(:))/2.0D0

    do j=1,SIZE(points)
       do i=1,SIZE(points)
          cart(i,j)%x = p(i)*p(j)*c(1)%x &
                      + q(i)*p(j)*c(2)%x &
                      + q(i)*q(j)*c(3)%x &
                      + p(i)*q(j)*c(4)%x
          cart(i,j)%y = p(i)*p(j)*c(1)%y &
                      + q(i)*p(j)*c(2)%y &
                      + q(i)*q(j)*c(3)%y &
                      + p(i)*q(j)*c(4)%y
          cart(i,j)%z = p(i)*p(j)*c(1)%z &
                      + q(i)*p(j)*c(2)%z &
                      + q(i)*q(j)*c(3)%z &
                      + p(i)*q(j)*c(4)%z

          ! project back to sphere:
          r = distance(cart(i,j))
          cart(i,j)%x = cart(i,j)%x/r
          cart(i,j)%y = cart(i,j)%y/r
          cart(i,j)%z = cart(i,j)%z/r
       end do
    end do
  end function element_var_coordinates3d

  !___________________________________________________________________
  subroutine allocate_element_desc(elem)

    type (element_t), intent(inout)   :: elem(:)
    integer                           :: num, j,i

    num = SIZE(elem)

    do j=1,num
       allocate(elem(j)%desc%putmapP(max_neigh_edges))
       allocate(elem(j)%desc%getmapP(max_neigh_edges))
       allocate(elem(j)%desc%putmapP_ghost(max_neigh_edges))
       allocate(elem(j)%desc%getmapP_ghost(max_neigh_edges))
       allocate(elem(j)%desc%putmapS(max_neigh_edges))
       allocate(elem(j)%desc%getmapS(max_neigh_edges))
       allocate(elem(j)%desc%reverse(max_neigh_edges))
       allocate(elem(j)%desc%globalID(max_neigh_edges))
       allocate(elem(j)%desc%loc2buf(max_neigh_edges))
       do i=1,max_neigh_edges
          elem(j)%desc%loc2buf(i)=i
          elem(j)%desc%globalID(i)=-1
       enddo

    end do
  end subroutine allocate_element_desc


! this should go in openACC's element_state.F90, but it cant because that
! module doesn't know about element_t.  
  subroutine setup_element_pointers(elem)
    use dimensions_mod, only: nelemd, qsize
#if USE_OPENACC
    use element_state, only : state_Qdp, derived_vn0, derived_divdp, derived_divdp_proj
#endif
    implicit none
    type(element_t), intent(inout) :: elem(:)
    integer :: ie
#if USE_OPENACC
    allocate( state_Qdp                (np,np,nlev,qsize,2,nelemd)            )
    allocate( derived_vn0              (np,np,2,nlev,nelemd)                  )
    allocate( derived_divdp            (np,np,nlev,nelemd)                    )
    allocate( derived_divdp_proj       (np,np,nlev,nelemd)                    )
    do ie = 1 , nelemd
      elem(ie)%state%Qdp                 => state_Qdp                (:,:,:,:,:,ie)
      elem(ie)%derived%vn0               => derived_vn0              (:,:,:,:,ie)  
      elem(ie)%derived%divdp             => derived_divdp            (:,:,:,ie)    
      elem(ie)%derived%divdp_proj        => derived_divdp_proj       (:,:,:,ie)    
    enddo
#endif
  end subroutine setup_element_pointers




end module element_mod
