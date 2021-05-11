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
  use control_mod, only : hypervis_scaling, cubed_sphere_map
  use spacecurve_mod, only : GilbertCurve

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

  ! ===============================
  ! Private methods
  ! ===============================
  private :: coordinates_atomic
  private :: metric_atomic
  private :: coriolis_init_atomic
  private :: Dmap

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

    call metric_atomic(elem,gll_points,alpha)

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
    use physical_constants, only: dx, dy, dx_ref, dy_ref, Sx, Sy, Lx, Ly

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

  end subroutine coordinates_atomic


  subroutine Dmap(D, a,b, corners3D, ref_map, cartp, facenum)
    use physical_constants, only: dx, dy
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

  end subroutine Dmap

  ! =========================================
  ! metric_atomic:
  !
  ! Initialize planar metric terms:
  ! initialize:
  !         metdet, rmetdet  (analytic)    = detD, 1/detD
  !         met                (analytic)    D^t D     (symmetric)
  !         metdet             (analytic)    = detD
  !         metinv             (analytic)    Dinv Dinv^t  (symmetic)
  !         D     (from subroutine vmap)
  !         Dinv  (computed directly from D)
  !
  ! ucontra = Dinv * u  =  metinv * ucov
  ! ucov    = D^t * u   =  met * ucontra
  !
  ! we also compute DE = D*E, where
  ! E = eigenvectors of metinv as a basis      metinv = E LAMBDA E^t
  !
  ! ueig = E^t ucov  = E^t D^t u =  (DE)^t u
  !
  !
  ! so if we want to tweak the mapping by a factor alpha (so the weights add up to domain size, for example)
  ! we take:
  !    NEW       OLD
  !       D = sqrt(alpha) D  and then rederive all quantities.
  !    detD = alpha detD
  !
  ! where alpha = domain size/SEMarea, SEMarea = global sum elem(ie)%mv(i,j)*elem(ie)%metdet(i,j)
  !
  ! =========================================

  subroutine metric_atomic(elem,gll_points,alpha)
    use element_mod, only : element_t
    use physical_constants, only: Lx, Ly, dx, dy

    type (element_t) :: elem
    real(kind=real_kind) :: alpha
    real (kind=longdouble_kind)      :: gll_points(np)
    ! Local variables
    integer ii
    integer i,j,nn
    integer iptr

    real (kind=real_kind) :: r         ! distance from origin for point on cube tangent to unit sphere

    real (kind=real_kind) :: const, norm
    real (kind=real_kind) :: detD      ! determinant of vector field mapping matrix.
    real (kind=real_kind) :: tmpD(2,2)
    real (kind=real_kind) :: M(2,2),E(2,2),eig(2),DE(2,2),DEL(2,2),V(2,2), nu1, nu2, lamStar1, lamStar2
    integer :: imaxM(2)
    real (kind=real_kind) :: l1, l2, sc,min_svd,max_svd,max_normDinv

    ! ==============================================
    ! Initialize differential mapping operator
    ! to and from vector fields on the physical domain to
    ! contravariant vector fields on the reference domain
    ! i.e. dM/dx^i in Sadourney (1972) and it's
    ! inverse
    ! ==============================================

    max_svd = 0.0d0
    max_normDinv = 0.0d0
    min_svd = 1d99
    do j=1,np
       do i=1,np

          call Dmap(elem%D(i,j,:,:),1.0D0,1.0D0,elem%corners3D,cubed_sphere_map)

          ! Numerical metric tensor based on analytic D: met = D^T times D
          ! (D maps between physical plane and reference element)
          elem%met(i,j,1,1) = elem%D(i,j,1,1)*elem%D(i,j,1,1) + &
                              elem%D(i,j,2,1)*elem%D(i,j,2,1)
          elem%met(i,j,1,2) = elem%D(i,j,1,1)*elem%D(i,j,1,2) + &
                              elem%D(i,j,2,1)*elem%D(i,j,2,2)
          elem%met(i,j,2,1) = elem%D(i,j,1,1)*elem%D(i,j,1,2) + &
                              elem%D(i,j,2,1)*elem%D(i,j,2,2)
          elem%met(i,j,2,2) = elem%D(i,j,1,2)*elem%D(i,j,1,2) + &
                              elem%D(i,j,2,2)*elem%D(i,j,2,2)

          ! compute D^-1...
          ! compute determinant of D mapping matrix... if not zero compute inverse

          detD = elem%D(i,j,1,1)*elem%D(i,j,2,2) - elem%D(i,j,1,2)*elem%D(i,j,2,1)

          elem%Dinv(i,j,1,1) =  elem%D(i,j,2,2)/detD
          elem%Dinv(i,j,1,2) = -elem%D(i,j,1,2)/detD
          elem%Dinv(i,j,2,1) = -elem%D(i,j,2,1)/detD
          elem%Dinv(i,j,2,2) =  elem%D(i,j,1,1)/detD

          ! L2 norm = sqrt max eigenvalue of metinv
          !         = 1/sqrt(min eigenvalue of met)
          ! l1 and l2 are eigenvalues of met
          ! (should both be positive, l1 > l2)
          l1 = (elem%met(i,j,1,1) + elem%met(i,j,2,2) + sqrt(4.0d0*elem%met(i,j,1,2)*elem%met(i,j,2,1) + &
              (elem%met(i,j,1,1) - elem%met(i,j,2,2))**2))/2.0d0
          l2 = (elem%met(i,j,1,1) + elem%met(i,j,2,2) - sqrt(4.0d0*elem%met(i,j,1,2)*elem%met(i,j,2,1) + &
              (elem%met(i,j,1,1) - elem%met(i,j,2,2))**2))/2.0d0
          ! Max L2 norm of Dinv is sqrt of max eigenvalue of metinv
          ! max eigenvalue of metinv is 1/min eigenvalue of met
          norm = 1.0d0/sqrt(min(abs(l1),abs(l2)))
          max_svd = max(norm, max_svd)
          ! Min L2 norm of Dinv is sqrt of min eigenvalue of metinv
          ! min eigenvalue of metinv is 1/max eigenvalue of met
          norm = 1.0d0/sqrt(max(abs(l1),abs(l2)))
          min_svd = min(norm, min_svd)

          ! some kind of pseudo-norm of Dinv
          ! C = 1/sqrt(2) sqrt( |g^x|^2 + |g^y|^2 + 2*|g^x dot g^y|)
          !   = 1/sqrt(2) sqrt( |g_x|^2 + |g_y|^2 + 2*|g_x dot g_y|) / J
          ! g^x = Dinv(:,1)    g_x = D(1,:)
          ! g^y = Dinv(:,2)    g_y = D(2,:)
          norm = (2*abs(sum(elem%Dinv(i,j,:,1)*elem%Dinv(i,j,:,2))) + sum(elem%Dinv(i,j,:,1)**2) + sum(elem%Dinv(i,j,:,2)**2))
          norm = sqrt(norm)
!          norm = (2*abs(sum(elem%D(1,:,i,j)*elem%D(2,:,i,j))) + sum(elem%D(1,:,i,j)**2) + sum(elem%D(2,:,i,j)**2))
!          norm = sqrt(norm)/detD
          max_normDinv = max(norm,max_normDinv)


          ! Need inverse of met if not calculated analytically
          elem%metdet(i,j) = abs(detD)
          elem%rmetdet(i,j) = 1.0D0/abs(detD)

          elem%metinv(i,j,1,1) =  elem%met(i,j,2,2)/(detD*detD)
          elem%metinv(i,j,1,2) = -elem%met(i,j,1,2)/(detD*detD)
          elem%metinv(i,j,2,1) = -elem%met(i,j,2,1)/(detD*detD)
          elem%metinv(i,j,2,2) =  elem%met(i,j,1,1)/(detD*detD)

          ! matricies for tensor hyper-viscosity
          ! compute eigenvectors of metinv (probably same as computed above)
          M = elem%metinv(i,j,:,:)

          eig(1) = (M(1,1) + M(2,2) + sqrt(4.0d0*M(1,2)*M(2,1) + &
              (M(1,1) - M(2,2))**2))/2.0d0
          eig(2) = (M(1,1) + M(2,2) - sqrt(4.0d0*M(1,2)*M(2,1) + &
              (M(1,1) - M(2,2))**2))/2.0d0

          ! use DE to store M - Lambda, to compute eigenvectors
          DE=M
          DE(1,1)=DE(1,1)-eig(1)
          DE(2,2)=DE(2,2)-eig(1)

          imaxM = maxloc(abs(DE))
          if (maxval(abs(DE))==0) then
             E(1,1)=1; E(2,1)=0;
          elseif ( imaxM(1)==1 .and. imaxM(2)==1 ) then
             E(2,1)=1; E(1,1) = -DE(2,1)/DE(1,1)
          else   if ( imaxM(1)==1 .and. imaxM(2)==2 ) then
             E(2,1)=1; E(1,1) = -DE(2,2)/DE(1,2)
          else   if ( imaxM(1)==2 .and. imaxM(2)==1 ) then
             E(1,1)=1; E(2,1) = -DE(1,1)/DE(2,1)
          else   if ( imaxM(1)==2 .and. imaxM(2)==2 ) then
             E(1,1)=1; E(2,1) = -DE(1,2)/DE(2,2)
          else
             call abortmp('Impossible error in planar_mod.F90::metric_atomic()')
          endif

          ! the other eigenvector is orthgonal:
          E(1,2)=-E(2,1)
          E(2,2)= E(1,1)

!normalize columns
	  E(:,1)=E(:,1)/sqrt(sum(E(:,1)*E(:,1)));
	  E(:,2)=E(:,2)/sqrt(sum(E(:,2)*E(:,2)));


! OBTAINING TENSOR FOR HV: follows same approach as in cube_mod, with spherical-specific scaling removed

!matrix D*E
          DE(1,1)=sum(elem%D(i,j,1,:)*E(:,1))
          DE(1,2)=sum(elem%D(i,j,1,:)*E(:,2))
          DE(2,1)=sum(elem%D(i,j,2,:)*E(:,1))
          DE(2,2)=sum(elem%D(i,j,2,:)*E(:,2))

	  lamStar1=1/(eig(1)**(hypervis_scaling/4.0d0))
	  lamStar2=1/(eig(2)**(hypervis_scaling/4.0d0))

!matrix (DE) * Lam^* * Lam , tensor HV when V is applied at each Laplace calculation
!          DEL(1:2,1) = lamStar1*eig(1)*DE(1:2,1)
!          DEL(1:2,2) = lamStar2*eig(2)*DE(1:2,2)

!matrix (DE) * (Lam^*)^2 * Lam, tensor HV when V is applied only once, at the last Laplace calculation
!will only work with hyperviscosity, not viscosity
          DEL(1:2,1) = (lamStar1**2) *eig(1)*DE(1:2,1)
          DEL(1:2,2) = (lamStar2**2) *eig(2)*DE(1:2,2)

!matrix (DE) * Lam^* * Lam  *E^t *D^t or (DE) * (Lam^*)^2 * Lam  *E^t *D^t
          V(1,1)=sum(DEL(1,:)*DE(1,:))
          V(1,2)=sum(DEL(1,:)*DE(2,:))
          V(2,1)=sum(DEL(2,:)*DE(1,:))
          V(2,2)=sum(DEL(2,:)*DE(2,:))

	  elem%tensorVisc(i,j,:,:)=V(:,:)

       end do
    end do

!    see Paul Ullrich writeup:
!    max_normDinv might be a tighter bound than max_svd for deformed elements
!    max_svd >= max_normDinv/sqrt(2), with equality holding if |g^x| = |g^y|
!    elem%normDinv=max_normDinv/sqrt(2)

    ! this norm is consistent with length scales defined below:
    elem%normDinv=max_svd


    ! compute element length scales, based on SVDs, in km:
    elem%dx_short = 1.0d0/(max_svd*0.5d0*dble(np-1)*1000.0d0)
    elem%dx_long  = 1.0d0/(min_svd*0.5d0*dble(np-1)*1000.0d0)

    ! Area correction: Bring numerical area from integration weights to
    ! agreement with geometric area.
    ! Three different cases:
    ! (1) alpha == 1, this means that cube_init_atomic wasn't
    ! called with alpha parameter there will be no correction,
    ! (2) alpha <> 1 and cubed_sphere_map=0 and correction is so-called
    ! 'alpha-correction',
    ! (3) alpha <> 1 and cubed_sphere_map=2 and it is 'epsilon-bubble'
    ! correction.

    if( cubed_sphere_map == 0 ) then
       ! alpha correction for cases (1) and (2).
       elem%D = elem%D * sqrt(alpha)
       elem%Dinv = elem%Dinv / sqrt(alpha)
       elem%metdet = elem%metdet * alpha
       ! replace "elem%rmetdet = elem%rmetdet / alpha" with the one below,
       ! to ensure that elem%rmetdet = 1/elem%metdet
       ! elem%rmetdet = elem%rmetdet / alpha
       elem%rmetdet = 1.0D0/elem%metdet
       elem%met = elem%met * alpha
       elem%metinv = elem%metinv / alpha
    elseif( cubed_sphere_map == 2 ) then
       ! eps bubble correction for case (3).
       do j=2,np-1
         do i=2,np-1
           elem%D(i,j,:,:) = elem%D(i,j,:,:) * sqrt(alpha)
           elem%Dinv(i,j,:,:) = elem%Dinv(i,j,:,:) / sqrt(alpha)
           elem%metdet(i,j) = elem%metdet(i,j) * alpha
           elem%rmetdet(i,j) = 1.0D0/elem%metdet(i,j)
           elem%met(i,j,:,:) = elem%met(i,j,:,:) * alpha
           elem%metinv(i,j,:,:) = elem%metinv(i,j,:,:) / alpha
         enddo
       enddo
    endif ! end of alpha/eps. bubble correction


  end subroutine metric_atomic


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
    use physical_constants, only : dx_ref, dy_ref

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
