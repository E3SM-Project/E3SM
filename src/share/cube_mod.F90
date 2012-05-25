#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#define _BEGIN_FACE 1
#define _END_FACE   4
#undef _FACE_6
#undef _FACE_5

module cube_mod
  use kinds, only : real_kind, long_kind, longdouble_kind
  use coordinate_systems_mod, only : spherical_polar_t, cartesian3D_t
  use physical_constants, only : dd_pi
  implicit none
  private

  integer,public, parameter :: nfaces = 6          ! number of faces on the cube
  integer,public, parameter :: nInnerElemEdge = 8  ! number of edges for an interior element
  integer,public, parameter :: nCornerElemEdge = 4 ! number of corner elements

  real(kind=real_kind), public, parameter :: cube_xstart = -0.25D0*DD_PI
  real(kind=real_kind), public, parameter :: cube_xend   =  0.25D0*DD_PI
  real(kind=real_kind), public, parameter :: cube_ystart = -0.25D0*DD_PI
  real(kind=real_kind), public, parameter :: cube_yend   =  0.25D0*DD_PI

  type, public :: face_t
     sequence
     type (spherical_polar_t) :: sphere0       ! tangent point of face on sphere
     type (spherical_polar_t) :: sw            ! sw corner of face on sphere
     type (spherical_polar_t) :: se            ! se corner of face on sphere
     type (spherical_polar_t) :: ne            ! ne corner of face on sphere
     type (spherical_polar_t) :: nw            ! nw corner of face on sphere
     type (cartesian3D_t)     :: P0
     type (cartesian3D_t)     :: X0
     type (cartesian3D_t)     :: Y0
     integer                  :: number
     integer                  :: padding       ! padd the struct
  end type face_t

  type, public :: cube_face_coord_t
     sequence
     real(real_kind) :: x             ! x coordinate
     real(real_kind) :: y             ! y coordinate
     type (face_t), pointer :: face     ! face
  end type cube_face_coord_t

  ! ==========================================
  ! Public Interfaces
  ! ==========================================

  public :: CubeTopology

  ! Rotate the North Pole:  used for JW baroclinic test case
  ! Settings this only changes Coriolis.  
  ! User must also rotate initial condition
  real (kind=real_kind), public :: rotate_grid = 0

  ! ===============================
  ! Public methods for cube
  ! ===============================

  public  :: cube_init_atomic
  public  :: convert_gbl_index
  public  :: cube_assemble
  public  :: vmap,dmap
  public  :: covariant_rot
  public  :: contravariant_rot
  public  :: set_corner_coordinates
  public  :: assign_node_numbers_to_elem


  public  :: CubeEdgeCount
  public  :: CubeElemCount
  public  :: CubeSetupEdgeIndex
  public  :: rotation_init_atomic

  ! ===============================
  ! Private methods
  ! ===============================


  private :: coordinates_atomic
  private :: metric_atomic
  private :: coreolis_init, coreolis_init_atomic
  private :: GetLatticeSpacing

contains



  ! =======================================
  !  cube_init_atomic:
  !
  ! Initialize element descriptors for 
  ! cube sphere case for each element ... 
  ! =======================================
  subroutine cube_init_atomic(elem,gll_points,alpha_in)
    use element_mod, only : element_t
    use dimensions_mod, only : np
    type (element_t),intent(inout) :: elem
    real (kind=real_kind),optional :: alpha_in
    real (kind=real_kind)          :: alpha=1
    real (kind=longdouble_kind)      :: gll_points(np)

    if(present(alpha_in)) alpha=alpha_in
    
    elem%FaceNum=elem%vertex%face_number
    call coordinates_atomic(elem,gll_points)

    call metric_atomic(elem,gll_points,alpha)

    call coreolis_init_atomic(elem)
    elem%desc%use_rotation= 0
    call solver_weights_atomic(elem)


  end subroutine cube_init_atomic

  ! =======================================
  ! coordinates_atomic:
  !
  ! Initialize element coordinates for
  ! cube-sphere case ... (atomic) 
  !
  ! =======================================

  subroutine coordinates_atomic(elem,gll_points)
    use element_mod, only : element_t, element_var_coordinates
    use coordinate_systems_mod, only : cartesian2d_t,ref2sphere,cubedsphere2cart, spherical_to_cart, sphere_tri_area
    use dimensions_mod, only : np
    use parallel_mod, only : abortmp
    type (element_t) :: elem
    real (kind=longdouble_kind)      :: gll_points(np)


    real (kind=real_kind)      :: area1,area2
    type (cartesian3d_t) :: quad(4)
    integer face_no,i,j

    ! =========================================
    ! compute coordinates of each GLL point
    ! =========================================
    face_no = elem%vertex%face_number
    do i=1,np
    do j=1,np
       elem%spherep(i,j)=ref2sphere(gll_points(i),gll_points(j),elem%corners,face_no)
    enddo
    enddo

    ! compute the corners in Cartesian coordinates
    do i=1,4
       elem%corners3D(i)=cubedsphere2cart(elem%corners(i),face_no)
    enddo

    ! also compute the [-pi/2,pi/2] cubed sphere coordinates:
    elem%cartp=element_var_coordinates(elem%corners,gll_points)

    ! Matrix describing vector conversion to cartesian
    ! Zonal direction
    elem%vec_sphere2cart(:,:,1,1) =                            -SIN(elem%spherep(:,:)%lon)
    elem%vec_sphere2cart(:,:,2,1) =                             COS(elem%spherep(:,:)%lon)
    elem%vec_sphere2cart(:,:,3,1) =                             0.0_real_kind
    ! Meridional direction
    elem%vec_sphere2cart(:,:,1,2) = -SIN(elem%spherep(:,:)%lat)*COS(elem%spherep(:,:)%lon)
    elem%vec_sphere2cart(:,:,2,2) = -SIN(elem%spherep(:,:)%lat)*SIN(elem%spherep(:,:)%lon)
    elem%vec_sphere2cart(:,:,3,2) =  COS(elem%spherep(:,:)%lat)

  end subroutine coordinates_atomic

  ! elem_jacobians:
  !
  ! Calculate Jacobian associated with mapping
  ! from arbitrary quadrilateral to [-1,1]^2
  ! along with its inverse and determinant
  ! ==========================================

  subroutine elem_jacobians(coords, unif2quadmap, nx)

    use coordinate_systems_mod, only : cartesian2d_t
!    use quadrature_mod, only : quadrature_t

    integer, intent(in) :: nx
    type (cartesian2D_t),  dimension(nx,nx), intent(in) :: coords
    ! unif2quadmap is the bilinear map from [-1,1]^2 -> arbitrary quadrilateral
    real (kind=real_kind), dimension(4,2), intent(out) :: unif2quadmap
!    type (quadrature_t), intent(in) :: gauss
!    real (kind=real_kind), dimension(2,2,nx,nx), intent(out) :: J, Jinv
!    real (kind=real_kind), dimension(nx,nx), intent(out) :: Jdet
    integer :: ii,jj

    unif2quadmap(1,1)=(coords(1,1)%x+coords(nx,1)%x+coords(nx,nx)%x+coords(1,nx)%x)/4.0d0
    unif2quadmap(1,2)=(coords(1,1)%y+coords(nx,1)%y+coords(nx,nx)%y+coords(1,nx)%y)/4.0d0
    unif2quadmap(2,1)=(-coords(1,1)%x+coords(nx,1)%x+coords(nx,nx)%x-coords(1,nx)%x)/4.0d0
    unif2quadmap(2,2)=(-coords(1,1)%y+coords(nx,1)%y+coords(nx,nx)%y-coords(1,nx)%y)/4.0d0
    unif2quadmap(3,1)=(-coords(1,1)%x-coords(nx,1)%x+coords(nx,nx)%x+coords(1,nx)%x)/4.0d0
    unif2quadmap(3,2)=(-coords(1,1)%y-coords(nx,1)%y+coords(nx,nx)%y+coords(1,nx)%y)/4.0d0
    unif2quadmap(4,1)=(coords(1,1)%x-coords(nx,1)%x+coords(nx,nx)%x-coords(1,nx)%x)/4.0d0
    unif2quadmap(4,2)=(coords(1,1)%y-coords(nx,1)%y+coords(nx,nx)%y-coords(1,nx)%y)/4.0d0

  end subroutine elem_jacobians

  ! =========================================
  ! metric_atomic:
  !
  ! Initialize cube-sphere metric terms:
  ! equal angular elements (atomic)
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
  !
  ! so if we want to tweak the mapping by a factor alpha (so he weights add up to 4pi, for example)
  ! we take:
  !    NEW       OLD     
  !       D = sqrt(alpha) D  and then rederive all quantities.  
  !    detD = alpha detD
  !    
  ! where alpha = 4pi/SEMarea, SEMarea = global sum elem(ie)%mv(i,j)*elem(ie)%metdet(i,j)
  ! 
  ! =========================================

  subroutine metric_atomic(elem,gll_points,alpha)
    use element_mod, only : element_t
    use dimensions_mod, only : np
    use physical_constants, only : rrearth
    use parallel_mod, only : abortmp

    type (element_t) :: elem
    real(kind=real_kind) :: alpha
    real (kind=longdouble_kind)      :: gll_points(np)
    ! Local variables
    integer ii,face_no
    integer i,j
    integer iptr

    real (kind=real_kind) :: r         ! distance from origin for point on cube tangent to unit sphere

    real (kind=real_kind) :: const, norm
    real (kind=real_kind) :: detD      ! determinant of vector field mapping matrix.  

    real (kind=real_kind) :: x1        ! 1st cube face coordinate
    real (kind=real_kind) :: x2        ! 2nd cube face coordinate
    real (kind=real_kind) :: tmpD(2,2)
    real (kind=real_kind) :: l1, l2     ! eigen values of met

    face_no = elem%vertex%face_number

    ! ==============================================
    ! Initialize differential mapping operator
    ! to and from vector fields on the sphere to 
    ! contravariant vector fields on the cube
    ! i.e. dM/dx^i in Sadourney (1972) and it's 
    ! inverse
    ! ==============================================

    ! MNL: Calculate Jacobians.  these must be computed before Dmap is used below
    call elem_jacobians(elem%cartp, elem%u2qmap, np)

    elem%max_eig = 0.0d0
    elem%min_eig = 1d99
    do j=1,np
       do i=1,np
          x1=gll_points(i)
          x2=gll_points(j)
          call Dmap(elem%D(:,:,i,j),elem,x1,x2)


          ! Numerical metric tensor based on analytic D: met = D^T times D
          ! (D maps between sphere and reference element)
          elem%met(1,1,i,j) = elem%D(1,1,i,j)*elem%D(1,1,i,j) + &
                              elem%D(2,1,i,j)*elem%D(2,1,i,j)
          elem%met(1,2,i,j) = elem%D(1,1,i,j)*elem%D(1,2,i,j) + &
                              elem%D(2,1,i,j)*elem%D(2,2,i,j)
          elem%met(2,1,i,j) = elem%D(1,1,i,j)*elem%D(1,2,i,j) + &
                              elem%D(2,1,i,j)*elem%D(2,2,i,j)
          elem%met(2,2,i,j) = elem%D(1,2,i,j)*elem%D(1,2,i,j) + &
                              elem%D(2,2,i,j)*elem%D(2,2,i,j)

          ! compute D^-1...
          ! compute determinant of D mapping matrix... if not zero compute inverse

          detD = elem%D(1,1,i,j)*elem%D(2,2,i,j) - elem%D(1,2,i,j)*elem%D(2,1,i,j)      

          elem%Dinv(1,1,i,j) =  elem%D(2,2,i,j)/detD
          elem%Dinv(1,2,i,j) = -elem%D(1,2,i,j)/detD
          elem%Dinv(2,1,i,j) = -elem%D(2,1,i,j)/detD
          elem%Dinv(2,2,i,j) =  elem%D(1,1,i,j)/detD

          ! L2 norm = sqrt max eigenvalue of metinv
          !         = 1/sqrt(min eigenvalue of met)
          ! l1 and l2 are eigenvalues of met
          ! (should both be positive, l1 > l2)
          l1 = (elem%met(1,1,i,j) + elem%met(2,2,i,j) + sqrt(4.0d0*elem%met(1,2,i,j)*elem%met(2,1,i,j) + &
              (elem%met(1,1,i,j) - elem%met(2,2,i,j))**2))/2.0d0
          l2 = (elem%met(1,1,i,j) + elem%met(2,2,i,j) - sqrt(4.0d0*elem%met(1,2,i,j)*elem%met(2,1,i,j) + &
              (elem%met(1,1,i,j) - elem%met(2,2,i,j))**2))/2.0d0
          ! Max L2 norm of Dinv is sqrt of max eigenvalue of metinv
          ! max eigenvalue of metinv is 1/min eigenvalue of met
          norm = 1.0d0/sqrt(min(abs(l1),abs(l2)))
          elem%max_eig = max(norm, elem%max_eig)
          ! Min L2 norm of Dinv is sqrt of min eigenvalue of metinv
          ! min eigenvalue of metinv is 1/max eigenvalue of met
          norm = 1.0d0/sqrt(max(abs(l1),abs(l2)))
          elem%min_eig = min(norm, elem%min_eig)

          ! Need inverse of met if not calculated analytically
          elem%metdet(i,j) = abs(detD)
          elem%rmetdet(i,j) = 1.0D0/abs(detD)

          elem%metinv(1,1,i,j) =  elem%met(2,2,i,j)/(detD*detD)
          elem%metinv(1,2,i,j) = -elem%met(1,2,i,j)/(detD*detD)
          elem%metinv(2,1,i,j) = -elem%met(2,1,i,j)/(detD*detD)
          elem%metinv(2,2,i,j) =  elem%met(1,1,i,j)/(detD*detD)
          
       end do
    end do

    elem%dx_short = 1.0d0/(elem%max_eig*0.5d0*dble(np-1)*rrearth*1000.0d0)
    elem%dx_long  = 1.0d0/(elem%min_eig*0.5d0*dble(np-1)*rrearth*1000.0d0)
    ! ===============================================
    !
    ! Initialize equal angular metric tensor on each 
    ! on velocity grid for unit sphere.
    !
    ! Initialize gdet = SQRT(ABS(DET(gij)))
    !
    ! These quantities are the same on every face
    ! of the cube.
    !
    ! =================================================

    ! mt: better might be to compute all these quantities directly from D
    ! for consistency?
    !
    ! MNL: done
    elem%D = elem%D * sqrt(alpha) 
    elem%Dinv = elem%Dinv / sqrt(alpha) 
    elem%metdet = elem%metdet * alpha
    elem%rmetdet = elem%rmetdet / alpha
    elem%met = elem%met * alpha
    elem%metinv = elem%metinv / alpha

  end subroutine metric_atomic

  ! =======================================
  ! solver_weights:
  !
  ! For nonstaggered GaussLobatto elements,
  ! compute weights for redundant points in 
  ! cg solver.
  !
  ! =======================================

  subroutine solver_weights_atomic(elem)
    use element_mod, only : element_t
    use dimensions_mod, only : np
    use parallel_mod, only : abortmp

    type (element_t) :: elem
    real (kind=real_kind) :: x 

    ! Local variables

    integer :: i, j
    ! =========================================
    ! compute cube face coordinates of element
    ! =========================================

    do i=1,np
      do j=1,np
        if (i==1) then
          if (j==1) then
             x = 1.0_real_kind/elem%node_multiplicity(1)
          else if (j==np) then
             x = 1.0_real_kind/elem%node_multiplicity(4)
          else
             x = 0.5_real_kind
          end if
        else if (i==np) then
          if (j==1) then
             x = 1.0_real_kind/elem%node_multiplicity(2)
          else if (j==np) then
             x = 1.0_real_kind/elem%node_multiplicity(3)
          else
             x = 0.5_real_kind
          end if
        else if (j==1 .or. j==np) then
           x = 0.5_real_kind
        else
           x = 1.0_real_kind
        end if
        elem%solver_wts(i,j) = x
      end do
    end do

  end subroutine solver_weights_atomic

#if 1
  ! ========================================
  ! covariant_rot:
  !
  ! 2 x 2 matrix multiply:  Db^T * Da^-T
  ! for edge rotations: maps face a to face b
  !
  ! ========================================

  function covariant_rot(Da,Db) result(R)

    real (kind=real_kind) :: Da(2,2)
    real (kind=real_kind) :: Db(2,2)
    real (kind=real_kind) :: R(2,2)

    real (kind=real_kind) :: detDa

    detDa = Da(2,2)*Da(1,1) - Da(1,2)*Da(2,1)

    R(1,1)=(Da(2,2)*Db(1,1) - Da(1,2)*Db(2,1))/detDa
    R(1,2)=(Da(1,1)*Db(2,1) - Da(2,1)*Db(1,1))/detDa
    R(2,1)=(Da(2,2)*Db(1,2) - Da(1,2)*Db(2,2))/detDa
    R(2,2)=(Da(1,1)*Db(2,2) - Da(2,1)*Db(1,2))/detDa

  end function covariant_rot
#else

  ! ========================================
  ! covariant_rot:
  !
  ! 2 x 2 matrix multiply:  Db * Da^-1
  ! for edge rotations: maps face a to face b
  !
  ! ========================================

  function covariant_rot(Da,Db) result(R)

    real (kind=real_kind) :: Da(2,2)
    real (kind=real_kind) :: Db(2,2)
    real (kind=real_kind) :: R(2,2)

    real (kind=real_kind) :: detDa

    detDa = Da(2,2)*Da(1,1) - Da(1,2)*Da(2,1)

    R(1,1)=(Da(2,2)*Db(1,1) - Da(2,1)*Db(1,2))/detDa
    R(1,2)=(Da(1,1)*Db(1,2) - Da(1,2)*Db(1,1))/detDa
    R(2,1)=(Da(2,2)*Db(2,1) - Da(2,1)*Db(2,2))/detDa
    R(2,2)=(Da(1,1)*Db(2,2) - Da(1,2)*Db(2,1))/detDa

  end function covariant_rot

#endif

  ! ========================================
  ! contravariant_rot:
  !
  ! 2 x 2 matrix multiply:  Db^-1 * Da
  ! that maps a contravariant vector field
  ! from an edge of cube face a to a contiguous 
  ! edge of cube face b.
  !
  ! ========================================

  function contravariant_rot(Da,Db) result(R)

    real (kind=real_kind) :: Da(2,2)
    real (kind=real_kind) :: Db(2,2)
    real (kind=real_kind) :: R(2,2)

    real (kind=real_kind) :: detDb

    detDb = Db(2,2)*Db(1,1) - Db(1,2)*Db(2,1)

    R(1,1)=(Da(1,1)*Db(2,2) - Da(2,1)*Db(1,2))/detDb
    R(1,2)=(Da(1,2)*Db(2,2) - Da(2,2)*Db(1,2))/detDb
    R(2,1)=(Da(2,1)*Db(1,1) - Da(1,1)*Db(2,1))/detDb
    R(2,2)=(Da(2,2)*Db(1,1) - Da(1,2)*Db(2,1))/detDb

  end function contravariant_rot

  ! ========================================================
  ! Dmap:
  !
  ! Initialize mapping that tranforms contravariant 
  ! vector fields on the reference element onto vector fields on
  ! the sphere. 
  ! For Gnomonic, followed by bilinear, this code uses the old vmap()
  ! for unstructured grids, this code uses the parametric map that
  ! maps quads on the sphere directly to the reference element
  ! ========================================================
  subroutine Dmap(D, elem, a,b)
    use coordinate_systems_mod, only : dist_threshold
    use element_mod, only : element_t
    type (element_t) :: elem
    real (kind=real_kind), intent(out)  :: D(2,2)
    real (kind=real_kind), intent(in)     :: a,b
    ! local
    real (kind=real_kind)  :: tmpD(2,2), Jp(2,2),x1,x2,pi,pj,qi,qj

    ! input (a,b) shold be a point in the reference element [-1,1]
    ! compute Jp(a,b)
    Jp(1,1) = elem%u2qmap(2,1) + elem%u2qmap(4,1)*b
    Jp(1,2) = elem%u2qmap(3,1) + elem%u2qmap(4,1)*a
    Jp(2,1) = elem%u2qmap(2,2) + elem%u2qmap(4,2)*b
    Jp(2,2) = elem%u2qmap(3,2) + elem%u2qmap(4,2)*a

    ! map (a,b) to the [-pi/2,pi/2] equi angular cube face:  x1,x2
    ! a = gp%points(i)
    ! b = gp%points(j)
    pi = (1-a)/2
    pj = (1-b)/2
    qi = (1+a)/2
    qj = (1+b)/2
    x1 = pi*pj*elem%corners(1)%x &
         + qi*pj*elem%corners(2)%x &
         + qi*qj*elem%corners(3)%x &
         + pi*qj*elem%corners(4)%x 
    x2 = pi*pj*elem%corners(1)%y &
         + qi*pj*elem%corners(2)%y &
         + qi*qj*elem%corners(3)%y &
         + pi*qj*elem%corners(4)%y 
    
    call vmap(tmpD,x1,x2,elem%vertex%face_number)

    ! Include map from element -> ref element in D
    D(1,1) = tmpD(1,1)*Jp(1,1) + tmpD(1,2)*Jp(2,1)
    D(1,2) = tmpD(1,1)*Jp(1,2) + tmpD(1,2)*Jp(2,2)
    D(2,1) = tmpD(2,1)*Jp(1,1) + tmpD(2,2)*Jp(2,1)
    D(2,2) = tmpD(2,1)*Jp(1,2) + tmpD(2,2)*Jp(2,2)
  end subroutine Dmap



  ! ========================================================
  ! vmap:
  !
  ! Initialize mapping that tranforms contravariant 
  ! vector fields on the cube onto vector fields on
  ! the sphere. This follows Taylor's D matrix 
  !
  !       | cos(theta)dlambda/dx1  cos(theta)dlambda/dx2 |
  !   D = |                                              |
  !       |     dtheta/dx1              dtheta/dx2       |
  !
  ! ========================================================

  subroutine vmap(D, x1, x2, face_no) 
    use coordinate_systems_mod, only : dist_threshold
    real (kind=real_kind), intent(inout)  :: D(2,2)
    real (kind=real_kind), intent(in)     :: x1
    real (kind=real_kind), intent(in)     :: x2
    integer              , intent(in)     :: face_no

    ! Local variables

    real (kind=real_kind) :: poledist  ! SQRT(TAN(x1)**2 +TAN(x2)**2)
    real (kind=real_kind) :: r         ! distance from cube point to center of sphere

    real (kind=real_kind) :: D11
    real (kind=real_kind) :: D12
    real (kind=real_kind) :: D21
    real (kind=real_kind) :: D22

    r=SQRT(1.0D0 + TAN(x1)**2 + TAN(x2)**2)

    if (face_no >= 1 .and. face_no <= 4) then

       D11 = 1.0D0/(r*COS(x1))
       D12 = 0.0D0
       D21 = -TAN(x1)*TAN(x2)/(COS(x1)*r*r)        
       D22 = 1.0D0/(r*r*COS(x1)*COS(x2)*COS(x2))

       D(1,1) =  D11
       D(1,2) =  D12
       D(2,1) =  D21
       D(2,2) =  D22


    else if (face_no==6) then
       poledist=SQRT( TAN(x1)**2 + TAN(x2)**2)
       if ( poledist <= DIST_THRESHOLD ) then

          ! we set the D transform to the identity matrix 
          ! which works ONLY for swtc1, phi starting at 
          ! 3*PI/2... assumes lon at pole == 0

          D(1,1) =  1.0D0
          D(1,2) =  0.0D0
          D(2,1) =  0.0D0
          D(2,2) =  1.0D0

       else

          D11 = -TAN(x2)/(poledist*COS(x1)*COS(x1)*r)
          D12 =  TAN(x1)/(poledist*COS(x2)*COS(x2)*r)
          D21 = -TAN(x1)/(poledist*COS(x1)*COS(x1)*r*r)
          D22 = -TAN(x2)/(poledist*COS(x2)*COS(x2)*r*r)

          D(1,1) =  D11
          D(1,2) =  D12
          D(2,1) =  D21
          D(2,2) =  D22

       end if
    else if (face_no==5) then
       poledist=SQRT( TAN(x1)**2 + TAN(x2)**2)
       if ( poledist <= DIST_THRESHOLD ) then

          ! we set the D transform to the identity matrix 
          ! which works ONLY for swtc1, phi starting at 
          ! 3*PI/2... assumes lon at pole == 0, i.e. very specific

          D(1,1) =  1.0D0
          D(1,2) =  0.0D0
          D(2,1) =  0.0D0
          D(2,2) =  1.0D0

       else

          D11 =  TAN(x2)/(poledist*COS(x1)*COS(x1)*r)
          D12 = -TAN(x1)/(poledist*COS(x2)*COS(x2)*r)
          D21 =  TAN(x1)/(poledist*COS(x1)*COS(x1)*r*r)
          D22 =  TAN(x2)/(poledist*COS(x2)*COS(x2)*r*r)

          D(1,1) =  D11
          D(1,2) =  D12
          D(2,1) =  D21
          D(2,2) =  D22

       end if
    end if

  end subroutine vmap

  ! ========================================
  ! coreolis_init:
  !
  ! Initialize coreolis term ...
  !
  ! ========================================

  subroutine coreolis_init(elem)
    use element_mod, only : element_t
    use dimensions_mod, only : np
    use physical_constants, only : omega
    type (element_t) :: elem(:)

    ! Local variables

    integer                  :: i,j
    integer                  :: ii

    do ii=1,SIZE(elem)
       call coreolis_init_atomic(elem(ii))
    end do

  end subroutine coreolis_init

  ! ========================================
  ! coreolis_init_atomic:
  !
  ! Initialize coreolis term ...
  !
  ! ========================================

  subroutine coreolis_init_atomic(elem)
    use element_mod, only : element_t
    use dimensions_mod, only : np
    use physical_constants, only : omega

    type (element_t) :: elem

    ! Local variables

    integer                  :: i,j
    real (kind=real_kind) :: lat,lon,rangle

    rangle = rotate_grid*DD_PI/180
    do j=1,np
       do i=1,np
             if ( rotate_grid /= 0) then
                lat = elem%spherep(i,j)%lat
                lon = elem%spherep(i,j)%lon
             	elem%fcor(i,j)= 2*omega* &
                     (-cos(lon)*cos(lat)*sin(rangle) + sin(lat)*cos(rangle))
             else
                elem%fcor(i,j) = 2.0D0*omega*SIN(elem%spherep(i,j)%lat)
             endif
       end do
    end do

  end subroutine coreolis_init_atomic

  ! =========================================
  ! rotation_init_atomic:
  !
  ! Initialize cube rotation terms resulting
  ! from changing cube face coordinate systems
  !
  ! =========================================

  subroutine rotation_init_atomic(elem, rot_type)
    use element_mod, only : element_t
    use dimensions_mod, only : np
    use control_mod, only : north, south, east, west, neast, seast, swest, nwest
    use parallel_mod, only : abortmp

    type (element_t) :: elem
    character(len=*) rot_type

    ! =======================================
    ! Local variables
    ! =======================================

    integer :: myface_no        ! current element face number
    integer :: nbrface_no       ! neighbor element face number
    integer :: inbr
    integer :: nrot,irot
    integer :: ii,i,j,k
    integer :: ir,jr

    real (kind=real_kind) :: Dloc(2,2,np)
    real (kind=real_kind) :: Drem(2,2,np)
    real (kind=real_kind) :: x1,x2


    myface_no = elem%vertex%face_number

    nrot   = 0

    do inbr=1,8
       if (ASSOCIATED(elem%vertex%nbrs(inbr)%n)) then
          do k=1,SIZE(elem%vertex%nbrs(inbr)%n)
             nbrface_no = elem%vertex%nbrs(inbr)%f(k)
             if (myface_no /= nbrface_no) nrot=nrot+1
          end do
       end if
    end do

    if(associated(elem%desc%rot)) then
       if (size(elem%desc%rot) > 0) then
          !         deallocate(elem%desc%rot)
          NULLIFY(elem%desc%rot)
       endif
    end if

    ! =====================================================
    ! If there are neighbors on other cube faces, allocate 
    ! an array of rotation matrix structs.
    ! =====================================================

    if (nrot > 0) then
       allocate(elem%desc%rot(nrot))
       elem%desc%use_rotation=1
       irot=0          
       do inbr=1,8
       if (ASSOCIATED(elem%vertex%nbrs(inbr)%n)) then
       do k=1,SIZE(elem%vertex%nbrs(inbr)%n)

          nbrface_no = elem%vertex%nbrs(inbr)%f(k)
          ! The cube edge (myface_no,nbrface_no) and inbr defines 
          ! a unique rotation given by (D^-1) on myface_no x (D) on nbrface_no

          if (myface_no /= nbrface_no .and. elem%vertex%nbrs(inbr)%n(k) /= -1 ) then           

             irot=irot+1

             if (inbr <= 4) then      
                allocate(elem%desc%rot(irot)%R(2,2,np))  ! edge
             else                     
                allocate(elem%desc%rot(irot)%R(2,2,1 ))   ! corner
             end if

             ! must compute Dloc on my face, Drem on neighbor face, 
             ! for each point on edge or corner.

             ! ==================================== 
             ! Equatorial belt east/west neighbors
             ! ==================================== 

             if (nbrface_no <= 4 .and. myface_no <= 4) then

                if (inbr == west) then
                   do j=1,np
                      x1 = elem%cartp(1,j)%x
                      x2 = elem%cartp(1,j)%y
                      call Vmap(Dloc(1,1,j), x1,x2,myface_no)
                      call Vmap(Drem(1,1,j),-x1,x2,nbrface_no)
                   end do
                else if (inbr == east) then
                   do j=1,np
                      x1 = elem%cartp(np,j)%x
                      x2 = elem%cartp(np,j)%y
                      call Vmap(Dloc(1,1,j), x1,x2,myface_no)
                      call Vmap(Drem(1,1,j),-x1,x2,nbrface_no)
                   end do
                else if (inbr == swest ) then
                   x1 = elem%cartp(1,1)%x
                   x2 = elem%cartp(1,1)%y
                   call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                   call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                else if (inbr == nwest ) then
                   x1 = elem%cartp(1,np)%x
                   x2 = elem%cartp(1,np)%y
                   call Vmap(Dloc(1,1,1), x1,x2,myface_no)
                   call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                else if (inbr == seast ) then
                   x1 = elem%cartp(np,1)%x
                   x2 = elem%cartp(np,1)%y
                   call Vmap(Dloc(1,1,1), x1,x2,myface_no)
                   call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                else if (inbr == neast ) then
                   x1 = elem%cartp(np,np)%x
                   x2 = elem%cartp(np,np)%y
                   call Vmap(Dloc(1,1,1), x1,x2,myface_no)
                   call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                end if

             end if

             ! Northern Neighbors of Equatorial Belt

             if ( myface_no <= 4 .and. nbrface_no == 6 ) then
                if (inbr == north) then
                   do i=1,np
                      ir=np+1-i
                      x1 = elem%cartp(i,np)%x
                      x2 = elem%cartp(i,np)%y
                      if ( myface_no == 1) then
                         call Vmap(Dloc(1,1,i), x1,x2,myface_no)
                         call Vmap(Drem(1,1,i),x1,-x2,nbrface_no)
                      end if
                      if ( myface_no == 2) then
                         call Vmap(Dloc(1,1,i),x1,x2,myface_no)
                         call Vmap(Drem(1,1,i),x2,x1,nbrface_no)

                      end if
                      if ( myface_no == 3) then
                         call Vmap(Dloc(1,1,ir), x1,x2,myface_no)
                         call Vmap(Drem(1,1,ir),-x1,x2,nbrface_no)
                      end if
                      if ( myface_no == 4) then
                         call Vmap(Dloc(1,1,ir), x1,x2,myface_no)
                         call Vmap(Drem(1,1,ir),-x2,-x1,nbrface_no)
                      end if
                   end do
                else if (inbr == nwest) then
                   x1 = elem%cartp(1,np)%x
                   x2 = elem%cartp(1,np)%y
                   call Vmap(Dloc(1,1,1), x1,x2,myface_no)
                   if ( myface_no == 1) call Vmap(Drem(1,1,1),x1,-x2,nbrface_no)
                   if ( myface_no == 2) call Vmap(Drem(1,1,1),x2, x1,nbrface_no)
                   if ( myface_no == 3) call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                   if ( myface_no == 4) call Vmap(Drem(1,1,1),-x2,-x1,nbrface_no)
                else if (inbr == neast) then
                   x1 = elem%cartp(np,np)%x
                   x2 = elem%cartp(np,np)%y
                   call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                   if ( myface_no == 1) call Vmap(Drem(1,1,1),x1,-x2,nbrface_no)
                   if ( myface_no == 2) call Vmap(Drem(1,1,1),x2, x1,nbrface_no)
                   if ( myface_no == 3) call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                   if ( myface_no == 4) call Vmap(Drem(1,1,1),-x2,-x1,nbrface_no)
                end if

             end if

             ! Southern Neighbors of Equatorial Belt

             if ( myface_no <= 4 .and. nbrface_no == 5 ) then
                if (inbr == south) then
                   do i=1,np
                      ir=np+1-i
                      x1 = elem%cartp(i,1)%x
                      x2 = elem%cartp(i,1)%y
                      if ( myface_no == 1) then
                         call Vmap(Dloc(1,1,i), x1, x2,myface_no)
                         call Vmap(Drem(1,1,i), x1,-x2,nbrface_no)
                      end if
                      if ( myface_no == 2) then
                         call Vmap(Dloc(1,1,ir),x1,x2,myface_no)
                         call Vmap(Drem(1,1,ir),-x2,-x1,nbrface_no)
                      end if
                      if ( myface_no == 3) then
                         call Vmap(Dloc(1,1,ir), x1,x2,myface_no)
                         call Vmap(Drem(1,1,ir),-x1,x2,nbrface_no)
                      end if
                      if ( myface_no == 4) then
                         call Vmap(Dloc(1,1,i), x1,x2,myface_no)
                         call Vmap(Drem(1,1,i), x2,x1,nbrface_no)
                      end if
                   end do
                else if (inbr == swest) then
                   x1 = elem%cartp(1,1)%x
                   x2 = elem%cartp(1,1)%y
                   call Vmap(Dloc(1,1,1),x1,x2,myface_no)


                   if ( myface_no == 1) call Vmap(Drem(1,1,1),x1,-x2,nbrface_no)
                   if ( myface_no == 2) call Vmap(Drem(1,1,1),-x2,-x1,nbrface_no)
                   if ( myface_no == 3) call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                   if ( myface_no == 4) call Vmap(Drem(1,1,1),x2,x1,nbrface_no)

                else if (inbr == seast) then
                   x1 = elem%cartp(np,1)%x
                   x2 = elem%cartp(np,1)%y
                   call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                   if ( myface_no == 1) call Vmap(Drem(1,1,1),x1,-x2,nbrface_no)
                   if ( myface_no == 2) call Vmap(Drem(1,1,1),-x2,-x1,nbrface_no)
                   if ( myface_no == 3) call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                   if ( myface_no == 4) call Vmap(Drem(1,1,1),x2,x1,nbrface_no)
                end if

             end if

             ! Neighbors of Northern Capping Face Number 6

             if ( myface_no == 6 ) then
                if (nbrface_no == 1) then
                   if (inbr == south) then
                      do i=1,np
                         x1 = elem%cartp(i,1)%x
                         x2 = elem%cartp(i,1)%y
                         call Vmap(Dloc(1,1,i),x1,x2,myface_no)
                         call Vmap(Drem(1,1,i),x1,-x2,nbrface_no)
                      end do
                   else if (inbr == swest) then
                      x1 = elem%cartp(1,1)%x
                      x2 = elem%cartp(1,1)%y
                      call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                      call Vmap(Drem(1,1,1),x1,-x2,nbrface_no)
                   else if (inbr == seast) then
                      x1 = elem%cartp(np,1)%x
                      x2 = elem%cartp(np,1)%y
                      call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                      call Vmap(Drem(1,1,1),x1,-x2,nbrface_no)
                   end if
                else if (nbrface_no == 2) then
                   if (inbr == east) then
                      do j=1,np
                         x1 = elem%cartp(np,j)%x
                         x2 = elem%cartp(np,j)%y
                         call Vmap(Dloc(1,1,j),x1,x2,myface_no)
                         call Vmap(Drem(1,1,j),x2,x1,nbrface_no)
                      end do
                   else if (inbr == seast) then
                      x1 = elem%cartp(np,1)%x
                      x2 = elem%cartp(np,1)%y
                      call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                      call Vmap(Drem(1,1,1),x2,x1,nbrface_no)
                   else if (inbr == neast) then
                      x1 = elem%cartp(np,np)%x
                      x2 = elem%cartp(np,np)%y
                      call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                      call Vmap(Drem(1,1,1),x2,x1,nbrface_no)
                   end if
                else if (nbrface_no == 3) then
                   if (inbr == north) then
                      do i=1,np
                         ir =np+1-i
                         x1 = elem%cartp(i,np)%x
                         x2 = elem%cartp(i,np)%y
                         call Vmap(Dloc(1,1,ir),x1,x2,myface_no)
                         call Vmap(Drem(1,1,ir),-x1,x2,nbrface_no)
                      end do
                   else if (inbr == nwest) then
                      x1 = elem%cartp(1,np)%x
                      x2 = elem%cartp(1,np)%y
                      call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                      call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                   else if (inbr == neast) then
                      x1 = elem%cartp(np,np)%x
                      x2 = elem%cartp(np,np)%y
                      call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                      call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                   end if
                else if (nbrface_no == 4) then
                   if (inbr == west) then
                      do j=1,np
                         jr=np+1-j
                         x1 = elem%cartp(1,j)%x
                         x2 = elem%cartp(1,j)%y
                         call Vmap(Dloc(1,1,jr), x1, x2,myface_no )
                         call Vmap(Drem(1,1,jr),-x2,-x1,nbrface_no)
                      end do
                   else if (inbr == swest) then
                      x1 = elem%cartp(1,1)%x
                      x2 = elem%cartp(1,1)%y
                      call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                      call Vmap(Drem(1,1,1),-x2,-x1,nbrface_no)
                   else if (inbr == nwest) then
                      x1 = elem%cartp(1,np)%x
                      x2 = elem%cartp(1,np)%y
                      call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                      call Vmap(Drem(1,1,1),-x2,-x1,nbrface_no)
                   end if
                end if
             end if

             ! Neighbors of South Capping Face Number 5

             if ( myface_no == 5 ) then
                if (nbrface_no == 1) then
                   if (inbr == north) then
                      do i=1,np
                         x1 = elem%cartp(i,np)%x
                         x2 = elem%cartp(i,np)%y
                         call Vmap(Dloc(1,1,i),x1,x2,myface_no)
                         call Vmap(Drem(1,1,i),x1,-x2,nbrface_no)
                      end do
                   else if (inbr == nwest) then
                      x1 = elem%cartp(1,np)%x
                      x2 = elem%cartp(1,np)%y
                      call Vmap(Dloc(:,:,1),x1,x2,myface_no)
                      call Vmap(Drem(:,:,1),x1,-x2,nbrface_no)
                   else if (inbr == neast) then
                      x1 = elem%cartp(np,np)%x
                      x2 = elem%cartp(np,np)%y
                      call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                      call Vmap(Drem(1,1,1),x1,-x2,nbrface_no)
                   end if
                else if (nbrface_no == 2) then
                   if (inbr == east) then
                      do j=1,np
                         jr=np+1-j
                         x1 = elem%cartp(np,j)%x
                         x2 = elem%cartp(np,j)%y
                         call Vmap(Dloc(1,1,jr),x1,  x2,myface_no)
                         call Vmap(Drem(1,1,jr),-x2,-x1,nbrface_no)
                      end do
                   else if (inbr == seast) then
                      x1 = elem%cartp(np,1)%x
                      x2 = elem%cartp(np,1)%y
                      call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                      call Vmap(Drem(1,1,1),-x2,-x1,nbrface_no)
                   else if (inbr == neast) then
                      x1 = elem%cartp(np,np)%x
                      x2 = elem%cartp(np,np)%y
                      call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                      call Vmap(Drem(1,1,1),-x2,-x1,nbrface_no)
                   end if
                else if (nbrface_no == 3) then
                   if (inbr == south) then
                      do i=1,np
                         ir=np+1-i
                         x1 = elem%cartp(i,1)%x
                         x2 = elem%cartp(i,1)%y
                         call Vmap(Dloc(1,1,ir),x1,x2,myface_no)
                         call Vmap(Drem(1,1,ir),-x1,x2,nbrface_no)
                      end do
                   else if (inbr == swest) then
                      x1 = elem%cartp(1,1)%x
                      x2 = elem%cartp(1,1)%y
                      call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                      call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                   else if (inbr == seast) then
                      x1 = elem%cartp(np,1)%x
                      x2 = elem%cartp(np,1)%y
                      call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                      call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                   end if
                else if (nbrface_no == 4) then
                   if (inbr == west) then
                      do j=1,np
                         x1 = elem%cartp(1,j)%x
                         x2 = elem%cartp(1,j)%y
                         call Vmap(Dloc(1,1,j),x1,x2,myface_no)
                         call Vmap(Drem(1,1,j),x2,x1,nbrface_no)
                      end do
                   else if (inbr == swest) then
                      x1 = elem%cartp(1,1)%x
                      x2 = elem%cartp(1,1)%y
                      call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                      call Vmap(Drem(1,1,1),x2,x1,nbrface_no)
                   else if (inbr == nwest) then
                      x1 = elem%cartp(1,np)%x
                      x2 = elem%cartp(1,np)%y
                      call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                      call Vmap(Drem(1,1,1),x2,x1,nbrface_no)
                   end if
                end if
             end if

             elem%desc%rot(irot)%nbr = inbr
             if (rot_type == "covariant") then
                do i=1,SIZE(elem%desc%rot(irot)%R(:,:,:),3)
                   elem%desc%rot(irot)%R(:,:,i)=covariant_rot(Dloc(:,:,i),Drem(:,:,i))
                end do
             else if (rot_type == "contravariant") then
                do i=1,SIZE(elem%desc%rot(irot)%R(:,:,:),3)
                   elem%desc%rot(irot)%R(:,:,i)=contravariant_rot(Dloc(:,:,i),Drem(:,:,i))
                end do
             end if

          endif
       end do
       end if
       end do
    end if

  end subroutine rotation_init_atomic

  subroutine set_corner_coordinates(elem)
    use element_mod,    only : element_t 
    use dimensions_mod, only : ne
    use parallel_mod, only : abortmp 
 
    type (element_t) :: elem 

    ! Local variables
    integer  i,ie,je,face_no,nn
    real (kind=real_kind)  :: dx,dy, startx, starty

    if (0==ne) call abortmp('Error in set_corner_coordinates: ne is zero')

    ! ========================================
    ! compute cube face coordinates of element
    ! =========================================

    call convert_gbl_index(elem%vertex%number,ie,je,face_no)

    elem%vertex%face_number = face_no 
    dx = (cube_xend-cube_xstart)/ne
    dy = (cube_yend-cube_ystart)/ne

    startx = cube_xstart+ie*dx
    starty = cube_ystart+je*dy

    elem%corners(1)%x = startx
    elem%corners(1)%y = starty
    elem%corners(2)%x = startx+dx
    elem%corners(2)%y = starty
    elem%corners(3)%x = startx+dx
    elem%corners(3)%y = starty+dy
    elem%corners(4)%x = startx   
    elem%corners(4)%y = starty+dy

    do i=1,4
       elem%node_multiplicity(i) = 4
    end do  
    ie = ie + 1
    je = je + 1
    if      (ie ==  1 .and. je ==  1) then 
       elem%node_multiplicity(1) = 3
    else if (ie == ne .and. je ==  1) then 
       elem%node_multiplicity(2) = 3
    else if (ie == ne .and. je == ne) then
       elem%node_multiplicity(3) = 3
    else if (ie ==  1 .and. je == ne) then
       elem%node_multiplicity(4) = 3
    end if  

  end subroutine set_corner_coordinates

  subroutine assign_node_numbers_to_elem(elements, GridVertex)
    use dimensions_mod, only : ne
    use element_mod,    only : element_t
    use control_mod,    only : north, south, east, west, neast, seast, swest, nwest
    use parallel_mod,   only : abortmp
    use gridgraph_mod,  only : GridVertex_t
    implicit none
    type (element_t), intent(inout)    :: elements(:)
    type (GridVertex_t), intent(in)    :: GridVertex(:)

    type (GridVertex_t)                :: vertex
    integer                            :: connectivity(6*ne*ne, 4)
    integer                            :: nn(4), en(4)
    integer el, i, n, direction
    integer current_node_num, tot_ne

    current_node_num = 0
    tot_ne = 6*ne*ne

    if (0==ne)      call abortmp('Error in assign_node_numbers_to_elem: ne is zero')
    if (tot_ne /= SIZE(GridVertex)) call abortmp('Error in assign_node_numbers_to_elem: GridVertex not correct length')

    connectivity = 0 

    do el = 1,tot_ne  
       vertex = GridVertex(el)
       en = 0 
       do direction = 1,8
          if (ASSOCIATED(vertex%nbrs(direction)%n)) then
             do i=1,SIZE(vertex%nbrs(direction)%n)
                n = vertex%nbrs(direction)%n(i)
                if (n /= -1) then
                   nn = connectivity(n,:)
                   select case (direction)
                   case (north)
                      if (nn(1)/=0) en(4) = nn(1)
                      if (nn(2)/=0) en(3) = nn(2)
                   case (south)
                      if (nn(4)/=0) en(1) = nn(4)
                      if (nn(3)/=0) en(2) = nn(3)
                   case (east)
                      if (nn(1)/=0) en(2) = nn(1)
                      if (nn(4)/=0) en(3) = nn(4)
                   case (west)
                      if (nn(2)/=0) en(1) = nn(2)
                      if (nn(3)/=0) en(4) = nn(3)
                   case (neast)
                      if (nn(1)/=0) en(3) = nn(1)
                   case (seast)
                      if (nn(4)/=0) en(2) = nn(4)
                   case (swest)
                      if (nn(3)/=0) en(1) = nn(3)
                   case (nwest)
                      if (nn(2)/=0) en(4) = nn(2)
                   end select
                end if
             end do
          end if
       end do
       do i=1,4
          if (en(i) == 0) then
             current_node_num = current_node_num + 1
             en(i) = current_node_num
          end if
       end do
       connectivity(el,:) = en
    end do
    if (current_node_num /= (6*ne*ne+2)) then
       call abortmp('Error in assignment of node numbers: Failed Euler test')
    end if
    do el = 1,SIZE(elements)
      elements(el)%node_numbers = connectivity(elements(el)%vertex%number, :)
    end do
  end subroutine assign_node_numbers_to_elem

  ! ================================================
  ! convert_gbl_index:
  !
  ! Convert global element index to cube index
  ! ================================================

  subroutine convert_gbl_index(number,ie,je,face_no)
    use dimensions_mod, only : ne
    use parallel_mod, only : abortmp
    integer, intent(in)  :: number
    integer, intent(out) :: ie,je,face_no

    if (0==ne) call abortmp('Error in cube_mod:convert_gbl_index: ne is zero')

    !  inverse of the function:      number = 1 + ie + ne*je + ne*ne*(face_no-1)
    face_no=((number-1)/(ne*ne))+1
    ie=MODULO(number-1,ne)
    je=(number-1)/ne - (face_no-1)*ne

  end subroutine convert_gbl_index
   
  subroutine CubeTopology(GridEdge, GridVertex)
    use params_mod, only : RECURSIVE, SFCURVE
    use control_mod, only: partmethod
    use gridgraph_mod, only : GridEdge_t, GridVertex_t, initgridedge
    use dimensions_mod, only : np, ne
    use spacecurve_mod, only :  IsFactorable, genspacecurve
    use control_mod, only : north, south, east, west, neast, seast, swest, nwest
    use parallel_mod, only : abortmp
    !-----------------------
    implicit none

    type (GridEdge_t),   intent(out),target     :: GridEdge(:)
    type (GridVertex_t), intent(out),target     :: GridVertex(:)


    integer,allocatable       :: Mesh(:,:)
    integer,allocatable       :: Mesh2(:,:),Mesh2_map(:,:,:),sfcij(:,:)
    type (GridVertex_t),allocatable        :: GridElem(:,:,:)
    integer                   :: i,j,k,l,number,irev,ne2,i2,j2,sfc_index
    integer                   :: EdgeWgtP,CornerWgt
    integer                   :: ielem,nedge
    integer                   :: offset, ierr

    if (0==ne) call abortmp('Error in CubeTopology: ne is zero')

    allocate(GridElem(ne,ne,nfaces),stat=ierr)
    if(ierr/=0) then
       call abortmp('error in allocation of GridElem structure')
    end if

    number=1
    EdgeWgtP   = np
    CornerWgt = 1
    do k=1,nfaces
       do j=1,ne
          do i=1,ne
             ! ====================================
             ! Number elements
             ! ====================================

             ! Do some initalization here
             do l=1,8
               NULLIFY(GridElem(i,j,k)%nbrs(l)%n)
               NULLIFY(GridElem(i,j,k)%nbrs(l)%f)
             end do
             GridElem(i,j,k)%wgtP(:)=0
             GridElem(i,j,k)%wgtP_ghost(:)=1  ! always this value
             GridElem(i,j,k)%SpaceCurve=0
             GridElem(i,j,k)%number=number 
             number=number+1

          end do
       end do
    end do

    !    print *,'CubeTopology: Ne, IsFactorable, IsLoadBalanced : ',ne,IsFactorable(ne),IsLoadBalanced(nelem,npart)

    allocate(Mesh(ne,ne))
    if(IsFactorable(ne)) then
       call GenspaceCurve(Mesh)
       !      call PrintCurve(Mesh) 
    else
       ! find the smallest ne2 which is a power of 2 and ne2>ne
       ne2=2**ceiling( log(real(ne))/log(2d0) )
       if (ne2<ne) call abortmp('Fatel SFC error')

       allocate(Mesh2(ne2,ne2))
       allocate(Mesh2_map(ne2,ne2,2))
       allocate(sfcij(0:ne2*ne2,2))

       call GenspaceCurve(Mesh2)  ! SFC partition for ne2

       ! associate every element on the ne x ne mesh (Mesh)
       ! with its closest element on the ne2 x ne2 mesh (Mesh2)
       ! Store this as a map from Mesh2 -> Mesh in Mesh2_map.
       ! elements in Mesh2 which are not mapped get assigned a value of 0
       Mesh2_map=0
       do j=1,ne
          do i=1,ne
             ! map this element to an (i2,j2) element
             ! [ (i-.5)/ne , (j-.5)/ne ]  = [ (i2-.5)/ne2 , (j2-.5)/ne2 ]
             i2=nint( ((i-.5)/ne)*ne2 + .5 )
             j2=nint( ((j-.5)/ne)*ne2 + .5 )
             if (i2<1) i2=1
             if (i2>ne2) i2=ne2
             if (j2<1) j2=1
             if (j2>ne2) j2=ne2
             Mesh2_map(i2,j2,1)=i
             Mesh2_map(i2,j2,2)=j
          enddo
       enddo

       ! create a reverse index array for Mesh2
       ! k = Mesh2(i,j) 
       ! (i,j) = (sfcij(k,1),sfci(k,2)) 
       do j=1,ne2
          do i=1,ne2
             k=Mesh2(i,j)
             sfcij(k,1)=i
             sfcij(k,2)=j
          enddo
       enddo

       ! generate a SFC for Mesh with the same ordering as the 
       ! elements in Mesh2 which map to Mesh.
       sfc_index=0
       do k=0,ne2*ne2-1
          i2=sfcij(k,1)
          j2=sfcij(k,2)
          i=Mesh2_map(i2,j2,1)
          j=Mesh2_map(i2,j2,2)
          if (i/=0) then
             ! (i2,j2) element maps to (i,j) element
             Mesh(i,j)=sfc_index
             sfc_index=sfc_index+1
          endif
       enddo
#if 0
       print *,'SFC Mapping to non powers of 2,3 used.  Mesh:'  
       do j=1,ne
          write(*,'(99i3)') (Mesh(i,j),i=1,ne)
       enddo
       call PrintCurve(Mesh2) 
#endif
       deallocate(Mesh2)
       deallocate(Mesh2_map)
       deallocate(sfcij)
    endif


    ! -------------------------------------------
    !  Setup the space-filling curve for face 1
    ! -------------------------------------------
    offset=0
    do j=1,ne
       do i=1,ne
          GridElem(i,j,1)%SpaceCurve = offset + Mesh(i,ne-j+1)
       enddo
    enddo

    ! -------------------------------------------
    !  Setup the space-filling curve for face 2
    ! -------------------------------------------
    offset = offset + ne*ne
    do j=1,ne
       do i=1,ne
          GridElem(i,j,2)%SpaceCurve = offset + Mesh(i,ne-j+1)
       enddo
    enddo

    ! -------------------------------------------
    !  Setup the space-filling curve for face 6
    ! -------------------------------------------
    offset = offset + ne*ne
    do j=1,ne
       do i=1,ne
          GridElem(i,j,6)%SpaceCurve = offset + Mesh(ne-i+1,ne-j+1)
       enddo
    enddo

    ! -------------------------------------------
    !  Setup the space-filling curve for face 4
    ! -------------------------------------------
    offset = offset + ne*ne
    do j=1,ne
       do i=1,ne
          GridElem(i,j,4)%SpaceCurve = offset + Mesh(ne-j+1,i)
       enddo
    enddo

    ! -------------------------------------------
    !  Setup the space-filling curve for face 5
    ! -------------------------------------------
    offset = offset + ne*ne
    do j=1,ne
       do i=1,ne
          GridElem(i,j,5)%SpaceCurve = offset + Mesh(i,j)
       enddo
    enddo


    ! -------------------------------------------
    !  Setup the space-filling curve for face 3
    ! -------------------------------------------
    offset = offset + ne*ne
    do j=1,ne
       do i=1,ne
          GridElem(i,j,3)%SpaceCurve = offset + Mesh(i,j)
       enddo
    enddo

    ! ==================
    ! face interiors
    ! ==================
    do k=1,6
       ! setup  SOUTH, WEST, SW neighbors
       do j=2,ne
          do i=2,ne
             ALLOCATE(GridElem(i,j,k)%nbrs(west)%n(1))
             ALLOCATE(GridElem(i,j,k)%nbrs(south)%n(1))
             ALLOCATE(GridElem(i,j,k)%nbrs(swest)%n(1))
             ALLOCATE(GridElem(i,j,k)%nbrs(west)%f(1))
             ALLOCATE(GridElem(i,j,k)%nbrs(south)%f(1))
             ALLOCATE(GridElem(i,j,k)%nbrs(swest)%f(1))
             GridElem(i,j,k)%nbrs(west)%n(1)  = GridElem(i-1,j,k)%number
             GridElem(i,j,k)%nbrs(west)%f(1)  = k
             GridElem(i,j,k)%wgtP(west)       = EdgeWgtP
             GridElem(i,j,k)%nbrs(south)%n(1) = GridElem(i,j-1,k)%number
             GridElem(i,j,k)%nbrs(south)%f(1) = k
             GridElem(i,j,k)%wgtP(south)      = EdgeWgtP
             GridElem(i,j,k)%nbrs(swest)%n(1) = GridElem(i-1,j-1,k)%number
             GridElem(i,j,k)%nbrs(swest)%f(1) = k
             GridElem(i,j,k)%wgtP(swest)      = CornerWgt
          end do
       end do

       !  setup EAST, NORTH, NE neighbors
       do j=1,ne-1
          do i=1,ne-1
             ALLOCATE(GridElem(i,j,k)%nbrs(east)%n(1))
             ALLOCATE(GridElem(i,j,k)%nbrs(north)%n(1))
             ALLOCATE(GridElem(i,j,k)%nbrs(neast)%n(1))
             ALLOCATE(GridElem(i,j,k)%nbrs(east)%f(1))
             ALLOCATE(GridElem(i,j,k)%nbrs(north)%f(1))
             ALLOCATE(GridElem(i,j,k)%nbrs(neast)%f(1))
             GridElem(i,j,k)%nbrs(east)%n(1)   = GridElem(i+1,j,k)%number
             GridElem(i,j,k)%nbrs(east)%f(1)   = k
             GridElem(i,j,k)%wgtP(east)        = EdgeWgtP
             GridElem(i,j,k)%nbrs(north)%n(1)  = GridElem(i,j+1,k)%number
             GridElem(i,j,k)%nbrs(north)%f(1)  = k
             GridElem(i,j,k)%wgtP(north)       = EdgeWgtP
             GridElem(i,j,k)%nbrs(neast)%n(1)  = GridElem(i+1,j+1,k)%number
             GridElem(i,j,k)%nbrs(neast)%f(1)  = k
             GridElem(i,j,k)%wgtP(neast)       = CornerWgt
          end do
       end do

       ! Setup the remaining SOUTH, EAST, and SE neighbors
       do j=2,ne
          do i=1,ne-1
             ALLOCATE(GridElem(i,j,k)%nbrs(south)%n(1))
             ALLOCATE(GridElem(i,j,k)%nbrs(east)%n(1))
             ALLOCATE(GridElem(i,j,k)%nbrs(seast)%n(1))
             ALLOCATE(GridElem(i,j,k)%nbrs(south)%f(1))
             ALLOCATE(GridElem(i,j,k)%nbrs(east)%f(1))
             ALLOCATE(GridElem(i,j,k)%nbrs(seast)%f(1))
             GridElem(i,j,k)%nbrs(south)%n(1)  = GridElem(i,j-1,k)%number
             GridElem(i,j,k)%nbrs(south)%f(1)  = k
             GridElem(i,j,k)%wgtP(south)       = EdgeWgtP
             GridElem(i,j,k)%nbrs(east)%n(1)   = GridElem(i+1,j,k)%number
             GridElem(i,j,k)%nbrs(east)%f(1)   = k
             GridElem(i,j,k)%wgtP(east)        = EdgeWgtP
             GridElem(i,j,k)%nbrs(seast)%n(1)  = GridElem(i+1,j-1,k)%number
             GridElem(i,j,k)%nbrs(seast)%f(1)  = k
             GridElem(i,j,k)%wgtP(seast)       = CornerWgt
          enddo
       enddo

       ! Setup the remaining NORTH, WEST, and NW neighbors
       do j=1,ne-1
          do i=2,ne
             ALLOCATE(GridElem(i,j,k)%nbrs(north)%n(1))
             ALLOCATE(GridElem(i,j,k)%nbrs(west)%n(1))
             ALLOCATE(GridElem(i,j,k)%nbrs(nwest)%n(1))
             ALLOCATE(GridElem(i,j,k)%nbrs(north)%f(1))
             ALLOCATE(GridElem(i,j,k)%nbrs(west)%f(1))
             ALLOCATE(GridElem(i,j,k)%nbrs(nwest)%f(1))
             GridElem(i,j,k)%nbrs(north)%n(1)  = GridElem(i,j+1,k)%number
             GridElem(i,j,k)%nbrs(north)%f(1)  = k
             GridElem(i,j,k)%wgtP(north)       = EdgeWgtP
             GridElem(i,j,k)%nbrs(west)%n(1)   = GridElem(i-1,j,k)%number
             GridElem(i,j,k)%nbrs(west)%f(1)   = k
             GridElem(i,j,k)%wgtP(west)        = EdgeWgtP
             GridElem(i,j,k)%nbrs(nwest)%n(1)  = GridElem(i-1,j+1,k)%number
             GridElem(i,j,k)%nbrs(nwest)%f(1)  = k
             GridElem(i,j,k)%wgtP(nwest)       = CornerWgt
          enddo
       enddo
    end do

    ! ======================
    ! west/east "belt" edges
    ! ======================

    do k=1,4
       do j=1,ne
          ALLOCATE(GridElem(1 ,j,k)%nbrs(west)%n(1))
          ALLOCATE(GridElem(ne,j,k)%nbrs(east)%n(1))
          ALLOCATE(GridElem(1 ,j,k)%nbrs(west)%f(1))
          ALLOCATE(GridElem(ne,j,k)%nbrs(east)%f(1))
          GridElem(1 ,j,k)%nbrs(west)%n(1) = GridElem(ne,j,MODULO(2+k,4)+1)%number
          GridElem(1 ,j,k)%nbrs(west)%f(1) = MODULO(2+k,4)+1
          GridElem(1 ,j,k)%wgtP(west)  = EdgeWgtP
          GridElem(ne,j,k)%nbrs(east)%n(1) = GridElem(1 ,j,MODULO(k  ,4)+1)%number
          GridElem(ne,j,k)%nbrs(east)%f(1) = MODULO(k  ,4)+1
          GridElem(ne,j,k)%wgtP(east)  = EdgeWgtP

          !  Special rules for corner 'edges'
          if( j /= 1) then
             ALLOCATE(GridElem(1 ,j,k)%nbrs(swest)%n(1))
             ALLOCATE(GridElem(ne,j,k)%nbrs(seast)%n(1))
             ALLOCATE(GridElem(1 ,j,k)%nbrs(swest)%f(1))
             ALLOCATE(GridElem(ne,j,k)%nbrs(seast)%f(1))
             GridElem(1 ,j,k)%nbrs(swest)%n(1)   = GridElem(ne,j-1,MODULO(2+k,4)+1)%number
             GridElem(1 ,j,k)%nbrs(swest)%f(1)   = MODULO(2+k,4)+1
             GridElem(1 ,j,k)%wgtP(swest)        = CornerWgt
             GridElem(ne,j,k)%nbrs(seast)%n(1)   = GridElem(1 ,j-1,MODULO(k  ,4)+1)%number
             GridElem(ne,j,k)%nbrs(seast)%f(1)   = MODULO(k  ,4)+1
             GridElem(ne,j,k)%wgtP(seast)        = CornerWgt
          endif
          if( j /= ne) then
             ALLOCATE(GridElem(1 ,j,k)%nbrs(nwest)%n(1))
             ALLOCATE(GridElem(ne,j,k)%nbrs(neast)%n(1))
             ALLOCATE(GridElem(1 ,j,k)%nbrs(nwest)%f(1))
             ALLOCATE(GridElem(ne,j,k)%nbrs(neast)%f(1))
             GridElem(1 ,j,k)%nbrs(nwest)%n(1)   = GridElem(ne,j+1,MODULO(2+k,4)+1)%number
             GridElem(1 ,j,k)%nbrs(nwest)%f(1)   = MODULO(2+k,4)+1
             GridElem(1 ,j,k)%wgtP(nwest)        = CornerWgt
             GridElem(ne,j,k)%nbrs(neast)%n(1)   = GridElem(1 ,j+1,MODULO(k  ,4)+1)%number
             GridElem(ne,j,k)%nbrs(neast)%f(1)   = MODULO(k  ,4)+1
             GridElem(ne,j,k)%wgtP(neast)        = CornerWgt
          endif
       end do
    end do


    ! ==================================
    ! south edge of 1 / north edge of 5
    ! ==================================

    do i=1,ne
       ALLOCATE(GridElem(i,1 ,1)%nbrs(south)%n(1))
       ALLOCATE(GridElem(i,ne,5)%nbrs(north)%n(1))
       ALLOCATE(GridElem(i,1 ,1)%nbrs(south)%f(1))
       ALLOCATE(GridElem(i,ne,5)%nbrs(north)%f(1))
       GridElem(i,1 ,1)%nbrs(south)%n(1) = GridElem(i,ne,5)%number
       GridElem(i,1 ,1)%nbrs(south)%f(1) = 5
       GridElem(i,1 ,1)%wgtP(south)      = EdgeWgtP
       GridElem(i,ne,5)%nbrs(north)%n(1) = GridElem(i,1 ,1)%number
       GridElem(i,ne,5)%nbrs(north)%f(1) = 1
       GridElem(i,ne,5)%wgtP(north)      = EdgeWgtP

       !  Special rules for corner 'edges'
       if( i /= 1) then
          ALLOCATE(GridElem(i,1 ,1)%nbrs(swest)%n(1))
          ALLOCATE(GridElem(i,ne,5)%nbrs(nwest)%n(1))
          ALLOCATE(GridElem(i,1 ,1)%nbrs(swest)%f(1))
          ALLOCATE(GridElem(i,ne,5)%nbrs(nwest)%f(1))
          GridElem(i,1 ,1)%nbrs(swest)%n(1)    = GridElem(i-1,ne,5)%number
          GridElem(i,1 ,1)%nbrs(swest)%f(1)    = 5
          GridElem(i,1 ,1)%wgtP(swest)         = CornerWgt
          GridElem(i,ne,5)%nbrs(nwest)%n(1)    = GridElem(i-1,1 ,1)%number
          GridElem(i,ne,5)%nbrs(nwest)%f(1)    = 1
          GridElem(i,ne,5)%wgtP(nwest)         = CornerWgt
       endif
       if( i /= ne) then
          ALLOCATE(GridElem(i,1 ,1)%nbrs(seast)%n(1))
          ALLOCATE(GridElem(i,ne,5)%nbrs(neast)%n(1))
          ALLOCATE(GridElem(i,1 ,1)%nbrs(seast)%f(1))
          ALLOCATE(GridElem(i,ne,5)%nbrs(neast)%f(1))
          GridElem(i,1 ,1)%nbrs(seast)%n(1)    = GridElem(i+1,ne,5)%number
          GridElem(i,1 ,1)%nbrs(seast)%f(1)    = 5
          GridElem(i,1 ,1)%wgtP(seast)         = CornerWgt
          GridElem(i,ne,5)%nbrs(neast)%n(1)    = GridElem(i+1,1 ,1)%number
          GridElem(i,ne,5)%nbrs(neast)%f(1)    = 1
          GridElem(i,ne,5)%wgtP(neast)         = CornerWgt
       endif

    end do

    ! ==================================
    ! south edge of 2 / east edge of 5
    ! ==================================

    do i=1,ne
       irev=ne+1-i
       ALLOCATE(GridElem(i, 1,2)%nbrs(south)%n(1))
       ALLOCATE(GridElem(ne,i,5)%nbrs(east)%n(1))
       ALLOCATE(GridElem(i, 1,2)%nbrs(south)%f(1))
       ALLOCATE(GridElem(ne,i,5)%nbrs(east)%f(1))
       GridElem(i,1 ,2)%nbrs(south)%n(1) = GridElem(ne,irev,5)%number
       GridElem(i,1 ,2)%nbrs(south)%f(1) = 5
       GridElem(i,1 ,2)%wgtP(south)      = EdgeWgtP
       GridElem(ne,i,5)%nbrs(east)%n(1)  = GridElem(irev,1 ,2)%number
       GridElem(ne,i,5)%nbrs(east)%f(1)  = 2
       GridElem(ne,i,5)%wgtP(east)       = EdgeWgtP

       !  Special rules for corner 'edges'
       if( i /= 1) then
          ALLOCATE(GridElem(i, 1,2)%nbrs(swest)%n(1))
          ALLOCATE(GridElem(ne,i,5)%nbrs(seast)%n(1))
          ALLOCATE(GridElem(i, 1,2)%nbrs(swest)%f(1))
          ALLOCATE(GridElem(ne,i,5)%nbrs(seast)%f(1))
          GridElem(i,1 ,2)%nbrs(swest)%n(1) = GridElem(ne,irev+1,5)%number
          GridElem(i,1 ,2)%nbrs(swest)%f(1) = 5
          GridElem(i,1 ,2)%wgtP(swest)      = CornerWgt
          GridElem(ne,i,5)%nbrs(seast)%n(1) = GridElem(irev+1,1 ,2)%number
          GridElem(ne,i,5)%nbrs(seast)%f(1) = 2
          GridElem(ne,i,5)%wgtP(seast)      = CornerWgt
       endif
       if(i /= ne) then
          ALLOCATE(GridElem(i, 1,2)%nbrs(seast)%n(1))
          ALLOCATE(GridElem(ne,i,5)%nbrs(neast)%n(1))
          ALLOCATE(GridElem(i, 1,2)%nbrs(seast)%f(1))
          ALLOCATE(GridElem(ne,i,5)%nbrs(neast)%f(1))
          GridElem(i,1 ,2)%nbrs(seast)%n(1)   = GridElem(ne,irev-1,5)%number
          GridElem(i,1 ,2)%nbrs(seast)%f(1)   = 5
          GridElem(i,1 ,2)%wgtP(seast)        = CornerWgt
          GridElem(ne,i,5)%nbrs(neast)%n(1)   = GridElem(irev-1,1 ,2)%number
          GridElem(ne,i,5)%nbrs(neast)%f(1)   = 2
          GridElem(ne,i,5)%wgtP(neast)        = CornerWgt
       endif
    enddo
    ! ==================================
    ! south edge of 3 / south edge of 5
    ! ==================================

    do i=1,ne
       irev=ne+1-i
       ALLOCATE(GridElem(i,1,3)%nbrs(south)%n(1))
       ALLOCATE(GridElem(i,1,5)%nbrs(south)%n(1))
       ALLOCATE(GridElem(i,1,3)%nbrs(south)%f(1))
       ALLOCATE(GridElem(i,1,5)%nbrs(south)%f(1))
       GridElem(i,1,3)%nbrs(south)%n(1) = GridElem(irev,1,5)%number
       GridElem(i,1,3)%nbrs(south)%f(1) = 5
       GridElem(i,1,3)%wgtP(south)      = EdgeWgtP
       GridElem(i,1,5)%nbrs(south)%n(1) = GridElem(irev,1,3)%number
       GridElem(i,1,5)%nbrs(south)%f(1) = 3
       GridElem(i,1,5)%wgtP(south)      = EdgeWgtP

       !  Special rules for corner 'edges'
       if( i /= 1) then
          ALLOCATE(GridElem(i,1,3)%nbrs(swest)%n(1))
          ALLOCATE(GridElem(i,1,5)%nbrs(swest)%n(1))
          ALLOCATE(GridElem(i,1,3)%nbrs(swest)%f(1))
          ALLOCATE(GridElem(i,1,5)%nbrs(swest)%f(1))
          GridElem(i,1,3)%nbrs(swest)%n(1) = GridElem(irev+1,1,5)%number
          GridElem(i,1,3)%nbrs(swest)%f(1) = 5
          GridElem(i,1,3)%wgtP(swest)      = CornerWgt
          GridElem(i,1,5)%nbrs(swest)%n(1) = GridElem(irev+1,1,3)%number
          GridElem(i,1,5)%nbrs(swest)%f(1) = 3
          GridElem(i,1,5)%wgtP(swest)      = CornerWgt
       endif
       if(i /= ne) then
          ALLOCATE(GridElem(i,1,3)%nbrs(seast)%n(1))
          ALLOCATE(GridElem(i,1,5)%nbrs(seast)%n(1))
          ALLOCATE(GridElem(i,1,3)%nbrs(seast)%f(1))
          ALLOCATE(GridElem(i,1,5)%nbrs(seast)%f(1))
          GridElem(i,1,3)%nbrs(seast)%n(1)    = GridElem(irev-1,1,5)%number
          GridElem(i,1,3)%nbrs(seast)%f(1)    = 5
          GridElem(i,1,3)%wgtP(seast)         = CornerWgt
          GridElem(i,1,5)%nbrs(seast)%n(1)    = GridElem(irev-1,1,3)%number
          GridElem(i,1,5)%nbrs(seast)%f(1)    = 3
          GridElem(i,1,5)%wgtP(seast)         = CornerWgt
       endif
    end do

    ! ==================================
    ! south edge of 4 / west edge of 5
    ! ==================================

    do i=1,ne
       irev=ne+1-i
       ALLOCATE(GridElem(i,1,4)%nbrs(south)%n(1))
       ALLOCATE(GridElem(1,i,5)%nbrs(west)%n(1))
       ALLOCATE(GridElem(i,1,4)%nbrs(south)%f(1))
       ALLOCATE(GridElem(1,i,5)%nbrs(west)%f(1))
       GridElem(i,1,4)%nbrs(south)%n(1) = GridElem(1,i,5)%number
       GridElem(i,1,4)%nbrs(south)%f(1) = 5
       GridElem(i,1,4)%wgtP(south)      = EdgeWgtP
       GridElem(1,i,5)%nbrs(west)%n(1)  = GridElem(i,1,4)%number
       GridElem(1,i,5)%nbrs(west)%f(1)  = 4
       GridElem(1,i,5)%wgtP(west)       = EdgeWgtP
       !  Special rules for corner 'edges'
       if( i /= 1) then
          ALLOCATE(GridElem(i,1,4)%nbrs(swest)%n(1))
          ALLOCATE(GridElem(1,i,5)%nbrs(swest)%n(1))
          ALLOCATE(GridElem(i,1,4)%nbrs(swest)%f(1))
          ALLOCATE(GridElem(1,i,5)%nbrs(swest)%f(1))
          GridElem(i,1,4)%nbrs(swest)%n(1)    = GridElem(1,i-1,5)%number
          GridElem(i,1,4)%nbrs(swest)%f(1)    = 5
          GridElem(i,1,4)%wgtP(swest)         = CornerWgt
          GridElem(1,i,5)%nbrs(swest)%n(1)    = GridElem(i-1,1,4)%number
          GridElem(1,i,5)%nbrs(swest)%f(1)    = 4
          GridElem(1,i,5)%wgtP(swest)         = CornerWgt
       endif
       if( i /= ne) then
          ALLOCATE(GridElem(i,1,4)%nbrs(seast)%n(1))
          ALLOCATE(GridElem(1,i,5)%nbrs(nwest)%n(1))
          ALLOCATE(GridElem(i,1,4)%nbrs(seast)%f(1))
          ALLOCATE(GridElem(1,i,5)%nbrs(nwest)%f(1))
          GridElem(i,1,4)%nbrs(seast)%n(1) = GridElem(1,i+1,5)%number
          GridElem(i,1,4)%nbrs(seast)%f(1) = 5
          GridElem(i,1,4)%wgtP(seast)      = CornerWgt
          GridElem(1,i,5)%nbrs(nwest)%n(1) = GridElem(i+1,1,4)%number
          GridElem(1,i,5)%nbrs(nwest)%f(1) = 4
          GridElem(1,i,5)%wgtP(nwest)      = CornerWgt
       endif
    end do

    ! ==================================
    ! north edge of 1 / south edge of 6
    ! ==================================

    do i=1,ne
       ALLOCATE(GridElem(i,ne,1)%nbrs(north)%n(1))
       ALLOCATE(GridElem(i,1 ,6)%nbrs(south)%n(1))
       ALLOCATE(GridElem(i,ne,1)%nbrs(north)%f(1))
       ALLOCATE(GridElem(i,1 ,6)%nbrs(south)%f(1))
       GridElem(i,ne,1)%nbrs(north)%n(1) = GridElem(i,1 ,6)%number
       GridElem(i,ne,1)%nbrs(north)%f(1) = 6
       GridElem(i,ne,1)%wgtP(north)      = EdgeWgtP
       GridElem(i,1 ,6)%nbrs(south)%n(1) = GridElem(i,ne,1)%number
       GridElem(i,1 ,6)%nbrs(south)%f(1) = 1
       GridElem(i,1 ,6)%wgtP(south)      = EdgeWgtP
       !  Special rules for corner 'edges'
       if( i /= 1) then
          ALLOCATE(GridElem(i,ne,1)%nbrs(nwest)%n(1))
          ALLOCATE(GridElem(i,1 ,6)%nbrs(swest)%n(1))
          ALLOCATE(GridElem(i,ne,1)%nbrs(nwest)%f(1))
          ALLOCATE(GridElem(i,1 ,6)%nbrs(swest)%f(1))
          GridElem(i,ne,1)%nbrs(nwest)%n(1) = GridElem(i-1,1 ,6)%number
          GridElem(i,ne,1)%nbrs(nwest)%f(1) = 6
          GridElem(i,ne,1)%wgtP(nwest)      = CornerWgt
          GridElem(i,1 ,6)%nbrs(swest)%n(1) = GridElem(i-1,ne,1)%number
          GridElem(i,1 ,6)%nbrs(swest)%f(1) = 1
          GridElem(i,1 ,6)%wgtP(swest)      = CornerWgt
       endif
       if( i /= ne) then
          ALLOCATE(GridElem(i,ne,1)%nbrs(neast)%n(1))
          ALLOCATE(GridElem(i,1 ,6)%nbrs(seast)%n(1))
          ALLOCATE(GridElem(i,ne,1)%nbrs(neast)%f(1))
          ALLOCATE(GridElem(i,1 ,6)%nbrs(seast)%f(1))
          GridElem(i,ne,1)%nbrs(neast)%n(1) = GridElem(i+1,1 ,6)%number
          GridElem(i,ne,1)%nbrs(neast)%f(1) = 6
          GridElem(i,ne,1)%wgtP(neast)      = CornerWgt
          GridElem(i,1 ,6)%nbrs(seast)%n(1) = GridElem(i+1,ne,1)%number
          GridElem(i,1 ,6)%nbrs(seast)%f(1) = 1
          GridElem(i,1 ,6)%wgtP(seast)      = CornerWgt
       endif
    end do

    ! ==================================
    ! north edge of 2 / east edge of 6
    ! ==================================

    do i=1,ne
       ALLOCATE(GridElem(i,ne,2)%nbrs(north)%n(1))
       ALLOCATE(GridElem(ne,i,6)%nbrs(east )%n(1))
       ALLOCATE(GridElem(i,ne,2)%nbrs(north)%f(1))
       ALLOCATE(GridElem(ne,i,6)%nbrs(east )%f(1))
       GridElem(i,ne,2)%nbrs(north)%n(1) = GridElem(ne,i,6)%number
       GridElem(i,ne,2)%nbrs(north)%f(1) = 6
       GridElem(i,ne,2)%wgtP(north)      = EdgeWgtP
       GridElem(ne,i,6)%nbrs(east)%n(1)  = GridElem(i,ne,2)%number
       GridElem(ne,i,6)%nbrs(east)%f(1)  = 2
       GridElem(ne,i,6)%wgtP(east)       = EdgeWgtP
       !  Special rules for corner 'edges'
       if( i /= 1) then
          ALLOCATE(GridElem(i,ne,2)%nbrs(nwest)%n(1))
          ALLOCATE(GridElem(ne,i,6)%nbrs(seast)%n(1))
          ALLOCATE(GridElem(i,ne,2)%nbrs(nwest)%f(1))
          ALLOCATE(GridElem(ne,i,6)%nbrs(seast)%f(1))
          GridElem(i,ne,2)%nbrs(nwest)%n(1)    = GridElem(ne,i-1,6)%number
          GridElem(i,ne,2)%nbrs(nwest)%f(1)    = 6
          GridElem(i,ne,2)%wgtP(nwest)         = CornerWgt
          GridElem(ne,i,6)%nbrs(seast)%n(1)    = GridElem(i-1,ne,2)%number
          GridElem(ne,i,6)%nbrs(seast)%f(1)    = 2
          GridElem(ne,i,6)%wgtP(seast)         = CornerWgt
       endif
       if( i /= ne) then
          ALLOCATE(GridElem(i,ne,2)%nbrs(neast)%n(1))
          ALLOCATE(GridElem(ne,i,6)%nbrs(neast)%n(1))
          ALLOCATE(GridElem(i,ne,2)%nbrs(neast)%f(1))
          ALLOCATE(GridElem(ne,i,6)%nbrs(neast)%f(1))
          GridElem(i,ne,2)%nbrs(neast)%n(1) = GridElem(ne,i+1,6)%number
          GridElem(i,ne,2)%nbrs(neast)%f(1) = 6
          GridElem(i,ne,2)%wgtP(neast)      = CornerWgt
          GridElem(ne,i,6)%nbrs(neast)%n(1) = GridElem(i+1,ne,2)%number
          GridElem(ne,i,6)%nbrs(neast)%f(1) = 2
          GridElem(ne,i,6)%wgtP(neast)      = CornerWgt
       endif
    end do

    ! ===================================
    ! north edge of 3 / north edge of 6
    ! ===================================

    do i=1,ne
       irev=ne+1-i
       ALLOCATE(GridElem(i,ne,3)%nbrs(north)%n(1))
       ALLOCATE(GridElem(i,ne,6)%nbrs(north)%n(1))
       ALLOCATE(GridElem(i,ne,3)%nbrs(north)%f(1))
       ALLOCATE(GridElem(i,ne,6)%nbrs(north)%f(1))
       GridElem(i,ne,3)%nbrs(north)%n(1) = GridElem(irev,ne,6)%number
       GridElem(i,ne,3)%nbrs(north)%f(1) = 6
       GridElem(i,ne,3)%wgtP(north)      = EdgeWgtP
       GridElem(i,ne,6)%nbrs(north)%n(1) = GridElem(irev,ne,3)%number
       GridElem(i,ne,6)%nbrs(north)%f(1) = 3
       GridElem(i,ne,6)%wgtP(north)      = EdgeWgtP
       !  Special rules for corner 'edges'
       if( i /= 1) then
          ALLOCATE(GridElem(i,ne,3)%nbrs(nwest)%n(1))
          ALLOCATE(GridElem(i,ne,6)%nbrs(nwest)%n(1))
          ALLOCATE(GridElem(i,ne,3)%nbrs(nwest)%f(1))
          ALLOCATE(GridElem(i,ne,6)%nbrs(nwest)%f(1))
          GridElem(i,ne,3)%nbrs(nwest)%n(1) = GridElem(irev+1,ne,6)%number
          GridElem(i,ne,3)%nbrs(nwest)%f(1) = 6
          GridElem(i,ne,3)%wgtP(nwest)      = CornerWgt
          GridElem(i,ne,6)%nbrs(nwest)%n(1) = GridElem(irev+1,ne,3)%number
          GridElem(i,ne,6)%nbrs(nwest)%f(1) = 3
          GridElem(i,ne,6)%wgtP(nwest)      = CornerWgt
       endif
       if( i /= ne) then
          ALLOCATE(GridElem(i,ne,3)%nbrs(neast)%n(1))
          ALLOCATE(GridElem(i,ne,6)%nbrs(neast)%n(1))
          ALLOCATE(GridElem(i,ne,3)%nbrs(neast)%f(1))
          ALLOCATE(GridElem(i,ne,6)%nbrs(neast)%f(1))
          GridElem(i,ne,3)%nbrs(neast)%n(1) = GridElem(irev-1,ne,6)%number
          GridElem(i,ne,3)%nbrs(neast)%f(1) = 6
          GridElem(i,ne,3)%wgtP(neast)      = CornerWgt
          GridElem(i,ne,6)%nbrs(neast)%n(1) = GridElem(irev-1,ne,3)%number
          GridElem(i,ne,6)%nbrs(neast)%f(1) = 3
          GridElem(i,ne,6)%wgtP(neast)      = CornerWgt
       endif
    end do

    ! ===================================
    ! north edge of 4 / west edge of 6
    ! ===================================

    do i=1,ne
       irev=ne+1-i
       ALLOCATE(GridElem(i,ne,4)%nbrs(north)%n(1))
       ALLOCATE(GridElem(1,i, 6)%nbrs(west)%n(1))
       ALLOCATE(GridElem(i,ne,4)%nbrs(north)%f(1))
       ALLOCATE(GridElem(1,i, 6)%nbrs(west)%f(1))
       GridElem(i,ne,4)%nbrs(north)%n(1) = GridElem(1,irev,6)%number
       GridElem(i,ne,4)%nbrs(north)%f(1) = 6
       GridElem(i,ne,4)%wgtP(north)      = EdgeWgtP
       GridElem(1,i,6)%nbrs(west)%n(1)   = GridElem(irev,ne,4)%number
       GridElem(1,i,6)%nbrs(west)%f(1)   = 4
       GridElem(1,i,6)%wgtP(west)        = EdgeWgtP
       !  Special rules for corner 'edges'
       if( i /= 1) then
          ALLOCATE(GridElem(i,ne,4)%nbrs(nwest)%n(1))
          ALLOCATE(GridElem(1,i, 6)%nbrs(swest)%n(1))
          ALLOCATE(GridElem(i,ne,4)%nbrs(nwest)%f(1))
          ALLOCATE(GridElem(1,i, 6)%nbrs(swest)%f(1))
          GridElem(i,ne,4)%nbrs(nwest)%n(1) = GridElem(1,irev+1,6)%number
          GridElem(i,ne,4)%nbrs(nwest)%f(1) = 6
          GridElem(i,ne,4)%wgtP(nwest)      = CornerWgt
          GridElem(1,i,6)%nbrs(swest)%n(1)  = GridElem(irev+1,ne,4)%number
          GridElem(1,i,6)%nbrs(swest)%f(1)  = 4
          GridElem(1,i,6)%wgtP(swest)       = CornerWgt
       endif
       if( i /= ne) then
          ALLOCATE(GridElem(i,ne,4)%nbrs(neast)%n(1))
          ALLOCATE(GridElem(1,i ,6)%nbrs(nwest)%n(1))
          ALLOCATE(GridElem(i,ne,4)%nbrs(neast)%f(1))
          ALLOCATE(GridElem(1,i ,6)%nbrs(nwest)%f(1))
          GridElem(i,ne,4)%nbrs(neast)%n(1) = GridElem(1,irev-1,6)%number
          GridElem(i,ne,4)%nbrs(neast)%f(1) = 6
          GridElem(i,ne,4)%wgtP(neast)      = CornerWgt
          GridElem(1,i,6)%nbrs(nwest)%n(1)  = GridElem(irev-1,ne,4)%number
          GridElem(1,i,6)%nbrs(nwest)%f(1)  = 4
          GridElem(1,i,6)%wgtP(nwest)       = CornerWgt
       endif
    end do
    

    ielem = 1                       ! Element counter
    do k=1,6
       do j=1,ne
          do i=1,ne
             do l=1,8
                if (ASSOCIATED(GridElem(i,j,k)%nbrs(l)%n)) then
                   ALLOCATE(GridVertex(ielem)%nbrs(l)%n(SIZE(GridElem(i,j,k)%nbrs(l)%n)))
                   ALLOCATE(GridVertex(ielem)%nbrs(l)%f(SIZE(GridElem(i,j,k)%nbrs(l)%n)))
                   GridVertex(ielem)%nbrs(l)%n(:)       = GridElem(i,j,k)%nbrs(l)%n(:)
                   GridVertex(ielem)%nbrs(l)%f(:)       = GridElem(i,j,k)%nbrs(l)%f(:)
                end if
             end do
             GridVertex(ielem)%wgtP       = GridElem(i,j,k)%wgtP
             GridVertex(ielem)%wgtP_ghost = GridElem(i,j,k)%wgtP_ghost
             GridVertex(ielem)%number     = GridElem(i,j,k)%number
             GridVertex(ielem)%processor_number  = 0
             GridVertex(ielem)%SpaceCurve = GridElem(i,j,k)%SpaceCurve
             ielem=ielem+1
          end do
       end do
    end do

    DEALLOCATE(Mesh)
    DEALLOCATE(GridElem)
#if 0
    if(OutputFiles) then
       close(7)
       close(8)
    endif
#endif

    ! =======================================
    ! Generate cube graph...
    ! =======================================

#if 0
    if(OutputFiles) then
       write(9,*)nelem,2*nelem      ! METIS requires this first line
    endif
#endif

    ! ============================================
    !  Setup the Grid edges (topology independent)
    ! ============================================
    call initgridedge(GridEdge,GridVertex)

    ! ============================================
    !  Setup the Grid edge Indirect addresses
    !          (topology dependent)
    ! ============================================
    nedge = SIZE(GridEdge)
    do i=1,nedge
       call CubeSetupEdgeIndex(GridEdge(i))
    enddo

  end subroutine CubeTopology

  ! =======================================
  ! cube_assemble:
  !
  ! Assemble the cube field element by element
  ! this routine is assumed to be single 
  ! threaded...
  ! =======================================

  function cube_assemble(gbl,fld,elem,par,nelemd,nelem,ielem) result(ierr)
    use element_mod, only : element_t
    use parallel_mod, only : abortmp
#ifdef _MPI
    use parallel_mod, only : parallel_t, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_STATUS_SIZE, MPI_REAL8,MPI_TAG
#else
    use parallel_mod, only : parallel_t
#endif
    real (kind=real_kind) :: gbl(:,:,:,:)    ! global output field 
    real (kind=real_kind) :: fld(:,:,:)      ! local model field  
    type (element_t)      :: elem            ! element to assemble 
    type (parallel_t)     :: par             ! parallel structure 
    integer               :: nelemd          ! number of elements on the node
    integer               :: nelem           ! number of elements on the node
    integer               :: ielem           ! local element ctr 
    integer               :: ierr            ! returned error code

    ! Local variables

    integer :: ie,je,face_no
    integer :: ibase,jbase
    integer :: i,j,k
    integer :: elem_number

    integer :: ne1,ne2    ! element dimensions
    integer :: n1,n2      ! gbl face dimensions
    integer :: nface      ! number of faces (must be 6)
    integer :: nlyr       ! number of layers

#if defined(_MPI)
    integer :: ectr       ! global element counter
    integer tag
    integer :: count      ! w/o "::", triggers PGI 3.1 F90 bug 
    integer pe
    integer status(MPI_STATUS_SIZE)
    integer mpi_err
#endif      

    call abortmp('Because convert_gbl_index is not used cube_assemble is broken. ')
    ne1   = SIZE(fld,1)
    ne2   = SIZE(fld,2)
    nlyr  = SIZE(fld,3)

    n1    = SIZE(gbl,1)
    n2    = SIZE(gbl,2)
    nface = SIZE(gbl,3)

    ! =========================
    ! Enforce certain rules...
    ! =========================

    ierr=0

    if (MODULO(n1,ne1) /= 0) then
       ierr=-1
       return
    end if

    if (MODULO(n2,ne2) /= 0) then 
       ierr=-2
       return
    end if

    if (nface /= 6) then
       ierr=-3
       return
    end if

    ! =========================================================
    ! Perform global assembly procedure element by element ...
    ! =========================================================

    if (par%rank==par%root) then

       if (ielem<=nelemd) then
          elem_number = elem%vertex%number

          call convert_gbl_index(elem_number,ie,je,face_no)
          if (face_no /= elem%vertex%face_number) call abortmp('Error in getting face number')

          ibase=ie*ne1
          jbase=je*ne2

          do k=1,nlyr
             do j=1,ne2
                do i=1,ne1
                   gbl(i+ibase,j+jbase,face_no,k)=fld(i,j,k)
                end do
             end do
          end do
       end if

#if defined(_MPI)
       if (ielem==nelemd) then
          ectr=nelemd
          do while(ectr<nelem)
             pe    = MPI_ANY_SOURCE
             tag   = MPI_ANY_TAG
             count = ne1*ne2*nlyr
             call MPI_RECV(fld(1,1,1),   &
                  count,        &
                  MPI_REAL8,    &
                  pe,           &
                  tag,          &  
                  par%comm,     &
                  status,       &
                  mpi_err) 

             elem_number = status(MPI_TAG)
             ! call convert_gbl_index(elem_number,ie,je,face_no)
             call abortmp('Because convert_gbl_index is not used for neghbors, the _MPI version needs to be fixed')

             ibase=ie*ne1
             jbase=je*ne2

             do k=1,nlyr
                do j=1,ne2
                   do i=1,ne1
                      gbl(i+ibase,j+jbase,face_no,k)=fld(i,j,k)
                   end do
                end do
             end do

             ectr=ectr+1
          end do
       end if

    else

       pe    = par%root
       tag   = elem%vertex%number
       count = ne1*ne2*nlyr
       call MPI_SEND(fld(1,1,1),    &
            count,         &
            MPI_REAL8,     &
            pe,            &
            tag,           &
            par%comm,      &
            mpi_err)
#endif
    end if

  end function cube_assemble

  ! ===================================================================
  ! CubeEdgeCount:
  !
  !  Determine the number of Grid Edges
  !
  ! ===================================================================

  function CubeEdgeCount()  result(nedge)
    use dimensions_mod, only     : ne
    use parallel_mod, only       : abortmp
    implicit none
    integer                     :: nedge

    if (0==ne) call abortmp('Error in CubeEdgeCount: ne is zero')
    nedge = nfaces*(ne*ne*nInnerElemEdge - nCornerElemEdge)

  end function CubeEdgeCount

  ! ===================================================================
  ! CubeElemCount:
  !
  !  Determine the number of Grid Elem
  !
  ! ===================================================================

  function CubeElemCount()  result(nelem)

    use dimensions_mod, only     : ne
    use parallel_mod, only       : abortmp

    implicit none
    integer                     :: nelem
    if (0==ne) call abortmp('Error in CubeElemCount: ne is zero')

    nelem = nfaces*ne*ne
  end function CubeElemCount

  subroutine CubeSetupEdgeIndex(Edge)
    use gridgraph_mod, only : gridedge_t
    use dimensions_mod, only : np
    use control_mod, only : north, south, east, west, neast, seast, swest, nwest
    type (GridEdge_t),target           :: Edge

    integer                            :: np0,sFace,dFace
    logical                            :: reverse
    integer,allocatable                :: forwardV(:), forwardP(:)
    integer,allocatable                :: backwardV(:), backwardP(:)
    integer                            :: i,ii

    ii=Edge%tail_face
    np0 = Edge%tail%wgtP(ii)

#ifdef TESTGRID

    allocate(forwardP(np0))
    allocate(backwardP(np0))

    do i=1,np0
       forwardP(i)  = i
       backwardP(i) = np0-i+1
    enddo
#endif

    sFace = Edge%tail_face
    dFace = Edge%head_face
    ! Do not reverse the indices
    reverse=.FALSE.

    ! Under special conditions use index reversal
    if(       (SFace == south .AND. dFace == east)  &
         .OR. (sFace == east  .AND. dFace == south) &
         .OR. (sFace == north .AND. dFace == west)  &
         .OR. (sFace == west  .AND. dFace == north) &
         .OR. (sFace == south .AND. dFace == south) &
         .OR. (sFace == north .AND. dFace == north) &
         .OR. (sFace == east  .AND. dFace == east ) &
         .OR. (sFace == west  .AND. dFace == west ) ) then
       reverse=.TRUE.
       Edge%reverse=.TRUE.
    endif

#ifdef TESTGRID
    !  Setup the destination indices
    select case(dFace)
    case(east)
       Edge%HeadIndex%ixV=nv
       Edge%HeadIndex%iyV=forwardV

       Edge%HeadIndex%ixP=np
       Edge%HeadIndex%iyP=forwardP
    case(west)
       Edge%HeadIndex%ixV=1
       Edge%HeadIndex%iyV=forwardV

       Edge%HeadIndex%ixP=1
       Edge%HeadIndex%iyP=forwardP
    case(north)
       Edge%HeadIndex%ixV=forwardV
       Edge%HeadIndex%iyV=nv

       Edge%HeadIndex%ixP=forwardP
       Edge%HeadIndex%iyP=np
    case(south)
       Edge%HeadIndex%ixV=forwardV
       Edge%HeadIndex%iyV=1

       Edge%HeadIndex%ixP=forwardP
       Edge%HeadIndex%iyP=1
    case(swest)
       Edge%HeadIndex%ixV=1
       Edge%HeadIndex%iyV=1

       Edge%HeadIndex%ixP=1
       Edge%HeadIndex%iyP=1
    case(seast)
       Edge%HeadIndex%ixV=nv
       Edge%HeadIndex%iyV=1

       Edge%HeadIndex%ixP=np
       Edge%HeadIndex%iyP=1
    case(nwest)
       Edge%HeadIndex%ixV=1
       Edge%HeadIndex%iyV=nv

       Edge%HeadIndex%ixP=1
       Edge%HeadIndex%iyP=np
    case(neast)
       Edge%HeadIndex%ixV=nv
       Edge%HeadIndex%iyV=nv

       Edge%HeadIndex%ixP=np
       Edge%HeadIndex%iyP=np
    case default
       write (*,*) 'SetupEdgeIndex: Error in dFace select statement'
    end select

    ! Setup the source indices
    select case(sFace)
    case(north)
       Edge%TailIndex%ixV=forwardV
       if(reverse) Edge%TailIndex%ixV=backwardV
       Edge%TailIndex%iyV=nv

       Edge%TailIndex%ixP=forwardP
       if(reverse) Edge%TailIndex%ixP=backwardP
       Edge%TailIndex%iyP=np
    case(south)
       Edge%TailIndex%ixV=forwardV
       if(reverse) Edge%TailIndex%ixV=backwardV
       Edge%TailIndex%iyV=1

       Edge%TailIndex%ixP=forwardP
       if(reverse) Edge%TailIndex%ixP=backwardP
       Edge%TailIndex%iyP=1
    case(east)
       Edge%TailIndex%ixV=nv
       Edge%TailIndex%iyV=forwardV
       if(reverse) Edge%TailIndex%iyV=backwardV

       Edge%TailIndex%ixP=np
       Edge%TailIndex%iyP=forwardP
       if(reverse) Edge%TailIndex%iyP=backwardP
    case(west)
       Edge%TailIndex%ixV=1
       Edge%TailIndex%iyV=forwardV
       if(reverse) Edge%TailIndex%iyV=backwardV

       Edge%TailIndex%ixP=1
       Edge%TailIndex%iyP=forwardP
       if(reverse) Edge%TailIndex%iyP=backwardP
    case(swest)
       Edge%TailIndex%ixV=1
       Edge%TailIndex%iyV=1

       Edge%TailIndex%ixP=1
       Edge%TailIndex%iyP=1
    case(seast)
       Edge%TailIndex%ixV=nv
       Edge%TailIndex%iyV=1

       Edge%TailIndex%ixP=np
       Edge%TailIndex%iyP=1
    case(nwest)
       Edge%TailIndex%ixV=1
       Edge%TailIndex%iyV=nv

       Edge%TailIndex%ixP=1
       Edge%TailIndex%iyP=np
    case(neast)
       Edge%TailIndex%ixV=nv
       Edge%TailIndex%iyV=nv

       Edge%TailIndex%ixP=np
       Edge%TailIndex%iyP=np
    case default
       write (*,*) 'SetupEdgeIndex: Error in sFace select statement'
    end select

    deallocate(forwardV)
    deallocate(forwardP)
    deallocate(backwardV)
    deallocate(backwardP)
#endif

  end subroutine CubeSetupEdgeIndex

  subroutine GetLatticeSpacing(spherev,dxv,spherep,dxp)
    use physical_constants, only : rearth
    use dimensions_mod, only : np

    type (spherical_polar_t), intent(in) :: spherev(np,np)
    real (kind=real_kind)                :: dxv
    type (spherical_polar_t), intent(in) :: spherep(np,np)
    real (kind=real_kind)                :: dxp

    real (kind=real_kind) xcorner,ycorner,zcorner
    real (kind=real_kind) x,y,z
    real (kind=real_kind) chord
    real (kind=real_kind) theta

    xcorner=COS(spherev(1,1)%lat)*COS(spherev(1,1)%lon)
    ycorner=COS(spherev(1,1)%lat)*SIN(spherev(1,1)%lon)
    zcorner=SIN(spherev(1,1)%lat)

    x=COS(spherev(2,1)%lat)*COS(spherev(2,1)%lon)
    y=COS(spherev(2,1)%lat)*SIN(spherev(2,1)%lon)
    z=SIN(spherev(2,1)%lat)

    chord = SQRT( (xcorner-x)**2 + &
         (ycorner-y)**2 + &
         (zcorner-z)**2 )

    theta = 2.0D0*ASIN(0.50D0*chord)

    dxv   = theta*rearth

    x=COS(spherep(1,1)%lat)*COS(spherep(1,1)%lon)
    y=COS(spherep(1,1)%lat)*SIN(spherep(1,1)%lon)
    z=SIN(spherep(1,1)%lat)

    chord = SQRT( (xcorner-x)**2 + &
         (ycorner-y)**2 + &
         (zcorner-z)**2 )

    theta = 2.0D0*ASIN(0.50D0*chord)

    dxp   = theta*rearth

  end subroutine GetLatticeSpacing

end module cube_mod



