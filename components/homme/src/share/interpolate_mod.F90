!---------------------------------------------------------------------------------------
! Interpolate_mod:
!
! 07/2016: O. Guba Changing interpolate_vector routine: 
!    (1) Instead of interpolating velocity in contravariant bases which are not
!    continuous across elements' edges, use interpolation in Cartesian basis.
!    (2) Affected routine is interpolate_vector3d. Its other version for 
!    velocities on 1 level only, interpolate_vector2d, is removed since it is 
!    not used.   
!    (3) Removing npts parameter from input params of interpolate_vector3d
!    because npts is always equal to np.
!
!

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module interpolate_mod
  use kinds, only : real_kind, iulog
  use element_mod, only : element_t
  use dimensions_mod, only : np, ne, nelemd
  use quadrature_mod, only : quadrature_t, legendre, quad_norm
  use coordinate_systems_mod, only : spherical_polar_t, cartesian2d_t, &
       cartesian3D_t, sphere2cubedsphere, spherical_to_cart, &
       cubedsphere2cart, distance, change_coordinates, projectpoint
  use physical_constants,     only : DD_PI
  use quadrature_mod,         only : quadrature_t, gauss, gausslobatto
  use parallel_mod,           only : abortmp, syncmp, parallel_t, MPIreal_t, MPIinteger_t
#ifdef _MPI
  use parallel_mod,           only : MPI_MAX, MPI_SUM, MPI_MIN
#endif
  use cube_mod,               only : convert_gbl_index, dmap, ref2sphere
  use mesh_mod,               only : MeshUseMeshFile
  use control_mod,            only : cubed_sphere_map

  implicit none
  private
  save

  logical   :: debug=.false.
#ifndef CAM
  integer, parameter, public :: MAX_VECVARS=25
  character(len=10), public :: vector_uvars(MAX_VECVARS), vector_vvars(MAX_VECVARS)
  logical, public :: replace_vec_by_vordiv(MAX_VECVARS)
#endif
! ^ ifndef CAM

  type, public :: interpolate_t
     real (kind=real_kind), dimension(:,:), pointer :: Imat  ! P_k(xj)*wj/gamma(k)
     real (kind=real_kind), dimension(:)  , pointer :: rk    ! 1/k
     real (kind=real_kind), dimension(:)  , pointer :: vtemp ! temp results
     real (kind=real_kind), dimension(:)  , pointer :: glp   ! GLL pts (nair)
  end type interpolate_t

  type, public :: interpdata_t
     ! Output Interpolation points.  Used to output data on lat-lon (or other grid)
     ! with native element interpolation.  Each element keeps a list of points from the
     ! interpolation grid that are in this element
     type (cartesian2D_t),pointer,dimension(:):: interp_xy      ! element coordinate
     integer, pointer,dimension(:)            :: ilat,ilon   ! position of interpolation point in lat-lon grid
     integer                                  :: n_interp
     integer                                  :: nlat
     integer                                  :: nlon
     logical                                  :: first_entry = .TRUE.
  end type interpdata_t

  real (kind=real_kind), private :: delta  = 1.0D-9  ! move tiny bit off center to
  ! avoid landing on element edges


  ! static data for interp_tracers
  logical                           :: interp_tracers_init=.false.
  real (kind=real_kind      )       :: interp_c(np,np)
  real (kind=real_kind      )       :: interp_gll(np)



  public :: interp_init
  public :: setup_latlon_interp
  public :: interpolate_scalar
  public :: interpolate_ce
  
  public :: interpolate_vector
  public :: set_interp_parameter
  public :: get_interp_parameter
  public :: get_interp_gweight
  public :: get_interp_lat
  public :: get_interp_lon
#ifndef CAM
  public :: var_is_vector_uvar, var_is_vector_vvar
#endif
  public :: cube_facepoint_ne
  public :: cube_facepoint_unstructured
  public :: parametric_coordinates

  public :: interpolate_tracers
  public :: interpolate_tracers_init
  public :: minmax_tracers
  public :: interpolate_2d
  public :: interpolate_create
  public :: point_inside_quad



  interface interpolate_scalar
     module procedure interpolate_scalar2d
     module procedure interpolate_scalar3d
  end interface
  interface interpolate_vector
     module procedure interpolate_vector3d
  end interface

  type (interpolate_t), target ::  interp_p

  ! store the  lat-lon grid
  ! gridtype = 1       equally spaced, including poles (FV scalars output grid)
  ! gridtype = 2       Gauss grid (CAM Eulerian)
  ! gridtype = 3       equally spaced, no poles (FV staggered velocity)
  ! Seven possible history files, last one is inithist and should be native grid
#ifndef CAM
  logical, public :: interpolate_analysis(8) = (/.true.,.false.,.false.,.false.,.false.,.false.,.false.,.false./)
#endif
  integer :: nlat,nlon
  real (kind=real_kind), pointer, public   :: lat(:)     => NULL()
  real (kind=real_kind), pointer, public   :: lon(:)     => NULL()
  real (kind=real_kind), pointer, public   :: gweight(:) => NULL()

  integer :: gridtype = 1        !
  integer :: itype = 1           ! 0 = native high order
                                 ! 1 = bilinear

  integer :: auto_grid = 0        ! 0 = interpolation grid set by namelist
                                  ! 1 = grid set via mesh resolution


contains


  subroutine set_interp_parameter(parm_name, value)
    character*(*), intent(in) :: parm_name
    character(len=80) :: msg
    integer :: value,power
    real (kind=real_kind) :: value_target

    if(parm_name .eq. 'itype') then
       itype=value
    else if(parm_name .eq. 'nlon') then
       nlon=value
    else if(parm_name .eq. 'nlat') then
       nlat=value
    else if(parm_name.eq. 'gridtype') then
       gridtype=value
    else if(parm_name.eq. 'auto') then
       auto_grid=1
       ! compute recommended nlat,nlon which has slightly higher
       ! resolution than the specifed number of points around equator given in "value"
       ! computed recommended lat-lon grid.
       ! nlon > peq   peq = points around equator cubed sphere grid
       ! take nlon power of 2, and at most 1 power of 3
       if (value.eq.0) then
           ! If reading in unstructured mesh, ne = 0
           ! This makes it hard to guess how many interpolation points to use
           ! So We'll set the default as 720 x 360
           ! BUT if you're running with an unstructured mesh, set interp_nlon and interp_nlat
           nlon = 1536
           nlat = 768
       else
           value_target=value*1.25
           power = nint(.5 +  log( value_target)/log(2d0) )
           power = max(power,7) ! min grid: 64x128
           if ( 3*2**(power-2) > value_target) then
               nlon=3*2**(power-2)   ! use 1 power of 3
           else
               nlon=2**power
           endif
       endif
       nlat=nlon/2
       if (gridtype==1) nlat=nlat+1
    else
       write(msg,*) 'Did not recognize parameter named ',parm_name,' in interpolate_mod:set_interp_parameter'
       call abortmp(msg)
    end if
  end subroutine set_interp_parameter
  function get_interp_parameter(parm_name) result(value)
    character*(*), intent(in) :: parm_name
    integer :: value
    character(len=80) :: msg
    if(parm_name .eq. 'itype') then
       value=itype
    else if(parm_name .eq. 'nlon') then
       value=nlon
    else if(parm_name .eq. 'nlat') then
       value=nlat
    else if(parm_name.eq. 'gridtype') then
       value=gridtype
    else if(parm_name.eq. 'auto_grid') then
       value=auto_grid
    else
       write(msg,*) 'Did not recognize parameter named ',parm_name,' in interpolate_mod:get_interp_parameter'
       value=-1
       call abortmp(msg)
    end if
    return
  end function get_interp_parameter
  function get_interp_gweight() result(gw)
    real(kind=real_kind) :: gw(nlat)
    gw=gweight
    return
  end function get_interp_gweight
  function get_interp_lat() result(thislat)
    use physical_constants, only : DD_PI
    real(kind=real_kind) :: thislat(nlat)
    thislat=lat*180.0D0/DD_PI
    return
  end function get_interp_lat
  function get_interp_lon() result(thislon)
    use physical_constants, only : DD_PI
    real(kind=real_kind) :: thislon(nlon)
    thislon=lon*180.0D0/DD_PI
    return
  end function get_interp_lon

  subroutine interpolate_create(gquad,interp)
    type (quadrature_t) , intent(in)   :: gquad
    type (interpolate_t), intent(out)  :: interp


    ! Local variables

    integer k,j
    integer npts
    real (kind=real_kind), dimension(:), allocatable :: gamma
    real (kind=real_kind), dimension(:), allocatable :: leg

    npts = size(gquad%points)

    allocate(interp%Imat(npts,npts))
    allocate(interp%rk(npts))
    allocate(interp%vtemp(npts))
    allocate(interp%glp(npts))
    allocate(gamma(npts))
    allocate(leg(npts))

    gamma = quad_norm(gquad,npts)

    do k=1,npts
       interp%rk(k) = 1.0D0/k
       interp%glp(k) = gquad%points(k)    !nair
    end do

    do j=1,npts
       leg=legendre(gquad%points(j),npts-1)
       do k=1,npts
          interp%Imat(j,k)=leg(k)*gquad%weights(j)/gamma(k)
       end do
    end do

    deallocate(gamma)
    deallocate(leg)

  end subroutine interpolate_create


  subroutine interpolate_tracers_init()
    use kinds,          only : longdouble_kind
    use dimensions_mod, only : np, qsize
    use quadrature_mod, only : quadrature_t, gausslobatto


    implicit none

    type (quadrature_t        )       :: gll        
    real (kind=real_kind      )       :: dp    (np)
    integer                           :: i,j

      gll=gausslobatto(np)
      dp = 1
      do i=1,np
        do j=1,np
          if (i /= j) then
            dp(i) = dp(i) * (gll%points(i) - gll%points(j))
          end if
        end do
      end do 
      do i=1,np
        do j=1,np
          interp_c(i,j) = 1/(dp(i)*dp(j))
        end do
      end do 
      interp_gll(:) = gll%points(:)
      interp_tracers_init = .true.

      deallocate(gll%points)
      deallocate(gll%weights)


  end subroutine interpolate_tracers_init




  subroutine interpolate_tracers(r, tracers, f)
    use kinds,          only : longdouble_kind
    use dimensions_mod, only : np, qsize


    implicit none
    type (cartesian2D_t), intent(in)  :: r
    real (kind=real_kind),intent(in)  :: tracers(np*np,qsize)
    real (kind=real_kind),intent(out) :: f(qsize)

    real (kind=real_kind      )       :: x     (np)
    real (kind=real_kind      )       :: y     (np)
    real (kind=real_kind      )       :: xy    (np*np)

    integer                           :: i,j

    
    if (.not. interp_tracers_init   ) then
       stop 'ERROR: interpolate_tracers() was not initialized'
    endif

    x = 1
    y = 1
    do i=1,np
      do j=1,np
        if (i /= j) then
          x(i) = x(i) * (r%x - interp_gll(j))
          y(i) = y(i) * (r%y - interp_gll(j))
        end if
      end do
    end do 

    do j=1,np  
      do i=1,np
        xy(i + (j-1)*np) = x(i)*y(j)*interp_c(i,j)
      end do
    end do 
    f = MATMUL(xy,tracers)
  end subroutine interpolate_tracers



  function linear_interpolate_2d(x,y,s) result(v)
    use dimensions_mod, only : np, qsize
    use kinds, only : longdouble_kind

    implicit none
    real (kind=longdouble_kind),intent(in)  :: x(np)
    real (kind=real_kind),intent(in)  :: y(np,np,qsize)
    type (cartesian2D_t), intent(in)  :: s        

    integer                           :: i,j,q
    real (kind=real_kind)  dx, dy(qsize), dydx(qsize), v(qsize)
    real (kind=real_kind)  y0(qsize), y1(qsize)
    type (cartesian2D_t)              :: r
 
    r = s
    if (r%x < -1) r%x = -1
    if (r%y < -1) r%y = -1
    if ( 1 < r%x) r%x =  1
    if ( 1 < r%y) r%y =  1
    do i=1,np  
      if (r%x < x(i)) exit
    end do 
    do j=1,np  
      if (r%y < x(j)) exit
    end do 
    if (1 < i) i = i-1
    if (1 < j) j = j-1
    if (np==i) i = i-1
    if (np==j) j = j-1

    dx = x(i+1)     - x(i)
    dy = y(i+1,j,:) - y(i,j,:)
    dydx = dy/dx
    y0 = y(i,j,:) + (r%x-x(i))*dydx 
    
    dy = y(i+1,j+1,:) - y(i,j+1,:)
    dydx = dy/dx
    y1 = y(i,j+1,:) + (r%x-x(i))*dydx 

    dx = x(j+1)     - x(j)
    dy = y1         - y0          
    dydx = dy/dx
    v  = y0         + (r%y-x(j))*dydx 

  end function linear_interpolate_2d

  subroutine minmax_tracers(r, tracers, mint, maxt) 
    use dimensions_mod, only : np, qsize
    use quadrature_mod, only : quadrature_t, gausslobatto


    implicit none

    type (cartesian2D_t), intent(in)  :: r
    real (kind=real_kind),intent(in)  :: tracers(np,np,qsize)
    real (kind=real_kind),intent(out) :: mint         (qsize)
    real (kind=real_kind),intent(out) :: maxt         (qsize)

    type (quadrature_t), save         :: gll        
    integer                           :: i,j
    logical            , save         :: first_time=.true.
    real (kind=real_kind)             :: y1           (qsize)
    real (kind=real_kind)             :: y2           (qsize)
    real (kind=real_kind)             :: q_interp     (4,qsize)
    type (cartesian2D_t)              :: s
    real (kind=real_kind)             :: delta
    integer :: q    

    do q=1,qsize
       mint(q) = minval(tracers(:,:,q))
       maxt(q) = maxval(tracers(:,:,q))
    enddo
    return

    delta = 1.D0/8.D0

    if (first_time) then
      first_time = .false.
      gll=gausslobatto(np)
    end if

    do i=1,np  
      if (r%x < gll%points(i)) exit
    end do 
    do j=1,np  
      if (r%y < gll%points(j)) exit
    end do 
    if (1 < i) i = i-1
    if (1 < j) j = j-1
    if (np==i) i = i-1
    if (np==j) j = j-1

!   mint(:) = minval(minval(tracers(i:i+1,j:j+1,:),1),1)
!   maxt(:) = maxval(maxval(tracers(i:i+1,j:j+1,:),1),1)

! Or check this out:
    s   = r
    s%x = s%x - delta
    s%y = s%y - delta
    q_interp(1,:) = linear_interpolate_2d(gll%points,tracers,s)
    s   = r
    s%x = s%x + delta
    s%y = s%y - delta
    q_interp(2,:) = linear_interpolate_2d(gll%points,tracers,s)
    s   = r
    s%x = s%x - delta
    s%y = s%y + delta
    q_interp(3,:) = linear_interpolate_2d(gll%points,tracers,s)
    s   = r
    s%x = s%x + delta
    s%y = s%y + delta
    q_interp(4,:) = linear_interpolate_2d(gll%points,tracers,s)

    mint(:) = minval(q_interp(:,:),1)
    maxt(:) = maxval(q_interp(:,:),1)
  end subroutine minmax_tracers

  function interpolate_2d(cart, f, interp, npts, fillvalue) result(fxy)
    integer, intent(in)               :: npts
    type (cartesian2D_t), intent(in)  :: cart
    real (kind=real_kind), intent(in) :: f(npts,npts)
    type (interpolate_t)              :: interp
    real (kind=real_kind)             :: fxy     ! value of f interpolated to (x,y)
    real (kind=real_kind), intent(in), optional :: fillvalue
    ! local variables

    real (kind=real_kind)             :: tmp_1,tmp_2
    real (kind=real_kind)             :: fk0,fk1
    real (kind=real_kind)             :: pk

    integer                           :: l,j,k

    if(present(fillvalue)) then
       if (any(f==fillvalue)) then
          fxy = fillvalue
          return
       endif
    endif


    do l=1,npts,2

       ! Compute Pk(cart%x) for Legendre order 0

       pk = 1.0D0

       fk0=0.0D0
       fk1=0.0D0
       do j=1,npts
          fk0 = fk0 + interp%Imat(j,1)*f(j,l  )
          fk1 = fk1 + interp%Imat(j,1)*f(j,l+1)
       end do
       interp%vtemp(l  ) = pk*fk0
       interp%vtemp(l+1) = pk*fk1

       ! Compute Pk(cart%x) for Legendre order 1

       tmp_2 = pk
       pk    = cart%x

       fk0=0.0D0
       fk1=0.0D0
       do j=1,npts
          fk0 = fk0 + interp%Imat(j,2)*f(j,l  )
          fk1 = fk1 + interp%Imat(j,2)*f(j,l+1)
       end do
       interp%vtemp(l  ) = interp%vtemp(l  ) + pk*fk0
       interp%vtemp(l+1) = interp%vtemp(l+1) + pk*fk1

       ! Compute Pk(cart%x) for Legendre order 2 to npts-1

       do k = 2,npts-1

          tmp_1  = tmp_2
          tmp_2  = pk
          pk = ( (2*k-1)*cart%x*tmp_2 - (k-1)*tmp_1 )*interp%rk(k)

          fk0=0.0D0
          fk1=0.0D0
          do j=1,npts
             fk0 = fk0 + interp%Imat(j,k+1)*f(j,l  )
             fk1 = fk1 + interp%Imat(j,k+1)*f(j,l+1)
          end do
          interp%vtemp(l  ) = interp%vtemp(l  ) + pk*fk0
          interp%vtemp(l+1) = interp%vtemp(l+1) + pk*fk1

       end do

    end do

    ! Compute Pk(cart%y) for Legendre order 0

    pk = 1.0

    fk0 = 0.0D0
    do j=1,npts
       fk0 = fk0 + interp%Imat(j,1)*interp%vtemp(j)
    end do
    fxy = pk*fk0

    ! Compute Pk(cart%y) for Legendre order 1

    tmp_2 = pk
    pk    = cart%y

    fk0=0.0D0
    do j=1,npts
       fk0 = fk0 + interp%Imat(j,2)*interp%vtemp(j)
    end do
    fxy = fxy + pk*fk0

    ! Compute Pk(cart%y) for Legendre order 2, npts-1

    do k = 2,npts-1
       tmp_1  = tmp_2
       tmp_2  = pk
       pk = ( (2*k-1)*cart%y*tmp_2 - (k-1)*tmp_1 )*interp%rk(k)

       fk0 = 0.0D0
       do j=1,npts
          fk0 = fk0 + interp%Imat(j,k+1)*interp%vtemp(j)
       end do

       fxy = fxy + pk*fk0

    end do

  end function interpolate_2d

  !===============================
  !(Nair) Bilinear interpolation for every GLL grid cell
  !===============================

  function interpol_bilinear(cart, f, interp, npts, fillvalue) result(fxy)

    integer, intent(in)               :: npts
    type (cartesian2D_t), intent(in)  :: cart
    real (kind=real_kind), intent(in) :: f(npts,npts)
    type (interpolate_t)              :: interp
    real (kind=real_kind)             :: fxy     ! value of f interpolated to (x,y)
    real (kind=real_kind), intent(in), optional :: fillvalue
    ! local variables

    real (kind=real_kind)             :: xoy(npts)
    real (kind=real_kind)             :: p,q,xp,yp ,y4(4)

    integer                           :: l,j,k, ii, jj, na,nb,nm

    xp = cart%x
    yp = cart%y

    xoy(:) = interp%glp(:)

    ! Search index along "x"  (bisection method)

    na = 1
    nb = npts
    do
       if  ((nb-na) <=  1)  exit
       nm = (nb + na)/2
       if (xp  >  xoy(nm)) then
          na = nm
       else
          nb = nm
       endif
    enddo
    ii = na

    ! Search index along "y"

    na = 1
    nb = npts
    do
       if  ((nb-na) <=  1)  exit
       nm = (nb + na)/2
       if (yp  >  xoy(nm)) then
          na = nm
       else
          nb = nm
       endif
    enddo
    jj = na

    ! GLL cell containing (xp,yp)

    y4(1) = f(ii,jj)
    y4(2) = f(ii+1,jj)
    y4(3) = f(ii+1,jj+1)
    y4(4) = f(ii,jj+1)

    if(present(fillvalue)) then
       if (any(y4==fillvalue)) then
          fxy = fillvalue
          return
       endif
    endif
       p = (xp - xoy(ii))/(xoy(ii+1) - xoy(ii))
       q = (yp - xoy(jj))/(xoy(jj+1) - xoy(jj))

       fxy = (1.0D0 - p)*(1.0D0 - q)* y4(1) + p*(1.0D0 - q) * y4(2)   &
            + p*q* y4(3) + (1.0D0 - p)*q * y4(4)

  end function interpol_bilinear

  function parametric_coordinates(sphere, corners3D,ref_map_in, corners,cartp,facenum) result (ref)

    implicit none
    type (spherical_polar_t), intent(in) :: sphere
    type (cartesian2D_t) :: ref

    type (cartesian3D_t)   :: corners3D(4)  !x,y,z coords of element corners
    integer,optional  :: ref_map_in    ! default is global variable 'cubed_sphere_map'
    ! optional arguments, only needed for ref_map=1 (equi-angle gnomonic projection):
    type (cartesian2D_t),optional   :: corners(4)    ! gnomonic coords of element corners
    type (cartesian2D_t),optional   :: cartp(np,np)  ! gnomonic coords of element data
    integer,optional  :: facenum


    ! local
    integer               :: i, MAX_NR_ITER=10
    real(kind=real_kind)  :: D(2,2),Dinv(2,2),detD,a,b,resa,resb,dela,delb,costh
    real(kind=real_kind)  :: tol_sq=1e-26
    type (spherical_polar_t) :: sphere1, sphere_tmp
    integer  :: ref_map

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! newton iteration on: ref=ref - df^-1 (ref2sphere(ref) - sphere)
    !
    ! Generic version written in terms of HOMME's 'ref2sphere' and 'Dmap' operaters,
    ! with no assumption as to the type of map (gnomonic, equi-angular, parametric)
    !
    ! Note that the coordinate increment from newton iterations is not a direction and thus 
    ! should not be converted into motion along a great circle arc - this routine
    ! correclty applies the increment by just adding it to the coordintes
    !
    ! f = ref2sphere(xvec) - sphere
    ! df = d(ref2sphere) 
    !
    ! D = diag(cos(theta),1) * d(ref2sphere)       d(ref2sphere) = diag(1/cos(theta),1)*D
    ! 
    ! df = diag(1/cos(theta),1)*D 
    ! df^-1 =  D^-1 *  diag(cos(theta),1)
    ! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (present(ref_map_in)) then
       ref_map=ref_map_in
    else
       ref_map=cubed_sphere_map
    endif
    costh=cos(sphere%lat)
    a=0
    b=0
    i=0
    do
       sphere1 = ref2sphere(a,b,corners3D,ref_map,corners,facenum)
       resa = sphere1%lon - sphere%lon
       if (resa>dd_pi) resa=resa-2*dd_pi
       if (resa<-dd_pi) resa=resa+2*dd_pi

       resb = sphere1%lat - sphere%lat 

       call Dmap(D,a,b,corners3D,ref_map,cartp,facenum)
       detD = D(1,1)*D(2,2) - D(1,2)*D(2,1)      
       Dinv(1,1) =  D(2,2)/detD
       Dinv(1,2) = -D(1,2)/detD
       Dinv(2,1) = -D(2,1)/detD
       Dinv(2,2) =  D(1,1)/detD
       
       dela =  Dinv(1,1)*costh*resa + Dinv(1,2)*resb 
       delb =  Dinv(2,1)*costh*resa + Dinv(2,2)*resb 
       a = a - dela
       b = b - delb
       i=i+1
       if ( (costh*resa)**2 + resb**2 < tol_sq .or. MAX_NR_ITER < i) exit
    end do
    ref%x=a
    ref%y=b

  end function parametric_coordinates




!
! find element containing given point, useing HOMME's standard
! equi-angular gnomonic map.
! note that with this map, only coordinate lines are great circle arcs
!
  function point_inside_equiangular(elem, sphere, sphere_xyz) result(inside)
    implicit none
    type (spherical_polar_t), intent(in)     :: sphere
    type (cartesian3D_t),     intent(in)    :: sphere_xyz
    type (element_t)        , intent(in)     :: elem
    logical                              :: inside, inside2
    integer               :: i,j
    type (cartesian2D_t) :: corners(4),sphere_xy,cart
    type (cartesian3D_t) :: corners_xyz(4),center,a,b,cross(4)
    real (kind=real_kind) :: yp(4), y, elem_diam,dotprod
    real (kind=real_kind) :: xp(4), x, xc,yc
    real (kind=real_kind) :: tol_inside
    real (kind=real_kind) :: d1,d2

    type (spherical_polar_t)    :: sphere_tmp

    inside = .false.


    ! first check if point is near the element:
    corners_xyz(:) = elem%corners3D(:)
    elem_diam = max( distance(corners_xyz(1),corners_xyz(3)), &
         distance(corners_xyz(2),corners_xyz(4)) )

    center%x = sum(corners_xyz(1:4)%x)/4
    center%y = sum(corners_xyz(1:4)%y)/4
    center%z = sum(corners_xyz(1:4)%z)/4
    if ( distance(center,sphere_xyz) > 1.0*elem_diam ) return

    tol_inside = 1e-10*elem_diam**2
    ! the point is close to the element, so project both to cubed sphere
    ! and perform contour integral
    sphere_xy=sphere2cubedsphere(sphere,elem%FaceNum)
    x = sphere_xy%x
    y = sphere_xy%y
    do i=1,4
      xp(i) = elem%corners(i)%x
      yp(i) = elem%corners(i)%y
    end do


    if (debug) then
       print *,'point: ',x,y,elem%FaceNum
       print *,'element:'
       write(*,'(a,4e16.8,a)') 'x=[',xp(1:4),']'
       write(*,'(a,4e16.8,a)') 'y=[',yp(1:4),']'

       ! first check if centroid is in this element (sanity check)
       sphere_tmp=change_coordinates(center)
       sphere_xy=sphere2cubedsphere(sphere_tmp,elem%FaceNum)
       xc=sphere_xy%x
       yc=sphere_xy%y
       print *,'cross product with centroid: all numbers should be negative'
       j = 4
       do i=1,4
          print *,i,(xc-xp(j))*(yp(i)-yp(j))  - (yc-yp(j))*(xp(i)-xp(j))
          j = i  ! within this loopk j = i-1
       end do

       print *,'cross product with search point'
       j = 4
       do i=1,4
          print *,i,(x-xp(j))*(yp(i)-yp(j))  - (y-yp(j))*(xp(i)-xp(j))
          j = i  ! within this loopk j = i-1
       end do
    endif

    
    j = 4
    do i=1,4
      ! a = x-xp(j), y-yp(j)
      ! b = xp(i)-xp(j), yp(i)-yp(j)
      ! compute a cross b:
      if ( -( (x-xp(j))*(yp(i)-yp(j))  - (y-yp(j))*(xp(i)-xp(j))) > tol_inside ) then
         return
      endif
      j = i  ! within this loopk j = i-1
    end do
    ! all cross products were negative, must be inside:
    inside=.true.
  end function point_inside_equiangular


!
! find if quad contains given point, with quad edges assumed to be great circle arcs
! this will work with any map where straight lines are mapped to great circle arcs.
! (thus it will fail on unstructured grids using the equi-angular gnomonic map)
!
  function point_inside_quad(corners_xyz, sphere_xyz) result(inside)
    implicit none
    type (cartesian3D_t),     intent(in)    :: sphere_xyz
    type (cartesian3D_t)    , intent(in)    :: corners_xyz(4)
    logical                              :: inside, inside2
    integer               :: i,j,ii
    type (cartesian2D_t) :: corners(4),sphere_xy,cart
    type (cartesian3D_t) :: center,a,b,cross(4)
    real (kind=real_kind) :: yp(4), y, elem_diam,dotprod
    real (kind=real_kind) :: xp(4), x
    real (kind=real_kind) :: d1,d2, tol_inside = 1e-12

    type (spherical_polar_t)   :: sphere  ! debug

    inside = .false.

    ! first check if point is near the corners:
    elem_diam = max( distance(corners_xyz(1),corners_xyz(3)), &
         distance(corners_xyz(2),corners_xyz(4)) )

    center%x = sum(corners_xyz(1:4)%x)/4
    center%y = sum(corners_xyz(1:4)%y)/4
    center%z = sum(corners_xyz(1:4)%z)/4
    if ( distance(center,sphere_xyz) > 1.0*elem_diam ) return

    j = 4
    do i=1,4
      ! outward normal to plane containing j->i edge:  corner(i) x corner(j)
      ! sphere dot (corner(i) x corner(j) ) = negative if inside
       cross(i)%x =  corners_xyz(i)%y*corners_xyz(j)%z - corners_xyz(i)%z*corners_xyz(j)%y
       cross(i)%y =-(corners_xyz(i)%x*corners_xyz(j)%z - corners_xyz(i)%z*corners_xyz(j)%x)
       cross(i)%z =  corners_xyz(i)%x*corners_xyz(j)%y - corners_xyz(i)%y*corners_xyz(j)%x
       dotprod = cross(i)%x*sphere_xyz%x + cross(i)%y*sphere_xyz%y +&
               cross(i)%z*sphere_xyz%z
       j = i  ! within this loopk j = i-1

       ! dot product is proportional to elem_diam. positive means outside,
       ! but allow machine precision tolorence: 
       if (dotprod > tol_inside*elem_diam) return 
       !if (dotprod > 0) return 
    end do
    inside=.true.
    return
  end function point_inside_quad

!
! find element containing given point, with element edges assumed to be great circle arcs
! this will work with any map where straight lines are mapped to great circle arcs.
! (thus it will fail on unstructured grids using the equi-angular gnomonic map)
!
  function point_inside_gc(elem, sphere_xyz) result(inside)
    implicit none
    type (cartesian3D_t),     intent(in)    :: sphere_xyz
    type (element_t)        , intent(in)     :: elem
    logical                              :: inside, inside2
    integer               :: i,j,ii
    type (cartesian2D_t) :: corners(4),sphere_xy,cart
    type (cartesian3D_t) :: corners_xyz(4),center,a,b,cross(4)
    real (kind=real_kind) :: yp(4), y, elem_diam,dotprod
    real (kind=real_kind) :: xp(4), x
    real (kind=real_kind) :: d1,d2, tol_inside = 1e-12

    type (spherical_polar_t)   :: sphere  ! debug

    inside = .false.

    ! first check if point is near the element:
    corners_xyz(:) = elem%corners3D(:)
    elem_diam = max( distance(corners_xyz(1),corners_xyz(3)), &
         distance(corners_xyz(2),corners_xyz(4)) )

    center%x = sum(corners_xyz(1:4)%x)/4
    center%y = sum(corners_xyz(1:4)%y)/4
    center%z = sum(corners_xyz(1:4)%z)/4
    if ( distance(center,sphere_xyz) > 1.0*elem_diam ) return


#if 0
    ! sanity check: see if center is inside element
    !       2 1
    !       3 4      (1) x (4) = outward normal
    !                 dot product with centroid should be negative (inside)
    j = 4
    do i=1,4
      ! outward normal to plane containing j->i edge:  corner(i) x corner(j)
      ! sphere dot (corner(i) x corner(j) ) = negative if on left (inside for counterclockwise)
      !                                       positive if on right 
      !
       cross(i)%x =  corners_xyz(i)%y*corners_xyz(j)%z - corners_xyz(i)%z*corners_xyz(j)%y
       cross(i)%y =-(corners_xyz(i)%x*corners_xyz(j)%z - corners_xyz(i)%z*corners_xyz(j)%x)
       cross(i)%z =  corners_xyz(i)%x*corners_xyz(j)%y - corners_xyz(i)%y*corners_xyz(j)%x
       dotprod = cross(i)%x*center%x + cross(i)%y*center%y +&
               cross(i)%z*center%z
       j = i  ! within this loopk j = i-1
       if (dotprod>0) then
          print *,'i,dotprod=',i,dotprod
          print *,'error: center is outside element. face=',elem%facenum
          print *,'corner (cubedsphere)'
          do ii=1,4
             print *,elem%corners(ii)%x,elem%corners(ii)%y
          enddo
          sphere=change_coordinates(center)
          print *,'center lon,lat: ',sphere%lon,sphere%lat
          do ii=1,4
             sphere=projectpoint(elem%corners(ii),elem%facenum)
             print *,sphere%lon,sphere%lat
          enddo
          print *,'test point'
          cart%x = -3.1459/4
          cart%y = 0
          sphere=projectpoint(cart,5)
        !  print *,'0,0 ->lon,lat: ',sphere%lon,sphere%lat
          cart%x = -3.1459/4 + .1
          cart%y = 0
          sphere=projectpoint(cart,5)
        !  print *,'+1,0 ->lon,lat: ',sphere%lon,sphere%lat
          cart%x = -3.1459/4 + .1
          cart%y  =0 + .1
          sphere=projectpoint(cart,5)
        !  print *,'+1,+1 ->lon,lat: ',sphere%lon,sphere%lat
          cart%x = -3.1459/4 
          cart%y  = 0 + .1
          sphere=projectpoint(cart,5)
        !  print *,'+1,+1 ->lon,lat: ',sphere%lon,sphere%lat

       endif
    end do
#endif

    j = 4
    do i=1,4
      ! outward normal to plane containing j->i edge:  corner(i) x corner(j)
      ! sphere dot (corner(i) x corner(j) ) = negative if inside
       cross(i)%x =  corners_xyz(i)%y*corners_xyz(j)%z - corners_xyz(i)%z*corners_xyz(j)%y
       cross(i)%y =-(corners_xyz(i)%x*corners_xyz(j)%z - corners_xyz(i)%z*corners_xyz(j)%x)
       cross(i)%z =  corners_xyz(i)%x*corners_xyz(j)%y - corners_xyz(i)%y*corners_xyz(j)%x
       dotprod = cross(i)%x*sphere_xyz%x + cross(i)%y*sphere_xyz%y +&
               cross(i)%z*sphere_xyz%z
       j = i  ! within this loopk j = i-1

       !if (dotprod>0 .and. dotprod/elem_diam < 1e-5) print *,dotprod/elem_diam

       ! dot product is proportional to elem_diam. positive means outside,
       ! but allow machine precision tolorence: 
       if (dotprod > tol_inside*elem_diam) return 
       !if (dotprod > 0) return 
    end do
    inside=.true.
    return
  end function point_inside_gc



  !================================================
  !  (Nair) Cube face index and local coordinates
  !================================================

  subroutine cube_facepoint_ne(sphere,ne,cart, number)
    use coordinate_systems_mod, only : cube_face_number_from_sphere, sphere2cubedsphere
    implicit none

    type (spherical_polar_t), intent (in) :: sphere
    integer             , intent(in)      :: ne
    type (cartesian2D_t), intent(out)     :: cart
    integer             , intent(out)     :: number

    real (kind=real_kind) :: xp,yp
    type (cartesian2D_t)  :: cube
    integer               :: ie, je, face_no
    real (kind=real_kind) :: x1,x2
    real (kind=real_kind) :: dx

    face_no = cube_face_number_from_sphere(sphere)
    cube    = sphere2cubedsphere(sphere, face_no)
    xp      = cube%x
    yp      = cube%y

    ! MNL: for uniform grids (on cube face), analytic solution is fine
    x1 = xp + 0.25D0*DD_PI
    x2 = yp + 0.25D0*DD_PI

    dx = (0.5D0*DD_PI)/ne
    ie = INT(ABS(x1)/dx)
    je = INT(ABS(x2)/dx)
    ! if we are exactly on an element edge, we can put the point in
    ! either the ie or ie+1 element, EXCEPT if ie==ne.
    if ( ABS(x1) < ne*dx  ) ie=ie+1
    if ( ABS(x2) < ne*dx  ) je=je+1
    if (ie>ne .or. je>ne) then
       write(iulog,*)'ERROR: ',ie,je,ne
       write(iulog,*)'lat,lon=',sphere%lat,sphere%lon
       write(iulog,*)'face no=',face_no
       write(iulog,*)x1,x2,x1/dx,x2/dx
       call abortmp('interpolate_mod: bad argument')
    endif

    ! bug fix MT 1/2009.  This was creating a plotting error at
    ! the row of elements in iface=2 at 50 degrees (NE=16 128x256 lat/lon grid)
    ! For point on element edge, we can have ie=2, but x1=dx
    ! but if ie>1, we must execute this statement.
    ! The only time we can skip this statement is if ie=1, but then
    ! the statement has no effect, so lets never skip it:
    !    if (x1 > dx ) then
    x1 = x1 - dble(ie-1)*dx
    !    endif

    x1 = 2.0D0*(x1/dx)-1.0D0

    !    if (x2 > dx ) then    ! removed MT 1/2009, see above
    x2 = x2 - dble(je-1)*dx
    !    endif

    x2 = 2.0D0*(x2/dx)-1.0D0

    ! coordinates within an element [-1,1]
    cart%x = x1
    cart%y = x2
    number = ie + (je-1)*ne + (face_no-1)*ne*ne
  end subroutine cube_facepoint_ne
  !================================================
  !  (Nair) Cube face index and local coordinates
  !================================================


  subroutine cube_facepoint_unstructured(sphere,cart, number, elem)
    use coordinate_systems_mod, only : cube_face_number_from_sphere, &
                                       sphere2cubedsphere,change_coordinates,cube_face_number_from_cart
    implicit none

    type (element_t)     , intent(in), target :: elem(:)
    type (spherical_polar_t), intent (in) :: sphere
    type (cartesian2D_t), intent(out)     :: cart
    integer             , intent(out)     :: number

    integer               :: ii
    Logical               :: found
    type (cartesian3D_t)       :: sphere_xyz
    type (cartesian2D_t)  :: cube
    sphere_xyz=spherical_to_cart(sphere)

    number=-1
!    print *,'WARNING: using GC map'
    do ii = 1,nelemd
       ! for equiangular gnomonic map:
       ! unstructed grid element edges are NOT great circles
       if (cubed_sphere_map==0) then
          found = point_inside_equiangular(elem(ii), sphere, sphere_xyz)
       else 
          ! assume element edges are great circle arcs:
          found = point_inside_gc(elem(ii), sphere_xyz)
       endif

       if (found) then
          number = ii
          cart = parametric_coordinates(sphere, elem(ii)%corners3D,&
               cubed_sphere_map,elem(ii)%corners,elem(ii)%cartp,elem(ii)%facenum)
          exit
       end if
    end do
  end subroutine cube_facepoint_unstructured


  subroutine interp_init()
    type (quadrature_t)   :: gp

    gp = gausslobatto(np)
    call interpolate_create(gp,interp_p)
  end subroutine interp_init


  subroutine setup_latlon_interp(elem,interpdata,par)
    !
    ! initialize interpolation data structures to interpolate to a lat-lon grid
    !
    !

    implicit none
    type (element_t)     , intent(in), target :: elem(:)
    type (parallel_t)      , intent(in)       :: par
    type (interpdata_t)  , intent(out)        :: interpdata(:)

    ! local
    integer i,j,ii,count_total,n_interp,count_max
    integer ngrid, number, elem_num, plat
    integer countx, missing_pts,ierr
    integer :: npts_mult_claims,max_claims

    real (kind=real_kind)    ::  dp,latdeg(nlat+1),clat(nlat+1),w(nlat+1),w_staggered(nlat)
    real (kind=real_kind)    ::  clat_staggered(nlat),latdeg_st(nlat),err,err2

    type (spherical_polar_t) :: sphere
    type (cartesian2D_t)     :: cart
    type (cartesian3D_t)     :: sphere_xyz,sphere2_xyz

    type (quadrature_t)       :: gp


    ! Array to make sure each interp point is on exactly one process
    type (cartesian2D_t),allocatable    :: cart_vec(:,:)
    integer :: k
    integer, allocatable :: global_elem_gid(:,:),local_elem_gid(:,:), local_elem_num(:,:)

    ! these arrays often are too large for stack, so lets make sure
    ! they go on the heap:
    allocate(local_elem_num(nlat,nlon))
    allocate(local_elem_gid(nlat,nlon))
    allocate(global_elem_gid(nlat,nlon))
    allocate(cart_vec(nlat,nlon))

    if (par%masterproc) then
       write(iulog,'(a,i4,a,i4,a)') 'Initializing ',nlat,' x ',nlon,' lat-lon interpolation grid: '
    endif

    do ii=1,nelemd
       interpdata(ii)%n_interp=0  ! reset counter
    enddo

    if (associated(lat))then
       deallocate(lat)
       nullify(lat)
    endif
    if (associated(gweight))then
       deallocate(gweight)
       nullify(gweight)
    endif

    if (associated(lon))then
       deallocate(lon)
       nullify(lon)
    endif

    allocate(lat(nlat))
    allocate(gweight(nlat))
    allocate(lon(nlon))
    call interp_init()
    gweight=0
    do i=1,nlon
       lon(i)=2*dd_pi*(i-1)/nlon
    enddo
    if (gridtype==1) then
       do j=1,nlat
          lat(j) = -dd_pi/2 + dd_pi*(j-1)/(nlat-1)
       end do
       plat=nlat
    endif
    if (gridtype==2) then
       gp=gauss(nlat)
       do j=1,nlat
          lat(j) = asin(gp%points(j))
          gweight(j) = gp%weights(j)
       end do
    endif
    if (gridtype==3) then
       do j=1,nlat
          lat(j) = -dd_pi/2 + dd_pi*(j-.5d0)/nlat
       end do
       plat=nlat+1
    endif

    if (gridtype==1 .or. gridtype==3) then
       ! gridtype=1    plat=nlat    gweight(1:nlat)=w(1:plat)
       ! gridtype=3    plat=nlat+1  gweight(1:nlat)=w_staggered(1:plat-1)

       ! L-R dynamics uses a regular latitude distribution (not gausian).
       ! The algorithm below is a bastardized version of LSM: map.F.
       dp = 180d0/(plat-1)
       do j = 1, plat
          latdeg(j) = -90d0 + (j-1)*dp
          clat(j) = latdeg(j)*dd_pi/180d0
       end do

       ! Calculate latitudes for the staggered grid

       do j = 1, plat-1
          clat_staggered(j) = (clat(j) + clat(j+1)) / 2
          latdeg_st     (j) = clat_staggered(j)*180d0/dd_pi
       end do

       ! Weights are defined as cos(phi)*(delta-phi)
       ! For a sanity check, the sum of w across all lats should be 2, or 1 across
       ! half of the latitudes.

       do j = 2, plat-1
          w(j) = sin(clat_staggered(j)) - sin(clat_staggered(j-1))
       end do
       w(1) = sin(clat_staggered(1)) + 1
       w(plat) = w(1)

       ! with nlat=2048, this error was 4e-16
       if (abs(sum(w(1:plat)) - 2) > 1e-8) then
          write(iulog,*) 'interpolate_mod: w weights do not sum to 2. sum=',sum(w(1:plat))
          call abortmp('interpolate_mod: weights do not sum to 2.')
       end if

       dp = dd_pi / (plat-1)
       do j = 1, plat-1
          w_staggered(j) = sin(clat(j+1)) - sin(clat(j))
       end do


       if (abs(sum(w_staggered(1:plat-1)) - 2) > 1e-8) then
          write(iulog,*) 'interpolate_mod: staggered weights do not sum to 2. sum=',sum(w_staggered(1:plat-1))
          call abortmp('interpolate_mod: weights do not sum to 2.')
       end if

       if (gridtype==1) then
          gweight(1:nlat)=w(1:plat)
       endif
       if (gridtype==3) then
          gweight(1:nlat)=w_staggered(1:plat-1)
       endif
    endif


    ! go through once, counting the number of points on each element
    sphere%r=1
    local_elem_num  = -1
    local_elem_gid  = -1
    global_elem_gid = -1
    err=0
    do j=1,nlat
       do i=1,nlon
          sphere%lat=lat(j)
          sphere%lon=lon(i)

          !debug=.false.
          !if (j==18 .and. i==95) debug=.true.
          number = -1
          if ( (cubed_sphere_map /= 0) .or. MeshUseMeshFile) then
             call cube_facepoint_unstructured(sphere, cart, number, elem)
             if (number /= -1) then
                ! If points are outside element but within tolerance, move to boundary
                if (cart%x + 1.0d0.le.0.0d0) cart%x = -1.0d0
                if (cart%x - 1.0d0.ge.0.0d0) cart%x = 1.0d0
                if (cart%y + 1.0d0.le.0.0d0) cart%y = -1.0d0
                if (cart%y - 1.0d0.ge.0.0d0) cart%y = 1.0d0

                local_elem_num(j,i) = number
                local_elem_gid(j,i) = elem(number)%vertex%number
                cart_vec(j,i)    = cart  ! local element coordiante of interpolation point
             endif
          else
             call cube_facepoint_ne(sphere, ne, cart, number)
             ! the sphere point belongs to the element number on face = face_no.
             ! do I own this element?
             if (number /= -1) then
                do ii=1,nelemd
                   if (number == elem(ii)%vertex%number) then
                      local_elem_gid(j,i) = number
                      local_elem_num(j,i) = ii
                      cart_vec(j,i)        = cart   ! local element coordinate found above
                      exit
                   endif
                enddo
             endif
          endif
          ii=local_elem_num(j,i)
          if (ii /= -1) then
             ! compute error: map 'cart' back to sphere and compare with original
             ! interpolation point:
             sphere2_xyz = spherical_to_cart( ref2sphere(cart%x,cart%y,     &
                  elem(ii)%corners3D,cubed_sphere_map,elem(ii)%corners,elem(ii)%facenum ))
             sphere_xyz = spherical_to_cart(sphere)
             err=max(err,distance(sphere2_xyz,sphere_xyz))
          endif
       enddo
       if (par%masterproc) then
          if ((MOD(j,64).eq.1).or.(j.eq.nlat)) then
             print *,'finished latitude ',j,' of ',nlat
          endif
       endif
    enddo
    err2=err
#ifdef _MPI
    call MPI_Allreduce(err,err2,1,MPIreal_t,MPI_MAX,par%comm,ierr)
#endif
    if (par%masterproc) then
       write(iulog,'(a,e12.4)') 'Max interpolation point search error: ',err2
    endif

    ! if multile elements claim a interpolation point, take the one with largest gid:
    global_elem_gid = local_elem_gid
#ifdef _MPI
    call MPI_Allreduce(local_elem_gid, global_elem_gid, nlat*nlon, MPIinteger_t, MPI_MAX, par%comm,ierr)
#endif

    missing_pts=0
    do j=1,nlat
       do i=1,nlon
          if (global_elem_gid(j,i) == -1 ) then
             missing_pts = missing_pts + 1
             if (par%masterproc) &
                  print *,'Error: point not claimed by any element j,i,lat(j),lon(i)=',j,i,lat(j),lon(i)
          else if (local_elem_gid(j,i) == global_elem_gid(j,i)  ) then
             ii = local_elem_num(j,i)
             interpdata(ii)%n_interp = interpdata(ii)%n_interp + 1
          endif
       end do
    end do

    countx=maxval(interpdata(1:nelemd)%n_interp)
    count_max = countx
#ifdef _MPI
    call MPI_Allreduce(countx,count_max,1,MPIinteger_t,MPI_MAX,par%comm,ierr)
#endif

    if (par%masterproc) then
       write(iulog,'(a,i6)') 'Maximum number of interpolation points claimed by an element: ',count_max
    endif

    ! allocate storage
    do ii=1,nelemd
       ngrid = interpdata(ii)%n_interp
       if(interpdata(ii)%first_entry)then
          NULLIFY(interpdata(ii)%interp_xy)
          NULLIFY(interpdata(ii)%ilat)
          NULLIFY(interpdata(ii)%ilon)

          interpdata(ii)%first_entry=.FALSE.
       endif
       if(associated(interpdata(ii)%interp_xy))then
          if(size(interpdata(ii)%interp_xy)>0)deallocate(interpdata(ii)%interp_xy)
       endif
       if(associated(interpdata(ii)%ilat))then
          if(size(interpdata(ii)%ilat)>0)deallocate(interpdata(ii)%ilat)
       endif

       if (associated(interpdata(ii)%ilon))then
          if(size(interpdata(ii)%ilon)>0)deallocate(interpdata(ii)%ilon)
       endif
       allocate(interpdata(ii)%interp_xy( ngrid ) )
       allocate(interpdata(ii)%ilat( ngrid ) )
       allocate(interpdata(ii)%ilon( ngrid ) )
       interpdata(ii)%n_interp=0  ! reset counter
    enddo
    do j=1,nlat
       do i=1,nlon
          if (local_elem_gid(j,i) == global_elem_gid(j,i) .and. &
               local_elem_gid(j,i) /= -1 ) then
             ii = local_elem_num(j,i)
             ngrid = interpdata(ii)%n_interp + 1
             interpdata(ii)%n_interp = ngrid
             interpdata(ii)%interp_xy( ngrid )   = cart_vec(j,i)
             interpdata(ii)%ilon( ngrid ) = i
             interpdata(ii)%ilat( ngrid ) = j
          endif
       enddo
    enddo

    ! now lets compute the number of points that were claimed by
    ! more than one element:
    do j=1,nlat
       do i=1,nlon
          if (local_elem_gid(j,i) == -1) then
             local_elem_gid(j,i)=0
          else
             local_elem_gid(j,i)=1
          endif
       enddo
    enddo
    global_elem_gid = local_elem_gid
#ifdef _MPI
    call MPI_Allreduce(local_elem_gid, global_elem_gid, nlat*nlon, MPIinteger_t, MPI_SUM, par%comm,ierr)
#endif
    if (par%masterproc) then
       countx=0
       do j=1,nlat
          do i=1,nlon
             if (global_elem_gid(j,i)>1) countx=countx+1
          enddo
       enddo
       npts_mult_claims=countx
       max_claims=maxval(global_elem_gid)
    endif

    if (par%masterproc) then
       print *,'Number of interpolation points claimed by more than one element: ',npts_mult_claims
       print *,'max number of elements which claimed the same interpolation point:',max_claims
    endif
    
    deallocate(global_elem_gid)
    deallocate(local_elem_num)
    deallocate(local_elem_gid)
    deallocate(cart_vec)

    ! check if every point in interpolation grid was claimed by an element:
    if (missing_pts>0) then
       count_total = nlat*nlon
       if(par%masterproc) then
          write(iulog,"(3A,I4,A,I7,a,i5)")"Error:",__FILE__," ",__LINE__," count_total:",count_total," missing:",missing_pts
       end if
       call syncmp(par)
       call abortmp('Error: interpolation points not claimed by any element')
    endif


  end subroutine setup_latlon_interp



! interpolate_scalar
!
! Interpolate a scalar field given in an element (fld_cube) to the points in
! interpdata%interp_xy(i), i=1 .. interpdata%n_interp.
!
! Note that it is possible the given element contains none of the interpolation points
! =======================================
subroutine interpolate_ce(cart,fld_cube,npts,fld, fillvalue)
  type (cartesian2D_t) :: cart
  integer                  ::  npts
  real (kind=real_kind)    ::  fld_cube(npts,npts) ! cube field
  real (kind=real_kind)    ::  fld          ! field at new grid lat,lon coordinates
  real (kind=real_kind), intent(in), optional :: fillvalue
  ! Local variables
  type (interpolate_t), pointer  ::  interp          ! interpolation structure

  integer :: ne

  integer :: i


  if (npts==np) then
     interp => interp_p
  else
     call abortmp('Error in interpolate_scalar(): must be called with p or v grid data')
  endif

  fld=interpolate_2d(cart,fld_cube,interp,npts,fillvalue)

end subroutine interpolate_ce



  ! =======================================
  ! interpolate_scalar
  !
  ! Interpolate a scalar field given in an element (fld_cube) to the points in
  ! interpdata%interp_xy(i), i=1 .. interpdata%n_interp.
  !
  ! Note that it is possible the given element contains none of the interpolation points
  ! =======================================
  subroutine interpolate_scalar2d(interpdata,fld_cube,npts,fld, fillvalue)
    integer                  ::  npts
    real (kind=real_kind)    ::  fld_cube(npts,npts) ! cube field
    real (kind=real_kind)    ::  fld(:)          ! field at new grid lat,lon coordinates
    type (interpdata_t)         ::  interpdata
    real (kind=real_kind), intent(in), optional :: fillvalue
    ! Local variables
    type (interpolate_t), pointer  ::  interp          ! interpolation structure

    integer :: ne

    integer :: i

    type (cartesian2D_t) :: cart

    if (npts==np) then
       interp => interp_p
    else
       call abortmp('Error in interpolate_scalar(): must be called with p or v grid data')
    endif

       ! Choice for Native (high-order) or Bilinear interpolations
    if(present(fillvalue)) then
       if (itype == 0) then
          do i=1,interpdata%n_interp
             fld(i)=interpolate_2d(interpdata%interp_xy(i),fld_cube,interp,npts,fillvalue)
          end do
       elseif (itype == 1) then
          do i=1,interpdata%n_interp
             fld(i)=interpol_bilinear(interpdata%interp_xy(i),fld_cube,interp,npts,fillvalue)
          end do
       end if
    else
       if (itype == 0) then
          do i=1,interpdata%n_interp
             fld(i)=interpolate_2d(interpdata%interp_xy(i),fld_cube,interp,npts)
          end do
       elseif (itype == 1) then
          do i=1,interpdata%n_interp
             fld(i)=interpol_bilinear(interpdata%interp_xy(i),fld_cube,interp,npts)
          end do
       end if
    endif


  end subroutine interpolate_scalar2d
  subroutine interpolate_scalar3d(interpdata,fld_cube,npts,nlev,fld, fillvalue)
    integer , intent(in)                 ::  npts, nlev
    real (kind=real_kind)    ::  fld_cube(npts,npts,nlev) ! cube field
    real (kind=real_kind)    ::  fld(:,:)          ! field at new grid lat,lon coordinates
    type (interpdata_t)         ::  interpdata
    real (kind=real_kind), intent(in), optional :: fillvalue
    ! Local variables
    type (interpolate_t), pointer  ::  interp          ! interpolation structure

    integer :: ne

    integer :: i, k

    type (cartesian2D_t) :: cart

    if (npts==np) then
       interp => interp_p
    else
       call abortmp('Error in interpolate_scalar(): must be called with p or v grid data')
    endif

    ! Choice for Native (high-order) or Bilinear interpolations
    if(present(fillvalue)) then
       if (itype == 0) then
          do k=1,nlev
             do i=1,interpdata%n_interp
                fld(i,k)=interpolate_2d(interpdata%interp_xy(i),fld_cube(:,:,k),interp,npts,fillvalue)
             end do
          end do
       elseif (itype == 1) then
          do k=1,nlev
             do i=1,interpdata%n_interp
                fld(i,k)=interpol_bilinear(interpdata%interp_xy(i),fld_cube(:,:,k),interp,npts,fillvalue)
             end do
          end do
       endif
    else
       if (itype == 0) then
          do k=1,nlev
             do i=1,interpdata%n_interp
                fld(i,k)=interpolate_2d(interpdata%interp_xy(i),fld_cube(:,:,k),interp,npts)
             end do
          end do
       elseif (itype == 1) then
          do k=1,nlev
             do i=1,interpdata%n_interp
                fld(i,k)=interpol_bilinear(interpdata%interp_xy(i),fld_cube(:,:,k),interp,npts)
             end do
          end do
       else
          write(iulog,*) itype
          call abortmp("wrong interpolation type")
       endif
    endif
  end subroutine interpolate_scalar3d



  ! =======================================
  ! interpolate_vector
  !
  ! Interpolate a vector field given in an element (fld_cube)
  ! to the points in interpdata%interp_xy(i), i=1 .. interpdata%n_interp.
  !
  ! input_coords = 0    fld_cube given in lat-lon
  ! input_coords = 1    fld_cube given in contravariant
  !
  ! Note that it is possible the given element contains none of the
  ! interpolation points
  ! =======================================
  subroutine interpolate_vector3d(interpdata,elem,fld_cube,nlev,fld,input_coords,fillvalue)
    implicit none
    type (interpdata_t),intent(in)       ::  interpdata
    type (element_t), intent(in)         ::  elem
    integer, intent(in)                  ::  nlev
    real (kind=real_kind), intent(in)    ::  fld_cube(np,np,2,nlev) ! vector field
    real (kind=real_kind), intent(out)   ::  fld(:,:,:)          ! field at new grid lat,lon coordinates
    real (kind=real_kind), intent(in),optional :: fillvalue
    integer, intent(in)                  ::  input_coords

    ! Local variables
    real (kind=real_kind)    ::  fld_contra(np,np,2,nlev) ! vector field
    real (kind=real_kind)    ::  fld_lonlat(np,np,2,nlev) ! vector field in lonlat
    real (kind=real_kind)    ::  fld_cart(np,np,3,nlev) ! vector field in cartesian
    real (kind=real_kind), allocatable    ::  fld_cart_interp(:,:,:) ! vector field in cartesian
    type (interpolate_t), pointer  ::  interp          ! interpolation structure

    real (kind=real_kind)    ::  v1,v2
    real (kind=real_kind)    ::  D(2,2)   ! derivative of gnomonic mapping
    real (kind=real_kind)    ::  JJ(2,2), tmpD(2,2)   ! derivative of gnomonic mapping
    real (kind=real_kind)    ::  Km(3,2) !transform from/to lonlat basis for vectors at interp. point
    real (kind=real_kind)    ::  llon, llat
    
    integer :: i,j,k 
    integer :: ilon, ilat
    type (cartesian2D_t) :: cart
    
    
    if(present(fillvalue)) then
       if (any(fld_cube==fillvalue)) then
          fld = fillvalue
          return
       end if
    end if

    interp => interp_p

    if (input_coords ==1 ) then
       ! convert to lonlat
       do k=1,nlev
          do j=1,np
             do i=1,np
                ! contra -> lonlat
                fld_lonlat(i,j,1,k) = elem%D(i,j,1,1)*fld_cube(i,j,1,k) + elem%D(i,j,1,2)*fld_cube(i,j,2,k)
                fld_lonlat(i,j,2,k) = elem%D(i,j,2,1)*fld_cube(i,j,1,k) + elem%D(i,j,2,2)*fld_cube(i,j,2,k)
                
                ! lonlat -> cart
                fld_cart(i,j,1,k) = elem%vec_sphere2cart(i,j,1,1)*fld_lonlat(i,j,1,k) + &
                                    elem%vec_sphere2cart(i,j,1,2)*fld_lonlat(i,j,2,k)
                fld_cart(i,j,2,k) = elem%vec_sphere2cart(i,j,2,1)*fld_lonlat(i,j,1,k) + &
                                    elem%vec_sphere2cart(i,j,2,2)*fld_lonlat(i,j,2,k)
                ! vec_sphere2cart(...,3,1) = 0
                ! fld_cart(i,j,3,k) =
                ! elem%vec_sphere2cart(i,j,3,1)*fld_lonlat(i,j,1,k) +
                ! elem%vec_sphere2cart(i,j,3,2)*fld_lonlat(i,j,2,k)
                fld_cart(i,j,3,k) = elem%vec_sphere2cart(i,j,3,2)*fld_lonlat(i,j,2,k)
             enddo
          enddo
       end do
    elseif(input_coords == 0) then 
       fld_lonlat = fld_cube
       do k=1,nlev
          do j=1,np
             do i=1,np               
                ! lonlat -> cart
                fld_cart(i,j,1,k) = elem%vec_sphere2cart(i,j,1,1)*fld_lonlat(i,j,1,k) + &
                                    elem%vec_sphere2cart(i,j,1,2)*fld_lonlat(i,j,2,k)
                fld_cart(i,j,2,k) = elem%vec_sphere2cart(i,j,2,1)*fld_lonlat(i,j,1,k) + &
                                    elem%vec_sphere2cart(i,j,2,2)*fld_lonlat(i,j,2,k)
                ! vec_sphere2cart(...,3,1) = 0
                ! fld_cart(i,j,3,k) =
                ! elem%vec_sphere2cart(i,j,3,1)*fld_lonlat(i,j,1,k) +
                ! elem%vec_sphere2cart(i,j,3,2)*fld_lonlat(i,j,2,k)
                fld_cart(i,j,3,k) = elem%vec_sphere2cart(i,j,3,2)*fld_lonlat(i,j,2,k)
             enddo
          enddo
       end do       
    else
       call abortmp('Error in interpolate_vector3d(): Unknown value of input_coords. Use 0 or 1.')
    endif

    allocate(fld_cart_interp(interpdata%n_interp,3,nlev))

    ! Choice for Native (high-order) or Bilinear interpolations
    if (itype == 0) then
       do k=1,nlev
          do i=1,interpdata%n_interp
             do j=1,3
                fld_cart_interp(i,j,k)=interpolate_2d(interpdata%interp_xy(i),fld_cart(:,:,j,k),interp,np)
             enddo
          end do
       end do
    elseif (itype == 1) then
       do k=1,nlev
          do i=1,interpdata%n_interp
             do j=1,3
                fld_cart_interp(i,j,k)=interpol_bilinear(interpdata%interp_xy(i),fld_cart(:,:,j,k),interp,np)
             end do
           enddo  
       end do
    else
       call abortmp("Error in interpolate_vector3d(): wrong interpolation type")
    endif


    do i=1,interpdata%n_interp
       ! compute D(:,:) at the point elem%interp_cube(i)
       call dmap(D,interpdata%interp_xy(i)%x,interpdata%interp_xy(i)%y,&
            elem%corners3D,cubed_sphere_map,elem%cartp,elem%facenum)
       do k=1,nlev
          ! convert fld from contra->latlon
          v1 = fld(i,k,1)
          v2 = fld(i,k,2)

          fld(i,k,1)=D(1,1)*v1 + D(1,2)*v2
          fld(i,k,2)=D(2,1)*v1 + D(2,2)*v2
       end do
    end do
    
    do i=1,interpdata%n_interp
       ! convert from cart to lonlat: we need to recover matrix K, for that we
       ! need lon,lat at the interp. point.
       ilon = interpdata%ilon(i)
       ilat = interpdata%ilat(i)
       llon = lon(ilon)
       llat = lat(ilat)
       Km(1,1) = -SIN(llon)
       Km(2,1) =  COS(llon)
       Km(3,1) =  0.0_real_kind
       Km(1,2) = -SIN(llat)*COS(llon)
       Km(2,2) = -SIN(llat)*SIN(llon)
       Km(3,2) =  COS(llat)
       do k=1,nlev
          fld(i,k,1) = Km(1,1)*fld_cart_interp(i,1,k) + Km(2,1)*fld_cart_interp(i,2,k) ! Km(3,1) = 0
          fld(i,k,2) = Km(1,2)*fld_cart_interp(i,1,k) + Km(2,2)*fld_cart_interp(i,2,k) + Km(3,2)*fld_cart_interp(i,3,k) 
       enddo
    end do

    deallocate(fld_cart_interp)

  end subroutine interpolate_vector3d



#ifndef CAM
  function var_is_vector_uvar(name)
    character(len=*), intent(in) :: name
    integer :: i, var_is_vector_uvar, null_index

    var_is_vector_uvar=0
    null_index=0
    do i=1,MAX_VECVARS
       if(trim(vector_uvars(i)).eq. '') then
          null_index=i
          exit
       endif
       if(trim(vector_uvars(i)).eq.name) then
          var_is_vector_uvar=i
          exit
       end if
    end do
  end function var_is_vector_uvar


  function var_is_vector_vvar(name)
    character(len=*), intent(in) :: name
    integer :: i, var_is_vector_vvar, null_index

    var_is_vector_vvar=0
    null_index=0
    do i=1,MAX_VECVARS
       if(trim(vector_vvars(i)).eq. '') then
          null_index=i
          exit
       endif
       if(trim(vector_vvars(i)).eq.name) then
          var_is_vector_vvar=i
          exit
       end if
    end do

  end function var_is_vector_vvar
#endif
! ^ ifndef CAM


end module interpolate_mod
