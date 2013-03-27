#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!#define _DBG_ print *,"file: ",__FILE__," line: ",__LINE__," ithr: ",hybrid%ithr
#define _DBG_

! number of points to use in numerical quadrature to compute control volume
! area in older, mass-centered control volumes for ASC's FV algorithm.
! Code that computes control volumes with the same area as the GLL weights
! (for SCRIP) uses analytic area formula.
! Looks like 6 gives enough precision in general
! for very low res that might not be the case...
#define _GLLPTS_ 6

module surfaces_mod
  use kinds, only : real_kind,longdouble_kind
  use coordinate_systems_mod, only : cartesian3d_t, cartesian2d_t, spherical_polar_t, change_coordinates
  use coordinate_systems_mod, only : cubedsphere2cart, cart2cubedsphere, cubedsphere2cart, &
                                     cube_face_number_from_cart, distance, sphere_tri_area

  use parallel_mod, only   : abortmp,  global_shared_buf, global_shared_sum
  use edge_mod, only       : EdgeBuffer_t, Ghostbuffer3d_t
  use dimensions_mod, only : np, ne, nelemd, max_elements_attached_to_node, s_nv
  use global_norms_mod, only: wrap_repro_sum
#if defined (_AIX) || defined (_BGL) || defined (_BGP) 
  use ieee_arithmetic, only: isnan => ieee_is_nan ! _EXTERNAL
#endif

  implicit none

  private

  !If this cade will be used with refined meshes, need to set to 14 (note: ok to use
  !14 for regular meshes - just wastes memory) - asking christoph about if this code
  !is used with refined meshes
  !TO DO: change code to dynamically allocate the right amount of memory
  
  !integer, parameter, public ::nv = 14 !2*max_elements_attached_to_node (refined mesh case)
  integer, parameter, public ::nv = 6 !2*max_elements_attached_to_node (regular case)

  type, public :: ctrlvol_t
     real (kind=real_kind)    :: vol(np,np)           ! area of the unit sphere covered (local)
     real (kind=real_kind)    :: totvol(np,np)        ! area of the unit sphere covered (local)
     real (kind=real_kind)    :: invvol(np,np)        ! inverse area (includes neigbors)
     type (cartesian3D_t)     :: vert(nv,np,np)       ! bounding box for the polygon
     type (cartesian2d_t)     :: cartp_dual(0:np,0:np)
     type (cartesian3d_t)     :: cart3d_dual(0:np,0:np)
     type (spherical_polar_t) :: vert_latlon(nv,np,np) ! bounding box for the polygon
     integer                  :: nvert(np,np)          ! number of vertex per polygon
     integer                  :: face_no(nv,np,np)     ! face_no of cv vertex.  0 if on cube edge

  end type ctrlvol_t

  ! three options:
  ! (1) cv_gnominc = .true.
  !     build control volumes out of lines which are
  !     always gnomonic coordinate lines.  results in hexagon control volumes
  !     at cube corners and edges.  control volumes at cube-sphere edges are
  !     non-convex, which breaks SCRIP.
  ! (2) cv_gnomonic = .false.  undef USE_PENTAGONS
  !     iterate to minimize difference between spherical area and GLL weight
  !     iteration will not converge for the 4 cvs in the center of some elements
  !     control volumes are all triangles or squares
  ! (3) cv_gnomonic = .false.  def USE_PENTAGONS
  !     iterate to minimize difference between spherical area and GLL weight
  !     introduce pentagons in the center of each element to make areas agree
  !     control volumes are triangles, squares or pentagons
  logical :: cv_gnomonic = .false.
#define USE_PENTAGONS

  type (ctrlvol_t),    public, allocatable, target  :: cvlist(:)
  type (EdgeBuffer_t), private,save  :: edge1
  type (Ghostbuffer3d_t),save                :: ghost_buf


  ! User interface
  public :: InitControlVolumesData   ! allocates internal data structure
  public :: InitControlVolumes       ! Inits all surfaces: vol,totvol, invvol

  private :: GetVolumeLocal
  private :: GetVolume
  logical, private :: initialized=.false.
contains

  ! Ok elemid is the local element id (in nets:nete)

  function GetVolume(elemid) result(vol)

    integer, intent(in) :: elemid
    real (kind=real_kind), dimension(:,:), pointer :: vol

    if(.not. initialized) call abortmp('Attempt to use volumes prior to initializing')
    vol => cvlist(elemid)%totvol

  end function GetVolume

  function GetVolumeLocal(elemid) result(vol)

    integer, intent(in) :: elemid
    real (kind=real_kind), dimension(:,:), pointer :: vol

    if(.not. initialized) call abortmp('Attempt to use volumes prior to initializing')
    vol => cvlist(elemid)%vol

  end function GetVolumeLocal
  subroutine InitControlVolumesData(nelemd)
    use edge_mod, only :   initedgebuffer, initGhostBuffer3d
    integer, intent(in) :: nelemd
    ! Cannot be done in a threaded region
    allocate(cvlist(nelemd))
    call initedgebuffer(edge1,3)
    call initGhostBuffer3d(ghost_buf, 3, np+1, 1)
  end subroutine InitControlVolumesData

  subroutine VerifyAreas(elem,hybrid,nets,nete)

    use element_mod,        only : element_t
    use hybrid_mod,         only : hybrid_t
    use reduction_mod, only : parallelmax

    integer,              intent(in)          :: nets,nete
    type (element_t),     intent(in), target  :: elem(:)
    type (hybrid_t),      intent(in)          :: hybrid

    integer                  :: i, j, ie, k, kptr, kmax
    real (kind=real_kind)    :: rspheremp(np,np)
    real (kind=real_kind)    :: invvol(np,np)
    real (kind=real_kind)    :: error, max_error, max_invvol, maxrsphere

    
    error = 0
    max_error = 0
    do ie=nets,nete
      rspheremp = elem(ie)%rspheremp
      invvol    = cvlist(ie)%invvol
      do j=1,np
        do i=1,np
          error = 100*ABS(rspheremp(i,j)-invvol(i,j))/invvol(i,j)
          if (max_error.lt.error) then
             max_error  = error
             max_invvol = invvol(i,j)
             maxrsphere = rspheremp(i,j)
          endif
        enddo
      enddo
    enddo
!    print '(A,F16.4 )',"Control Volume Stats: Max error percent:", max_error 
!    print '(A,F16.12)',"                     Value From Element:",1/maxrsphere
!    print '(A,F16.12)',"              Value From Control Volume:",1/max_invvol
    max_error=parallelmax(max_error,hybrid)
    if (hybrid%masterthread) print '(a,f16.4)',&
         "Control volume area vs. gll area: max error (percent):",max_error
         
  end subroutine VerifyAreas


  subroutine InitControlVolumes(elem, hybrid,nets,nete)
    use element_mod,  only : element_t
    use hybrid_mod,   only : hybrid_t
    integer,              intent(in)          :: nets,nete
    type (element_t),     intent(in), target  :: elem(:)
    type (hybrid_t),      intent(in)          :: hybrid
    if (ne == 0 ) then
       call InitControlVolumes_duel(elem, hybrid,nets,nete)
    else
       call InitControlVolumes_gll(elem, hybrid,nets,nete)
       call VerifVolumes(elem, hybrid,nets,nete)
    endif
  end subroutine InitControlVolumes

  subroutine InitControlVolumes_duel(elem, hybrid,nets,nete)
    use bndry_mod,          only : bndry_exchangev
    use edge_mod,           only : edgeVpack, edgeVunpack, freeedgebuffer, freeghostbuffer3d
    use element_mod,        only : element_t, element_var_coordinates, element_var_coordinates3d
    use hybrid_mod,         only : hybrid_t

    use quadrature_mod,     only : quadrature_t, gausslobatto
    use coordinate_systems_mod, only : cube_face_number_from_sphere

    integer,              intent(in)          :: nets,nete
    type (element_t),     intent(in), target  :: elem(:)
    type (hybrid_t),      intent(in)          :: hybrid

    type (quadrature_t)              :: gll_pts
    type (cartesian2d_t)             :: cv_cartp(0:np,0:np)
    type (cartesian3d_t)             :: quad(4),corners3d(4)
    real (kind=longdouble_kind)      :: cv_pts(0:np)
    real (kind=real_kind)            :: rvert
    real (kind=real_kind)            :: test(np,np,1)

    integer                          :: i, j, ie, k, kmax
    logical                          :: Debug=.FALSE.,keep

    gll_pts = gausslobatto(np)
    ! gll points
    cv_pts(0)=-1
    do i=1,np
       cv_pts(i) = cv_pts(i-1) + gll_pts%weights(i)
    enddo
    cv_pts(np)=1
    do i=1,np-1
       if (gll_pts%points(i) > cv_pts(i) .or. cv_pts(i) > gll_pts%points(i+1)) then
          call abortmp("Error: CV and GLL points not interleaved")
       endif
    enddo


    ! intialize local element areas
    test = 0
    do ie=nets,nete
!!$       cv_cartp(0:np,0:np) = element_var_coordinates(elem(ie)%corners, cv_pts)
!!$       do i=0,np
!!$          do j=0,np
!!$             cvlist(ie)%cart3d_dual(i,j) = cubedsphere2cart(cv_cartp(i,j), elem(ie)%FaceNum)
!!$          enddo
!!$       enddo
!       do i=1,4
!          corners3d(i)=cubedsphere2cart(elem(ie)%corners(i),elem(ie)%FaceNum)
!       enddo
       cvlist(ie)%cart3d_dual(0:np,0:np) = element_var_coordinates3D(elem(ie)%corners3D, cv_pts)



       ! compute true area of element and SEM area
       cvlist(ie)%vol=0
       do i=1,np
          do j=1,np
             ! (gnomonic coordinate lines only), more accurate
             quad(1) = cvlist(ie)%cart3d_dual(i-1,j-1)
             quad(2) = cvlist(ie)%cart3d_dual(i,j-1)
             quad(3) = cvlist(ie)%cart3d_dual(i,j)
             quad(4) = cvlist(ie)%cart3d_dual(i-1,j)
             cvlist(ie)%vol(i,j) = surfarea(quad,4)
          enddo
       enddo
       test(:,:,1) = cvlist(ie)%vol(:,:)
       call edgeVpack(edge1,test,1,0,elem(ie)%desc)
    enddo

    _DBG_
    call bndry_exchangeV(hybrid, edge1)
    _DBG_

    test = 0
    do ie=nets,nete
       test(:,:,1) = cvlist(ie)%vol(:,:)
       call edgeVunpack(edge1, test, 1, 0, elem(ie)%desc)
       cvlist(ie)%totvol(:,:) = test(:,:,1) 
       cvlist(ie)%invvol(:,:)=1.0_real_kind/cvlist(ie)%totvol(:,:)
    enddo


    call VerifyAreas(elem,hybrid,nets,nete)

    ! construct the global CV grid and global CV areas from the
    ! local dual grid (cvlist()%cart_dual) and local areas (cvlist()%vol)
    call construct_cv_duel(elem, hybrid, nets, nete)
    ! compute output needed for SCRIP:  lat/lon coordinates, and for the
    ! control volume with only 3 corners, repeat the last point to make a
    ! degenerate quad.
    kmax=0
    do ie=nets,nete
       kmax = MAX(kmax,MAXVAL(cvlist(ie)%nvert))
    enddo
    do ie=nets,nete
       do j=1,np
          do i=1,np
             cvlist(ie)%vert_latlon(:,i,j)%lat=0
             cvlist(ie)%vert_latlon(:,i,j)%lon=0
             k = cvlist(ie)%nvert(i,j)
             cvlist(ie)%vert        (k+1:kmax, i, j) = cvlist(ie)%vert(k,i,j)
             do k=1,kmax
               cvlist(ie)%vert_latlon(k, i, j) = change_coordinates(cvlist(ie)%vert(k, i, j))
               cvlist(ie)%face_no(k, i, j) = cube_face_number_from_sphere(cvlist(ie)%vert_latlon(k, i, j))
            enddo
          enddo
       enddo
    enddo
    ! Release memory
    if(hybrid%masterthread) then
       call freeedgebuffer(edge1)
       call FreeGhostBuffer3D(ghost_buf)
    end if

    initialized=.true.
  end subroutine  InitControlVolumes_duel




function  average(t, n) result(a)
    integer             , intent(in)       :: n
    type (cartesian3d_t), intent(in)       :: t(n)
    type (cartesian3d_t)                   :: a 
    integer                                :: i
    a%x = 0
    a%y = 0
    a%z = 0
    do i=1,n
      a%x = a%x + t(i)%x
      a%y = a%y + t(i)%y
      a%z = a%z + t(i)%z
    enddo
    a%x = a%x/n
    a%y = a%y/n
    a%z = a%z/n
    return
end function  average

function  make_unique(a, n) result(m)
    use physical_constants,      only : dd_pi
    integer             , intent(in)       :: n
    real (kind=real_kind),intent(inout)    :: a(n) 
    integer                                :: m
    integer                                :: i,j
    do i=1,n-1
      do j=i+1,n
        if (ABS(a(j)-a(i)).lt. 1e-6)  a(j) = 9999
      enddo
    enddo
    m = 0
    do i=1,n
      if (a(i).lt.9000) m = m + 1
    enddo
    if (mod(m,2).ne.0) then
       do i=1,n
          print *,'angle with centroid: ',i,mod(a(i)+dd_pi,dd_pi)
       enddo
       call abortmp("Error: Found an odd number or nodes for cv element. Should be even.")
    endif
    return
end function  make_unique


function SortNodes(t3, n) result(m)
    use coordinate_systems_mod,  only : cube_face_number_from_cart, cart2cubedsphere, change_coordinates
    use physical_constants,      only : dd_pi


    integer             , intent(in)       :: n
    type (cartesian3d_t), intent(inout)    :: t3(n)

    type (cartesian3d_t)                   :: c3, t(n)
    type (cartesian2d_t)                   :: c2, t2
    real (kind=real_kind)                  :: angle(n) 
    integer                                :: i,j,k,m,f
    integer                                :: ip(n)  

    c3 = average(t3, n)
    f  = cube_face_number_from_cart(c3)
    c2 = cart2cubedsphere(c3, f)

    do i=1,n
      t2 = cart2cubedsphere(t3(i), f)
      t2%x = t2%x - c2%x
      t2%y = t2%y - c2%y
      angle(i) = atan2(t2%y, t2%x)
    enddo
    m = make_unique(angle,n)
    do i=1,m
      k = 1
      do j=2,n
        if (angle(j)<angle(k)) k=j
      enddo
      angle(k) = 9999 ! greater than pi
      ip(i)=k
    enddo
    t(1:m) = t3(ip(1:m))
    t3(1:m) = t(1:m)
    return
end function SortNodes




subroutine construct_cv_duel(elem,hybrid,nets,nete)
!
! construct global dual grid from local element dual grid cvlist(ie)%cart3d_dual(:,:)
! all control volumes will be squares or triangles (at cube corners)
!
! 10/2009: added option to make hexagon control volumes at cube edges and corners
!
    use bndry_mod,    only : ghost_exchangev3d
    use element_mod,  only : element_t
    use hybrid_mod,   only : hybrid_t
    use edge_mod,     only : ghostVpack3d, ghostVunpack3d
    use parallel_mod, only : abortmp
    use dimensions_mod, only : max_corner_elem
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    use parallel_mod, only : abortmp

    implicit none
    type (element_t),     intent(in), target    :: elem(:)
    type (hybrid_t),      intent(in) :: hybrid
    integer                          :: nets,nete
!   local
    integer :: i,j,k,m,n,o,p,ie,m2
    real (kind=real_kind)    :: vertpack  (    0:np,       0:np,    3)
    real (kind=real_kind)    :: vertpack2 (    1:np+1,     1:np+1,  3)
    real (kind=real_kind)    :: vertunpack(   -1:np+1,    -1:np+1,  3)
    real (kind=real_kind)    :: vertunpack2(   0:np+2,     0:np+2,  3)
    real (kind=real_kind)    ::         sw(   -1:   0,    -1:   0,  3, max_corner_elem-1)
    real (kind=real_kind)    ::         se(   np:np+1,    -1:   0,  3, max_corner_elem-1)
    real (kind=real_kind)    ::         ne(   np:np+1,    np:np+1,  3, max_corner_elem-1)
    real (kind=real_kind)    ::         nw(   -1:   0,    np:np+1,  3, max_corner_elem-1)
    real (kind=real_kind)    ::        sw2(    0:   1,     0:   1,  3, max_corner_elem-1)
    real (kind=real_kind)    ::        se2( np+1:np+2,     0:   1,  3, max_corner_elem-1)
    real (kind=real_kind)    ::        ne2( np+1:np+2,  np+1:np+2,  3, max_corner_elem-1)
    real (kind=real_kind)    ::        nw2(    0:   1,  np+1:np+2,  3, max_corner_elem-1)

    type (cartesian3d_t)    :: vert(2*nv)
    type (cartesian3d_t)    :: pt_3d
    type (cartesian3d_t)    :: cv   (-1:np+1,   -1:np+1)
    type (cartesian3d_t)    :: cv_sw(-1:   0,   -1:   0,  max_corner_elem-1)
    type (cartesian3d_t)    :: cv_se(np:np+1,   -1:   0,  max_corner_elem-1)
    type (cartesian3d_t)    :: cv_ne(np:np+1,   np:np+1,  max_corner_elem-1)
    type (cartesian3d_t)    :: cv_nw(-1:   0,   np:np+1,  max_corner_elem-1)
    integer                 :: mlt(5:8)


    vertpack   = 0
    vertunpack = 0
    vertpack2  = 0
    vertunpack2= 0
    sw         = 0
    se         = 0
    ne         = 0
    nw         = 0
    do ie=nets,nete
       do j=0,np
          do i=0,np
             pt_3d = cvlist(ie)%cart3d_dual(i,j)
             vertpack(i,j,1) = pt_3d%x
             vertpack(i,j,2) = pt_3d%y
             vertpack(i,j,3) = pt_3d%z
          enddo
       enddo
       do j=0,np
          do i=0,np
             do k=1,3
                vertpack2(i+1,j+1,k) = vertpack(i,j,k)
             enddo
          enddo
       enddo
       call ghostVpack3d(ghost_buf, vertpack2, 3, 0, elem(ie)%desc)
    enddo

    call ghost_exchangev3d (hybrid, ghost_buf)


    do ie=nets,nete
       do j=0,np
          do i=0,np
             pt_3d = cvlist(ie)%cart3d_dual(i,j)
             vertunpack(i,j,1) = pt_3d%x
             vertunpack(i,j,2) = pt_3d%y
             vertunpack(i,j,3) = pt_3d%z
          enddo
       enddo
       do j=0,np
          do i=0,np
             do k=1,3
                vertunpack2(i+1,j+1,k) = vertunpack(i,j,k)
             enddo
          enddo
       enddo
       sw2=0
       se2=0
       nw2=0
       ne2=0
       ! call ghostVunpack3d(ghost_buf, vertunpack2, 3, 0, elem(ie)%desc, sw(-1,-1,:,:),se(np,-1,:,:),nw(-1,np,:,:),ne(np,np,:,:),mlt)
       call ghostVunpack3d(ghost_buf, vertunpack2, 3, 0, elem(ie)%desc, sw2,se2,nw2,ne2,mlt)
       do j=0,np+2
          do i=0,np+2
             do k=1,3
                vertunpack(i-1,j-1,k) = vertunpack2(i,j,k)
             enddo
          enddo
       enddo
       sw=0
       se=0
       nw=0
       ne=0
       do m=1,mlt(swest)-1
         do k=1,3
           do j=0,1    
             do i=0,1    
               sw(i-1,j-1,k,m) = sw2(i,j,k,m)
             enddo
           enddo
         enddo
       enddo
       do m=1,mlt(seast)-1
         do k=1,3
           do j=0,1    
             do i=np+1,np+2
               se(i-1,j-1,k,m) = se2(i,j,k,m)
             enddo
           enddo
         enddo
       enddo
       do m=1,mlt(nwest)-1
         do k=1,3
           do j=np+1,np+2
             do i=0,1
               nw(i-1,j-1,k,m) = nw2(i,j,k,m)
             enddo
           enddo
         enddo
       enddo
       do m=1,mlt(neast)-1
         do k=1,3
           do j=np+1,np+2
             do i=np+1,np+2
               ne(i-1,j-1,k,m) = ne2(i,j,k,m)
             enddo
           enddo
         enddo
       enddo
       ! Count and orient vert array
       ! Positive: 1->2->3->4 is counter clockwise on the sphere
       ! Negative: clockwise orientation

       do j=1,np
          do i=1,np
             cvlist(ie)%vert(:,i,j)%x = 0_real_kind
             cvlist(ie)%vert(:,i,j)%y = 0_real_kind
             cvlist(ie)%vert(:,i,j)%z = 0_real_kind
          enddo
       enddo

       do j=-1,np+1
          do i=-1,np+1
             cv(i,j)%x = vertunpack(i,j,1)
             cv(i,j)%y = vertunpack(i,j,2)
             cv(i,j)%z = vertunpack(i,j,3)
          enddo
       enddo

       do j=-1,0
          do i=-1,0
             do k=1,mlt(swest)-1
                cv_sw(i,j,k)%x = sw(i,j,1,k)
                cv_sw(i,j,k)%y = sw(i,j,2,k)
                cv_sw(i,j,k)%z = sw(i,j,3,k)
             enddo
          enddo
       enddo
       do j=-1,0
          do i=np,np+1
             do k=1,mlt(seast)-1
                cv_se(i,j,k)%x = se(i,j,1,k)
                cv_se(i,j,k)%y = se(i,j,2,k)
                cv_se(i,j,k)%z = se(i,j,3,k)
             enddo
          enddo
       enddo
       do j=np,np+1
          do i=-1,0
             do k=1,mlt(nwest)-1
                cv_nw(i,j,k)%x = nw(i,j,1,k)
                cv_nw(i,j,k)%y = nw(i,j,2,k)
                cv_nw(i,j,k)%z = nw(i,j,3,k)
             enddo
          enddo
       enddo
       do j=np,np+1
          do i=np,np+1
             do k=1,mlt(neast)-1
                cv_ne(i,j,k)%x = ne(i,j,1,k)
                cv_ne(i,j,k)%y = ne(i,j,2,k)
                cv_ne(i,j,k)%z = ne(i,j,3,k)
             enddo
          enddo
       enddo

       do j=2,np-1
          do i=2,np-1
             ! internal vertex on Cubed sphere
             ! Here is the order:
             !
             ! 4NW <- 3NE
             !  |      ^
             !  v      |
             ! 1SW -> 2SE
             vert(1) = cv(i-1, j-1)
             vert(2) = cv(i  , j-1)
             vert(3) = cv(i  , j  )
             vert(4) = cv(i-1, j  )
             cvlist(ie)%vert(1:4,i,j) = vert(1:4)
             cvlist(ie)%nvert(i,j) = 4
          end do
       enddo

       do j=0,np,np
          do i=2,np-1
             vert(1) = cv(i-1, j-1)
             vert(2) = cv(i  , j-1)
             vert(3) = cv(i  , j  )
             vert(4) = cv(i  , j+1)
             vert(5) = cv(i-1, j+1)
             vert(6) = cv(i-1, j  )
             p = j
             if (p.eq.0) p=1
             cvlist(ie)%vert(1:6,i,p) = vert(1:6)
             cvlist(ie)%nvert(i,p) = 6
          enddo
       enddo

       do j=2,np-1
          do i=0,np,np
             vert(1) = cv(i-1, j-1)
             vert(2) = cv(i  , j-1)
             vert(3) = cv(i+1, j-1)
             vert(4) = cv(i+1, j  )
             vert(5) = cv(i  , j  )
             vert(6) = cv(i-1, j  )
             o = i
             if (o.eq.0) o=1
             cvlist(ie)%vert(1:6,o,j) = vert(1:6)
             cvlist(ie)%nvert(o,j) = 6
          enddo
       enddo
       do j=0,np,np
          do i=0,np,np
             m = 0
             vert(:)%x = 0
             vert(:)%y = 0
             vert(:)%z = 0
             if (i.eq.0.and.j.eq.0) then 
                                          ! counterclockwise from lower right
                vert(m+1) = cv(i+1, j-1)  !     5       4    
                vert(m+2) = cv(i+1, j  )  !  (-1,+1)  (0,+1)  (+1,+1)  3
                vert(m+3) = cv(i+1, j+1)  !
                vert(m+4) = cv(i  , j+1)  !  (-1, 0)  (i, j)  (+1, 0)  2
                vert(m+5) = cv(i-1, j+1)  !                  
                vert(m+6) = cv(i-1, j  )  !     X       X     (+1,-1)  1 
                m = m + 6
                if (mlt(swest).ne.0) then
                   vert(m+1) = cv(i-1, j-1)
                   vert(m+2) = cv(i  , j-1)
                   m = m+2
                   do k=1,mlt(swest)-1   ! Bummer, toss in (-1,0) because transpose is undetectable
                      vert(m+1) = cv_sw(i-1, j  , k)
                      vert(m+2) = cv_sw(i-1, j-1, k)
                      vert(m+3) = cv_sw(i  , j-1, k)
                      m=m+3
                   enddo
                endif
             endif
             if (i.eq.np.and.j.eq.0) then
                if (mlt(seast).ne.0) then
                   vert(m+1) = cv(i+1, j-1)
                   vert(m+2) = cv(i+1, j  )
                   m = m+2
                   do k=1,mlt(seast)-1
                      vert(m+1) = cv_se(i  , j-1, k)
                      vert(m+2) = cv_se(i+1, j-1, k)
                      vert(m+3) = cv_se(i+1, j  , k)
                      m=m+3
                   enddo
                endif
                vert(m+1) = cv(i+1, j+1)  
                vert(m+2) = cv(i  , j+1)  
                vert(m+3) = cv(i-1, j+1)  
                vert(m+4) = cv(i-1, j  )   
                vert(m+5) = cv(i-1, j-1)
                vert(m+6) = cv(i  , j-1)
                m = m + 6
             endif
             if (i.eq.np.and.j.eq.np) then
                vert(1) = cv(i+1, j-1)  
                vert(2) = cv(i+1, j  )  
                m = m + 2
                if (mlt(neast).ne.0) then
                   vert(m+1) = cv(i+1, j+1)
                   vert(m+2) = cv(i  , j+1)
                   m = m+2
                   do k=1,mlt(neast)-1
                      vert(m+1) = cv_ne(i+1, j  , k)
                      vert(m+2) = cv_ne(i+1, j+1, k)
                      vert(m+3) = cv_ne(i  , j+1, k)
                      m=m+3
                   enddo
                endif
                vert(m+1) = cv(i-1, j+1)  
                vert(m+2) = cv(i-1, j  )   
                vert(m+3) = cv(i-1, j-1)
                vert(m+4) = cv(i  , j-1)
                m = m + 4
             endif
             if (i.eq.0.and.j.eq.np) then
                vert(m+1) = cv(i+1, j-1)  
                vert(m+2) = cv(i+1, j  )  
                vert(m+3) = cv(i+1, j+1)  
                vert(m+4) = cv(i  , j+1)  
                m = m + 4
                if (mlt(nwest).ne.0) then
                   vert(m+1) = cv(i-1, j+1)
                   vert(m+2) = cv(i-1, j  )
                   m = m+2
                   do k=1,mlt(nwest)-1
                      vert(m+1) = cv_nw(i  , j+1, k)
                      vert(m+2) = cv_nw(i-1, j+1, k)
                      vert(m+3) = cv_nw(i-1, j  , k)
                      m=m+3
                   enddo
                endif
                vert(m+1) = cv(i-1, j-1)
                vert(m+2) = cv(i  , j-1)
                m = m + 2
             endif
             o = i
             p = j
             if (o.eq.0) o=1
             if (p.eq.0) p=1
             m2=m
             if (8.lt.m) m=SortNodes(vert, m2)
             if (m>nv) call abortmp("error: vert dimensioned too small")
             cvlist(ie)%vert(1:m,o,p) = vert(1:m)
             cvlist(ie)%nvert(o,p) = m
          enddo
       enddo
    enddo
  end subroutine construct_cv_duel

  function SurfArea( cv, nvert ) result(area)
  implicit none
  real(kind=real_kind) :: area,area1,area2,area3
  type (cartesian3D_t) :: cv(:)
  integer :: nvert
  if (abs(nvert) == 3 ) then
     area2=0
     if (cv(1)%x==0) then
        call sphere_tri_area(cv(2),cv(3),cv(4),area1)
     else if (cv(2)%x==0) then
        call sphere_tri_area(cv(1),cv(3),cv(4),area1)
     else if (cv(3)%x==0) then
        call sphere_tri_area(cv(1),cv(2),cv(4),area1)
     else if (cv(4)%x==0) then
        call sphere_tri_area(cv(1),cv(2),cv(3),area1)
     else
        print *,cv(1)%x,cv(1)%y
        print *,cv(2)%x,cv(2)%y
        print *,cv(3)%x,cv(3)%y
        print *,cv(4)%x,cv(4)%y
        print *,'surfarea error: should never happen'
        stop
     endif
  else if (abs(nvert)==4) then
     call sphere_tri_area(cv(1),cv(2),cv(3),area1)
     call sphere_tri_area(cv(1),cv(3),cv(4),area2)
     area3=0

  else if (abs(nvert)==5) then
     call sphere_tri_area(cv(1),cv(2),cv(3),area1)
     call sphere_tri_area(cv(1),cv(3),cv(4),area2)
     call sphere_tri_area(cv(1),cv(4),cv(5),area3)
  else
     stop 'error: surfarea nvert>5 not yet coded'
  endif
  area=area1+area2+area3
  end function SurfArea
  










  !   ^
  !   |dy  o
  !   |
  ! (x,y) ---->dx
  function  SurfArea_dxdy(dx,dy,corner) result(integral)
    use quadrature_mod, only : quadrature_t
    use physical_constants, only : dd_pi

    type (cartesian2d_t)  :: corner
    real (kind=real_kind) :: dx,dy,da,integral,x,y

    integer               :: i,j


    ! cubed-sphere cell area, from Lauritzen & Nair MWR 2008
    ! central angles:
    ! cube face: -pi/4,-pi/4 -> pi/4,pi/4
    ! this formula gives 2   so normalize by 4pi/6 / 2 = pi/3
    real (kind=real_kind) :: alpha,beta,a1,a2,a3,a4
    alpha = corner%x ; beta  = corner%y
    a1 = acos(-sin(alpha)*sin(beta))             !  2.094
    a2 =-acos(-sin(alpha+dx)*sin(beta) )         ! -1.047
    a3 =-acos(-sin(alpha)*sin(beta+dy) )         ! -1.047
    a4 = acos(-sin(alpha+dx)*sin(beta+dy) )      !  2.094 
    integral = (a1+a2+a3+a4)
    return
  end function SurfArea_dxdy



  function find_intersect( x1in,x2in,y1in,y2in ) result(sect)
  implicit none
  type (cartesian2D_t) :: x1in,x2in,y1in,y2in,sect,x,y,b
  type (cartesian2D_t) :: x1,x2,y1,y2
!  integer :: nvert
  real(kind=real_kind) :: s1,s2,detA

!  x1 + (x2-x1)*s1  = y1 + (y2-y1)*s2
! b = y1-x1
! x=x2-x1
! y=y2-y1
!  x s1 - y s2 = b
!  x(1) s1 - y(1) s2 = b(1)
!  x(2) s1 - y(2) s2 = b(2)
!
!  x(1) -y(1)   s1   =  b(1)        A s = b
!  x(2) -y(2)   s2   =  b(2)
!
!  A2=  -y(2)  y(1)  
!       -x(2)  x(1)                s = A2 * b /detA

  ! convert to gnomonic
  x1%x = tan(x1in%x)
  x2%x = tan(x2in%x)
  y1%x = tan(y1in%x)
  y2%x = tan(y2in%x)
  x1%y = tan(x1in%y)
  x2%y = tan(x2in%y)
  y1%y = tan(y1in%y)
  y2%y = tan(y2in%y)

  x%x=x2%x-x1%x
  x%y=x2%y-x1%y
  y%x=y2%x-y1%x
  y%y=y2%y-y1%y
  b%x=y1%x-x1%x 
  b%y=y1%y-x1%y 

  detA = -x%x*y%y + x%y*y%x
   
  s1 =  (-y%y*b%x + y%x*b%y )/detA
  s2 =  (-x%y*b%x + x%x*b%y )/detA

  sect%x = x1%x + (x2%x-x1%x)*s1
  sect%y = x1%y + (x2%y-x1%y)*s1

  sect%x = (sect%x + y1%x + (y2%x-y1%x)*s2)/2
  sect%y = (sect%y + y1%y + (y2%y-y1%y)*s2)/2

  if (s1<0 .or. s1>1) then
     print *,'failed: intersection: ',s1,s2
     stop
  endif

  ! convert back to equal angle:
  sect%x = atan(sect%x)
  sect%y = atan(sect%y)
  end function find_intersect





   subroutine pentagon_iteration(sq1,sq2,pent,asq1,asq2,apent,faceno,anorm)
!               sq2 
!              4  3                
!              1  2
!
!    sq1       4  3 
!    2  1    5  pent  
!    3  4    1    2
!
!
!   d/dt sq1(1)  =    (area(sq1)-asq1) * [ com(sq1)-sq1(1) ]  
!                    +(area(sq2)-asq2) * [ com(sq2)-sq1(1) ]  
!                    +(area(pent)-apent) * [ com(pent)-sq1(1) ]  
! 
!  
!   
    type (cartesian2d_t)             :: sq1(4),sq2(4),pent(5)
    real (kind=real_kind)            :: asq1,asq2,apent,anorm
    integer :: faceno

    type (cartesian3D_t)             :: sq1_3d(4),sq2_3d(4),pent_3d(5)
    real (kind=real_kind)            :: isq1,isq2,ipent,dt,diff1,diff2,diffp,err
    real (kind=real_kind)            :: tol_pentagon_iteration=1e-10
    type (cartesian2d_t)             :: sq1com,sq2com,pentcom,ds1,ds2
    integer :: i,iter,iter_max
    
    ! compute center of mass:
    sq1com%x = sum(sq1(:)%x)/4
    sq1com%y = sum(sq1(:)%y)/4
    sq2com%x = sum(sq2(:)%x)/4
    sq2com%y = sum(sq2(:)%y)/4
    pentcom%x = sum(pent(:)%x)/5
    pentcom%y = sum(pent(:)%y)/5

    do i=1,4
       sq1_3d(i)=cubedsphere2cart(sq1(i),faceno  )
       sq2_3d(i)=cubedsphere2cart(sq2(i),faceno  )
       pent_3d(i)=cubedsphere2cart(pent(i),faceno  )
    enddo
    pent_3d(5)=cubedsphere2cart(pent(5),faceno  )

    dt=.5
    iter_max = 10000
    do iter=1,iter_max
       isq1 = SurfArea(sq1_3d,4)
       isq2 = SurfArea(sq2_3d,4)
       ipent = SurfArea(pent_3d,5)

!   d/dt sq1(1)  =    (area(sq1)-asq1) * [ com(sq1)-sq1(1) ]  
!                    +(area(sq2)-asq2) * [ com(sq2)-sq1(1) ]  
!                    +(area(pent)-apent) * [ com(pent)-sq1(1) ]  
! 
       diff1 = (isq1-asq1)/anorm
       diff2 = (isq2-asq2)/anorm
       diffp = (ipent-apent)/anorm

       err = abs(diff1)+abs(diff2)+abs(diffp)
       if (err< tol_pentagon_iteration) exit
       !write(*,'(i5,3e18.5)') iter,diff1,diff2,diffp
       if (mod(iter,1000).eq.0) write(*,'(i5,3e18.5)') iter,err

       ds1%x = diff1* ( sq1com%x - sq1(1)%x )
       ds1%y = diff1* ( sq1com%y - sq1(1)%y )
       ds1%x = ds1%x + diffp* ( pentcom%x - sq1(1)%x )
       ds1%y = ds1%y + diffp* ( pentcom%y - sq1(1)%y )

       ds2%x = diff2* ( sq2com%x - sq2(1)%x )
       ds2%y = diff2* ( sq2com%y - sq2(1)%y )
       ds2%x = ds2%x + diffp* ( pentcom%x - sq2(1)%x )
       ds2%y = ds2%y + diffp* ( pentcom%y - sq2(1)%y )

       sq1(1)%x = sq1(1)%x + dt*ds1%x
       sq1(1)%y = sq1(1)%y + dt*ds1%y
       sq2(1)%x = sq2(1)%x + dt*ds2%x
       sq2(1)%y = sq2(1)%y + dt*ds2%y
       pent(4)=sq2(1)
       pent(5)=sq1(1)
       sq1_3d(1)=cubedsphere2cart(sq1(1),faceno  )
       sq2_3d(1)=cubedsphere2cart(sq2(1),faceno  )
       pent_3d(4)=sq2_3d(1)
       pent_3d(5)=sq1_3d(1)
    enddo
    if (iter >= iter_max) then
       print *,'pentagon iteration did not converge err=',err
    endif
   end subroutine pentagon_iteration






  



  subroutine InitControlVolumes_gll(elem, hybrid,nets,nete)    
    use bndry_mod, only : bndry_exchangev
    use edge_mod, only : freeedgebuffer
    use element_mod, only : element_t,element_coordinates
    use hybrid_mod, only : hybrid_t

    use quadrature_mod, only : quadrature_t, gausslobatto
    use dimensions_mod, only : nlev
    use physical_constants, only : dd_pi
    use cube_mod, only : convert_gbl_index
    use reduction_mod, only : parallelmax

    integer,              intent(in)    :: nets,nete
    type (element_t),     intent(in), target    :: elem(:)
    type (hybrid_t),      intent(in) :: hybrid

    type (cartesian2d_t)             :: cartp_com(np,np)  ! center of mass
    type (cartesian2d_t)             :: cartp_nm1(0:np,0:np)
    real (kind=real_kind)            :: delx_k,dely_k,sum_dbg,r
    integer                          :: i,j,ie,k,kptr,gllpts,nvert,k2,ie1,je1,face_no,kinsert
    integer                          ::  iter,iter_max,i1,j1
    real (kind=real_kind)            :: diff(np,np),diffy(np-1,np-1),diffx(np-1,np-1)
    real (kind=real_kind)            :: dx,dy,a1(nets:nete),a2(nets:nete),d1(nets:nete),d1mid(nets:nete)
    real (kind=real_kind)            :: d2,d1_global,d1_global_mid,sphere1,sphere2,diff2,diff3
    real (kind=real_kind)            :: diff23,diff32,diff33,diff22
    real (kind=longdouble_kind)      :: gllnm1(0:np)
    type (cartesian2d_t)             :: corner,start,end,cv_loc_2d(4,np,np),cvnew_loc_2d(4,np,np)
    type (cartesian3D_t)             :: cart,cv_loc_3d(nv,np,np)
    type (cartesian3D_t)             :: temp3d(nv)
    type (cartesian2d_t)             :: cartp2d(np,np)
    type (cartesian2d_t)             :: x1,x2,x3,x
    type (cartesian2d_t)             :: sq1(4),sq2(4),pent(5)
    type (cartesian3D_t)             :: x1_3d,x2_3d,x3_3d
    type (quadrature_t)              :: gll
    type (cartesian2d_t)             :: dir,dirsum
    type (spherical_polar_t)         :: polar_tmp(0:np,0:np)
    real (kind=real_kind)            :: rvert,area1,area2,ave,lat(4),lon(4)
    real (kind=real_kind)            :: s,ds,triarea,triarea_target
    real (kind=real_kind)            :: xp1,xm1,yp1,ym1,sumdiff
    real (kind=real_kind)            :: tiny=1e-11,norm
    real (kind=real_kind)            :: tol=2e-11  ! convergece outer iteration
    real (kind=real_kind)            :: tol_pentagons=1e-13  ! convergece pentagon iteration

    ! area difference to trigger pentagons.  
    ! if it is too small, we will have pentagons with 1 very short edges
    ! accuracy of surfarea() with very thin triangles seems poor (1e-11) 
    ! ne=30  1e-3:  add 648 pentagons.  area ratio:  1.003
    ! ne=30  1e-4:  add 696 pentagons.  area ratio:  1.000004102
    ! ne=30  1e-5:  add 696 pentagons.  area ratio:  1.000004102
    ! ne=240 1e-4:  add 5688/ 345600 pentagons, area ratio: 1.0004
    ! ne=240 1e-5:  add 5736/ 345600 pentagons, area ratio: 1.000000078
    real (kind=real_kind)            :: tol_use_pentagons=1e-5
    logical                          :: Debug=.FALSE.,keep

    integer :: face1,face2,found,ie_max,movex,movey,moved,ii,kmax
    integer :: nskip,npent
    integer :: nskipie(nets:nete), npentie(nets:nete)
    type (cartesian2d_t)             :: vert1_2d, vert_2d,vert2_2d
    type (cartesian3D_t)             :: vert1,vert2,vert_inserted(7)


    kmax=4


    gll = gausslobatto(np)
    ! mid point rule: 
    do i=1,np-1
       gllnm1(i) = ( gll%points(i) + gll%points(i+1) ) /2
    enddo
    ! check that gll(i) < gllnm1(i) < gll(i+1)
    do i=1,np-1
       if (gll%points(i) > gllnm1(i) .or. gllnm1(i) > gll%points(i+1)) then
          call abortmp("Error: CV and GLL points not interleaved") 
       endif
    enddo
    gllnm1(0)=-1
    gllnm1(np)=1

    ! MNL: dx and dy are no longer part of element_t
    !      but they are easily computed for the
    !      uniform case
    dx = DD_PI/(2.0d0*dble(ne))
    dy = dx

    ! intialize local element dual grid, local element areas

    do ie=nets,nete

       call convert_gbl_index(elem(ie)%vertex%number,ie1,je1,face_no)
       start%x=-DD_PI/4 + ie1*dx
       start%y=-DD_PI/4 + je1*dy
       end%x  =start%x + dx
       end%y  =start%y + dy
       cartp_nm1(0:np,0:np) = element_coordinates(start,end,gllnm1)
       cvlist(ie)%cartp_dual = cartp_nm1  

       ! compute true area of element and SEM area
       a1(ie) = SurfArea_dxdy(dx,dy,elem(ie)%cartp(1,1))
       a2(ie) = sum(elem(ie)%spheremp(:,:))
       do i=1,np
          do j=1,np
             ! (gnomonic coordinate lines only), more accurate
             delx_k = cartp_nm1(i,j-1)%x - cartp_nm1(i-1,j-1)%x
             dely_k = cartp_nm1(i-1,j)%y - cartp_nm1(i-1,j-1)%y
             cvlist(ie)%vol(i,j) = SurfArea_dxdy(delx_k,dely_k,cartp_nm1(i-1,j-1))
          enddo
       enddo
       global_shared_buf(ie,1) = a1(ie)
       global_shared_buf(ie,2) = a2(ie)
    enddo
    call wrap_repro_sum(nvars=2, comm=hybrid%par%comm)
    sphere1 = global_shared_sum(1)
    sphere2 = global_shared_sum(2)

    ! construct the global CV grid and global CV areas from the 
    ! local dual grid (cvlist()%cart_dual) and local areas (cvlist()%vol)
    call construct_cv_gll(elem,hybrid,nets,nete)

    iter_max=2000
    if (iter_max>0) then
       ! areas computed from eleemnts on boundaries are from hexagons and pentagons
       ! compute new areas where all CVs are squares or triangles
       do ie=nets,nete
          do i=1,np
             do j=1,np
                ! ifort bug if we try this: 
                ! area2 = surfarea(cvlist(ie)%vert(1:4,i,j),cvlist(ie)%nvert(i,j))
                cv_loc_3d(:,i,j)=cvlist(ie)%vert(:,i,j)
                area2 = surfarea(cv_loc_3d(:,i,j),cvlist(ie)%nvert(i,j))
                cvlist(ie)%totvol(i,j)=area2
             enddo
          enddo
       enddo
    endif
    ! iteration over cvlist(ie)%totvol
    d1_global=0
    do iter=1,iter_max
       ie_max=-1
       do ie=nets,nete
          ! we want at each point, the gll_area = true_area
          ! but sum(gll_area) = a2 and sum(true_area)=a1
          ! so normalize so that: gll_area/a2 = true_area/a1, or gll_area = area*a2/a1
!          diff(:,:) = ( cvlist(ie)%totvol(:,:)*sphere2/sphere1 - 1/elem(ie)%rspheremp(:,:) )
!          sumdiff=sum(diff)
!          diff(:,:) = diff(:,:)/(a1(ie)/(np*np))

          ! requires more iterations, but the total volume within an 
          ! element is always correct
          diff(:,:) = ( cvlist(ie)%vol(:,:) - elem(ie)%spheremp(:,:)*a1(ie)/a2(ie) )
          sumdiff=sum( cvlist(ie)%vol(:,:)) - a1(ie)
          diff(:,:) = diff(:,:)/(a1(ie)/(np*np))



          ! set boundary values (actually not used)
          cartp_nm1 = cvlist(ie)%cartp_dual(0:np,0:np)
          ! convert 9 cv corners in this element into cart_nm1 cubed-sphere coordiantes
          do i=1,np-1
             do j=1,np-1
                cartp_nm1(i,j) = cart2cubedsphere( cvlist(ie)%vert(3,i,j),elem(ie)%FaceNum  )
             enddo
          enddo
          ! compute center of mass of control volumes:
          ! todo: move points towards GLL points, not center of mass
          ! center of mass could send up a feedback with CV points!
          do i=1,np
             do j=1,np
                cart%x = sum( cvlist(ie)%vert(:,i,j)%x )/abs(cvlist(ie)%nvert(i,j))
                cart%y = sum( cvlist(ie)%vert(:,i,j)%y )/abs(cvlist(ie)%nvert(i,j))
                cart%z = sum( cvlist(ie)%vert(:,i,j)%z )/abs(cvlist(ie)%nvert(i,j))
                cartp_com(i,j) = cart2cubedsphere( cart,elem(ie)%FaceNum  )
             enddo
          enddo
          d2=0
          do i=1,np-1
             do j=1,np-1
                dirsum%x=0
                dirsum%y=0
                movex=1
                movey=1
                moved=0

                do i1=0,1
                do j1=0,1
                    ! keep=.true. : .85/1.05 
                    ! corners only: .93/1.07
                    ! corners and edges:  .89/1.11
#ifdef USE_PENTAGONS
                   keep=.false.
                   ! corner volumes 
                   if (i==1 .and. j==1) then  
                      if (i1==0 .and. j1==0) keep=.true.
                      moved=1
                   else if (i==np-1 .and. j==1) then
                      if (i1==1 .and. j1==0) keep=.true.
                      moved=-1
                   else if (i==1 .and. j==np-1) then
                      if (i1==0 .and. j1==1) keep=.true.
                      moved=-1
                   else if (i==np-1 .and. j==np-1) then
                      if (i1==1 .and. j1==1) keep=.true.
                      moved=1
                   ! edge volumes


                   else if (i==1) then
                      if (i1==0) keep=.true.
                   else if (i==np-1) then 
                      if (i1==1) keep=.true.
                   else if (j==1) then
                      if (j1==0) keep=.true.
                   else if (j==np-1) then
                      if (j1==1) keep=.true.
                   else
                      keep=.true.
                   endif
#else                      
                   keep=.false.
                   ! corner volumes 
                   if (i==1 .and. j==1) then  
                      if (i1==0 .and. j1==0) keep=.true.
                      moved=1
                   else if (i==np-1 .and. j==1) then
                      if (i1==1 .and. j1==0) keep=.true.
                      moved=-1
                   else if (i==1 .and. j==np-1) then
                      if (i1==0 .and. j1==1) keep=.true.
                      moved=-1
                   else if (i==np-1 .and. j==np-1) then
                      if (i1==1 .and. j1==1) keep=.true.
                      moved=1
                   ! edge volumes

                   else if (i==1) then
                      if (i1==0) keep=.true.
                      keep=.true.
                   else if (i==np-1) then 
                      if (i1==1) keep=.true.
                      keep=.true.
                   else if (j==1) then
                      if (j1==0) keep=.true.
                      keep=.true.
                   else if (j==np-1) then
                      if (j1==1) keep=.true.
                      keep=.true.
                   else
                      !if (j1==0) keep=.true.
                      if (i1==0 .and. j1==0) keep=.true.
                      keep=.true.
                   endif
#endif
                   if (keep) then
                      ! error weighted direction towards center of mass of area
                      !dir%x =  (cartp_com(i+i1,j+j1)%x - cartp_nm1(i,j)%x )*(abs(diff(i+i1,j+j1)))
                      !dir%y =  (cartp_com(i+i1,j+j1)%y - cartp_nm1(i,j)%y )*(abs(diff(i+i1,j+j1)))
                      ! move towards grid point
                      dir%x =  (elem(ie)%cartp(i+i1,j+j1)%x - cartp_nm1(i,j)%x )*(abs(diff(i+i1,j+j1)))
                      dir%y =  (elem(ie)%cartp(i+i1,j+j1)%y - cartp_nm1(i,j)%y )*(abs(diff(i+i1,j+j1)))
                      if (moved==1) then
                         ! project onto (1,1)/sqrt(2)
                         dir%x = dir%x/sqrt(2d0) + dir%y/sqrt(2d0)
                         dir%y = dir%x
                      endif
                      if (moved==-1) then
                         ! project onto (-1,1)/sqrt(2)
                         dir%y = -dir%x/sqrt(2d0) + dir%y/sqrt(2d0)
                         dir%x = -dir%y
                      endif
                      

                      if ( diff(i+i1,j+j1) > 0 ) then
                         ! this volume is too big, so move cv point towards grid center
                         ! weighted by length error
                         dirsum%x = dirsum%x + movex*dir%x
                         dirsum%y = dirsum%y + movey*dir%y
                      else
                         dirsum%x = dirsum%x - movex*dir%x
                         dirsum%y = dirsum%y - movey*dir%y
                      endif
                   endif
                enddo
                enddo
                d2 = d2 + dirsum%x**2 + dirsum%y**2
                cartp_nm1(i,j)%x = cartp_nm1(i,j)%x + 0.25*dirsum%x
                cartp_nm1(i,j)%y = cartp_nm1(i,j)%y + 0.25*dirsum%y

             enddo
          enddo
          cvlist(ie)%cartp_dual(0:np,0:np) = cartp_nm1  
          d2=sqrt(d2) 


#if 0
          if (hybrid%masterthrea) then
          if (mod(iter-1,250).eq.0  ) then
             if (d1_global == d1(ie) ) then  ! max error, last iteration
                ie_max=ie
                !print *,"max element error grad=",d2," ie=",ie
                print *,iter,"sum diff",sumdiff
                do i=4,1,-1
                   write(*,'(a,4e13.4)') 'diff(:,i)=',diff(:,i)
                enddo
             endif
          endif
          endif
#endif
          d1(ie)=sqrt(sum(diff**2))

          d1mid(ie)=d1(ie)
#ifdef USE_PENTAGONS
          ! ignore center cv's:
          diff(2:3,2:3)=0
          d1mid(ie)=sqrt(sum(diff**2))
#endif


       enddo  ! ie loop
       dx=maxval(d1)
       d1_global = ParallelMax(dx,hybrid)
       dx=maxval(d1mid)
       d1_global_mid = ParallelMax(dx,hybrid)
       if (mod(iter-1,250).eq.0) then
          if (hybrid%masterthread) print *,iter,"max d1=",d1_global,d1_global_mid
       endif
       ! compute new global CV  (cvlist(ie)%vert from cvlist(ie)%cartp_dual).  
       ! cvlist()%totarea incorrect since local volumes not computed above
       call construct_cv_gll(elem,hybrid,nets,nete)

       ! update totvol (area of multi-element cv)
       do ie=nets,nete
          do i=1,np
             do j=1,np
                ! ifort bug if we try this: 
                ! area2 = surfarea(cvlist(ie)%vert(1:4,i,j),cvlist(ie)%nvert(i,j))
                cv_loc_3d(:,i,j)=cvlist(ie)%vert(:,i,j)
                area2 = surfarea(cv_loc_3d(:,i,j),cvlist(ie)%nvert(i,j))
                cvlist(ie)%totvol(i,j) = area2
                if (isnan(area2)) then
                   print *,'ie,i,j',ie,i,j
                   print *,cvlist(ie)%nvert(i,j)
                   print *,cv_loc_3d(1,i,j)
                   print *,cv_loc_3d(2,i,j)
                   print *,cv_loc_3d(3,i,j)
                   print *,cv_loc_3d(4,i,j)
                   stop 'area = NaN'
                endif
             enddo
          enddo
       enddo



       ! update %vol (local control volume within each element)
       do ie=nets,nete
          cartp2d = elem(ie)%cartp
          do i=1,np
             do j=1,np
                ! ifort bug if we try this: 
                ! area2 = surfarea(cvlist(ie)%vert(1:4,i,j),cvlist(ie)%nvert(i,j))
                do ii=1,4
                   cv_loc_2d(ii,i,j) = cart2cubedsphere( cvlist(ie)%vert(ii,i,j),elem(ie)%FaceNum  )
                enddo
                
                !print *,'intersection loop',i,j
                if (i==1 .and. j==1) then
                   cv_loc_2d(1,i,j)=cartp2d(i,j)
                endif
                if (i==np .and. j==1) then
                   cv_loc_2d(2,i,j)=cartp2d(i,j)
                endif
                if (i==1 .and. j==np) then
                   cv_loc_2d(4,i,j)=cartp2d(i,j)
                endif
                if (i==np .and. j==np) then
                   cv_loc_2d(3,i,j)=cartp2d(i,j)
                endif
                cvnew_loc_2d(:,i,j)=cv_loc_2d(:,i,j)
                
                !
                ! 4NW <- 3NE
                !  |      ^
                !  v      |
                ! 1SW -> 2SE
                if (i==1) then
                   ! replace points with x< elem(ie)%vert(i,j)%x
                   if (cv_loc_2d(1,i,j)%x < cartp2d(i,j)%x) then
                      cvnew_loc_2d(1,i,j) = find_intersect(&
                           cv_loc_2d(1,i,j), cv_loc_2d(2,i,j),&
                           elem(ie)%cartp(i,1),elem(ie)%cartp(i,np))
                   endif
                   if (cv_loc_2d(4,i,j)%x < cartp2d(i,j)%x) then
                      cvnew_loc_2d(4,i,j) = find_intersect(&
                           cv_loc_2d(4,i,j), cv_loc_2d(3,i,j),&
                           elem(ie)%cartp(i,1),elem(ie)%cartp(i,np))
                   endif
                endif

                if (i==np) then
                   ! replace points with x> elem(ie)%vert(i,j)%x
                   if (cv_loc_2d(2,i,j)%x > cartp2d(i,j)%x) then
                      cvnew_loc_2d(2,i,j) = find_intersect(&
                           cv_loc_2d(1,i,j), cv_loc_2d(2,i,j),&
                           elem(ie)%cartp(i,1),elem(ie)%cartp(i,np))
                   endif
                   if (cv_loc_2d(3,i,j)%x > cartp2d(i,j)%x) then
                      cvnew_loc_2d(3,i,j) = find_intersect(&
                           cv_loc_2d(4,i,j), cv_loc_2d(3,i,j),&
                           elem(ie)%cartp(i,1),elem(ie)%cartp(i,np))
                   endif
                endif
                !
                ! 4NW <- 3NE
                !  |      ^
                !  v      |
                ! 1SW -> 2SE
                if (j==1) then
                   ! replace points with y < elem(ie)%vert(i,j)%y
                   if (cv_loc_2d(1,i,j)%y < cartp2d(i,j)%y) then
                      cvnew_loc_2d(1,i,j) = find_intersect(&
                           cv_loc_2d(1,i,j), cv_loc_2d(4,i,j),&
                           elem(ie)%cartp(1,j),elem(ie)%cartp(np,j))
                   endif
                   if (cv_loc_2d(2,i,j)%y < cartp2d(i,j)%y) then
                      cvnew_loc_2d(2,i,j) = find_intersect(&
                           cv_loc_2d(2,i,j), cv_loc_2d(3,i,j),&
                           elem(ie)%cartp(1,j),elem(ie)%cartp(np,j))
                   endif
                endif
                if (j==np) then
                   ! replace points with y > elem(ie)%vert(i,j)%y
                   if (cv_loc_2d(4,i,j)%y > cartp2d(i,j)%y) then
                      cvnew_loc_2d(4,i,j) = find_intersect(&
                           cv_loc_2d(1,i,j), cv_loc_2d(4,i,j),&
                           elem(ie)%cartp(1,j),elem(ie)%cartp(np,j))
                   endif
                   if (cv_loc_2d(3,i,j)%y > cartp2d(i,j)%y) then
                      cvnew_loc_2d(3,i,j) = find_intersect(&
                           cv_loc_2d(2,i,j), cv_loc_2d(3,i,j),&
                           elem(ie)%cartp(1,j),elem(ie)%cartp(np,j))
                   endif
                endif
                do ii=1,4
                   cv_loc_3d(ii,i,j)=cubedsphere2cart(cvnew_loc_2d(ii,i,j),elem(ie)%FaceNum  )
                enddo
                area2 = surfarea(cv_loc_3d(:,i,j),4)
                cvlist(ie)%vol(i,j) = area2
                if (isnan(area2)) then
                   print *,'ie,i,j',ie,i,j
                   print *,cvlist(ie)%nvert(i,j)
                   print *,cv_loc_3d(1,i,j)
                   print *,cv_loc_3d(2,i,j)
                   print *,cv_loc_3d(3,i,j)
                   print *,cv_loc_3d(4,i,j)
                   stop 'area = NaN'
                endif
             enddo
          enddo
       enddo
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

       
       if ( d1_global > 10.0 .or. d1_global_mid < tol) then
          if (hybrid%masterthread) print *,'first iteration stopping:'
          if (hybrid%masterthread) print *,iter,"max error=",d1_global_mid
          exit
       endif
    enddo  ! iteration loop


#ifdef USE_PENTAGONS
    kmax=5

    nskip=0
    npent=0
    nskipie(:) = 0
    npentie(:) = 0
    do ie=nets,nete
    diff = ( cvlist(ie)%vol(:,:) - elem(ie)%spheremp(:,:)*a1(ie)/a2(ie) )
    if ( maxval(abs(diff(2:3,2:3)))/a1(ie)  > tol_use_pentagons ) then
       npent=npent+1
       npentie(ie) = npentie(ie) + 1
                !
                ! 4NW <- 3NE
                !  |      ^
                !  v      |             23   33
                ! 1SW -> 2SE            22   32
       if (diff(2,2)>0 .and. diff(3,3)>0) then
          x1 = cart2cubedsphere( cvlist(ie)%vert(3,2,2),elem(ie)%FaceNum  )
          x2 = cart2cubedsphere( cvlist(ie)%vert(1,2,2),elem(ie)%FaceNum  )
          s=.99
          x3%x = x2%x + (x1%x-x2%x)*s
          x3%y = x2%y + (x1%y-x2%y)*s

          sq1(1) = x3
          sq1(2) = cart2cubedsphere( cvlist(ie)%vert(4,2,2),elem(ie)%FaceNum  ) 
          sq1(3) = cart2cubedsphere( cvlist(ie)%vert(1,2,2),elem(ie)%FaceNum  )
          sq1(4) = cart2cubedsphere( cvlist(ie)%vert(2,2,2),elem(ie)%FaceNum  )

          x2 = cart2cubedsphere( cvlist(ie)%vert(3,3,3),elem(ie)%FaceNum  )          
          s=.99
          x3%x = x2%x + (x1%x-x2%x)*s
          x3%y = x2%y + (x1%y-x2%y)*s

          sq2(1) = x3
          sq2(2) = cart2cubedsphere( cvlist(ie)%vert(2,3,3),elem(ie)%FaceNum  )
          sq2(3) = cart2cubedsphere( cvlist(ie)%vert(3,3,3),elem(ie)%FaceNum  )
          sq2(4) = cart2cubedsphere( cvlist(ie)%vert(4,3,3),elem(ie)%FaceNum  )

          pent(1) = cart2cubedsphere( cvlist(ie)%vert(1,3,2),elem(ie)%FaceNum  )
          pent(2) = cart2cubedsphere( cvlist(ie)%vert(2,3,2),elem(ie)%FaceNum  )
          pent(3) = cart2cubedsphere( cvlist(ie)%vert(3,3,2),elem(ie)%FaceNum  )
          pent(4) = sq2(1)
          pent(5) = sq1(1)

          call pentagon_iteration(sq1,sq2,pent,&
               elem(ie)%spheremp(2,2)*a1(ie)/a2(ie), &
               elem(ie)%spheremp(3,3)*a1(ie)/a2(ie), &
               elem(ie)%spheremp(3,2)*a1(ie)/a2(ie),elem(ie)%FaceNum,a1(ie))

          x2_3d=cubedsphere2cart(sq1(1),elem(ie)%FaceNum  )
          x3_3d=cubedsphere2cart(sq2(1),elem(ie)%FaceNum  )

          cvlist(ie)%vert(3,2,2)=x2_3d
          cvlist(ie)%vert(1,3,3)=x3_3d

          cvlist(ie)%vert(5,2,3)=cvlist(ie)%vert(4,2,3)
          cvlist(ie)%vert(4,2,3)=cvlist(ie)%vert(3,2,3)
          cvlist(ie)%vert(2,2,3)=x2_3d
          cvlist(ie)%vert(3,2,3)=x3_3d

          cvlist(ie)%vert(5,3,2)=x2_3d
          cvlist(ie)%vert(4,3,2)=x3_3d

          cvlist(ie)%nvert(2,3)=sign(5,cvlist(ie)%nvert(2,3))
          cvlist(ie)%nvert(3,2)=sign(5,cvlist(ie)%nvert(3,2))
       else if (diff(2,3) >0 .and. diff(3,2)>0) then
                !
                ! 4NW <- 3NE
                !  |      ^
                !  v      |             23   33
                ! 1SW -> 2SE            22   32
          x1 = cart2cubedsphere( cvlist(ie)%vert(2,2,3),elem(ie)%FaceNum  )
          x2 = cart2cubedsphere( cvlist(ie)%vert(4,2,3),elem(ie)%FaceNum  )
          s=.99
          x3%x = x2%x + (x1%x-x2%x)*s
          x3%y = x2%y + (x1%y-x2%y)*s

          sq1(1) = x3
          sq1(2) = cart2cubedsphere( cvlist(ie)%vert(3,2,3),elem(ie)%FaceNum  ) 
          sq1(3) = cart2cubedsphere( cvlist(ie)%vert(4,2,3),elem(ie)%FaceNum  )
          sq1(4) = cart2cubedsphere( cvlist(ie)%vert(1,2,3),elem(ie)%FaceNum  )

          x2 = cart2cubedsphere( cvlist(ie)%vert(2,3,2),elem(ie)%FaceNum  )          
          s=.99
          x3%x = x2%x + (x1%x-x2%x)*s
          x3%y = x2%y + (x1%y-x2%y)*s

          sq2(1) = x3
          sq2(2) = cart2cubedsphere( cvlist(ie)%vert(1,3,2),elem(ie)%FaceNum  )
          sq2(3) = cart2cubedsphere( cvlist(ie)%vert(2,3,2),elem(ie)%FaceNum  )
          sq2(4) = cart2cubedsphere( cvlist(ie)%vert(3,3,2),elem(ie)%FaceNum  )

          pent(1) = cart2cubedsphere( cvlist(ie)%vert(4,2,2),elem(ie)%FaceNum  )
          pent(2) = cart2cubedsphere( cvlist(ie)%vert(1,2,2),elem(ie)%FaceNum  )
          pent(3) = cart2cubedsphere( cvlist(ie)%vert(2,2,2),elem(ie)%FaceNum  )
          pent(4) = sq2(1)
          pent(5) = sq1(1)

          call pentagon_iteration(sq1,sq2,pent,&
               elem(ie)%spheremp(2,3)*a1(ie)/a2(ie), &
               elem(ie)%spheremp(3,2)*a1(ie)/a2(ie), &
               elem(ie)%spheremp(2,2)*a1(ie)/a2(ie),elem(ie)%FaceNum,a1(ie))

          x2_3d=cubedsphere2cart(sq1(1),elem(ie)%FaceNum  )
          x3_3d=cubedsphere2cart(sq2(1),elem(ie)%FaceNum  )

          cvlist(ie)%vert(2,2,3)=x2_3d

          cvlist(ie)%vert(4,3,2)=x3_3d

          cvlist(ie)%vert(5,2,2)=cvlist(ie)%vert(4,2,2)
          cvlist(ie)%vert(3,2,2)=x3_3d
          cvlist(ie)%vert(4,2,2)=x2_3d


          cvlist(ie)%vert(1,3,3)=x3_3d
          cvlist(ie)%vert(5,3,3)=x2_3d

          cvlist(ie)%nvert(3,3)=sign(5,cvlist(ie)%nvert(3,3))
          cvlist(ie)%nvert(2,2)=sign(5,cvlist(ie)%nvert(2,2))
       else
          print *,ie,'bad type'
          stop
       endif
       ! recompute areas:
       do i=2,3
          do j=2,3
             nvert=abs(cvlist(ie)%nvert(i,j))
             temp3d(1:nvert)=cvlist(ie)%vert(1:nvert,i,j)
             cvlist(ie)%vol(i,j)=surfarea(temp3d,nvert)
             cvlist(ie)%totvol(i,j)=cvlist(ie)%vol(i,j)
          enddo
       enddo
    else
       !print *,'skipping pentagon procedure ie=',ie
       !print *,'maxval diff: ',maxval(abs(diff(2:3,2:3)))/a1(ie) 
       nskip=nskip+1
       nskipie(ie) = nskipie(ie) + 1
    endif
#if 0
    diff = ( cvlist(ie)%vol(:,:) - elem(ie)%spheremp(:,:)*a1(ie)/a2(ie) )
    if (maxval(abs(diff))>1e-8) then
       print *,'after pentagon iteration ie=',ie
       do i=4,1,-1
          write(*,'(a,4e13.4)') 'diff(:,i)=',diff(:,i)
       enddo
       stop
    endif
#endif
     global_shared_buf(ie,1) = nskipie(ie)
     global_shared_buf(ie,2) = npentie(ie)
    enddo
    call wrap_repro_sum(nvars=2, comm=hybrid%par%comm)
    nskip = global_shared_sum(1)
    npent = global_shared_sum(2)
    if (hybrid%masterthread) then
       write(*,'(a,i7,a,i7)') "no. elements where pentagons were added: ",npent,"/",npent+nskip
    endif

#endif

    ! compute output needed for SCRIP:  lat/lon coordinates, and for the 
    ! control volume with only 3 corners, repeat the last point to make a 
    ! degenerate quad.  
    do ie=nets,nete
       do j=1,np
          do i=1,np
             cvlist(ie)%vert_latlon(:,i,j)%lat=0
             cvlist(ie)%vert_latlon(:,i,j)%lon=0
             do k=1,kmax
                rvert = cvlist(ie)%vert(k,i,j)%x**2+cvlist(ie)%vert(k,i,j)%y**2+cvlist(ie)%vert(k,i,j)%z**2
                if(rvert>0.9) then
                   cvlist(ie)%vert_latlon(k,i,j) = change_coordinates(cvlist(ie)%vert(k,i,j))
                else
                   ! coordinates = 0, this corner was not set above because this point
                   ! only has 3 cells (corner point) pick either neighbor to make a degenerate quad
                   k2=k-1
                   if (k2==0) k2=2  ! can only happen for corner point with data in 2,3,4
                   cvlist(ie)%vert_latlon(k,i,j) = change_coordinates(cvlist(ie)%vert(k2,i,j))
                   cvlist(ie)%vert(k,i,j) = cvlist(ie)%vert(k2,i,j)
                endif
             enddo
          enddo
       enddo
    enddo
    ! Release memory
    if(hybrid%masterthread) then
       call freeedgebuffer(edge1)
    end if
    
    initialized=.true.
  end subroutine  InitControlVolumes_gll





subroutine construct_cv_gll(elem,hybrid,nets,nete)
!
! construct global dual grid from local element dual grid cvlist(ie)%cartp_dual(:,:)
! all control volumes will be squares or triangles (at cube corners)
!
! 10/2009: added option to make hexagon control volumes at cube edges and corners
!
    use bndry_mod, only : bndry_exchangev
    use element_mod, only : element_t
    use hybrid_mod, only : hybrid_t
    use edge_mod, only : edgeVpack, edgeVunpack, edgeVunpackVert
    use parallel_mod, only : abortmp	

    implicit none
    type (element_t),     intent(in), target    :: elem(:)
    type (hybrid_t),      intent(in) :: hybrid
    integer :: nets,nete
!   local      
    integer :: i,j,k,ie,kptr,nvert,ie2
    logical :: corner
    real (kind=real_kind)  :: test(np,np,1),vertpack(np,np,3),rvert
    type (cartesian2d_t)             :: vert(4)
    type (cartesian2d_t)    :: cartp_nm1(0:np,0:np)
    
    test(:,:,:) = 0

    do ie=nets,nete
       ! now construct the dual grid

       cartp_nm1 = cvlist(ie)%cartp_dual

       do j=1,np
          do i=1,np
             cvlist(ie)%vert(:,i,j)%x = 0_real_kind
             cvlist(ie)%vert(:,i,j)%y = 0_real_kind
             cvlist(ie)%vert(:,i,j)%z = 0_real_kind
          enddo
       enddo
       
       ! interior
       
       do j=2,np-1       
          do i=2,np-1
             
             ! internal vertex on Cubed sphere
             ! Here is the order:
             ! 
             ! 4NW <- 3NE
             !  |      ^
             !  v      |
             ! 1SW -> 2SE
             vert(1)%x = cartp_nm1(i-1,j-1)%x
             vert(1)%y = cartp_nm1(i-1,j-1)%y
             vert(2)%x = cartp_nm1(i  ,j-1)%x
             vert(2)%y = cartp_nm1(i  ,j-1)%y
             vert(3)%x = cartp_nm1(i  ,j  )%x
             vert(3)%y = cartp_nm1(i  ,j  )%y
             vert(4)%x = cartp_nm1(i-1,j  )%x
             vert(4)%y = cartp_nm1(i-1,j  )%y

             do k=1,4
                cvlist(ie)%vert(k,i,j) = cubedsphere2cart(vert(k),elem(ie)%FaceNum)
             enddo
             cvlist(ie)%nvert(i,j) = 4

          end do
       enddo

       ! Compute everything on the edges and then sum
       do i=2,np-1
          j=1
          !              
          ! 4NW <- 3NE   
          !  |      ^    
          !  v      |    
          ! 1SW -> 2SE   
          !              
          !
          ! only pack top two nodes.  
          ! leave other data zero, filled in by edgeexchange
          cvlist(ie)%vert(4,i,j)%x = cvlist(ie)%vert(1,i,j+1)%x
          cvlist(ie)%vert(4,i,j)%y = cvlist(ie)%vert(1,i,j+1)%y
          cvlist(ie)%vert(4,i,j)%z = cvlist(ie)%vert(1,i,j+1)%z
          cvlist(ie)%vert(3,i,j)%x = cvlist(ie)%vert(2,i,j+1)%x
          cvlist(ie)%vert(3,i,j)%y = cvlist(ie)%vert(2,i,j+1)%y
          cvlist(ie)%vert(3,i,j)%z = cvlist(ie)%vert(2,i,j+1)%z


          j=np

          cvlist(ie)%vert(1,i,j)%x = cvlist(ie)%vert(4,i,j-1)%x
          cvlist(ie)%vert(1,i,j)%y = cvlist(ie)%vert(4,i,j-1)%y
          cvlist(ie)%vert(1,i,j)%z = cvlist(ie)%vert(4,i,j-1)%z
          cvlist(ie)%vert(2,i,j)%x = cvlist(ie)%vert(3,i,j-1)%x
          cvlist(ie)%vert(2,i,j)%y = cvlist(ie)%vert(3,i,j-1)%y
          cvlist(ie)%vert(2,i,j)%z = cvlist(ie)%vert(3,i,j-1)%z

       enddo

       do j=2,np-1
          i=1

          cvlist(ie)%vert(2,i,j)%x = cvlist(ie)%vert(1,i+1,j)%x
          cvlist(ie)%vert(2,i,j)%y = cvlist(ie)%vert(1,i+1,j)%y
          cvlist(ie)%vert(2,i,j)%z = cvlist(ie)%vert(1,i+1,j)%z
          cvlist(ie)%vert(3,i,j)%x = cvlist(ie)%vert(4,i+1,j)%x
          cvlist(ie)%vert(3,i,j)%y = cvlist(ie)%vert(4,i+1,j)%y
          cvlist(ie)%vert(3,i,j)%z = cvlist(ie)%vert(4,i+1,j)%z

          i=np

          cvlist(ie)%vert(4,i,j)%x = cvlist(ie)%vert(3,i-1,j)%x
          cvlist(ie)%vert(4,i,j)%y = cvlist(ie)%vert(3,i-1,j)%y
          cvlist(ie)%vert(4,i,j)%z = cvlist(ie)%vert(3,i-1,j)%z
          cvlist(ie)%vert(1,i,j)%x = cvlist(ie)%vert(2,i-1,j)%x
          cvlist(ie)%vert(1,i,j)%y = cvlist(ie)%vert(2,i-1,j)%y
          cvlist(ie)%vert(1,i,j)%z = cvlist(ie)%vert(2,i-1,j)%z

       enddo

       ! Corners       
       ! SW
       cvlist(ie)%vert(3,1,1)%x = cvlist(ie)%vert(1,2,2)%x
       cvlist(ie)%vert(3,1,1)%y = cvlist(ie)%vert(1,2,2)%y
       cvlist(ie)%vert(3,1,1)%z = cvlist(ie)%vert(1,2,2)%z

       ! SE
       cvlist(ie)%vert(4,np,1)%x = cvlist(ie)%vert(2,np-1,2)%x
       cvlist(ie)%vert(4,np,1)%y = cvlist(ie)%vert(2,np-1,2)%y
       cvlist(ie)%vert(4,np,1)%z = cvlist(ie)%vert(2,np-1,2)%z

       ! NE
       cvlist(ie)%vert(1,np,np)%x = cvlist(ie)%vert(3,np-1,np-1)%x
       cvlist(ie)%vert(1,np,np)%y = cvlist(ie)%vert(3,np-1,np-1)%y
       cvlist(ie)%vert(1,np,np)%z = cvlist(ie)%vert(3,np-1,np-1)%z

       ! NW
       cvlist(ie)%vert(2,1,np)%x = cvlist(ie)%vert(4,2,np-1)%x
       cvlist(ie)%vert(2,1,np)%y = cvlist(ie)%vert(4,2,np-1)%y
       cvlist(ie)%vert(2,1,np)%z = cvlist(ie)%vert(4,2,np-1)%z





       kptr=0
       test(:,:,1) = cvlist(ie)%vol(:,:)
       call edgeVpack(edge1,test(1,1,1),1,kptr,elem(ie)%desc)

       cvlist(ie)%invvol(:,:) = cvlist(ie)%vol(:,:)

    enddo  ! loop over NE


    _DBG_
    call bndry_exchangeV(hybrid,edge1)
    _DBG_

    do ie=nets,nete
       kptr=0
       call edgeVunpack(edge1, cvlist(ie)%invvol(1,1),1, kptr, elem(ie)%desc)
       cvlist(ie)%totvol(:,:)=cvlist(ie)%invvol(:,:)
       cvlist(ie)%invvol(:,:)=1.0_real_kind/cvlist(ie)%invvol(:,:)
    enddo

    ! Create the polygon at the edges of the element


    if(.NOT.(MODULO(np,2)==0)) then
       call abortmp("surfaces_mod: NV odd not implemented")
    endif
    vertpack = 0
    do ie=nets,nete
       ! Special messed up copy
       ! 
       !ASC should be replaced by a edgepack
       ! S+N
       do i=1,np/2
          j=1
          vertpack(i,j,1) = cvlist(ie)%vert(3,i,j)%x
          vertpack(i,j,2) = cvlist(ie)%vert(3,i,j)%y
          vertpack(i,j,3) = cvlist(ie)%vert(3,i,j)%z
          j=np
          vertpack(i,j,1) = cvlist(ie)%vert(2,i,j)%x
          vertpack(i,j,2) = cvlist(ie)%vert(2,i,j)%y
          vertpack(i,j,3) = cvlist(ie)%vert(2,i,j)%z
       enddo

       do i=np/2+1,np
          j=1
          vertpack(i,j,1) = cvlist(ie)%vert(4,i,j)%x
          vertpack(i,j,2) = cvlist(ie)%vert(4,i,j)%y
          vertpack(i,j,3) = cvlist(ie)%vert(4,i,j)%z
          j=np
          vertpack(i,j,1) = cvlist(ie)%vert(1,i,j)%x
          vertpack(i,j,2) = cvlist(ie)%vert(1,i,j)%y
          vertpack(i,j,3) = cvlist(ie)%vert(1,i,j)%z
       enddo

       ! E+W
       do j=2,np/2
          i=1
          vertpack(i,j,1) = cvlist(ie)%vert(3,i,j)%x
          vertpack(i,j,2) = cvlist(ie)%vert(3,i,j)%y
          vertpack(i,j,3) = cvlist(ie)%vert(3,i,j)%z
          i=np
          vertpack(i,j,1) = cvlist(ie)%vert(4,i,j)%x
          vertpack(i,j,2) = cvlist(ie)%vert(4,i,j)%y
          vertpack(i,j,3) = cvlist(ie)%vert(4,i,j)%z
       enddo

       do j=np/2+1,np-1
          i=1
          vertpack(i,j,1) = cvlist(ie)%vert(2,i,j)%x
          vertpack(i,j,2) = cvlist(ie)%vert(2,i,j)%y
          vertpack(i,j,3) = cvlist(ie)%vert(2,i,j)%z
          i=np
          vertpack(i,j,1) = cvlist(ie)%vert(1,i,j)%x
          vertpack(i,j,2) = cvlist(ie)%vert(1,i,j)%y
          vertpack(i,j,3) = cvlist(ie)%vert(1,i,j)%z
       enddo

       do j=2,np-1
          do i=2,np-1
             vertpack(i,j,1) =0_real_kind
             vertpack(i,j,2) =0_real_kind
             vertpack(i,j,3) =0_real_kind
          enddo
       enddo

       kptr=0
       call edgeVpack(edge1,vertpack,3,kptr,elem(ie)%desc)
    enddo

    call bndry_exchangeV(hybrid,edge1)

    do ie=nets,nete
       kptr=0
       call edgeVunpackVert(edge1, cvlist(ie)%vert,elem(ie)%desc)       
       ! Count and orient vert array
       ! nvert is an integer: -4,-3,3,4
       ! Positive: 1->2->3->4 is counter clockwise on the sphere
       ! Negative: clockwise orientation
       do j=1,np
          do i=1,np
             nvert=0
             do k=1,4
                rvert = cvlist(ie)%vert(k,i,j)%x**2+cvlist(ie)%vert(k,i,j)%y**2+cvlist(ie)%vert(k,i,j)%z**2
                if(rvert>0.9_real_kind)nvert=nvert+1
             enddo
             if(.NOT.Orientation(cvlist(ie)%vert(:,i,j),elem(ie)%FaceNum))nvert=-nvert
             cvlist(ie)%nvert(i,j) = nvert
             corner = ( ((i==1) .and. (j==1)) .or. &
                  ((i==1) .and. (j==np)) .or. &
                  ((i==np) .and. (j==1)) .or. &
                  ((i==np) .and. (j==np)) )
             if (abs(nvert)/=4) then
                if (abs(nvert)/=3) then
                   print *,'i,j,nvert=',i,j,nvert
                   stop 'bad value of nvert'
                endif
                if (.not. corner) then
                   print *,'ie,i,j,nvert,corner=',ie,i,j,nvert,corner
                   print *,cvlist(ie)%vert(1,i,j)%x
                   print *,cvlist(ie)%vert(2,i,j)%x
                   print *,cvlist(ie)%vert(3,i,j)%x
                   print *,cvlist(ie)%vert(4,i,j)%x
                   !print *,'dual:'
                   !do ie2=nets,nete
                   !   print *,ie2,maxval(cvlist(ie2)%cartp_dual(:,:)%x)
                   !   print *,ie2,maxval(cvlist(ie2)%cartp_dual(:,:)%y)
                   !enddo
                   stop 'ERROR: corner point should have nvert=3'
                endif
             endif
          enddo
       enddo
    enddo
end subroutine  construct_cv_gll




  function Orientation(v,FaceNum) result(orient)

    type (cartesian3d_t)  :: v(3)
    integer               :: FaceNum
    type (cartesian3D_t)  :: v12,v23
    real (kind=real_kind) :: test,cart(3,3)

    logical               :: orient

    orient=.FALSE.

    if ((FaceNum == 5).OR.(FaceNum == 6)) then

       cart(1,1) = v(1)%x
       cart(2,1) = v(1)%y
       cart(3,1) = v(1)%z

       cart(1,2) = v(2)%x
       cart(2,2) = v(2)%y
       cart(3,2) = v(2)%z

       cart(1,3) = v(3)%x
       cart(2,3) = v(3)%y
       cart(3,3) = v(3)%z

       v12%x = cart(1,2) - cart(1,1)
       v12%y = cart(2,2) - cart(2,1)
       v12%z = cart(3,2) - cart(3,1)

       v23%x = cart(1,3) - cart(1,2)
       v23%y = cart(2,3) - cart(2,2)
       v23%z = cart(3,3) - cart(3,2)

       test = (v12%y*v23%z - v12%z*v23%y)*v12%x &
            - (v12%x*v23%z - v12%z*v23%x)*v12%y &
            + (v12%x*v23%y - v12%y*v23%x)*v12%z

       if (test > 0_real_kind)then
          orient=.TRUE.
       endif

    else
       orient=.TRUE.
    endif

  end function Orientation




  subroutine VerifVolumes(elem, hybrid,nets,nete)
    use reduction_mod, only : red_sum, parallelmin, parallelmax
    use hybrid_mod, only : hybrid_t
    use physical_constants, only : DD_PI
    use element_mod, only : element_t
    type(element_t), intent(in) :: elem(:)
    integer,              intent(in) :: nets,nete
    type (hybrid_t),      intent(in) :: hybrid

    real (kind=real_kind)            :: psum,ptot,Vol_tmp(1),corr,maxelem_variation
    real (kind=real_kind)            :: vol(np,np,nets:nete),r,rmin,rmax,a1,a2,locmin,locmax,emin,emax,dx,dy
    integer                          :: i,j,ie,kptr,face

    real (kind=real_kind), dimension(:,:), pointer :: locvol

    dx = DD_PI/(2.0d0*dble(ne))
    dy = dx

    if(.not. initialized) call abortmp('Attempt to use volumes prior to initializing')
    rmin=2
    rmax=0
    maxelem_variation=0
    do ie=nets,nete
       locvol => GetVolume(ie)
       locmin = minval(locvol(:,:)*elem(ie)%rspheremp(:,:))
       locmax = maxval(locvol(:,:)*elem(ie)%rspheremp(:,:))
       rmin = min(rmin,locmin)
       rmax = max(rmax,locmax)

       if (locmax>1.01) then
!       do i=4,1,-1
          print *,'locmin(:,i)=',ie,locvol(1,1),1/elem(ie)%rspheremp(1,1)
!       enddo
       endif


       if (locmax-locmin > maxelem_variation) then
          maxelem_variation = locmax-locmin
          emin=locmin
          emax=locmax
       endif
    enddo
    rmin = ParallelMin(rmin,hybrid)
    rmax = ParallelMax(rmax,hybrid)
    if(hybrid%masterthread) then
       write(*,'(a,2e14.7)') "Min/max ratio between spherical and GLL area:",rmin,rmax
    endif
    if (maxelem_variation == ParallelMax(maxelem_variation,hybrid) ) then
       write(*,'(a,2e14.7)') "Min/max ratio element with largest variation:",emin,emax
    endif


    rmin=2
    rmax=0
    do ie=nets,nete
       a1 = SurfArea_dxdy(dx,dy,elem(ie)%cartp(1,1))
       a2 = sum(elem(ie)%spheremp(:,:))
       r=a1/a2
       rmin = min(r,rmin)
       rmax = max(r,rmax)
    enddo
    rmin = ParallelMin(rmin,hybrid)
    rmax = ParallelMax(rmax,hybrid)
    if(hybrid%masterthread) then
       write(*,'(a,2f12.9)') "Min/max ratio spherical and GLL element area:",rmin,rmax
    endif

    do ie=nets,nete
       global_shared_buf(ie,1:6) = 0.d0
       face = elem(ie)%FaceNum
       locvol => GetVolumeLocal(ie)
       do j=1,np
          do i=1,np
             global_shared_buf(ie,face) = global_shared_buf(ie,face) + locvol(i,j)
          enddo
       enddo
    enddo
    call wrap_repro_sum(nvars=6, comm=hybrid%par%comm)

    ptot=0_real_kind
    _DBG_
    do face=1,6
       red_sum%buf(1) = global_shared_sum(face)
       psum = red_sum%buf(1)

       ptot = ptot + psum

       if(hybrid%masterthread) then
          write(*,'(a,i2,a,2e23.15)') "cube face:",face," : SURFACE FV =",&
               6_real_kind*psum/(4_real_kind * DD_PI), &
               6_real_kind*psum/(4_real_kind * DD_PI)-1
       endif
    enddo
    _DBG_

    if(hybrid%masterthread) then
       print *,"SURFACE FV (total)= ", ptot/(4_real_kind * DD_PI)
       print *
    endif


100 format (A13,1(E24.15))

  end subroutine VerifVolumes



#if defined (_AIX) || defined (_BGL) || defined (_BGP) 
   ! for these machines, we use (see above):
   ! use ieee_arithmetic, only: isnan => ieee_is_nan 
#else
! INTEL compiler does have a isnan() function 
! PGI compiler does not have isnan() 
logical function isnan(a)
real(8), intent(in) :: a
if (a.ne.a) then
   isnan = .true.
else
   isnan = .false.
end if
return
end function isnan
#endif

end module surfaces_mod


