!MODULE FVM_ANALYTIC_MOD--------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 17.October 2011                                    !
! This module contains all analytical terms for fvm                           !
!-----------------------------------------------------------------------------!
module fvm_analytic_mod
  use shr_kind_mod,   only: r8=>shr_kind_r8
  use control_mod,    only : north, south, east, west, neast, nwest, seast, swest
  use cam_abortutils, only: endrun

  implicit none
  private

  public :: get_high_order_weights_over_areas,  compute_reconstruct_matrix
  public :: compute_halo_vars, init_flux_orient
  public :: I_00, I_10, I_01, I_20, I_02, I_11, gauss_points
  public :: F_00, F_10, F_01, F_20, F_02, F_11
  public :: create_interpolation_points, compute_basic_coordinate_vars

CONTAINS

  subroutine compute_basic_coordinate_vars(elem,&
       nc,irecons,dalpha,dbeta,vtx_cart,center_cart,area_sphere,spherecentroid)
    use coordinate_systems_mod, only: cart2spherical
    use element_mod,            only: element_t
    use coordinate_systems_mod, only: spherical_polar_t

    type (element_t),         intent(in ) :: elem
    integer,                  intent(in)  :: nc,irecons

    real (kind=r8),           intent(out) :: dalpha, dbeta
    real (kind=r8),           intent(out) :: vtx_cart   (4,2,nc,nc)
    real (kind=r8),           intent(out) :: area_sphere(nc,nc)
    real (kind=r8),           intent(out) :: spherecentroid(irecons-1,nc,nc)
    type (spherical_polar_t), intent(out) :: center_cart(nc,nc) ! Spherical coordinates of fvm grid

    integer        :: i,j
    real (kind=r8) :: centerx,centery
    real (kind=r8) :: acartx(nc+1), acarty(nc+1)

    dalpha=abs(elem%corners(1)%x-elem%corners(2)%x)/nc
    dbeta =abs(elem%corners(1)%y-elem%corners(4)%y)/nc

    do i=1,nc+1
      acartx(i) = tan(elem%corners(1)%x+(i-1)*dalpha)
      acarty(i) = tan(elem%corners(1)%y+(i-1)*dbeta)
    end do

    do j=1,nc
      do i=1,nc
        centerx = tan(elem%corners(1)%x+(i-0.5_r8)*dalpha)
        centery = tan(elem%corners(1)%y+(j-0.5_r8)*dbeta)
        center_cart(i,j) = cart2spherical(centerx,centery,elem%FaceNum)
      enddo
    enddo

    vtx_cart = -9D9
    do j=1,nc
      do i=1,nc
        vtx_cart(1,1,i,j) = acartx(i  )
        vtx_cart(1,2,i,j) = acarty(j  )
 
        vtx_cart(2,1,i,j) = acartx(i+1)
        vtx_cart(2,2,i,j) = acarty(j  )
 
        vtx_cart(3,1,i,j) = acartx(i+1)
        vtx_cart(3,2,i,j) = acarty(j+1)

        vtx_cart(4,1,i,j) = acartx(i  )
        vtx_cart(4,2,i,j) = acarty(j+1)
      end do
    end do
    ! compute area and centroid for the interior and halo zone of interior elements
    call moment_onsphere(nc,irecons,area_sphere,vtx_cart,.true.,spherecentroid)
  end subroutine compute_basic_coordinate_vars

  subroutine compute_halo_vars(faceno,cubeboundary,nc,nhc,nhe,&
       jx_min,jx_max,jy_min,jy_max,flux_orient, ifct, rot_matrix)
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest

    integer,          intent(in)  :: faceno,nc,nhc,nhe,cubeboundary

    integer,          intent(out) :: jx_min(3),jx_max(3),jy_min(3),jy_max(3)
    real (kind=r8),   intent(out) :: flux_orient(2, 1-nhc:nc+nhc,1-nhc:nc+nhc)
    integer,          intent(out) :: ifct          (1-nhc:nc+nhc,1-nhc:nc+nhc)
    integer,          intent(out) :: rot_matrix(2,2,1-nhc:nc+nhc,1-nhc:nc+nhc)

    integer :: i,j
    integer :: rot90_matrix(2,2)
    integer :: ishft


    jx_min(2) = 0; jx_max(2) = -1; jy_min(2) = 0; jy_max(2) = -1
    jx_min(3) = 0; jx_max(3) = -1; jy_min(3) = 0; jy_max(3) = -1

    select case (cubeboundary)
    case (0)
      jx_min(1)=1-nhe; jx_max(1)=nc+1+nhe; jy_min(1)=1-nhe; jy_max(1)=nc+1+nhe
    case (west)
      jx_min(1)=1    ; jx_max(1)=nc+1+nhe; jy_min(1)=1-nhe; jy_max(1)=nc+1+nhe
      jx_min(2)=1-nhe; jx_max(2)=1       ; jy_min(2)=1-nhe; jy_max(2)=nc+1+nhe
    case(east)
      jx_min(1)=1-nhe; jx_max(1)=nc+1    ; jy_min(1)=1-nhe; jy_max(1)=nc+1+nhe
      jx_min(2)=nc+1 ; jx_max(2)=nc+1+nhe; jy_min(2)=1-nhe; jy_max(2)=nc+1+nhe
    case(north)
      jx_min(1)=1-nhe; jx_max(1)=nc+1+nhe; jy_min(1)=1-nhe; jy_max(1)=nc+1
      jx_min(2)=1-nhe; jx_max(2)=nc+1+nhe; jy_min(2)=nc+1 ; jy_max(2)=nc+1+nhe
    case(south)
      jx_min(1)=1-nhe; jx_max(1)=nc+1+nhe; jy_min(1)=1    ; jy_max(1)=nc+1+nhe
      jx_min(2)=1-nhe; jx_max(2)=nc+1+nhe; jy_min(2)=1-nhe; jy_max(2)=1
    case(swest)
      jx_min(1)=1    ; jx_max(1)=nc+1+nhe; jy_min(1)=1    ; jy_max(1)=nc+1+nhe
      jx_min(2)=1    ; jx_max(2)=nc+1+nhe; jy_min(2)=1-nhe; jy_max(2)=1
      jx_min(3)=1-nhe; jx_max(3)=1       ; jy_min(3)=1    ; jy_max(3)=nc+1+nhe
    case(seast)
      jx_min(1)=1-nhe; jx_max(1)=nc+1    ; jy_min(1)=1    ; jy_max(1)=nc+1+nhe
      jx_min(2)=1-nhe; jx_max(2)=nc+1    ; jy_min(2)=1-nhe; jy_max(2)=1
      jx_min(3)=nc+1 ; jx_max(3)=nc+1+nhe; jy_min(3)=1    ; jy_max(3)=nc+1+nhe
    case(neast)
      jx_min(1)=1-nhe; jx_max(1)=nc+1    ; jy_min(1)=1-nhe; jy_max(1)=nc+1
      jx_min(2)=1-nhe; jx_max(2)=nc+1    ; jy_min(2)=nc+1 ; jy_max(2)=nc+1+nhe
      jx_min(3)=nc+1 ; jx_max(3)=nc+1+nhe; jy_min(3)=1-nhe; jy_max(3)=nc+1
    case(nwest)
      jx_min(1)=1    ; jx_max(1)=nc+1+nhe; jy_min(1)=1-nhe; jy_max(1)=nc+1
      jx_min(2)=1    ; jx_max(2)=nc+1+nhe; jy_min(2)=nc+1 ; jy_max(2)=nc+1+nhe
      jx_min(3)=1-nhe; jx_max(3)=1       ; jy_min(3)=1-nhe; jy_max(3)=nc+1

    case default
      print *, 'Fatal Error in fvm_line_integrals_mod.F90.'
      call endrun('Selected case for cubeboundary does not exists!')
    end select
    !
    ! init location of flux-sides
    !
    call init_flux_orient(flux_orient,ifct,nc,nhc,cubeboundary,faceno)
    rot_matrix(1,1,:,:) = 1; rot_matrix(1,2,:,:) = 0;
    rot_matrix(2,1,:,:) = 0; rot_matrix(2,2,:,:) = 1;

    if (cubeboundary>0) then
      !
      ! clockwise 90 rotation of vectors
      !
      rot90_matrix(1,1) = 0; rot90_matrix(2,1) = -1;
      rot90_matrix(1,2) = 1; rot90_matrix(2,2) =  0;
      do j=1-nhc,nc+nhc
        do i=1-nhc,nc+nhc
          do ishft=1,4-nint(flux_orient(2,i,j))
            rot_matrix(:,:,i,j) = MATMUL(rot90_matrix,rot_matrix(:,:,i,j))
          end do
        enddo
      enddo
    end if
  end subroutine compute_halo_vars


  ! ----------------------------------------------------------------------------------!
  !SUBROUTINE MOMENT_ONSPHERE-----------------------------------------------CE-for FVM!
  ! AUTHOR: CHRISTOPH ERATH, 20.July 2011                                             !
  ! DESCRIPTION: Compute area and centroids/moments via line integrals                !
  !                                                                                   !
  ! INPUT:  x  ...  x cartesian coordinats of the arrival grid on the cube            !
  !         y  ...  y cartesian coordinats of the arrival grid on the cube            !
  !            ... cell boundaries in x and y directions                              !
  ! INPUT/OUTPUT:                                                                     !
  !         area      ... area of cells on the sphere                                 !
  !         centroid  ... x,y,x^2,y^2,xy                                              !
  !-----------------------------------------------------------------------------------!
  subroutine moment_onsphere(nc,irecons,area,vtx_cart,lanalytic,spherecentroid)
    use dimensions_mod, only: ngpc

    integer, intent(in)  :: nc,irecons
    real (kind=r8), dimension(nc,nc)                              , intent(out) :: area
    real (kind=r8), dimension(irecons-1,nc,nc), intent(out) :: spherecentroid
    real (kind=r8), dimension(4,2,nc,nc)      , intent(in)  :: vtx_cart
    logical, optional, intent(in) :: lanalytic
    integer                                              :: i,j
    !
    ! variables for call to get_high_order_weights_over_areas
    !
    integer, parameter :: num_area=1, num_seg_max=2
    REAL(KIND=r8), dimension(2,num_seg_max,num_area) :: xx, dxx
    integer             , dimension(num_area               ), parameter :: num_seg=2
    REAL(KIND=r8), dimension(irecons,num_area):: weights
    real (kind=r8), dimension(nc+1) :: x, y


    real (kind=r8), dimension(ngpc):: gsweights, gspts
    !
    ! initialize quadrature weights for get_high_order_weights_over_areas
    !
    call gauss_points(ngpc,gsweights,gspts) !set gauss points/weights
    gspts = 0.5_r8*(gspts+1.0_r8) !shift location so in [0:1] instead of [-1:1]

    x(1:nc) = vtx_cart(1,1,1:nc,1   )
    y(1:nc) = vtx_cart(1,2,1   ,1:nc)
    x(nc+1) = vtx_cart(2,1,  nc,1   )
    y(nc+1) = vtx_cart(3,2,1   ,nc  )

    select case (irecons)
    case(1)
      if (present(lanalytic)) then
        do j=1,nc
          do i=1,nc
            area(i,j) = (I_00(x(i+1),y(j+1)) - I_00(x(i),y(j+1)) + &
                 I_00(x(i),y(j)) - I_00(x(i+1),y(j)))
          end do
        end do
      else
        call endrun("non-analytic moments not coded for irecons=1")
      end if

    case(3)
      if (present(lanalytic)) then
        do j=1,nc
          do i=1,nc
            area(i,j) = (I_00(x(i+1),y(j+1)) - I_00(x(i),y(j+1)) + &
                 I_00(x(i),y(j)) - I_00(x(i+1),y(j)))
            ! Compute centroids via line integrals
            spherecentroid(1,i,j) = (I_10(x(i+1),y(j+1)) - I_10(x(i),y(j+1)) + &
                 I_10(x(i),y(j)) - I_10(x(i+1),y(j))) / area(i,j)
            spherecentroid(2,i,j) = (I_01(x(i+1),y(j+1)) - I_01(x(i),y(j+1)) + &
                 I_01(x(i),y(j)) - I_01(x(i+1),y(j))) / area(i,j)
          end do
        end do
      else
        call endrun("non-analytic moments not coded for irecons=3")
      end if


    case(6)
      if (present(lanalytic)) then
        do j=1,nc
          do i=1,nc
            !       area(i,j) = surfareaxy(x(i),x(i+1),y(j),y(j+1))
            area(i,j) = (I_00(x(i+1),y(j+1)) - I_00(x(i),y(j+1)) + &
                 I_00(x(i),y(j)) - I_00(x(i+1),y(j)))
            ! Compute centroids via line integrals
            spherecentroid(1,i,j) = (I_10(x(i+1),y(j+1)) - I_10(x(i),y(j+1)) + &
                 I_10(x(i),y(j)) - I_10(x(i+1),y(j))) / area(i,j)
            spherecentroid(2,i,j) = (I_01(x(i+1),y(j+1)) - I_01(x(i),y(j+1)) + &
                 I_01(x(i),y(j)) - I_01(x(i+1),y(j))) / area(i,j)
            ! TAN(alpha)^2 component
            spherecentroid(3,i,j) = (I_20(x(i+1),y(j+1)) - I_20(x(i),y(j+1)) + &
                 I_20(x(i),y(j)) - I_20(x(i+1),y(j))) / area(i,j)
            ! TAN(beta)^2 component
            spherecentroid(4,i,j) = (I_02(x(i+1),y(j+1)) - I_02(x(i),y(j+1)) + &
                 I_02(x(i),y(j)) - I_02(x(i+1),y(j))) / area(i,j)
            ! TAN(alpha) TAN(beta) component
            spherecentroid(5,i,j) = (I_11(x(i+1),y(j+1)) - I_11(x(i),y(j+1)) + &
                 I_11(x(i),y(j)) - I_11(x(i+1),y(j))) / area(i,j)
          end do
        end do
      else
        do j=1,nc
          do i=1,nc

            xx (1,1,1) = x(i)       ; xx (2,1,1) = y(j+1);
            dxx(1,1,1) = x(i+1)-x(i); dxx(2,1,1) = 0.0_r8 ;

            xx (1,2,1) = x(i+1)     ; xx (2,2,1) = y(j)  ;
            dxx(1,2,1) = x(i)-x(i+1); dxx(2,2,1) = 0.0_r8 ;

            call get_high_order_weights_over_areas(xx,dxx,num_seg,num_seg_max,num_area,weights,ngpc,gsweights,gspts,irecons)

            area(i,j)         = weights(1,1)

            spherecentroid(1:5,i,j) = weights(2:6,1)/area(i,j)
          end do
        end do
      end if
    case default
      call endrun('SUBROUTINE moment_on_sphere: irecons out of range')
    end select
  end subroutine moment_onsphere


  ! ----------------------------------------------------------------------------------!
  !SUBROUTINES I_00, I_01, I_20, I_02, I11----------------------------------CE-for FVM!
  ! AUTHOR: CHRISTOPH ERATH, 17.October 2011                                          !
  ! DESCRIPTION: calculates the exact integrals                                       !
  !                                                                                   !
  ! CALLS: none                                                                       !
  ! INPUT: x    ... x coordinate of the evaluation point (Cartesian on the cube)      !
  !        y    ... y coordinate of the evaluation point (Cartesian on the cube)      !
  ! OUTPUT: I_00, I_01, I_20, I_02, I11                                               !
  !-----------------------------------------------------------------------------------!
  function I_00(x,y)
    implicit none
    real (kind=r8)                 :: I_00
    real (kind=r8), intent(in)     :: x,y

    I_00 = ATAN(x*y/SQRT(1.0_r8+x*x+y*y))
  end function I_00

  function I_10(x,y)
    implicit none
    real (kind=r8)                 :: I_10
    real (kind=r8), intent(in)     :: x,y
    real (kind=r8)                 :: tmp

    !   tmp = ATAN(x)
    !   I_10 = -ASINH(y*COS(tmp))
    tmp = y*COS(ATAN(x))
    I_10 = -log(tmp+sqrt(tmp*tmp+1))
  end function I_10


  function I_01(x,y)
    implicit none
    real (kind=r8)                 :: I_01
    real (kind=r8), intent(in)     :: x,y
    real (kind=r8)                 :: tmp

    !   I_01 = -ASINH(x/SQRT(1+y*y))
    tmp=x/SQRT(1+y*y)
    I_01 = -log(tmp+sqrt(tmp*tmp+1))
  end function I_01

  function I_20(x,y)
    implicit none
    real (kind=r8)                 :: I_20
    real (kind=r8), intent(in)     :: x,y
    real (kind=r8)                 :: tmp,tmp1

    tmp = 1.0_r8+y*y
    tmp1=x/SQRT(tmp)
    I_20 = y*log(tmp1+sqrt(tmp1*tmp1+1))+ACOS(x*y/(SQRT((1.0_r8+x*x)*tmp)))
  end function I_20

  function I_02(x,y)
    implicit none
    real (kind=r8)                  :: I_02
    real (kind=r8), intent(in)      :: x,y
    real (kind=r8)                  :: tmp,tmp1

    !   tmp=1.0_r8+x*x
    !   I_02 = x*ASINH(y/SQRT(tmp))+ACOS(x*y/SQRT(tmp*(1+y*y)))
    tmp=1.0_r8+x*x
    tmp1=y/SQRT(tmp)

    I_02 = x*log(tmp1+sqrt(tmp1*tmp1+1))+ACOS(x*y/SQRT(tmp*(1+y*y)))

  end function I_02

  function I_11(x,y)
    implicit none
    real (kind=r8)                   :: I_11
    real (kind=r8), intent(in)       :: x,y

    I_11 = -SQRT(1+x*x+y*y)
  end function I_11
  !END SUBROUTINES I_00, I_01, I_20, I_02, I11------------------------------CE-for FVM!


  real (kind=r8) function F_00(x_in,y_in)
    implicit none
    real (kind=r8), intent(in) :: x_in,y_in
    real (kind=r8) :: x,y
    !
    x = x_in
    y = y_in
    F_00 =y/((1.0_r8+x*x)*SQRT(1.0_r8+x*x+y*y))
  end function F_00

  real (kind=r8) function F_10(x_in,y_in)
    implicit none
    real (kind=r8), intent(in) :: x_in,y_in
    real (kind=r8) :: x,y

    x = x_in
    y = y_in

    F_10 =x*y/((1.0_r8+x*x)*SQRT(1.0_r8+x*x+y*y))
  end function F_10

  real (kind=r8) function F_01(x_in,y_in)
    implicit none
    real (kind=r8), intent(in) :: x_in,y_in
    real (kind=r8) :: x,y

    x = x_in
    y = y_in

    F_01 =-1.0_r8/(SQRT(1.0_r8+x*x+y*y))
  end function F_01

  real (kind=r8) function F_20(x_in,y_in)
    implicit none
    real (kind=r8), intent(in) :: x_in,y_in
    real (kind=r8) :: x,y

    x = x_in
    y = y_in

    F_20 =x*x*y/((1.0_r8+x*x)*SQRT(1.0_r8+x*x+y*y))
  end function F_20

  real (kind=r8) function F_02(x_in,y_in)
    implicit none
    real (kind=r8), intent(in) :: x_in,y_in
    real (kind=r8) :: x,y,alpha,tmp

    x = x_in
    y = y_in

    alpha = ATAN(x)
!     F_02 =-y/SQRT(1.0_r8+x*x+y*y)+ASINH(y*COS(alpha))
    tmp=y*COS(alpha)
    F_02 =-y/SQRT(1.0_r8+x*x+y*y)+log(tmp+sqrt(tmp*tmp+1))

    !
    ! cos(alpha) = 1/sqrt(1+x*x)
    !
  end function F_02

  real (kind=r8) function F_11(x_in,y_in)
    implicit none
    real (kind=r8), intent(in) :: x_in,y_in
    real (kind=r8) :: x,y

    x = x_in
    y = y_in

    F_11 =-x/(SQRT(1.0_r8+x*x+y*y))
  end function F_11



  !
  ! matrix version of reconstruct_cubic_onface
  !
  subroutine compute_reconstruct_matrix(nc,nhe,nhc,irecons,dalpha,dbeta,spherecentroid,vtx_cart,&
       centroid_stretch,vertex_recons_weights,recons_metrics,recons_metrics_integral)
    implicit none
    integer              , intent(in) :: nc,nhe,irecons,nhc
    real (kind=r8), intent(in) :: dalpha,dbeta
    real (kind=r8), dimension(irecons-1,1-nhc:nc+nhc,1-nhc:nc+nhc), intent(in) :: spherecentroid
    real (kind=r8), dimension(4,2,1-nhc:nc+nhc,1-nhc:nc+nhc)      , intent(in) :: vtx_cart

    real (kind=r8), dimension(7,1-nhe:nc+nhe,1-nhe:nc+nhe)        , intent(out):: centroid_stretch
    real (kind=r8), dimension(1:irecons-1,4,1-nhe:nc+nhe,1-nhe:nc+nhe), intent(out):: vertex_recons_weights
    real (kind=r8), dimension(3,1-nhe:nc+nhe,1-nhe:nc+nhe)        , intent(out):: recons_metrics
    real (kind=r8), dimension(3,1-nhe:nc+nhe,1-nhe:nc+nhe)        , intent(out):: recons_metrics_integral

    !
    integer  :: i, j, count, m, n
    real (kind=r8) :: coef,tmp,cartx,carty
    !
    ! pre-compute variables for reconstruction
    !
    select case (irecons)
    case(3)
      do j= 1-nhe,nc+nhe
        do i=1-nhe,nc+nhe
          count = 1
          do n = j, j+1
            do m = i, i+1
              cartx = vtx_cart(count,1,i,j); carty = vtx_cart(count,2,i,j);

              vertex_recons_weights(1,count,i,j) = cartx - spherecentroid(1,i,j)
              vertex_recons_weights(2,count,i,j) = carty - spherecentroid(2,i,j)

              count=count+1
            end do
          enddo
        end do
      end do
      call endrun("recons_metrics and recons_metrics_integral not initialize")
      !
      ! for reconstruction
      !
      do j= 1-nhe,nc+nhe
        do i=1-nhe,nc+nhe
          !
          !***************
          !*    dfdx     *
          !***************
          !
          coef = 1.0_r8/(12.0_r8 * dalpha)                   !finite difference coefficient
          coef = coef /( 1.0_r8 + spherecentroid(1,i,j)**2) !stretching coefficient

          centroid_stretch(1,i,j) = coef
          !
          !***************
          !*    dfdy     *
          !***************
          !
          coef = 1.0_r8/(12.0_r8 * dbeta)                    !finite difference coefficient
          coef = coef /( 1.0_r8 + spherecentroid(2,i,j)**2) !stretching coefficient

          centroid_stretch(2,i,j) = coef
        end do
      end do
    case(6)
      do j= 1-nhe,nc+nhe
        do i=1-nhe,nc+nhe
          do count=1,4
            cartx = vtx_cart(count,1,i,j); carty = vtx_cart(count,2,i,j);

            vertex_recons_weights(1,count,i,j) = cartx - spherecentroid(1,i,j)
            vertex_recons_weights(2,count,i,j) = carty - spherecentroid(2,i,j)

            vertex_recons_weights(3,count,i,j) = (spherecentroid(1,i,j)**2 - &
                 spherecentroid(3,i,j))   + &
                 (cartx - spherecentroid(1,i,j))**2
            vertex_recons_weights(4,count,i,j) = (spherecentroid(2,i,j)**2 - &
                 spherecentroid(4,i,j)) + &
                 (carty - spherecentroid(2,i,j))**2

            vertex_recons_weights(5,count,i,j) = (cartx - spherecentroid(1,i,j))*     &
                 (carty - spherecentroid(2,i,j))+     &
                 (spherecentroid(1,i,j) *             &
                 spherecentroid(2,i,j) - &
                 spherecentroid(5,i,j))
          end do
        end do
      end do

      do j= 1-nhe,nc+nhe
        do i=1-nhe,nc+nhe
          recons_metrics(1,i,j) = spherecentroid(1,i,j)**2 -spherecentroid(3,i,j)
          recons_metrics(2,i,j) = spherecentroid(2,i,j)**2 -spherecentroid(4,i,j)
          recons_metrics(3,i,j) = spherecentroid(1,i,j)*spherecentroid(2,i,j)-&
               spherecentroid(5,i,j)

          recons_metrics_integral(1,i,j) = &
               2.0_r8*spherecentroid(1,i,j)**2 -spherecentroid(3,i,j)
          recons_metrics_integral(2,i,j) = &
               2.0_r8*spherecentroid(2,i,j)**2 -spherecentroid(4,i,j)
          recons_metrics_integral(3,i,j) = &
               2.0_r8*spherecentroid(1,i,j)*spherecentroid(2,i,j)-&
               spherecentroid(5,i,j)
        end do
      end do



      !
      ! pre-compute variables for reconstruction
      !
      do j= 1-nhe,nc+nhe
        do i=1-nhe,nc+nhe
          !
          !***************
          !*    dfdx     *
          !***************
          !
          coef = 1.0_r8/(12.0_r8 * dalpha)                   !finite difference coefficient
          coef = coef /( 1.0_r8 + spherecentroid(1,i,j)**2) !stretching coefficient

          centroid_stretch(1,i,j) = coef
          !
          !***************
          !*    dfdy     *
          !***************
          !
          coef = 1.0_r8/(12.0_r8 * dbeta)                    !finite difference coefficient
          coef = coef /( 1.0_r8 + spherecentroid(2,i,j)**2) !stretching coefficient

          centroid_stretch(2,i,j) = coef

          !*****************
          !*    d2fdx2     *
          !*****************
          !
          coef = 1.0_r8 / (12.0_r8 * dalpha**2)                  !finite difference coefficient
          !
          ! stretching coefficient part 2
          !      recons(3,i,j) = (a * recons(1,i,j)+ recons(3,i,j))*b
          !
          tmp = 0.5_r8/((1.0_r8 + spherecentroid(1,i,j)**2)**2)

          centroid_stretch(3,i,j) = coef*tmp
          centroid_stretch(6,i,j) = -spherecentroid(1,i,j)/(1.0_r8 + spherecentroid(1,i,j)**2)

          !
          !*****************
          !*    d2fdy2     *
          !*****************
          !
          !
          coef = 1.0_r8 / (12.0_r8 * dbeta**2)                     !finite difference coefficient
          !
          ! stretching coefficient part 2
          !
          !      recons(4,i,j) = (a * recons(1,i,j)+ recons(4,i,j))*b
          !
          tmp =0.5_r8/((1.0_r8 + spherecentroid(2,i,j)**2)**2)

          centroid_stretch(4,i,j) = coef*tmp
          centroid_stretch(7,i,j) = -spherecentroid(2,i,j)/(1.0_r8 + spherecentroid(2,i,j)**2)
          !
          !*****************
          !*    d2fdxdy    *
          !*****************
          !
          !
          coef = 1.0_r8 / (4.0_r8 * dalpha * dbeta)            !finite difference coefficient
          coef = coef  / ((1.0_r8 + spherecentroid(1,i,j)**2) * &
               (1.0_r8 + spherecentroid(2,i,j)**2))    !stretching coefficient

          centroid_stretch(5,i,j) = coef
        enddo
      enddo
    case default
      call endrun('SUBROUTINE compute_reconstruct_matrix: irecons out of range')
    end select
  end subroutine compute_reconstruct_matrix


  subroutine get_high_order_weights_over_areas(x,dx,num_seg,num_seg_max,num_area,weights,ngpc,gsweights, gspts,irecons)
    implicit none
    integer                                                 , intent(in)    :: num_area, num_seg_max, irecons
    REAL(KIND=r8), dimension(2,num_seg_max,num_area ), intent(inout) :: x, dx
    integer                                                 , intent(in)    :: ngpc
    integer             , dimension(num_area               ), intent(in)    :: num_seg
    REAL(KIND=r8), dimension(irecons,num_area), intent(out)   :: weights

    real (kind=r8), dimension(ngpc,num_seg_max               ) :: xq,yq        !quadrature points along line segments
    real (kind=r8), dimension(ngpc,num_seg_max,irecons) :: F            !potentials
    real (kind=r8), dimension(                 irecons) :: weights_area
    real (kind=r8), dimension(ngpc,num_seg_max) :: xq2, yrh, rho, tmp !intermediate variables for optimization
    REAL(KIND=r8) , dimension(ngpc,num_seg_max) :: xq2ir, xq2i, rhoi  !intermediate variables for optimization

    integer :: iseg,iarea,i,j,k

    real (kind=r8), dimension(ngpc) :: gsweights, gspts

    weights(1:irecons,1:num_area) = 0.0_r8 !may not be necessary dbgxxx
    do iarea=1,num_area
      do iseg=1,num_seg(iarea)
        xq(:,iseg) = x(1,iseg,iarea)+dx(1,iseg,iarea)*gspts(:)
        yq(:,iseg) = x(2,iseg,iarea)+dx(2,iseg,iarea)*gspts(:)
      end do
      !
      ! potentials (equation's 23-28 in CSLAM paper; Lauritzen et al., 2010):
      !
      ! (Rory Kelly optimization)
      !
      do j=1,num_seg(iarea)
!DIR$ SIMD
        do i=1,ngpc
          xq2(i,j)   =  xq(i,j)*xq(i,j)
          xq2i(i,j)  =  1.0_r8/(1.0_r8+xq2(i,j))
          xq2ir(i,j) =  SQRT(xq2i(i,j))
          rho(i,j)   =  SQRT(1.0_r8+xq2(i,j)+yq(i,j)*yq(i,j))
          rhoi(i,j)  =  1.0_r8/rho(i,j)
          yrh(i,j)   =  yq(i,j)*rhoi(i,j)
          tmp(i,j)   =  yq(i,j)*xq2ir(i,j)
          F(i,j,1)   =  yrh(i,j)*xq2i(i,j)                 !F_00 !F_00
          F(i,j,2)   =  xq(i,j)*yrh(i,j)*xq2i(i,j)         !F_10 !F_10
          F(i,j,3)   = -1.0_r8*rhoi(i,j)                    !F_01 !F_01
          F(i,j,4)   =  xq2(i,j)*yrh(i,j)*xq2i(i,j)        !F_20 !F_20
          F(i,j,6)   = -xq(i,j)*rhoi(i,j)                  !F_11 !F_11
        enddo
        !
        ! take F(i,j,5) out of loop above since it prevents vectorization
        !
        do i=1,ngpc
          F(i,j,5)   = -yq(i,j)*rhoi(i,j)+log(tmp(i,j)+rho(i,j)*xq2ir(i,j))  !F_02 !F_02
        end do
      enddo
      weights_area = 0.0_r8
      do k=1,irecons
        do iseg=1,num_seg(iarea)
          weights_area(k) = weights_area(k) + sum(gsweights(:)*F(:,iseg,k))*0.5_r8*dx(1,iseg,iarea)
        end do
      end do
      weights(1:irecons,iarea) = weights_area(1:irecons)
    end do
  end subroutine get_high_order_weights_over_areas


  !********************************************************************************
  !
  ! Gauss-Legendre quadrature
  !
  ! Tabulated values
  !
  !********************************************************************************
  subroutine gauss_points(n,weights,points)
    implicit none
    integer, intent(in ) :: n
    real (kind=r8), dimension(:), intent(out) :: weights, points !dimension(n)

    select case (n)
      !    CASE(1)
      !       abscissae(1) = 0.0_r8
      !       weights(1)   = 2.0_r8
    case(2)
      points(1)    = -sqrt(1.0_r8/3.0_r8)
      points(2)    =  sqrt(1.0_r8/3.0_r8)
      weights(1)   =  1.0_r8
      weights(2)   =  1.0_r8
    case(3)
      points(1)    = -0.774596669241483377035853079956_r8
      points(2)    =  0.0_r8
      points(3)    =  0.774596669241483377035853079956_r8
      weights(1)   =  0.555555555555555555555555555556_r8
      weights(2)   =  0.888888888888888888888888888889_r8
      weights(3)   =  0.555555555555555555555555555556_r8
    case(4)
      points(1)    = -0.861136311594052575223946488893_r8
      points(2)    = -0.339981043584856264802665659103_r8
      points(3)    =  0.339981043584856264802665659103_r8
      points(4)    =  0.861136311594052575223946488893_r8
      weights(1)   =  0.347854845137453857373063949222_r8
      weights(2)   =  0.652145154862546142626936050778_r8
      weights(3)   =  0.652145154862546142626936050778_r8
      weights(4)   =  0.347854845137453857373063949222_r8
    case(5)
      points(1)    = -(1.0_r8/3.0_r8)*sqrt(5.0_r8+2.0_r8*sqrt(10.0_r8/7.0_r8))
      points(2)    = -(1.0_r8/3.0_r8)*sqrt(5.0_r8-2.0_r8*sqrt(10.0_r8/7.0_r8))
      points(3)    =  0.0_r8
      points(4)    =  (1.0_r8/3.0_r8)*sqrt(5.0_r8-2.0_r8*sqrt(10.0_r8/7.0_r8))
      points(5)    =  (1.0_r8/3.0_r8)*sqrt(5.0_r8+2.0_r8*sqrt(10.0_r8/7.0_r8))
      weights(1)   = (322.0_r8-13.0_r8*sqrt(70.0_r8))/900.0_r8
      weights(2)   = (322.0_r8+13.0_r8*sqrt(70.0_r8))/900.0_r8
      weights(3)   = 128.0_r8/225.0_r8
      weights(4)   = (322.0_r8+13.0_r8*sqrt(70.0_r8))/900.0_r8
      weights(5)   = (322.0_r8-13.0_r8*sqrt(70.0_r8))/900.0_r8
    case default
      call endrun('SUBROUTINE gauss_points: n out of range in (0<n<5)')
    end select

  end subroutine gauss_points



subroutine init_flux_orient(flux_orient,ifct,nc,nhc,cubeboundary,faceno)
  implicit none
  integer              , intent(in)  :: cubeboundary,faceno,nc,nhc
  real (kind=r8), intent(out) :: flux_orient(2  ,1-nhc:nc+nhc,1-nhc:nc+nhc)
  integer              , intent(out) :: ifct           (1-nhc:nc+nhc,1-nhc:nc+nhc)

  integer :: ib

  flux_orient      = 99.9D9 !for debugging
  !
  ! halo of flux_orient will be filled through boundary exchange
  !
  flux_orient (1,1:nc,1:nc) = dble(faceno)
  flux_orient (2,:,:) = 0.0_r8
  ifct(:,:) = 1
  if (cubeboundary>0) then

     !
     ! cshift (permute) value needed to be applied to vertex number so that they match orientation
     ! of the interior of the panel
     !
     !
     ib = cubeboundary
     if (faceno==2) then
        if (ib==north.or.ib==nwest.or.ib==neast) flux_orient(2,1-nhc:nc+nhc,nc+1 :nc+nhc) = 1
        if (ib==south.or.ib==swest.or.ib==seast) flux_orient(2,1-nhc:nc+nhc,1-nhc:0     ) = 3
     end if
     if (faceno==3) then
        if (ib==north.or.ib==nwest.or.ib==neast) flux_orient (2,1-nhc:nc+nhc,nc+1 :nc+nhc) = 2
        if (ib==south.or.ib==swest.or.ib==seast) flux_orient (2,1-nhc:nc+nhc,1-nhc:0     ) = 2
     end if
     if (faceno==4) then
        if (ib==north.or.ib==nwest.or.ib==neast) flux_orient (2,1-nhc:nc+nhc,nc+1 :nc+nhc) = 3
        if (ib==south.or.ib==swest.or.ib==seast) flux_orient (2,1-nhc:nc+nhc,1-nhc:0     ) = 1
     end if
     if (faceno==5) then
        if (ib==south.or.ib==swest.or.ib==seast) flux_orient (2,1-nhc:nc+nhc,1-nhc:0     ) = 2
        if (ib== west.or.ib==swest.or.ib==nwest) flux_orient (2,1-nhc:0     ,1-nhc:nc+nhc) = 3
        if (ib== east.or.ib==seast.or.ib==neast) flux_orient (2, nc+1:nc+nhc,1-nhc:nc+nhc) = 1
     end if

     if (faceno==6) then
        if (ib==north.or.ib==nwest.or.ib==neast ) flux_orient (2,1-nhc:nc+nhc,nc+1 :nc+nhc) = 2
        if (ib==west .or.ib==swest.or.ib==nwest ) flux_orient (2,1-nhc:0    ,1-nhc:nc+nhc)  = 1
        if (ib==east .or.ib==seast.or.ib==neast ) flux_orient (2,nc+1:nc+nhc,1-nhc:nc+nhc)  = 3
     end if
     !
     ! non-existent cells in physical space
     !
     if (cubeboundary==nwest) then
        flux_orient(2,1-nhc:0     ,nc+1 :nc+nhc) = 0
        ifct       (  1-nhc:0     ,nc+1 :nc+nhc) = 0
     else if (cubeboundary==swest) then
        flux_orient (2,1-nhc:0     ,1-nhc:0     ) = 0
        ifct        (  1-nhc:0     ,1-nhc:0     ) = 0
     else if (cubeboundary==neast) then
        flux_orient (2,nc+1 :nc+nhc,nc+1 :nc+nhc) = 0
        ifct        (  nc+1 :nc+nhc,nc+1 :nc+nhc) = 0
     else if (cubeboundary==seast) then
        flux_orient (2,nc+1 :nc+nhc,1-nhc:0     ) = 0
        ifct        (  nc+1 :nc+nhc,1-nhc:0     ) = 0
     end if
   end if

 end subroutine init_flux_orient

!
!
!

! ----------------------------------------------------------------------------------!
!SUBROUTINE CREATE_INTERPOLATIION_POINTS----------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 17.October 2011                                          !
! DESCRIPTION: for elements, which share a cube edge, we have to do some            !
!        interpolation on different cubic faces, also in the halo region:           !
!        because we also need the reconstruction coefficients in the halo zone,     !
!        which is basically calculated twice, on the original cell of an element    !
!        on face A and on a cell in the halo region of an element of face B         !
!        The crux is, that the interpolation has to be the same to ensure           !
!        conservation of the scheme                                                 !
!        SYMMETRY of the CUBE is used for calucaltion the interpolation_point       !
!                                                                                   !
! CALLS: interpolation_point                                                        !
! INPUT/OUTPUT:                                                                     !
!        elem     ...  element structure from HOMME                                 !
!        fvm   ...  structure                                                       !
!-----------------------------------------------------------------------------------!
  subroutine create_interpolation_points(elem,&
       nc,nhc,nhr,ns,nh,cubeboundary,&
       dalpha,dbeta,ibase,halo_interp_weight)
    use element_mod           , only: element_t
    use coordinate_systems_mod, only: cartesian2D_t
    use control_mod           , only: north, south, east, west, neast, nwest, seast, swest
    use cube_mod              , only: cube_xstart, cube_xend, cube_ystart, cube_yend

    implicit none
    type (element_t), intent(in)     :: elem

    integer              , intent(in)   :: nc,nhc,nhr,ns,nh,cubeboundary
    integer              , intent(out) :: ibase(1-nh:nc+nh,1:nhr,2)
    real (kind=r8), intent(out) :: halo_interp_weight(1:ns,1-nh:nc+nh,1:nhr,2)
    !
    ! pre-compute weight/index matrices
    !
    integer :: imin,imax,jmin,jmax,iinterp
    real (kind=r8), intent(in)   :: dalpha,dbeta

    real    (kind=r8), dimension(1-nhc:nc+nhc) :: gnomxstart, gnomxend, gnomystart, gnomyend
    integer                                         :: i, halo, ida, ide, iref1
    type (cartesian2D_t)                            :: tmpgnom
    real (kind=r8)  :: interp(1-nh:nc+nh,1:nhr,2)
    integer                              ::ibaseref
    integer                             :: ibase_tmp(1-nh:nc+nh,1:nhr,2)

  ibase = 99999 !dbg
  halo_interp_weight(:,:,:,:) = 9.99E9_r8 !dbg

  ! element is not on a corner, but shares a cube edge (call of subroutine)
  if(cubeboundary <= 4) then
    gnomxstart(1-nhc)=elem%corners(1)%x-(nhc-0.5_r8)*dalpha
    gnomystart(1-nhc)=elem%corners(1)%y-(nhc-0.5_r8)*dbeta
    do i=2-nhc,nc+nhc
      gnomxstart(i)=gnomxstart(i-1)+dalpha
      gnomystart(i)=gnomystart(i-1)+dbeta
    end do
    ida=1-nhc  !lower bound
    ide=nc+nhc !upper bound
    select case (cubeboundary)
      !INTERIOR element
      case(0)
        ! nothing to do!
      !CASE WEST
      case(west)
        do halo=1,nhr
!          iref1=ida
          tmpgnom%x=cube_xstart-(halo-0.5_r8)*dalpha
          do i=halo-nh,nc+nh-(halo-1) !see fvm_reconstruction to understand these boundaries
            iref1=ida
            tmpgnom%y=gnomystart(i)
            call interpolation_point(nc,ns,tmpgnom,gnomystart,1,4,1,interp(i,halo,1),&
                                     ida,ide,iref1,ibase_tmp(i,halo,1))
          end do
        end do

      !CASE EAST
      case(east)
        ! east zone
        do halo=1,nhr
          iref1=ida
          tmpgnom%x=cube_xend+(halo-0.5_r8)*dalpha
          do i=halo-nh,nc+nh-(halo-1)
            tmpgnom%y=gnomystart(i)
            call interpolation_point(nc,ns,tmpgnom,gnomystart,1,2,1,interp(i,halo,1),&
                                     ida,ide,iref1,ibase_tmp(i,halo,1))
          end do
        end do

      !CASE NORTH
      case(north)
        ! north zone
        do halo=1,nhr
          tmpgnom%y=cube_yend+(halo-0.5_r8)*dbeta
          iref1=ida
          do i=halo-nh,nc+nh-(halo-1)
            tmpgnom%x=gnomxstart(i)
            !
            ! dbg - change to interp(i,halo,1) instead of interp(i,halo,2)
            !       so that I can get rid of iinterp = 1 in fvm_reconstruction_mod
            !
            call interpolation_point(nc,ns,tmpgnom,gnomxstart,1,6,0,interp(i,halo,2),&
                                     ida,ide,iref1,ibase_tmp(i,halo,2))
          end do
        end do
      !CASE SOUTH
      case(south)
       !south zone
       do halo=1,nhr
          iref1=ida
          tmpgnom%y=cube_ystart-(halo-0.5_r8)*dbeta
          do i=halo-nh,nc+nh-(halo-1)
            tmpgnom%x=gnomxstart(i)
            call interpolation_point(nc,ns,tmpgnom,gnomxstart,1,5,0,interp(i,halo,2),&
                                     ida,ide,iref1,ibase_tmp(i,halo,2))
          end do
        end do

        !
        !THIS CASE SHOULD NOT HAPPEN!
     case default
           print *,'Fatal Error in first select statement:'
           call endrun('fvm_reconstruction_mod.F90 subroutine fillhalo_cubic!' )
    end select
  !CORNER TREATMENT
  else
    gnomxstart(1-nhc)=cube_xstart-(nhc-0.5_r8)*dalpha
    gnomxend(nc+nhc)=cube_xend+(nhc-0.5_r8)*dalpha
    gnomystart(1-nhc)=cube_ystart-(nhc-0.5_r8)*dbeta
    gnomyend(nc+nhc)=cube_yend+(nhc-0.5_r8)*dbeta
    do i=2-nhc,nc+nhc
      gnomxstart(i)=gnomxstart(i-1)+dalpha
      gnomxend(nc+1-i)=gnomxend(nc+2-i)-dalpha
      gnomystart(i)=gnomystart(i-1)+dbeta
      gnomyend(nc+1-i)=gnomyend(nc+2-i)-dbeta
    end do

    select case (cubeboundary)
      !CASE SOUTH WEST
      case(swest)
        ! west zone
        do halo=1,nhr
          tmpgnom%x=cube_xstart-(halo-0.5_r8)*dalpha
          ida=1
          ide=nc+nc
          iref1=ida
          do i=0,nc+nh-(halo-1)
            tmpgnom%y=gnomystart(i)
            call interpolation_point(nc,ns,tmpgnom,gnomystart,1,4,1,interp(i,halo,1),&
                                     ida,ide,iref1,ibase_tmp(i,halo,1))
          end do
       end do
      !CASE SOUTH EAST
      case(seast)
        ! east zone
        do halo=1,nhr
          tmpgnom%x=cube_xend+(halo-0.5_r8)*dalpha
          ida=1
          ide=nc+nc
          iref1=ida
          do i=0,nc+nh-(halo-1)
            tmpgnom%y=gnomystart(i)
            call interpolation_point(nc,ns,tmpgnom,gnomystart,1,2,1, interp(i,halo,1),&
                                     ida,ide,iref1,ibase_tmp(i,halo,1))
          end do
       end do
      !CASE NORTH EAST
      case(neast)
        ! east zone
        do halo=1,nhr
          tmpgnom%x=cube_xend+(halo-0.5_r8)*dalpha
          ida=1-nc
          ide=nc
          iref1=ida
          do i=halo-nh,nc+1
            tmpgnom%y=gnomyend(i)
            call interpolation_point(nc,ns,tmpgnom,gnomyend,1,2,1, interp(i,halo,1),&
                                     ida,ide,iref1,ibase_tmp(i,halo,1))
          end do
       end do
      !CASE NORTH WEST
      case(nwest)
        ! west zone
        do halo=1,2
          tmpgnom%x=cube_xstart-(halo-0.5_r8)*dalpha
          ida=1-nc
          ide=nc
          iref1=ida
          do i=halo-nh,nc+1
            tmpgnom%y=gnomyend(i)
            call interpolation_point(nc,ns,tmpgnom,gnomyend,1,4,1, interp(i,halo,1),&
                                     ida,ide,iref1,ibase_tmp(i,halo,1))
          end do
       end do
        !THIS CASE SHOULD NOT HAPPEN!
          case default
            print *,'Fatal Error in second select statement:'
            call endrun('fvm_reconstruction_mod.F90 subroutine create_interpolationpoint!')
         end select
  endif

  !**************************
  !
  ! compute haloe weights and indices
  !
  if (cubeboundary>0) then
     if (cubeboundary<5) then
        !
        ! element is located at a panel side but is not a corner element
        ! (west,east,south,north) = (1,2,3,4)
        !
        if (cubeboundary==west .or.cubeboundary==east ) then
           iinterp = 1
        end if
        if (cubeboundary==north.or.cubeboundary==south) iinterp = 2
        do halo=1,nhr
           do i=halo-nh,nc+nh-(halo-1)
              ibaseref=ibase_tmp(i,halo,iinterp)
              ibase(i,halo,1) = ibaseref
              call get_equispace_weights(dbeta, interp(i,halo,iinterp),&
                   halo_interp_weight(:,i,halo,1),ns)
           end do
        end do
     else
        !
        ! element is located at a cube corner
        ! (swest,seast,nwest,neast)=(5,6,7,8)
        !
        do halo=1,nhr
           if (cubeboundary==swest .or.cubeboundary==seast) then
              imin = 0      ; imax = nc+nh-(halo-1);
              jmin = halo-nh; jmax = nc+1;
           else
              jmin = 0      ; jmax = nc+nh-(halo-1);
              imin = halo-nh; imax = nc+1;
           end if
           do i=imin,imax
              ibaseref=ibase_tmp(i,halo,1)
              ibase(i,halo,1) = ibaseref
              call get_equispace_weights(dbeta, interp(i,halo,1),halo_interp_weight(:,i,halo,1),ns)
           end do
           !
           ! reverse weights/indices for fotherpanel (see details on reconstruct_matrix)
           !
           halo_interp_weight(1:ns,jmin:jmax,halo,2) = halo_interp_weight(ns:1:-1,imax:imin:-1,halo,1)
           ibase       (jmin:jmax,halo     ,2) = nc+1-(ns-1)-ibase(imax:imin:-1,halo        ,1)
        end do
     end if

  end if


end subroutine create_interpolation_points




!END SUBROUTINE CREATE_INTERPOLATION_POINTS-------------------------------CE-for FVM!



! ----------------------------------------------------------------------------------!
!SUBROUTINE INTERPOLATION_POINT-------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 14.November 2011                                         !
! DESCRIPTION: calculates the interpolation point on from face 1 in face 2 in       !
!        alpha/beta coordinates, only 1D                                            !
!                                                                                   !
! CALLS: cubedsphere2cart, cart2cubedsphere                                         !
! INPUT: gnom... 1D coordinates                                                     !
!        gnom1d... 1d coordinates                                                   !
!        face1... orginal face                                                      !
!        face2... target face (where the interpolation has to be done)              !
!        xy ... 0 for alpha coordinate, any other for beta                          !
!        except.which type, interior, left edge (-1), right edge (1)                !
!        point... interpolation point                                               !
!        ida  ... begin of interpval                                                !
!        ide  ... end of interpval                                                  !


! INPUT/OUTPUT/RETURN:                                                              !
!        iref ... where we start the search, is also an OUTPUT, so we know for the  !
!                 next point where to start                                         !
!-----------------------------------------------------------------------------------!
                                      !
! DESCRIPTION: searchs where the interpolation point has to be (iref), two values   !
!        of interpval on the left and on the right, except if we are out of range   !
!        which is indicated through ia and ie, respectively                         !
!        It is a 1D interpolation, use alpha/beta coordinates!!!                    !
!                                                                                   !
! CALLS: cubic_equispace_interp                                                     !
! INPUT: iref ... where we start the search, is also an OUTPUT, so we know for the  !
!                 next point where to start                                         !
!        ibaseref ... startindex of the four tracer value for the reconstruction    !
!        point    ... provides the difference of the interpolation point to use it  !
!                     directly in CUBIC_EQUISPACE_INTERP                            !
!-----------------------------------------------------------------------------------!
function get_gno_point(gnom,face1,face2,xy) result(point)
  use coordinate_systems_mod, only : cubedsphere2cart, cart2cubedsphere, &
                                     cartesian2D_t,cartesian3D_t
  implicit none
  type (cartesian2D_t), intent(in)                     :: gnom
  integer, intent(in)                                  :: face1, face2, xy
  real (kind=r8)                                :: point

  type(cartesian3D_t)                 :: tmpcart3d
  type (cartesian2D_t)                :: tmpgnom

  tmpcart3d=cubedsphere2cart(gnom,face1)
  tmpgnom=cart2cubedsphere(tmpcart3d,face2)
  if(xy==0) then
    point=tmpgnom%x
  else
    point=tmpgnom%y
  end if
end function get_gno_point

subroutine interpolation_point(nc,ns,gnom,gnom1d,face1,face2,xy,point,ida,ide,iref,ibaseref)
  use coordinate_systems_mod, only : cartesian2D_t
  implicit none
  integer                                , intent(in) :: nc,ns
  type (cartesian2D_t), intent(in)                    :: gnom
  real (kind=r8), dimension(1-nc:), intent(in) :: gnom1d  !dimension(1-nhc:nc+nhc)
  integer, intent(in)                                 :: face1, face2, xy
  integer,intent(in)                                  :: ida, ide
  integer,intent(inout)                               :: iref,ibaseref
  real (kind=r8), intent(inout)                :: point

!  type(cartesian3D_t)                 :: tmpcart3d
!  type (cartesian2D_t)                :: tmpgnom

  point = get_gno_point(gnom,face1,face2,xy)

!  tmpcart3d=cubedsphere2cart(gnom,face1)
!  tmpgnom=cart2cubedsphere(tmpcart3d,face2)
!  if(xy==0) then
!    point=tmpgnom%x
!  else
!    point=tmpgnom%y
!  end if
  !
  ! in which cell is interpolation point located? gno(iref) is location of point to the right that is closest
  !
  ! |----------|---------|------x---|----------|------|------
  !                 gno(iref-1)  gno(iref)
  !
  iref=ida
  do while (point>gnom1d(iref))
    iref = iref + 1
    if (iref>ide+1) then
      call endrun("error in search - ABORT; probably invalid ns-nc combination")
    end if
    if (iref>ide) then
!       write(*,*) "extrapolation in interpolation_point",iref,ide
       iref=ide
       exit
    endif
  end do
  !
  ! this routine works for ns=1 and ns even
  !
  if (MOD(ns,2)==1) then
     iref = max(iref,ida+1)!make sure gnom1d does not go out of bounds for extrapolation
     if (gnom1d(iref)-point>point-gnom1d(iref-1)) iref=iref-1
     iref=iref-((ns-1)/2)
     ibaseref = min(max(iref,ida),ide-(ns-1))!extrapolation
     point=point-gnom1d(ibaseref)
  else if (MOD(ns, 2)==0) then
     !
     ! this code is only coded for ns even
     !
     ! ibaseref is the left most index used for 1D interpolation
     ! (hence iref = iref-ns/2 except near corners)
     !
     iref = iref-ns/2
     ibaseref = min(max(iref,ida),ide-(ns-1))
     point=point-gnom1d(ibaseref)
  end if
end subroutine interpolation_point
!END SUBROUTINE INTERPOLATION_POINT---------------------------------------CE-for FVM!
! ---------------------------------------------------------------------!
!                                                                      !
! Precompute weights for Lagrange interpolation                        !
! for equi-distant source grid values                                  !
!                                                                      !
!----------------------------------------------------------------------!

subroutine get_equispace_weights(dx, x, w,ns)
  !
  ! Coordinate system for Lagrange interpolation:
  !
  ! |------|------|------|------|
  ! 0     dx    2*dx   3*dx   ns*dx
  !
  implicit none
  real (kind=r8),intent(in)                  :: dx  ! spacing of points, alpha/beta
  real (kind=r8),intent(in)                  :: x   ! X coordinate where interpolation is to be applied
  real (kind=r8),dimension(:),intent(out)    :: w   ! dimension(ns)
  integer              ,intent(in)                  :: ns
  !
  integer :: j,k
  !
  ! use Lagrange interpolation formulae, e.g.,:
  !
  !                http://mathworld.wolfram.com/LagrangeInterpolatingPolynomial.html
  !
  w = 1.0_r8
  if (ns.ne.1) then
     do j=1,ns
        do k=1,ns
           if (k.ne.j) then
              w(j)=w(j)*(x-dble(k-1)*dx)/(dble(j-1)*dx-dble(k-1)*dx)
           end if
        end do
     end do
  end if
end subroutine get_equispace_weights

end module fvm_analytic_mod
