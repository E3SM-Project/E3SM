#undef _GAUSS_TABLE
module quadrature_mod
  use shr_kind_mod,   only: r8=>shr_kind_r8

  implicit none
  private

  type, public :: quadrature_t
    real (kind=r8), dimension(:), pointer :: points
    real (kind=r8), dimension(:), pointer :: weights
  end type quadrature_t

  public  :: gausslobatto
  public  :: test_gausslobatto
  public  :: gauss
  public  :: test_gauss
  public  :: legendre
  public  :: jacobi
  public  :: quad_norm

  public  :: trapezoid
  private :: trapN
  public  :: simpsons
  public  :: gaussian_int

  private :: gausslobatto_pts
  private :: gausslobatto_wts
  private :: gauss_pts
  private :: gauss_wts
  private :: jacobi_polynomials
  private :: jacobi_derivatives


contains

  ! ==============================================================
  ! gauss:
  !
  ! Find the Gauss collocation points and the corresponding weights.
  !
  ! ==============================================================

  function gauss(npts) result(gs)
    integer, intent(in) :: npts
    type (quadrature_t) :: gs

    allocate(gs%points(npts))
    allocate(gs%weights(npts))

    gs%points=gauss_pts(npts)
    gs%weights=gauss_wts(npts,gs%points)

  end function gauss

#if defined(_GAUSS_TABLE)
  function gauss_pts(npts) result(pts)

    integer, intent(in) :: npts
    real (kind=r8) :: pts(npts)

    pts(1) = -0.93246951420315202781_r8
    pts(2) = -0.66120938646626451366_r8
    pts(3) = -0.23861918608319690863_r8
    pts(4) = -pts(3)
    pts(5) = -pts(2)
    pts(6) = -pts(1)

  end function gauss_pts


  function gauss_wts(npts,pts) result(wts)

    integer, intent(in) :: npts
    real (kind=r8) :: pts(npts)
    real (kind=r8) :: wts(npts)

    wts(1)  =  0.17132449237917034504_r8
    wts(2)  =  0.36076157304813860756_r8
    wts(3)  =  0.46791393457269104738_r8
    wts(4)  =  wts(3)
    wts(5)  =  wts(2)
    wts(6)  =  wts(1)

  end function gauss_wts
#else

  ! ==============================================================
  ! gauss_pts:
  !
  ! Compute the Gauss Collocation points
  ! for Jacobi Polynomials
  !
  ! ==============================================================

  function gauss_pts(np1) result(pts)
    use physconst, only: pi

    integer, intent(in)     :: np1        ! Number of velocity grid points
    real (kind=r8) :: pts(np1)

    ! Local variables

    real (kind=r8) :: alpha,beta
    real (kind=r8) :: xjac(0:np1-1)
    real (kind=r8) :: jac(0:np1)
    real (kind=r8) :: djac(0:np1)

    integer  prec                    ! number of mantissa bits
    real (kind=r8) eps      ! machine epsilon
    real (kind=r8), parameter :: convthresh = 10  ! convergence threshold relative\

    ! to machine epsilon
    integer, parameter :: kstop = 30 ! max iterations for polynomial deflation

    real (kind=r8) :: poly
    real (kind=r8) :: pder
    real (kind=r8) :: recsum,thresh
    real (kind=r8) :: dth

    real (kind=r8) :: x
    real (kind=r8) :: delx
    real (kind=r8) :: c0,c1,c2,c10

    integer i,j,k
    integer n, nh

    n  = np1 - 1
    c0 = 0.0_r8
    c1 = 1.0_r8
    c2 = 2.0_r8
    c10 = 10.0_r8
    alpha = c0
    beta  = c0

    ! =========================================================
    ! compute machine precision and set the convergence
    ! threshold thresh to 10 times that level
    ! =========================================================

    prec   = precision(c10)
    eps    = c10**(-prec)
    thresh = convthresh*eps

    ! ============================================================
    ! Compute first half of the roots by "polynomial deflation".
    ! ============================================================

    dth = PI/(2*n+2)

    nh  = (n+1)/2

    do j=0,nh-1
       x=COS((c2*j+1)*dth)          ! first guess at root
       k=0
       delx=c1
       do while(k<kstop .and. ABS(delx) > thresh)
          call jacobi(n+1,x,alpha,beta,jac(0:n+1),djac(0:n+1))
          poly =  jac(n+1)
          pder = djac(n+1)
          recsum=c0
          do i=0,j-1
             recsum = recsum + c1/(x-xjac(i))
          end do
          delx = -poly/(pder-recsum*poly)
          x = x + delx
          k = k + 1
       end do

       xjac(j)=x

    end do

    ! ================================================
    ! compute the second half of the roots by symmetry
    ! ================================================

    do j=0,nh
       xjac(n-j) = -xjac(j)
    end do

    if (MODULO(n,2)==0) xjac(nh)=c0

    ! ====================================================
    ! Reverse the sign of everything so that indexing
    ! increases with position
    ! ====================================================

    do j=0,n
       pts(j+1) = -xjac(j)
    end do

  end function gauss_pts

  ! ================================================
  ! gauss_wts:
  !
  ! Gauss Legendre Weights
  ! ================================================

  function gauss_wts(np1, gpts) result(wts)

    integer, intent(in)                 :: np1
    real (kind=r8), intent(in) :: gpts(np1)  ! Gauss-Legendre points
    real (kind=r8)             :: wts(np1)   ! Gauss-Legendre weights

    ! Local variables

    real (kind=r8) :: c0,c1,c2
    real (kind=r8) :: alpha
    real (kind=r8) :: beta
    real (kind=r8) :: djac(np1)
    integer i,n

    c0    = 0.0_r8
    c1    = 1.0_r8
    c2    = 2.0_r8

    alpha = c0
    beta  = c0
    n     = np1-1

    djac=jacobi_derivatives(np1,alpha,beta,np1,gpts)     

    do i=1,np1
       wts(i)=c2/((c1-gpts(i)**2)*djac(i)*djac(i))
    end do

  end function gauss_wts

#endif

  ! ==============================================================
  ! test_gauss:
  !
  ! Unit Tester for Gaussian Points, Weights
  ! ==============================================================

  subroutine test_gauss(npts)

    integer, intent(in) :: npts
    type (quadrature_t) :: gs

    integer i
    real (kind=r8) :: gssum
    gs=gauss(npts)

    print *
    print *,"============================================"
    print *,"        Testing Gaussian Quadrature..."
    print *
    print *,"          points              weights"
    print *,"============================================"
    do i=1,npts
       print *,i,gs%points(i),gs%weights(i)
    end do
    print *,"============================================"
    gssum=SUM(gs%weights(:))
    print *,"sum of Gaussian weights=",gssum
    print *,"============================================"

    deallocate(gs%points)
    deallocate(gs%weights)

  end subroutine test_gauss

  ! ==============================================================
  ! gausslobatto:
  !
  ! Find the Gauss-Lobatto Legendre collocation points xgl(i) and the
  ! corresponding weights.
  !
  ! ==============================================================

  function gausslobatto(npts) result(gll)

    integer, intent(in) :: npts
    type (quadrature_t) :: gll

    allocate(gll%points(npts))
    allocate(gll%weights(npts))

    gll%points=gausslobatto_pts(npts)
    gll%weights=gausslobatto_wts(npts,gll%points)

  end function gausslobatto

  ! ==============================================================
  ! gausslobatto_pts:
  !
  ! Compute the Gauss-Lobatto Collocation points 
  ! for Jacobi Polynomials
  !
  ! ==============================================================

  function gausslobatto_pts(np1) result(pts)
    use physconst, only: pi

    integer, intent(in)     :: np1        ! Number of velocity grid points
    real (kind=r8) :: pts(np1)

    ! Local variables

    real (kind=r8) :: alpha,beta
    real (kind=r8) :: xjac(0:np1-1)
    real (kind=r8) :: jac(0:np1)
    real (kind=r8) :: jacm1(0:np1)
    real (kind=r8) :: djac(0:np1)

    integer  prec                    ! number of mantissa bits 
    real (kind=r8) eps      ! machine epsilon
    real (kind=r8), parameter :: convthresh = 10  ! convergence threshold relative 
    ! to machine epsilon 
    integer, parameter :: kstop = 30 ! max iterations for polynomial deflation

    real (kind=r8) :: a,b,det
    real (kind=r8) :: poly
    real (kind=r8) :: pder
    real (kind=r8) :: recsum,thresh
    real (kind=r8) :: dth,cd,sd,cs,ss,cstmp

    real (kind=r8) :: x
    real (kind=r8) :: delx
    real (kind=r8) :: c0,c1,c2,c10

    integer i,j,k
    integer n, nh

    n  = np1 - 1
    c0 = 0.0_r8
    c1 = 1.0_r8
    c2 = 2.0_r8
    c10 = 10.0_r8

    alpha = c0
    beta  = c0

    ! =========================================================
    ! compute machine precision and set the convergence
    ! threshold thresh to 10 times that level 
    ! =========================================================

    prec   = PRECISION(c10)
    eps    = c10**(-prec)
    thresh = convthresh*eps

    ! =====================================================
    ! initialize the end points
    ! =====================================================

    xjac(0) =  c1
    xjac(n) = -c1

    ! ============================================================
    ! Compute first half of the roots by "polynomial deflation".
    ! ============================================================

    ! ============================================================
    ! compute the parameters in the polynomial whose 
    ! roots are desired...
    ! ============================================================

    call jacobi(n+1, c1,alpha,beta,jac(0:n+1),djac(0:n+1))
    call jacobi(n+1,-c1,alpha,beta,jacm1(0:n+1),djac(0:n+1))

    det =   jac(n  )*jacm1(n-1)-jacm1(n  )*jac(n-1)
    a   = -(jac(n+1)*jacm1(n-1)-jacm1(n+1)*jac(n-1))/det
    b   = -(jac(n  )*jacm1(n+1)-jacm1(n  )*jac(n+1))/det

    dth = PI/(2*n+1)
    cd  = COS(c2*dth)
    sd  = SIN(c2*dth)
    cs  = COS(dth)
    ss  = SIN(dth)

    nh  = (n+1)/2

    do j=1,nh-1
       x=cs          ! first guess at root 
       k=0
       delx=c1
       do while(k<kstop .and. ABS(delx) > thresh)
          call jacobi(n+1,x,alpha,beta,jac(0:n+1),djac(0:n+1))
          poly =  jac(n+1)+a* jac(n)+b* jac(n-1)
          pder = djac(n+1)+a*djac(n)+b*djac(n-1)
          recsum=c0
          do i=0,j-1
             recsum = recsum + c1/(x-xjac(i))
          end do
          delx = -poly/(pder-recsum*poly)
          x = x + delx
          k = k + 1
       end do

       xjac(j)=x

       ! =====================================================
       !  compute the guesses for the roots
       !  for the next points, i.e :
       !
       !  ss = sn(theta) => sin(theta+2*dth)
       !  cs = cs(theta) => cs(theta+2*dth)
       ! =====================================================

       cstmp=cs*cd-ss*sd      
       ss=cs*sd+ss*cd    
       cs=cstmp          
    end do

    ! ================================================
    ! compute the second half of the roots by symmetry
    ! ================================================

    do j=1,nh 
       xjac(n-j) = -xjac(j)
    end do

    if (MODULO(n,2)==0) xjac(nh)=c0

    ! ====================================================
    ! Reverse the sign of everything so that indexing
    ! increases with position          
    ! ====================================================

    do j=0,n
       pts(j+1) = -xjac(j)
    end do

  end function gausslobatto_pts

  ! ================================================
  ! Gauss Lobatto Legendre Weights   
  ! ================================================ 

  function gausslobatto_wts(np1, glpts) result(wts)

    integer, intent(in)                 :: np1
    real (kind=r8), intent(in) :: glpts(np1)
    real (kind=r8)             :: wts(np1)

    ! Local variables

    real (kind=r8) :: c0,c2
    real (kind=r8) :: alpha
    real (kind=r8) :: beta
    real (kind=r8) :: jac(np1)
    integer i,n

    c0    = 0.0_r8
    c2    = 2.0_r8
    alpha = c0
    beta  = c0
    n     = np1-1

    jac=jacobi_polynomials(n,alpha,beta,np1,glpts)

    do i=1,np1
       wts(i)=c2/(n*(n+1)*jac(i)*jac(i))
    end do

  end function gausslobatto_wts

  ! ==============================================================
  ! test_gausslobatto:
  !
  ! Unit Tester for Gaussian Lobatto Quadrature...
  ! ==============================================================

  subroutine test_gausslobatto(npts)
    integer, intent(in) :: npts
    type (quadrature_t) :: gll

    integer i
    real (kind=r8) :: gllsum
    gll=gausslobatto(npts)

    print *
    print *,"============================================"
    print *,"      Testing Gauss-Lobatto Quadrature..."
    print *
    print *,"          points              weights"
    print *,"============================================"
    do i=1,npts
       print *,i,gll%points(i),gll%weights(i)
    end do
    print *,"============================================"
    gllsum=SUM(gll%weights(:))
    print *,"sum of Gauss-Lobatto weights=",gllsum
    print *,"============================================"

    deallocate(gll%points)
    deallocate(gll%weights)

  end subroutine test_gausslobatto

  ! ================================================
  !
  ! subroutine jacobi:
  !
  !  Computes the Jacobi Polynomials (jac) and their
  !  first derivatives up to and including degree n 
  !  at point x on the interval (-1,1).
  !
  !    See for example the recurrence relations 
  !    in equation 2.5.4 (page 70) in 
  !
  !    "Spectral Methods in Fluid Dynamics",
  !    by C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A.Zang
  !    Springer-Verlag, 1988.
  ! ================================================

  subroutine jacobi(n, x, alpha, beta, jac, djac)

    integer, intent(in)                 :: n
    real (kind=r8), intent(in) :: x
    real (kind=r8), intent(in) :: alpha
    real (kind=r8), intent(in) :: beta
    real (kind=r8)             :: jac(0:n)
    real (kind=r8)             :: djac(0:n)

    ! Local variables

    real (kind=r8) :: a1k
    real (kind=r8) :: a2k
    real (kind=r8) :: a3k
    real (kind=r8) :: da2kdx

    real (kind=r8) :: c2,c1,c0

    integer ::  k

    c0 = 0.0_r8
    c1 = 1.0_r8
    c2 = 2.0_r8

    jac(0)=c1
    jac(1)=(c1 + alpha)*x

    djac(0)=c0
    djac(1)=(c1 + alpha)

    do k=1,n-1
       a1k      = c2*( k + c1 )*( k + alpha + beta + c1 )*( c2*k + alpha + beta )
       da2kdx   = ( c2*( k + c1 ) + alpha + beta )*( c2*k + alpha + beta + c1 )*( c2*k + alpha + beta )
       a2k      = ( c2*k + alpha + beta + c1 )*( alpha*alpha - beta*beta ) + x*da2kdx
       a3k      = c2*(k + alpha)*( k + beta )*( c2*k + alpha + beta + c2 )
       jac(k+1) = ( a2k*jac(k)-a3k*jac(k-1) )/a1k
       djac(k+1)= ( a2k*djac(k) + da2kdx*jac(k) - a3k*djac(k-1) )/a1k          
    end do

  end subroutine jacobi


  ! ==========================================================
  ! This routine computes the Nth order Jacobi Polynomials 
  ! (jac) for a vector of positions x on the interval (-1,1),
  ! of length npoints.
  !
  !    See for example the recurrence relations 
  !    in equation 2.5.4 (page 70) in 
  !
  !     "Spectral Methods in Fluid Dynamics",
  !     by C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A.Zang
  !     Springer-Verlag, 1988.
  !
  ! ===========================================================

  function jacobi_polynomials(n, alpha, beta, npoints, x) result(jac)

    integer, intent(in)     :: n         ! order of the Jacobi Polynomial
    real (kind=r8) :: alpha 
    real (kind=r8) :: beta
    integer, intent(in)     :: npoints
    real (kind=r8) :: x(npoints)
    real (kind=r8) :: jac(npoints)

    ! Local variables

    real (kind=r8) :: a1k
    real (kind=r8) :: a2k
    real (kind=r8) :: a3k
    real (kind=r8) :: da2kdx

    real (kind=r8) :: jacp1
    real (kind=r8) :: jacm1
    real (kind=r8) :: jac0
    real (kind=r8) :: xtmp

    real (kind=r8) :: c2,c1,c0
    integer j,k

    c0 = 0.0_r8
    c1 = 1.0_r8
    c2 = 2.0_r8

    do j = 1,npoints

       xtmp=x(j)

       jacm1=c1
       jac0 =(c1+alpha)*xtmp

       do k=1,n-1
          a1k=c2*(k+c1)*(k+alpha+beta+c1)*(c2*k+alpha+beta)
          da2kdx=(c2*k+alpha+beta+c2)*(c2*k+alpha+beta+c1)*(c2*k+alpha+beta)
          a2k=(c2*k+alpha+beta+c1)*(alpha*alpha-beta*beta) + xtmp*da2kdx
          a3k=c2*(k+alpha)*(k+beta)*(c2*k+alpha+beta+c2)
          jacp1=(a2k*jac0-a3k*jacm1)/a1k
          jacm1=jac0
          jac0 =jacp1
       end do

       if (n==0)jac0=jacm1
       jac(j)=jac0
    end do

  end function jacobi_polynomials

  ! ================================================
  ! This routine computes the first derivatives of Nth
  ! order Jacobi Polynomials (djac) for a vector of 
  ! positions x on the interval (-1,1), of length npoints.
  !
  ! See for example the recurrence relations 
  ! in equation 2.5.4 (page 70) in 
  !
  ! "Spectral Methods in Fluid Dynamics",
  ! by C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A.Zang
  ! Springer-Verlag, 1988.
  !
  ! ================================================

  function jacobi_derivatives(n, alpha, beta, npoints, x) result(djac)

    integer                , intent(in) :: n         ! order of the Jacobi Polynomial
    real (kind=r8), intent(in) :: alpha 
    real (kind=r8), intent(in) :: beta
    integer                , intent(in) :: npoints
    real (kind=r8), intent(in) :: x(npoints)

    real (kind=r8)             :: djac(npoints)

    ! Local variables

    ! Local variables

    real (kind=r8) :: a1k
    real (kind=r8) :: a2k
    real (kind=r8) :: a3k
    real (kind=r8) :: da2kdx

    real (kind=r8) :: jacp1
    real (kind=r8) :: jacm1
    real (kind=r8) :: jac0
    real (kind=r8) :: djacp1
    real (kind=r8) :: djacm1
    real (kind=r8) :: djac0

    real (kind=r8) :: xtmp

    real (kind=r8) :: c2,c1,c0
    integer j,k

    c0 = 0.0_r8
    c1 = 1.0_r8
    c2 = 2.0_r8

    do j = 1,npoints

       xtmp=x(j)

       jacm1=c1
       jac0 =(c1+alpha)*xtmp

       djacm1 = c0
       djac0  = (c1+alpha)

       do k=1,n-1
          a1k=c2*(k+c1)*(k+alpha+beta+c1)*(c2*k+alpha+beta)
          da2kdx=(c2*k+alpha+beta+c2)*(c2*k+alpha+beta+c1)*(c2*k+alpha+beta)
          a2k=(c2*k+alpha+beta+c1)*(alpha*alpha-beta*beta) + xtmp*da2kdx
          a3k=c2*(k+alpha)*(k+beta)*(c2*k+alpha+beta+c2)

          jacp1=(a2k*jac0-a3k*jacm1)/a1k
          djacp1=(a2k*djac0+da2kdx*jac0-a3k*djacm1)/a1k

          jacm1=jac0
          jac0=jacp1

          djacm1=djac0
          djac0=djacp1

       end do

       if (n==0)djac0=djacm1
       djac(j)=djac0

    end do

  end function jacobi_derivatives

  ! ===================================================
  !
  ! legendre:
  !
  ! Compute the legendre polynomials using
  ! the recurrence relationship.
  ! return leg(m+1) = P_N(x) for m=0..N
  ! p_3 = Legendre polynomial of degree N
  ! p_2 = Legendre polynomial of degree N-1 at x
  ! p_1 = Legendre polynomial of degree N-2 at x
  !
  ! ===================================================

  function legendre(x,N) result(leg)

    integer   :: N
    real (kind=r8) :: x
    real (kind=r8) :: leg(N+1)

    real (kind=r8) ::  p_1, p_2, p_3
    integer   :: k

    p_3 = 1.0_r8
    leg(1)=p_3
    if (n.ne.0) then
       p_2 = p_3
       p_3 = x
       leg(2)=p_3
       do k = 2,N
          p_1 = p_2
          p_2 = p_3
          p_3 = ( (2*k-1)*x*p_2 - (k-1)*p_1 ) / k
          leg(k+1)=p_3
       end do
    end if

  end function legendre


  ! ===========================================
  ! quad_norm:
  !
  ! compute normalization constants 
  ! for k=1,N order Legendre polynomials
  !
  ! e.g. gamma(k) in Canuto, page 58.
  !
  ! ===========================================

  function quad_norm(gquad,N) result(gamma)
    type (quadrature_t), intent(in) :: gquad
    integer            , intent(in) :: N

    real (kind=r8) :: gamma(N)

    ! Local variables
    real (kind=r8) :: leg(N)
    integer               :: i,k

    gamma(:)=0.0_r8

    do i=1,N
       leg=legendre(gquad%points(i),N-1)
       do k=1,N
          gamma(k)= gamma(k)+leg(k)*leg(k)*gquad%weights(i)
       end do
    end do

  end function quad_norm

  ! =======================
  ! TrapN:
  ! Numerical recipes
  ! =======================

  subroutine trapN(f,a,b,N,it,s)
    INTERFACE
       FUNCTION f(x) RESULT(f_x)   ! Function to be integrated
         use shr_kind_mod, only: r8=>shr_kind_r8
         real(kind=r8), INTENT(IN) :: x
         real(kind=r8) :: f_x
       END FUNCTION f
    END INTERFACE

    real(kind=r8),intent(in) :: a,b
    integer, intent(in)             :: N
    integer, intent(inout)          :: it
    real(kind=r8), intent(inout) :: s

    real(kind=r8) :: ssum
    real(kind=r8) :: del
    real(kind=r8) :: rtnm
    real(kind=r8) :: x

    integer :: j

    if (N==1) then
       s = 0.5_r8*(b-a)*(f(a) + f(b))
       it =1
    else
       ssum = 0.0_r8
       rtnm =1.0_r8/it
       del = (b-a)*rtnm
       x=a+0.5_r8*del
       do j=1,it
          ssum = ssum + f(x)
          x=x+del
       end do
       s=0.5_r8*(s + del*ssum)
       it=2*it  
    end if

  end subroutine trapN

  ! ==========================================
  ! Trapezoid Rule for integrating functions 
  ! from a to b with residual error eps
  ! ==========================================

  function trapezoid(f,a,b,eps) result(Integral)

    integer, parameter :: Nmax = 25  ! At most 2^Nmax + 1 points in integral

    INTERFACE
       FUNCTION f(x) RESULT(f_x)   ! Function to be integrated
         use shr_kind_mod, only: r8=>shr_kind_r8
         real(kind=r8), INTENT(IN) :: x
         real(kind=r8) :: f_x
       END FUNCTION f
    END INTERFACE

    real(kind=r8), intent(in) :: a,b       ! The integral bounds
    real(kind=r8), intent(in) :: eps       ! relative error bound for integral
    real(kind=r8)             :: Integral  ! the integral result (within eps)
    real(kind=r8)             :: s         ! Integral approximation
    real(kind=r8)             :: sold      ! previous integral approx

    integer                          :: N
    integer                          :: it

    ! ==============================================================
    ! Calculate I here using trapezoid rule using f and a DO loop...
    ! ==============================================================

    s    = 1.0e30_r8
    sold = 0.0_r8
    N=1
    it=0
    do while(N<=Nmax .and. ABS(s-sold)>eps*ABS(sold))
       sold=s
       call trapN(f,a,b,N,it,s)
       N=N+1
    end do

    Integral = s

  end function trapezoid

  ! ==========================================
  ! Simpsons Rule for integrating functions 
  ! from a to b with residual error eps
  ! ==========================================

  function simpsons(f,a,b,eps) result(Integral)

    integer, parameter :: Nmax = 25  ! At most 2^Nmax + 1 points in integral

    INTERFACE
       FUNCTION f(x) RESULT(f_x)   ! Function to be integrated
         use shr_kind_mod, only: r8=>shr_kind_r8
         real(kind=r8), INTENT(IN) :: x
         real(kind=r8) :: f_x
       END FUNCTION f
    END INTERFACE

    real(kind=r8), intent(in) :: a,b       ! The integral bounds
    real(kind=r8), intent(in) :: eps       ! relative error bound for integral
    real(kind=r8)             :: Integral  ! the integral result (within eps)
    real(kind=r8)             :: s         ! Integral approximation
    real(kind=r8)             :: os        ! previous integral approx
    real(kind=r8)             :: st        ! Integral approximation
    real(kind=r8)             :: ost       ! previous integral approx

    integer                          :: N
    integer                          :: it

    ! ==============================================================
    ! Calculate I here using trapezoid rule using f and a DO loop...
    ! ==============================================================

    ost= 0.0_r8
    s  = 1.0e30_r8
    os = 0.0_r8

    N=1
    it=0
    do while ((N<=Nmax .and. ABS(s-os)>eps*ABS(os) ) .or. N<=2)
       os = s
       call trapN(f,a,b,N,it,st)
       s=(4.0_r8*st-ost)/3.0_r8
       ost=st
       N=N+1
    end do

    Integral = s

  end function simpsons


  ! ==========================================
  ! gaussian_int:
  !
  ! Gaussian Quadrature Rule for integrating 
  ! function f from a to b  with gs weights and 
  ! points with precomputed gaussian quadrature 
  ! and weights.
  ! ==========================================

  function gaussian_int(f,a,b,gs) result(Integral)

    integer, parameter :: Nmax = 10  ! At most 2^Nmax + 1 points in integral

    INTERFACE
       FUNCTION f(x) RESULT(f_x)   ! Function to be integrated
         use shr_kind_mod, only: r8=>shr_kind_r8
         real(kind=r8), INTENT(IN) :: x
         real(kind=r8) :: f_x
       END FUNCTION f
    END INTERFACE

    real(kind=r8), intent(in) :: a,b       ! The integral bounds
    type(quadrature_t), intent(in)   :: gs        ! gaussian points/wts
    real(kind=r8)             :: Integral  ! the integral result (within eps)

    integer                          :: i
    real (kind=r8)            :: s,x
    ! ==============================================================
    ! Calculate I = S f(x)dx here using gaussian quadrature
    ! ==============================================================

    s = 0.0_r8
    do i=1,SIZE(gs%points)
       x = 0.50_r8*((b-a)*gs%points(i) + (b+a))
       s = s + gs%weights(i)*f(x)
    end do
    Integral = s*(0.5_r8*(b-a))

  end function gaussian_int

end module quadrature_mod





