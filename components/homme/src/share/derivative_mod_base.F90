#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


module derivative_mod_base

  use kinds,          only : real_kind, longdouble_kind
  use dimensions_mod, only : np, nc, nelemd, nlev
  use quadrature_mod, only : quadrature_t, gauss, gausslobatto,legendre, jacobi
  use parallel_mod,   only : abortmp
  use element_mod,    only : element_t
  use control_mod,    only : hypervis_scaling, hypervis_power
  use physical_constants, only : rrearth

implicit none
private

  type, public :: derivative_t
     real (kind=real_kind) :: Dvv(np,np)
     real (kind=real_kind) :: Dvv_diag(np,np)
     real (kind=real_kind) :: Dvv_twt(np,np)
     real (kind=real_kind) :: Mvv_twt(np,np)  ! diagonal matrix of GLL weights
     real (kind=real_kind) :: Mfvm(np,nc+1)
     real (kind=real_kind) :: Cfvm(np,nc)
     real (kind=real_kind) :: legdg(np,np)
  end type derivative_t

  type, public :: derivative_stag_t
     real (kind=real_kind) :: D(np,np)
     real (kind=real_kind) :: M(np,np)
     real (kind=real_kind) :: Dpv(np,np)
     real (kind=real_kind) :: D_twt(np,np)
     real (kind=real_kind) :: M_twt(np,np)
     real (kind=real_kind) :: M_t(np,np)
  end type derivative_stag_t

  real (kind=real_kind), allocatable :: integration_matrix(:,:)
  real (kind=real_kind), allocatable :: boundary_interp_matrix(:,:,:)

! ======================================
! Public Interfaces
! ======================================

  public :: subcell_integration
  public :: subcell_dss_fluxes
  public :: subcell_div_fluxes
  public :: subcell_Laplace_fluxes
  public :: allocate_subcell_integration_matrix

  public :: derivinit
  public :: deriv_print

  public :: gradient
  public :: gradient_wk
  public :: vorticity
  public :: divergence

  public :: interpolate_gll2fvm_corners
  public :: interpolate_gll2fvm_points
  public :: remap_phys2gll


  interface divergence
      module procedure divergence_nonstag
      module procedure divergence_stag
  end interface

  interface gradient
      module procedure gradient_str_nonstag
      module procedure gradient_str_stag
  end interface

  interface gradient_wk
      module procedure gradient_wk_nonstag
      module procedure gradient_wk_stag
  end interface

  public :: v2pinit

  private :: dmatinit
  private :: dvvinit
  private :: dpvinit

! these routines compute spherical differential operators as opposed to
! the gnomonic coordinate operators above.  Vectors (input or output)
! are always expressed in lat-lon coordinates
!
! note that weak derivatives (integrated by parts form) can be defined using
! contra or co-variant test functions, so 
!
  public  :: gradient_sphere
  public  :: gradient_sphere_wk_testcov
  public  :: gradient_sphere_wk_testcontra   ! only used for debugging
  public  :: ugradv_sphere
  public  :: vorticity_sphere
  public  :: vorticity_sphere_diag
  public  :: divergence_sphere
  public  :: curl_sphere
  public  :: curl_sphere_wk_testcov
! public  :: curl_sphere_wk_testcontra  ! not coded
  public  :: divergence_sphere_wk
  public  :: laplace_sphere_wk
  public  :: vlaplace_sphere_wk
  public  :: element_boundary_integral
  public  :: edge_flux_u_cg
  public  :: limiter_optim_iter_full

contains

! ==========================================
! derivinit:
!
! Initialize the matrices for taking 
! derivatives and interpolating
! ==========================================

  subroutine derivinit(deriv,fvm_corners, fvm_points)
    type (derivative_t)      :: deriv
!    real (kind=longdouble_kind),optional :: phys_points(:)
    real (kind=longdouble_kind),optional :: fvm_corners(nc+1)
    real (kind=longdouble_kind),optional :: fvm_points(nc)

    ! Local variables
    type (quadrature_t) :: gp   ! Quadrature points and weights on pressure grid
    
    real (kind=longdouble_kind) :: dmat(np,np)
    real (kind=longdouble_kind) :: dpv(np,np)
    real (kind=longdouble_kind) :: v2p(np,np)
    real (kind=longdouble_kind) :: p2v(np,np)
    real (kind=longdouble_kind) :: dvv(np,np)
    real (kind=longdouble_kind) :: dvv_diag(np,np)
    real (kind=longdouble_kind) :: v2v(np,np)
    real (kind=longdouble_kind) :: xnorm
    integer i,j

    ! ============================================
    ! initialize matrices in longdouble_kind precision
    ! and transfer results into real_kind
    ! floating point precision
    ! ============================================

    gp=gausslobatto(np)

    call dvvinit(dvv,gp)
    deriv%Dvv(:,:)   = dvv(:,:)

    do i=1,np
       do j=1,np
          if (i.eq.j) then
             deriv%dvv_diag(i,j)   = dvv(i,j)
          else
             deriv%dvv_diag(i,j) = 0.0D0
          endif 
        end do
     end do


    v2v = 0.0D0
    do i=1,np
       v2v(i,i) = gp%weights(i)
    end do

    do i=1,np
       do j=1,np
          dvv(j,i) = dvv(j,i)*gp%weights(i)
       end do
    end do

    deriv%Dvv_twt = TRANSPOSE(dvv)
    deriv%Mvv_twt = v2v
    if (present(fvm_corners)) &
         call v2pinit(deriv%Mfvm,gp%points,fvm_corners,np,nc+1)

    if (present(fvm_points)) &
         call v2pinit(deriv%Cfvm,gp%points,fvm_points,np,nc)

    ! notice we deallocate this memory here even though it was allocated 
    ! by the call to gausslobatto.
    deallocate(gp%points)
    deallocate(gp%weights)

  end subroutine derivinit

  subroutine deriv_print(deriv)
    type (derivative_t), intent(in) :: deriv
    
    ! Local variables

    integer j
    print *,"Derivative Matrix Dvv"
    do j=1,np
       write(6,*)deriv%Dvv(:,j)
    end do

    print *,"Weak Derivative Matrix Dvv_twt"
    do j=1,np
       write(6,*)deriv%Dvv_twt(:,j)
    end do


  end subroutine deriv_print

! =======================================
! dmatinit:
!
! Compute rectangular v->p 
! derivative matrix (dmat)
! =======================================

  subroutine dmatinit(dmat)

    real (kind=longdouble_kind) :: dmat(np,np)

    ! Local variables

    type (quadrature_t) :: gll
    type (quadrature_t) :: gs

    integer i,j
    real(kind=longdouble_kind)  fact,f1,f2
    real(kind=longdouble_kind)  func0,func1
    real(kind=longdouble_kind)  dis,c0,c1

    real(kind=longdouble_kind)  :: leg(np,np)
    real(kind=longdouble_kind)  ::  jac(0:np-1)
    real(kind=longdouble_kind)  :: djac(0:np-1)

    c0 = 0.0_longdouble_kind
    c1 = 1.0_longdouble_kind

    gll= gausslobatto(np)
    gs = gauss(np)

    ! =============================================================
    ! Compute Legendre polynomials on Gauss-Lobatto grid (velocity)
    ! =============================================================

    do i=1,np
       leg(:,i) = legendre(gll%points(i),np-1)
    end do

    ! ================================================================
    !  Derivatives of velocity cardinal functions on pressure grid
    !  d(i,j) = D(j,i) = D' (D-transpose) since D(i,j) = dh_j(x_i)/dx
    ! ================================================================

    fact = np*(np-1)

    do j=1,np
       call jacobi(np-1,gs%points(j),c0,c0,jac(0:np-1),djac(0:np-1))
       func0 =  jac(np-1)
       func1 = djac(np-1)
       f1 = fact*func0
       f2 = (c1 - gs%points(j))*(c1 + gs%points(j)) * func1
       do i = 1, np
          if ( gs%points(j) /= gll%points(i) ) then
             dis = gs%points(j) - gll%points(i)
             dmat(i,j) = func0 / ( leg(np,i)*dis ) + f2 / (fact*leg(np,i)*dis*dis)
!!! OTHER             dmat(i,j) = (1.0D0/(fact*leg(np,i)*dis*dis))* (func0*fact*dis + f2)
          else
             dmat(i,j) = c0
          endif
       end do
    end do

    deallocate(gll%points)
    deallocate(gll%weights)

    deallocate(gs%points)
    deallocate(gs%weights)

end subroutine dmatinit

! =======================================
! dpvinit:
!
! Compute rectangular p->v
! derivative matrix (dmat) 
! for strong gradients
! =======================================

subroutine dpvinit(dmat)

real (kind=longdouble_kind) :: dmat(np,np)

! Local variables

type (quadrature_t) :: gll
type (quadrature_t) :: gs

integer i,j
real(kind=longdouble_kind)  dis,c0,c1

real(kind=longdouble_kind)  :: legv(0:np,np)
real(kind=longdouble_kind)  :: dlegv(0:np,np)

real(kind=longdouble_kind)  :: leg(0:np)
real(kind=longdouble_kind)  :: dleg(0:np)

c0 = 0.0_longdouble_kind
c1 = 1.0_longdouble_kind

gll= gausslobatto(np)
gs = gauss(np)

! =============================================================
! Compute Legendre polynomials on Gauss-Lobatto grid (velocity)
! =============================================================

do i=1,np
call jacobi(np,gll%points(i),c0,c0,legv(0:np,i),dlegv(0:np,i))
end do

! ================================================================
!  Derivatives of velocity cardinal functions on pressure grid
    !  d(i,j) = D(j,i) = D' (D-transpose) since D(i,j) = dh_j(x_i)/dx
    ! ================================================================

    do j=1,np
       call jacobi(np,gs%points(j),c0,c0,leg(0:np),dleg(0:np))
       do i = 1, np
          if ( gs%points(j) /= gll%points(i) ) then
             dis = gll%points(i) - gs%points(j)
             dmat(j,i) = dlegv(np,i)/( dleg(np)*dis ) -  legv(np,i)/ (dleg(np)*dis*dis)
          else
             dmat(j,i) = c0
          endif
       end do
    end do

    deallocate(gll%points)
    deallocate(gll%weights)

    deallocate(gs%points)
    deallocate(gs%weights)

  end subroutine dpvinit

! =======================================
! v2pinit:
! Compute interpolation matrix from gll(1:n1) -> gs(1:n2)
! =======================================
  subroutine v2pinit(v2p,gll,gs,n1,n2)
    integer :: n1,n2
    real(kind=real_kind)  ::  v2p(n1,n2)
    real(kind=real_kind)  ::  v2p_new(n1,n2)
    real(kind=longdouble_kind)  ::  gll(n1),gs(n2)
    ! Local variables

    integer i,j,k,m,l
    real(kind=longdouble_kind)  fact,f1, sum
    real(kind=longdouble_kind)  func0,func1

    real(kind=longdouble_kind)  :: leg(n1,n1)
    real(kind=longdouble_kind)  ::  jac(0:n1-1)
    real(kind=longdouble_kind)  :: djac(0:n1-1)
    real(kind=longdouble_kind)  :: c0,c1

    type(quadrature_t) :: gll_pts
    real(kind=longdouble_kind)  :: leg_out(n1,n2)
    real(kind=longdouble_kind)  :: gamma(n1)

    c0 = 0.0_longdouble_kind
    c1 = 1.0_longdouble_kind

    ! ==============================================================
    ! Compute Legendre polynomials on Gauss-Lobatto grid (velocity)
    ! ==============================================================

    fact = -n1*(n1-1)
    do i=1,n1
       leg(:,i) = legendre(gll(i),n1-1)
       leg(n1,i) = fact * leg(n1,i)
    end do

    ! ===================================================
    !  Velocity cardinal functions on pressure grid
    ! ===================================================
#if 0
    do j=1,n2
       call jacobi(n1-1,gs(j),c0,c0,jac(0:n1-1),djac(0:n1-1))
       func0 = jac(n1-1)
       func1 = djac(n1-1)
       f1 = (c1 - gs(j)**2) * func1
       do i = 1, n1
          if ( gs(j) /= gll(i) ) then
             v2p(i,j) = f1 / ( leg(n1,i) * (gs(j)-gll(i)))
          else
             v2p(i,j) = c1
          endif
       end do
    end do
#endif

    ! NEW VERSION, with no division by (gs(j)-gll(i)):

    ! compute legendre polynomials at output points:
    gll_pts = gausslobatto(n1)

    fact = -n1*(n1-1)
    do i=1,n2
       leg_out(:,i) = legendre(gs(i),n1-1)
       leg_out(n1,i) = fact * leg_out(n1,i)
    end do


    ! compute gamma: (normalization factor for inv(leg)
    do m=1,n1
       gamma(m)=0
       do i=1,n1
          gamma(m)=gamma(m)+leg(m,i)*leg(m,i)*gll_pts%weights(i) 
       enddo
       gamma(m)=1/gamma(m)
    enddo

    ! compute product of leg_out * inv(leg):
    do j=1,n2   ! this should be fvm points
       do l=1,n1   ! this should be GLL points
          sum=0
          do k=1,n1  ! number of polynomials = number of GLL points
             sum=sum + leg_out(k,j)*gamma(k)*leg(k,l)
          enddo
          v2p_new(l,j) = gll_pts%weights(l)*sum
       enddo
    enddo
    deallocate(gll_pts%points)
    deallocate(gll_pts%weights)

#if 0
    do j=1,n2   ! this should be fvm points
       do l=1,n1   ! this should be GLL points
          print *,l,j,v2p_new(l,j),v2p(l,j)
       enddo
    enddo
    print *,'max error: ',maxval(abs(v2p_new-v2p))
#endif

    v2p=v2p_new
  end subroutine v2pinit



! =======================================
! dvvinit:
!
! Compute rectangular v->v
! derivative matrix (dvv)
! =======================================

  subroutine dvvinit(dvv,gll)

    real(kind=longdouble_kind)  ::  dvv(np,np)
    type (quadrature_t)   :: gll

    ! Local variables

    real(kind=longdouble_kind)  :: leg(np,np)
    real(kind=longdouble_kind)  :: c0,c1,c4

    integer i,j

    c0 = 0.0_longdouble_kind
    c1 = 1.0_longdouble_kind
    c4 = 4.0_longdouble_kind

    do i=1,np
       leg(:,i) = legendre(gll%points(i),np-1)
    end do

    dvv(:,:) = c0
    do j=1,np
       do i=1,j-1
          dvv(j,i) = (c1/(gll%points(i)-gll%points(j)))*leg(np,i)/leg(np,j)
       end do
       dvv(j,j) = c0
       do i=j+1,np
          dvv(j,i) = (c1/(gll%points(i)-gll%points(j)))*leg(np,i)/leg(np,j)
       end do
    end do


    dvv(np,np) = + np*(np-1)/c4
    dvv(1,1)   = - np*(np-1)/c4

  end subroutine dvvinit

!  ================================================
!  divergence_stag: 
!
!  Compute divergence (maps v grid -> p grid)
!  ================================================

  function divergence_stag(v,deriv) result(div)

    real(kind=real_kind), intent(in) :: v(np,np,2)
    type (derivative_stag_t), intent(in) :: deriv

    real(kind=real_kind) :: div(np,np)

    ! Local

    integer i
    integer j
    integer l

    real(kind=real_kind)  sumx00
    real(kind=real_kind)  sumy00

    real(kind=real_kind) :: vtemp(np,np,2)
    

#ifdef DEBUG
    print *, "divergence_stag"
#endif
     do j=1,np
        do l=1,np
 
           sumx00=0.0d0
           sumy00=0.0d0
!DIR$ UNROLL(NP)
           do i=1,np
              sumx00 = sumx00 + deriv%D(i,l  )*v(i,j  ,1)
              sumy00 = sumy00 + deriv%M(i,l  )*v(i,j  ,2)
           enddo
          vtemp(j  ,l  ,1) = sumx00
          vtemp(j  ,l  ,2) = sumy00
       enddo
    enddo
    do j=1,np
       do i=1,np
          sumx00=0.0d0
          sumy00=0.0d0
!DIR$ UNROLL(NP)
          do l=1,np
             sumx00 = sumx00 +  deriv%M(l,j  )*vtemp(l,i  ,1)
             sumy00 = sumy00 +  deriv%D(l,j  )*vtemp(l,i  ,2)
          enddo
          div(i  ,j  ) = sumx00 + sumy00

       enddo
    enddo

  end function divergence_stag

!  ================================================
!  divergence_nonstag: 
!
!  Compute divergence (maps v->v)
!  ================================================

  function divergence_nonstag(v,deriv) result(div)

    real(kind=real_kind), intent(in) :: v(np,np,2)
    type (derivative_t), intent(in) :: deriv

    real(kind=real_kind) :: div(np,np)

    ! Local

    integer i
    integer j
    integer l

    real(kind=real_kind) ::  dudx00
    real(kind=real_kind) ::  dvdy00

    real(kind=real_kind) ::  vvtemp(np,np)

    !write(*,*) "divergence_nonstag"

     do j=1,np
        do l=1,np
           dudx00=0.0d0
           dvdy00=0.0d0
!DIR$ UNROLL(NP)
           do i=1,np
              dudx00 = dudx00 + deriv%Dvv(i,l  )*v(i,j  ,1)
              dvdy00 = dvdy00 + deriv%Dvv(i,l  )*v(j  ,i,2)
           end do

           div(l  ,j  ) = dudx00
           vvtemp(j  ,l  ) = dvdy00
        end do
    end do
    do j=1,np
       do i=1,np
          div(i,j)=div(i,j)+vvtemp(i,j)
       end do
    end do

  end function divergence_nonstag

!  ================================================
!  gradient_wk_stag:
! 
!  Compute the weak form gradient:
!  maps scalar field on the pressure grid to the
!  velocity grid
!  ================================================

  function gradient_wk_stag(p,deriv) result(dp)

    type (derivative_stag_t), intent(in) :: deriv
    real(kind=real_kind), intent(in) :: p(np,np)

    real(kind=real_kind)             :: dp(np,np,2)

    ! Local
      
    integer i
    integer j
    integer l

    real(kind=real_kind)  sumx00,sumx01
    real(kind=real_kind)  sumy00,sumy01

    real(kind=real_kind)  :: vtempt(np,np,2)

#ifdef DEBUG
    print *, "gradient_wk_stag"
#endif
    !JMD ================================
    !JMD 2*np*np*np Flops 
    !JMD ================================

    do j=1,np
       do l=1,np
          sumx00=0.0d0
          sumy00=0.0d0
!DIR$ UNROLL(NP)
           do i=1,np
              sumx00 = sumx00 + deriv%D_twt(i,l  )*p(i,j  )
              sumy00 = sumy00 + deriv%M_twt(i,l  )*p(i,j  )
           enddo
           vtempt(j  ,l  ,1) = sumx00
           vtempt(j  ,l  ,2) = sumy00
        enddo
    enddo
    do j=1,np
       do i=1,np
          sumx00=0.0d0
          sumy00=0.0d0
!DIR$ UNROLL(NP)
          do l=1,np
             sumx00 = sumx00 +  deriv%M_twt(l,j  )*vtempt(l,i  ,1)
             sumy00 = sumy00 +  deriv%D_twt(l,j  )*vtempt(l,i  ,2)
          enddo
          dp(i  ,j  ,1) = sumx00
          dp(i  ,j  ,2) = sumy00
      enddo
    enddo


  end function gradient_wk_stag

!  ================================================
!  gradient_wk_nonstag:
! 
!  Compute the weak form gradient:
!  maps scalar field on the Gauss-Lobatto grid to the
!  weak gradient on the Gauss-Lobbatto grid
!  ================================================

  function gradient_wk_nonstag(p,deriv) result(dp)

    type (derivative_t), intent(in) :: deriv
    real(kind=real_kind), intent(in) :: p(np,np)

    real(kind=real_kind)             :: dp(np,np,2)

    ! Local
      
    integer i
    integer j
    integer l

    real(kind=real_kind)  sumx00
    real(kind=real_kind)  sumy00

    real(kind=real_kind)  :: vvtempt(np,np,2)

    !JMD ================================
    !JMD 2*np*np*np Flops 
    !JMD ================================

!   print *, "gradient_wk_nonstag"
       do j=1,np
          do l=1,np
             sumx00=0.0d0
             sumy00=0.0d0
!DIR$ UNROLL(NP)
             do i=1,np
                sumx00 = sumx00 + deriv%Dvv_twt(i,l  )*p(i,j  )
                sumy00 = sumy00 + deriv%Mvv_twt(i,l  )*p(i,j  )
             end do
             vvtempt(j  ,l  ,1) = sumx00
             vvtempt(j  ,l  ,2) = sumy00
          end do
       end do

       !JMD ================================
       !JMD 2*np*np*np Flops 
       !JMD ================================

       do j=1,np
          do i=1,np
             sumx00=0.0d0
             sumy00=0.0d0
!DIR$ UNROLL(NP)
             do l=1,np
                sumx00 = sumx00 +  deriv%Mvv_twt(l,j  )*vvtempt(l,i  ,1)
                sumy00 = sumy00 +  deriv%Dvv_twt(l,j  )*vvtempt(l,i  ,2)
             end do
             dp(i  ,j  ,1) = sumx00
             dp(i  ,j  ,2) = sumy00
          end do
       end do
  end function gradient_wk_nonstag

!  ================================================
!  gradient_str_stag:
! 
!  Compute the *strong* form gradient:
!  maps scalar field on the pressure grid to the
!  velocity grid
!  ================================================

  function gradient_str_stag(p,deriv) result(dp)

    type (derivative_stag_t), intent(in) :: deriv
    real(kind=real_kind), intent(in) :: p(np,np)

    real(kind=real_kind)             :: dp(np,np,2)

    ! Local
      
    integer i
    integer j
    integer l

    real(kind=real_kind)  sumx00
    real(kind=real_kind)  sumy00

    real(kind=real_kind)  :: vtempt(np,np,2)
#ifdef DEBUG
    print *, "gradient_str_stag"
#endif
    do j=1,np
       do l=1,np
          sumx00=0.0d0
          sumy00=0.0d0
!DIR$ UNROLL(NP)
          do i=1,np
             sumx00 = sumx00 + deriv%Dpv(i,l  )*p(i,j  )
             sumy00 = sumy00 + deriv%M_t(i,l  )*p(i,j  )
          enddo
          vtempt(j  ,l  ,1) = sumx00
          vtempt(j  ,l  ,2) = sumy00
       enddo
    enddo
    do j=1,np
       do i=1,np
          sumx00=0.0d0
          sumy00=0.0d0
!DIR$ UNROLL(NP)
          do l=1,np
             sumx00 = sumx00 +  deriv%M_t(l,j  )*vtempt(l,i  ,1)
             sumy00 = sumy00 +  deriv%Dpv(l,j  )*vtempt(l,i  ,2)
          enddo
          dp(i  ,j  ,1) = sumx00
          dp(i  ,j  ,2) = sumy00
       enddo
    enddo

  end function gradient_str_stag

!  ================================================
!  gradient_str_nonstag:
!
!  Compute the *strong* gradient on the velocity grid
!  of a scalar field on the velocity grid
!  ================================================

  function gradient_str_nonstag(s,deriv) result(ds)

    type (derivative_t), intent(in) :: deriv
    real(kind=real_kind), intent(in) :: s(np,np)

    real(kind=real_kind) :: ds(np,np,2)

    integer i
    integer j
    integer l
    real(kind=real_kind) ::  dsdx00,dsdx01
    real(kind=real_kind) ::  dsdy00,dsdy01
#ifdef DEBUG
    print *, "gradient_str_nonstag"
!   write(17) np,s,deriv
#endif
       do j=1,np
          do l=1,np
             dsdx00=0.0d0
             dsdy00=0.0d0
!DIR$ UNROLL(NP)
             do i=1,np
                dsdx00 = dsdx00 + deriv%Dvv(i,l  )*s(i,j  )
                dsdy00 = dsdy00 + deriv%Dvv(i,l  )*s(j  ,i)
             end do
             ds(l  ,j  ,1) = dsdx00
             ds(j  ,l  ,2) = dsdy00
          end do
       end do
  end function gradient_str_nonstag

!  ================================================
!  vorticity:
!
!  Compute the vorticity of the velocity field on the
!  velocity grid
!  ================================================

  function vorticity(v,deriv) result(vort)

    type (derivative_t), intent(in) :: deriv
    real(kind=real_kind), intent(in) :: v(np,np,2)

    real(kind=real_kind) :: vort(np,np)

    integer i
    integer j
    integer l
    
    real(kind=real_kind) ::  dvdx00,dvdx01
    real(kind=real_kind) ::  dudy00,dudy01

    real(kind=real_kind)  :: vvtemp(np,np)
    do j=1,np
       do l=1,np
          dudy00=0.0d0
          dvdx00=0.0d0
!DIR$ UNROLL(NP)
          do i=1,np
             dvdx00 = dvdx00 + deriv%Dvv(i,l  )*v(i,j  ,2)
             dudy00 = dudy00 + deriv%Dvv(i,l  )*v(j  ,i,1)
          enddo
          vort(l  ,j  ) = dvdx00
          vvtemp(j  ,l  ) = dudy00
       enddo
    enddo
    do j=1,np
       do i=1,np
          vort(i,j)=vort(i,j)-vvtemp(i,j)
       end do
    end do

  end function vorticity




!  ================================================
!  interpolate_gll2fvm_points:
!
!  shape funtion interpolation from data on GLL grid to cellcenters on physics grid
!  Author: Christoph Erath
!  ================================================
  function interpolate_gll2fvm_points(v,deriv) result(p)

    real(kind=real_kind), intent(in) :: v(np,np)
    type (derivative_t)         :: deriv
    real(kind=real_kind) :: p(nc,nc)

    ! Local
    integer i
    integer j
    integer l

    real(kind=real_kind)  sumx00,sumx01
    real(kind=real_kind)  sumx10,sumx11
    real(kind=real_kind)  vtemp(np,nc)

    do j=1,np
       do l=1,nc
          sumx00=0.0d0
!DIR$ UNROLL(NP)
          do i=1,np
             sumx00 = sumx00 + deriv%Cfvm(i,l  )*v(i,j  )
          enddo
          vtemp(j  ,l) = sumx00
        enddo
    enddo
    do j=1,nc
       do i=1,nc
          sumx00=0.0d0
!DIR$ UNROLL(NP)
          do l=1,np
             sumx00 = sumx00 + deriv%Cfvm(l,j  )*vtemp(l,i)
          enddo
          p(i  ,j  ) = sumx00
       enddo
    enddo
  end function interpolate_gll2fvm_points

!  ================================================
!  interpolate_gll2fvm_corners:
!
!  shape funtion interpolation from data on GLL grid to physics grid
!
!  ================================================
  function interpolate_gll2fvm_corners(v,deriv) result(p)

    real(kind=real_kind), intent(in) :: v(np,np)
    type (derivative_t), intent(in) :: deriv
    real(kind=real_kind) :: p(nc+1,nc+1)

    ! Local
    integer i
    integer j
    integer l

    real(kind=real_kind)  sumx00,sumx01
    real(kind=real_kind)  sumx10,sumx11
    real(kind=real_kind)  vtemp(np,nc+1)

    do j=1,np
       do l=1,nc+1
          sumx00=0.0d0
!DIR$ UNROLL(NP)
          do i=1,np
             sumx00 = sumx00 + deriv%Mfvm(i,l  )*v(i,j  )
          enddo
          vtemp(j  ,l) = sumx00
        enddo
    enddo
    do j=1,nc+1
       do i=1,nc+1
          sumx00=0.0d0
!DIR$ UNROLL(NP)
          do l=1,np
             sumx00 = sumx00 + deriv%Mfvm(l,j  )*vtemp(l,i)
          enddo
          p(i  ,j  ) = sumx00
       enddo
    enddo
  end function interpolate_gll2fvm_corners


#if 0
!  ================================================
!  bilin_phys2gll:
!
!  interpolate to an equally spaced (in reference element coordinate system)
!  "physics" grid to the GLL grid
!
!  For edge/corner nodes:  compute only contribution from this element,
!  assuming there is a corresponding (symmetric) contribution from the
!  neighboring element which will be incorporated after DSS'ing the results
!  Note: this symmetry assuption is false at cube panel edges.  
!
!  MT initial version 3/2014
!  2nd order at all gll points (including cube corners) *except*
!  1st order at edge points on cube panel edges
!  ================================================
  function bilin_phys2gll(pin,nphys) result(pout)
    integer :: nphys
    real(kind=real_kind), intent(in) :: pin(nphys,nphys)
    real(kind=real_kind) :: pout(np,np)
    
    ! static data, shared by all threads
    integer, save  :: nphys_init=0
    integer, save  :: index_l(np),index_r(np)
    real(kind=real_kind),save :: w_l(np),w_r(np)

    ! local
    real(kind=real_kind) :: px(np,nphys)  ! interpolate in x to this array
                                          ! then interpolate in y to pout
    integer i,j,i1,i2
    real(kind=real_kind) :: phys_centers(nphys)
    real(kind=real_kind) :: dx,d_l,d_r,gll
    type(quadrature_t) :: gll_pts


    if (nphys_init/=nphys) then
    ! setup (most be done on masterthread only) since all data is static
    ! MT: move inside if - we dont want a barrier every time this is called       
#if (defined HORIZ_OPENMP)
!OMP BARRIER
!OMP MASTER
#endif
       nphys_init=nphys

       ! compute phys grid cell edges on [-1,1]
       do i=1,nphys
          dx = 2d0/nphys
          phys_centers(i)=-1 + (dx/2) + (i-1)*dx
       enddo

       ! compute 1D interpolation weights
       ! for every GLL point, find the index of the phys point to the left and right:
       !    index_l, index_r
       ! Then compute the weights
       !    w_l, w_r
       !
       ! for points on the edge, we just take data from a single point
       ! which we'll fake by setting index_l=index_r and w_l = w_r = 0.5

       ! first GLL point is just a copy from first PHYS point:    
       w_l(1)=0.5d0
       w_r(1)=0.5d0
       index_l(1)=1
       index_r(1)=1

       w_l(np)=0.5d0
       w_r(np)=0.5d0
       index_l(np)=nphys
       index_r(np)=nphys

       ! compute GLL cell edges on [-1,1]
       gll_pts = gausslobatto(np)

       do i=2,np-1
          ! find: i1,i2 such that:
          ! phys_centers(i1) <=  gll(i)  <= phys_centers(i2)
          gll = gll_pts%points(i)
          i1=0
          do j=1,nphys-1
             if (phys_centers(j) <= gll .and. phys_centers(j+1)>gll) then
                i1=j
                i2=j+1
             endif
          enddo
          if (i1==0) call abortmp('ERROR: bilin_phys2gll() bad logic0')

          d_l = gll-phys_centers(i1) 
          d_r = phys_centers(i2) - gll

          if (d_l<0) call abortmp('ERROR: bilin_phys2gll() bad logic1')
          if (d_r<0) call abortmp('ERROR: bilin_phys2gll() bad logic2')

          w_l(i) = d_r/(d_r + d_l) 
          w_r(i) = d_l/(d_r + d_l) 
          index_l(i)=i1
          index_r(i)=i2  
       enddo

       deallocate(gll_pts%points)
       deallocate(gll_pts%weights)

#if (defined HORIZ_OPENMP)
!OMP END MASTER
!OMP BARRIER
#endif
    endif

    ! interpolate i dimension
    do j=1,nphys
       do i=1,np
          ! pin(nphys,nphys) -> px(np,nphys)
          i1=index_l(i)
          i2=index_r(i)
          px(i,j) = w_l(i)*pin(i1,j) + w_r(i)*pin(i2,j)
       enddo
    enddo

    ! interpolate j dimension
    do j=1,np
       do i=1,np
          ! px(np,nphys) -> pout(np,np)
          i1=index_l(j)
          i2=index_r(j)
          pout(i,j) = w_l(j)*px(i,i1) + w_r(j)*px(i,i2)
       enddo
    enddo
    end function bilin_phys2gll
!----------------------------------------------------------------
#endif

!  ================================================
!  remap_phys2gll:
!
!  interpolate to an equally spaced (in reference element coordinate system)
!  "physics" grid to the GLL grid
!
!  1st order, monotone, conservative
!  MT initial version 2013
!  ================================================
  function remap_phys2gll(pin,nphys) result(pout)
    integer :: nphys
    real(kind=real_kind), intent(in) :: pin(nphys*nphys)
    real(kind=real_kind) :: pout(np,np)
    
    ! Local
    integer, save  :: nphys_init=0
    integer, save  :: nintersect
    real(kind=real_kind),save,pointer :: acell(:)  ! arrivial cell index of i'th intersection
    real(kind=real_kind),save,pointer :: dcell(:)  ! departure cell index of i'th intersection
    real(kind=real_kind),save,pointer :: delta(:)  ! length of i'th intersection
    real(kind=real_kind),save,pointer :: delta_a(:)  ! length of arrival cells
    integer in_i,in_j,ia,ja,id,jd,count,i,j
    logical :: found

    real(kind=real_kind) :: tol=1e-13
    real(kind=real_kind) :: weight,x1,x2,dx
!    real(kind=longdouble_kind) :: gll_edges(np+1),phys_edges(nphys+1)
    real(kind=real_kind) :: gll_edges(np+1),phys_edges(nphys+1)
    type(quadrature_t) :: gll_pts
    if (nphys_init/=nphys) then
       ! setup (most be done on masterthread only) since all data is static
       ! MT: move barrier inside if loop - we dont want a barrier every regular call
#if (defined HORIZ_OPENMP)
!OMP BARRIER
!OMP MASTER
#endif
       nphys_init=nphys
       ! find number of intersections
       nintersect = np+nphys-1  ! max number of possible intersections
       allocate(acell(nintersect))
       allocate(dcell(nintersect))
       allocate(delta(nintersect))
       allocate(delta_a(np))

       ! compute phys grid cell edges on [-1,1]
       do i=1,nphys+1
          dx = 2d0/nphys
          phys_edges(i)=-1 + (i-1)*dx
       enddo

       ! compute GLL cell edges on [-1,1]
       gll_pts = gausslobatto(np)
       gll_edges(1)=-1
       do i=2,np
          gll_edges(i) = gll_edges(i-1) + gll_pts%weights(i-1)
       enddo
       gll_edges(np+1)=1
       delta_a=gll_pts%weights
       deallocate(gll_pts%points)
       deallocate(gll_pts%weights)

       count=0
       x1=-1
       do while ( abs(x1-1) > tol )
          ! find point x2 closet to x1 and x2>x1:
          x2=1.1
          do ia=2,np+1
             if (gll_edges(ia)>x1) then
                if ( ( gll_edges(ia)-x1) < (x2-x1) ) then
                   x2=gll_edges(ia)
                endif
             endif
          enddo
          do id=2,nphys+1
             if (phys_edges(id)>x1) then
                if ( ( phys_edges(id)-x1) < (x2-x1) ) then
                   x2=phys_edges(id)
                endif
             endif
          enddo
          print *,'x2=',x2
          if (x2>1+tol) call abortmp('ERROR: did not find next intersection point')
          if (x2<=x1) call abortmp('ERROR: next intersection point did not advance')
          count=count+1
          if (count>nintersect) call abortmp('ERROR: search failuer: nintersect was too small')
          delta(count)=x2-x1
          
          found=.false.
          do ia=1,np
             if (gll_edges(ia) <= x1+tol  .and.  x2-tol <= gll_edges(ia+1)) then
                found=.true.
                acell(count)=ia
             endif
          enddo
          if (.not. found) call abortmp('ERROR: interval search problem')
       
          found=.false.
          do id=1,nphys
             if (phys_edges(id) <= x1+tol .and.  x2-tol <= phys_edges(id+1)) then
                found=.true.
                dcell(count)=id
             endif
          enddo
          if (.not. found) call abortmp('ERROR: interval search problem')
          x1=x2
       enddo
       ! reset to actual number of intersections
       nintersect=count
#if 0
       print *,'gll->phys conservative monotone remap algorithm:'
       print *,'np,nphys,nintersect',np,nphys,nintersect
       print *,'i   [x1,x2]   [acell]   [dcell]'
       x1=-1
       do in_i=1,nintersect
          ia=acell(in_i)
          id=dcell(in_i)
          write(*,'(i3,a,2f10.6,a,a,2f10.6,a,a,2f10.6,a)') in_i,&
               '[',x1,x1+delta(in_i),']',&
               '[',gll_edges(ia),gll_edges(ia+1),']',&
               '[',phys_edges(id),phys_edges(id+1),']'
          x1=x1+delta(in_i)
       enddo

    pout=0
    do in_i = 1,nintersect
       do in_j = 1,nintersect
          ia = acell(in_i)
          ja = acell(in_j)
          id = dcell(in_i)
          jd = dcell(in_j)
          weight = (  delta(in_i)*delta(in_j) ) / ( delta_a(ia)*delta_a(ja))
          pout(ia,ja) = pout(ia,ja) + weight
       enddo
    enddo
    print *,'sum of weights: ',pout(:,:)
    call abortmp(__FILE__)
#endif
#if (defined HORIZ_OPENMP)
!OMP END MASTER
!OMP BARRIER
#endif
    endif

    pout=0
    do in_i = 1,nintersect
       do in_j = 1,nintersect
          ia = acell(in_i)
          ja = acell(in_j)
          id = dcell(in_i)
          jd = dcell(in_j)
          ! mass in intersection region:  value*area_intersect
          ! value_arrival = value*area_intersect/area_arrival
          weight = (  delta(in_i)*delta(in_j) ) / ( delta_a(ia)*delta_a(ja))
          ! accumulate contribution from each intersection region:
          pout(ia,ja) = pout(ia,ja) + weight*pin(id+(jd-1)*nphys)
       enddo
    enddo
    
    end function remap_phys2gll
    
!----------------------------------------------------------------

!DIR$ ATTRIBUTES FORCEINLINE :: gradient_sphere
  function gradient_sphere(s,deriv,Dinv) result(ds)
!
!   input s:  scalar
!   output  ds: spherical gradient of s, lat-lon coordinates
!

    type (derivative_t), intent(in) :: deriv
    real(kind=real_kind), intent(in), dimension(np,np,2,2) :: Dinv
    real(kind=real_kind), intent(in) :: s(np,np)

    real(kind=real_kind) :: ds(np,np,2)

    integer i
    integer j
    integer l

    real(kind=real_kind) ::  dsdx00, dsdy00
    real(kind=real_kind) ::  v1(np,np),v2(np,np)

    do j=1,np
       do l=1,np
          dsdx00=0.0d0
          dsdy00=0.0d0
!DIR$ UNROLL(NP)
          do i=1,np
             dsdx00 = dsdx00 + deriv%Dvv(i,l  )*s(i,j  )
             dsdy00 = dsdy00 + deriv%Dvv(i,l  )*s(j  ,i)
          end do
          v1(l  ,j  ) = dsdx00*rrearth
          v2(j  ,l  ) = dsdy00*rrearth
       end do
    end do
    ! convert covarient to latlon
    do j=1,np
       do i=1,np
          ds(i,j,1)=Dinv(i,j,1,1)*v1(i,j) + Dinv(i,j,2,1)*v2(i,j)
          ds(i,j,2)=Dinv(i,j,1,2)*v1(i,j) + Dinv(i,j,2,2)*v2(i,j)
       enddo
    enddo

    end function gradient_sphere


  function curl_sphere_wk_testcov(s,deriv,elem) result(ds)
!
!   integrated-by-parts gradient, w.r.t. COVARIANT test functions
!   input s:  scalar  (assumed to be s*khat)
!   output  ds: weak curl, lat/lon coordinates
!   
! starting with: 
!   PHIcov1 = (PHI,0)  covariant vector 
!   PHIcov2 = (0,PHI)  covariant vector 
!
!   ds1 = integral[ PHIcov1 dot curl(s*khat) ] 
!   ds2 = integral[ PHIcov2 dot curl(s*khat) ] 
! integrate by parts: 
!   ds1 = integral[ vor(PHIcov1) * s ]       
!   ds2 = integral[ vor(PHIcov1) * s ]
!
!     PHIcov1 = (PHI^mn,0)   
!     PHIcov2 = (0,PHI^mn)
!  vorticity() acts on covariant vectors:
!   ds1 = sum wij g  s_ij 1/g (  (PHIcov1_2)_x  - (PHIcov1_1)_y ) 
!       = -sum wij s_ij  d/dy (PHI^mn )
! for d/dy component, only sum over i=m
!       = -sum  w_mj s_mj   d( PHI^n)(j)
!           j
!
!   ds2 = sum wij g  s_ij 1/g (  (PHIcov2_2)_x  - (PHIcov2_1)_y ) 
!       = +sum wij s_ij  d/dx (PHI^mn )
! for d/dx component, only sum over j=n
!       = +sum  w_in s_in  d( PHI^m)(i)
!           i
!
    type (derivative_t), intent(in) :: deriv
    type (element_t), intent(in) :: elem
    real(kind=real_kind), intent(in) :: s(np,np)

    real(kind=real_kind) :: ds(np,np,2)

    integer i,j,l,m,n
    real(kind=real_kind) ::  dscontra(np,np,2)

    dscontra=0
    do n=1,np
       do m=1,np
!DIR$ UNROLL(NP)
          do j=1,np
             ! phi(n)_y  sum over second index, 1st index fixed at m
             dscontra(m,n,1)=dscontra(m,n,1)-(elem%mp(m,j)*s(m,j)*deriv%Dvv(n,j) )*rrearth
             ! phi(m)_x  sum over first index, second index fixed at n
             dscontra(m,n,2)=dscontra(m,n,2)+(elem%mp(j,n)*s(j,n)*deriv%Dvv(m,j) )*rrearth
          enddo
       enddo
    enddo

    ! convert contra -> latlon 
    do j=1,np
       do i=1,np
          ds(i,j,1)=(elem%D(i,j,1,1)*dscontra(i,j,1) + elem%D(i,j,1,2)*dscontra(i,j,2))
          ds(i,j,2)=(elem%D(i,j,2,1)*dscontra(i,j,1) + elem%D(i,j,2,2)*dscontra(i,j,2))
       enddo
    enddo
    end function curl_sphere_wk_testcov


  function gradient_sphere_wk_testcov(s,deriv,elem) result(ds)
!
!   integrated-by-parts gradient, w.r.t. COVARIANT test functions
!   input s:  scalar
!   output  ds: weak gradient, lat/lon coordinates
!   ds = - integral[ div(PHIcov) s ]
!
!     PHIcov1 = (PHI^mn,0)   
!     PHIcov2 = (0,PHI^mn)
!   div() acts on contra components, so convert test function to contra: 
!     PHIcontra1 =  metinv PHIcov1  = (a^mn,b^mn)*PHI^mn   
!                                     a = metinv(1,1)  b=metinv(2,1)
!
!   ds1 = sum wij g  s_ij 1/g ( g a PHI^mn)_x  + ( g b PHI^mn)_y ) 
!       = sum  wij s_ij  ag(m,n)  d/dx( PHI^mn ) + bg(m,n) d/dy( PHI^mn)
!          i,j 
! for d/dx component, only sum over j=n
!       = sum  w_in s_in  ag(m,n)  d( PHI^m)(i)
!          i
! for d/dy component, only sum over i=m
!       = sum  w_mj s_mj  bg(m,n)  d( PHI^n)(j)
!          j
!  
!
! This formula is identical to gradient_sphere_wk_testcontra, except that
!    g(m,n) is replaced by a(m,n)*g(m,n)   
!  and we have two terms for each componet of ds 
!
!
    type (derivative_t), intent(in) :: deriv
    type (element_t), intent(in) :: elem
    real(kind=real_kind), intent(in) :: s(np,np)

    real(kind=real_kind) :: ds(np,np,2)

    integer i,j,l,m,n
    real(kind=real_kind) ::  dscontra(np,np,2)


    dscontra=0
    do n=1,np
       do m=1,np
!DIR$ UNROLL(NP)
          do j=1,np
             dscontra(m,n,1)=dscontra(m,n,1)-(&
                  (elem%mp(j,n)*elem%metinv(m,n,1,1)*elem%metdet(m,n)*s(j,n)*deriv%Dvv(m,j) ) +&
                  (elem%mp(m,j)*elem%metinv(m,n,2,1)*elem%metdet(m,n)*s(m,j)*deriv%Dvv(n,j) ) &
                  ) *rrearth

             dscontra(m,n,2)=dscontra(m,n,2)-(&
                  (elem%mp(j,n)*elem%metinv(m,n,1,2)*elem%metdet(m,n)*s(j,n)*deriv%Dvv(m,j) ) +&
                  (elem%mp(m,j)*elem%metinv(m,n,2,2)*elem%metdet(m,n)*s(m,j)*deriv%Dvv(n,j) ) &
                  ) *rrearth
          enddo
       enddo
    enddo
    ! convert contra -> latlon 
    do j=1,np
       do i=1,np
          ds(i,j,1)=(elem%D(i,j,1,1)*dscontra(i,j,1) + elem%D(i,j,1,2)*dscontra(i,j,2))
          ds(i,j,2)=(elem%D(i,j,2,1)*dscontra(i,j,1) + elem%D(i,j,2,2)*dscontra(i,j,2))
       enddo
    enddo

    end function gradient_sphere_wk_testcov


  function gradient_sphere_wk_testcontra(s,deriv,elem) result(ds)
!
!   integrated-by-parts gradient, w.r.t. CONTRA test functions
!   input s:  scalar
!   output  ds: weak gradient, lat/lon coordinates
!
!   integral[ div(phivec) s ] = sum  spheremp()* divergence_sphere(phivec) *s
!   ds1 = above formual with phivec=(PHI,0) in CONTRA coordinates
!   ds2 = above formual with phivec=(0,PHI) in CONTRA coordinates
!   
! PHI = (phi,0)
!   s1 =  sum w_ij s_ij g_ij 1/g_ij ( g_ij PHI^mn )x  
!      =  sum w_ij s_ij g_mn dx(PHI^mn)_ij 
!         ij
! because x derivative is zero for j<>n, only have to sum over j=n
!   s1(m,n)  =  sum w_i,n g_mn dx(PHI^m)_i,n s_i,n
!                i
!
    type (derivative_t), intent(in) :: deriv
    type (element_t), intent(in) :: elem
    real(kind=real_kind), intent(in) :: s(np,np)

    real(kind=real_kind) :: ds(np,np,2)

    integer i,j,l,m,n
    real(kind=real_kind) ::  dscov(np,np,2)

    ! debug: 
    real(kind=real_kind) ::  vcontra(np,np,2)
    real(kind=real_kind) ::  v(np,np,2)
    real(kind=real_kind) ::  div(np,np)



    dscov=0
    do n=1,np
       do m=1,np
!DIR$ UNROLL(NP)
          do j=1,np
             ! phi(m)_x  sum over first index, second index fixed at n
             dscov(m,n,1)=dscov(m,n,1)-(elem%mp(j,n)*elem%metdet(m,n)*s(j,n)*deriv%Dvv(m,j) )*rrearth
             ! phi(n)_y  sum over second index, 1st index fixed at m
             dscov(m,n,2)=dscov(m,n,2)-(elem%mp(m,j)*elem%metdet(m,n)*s(m,j)*deriv%Dvv(n,j) )*rrearth
          enddo
       enddo
    enddo

#if 0
    ! slow form, for debugging
    do m=1,np
       do n=1,np
          vcontra=0
          vcontra(m,n,1)=1

          ! contra->latlon:
          v(:,:,1)=(elem%D(:,:,1,1)*vcontra(:,:,1) + elem%D(:,:,1,2)*vcontra(:,:,2))
          v(:,:,2)=(elem%D(:,:,2,1)*vcontra(:,:,1) + elem%D(:,:,2,2)*vcontra(:,:,2))


          ! compute div(metdet phivec) * s
          div = divergence_sphere(v,deriv,elem)
          ! compute integral[ div(phi) * s ]
          ds(m,n,1)=0
          do i=1,np
             do j=1,np
                ds(m,n,1)=ds(m,n,1) + div(i,j)*s(i,j)*elem%spheremp(i,j)
             enddo
          enddo

          vcontra=0
          vcontra(m,n,2)=1

          ! contra->latlon:
          v(:,:,1)=(elem%D(:,:,1,1)*vcontra(:,:,1) + elem%D(:,:,1,2)*vcontra(:,:,2))
          v(:,:,2)=(elem%D(:,:,2,1)*vcontra(:,:,1) + elem%D(:,:,2,2)*vcontra(:,:,2))

          ! compute div(metdet phivec) * s
          div = divergence_sphere(v,deriv,elem)
          ! compute integral[ div(phi) * s ]
          ds(m,n,2)=0
          do i=1,np
             do j=1,np
                ds(m,n,2)=ds(m,n,2) + div(i,j)*s(i,j)*elem%spheremp(i,j)
             enddo
          enddo
       enddo
    enddo
    ! change sign 
    ds=-ds
    print *,'ds,dscov:1 ',ds(1,1,1),dscov(1,1,1),ds(1,1,1)/dscov(1,1,1)
    print *,'ds,dscov:2 ',ds(1,1,2),dscov(1,1,2),ds(1,1,2)/dscov(1,1,2)

    dscov=ds
#endif
    ! convert covariant -> latlon 
    ds(:,:,1)=elem%Dinv(:,:,1,1)*dscov(:,:,1) + elem%Dinv(:,:,2,1)*dscov(:,:,2)
    ds(:,:,2)=elem%Dinv(:,:,1,2)*dscov(:,:,1) + elem%Dinv(:,:,2,2)*dscov(:,:,2)

    end function gradient_sphere_wk_testcontra

  function ugradv_sphere(u,v,deriv,elem) result(ugradv)
!
!   input:  vectors u and v  (latlon coordinates)
!   output: vector  [ u dot grad ] v  (latlon coordinates)
!
    type (derivative_t), intent(in) :: deriv
    type (element_t), intent(in) :: elem
    real(kind=real_kind), intent(in) :: u(np,np,2)
    real(kind=real_kind), intent(in) :: v(np,np,2)

    real(kind=real_kind) :: ugradv(np,np,2)
    real(kind=real_kind) :: dum_cart(np,np,3)

    integer :: component

    ! latlon -> cartesian
    do component=1,3
       ! Summing along the third dimension is a sum over components for each point.
       ! (This is just a faster way of doing a dot product for each grid point,
       ! since reindexing the inputs to use the intrinsic effectively would be
       ! just asking for trouble.)
       dum_cart(:,:,component)=sum( elem%vec_sphere2cart(:,:,component,:)*v(:,:,:) ,3)
!       dum_cart(:,:,component)= elem%vec_sphere2cart(:,:,component,1)*v(:,:,1) + &
!                                elem%vec_sphere2cart(:,:,component,2)*v(:,:,2)
    end do

    ! Do ugradv on the cartesian components.
    do component=1,3
       ! Dot u with the gradient of each component
       dum_cart(:,:,component) = sum( u(:,:,:) * &
            gradient_sphere(dum_cart(:,:,component),deriv,elem%Dinv) ,3)
    enddo

    ! cartesian -> latlon
    do component=1,2
       ! vec_sphere2cart is its own pseudoinverse.
       ugradv(:,:,component) = sum(dum_cart(:,:,:)*elem%vec_sphere2cart(:,:,:,component), 3)
!       ugradv(:,:,component) = dum_cart(:,:,1)*elem%vec_sphere2cart(:,:,1,component) + &
!                               dum_cart(:,:,2)*elem%vec_sphere2cart(:,:,2,component) + &
!                               dum_cart(:,:,3)*elem%vec_sphere2cart(:,:,3,component)
    end do

  end function ugradv_sphere



  function curl_sphere(s,deriv,elem) result(ds)
!
!   input s:  scalar  (assumed to be  s khat)
!   output  curl(s khat) vector in lat-lon coordinates
! 
!   This subroutine can be used to compute divergence free velocity fields,
!   since div(ds)=0
!
!    first compute:  
!    curl(s khat) = (1/jacobian) ( ds/dy, -ds/dx ) in contra-variant coordinates
!    then map to lat-lon
!
    type (derivative_t), intent(in) :: deriv
    type (element_t), intent(in) :: elem
    real(kind=real_kind), intent(in) :: s(np,np)

    real(kind=real_kind) :: ds(np,np,2)

    integer i
    integer j
    integer l

    real(kind=real_kind) ::  dsdx00
    real(kind=real_kind) ::  dsdy00
    real(kind=real_kind) ::  v1(np,np),v2(np,np)
    
    do j=1,np
       do l=1,np
          dsdx00=0.0d0
          dsdy00=0.0d0
!DIR$ UNROLL(NP)
          do i=1,np
             dsdx00 = dsdx00 + deriv%Dvv(i,l  )*s(i,j  )
             dsdy00 = dsdy00 + deriv%Dvv(i,l  )*s(j  ,i)
          end do
          v2(l  ,j  ) = -dsdx00*rrearth
          v1(j  ,l  ) =  dsdy00*rrearth
       end do
    end do
    ! convert contra -> latlon *and* divide by jacobian
    do j=1,np
       do i=1,np
          ds(i,j,1)=(elem%D(i,j,1,1)*v1(i,j) + elem%D(i,j,1,2)*v2(i,j))/elem%metdet(i,j)
          ds(i,j,2)= (elem%D(i,j,2,1)*v1(i,j) + elem%D(i,j,2,2)*v2(i,j))/elem%metdet(i,j)
       enddo
    enddo
 
    end function curl_sphere


!--------------------------------------------------------------------------



  function divergence_sphere_wk(v,deriv,elem) result(div)
!
!   input:  v = velocity in lat-lon coordinates
!   ouput:  div(v)  spherical divergence of v, integrated by parts
!
!   Computes  -< grad(psi) dot v > 
!   (the integrated by parts version of < psi div(v) > )
!
!   note: after DSS, divergence_sphere() and divergence_sphere_wk() 
!   are identical to roundoff, as theory predicts.
!
    real(kind=real_kind), intent(in) :: v(np,np,2)  ! in lat-lon coordinates
    type (derivative_t), intent(in) :: deriv
    type (element_t), intent(in) :: elem
    real(kind=real_kind) :: div(np,np)

    ! Local

    integer i,j,m,n

    real(kind=real_kind) ::  vtemp(np,np,2)
    real(kind=real_kind) ::  ggtemp(np,np,2)
    real(kind=real_kind) ::  gtemp(np,np,2)
    real(kind=real_kind) ::  psi(np,np)
    real(kind=real_kind) :: xtmp

    ! latlon- > contra
    do j=1,np
       do i=1,np
          vtemp(i,j,1)=(elem%Dinv(i,j,1,1)*v(i,j,1) + elem%Dinv(i,j,1,2)*v(i,j,2))
          vtemp(i,j,2)=(elem%Dinv(i,j,2,1)*v(i,j,1) + elem%Dinv(i,j,2,2)*v(i,j,2))
       enddo
    enddo

    do n=1,np
       do m=1,np

          div(m,n)=0
!DIR$ UNROLL(NP)
          do j=1,np
             div(m,n)=div(m,n)-(elem%spheremp(j,n)*vtemp(j,n,1)*deriv%Dvv(m,j) &
                              +elem%spheremp(m,j)*vtemp(m,j,2)*deriv%Dvv(n,j)) &
                              * rrearth
          enddo

#if 0
! debug the above formula using the N^4 slow formulation:
          psi=0
          psi(m,n)=1
          ggtemp=gradient_sphere(psi,deriv,elem%Dinv)
          ! latlon -> covarient
          do j=1,np
             do i=1,np
                gtemp(i,j,1)=(elem%D(i,j,1,1)*ggtemp(i,j,1) + elem%D(i,j,2,1)*ggtemp(i,j,2))
                gtemp(i,j,2)=(elem%D(i,j,1,2)*ggtemp(i,j,1) + elem%D(i,j,2,2)*ggtemp(i,j,2))
             enddo
          enddo
! grad(psi) dot v:
          xtmp=0
          do j=1,np
          do i=1,np
             xtmp=xtmp-elem%spheremv(i,j)*(vtemp(i,j,1)*gtemp(i,j,1)+vtemp(i,j,2)*gtemp(i,j,2))
          enddo
          enddo
          if (abs(xtmp-div(m,n)) > 3e-17) then
             print *,m,n,xtmp,div(m,n),xtmp-div(m,n)
          endif
#endif          
       end do
    end do
    
  end function divergence_sphere_wk



  function element_boundary_integral(v,deriv,elem) result(result)
!
!   input:  v = velocity in lat-lon coordinates
!   ouput:  result(i,j) = contour integral of PHI_ij * v dot normal
!           where PHI_ij = cardinal function at i,j GLL point 
!
!   this routine is used just to check spectral element integration by parts identities
!
    real(kind=real_kind), intent(in) :: v(np,np,2)  ! in lat-lon coordinates
    type (derivative_t), intent(in) :: deriv
    type (element_t), intent(in) :: elem
    real(kind=real_kind) :: result(np,np)

    ! Local
    real(kind=real_kind) :: ucontra(np,np,2)  ! in lat-lon coordinates
    integer i,j

    ! latlon->contra
    do j=1,np
       do i=1,np
          ucontra(i,j,1)=(elem%Dinv(i,j,1,1)*v(i,j,1) + elem%Dinv(i,j,1,2)*v(i,j,2))
          ucontra(i,j,2)=(elem%Dinv(i,j,2,1)*v(i,j,1) + elem%Dinv(i,j,2,2)*v(i,j,2))
       enddo
    enddo

    ! note: GLL weights  weight(i) = Mvv_twt(i,i)
    result=0
    j=1
    do i=1,np
       result(i,j)=result(i,j)-deriv%Mvv_twt(i,i)*elem%metdet(i,j)*ucontra(i,j,2)*rrearth
    enddo
    
    j=np
    do i=1,np
       result(i,j)=result(i,j)+deriv%Mvv_twt(i,i)*elem%metdet(i,j)*ucontra(i,j,2)*rrearth
    enddo
    
    i=1
    do j=1,np
       result(i,j)=result(i,j)-deriv%Mvv_twt(j,j)*elem%metdet(i,j)*ucontra(i,j,1)*rrearth
    enddo
    
    i=np
    do j=1,np
       result(i,j)=result(i,j)+deriv%Mvv_twt(j,j)*elem%metdet(i,j)*ucontra(i,j,1)*rrearth
    enddo
  end function element_boundary_integral



  function edge_flux_u_cg( v,p,pedges, deriv, elem, u_is_contra) result(result)
!
!
!   input:  v = velocity in contra or lat-lon coordinates (CONTINUIOUS)
!           p      = scalar on this element
!           pedges = scalar edge data from neighbor elements
!
!   ouput:  result(i,j) = contour integral of PHI_ij * pstar * v dot normal
!           where PHI_ij = cardinal function at i,j GLL point 
!           pstar = centered or other flux
!
    real(kind=real_kind), intent(in) :: v(np,np,2) 
    real(kind=real_kind), intent(in) :: p(np,np) 
    real(kind=real_kind), intent(in) :: pedges(0:np+1,0:np+1) 
    type (derivative_t), intent(in) :: deriv
    type (element_t), intent(in) :: elem
    real(kind=real_kind) :: result(np,np)
    logical :: u_is_contra

    ! Local
    real(kind=real_kind) :: ucontra(np,np,2)  ! in lat-lon coordinates
    real(kind=real_kind) :: flux,pstar
    integer i,j


    result=0


    if (u_is_contra) then
       ucontra=v
    else
       ! latlon->contra
       do j=1,np
          do i=1,np
             ucontra(i,j,1)=(elem%Dinv(i,j,1,1)*v(i,j,1) + elem%Dinv(i,j,1,2)*v(i,j,2))
             ucontra(i,j,2)=(elem%Dinv(i,j,2,1)*v(i,j,1) + elem%Dinv(i,j,2,2)*v(i,j,2))
          enddo
       enddo
    endif
#if 0
    ! centered
    do i=1,np
       j=1
       pstar=(pedges(i,0) + p(i,j) ) /2
       flux = -pstar*ucontra(i,j,2)*( deriv%Mvv_twt(i,i)*elem%metdet(i,j)*rrearth)
       result(i,j)=result(i,j)+flux
       
       j=np
       pstar=(pedges(i   ,np+1) + p(i,j) ) /2
       flux = pstar*ucontra(i,j,2)* ( deriv%Mvv_twt(i,i)*elem%metdet(i,j)*rrearth)
       result(i,j)=result(i,j)+flux
    enddo
    
    do j=1,np
       i=1
       pstar=(pedges(0   ,j   ) + p(i,j) )/2
       flux = -pstar*ucontra(i,j,1)* ( deriv%Mvv_twt(j,j)*elem%metdet(i,j)*rrearth)
       result(i,j)=result(i,j)+flux
       
       i=np  
       pstar=(pedges(np+1,j   ) + p(i,j) ) /2
       flux = pstar*ucontra(i,j,1)* ( deriv%Mvv_twt(j,j)*elem%metdet(i,j)*rrearth)
       result(i,j)=result(i,j)+flux
    end do
#else
    ! upwind
    do i=1,np
       j=1
       pstar=p(i,j)
       if (ucontra(i,j,2)>0) pstar=pedges(i,0)
       flux = -pstar*ucontra(i,j,2)*( deriv%Mvv_twt(i,i)*elem%metdet(i,j)*rrearth)
       result(i,j)=result(i,j)+flux
       
       j=np
       pstar=p(i,j)
       if (ucontra(i,j,2)<0) pstar=pedges(i,np+1)
       flux = pstar*ucontra(i,j,2)* ( deriv%Mvv_twt(i,i)*elem%metdet(i,j)*rrearth)
       result(i,j)=result(i,j)+flux
    enddo
    
    do j=1,np
       i=1
       pstar=p(i,j)
       if (ucontra(i,j,1)>0) pstar=pedges(0,j)
       flux = -pstar*ucontra(i,j,1)* ( deriv%Mvv_twt(j,j)*elem%metdet(i,j)*rrearth)
       result(i,j)=result(i,j)+flux
       
       i=np  
       pstar=p(i,j)
       if (ucontra(i,j,1)<0) pstar=pedges(np+1,j)
       flux = pstar*ucontra(i,j,1)* ( deriv%Mvv_twt(j,j)*elem%metdet(i,j)*rrearth)
       result(i,j)=result(i,j)+flux
    end do
#endif    

  end function edge_flux_u_cg

    
!DIR$ ATTRIBUTES FORCEINLINE :: vorticity_sphere
  function vorticity_sphere(v,deriv,elem) result(vort)
!
!   input:  v = velocity in lat-lon coordinates
!   ouput:  spherical vorticity of v
!

    type (derivative_t), intent(in) :: deriv
    type (element_t), intent(in) :: elem
    real(kind=real_kind), intent(in) :: v(np,np,2)

    real(kind=real_kind) :: vort(np,np)

    integer i
    integer j
    integer l
    
    real(kind=real_kind) ::  dvdx00,dudy00
    real(kind=real_kind) ::  vco(np,np,2)
    real(kind=real_kind) ::  vtemp(np,np)

    ! convert to covariant form
    do j=1,np
       do i=1,np
          vco(i,j,1)=(elem%D(i,j,1,1)*v(i,j,1) + elem%D(i,j,2,1)*v(i,j,2))
          vco(i,j,2)=(elem%D(i,j,1,2)*v(i,j,1) + elem%D(i,j,2,2)*v(i,j,2))
       enddo
    enddo

    do j=1,np
       do l=1,np

          dudy00=0.0d0
          dvdx00=0.0d0

!DIR$ UNROLL(NP)
          do i=1,np
             dvdx00 = dvdx00 + deriv%Dvv(i,l  )*vco(i,j  ,2)
             dudy00 = dudy00 + deriv%Dvv(i,l  )*vco(j  ,i,1)
          enddo
 
          vort(l  ,j  ) = dvdx00
          vtemp(j  ,l  ) = dudy00
       enddo
    enddo

    do j=1,np
       do i=1,np
          vort(i,j)=(vort(i,j)-vtemp(i,j))*(elem%rmetdet(i,j)*rrearth)
       end do
    end do

  end function vorticity_sphere

  function vorticity_sphere_diag(v,deriv,elem) result(vort)
  !
  !   input:  v = velocity in lat-lon coordinates
  !   ouput:  diagonal component of spherical vorticity of v
  !

      type (derivative_t), intent(in) :: deriv
      type (element_t), intent(in) :: elem
      real(kind=real_kind), intent(in) :: v(np,np,2)

      real(kind=real_kind) :: vort(np,np)

      integer i
      integer j
      integer l

      real(kind=real_kind) :: dvdx00,dudy00
      real(kind=real_kind) :: vco(np,np,2)
      real(kind=real_kind) :: vtemp(np,np)
      real(kind=real_kind) :: rdx
      real(kind=real_kind) :: rdy

      ! convert to covariant form
                                                                    
      do j=1,np
         do i=1,np
            vco(i,j,1)=(elem%D(i,j,1,1)*v(i,j,1) + elem%D(i,j,2,1)*v(i,j,2))
            vco(i,j,2)=(elem%D(i,j,1,2)*v(i,j,1) + elem%D(i,j,2,2)*v(i,j,2))
         enddo
      enddo

                                                                                                               
      do j=1,np
         do l=1,np
            dudy00=0.0d0
            dvdx00=0.0d0
!DIR$ UNROLL(NP)
            do i=1,np
               dvdx00 = dvdx00 + deriv%Dvv_diag(i,l)*vco(i,j ,2)
               dudy00 = dudy00 + deriv%Dvv_diag(i,l)*vco(j ,i,1)
            enddo 
            vort(l ,j) = dvdx00 
            vtemp(j ,l) = dudy00
         enddo
      enddo

      do j=1,np
         do i=1,np 
          vort(i,j)=(vort(i,j)-vtemp(i,j))*(elem%rmetdet(i,j)*rrearth)
         end do 
      end do 
     
  end function vorticity_sphere_diag

!DIR$ ATTRIBUTES FORCEINLINE :: divergence_sphere
  function divergence_sphere(v,deriv,elem) result(div)
!
!   input:  v = velocity in lat-lon coordinates
!   ouput:  div(v)  spherical divergence of v
!


    real(kind=real_kind), intent(in) :: v(np,np,2)  ! in lat-lon coordinates
    type (derivative_t), intent(in) :: deriv
    type (element_t), intent(in) :: elem
    real(kind=real_kind) :: div(np,np)

    ! Local

    integer i
    integer j
    integer l

    real(kind=real_kind) ::  dudx00
    real(kind=real_kind) ::  dvdy00
    real(kind=real_kind) ::  gv(np,np,2),vvtemp(np,np)

    ! convert to contra variant form and multiply by g
    do j=1,np
       do i=1,np
          gv(i,j,1)=elem%metdet(i,j)*(elem%Dinv(i,j,1,1)*v(i,j,1) + elem%Dinv(i,j,1,2)*v(i,j,2))
          gv(i,j,2)=elem%metdet(i,j)*(elem%Dinv(i,j,2,1)*v(i,j,1) + elem%Dinv(i,j,2,2)*v(i,j,2))
       enddo
    enddo

    ! compute d/dx and d/dy         
    do j=1,np
       do l=1,np
          dudx00=0.0d0
          dvdy00=0.0d0
!DIR$ UNROLL(NP)
          do i=1,np
             dudx00 = dudx00 + deriv%Dvv(i,l  )*gv(i,j  ,1)
             dvdy00 = dvdy00 + deriv%Dvv(i,l  )*gv(j  ,i,2)
          end do
          div(l  ,j  ) = dudx00
          vvtemp(j  ,l  ) = dvdy00
       end do
    end do

!dir$ simd
    div(:,:)=(div(:,:)+vvtemp(:,:))*(elem%rmetdet(:,:)*rrearth)
    
  end function divergence_sphere


!DIR$ ATTRIBUTES FORCEINLINE :: laplace_sphere_wk
  function laplace_sphere_wk(s,deriv,elem,var_coef) result(laplace)
!
!   input:  s = scalar
!   ouput:  -< grad(PHI), grad(s) >   = weak divergence of grad(s)
!     note: for this form of the operator, grad(s) does not need to be made C0
!            
    real(kind=real_kind), intent(in) :: s(np,np) 
    logical, intent(in) :: var_coef
    type (derivative_t), intent(in) :: deriv
    type (element_t), intent(in) :: elem
    real(kind=real_kind)             :: laplace(np,np)
    real(kind=real_kind)             :: laplace2(np,np)
    integer i,j

    ! Local
    real(kind=real_kind) :: grads(np,np,2), oldgrads(np,np,2)

    grads=gradient_sphere(s,deriv,elem%Dinv)
 
    if (var_coef) then
       if (hypervis_power/=0 ) then
          ! scalar viscosity with variable coefficient
          grads(:,:,1) = grads(:,:,1)*elem%variable_hyperviscosity(:,:)
          grads(:,:,2) = grads(:,:,2)*elem%variable_hyperviscosity(:,:)
       else if (hypervis_scaling /=0 ) then
          ! tensor hv, (3)
          oldgrads=grads
          do j=1,np
             do i=1,np
!JMD                grads(i,j,1) = sum(oldgrads(i,j,:)*elem%tensorVisc(i,j,1,:))
!JMD                grads(i,j,2) = sum(oldgrads(i,j,:)*elem%tensorVisc(i,j,2,:))
                grads(i,j,1) = oldgrads(i,j,1)*elem%tensorVisc(i,j,1,1) + &
                               oldgrads(i,j,2)*elem%tensorVisc(i,j,1,2)
                grads(i,j,2) = oldgrads(i,j,1)*elem%tensorVisc(i,j,2,1) + &
                               oldgrads(i,j,2)*elem%tensorVisc(i,j,2,2)
             end do
          end do
       else
          ! do nothing: constant coefficient viscsoity
       endif
    endif

    ! note: divergnece_sphere and divergence_sphere_wk are identical *after* bndry_exchange
    ! if input is C_0.  Here input is not C_0, so we should use divergence_sphere_wk().  
    laplace=divergence_sphere_wk(grads,deriv,elem)

  end function laplace_sphere_wk

!DIR$ ATTRIBUTES FORCEINLINE :: vlaplace_sphere_wk
  function vlaplace_sphere_wk(v,deriv,elem,var_coef,nu_ratio) result(laplace)
!
!   input:  v = vector in lat-lon coordinates
!   ouput:  weak laplacian of v, in lat-lon coordinates
!
!   logic:
!      tensorHV:     requires cartesian
!      nu_div/=nu:   requires contra formulatino
!
!   One combination NOT supported:  tensorHV and nu_div/=nu then abort
!
    real(kind=real_kind), intent(in) :: v(np,np,2) 
    logical, intent(in) :: var_coef
    type (derivative_t), intent(in) :: deriv
    type (element_t), intent(in) :: elem
    real(kind=real_kind), optional :: nu_ratio
    real(kind=real_kind) :: laplace(np,np,2)


    if (hypervis_scaling/=0 .and. var_coef) then
       ! tensorHV is turned on - requires cartesian formulation
       if (present(nu_ratio)) then
          if (nu_ratio /= 1) then
             call abortmp('ERROR: tensorHV can not be used with nu_div/=nu')
          endif
       endif
       laplace=vlaplace_sphere_wk_cartesian(v,deriv,elem,var_coef)
    else  
       ! all other cases, use contra formulation:
       laplace=vlaplace_sphere_wk_contra(v,deriv,elem,var_coef,nu_ratio)
    endif

  end function vlaplace_sphere_wk



  function vlaplace_sphere_wk_cartesian(v,deriv,elem,var_coef) result(laplace)
!
!   input:  v = vector in lat-lon coordinates
!   ouput:  weak laplacian of v, in lat-lon coordinates

    real(kind=real_kind), intent(in) :: v(np,np,2) 
    logical :: var_coef
    type (derivative_t), intent(in) :: deriv
    type (element_t), intent(in) :: elem
    real(kind=real_kind) :: laplace(np,np,2)
    ! Local

    integer component
    real(kind=real_kind) :: dum_cart(np,np,3)


    ! latlon -> cartesian
    do component=1,3
!JMD       dum_cart(:,:,component)=sum( elem%vec_sphere2cart(:,:,component,:)*v(:,:,:) ,3)
       dum_cart(:,:,component) = elem%vec_sphere2cart(:,:,component,1)*v(:,:,1) + &
                                elem%vec_sphere2cart(:,:,component,2)*v(:,:,2)
    end do

    ! Do laplace on cartesian comps
    do component=1,3
       dum_cart(:,:,component) = laplace_sphere_wk(dum_cart(:,:,component),deriv,elem,var_coef)
    enddo

    ! cartesian -> latlon
    do component=1,2
       ! vec_sphere2cart is its own pseudoinverse.
!JMD       laplace(:,:,component)=sum( dum_cart(:,:,:)*elem%vec_sphere2cart(:,:,:,component) ,3)
       laplace(:,:,component) = dum_cart(:,:,1)*elem%vec_sphere2cart(:,:,1,component) + &
                                dum_cart(:,:,2)*elem%vec_sphere2cart(:,:,2,component) + &
                                dum_cart(:,:,3)*elem%vec_sphere2cart(:,:,3,component) 
    end do 

#define UNDAMPRRCART
#ifdef UNDAMPRRCART
    ! add in correction so we dont damp rigid rotation
    laplace(:,:,1)=laplace(:,:,1) + 2*elem%spheremp(:,:)*v(:,:,1)*(rrearth**2)
    laplace(:,:,2)=laplace(:,:,2) + 2*elem%spheremp(:,:)*v(:,:,2)*(rrearth**2)
#endif
  end function vlaplace_sphere_wk_cartesian



  function vlaplace_sphere_wk_contra(v,deriv,elem,var_coef,nu_ratio) result(laplace)
!
!   input:  v = vector in lat-lon coordinates
!   ouput:  weak laplacian of v, in lat-lon coordinates
!
!   du/dt = laplace(u) = grad(div) - curl(vor)
!   < PHI du/dt > = < PHI laplace(u) >        PHI = covariant, u = contravariant
!                 = < PHI grad(div) >  - < PHI curl(vor) >
!                 = grad_wk(div) - curl_wk(vor)               
!
    real(kind=real_kind), intent(in) :: v(np,np,2) 
    logical, intent(in) :: var_coef
    type (derivative_t), intent(in) :: deriv
    type (element_t), intent(in) :: elem
    real(kind=real_kind) :: laplace(np,np,2)
    real(kind=real_kind), optional :: nu_ratio
    ! Local

    integer i,j,l,m,n
    real(kind=real_kind) :: vor(np,np),div(np,np)
    real(kind=real_kind) :: v1,v2,div1,div2,vor1,vor2,phi_x,phi_y

    div=divergence_sphere(v,deriv,elem)
    vor=vorticity_sphere(v,deriv,elem)

    if (var_coef .and. hypervis_power/=0 ) then
          ! scalar viscosity with variable coefficient
          div = div*elem%variable_hyperviscosity(:,:)
          vor = vor*elem%variable_hyperviscosity(:,:)
    endif

    if (present(nu_ratio)) div = nu_ratio*div

    laplace = gradient_sphere_wk_testcov(div,deriv,elem) - &
         curl_sphere_wk_testcov(vor,deriv,elem)

    do n=1,np
       do m=1,np
          ! add in correction so we dont damp rigid rotation
#define UNDAMPRR
#ifdef UNDAMPRR
          laplace(m,n,1)=laplace(m,n,1) + 2*elem%spheremp(m,n)*v(m,n,1)*(rrearth**2)
          laplace(m,n,2)=laplace(m,n,2) + 2*elem%spheremp(m,n)*v(m,n,2)*(rrearth**2)
#endif
       enddo
    enddo
  end function vlaplace_sphere_wk_contra


  function subcell_dss_fluxes(dss, p, n, metdet, C) result(fluxes)

    implicit none

    integer              , intent(in)  :: p
    integer              , intent(in)  :: n
    real (kind=real_kind), intent(in)  :: dss     (p,p)
    real (kind=real_kind), intent(in)  :: metdet  (p,p)
    real (kind=real_kind), intent(in)  :: C       (2,2,2)  

    real (kind=real_kind)              :: fluxes  (n,n,4)

    real (kind=real_kind)              :: Bp(p,p)
    real (kind=real_kind)              :: Tp(p,p)
    real (kind=real_kind)              :: Lp(p,p)
    real (kind=real_kind)              :: Rp(p,p)

    real (kind=real_kind)              :: B(n,n)
    real (kind=real_kind)              :: T(n,n)
    real (kind=real_kind)              :: L(n,n)
    real (kind=real_kind)              :: R(n,n)

    integer            :: i,j

    fluxes  = 0

    Bp = 0
    Tp = 0
    Rp = 0
    Lp = 0

    Bp(:,1)  = dss(:,1)  ! bottom
    Tp(:,p)  = dss(:,p)  ! top
    Rp(p,:)  = dss(p,:)  ! right
    Lp(1,:)  = dss(1,:)  ! left

    Bp(1,1)  = C(1,1,2)
    Lp(1,1)  = C(1,1,1)
    Bp(p,1)  = C(2,1,2)
    Rp(p,1)  = C(2,1,1) 

    Tp(1,p)  = C(1,2,2)
    Lp(1,p)  = C(1,2,1) 
    Tp(p,p)  = C(2,2,2)
    Rp(p,p)  = C(2,2,1) 

    B = subcell_integration(Bp, p, n, metdet)
    T = subcell_integration(Tp, p, n, metdet)
    L = subcell_integration(Lp, p, n, metdet)
    R = subcell_integration(Rp, p, n, metdet)

    do i = 1,n
    do j = 1,n
      if (1<j) T(i,j) = T(i,j) + T(i,j-1) 
      if (1<i) R(i,j) = R(i,j) + R(i-1,j) 
    end do
    end do

    do i = n,1,-1
    do j = n,1,-1
      if (j<n) B(i,j) = B(i,j) + B(i,j+1) 
      if (i<n) L(i,j) = L(i,j) + L(i+1,j) 
    end do
    end do

    do i = 1,n
      do j = 1,n
        if (1==j) fluxes(i,j,1) =  B(i,j)
        if (n==i) fluxes(i,j,2) =  R(i,j)
        if (j==n) fluxes(i,j,3) =  T(i,j)
        if (1==i) fluxes(i,j,4) =  L(i,j)

        if (1< j) fluxes(i,j,1) =   B(i,j) - T(i,j-1)
        if (i< n) fluxes(i,j,2) =   R(i,j) - L(i+1,j)
        if (j< n) fluxes(i,j,3) =   T(i,j) - B(i,j+1)
        if (1< i) fluxes(i,j,4) =   L(i,j) - R(i-1,j)
      end do
    end do

  end function subcell_dss_fluxes

  function subcell_div_fluxes(u, p, n, metdet) result(fluxes)

    implicit none

    integer              , intent(in)  :: p
    integer              , intent(in)  :: n
    real (kind=real_kind), intent(in)  :: u(p,p,2)
    real (kind=real_kind), intent(in)  :: metdet(p,p)

    real (kind=real_kind)              :: v(p,p,2)
    real (kind=real_kind)              :: fluxes(n,n,4)
    real (kind=real_kind)              :: tb(n,p)
    real (kind=real_kind)              :: lr(p,n)
    real (kind=real_kind)              :: flux_l(n,n)
    real (kind=real_kind)              :: flux_r(n,n)
    real (kind=real_kind)              :: flux_b(n,n)
    real (kind=real_kind)              :: flux_t(n,n)
    
    integer i,j

    if (.not.ALLOCATED(integration_matrix)      .or. &
        SIZE(integration_matrix,1).ne.n .or. &
        SIZE(integration_matrix,2).ne.p) then
      call abortmp( 'FATAL ERROR: allocate_subcell_integration_matrix not called')
    end if

    v(:,:,1) = u(:,:,1)*metdet(:,:)
    v(:,:,2) = u(:,:,2)*metdet(:,:)

    tb = MATMUL(integration_matrix, v(:,:,2))
    flux_b(:,:) = MATMUL(tb,TRANSPOSE(boundary_interp_matrix(:,1,:)))
    flux_t(:,:) = MATMUL(tb,TRANSPOSE(boundary_interp_matrix(:,2,:)))

    lr = MATMUL(v(:,:,1),TRANSPOSE(integration_matrix))
    flux_l(:,:) = MATMUL(boundary_interp_matrix(:,1,:),lr)
    flux_r(:,:) = MATMUL(boundary_interp_matrix(:,2,:),lr)

    fluxes(:,:,1) = -flux_b(:,:)*rrearth
    fluxes(:,:,2) =  flux_r(:,:)*rrearth
    fluxes(:,:,3) =  flux_t(:,:)*rrearth
    fluxes(:,:,4) = -flux_l(:,:)*rrearth

  end function subcell_div_fluxes

  function subcell_Laplace_fluxes(u, deriv, elem, p, n) result(fluxes)

    implicit none

    integer              , intent(in)  :: p
    integer              , intent(in)  :: n
    type (derivative_t)  , intent(in)  :: deriv
    type (element_t)     , intent(in)  :: elem
    real (kind=real_kind), intent(in)  :: u(p,p)

    real (kind=real_kind)              :: g(p,p,2)
    real (kind=real_kind)              :: v(p,p,2)
    real (kind=real_kind)              :: div(p,p,2)
    real (kind=real_kind)              :: sub_int(n,n,2)
    real (kind=real_kind)              :: fluxes(n,n,4)
    
    integer i,j

    g = gradient_sphere(u,deriv,elem%Dinv)

    v(:,:,1) = elem%Dinv(:,:,1,1)*g(:,:,1) + elem%Dinv(:,:,1,2)*g(:,:,2)
    v(:,:,2) = elem%Dinv(:,:,2,1)*g(:,:,1) + elem%Dinv(:,:,2,2)*g(:,:,2)
    do j=1,p
    do i=1,p
       div(i,j,1) = -SUM(elem%spheremp(:,j)*v(:,j,1)*deriv%Dvv(i,:))
       div(i,j,2) = -SUM(elem%spheremp(i,:)*v(i,:,2)*deriv%Dvv(j,:))
    end do
    end do
    div = div * rrearth

    div(:,:,1) = div(:,:,1) / elem%spheremp(:,:)
    div(:,:,2) = div(:,:,2) / elem%spheremp(:,:)

    sub_int(:,:,1)  = subcell_integration(div(:,:,1), p, n, elem%metdet)
    sub_int(:,:,2)  = subcell_integration(div(:,:,2), p, n, elem%metdet)

    do i=1,n
    do j=2,n
      sub_int(j,i,1) = sub_int(j,i,1) + sub_int(j-1,i,1) 
      sub_int(i,j,2) = sub_int(i,j,2) + sub_int(i,j-1,2) 
    end do
    end do

    fluxes = 0
    do i=1,n
    do j=1,n
      if (1.lt.j) fluxes(i,j,1) = -sub_int(i,  j-1,2)
      if (i.lt.n) fluxes(i,j,2) =  sub_int(i,  j,  1)
      if (j.lt.n) fluxes(i,j,3) =  sub_int(i,  j,  2)
      if (1.lt.i) fluxes(i,j,4) = -sub_int(i-1,j,  1)
    end do
    end do

  end function subcell_Laplace_fluxes


  ! Given a field defined on the unit element, [-1,1]x[-1,1]
  ! sample values, sampled_val, premultiplied by integration weights,
  ! and a number, np, of Gauss-Lobatto-Legendre points. Divide
  ! the square up into intervals by intervals sub-squares so that
  ! there are now intervals**2 sub-cells.  Integrate the 
  ! function defined by sampled_val over each of these
  ! sub-cells and return the integrated values as an 
  ! intervals by intervals matrix.
  !
  ! Efficiency is obtained by computing and caching the appropriate
  ! integration matrix the first time the function is called.
  function subcell_integration(sampled_val, np, intervals, metdet) result(values)

    implicit none

    integer              , intent(in)  :: np
    integer              , intent(in)  :: intervals
    real (kind=real_kind), intent(in)  :: sampled_val(np,np)
    real (kind=real_kind), intent(in)  :: metdet(np,np)
    real (kind=real_kind)              :: val(np,np)
    real (kind=real_kind)              :: values(intervals,intervals)

    integer i,j

    if (.not.ALLOCATED(integration_matrix)      .or. &
        SIZE(integration_matrix,1).ne.intervals .or. &
        SIZE(integration_matrix,2).ne.np) then
      call abortmp( 'FATAL ERROR: allocate_subcell_integration_matrix not called')
    end if

    ! Multiply sampled values by spectral element weights
    val = sampled_val * metdet 

    ! Multiply the sampled values by the weighted jacobians.  
    ! Symmetry allows us to write this as J^t V J
    ! where J is a vector.  

    values = MATMUL(integration_matrix, &
             MATMUL(val,TRANSPOSE(integration_matrix)))

  end function subcell_integration


  ! Helper subroutine that will fill in a matrix needed to 
  ! integrate a function defined on the GLL points of a unit
  ! square on sub-cells.  So np is the number of integration
  ! GLL points defined on the unit square (actually [-1,1]x[-1,1])
  ! and intervals is the number to cut it up into, say a 3 by 3
  ! set of uniform sub-cells.  This function will fill the 
  ! subcell_integration matrix with the correct coefficients
  ! to integrate over each subcell.  
  subroutine allocate_subcell_integration_matrix(np, intervals)
    !-----------------
    !-----------------
    use quadrature_mod, only : gausslobatto, quadrature_t
    
    implicit none

    integer              , intent(in)  :: np
    integer              , intent(in)  :: intervals
    real (kind=real_kind)              :: values(intervals,intervals)


    real(kind=real_kind), parameter :: zero = 0.0D0, one=1.0D0, two=2.0D0


    real (kind=real_kind) :: sub_gll        (intervals,np)

    real (kind=real_kind) :: Lagrange_interp(intervals,np,np)
    type (quadrature_t)   :: gll 

    real (kind=real_kind) :: legrange_div(np)
    real (kind=real_kind) :: a,b,x,y, x_j, x_i 
    real (kind=real_kind) :: r(1) 
    integer i,j,n,m

    if (ALLOCATED(integration_matrix)) deallocate(integration_matrix)
    allocate(integration_matrix(intervals,np))
    if (ALLOCATED(boundary_interp_matrix)) deallocate(boundary_interp_matrix)
    allocate(boundary_interp_matrix(intervals,2,np))

    gll = gausslobatto(np)
 
    ! The GLL (Gauss-Lobatto-Legendre) points are from [-1,1], 
    ! we have a bunch of sub-intervals defined by intervals that 
    ! go from [a,b] so we need to linearly map [-1,1] -> [a,b] 
    ! all the  GLL points by  y = (a/2)(1-x) + (b/2)(1+x)
    do i=1,intervals
      a = -one + (i-one)*two/intervals   
      b = -one +  i     *two/intervals  
      sub_gll(i,:) = (a+b)/two + gll%points(:)/intervals
    end do

    ! Now to interpolate from the values at the input GLL
    ! points to the sub-GLL points.  Do this by Lagrange
    ! interpolation.  The jth Lagrange interpolating polynomial
    ! for points x_i is 
    !              \prod_{i\ne j} (x-x_i)/(x_j-x_i)
    ! These are then multiplied by the sampled values y_i 
    ! and summed. 
    
    ! Save some time by pre-computing the denominitor. I think 
    ! this is OK since all the points are of order 1 so should
    ! be well behaved.
    do n = 1,np
      x_j = gll%points(n)
      x   = one 
      do m = 1,np 
        if (m.ne.n) then
          x_i = gll%points(m)
          x = x * (x_j-x_i)
        endif
      end do
      legrange_div(n)= x
    end do 
    do i=1,intervals
      do n=1,np
        x = sub_gll(i,n)
        do j = 1,np
          y = one
          do m = 1,np
            if (m.ne.j) then
              x_i = gll%points(m)
              y = y * (x-x_i)
            end if
          end do
          Lagrange_interp(i,n,j) = y/legrange_div(j)
        end do
      end do
    end do

    ! Integration is the GLL weights times Jacobians times
    ! the interpolated values:
    !                   w^t I Y I^t w 
    ! where  
    ! w is GLL weights and Jacobians, 
    ! I is the Lagrange_interp matrix, and
    ! Y is the coefficient matrix, sampled_val.
    ! This can be written  J Y J^t where
    !                       J = w^t I
    ! J is integration_matrix
    do i=1,intervals
      integration_matrix(i,:) = MATMUL(gll%weights(:),Lagrange_interp(i,:,:))
    end do

    ! There is still the Jacobian to consider.  We are 
    ! integrating over [a,b] x [c,d] where 
    !        |b-a| = |d-c| = 2/Intervals
    ! Multiply the weights appropriately given that 
    ! they are defined for a 2x2 square
    integration_matrix = integration_matrix/intervals

    boundary_interp_matrix(:,:,:) = Lagrange_interp(:,(/1,np/),:)
  end subroutine allocate_subcell_integration_matrix



  subroutine limiter_optim_iter_full(ptens,sphweights,minp,maxp,dpmass)
    ! 
    !The idea here is the following: We need to find a grid field which is closest
    !to the initial field (in terms of weighted sum), but satisfies the min/max constraints.
    !So, first we find values which do not satisfy constraints and bring these values
    !to a closest constraint. This way we introduce some mass change (addmass),
    !so, we redistribute addmass in the way that l2 error is smallest.
    !This redistribution might violate constraints thus, we do a few iterations.
    !
    ! O. Guba ~2012                    Documented in Guba, Taylor & St-Cyr, JCP 2014
    ! I. Demeshko & M. Taylor 7/2015:  Removed indirect addressing.  
    ! N. Lopez & M. Taylor 8/2015:     Mass redistributon tweak which is better at 
    !                                  linear coorelation preservation
    !
    use kinds         , only : real_kind
    use dimensions_mod, only : np, np, nlev

    real (kind=real_kind), dimension(nlev), intent(inout)   :: minp, maxp
    real (kind=real_kind), dimension(np,np,nlev), intent(inout)   :: ptens
    real (kind=real_kind), dimension(np,np,nlev), intent(in), optional  :: dpmass
    real (kind=real_kind), dimension(np,np), intent(in)   :: sphweights

    real (kind=real_kind), dimension(np,np) :: ptens_mass
    integer  k1, k, i, j, iter, weightsnum
    real (kind=real_kind) :: addmass, weightssum, mass, sumc
    real (kind=real_kind) :: x(np*np),c(np*np)
    integer :: maxiter = np*np-1
    real (kind=real_kind) :: tol_limiter = 5e-14

 
    do k=1,nlev

     k1=1
     do i=1,np
     do j=1,np
       c(k1)=sphweights(i,j)*dpmass(i,j,k)
       x(k1)=ptens(i,j,k)/dpmass(i,j,k)
       k1=k1+1
      enddo
     enddo

     sumc=sum(c)
     if (sumc <= 0 ) CYCLE   ! this should never happen, but if it does, dont limit
     mass=sum(c*x)

    

      ! relax constraints to ensure limiter has a solution:
      ! This is only needed if runnign with the SSP CFL>1 or
      ! due to roundoff errors
      if( mass < minp(k)*sumc ) then
        minp(k) = mass / sumc
      endif
      if( mass > maxp(k)*sumc ) then
        maxp(k) = mass / sumc
      endif

    

     do iter=1,maxiter

      addmass=0.0d0

       do k1=1,np*np
         if((x(k1)>maxp(k))) then
           addmass=addmass+(x(k1)-maxp(k))*c(k1)
           x(k1)=maxp(k)
         endif
         if((x(k1)<minp(k))) then
           addmass=addmass-(minp(k)-x(k1))*c(k1)
           x(k1)=minp(k)
         endif
       enddo !k1

       if(abs(addmass)<=tol_limiter*abs(mass)) exit

       weightssum=0.0d0
!       weightsnum=0
       if(addmass>0)then
        do k1=1,np*np
          if(x(k1)<maxp(k))then
            weightssum=weightssum+c(k1)
!            weightsnum=weightsnum+1
          endif
        enddo !k1
        do k1=1,np*np
          if(x(k1)<maxp(k))then
              x(k1)=x(k1)+addmass/weightssum
!              x(k1)=x(k1)+addmass/(c(k1)*weightsnum)
          endif
        enddo
      else
        do k1=1,np*np
          if(x(k1)>minp(k))then
            weightssum=weightssum+c(k1)
!            weightsnum=weightsnum+1
          endif
        enddo
        do k1=1,np*np
          if(x(k1)>minp(k))then
            x(k1)=x(k1)+addmass/weightssum
!           x(k1)=x(k1)+addmass/(c(k1)*weightsnum)
          endif
        enddo
      endif


   enddo!end of iteration

   k1=1
   do i=1,np
    do j=1,np
      ptens(i,j,k)=x(k1)
      k1=k1+1
    enddo
   enddo

  enddo

  do k=1,nlev
    ptens(:,:,k)=ptens(:,:,k)*dpmass(:,:,k)
  enddo
 
  end subroutine limiter_optim_iter_full





end module derivative_mod_base
