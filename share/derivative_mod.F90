#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


module derivative_mod
  use kinds, only : real_kind, longdouble_kind
  use dimensions_mod, only : np, nc, npdg, nep
  use quadrature_mod, only : quadrature_t, gauss, gausslobatto,legendre, jacobi
  use parallel_mod, only : abortmp
  ! needed for spherical differential operators:
  use physical_constants, only : rrearth 
  use element_mod, only : element_t
  use control_mod, only : which_vlaplace 

implicit none
private

  type, public :: derivative_t
     real (kind=real_kind) :: Dvv(np,np)
     real (kind=real_kind) :: Dvv_diag(np,np)
     real (kind=real_kind) :: Dvv_twt(np,np)
     real (kind=real_kind) :: Mvv_twt(np,np)  ! diagonal matrix of GLL weights
     real (kind=real_kind) :: vvtemp(np,np)
     real (kind=real_kind) :: vvtempt(np,np,2)
     real (kind=real_kind) :: Mfvm(np,nc+1)
     real (kind=real_kind) :: Cfvm(np,nc)
     real (kind=real_kind) :: Sfvm(np,nep)
     real (kind=real_kind) :: legdg(np,np)
  end type derivative_t

  type, public :: derivative_stag_t
     real (kind=real_kind) :: D(np,np)
     real (kind=real_kind) :: M(np,np)
     real (kind=real_kind) :: Dpv(np,np)
     real (kind=real_kind) :: D_twt(np,np)
     real (kind=real_kind) :: M_twt(np,np)
     real (kind=real_kind) :: M_t(np,np)
     real (kind=real_kind) :: vtemp(np,np,2)
     real (kind=real_kind) :: vtempt(np,np,2)
  end type derivative_stag_t

  real (kind=real_kind), allocatable :: integration_matrix(:,:)
  private :: allocate_subcell_integration_matrix

! ======================================
! Public Interfaces
! ======================================

  public :: subcell_integration

  public :: derivinit
  public :: deriv_print

  public :: gradient
  public :: gradient_wk
  public :: vorticity
  public :: divergence

  public :: interpolate_gll2fvm_corners
  public :: interpolate_gll2fvm_points
  public :: interpolate_gll2spelt_points
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
  public  :: gradient_sphere
  public  :: gradient_sphere_wk
  public  :: ugradv_sphere
  public  :: vorticity_sphere
  public  :: vorticity_sphere_diag
  public  :: divergence_sphere
  public  :: curl_sphere
  public  :: divergence_sphere_wk
  public  :: laplace_sphere_wk
  public  :: vlaplace_sphere_wk
  public  :: element_boundary_integral
  public  :: edge_flux_u_cg
  public  :: gll_to_dgmodal
  public  :: dgmodal_to_gll



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

    ! Legendre polynomials of degree npdg-1, on the np GLL grid:
    if (npdg>np) call abortmp( 'FATAL ERROR: npdg>np')
    if (npdg>0 .and. npdg<np) then
       ! in this case, we use a DG basis of Legendre polynomials
       ! stored at the GLL points.  
       do i=1,np
          deriv%legdg(1:npdg,i) = legendre(gp%points(i),npdg-1)
       end do
       ! normalize
       do j=1,npdg
          xnorm=sqrt(sum(deriv%legdg(j,:)*deriv%legdg(j,:)*gp%weights(:)))
          deriv%legdg(j,:)=deriv%legdg(j,:)/xnorm
       enddo
    endif
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
    type (derivative_t) :: deriv
    
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
    real(kind=real_kind)  ::  v2p(n1,n2)
    real(kind=real_kind)  ::  v2p_new(n1,n2)
    real(kind=longdouble_kind)  ::  gll(n1),gs(n2)
    integer :: n1,n2
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
    type (derivative_stag_t)         :: deriv

    real(kind=real_kind) :: div(np,np)

    ! Local

    integer i
    integer j
    integer l
    logical, parameter :: UseUnroll = .TRUE.

    real(kind=real_kind)  sumx00,sumx01
    real(kind=real_kind)  sumy00,sumy01
    real(kind=real_kind)  sumx10,sumx11
    real(kind=real_kind)  sumy10,sumy11

#ifdef DEBUG
    print *, "divergence_stag"
#endif
if(MODULO(np,2) == 0 .and. UseUnroll) then 
    !JMD====================================
    !JMD  2*np*np*np Flops
    !JMD====================================
    do j=1,np,2
       do l=1,np,2

          sumx00=0.0d0
          sumx01=0.0d0
          sumx10=0.0d0
          sumx11=0.0d0

          sumy00=0.0d0
          sumy01=0.0d0
          sumy10=0.0d0
          sumy11=0.0d0

          do i=1,np
             sumx00 = sumx00 + deriv%D(i,l  )*v(i,j  ,1)
             sumx01 = sumx01 + deriv%D(i,l+1)*v(i,j  ,1)
             sumx10 = sumx10 + deriv%D(i,l  )*v(i,j+1,1)
             sumx11 = sumx11 + deriv%D(i,l+1)*v(i,j+1,1)

             sumy00 = sumy00 + deriv%M(i,l  )*v(i,j  ,2)
             sumy01 = sumy01 + deriv%M(i,l+1)*v(i,j  ,2)
             sumy10 = sumy10 + deriv%M(i,l  )*v(i,j+1,2)
             sumy11 = sumy11 + deriv%M(i,l+1)*v(i,j+1,2)
          end do

          deriv%vtemp(j  ,l  ,1) = sumx00
          deriv%vtemp(j  ,l+1,1) = sumx01
          deriv%vtemp(j+1,l  ,1) = sumx10
          deriv%vtemp(j+1,l+1,1) = sumx11

          deriv%vtemp(j  ,l  ,2) = sumy00
          deriv%vtemp(j  ,l+1,2) = sumy01
          deriv%vtemp(j+1,l  ,2) = sumy10
          deriv%vtemp(j+1,l+1,2) = sumy11

       end do
    end do


    !JMD====================================
    !JMD  2*np*np*np Flops
    !JMD====================================
    do j=1,np,2
       do i=1,np,2
          sumx00=0.0d0
          sumx01=0.0d0
          sumx10=0.0d0
          sumx11=0.0d0

          sumy00=0.0d0
          sumy01=0.0d0
          sumy10=0.0d0
          sumy11=0.0d0

          do l=1,np
             sumx00 = sumx00 +  deriv%M(l,j  )*deriv%vtemp(l,i  ,1)
             sumx01 = sumx01 +  deriv%M(l,j+1)*deriv%vtemp(l,i  ,1)
             sumx10 = sumx10 +  deriv%M(l,j  )*deriv%vtemp(l,i+1,1)
             sumx11 = sumx11 +  deriv%M(l,j+1)*deriv%vtemp(l,i+1,1)

             sumy00 = sumy00 +  deriv%D(l,j  )*deriv%vtemp(l,i  ,2)
             sumy01 = sumy01 +  deriv%D(l,j+1)*deriv%vtemp(l,i  ,2)
             sumy10 = sumy10 +  deriv%D(l,j  )*deriv%vtemp(l,i+1,2)
             sumy11 = sumy11 +  deriv%D(l,j+1)*deriv%vtemp(l,i+1,2)
          end do

          div(i  ,j  ) = sumx00 + sumy00
          div(i  ,j+1) = sumx01 + sumy01
          div(i+1,j  ) = sumx10 + sumy10
          div(i+1,j+1) = sumx11 + sumy11

       end do
    end do
else
     do j=1,np
        do l=1,np
 
           sumx00=0.0d0
           sumy00=0.0d0
           do i=1,np
              sumx00 = sumx00 + deriv%D(i,l  )*v(i,j  ,1)
              sumy00 = sumy00 + deriv%M(i,l  )*v(i,j  ,2)
 	   enddo
          deriv%vtemp(j  ,l  ,1) = sumx00
          deriv%vtemp(j  ,l  ,2) = sumy00
       enddo
    enddo
    do j=1,np
       do i=1,np
          sumx00=0.0d0
	  sumy00=0.0d0
          do l=1,np
             sumx00 = sumx00 +  deriv%M(l,j  )*deriv%vtemp(l,i  ,1)
	     sumy00 = sumy00 +  deriv%D(l,j  )*deriv%vtemp(l,i  ,2)
	  enddo
          div(i  ,j  ) = sumx00 + sumy00

	enddo
    enddo
endif

  end function divergence_stag

!  ================================================
!  divergence_nonstag: 
!
!  Compute divergence (maps v->v)
!  ================================================

  function divergence_nonstag(v,deriv) result(div)

    real(kind=real_kind), intent(in) :: v(np,np,2)
    type (derivative_t)              :: deriv

    real(kind=real_kind) :: div(np,np)

    ! Local

    integer i
    integer j
    integer l

    logical, parameter :: UseUnroll = .TRUE.

    real(kind=real_kind) ::  dudx00,dudx01
    real(kind=real_kind) ::  dudx10,dudx11

    real(kind=real_kind) ::  dvdy00,dvdy01
    real(kind=real_kind) ::  dvdy10,dvdy11

if(modulo(np,2) .eq. 0 .and. UseUnroll) then
! this is just loop unrolling - a good compiler should do it for you jpe
       do j=1,np,2
          do l=1,np,2

             dudx00=0.0d0
             dudx01=0.0d0
             dudx10=0.0d0
             dudx11=0.0d0

             dvdy00=0.0d0
             dvdy01=0.0d0
             dvdy10=0.0d0
             dvdy11=0.0d0

             do i=1,np
                
                dudx00 = dudx00 + deriv%Dvv(i,l  )*v(i,j  ,1)
                dudx01 = dudx01 + deriv%Dvv(i,l+1)*v(i,j  ,1)
                dudx10 = dudx10 + deriv%Dvv(i,l  )*v(i,j+1,1)
                dudx11 = dudx11 + deriv%Dvv(i,l+1)*v(i,j+1,1)
                
                dvdy00 = dvdy00 + deriv%Dvv(i,l  )*v(j  ,i,2)
                dvdy01 = dvdy01 + deriv%Dvv(i,l+1)*v(j  ,i,2)
                dvdy10 = dvdy10 + deriv%Dvv(i,l  )*v(j+1,i,2)
                dvdy11 = dvdy11 + deriv%Dvv(i,l+1)*v(j+1,i,2)

             end do

             div(l  ,j  ) = dudx00
             div(l+1,j  ) = dudx01
             div(l  ,j+1) = dudx10
             div(l+1,j+1) = dudx11

             deriv%vvtemp(j  ,l  ) = dvdy00
             deriv%vvtemp(j  ,l+1) = dvdy01
             deriv%vvtemp(j+1,l  ) = dvdy10
             deriv%vvtemp(j+1,l+1) = dvdy11

          end do
       end do
    else

       do j=1,np
          do l=1,np
             dudx00=0.0d0
             dvdy00=0.0d0

             do i=1,np
                dudx00 = dudx00 + deriv%Dvv(i,l  )*v(i,j  ,1)
                dvdy00 = dvdy00 + deriv%Dvv(i,l  )*v(j  ,i,2)
             end do

             div(l  ,j  ) = dudx00
             deriv%vvtemp(j  ,l  ) = dvdy00


          end do
       end do
    end if
    do j=1,np
       do i=1,np
          div(i,j)=div(i,j)+deriv%vvtemp(i,j)
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

    type (derivative_stag_t)         :: deriv
    real(kind=real_kind), intent(in) :: p(np,np)

    real(kind=real_kind)             :: dp(np,np,2)

    ! Local
      
    integer i
    integer j
    integer l
    logical, parameter :: UseUnroll = .TRUE.

    real(kind=real_kind)  sumx00,sumx01
    real(kind=real_kind)  sumy00,sumy01
    real(kind=real_kind)  sumx10,sumx11
    real(kind=real_kind)  sumy10,sumy11

#ifdef DEBUG
    print *, "gradient_wk_stag"
#endif
    !JMD ================================
    !JMD 2*np*np*np Flops 
    !JMD ================================

if(MODULO(np,2) == 0 .and. UseUnroll) then 


    do j=1,np,2
       do l=1,np,2
          sumx00=0.0d0
          sumx01=0.0d0
          sumx10=0.0d0
          sumx11=0.0d0

          sumy00=0.0d0
          sumy01=0.0d0
          sumy10=0.0d0
          sumy11=0.0d0

          do i=1,np
             sumx00 = sumx00 + deriv%D_twt(i,l  )*p(i,j  )
             sumx01 = sumx01 + deriv%D_twt(i,l+1)*p(i,j  )
             sumx10 = sumx10 + deriv%D_twt(i,l  )*p(i,j+1)
             sumx11 = sumx11 + deriv%D_twt(i,l+1)*p(i,j+1)

             sumy00 = sumy00 + deriv%M_twt(i,l  )*p(i,j  )
             sumy01 = sumy01 + deriv%M_twt(i,l+1)*p(i,j  )
             sumy10 = sumy10 + deriv%M_twt(i,l  )*p(i,j+1)
             sumy11 = sumy11 + deriv%M_twt(i,l+1)*p(i,j+1)
          end do

          deriv%vtempt(j  ,l  ,1) = sumx00
          deriv%vtempt(j  ,l+1,1) = sumx01
          deriv%vtempt(j+1,l  ,1) = sumx10
          deriv%vtempt(j+1,l+1,1) = sumx11

          deriv%vtempt(j  ,l  ,2) = sumy00
          deriv%vtempt(j  ,l+1,2) = sumy01
          deriv%vtempt(j+1,l  ,2) = sumy10
          deriv%vtempt(j+1,l+1,2) = sumy11
       end do
    end do

	
    !JMD ================================
    !JMD 2*np*np*np Flops 
    !JMD ================================

    do j=1,np,2
       do i=1,np,2
          sumx00=0.0d0
          sumx01=0.0d0
          sumx10=0.0d0
          sumx11=0.0d0

          sumy00=0.0d0
          sumy01=0.0d0
          sumy10=0.0d0
          sumy11=0.0d0

          do l=1,np
             sumx00 = sumx00 +  deriv%M_twt(l,j  )*deriv%vtempt(l,i  ,1)
             sumx01 = sumx01 +  deriv%M_twt(l,j+1)*deriv%vtempt(l,i  ,1)
             sumx10 = sumx10 +  deriv%M_twt(l,j  )*deriv%vtempt(l,i+1,1)
             sumx11 = sumx11 +  deriv%M_twt(l,j+1)*deriv%vtempt(l,i+1,1)

             sumy00 = sumy00 +  deriv%D_twt(l,j  )*deriv%vtempt(l,i  ,2)
             sumy01 = sumy01 +  deriv%D_twt(l,j+1)*deriv%vtempt(l,i  ,2)
             sumy10 = sumy10 +  deriv%D_twt(l,j  )*deriv%vtempt(l,i+1,2)
             sumy11 = sumy11 +  deriv%D_twt(l,j+1)*deriv%vtempt(l,i+1,2)
          end do

          dp(i  ,j  ,1) = sumx00
          dp(i  ,j+1,1) = sumx01
          dp(i+1,j  ,1) = sumx10
          dp(i+1,j+1,1) = sumx11

          dp(i  ,j  ,2) = sumy00
          dp(i  ,j+1,2) = sumy01
          dp(i+1,j  ,2) = sumy10
          dp(i+1,j+1,2) = sumy11

       end do
    end do
else
    do j=1,np
       do l=1,np
 	  sumx00=0.0d0
          sumy00=0.0d0
           do i=1,np
              sumx00 = sumx00 + deriv%D_twt(i,l  )*p(i,j  )
              sumy00 = sumy00 + deriv%M_twt(i,l  )*p(i,j  )
 	  enddo
           deriv%vtempt(j  ,l  ,1) = sumx00
           deriv%vtempt(j  ,l  ,2) = sumy00
        enddo
    enddo
    do j=1,np
       do i=1,np
          sumx00=0.0d0
	  sumy00=0.0d0
          do l=1,np
             sumx00 = sumx00 +  deriv%M_twt(l,j  )*deriv%vtempt(l,i  ,1)
	     sumy00 = sumy00 +  deriv%D_twt(l,j  )*deriv%vtempt(l,i  ,2)
	  enddo
	  dp(i  ,j  ,1) = sumx00
          dp(i  ,j  ,2) = sumy00
      enddo
    enddo
endif


  end function gradient_wk_stag

!  ================================================
!  gradient_wk_nonstag:
! 
!  Compute the weak form gradient:
!  maps scalar field on the Gauss-Lobatto grid to the
!  weak gradient on the Gauss-Lobbatto grid
!  ================================================

  function gradient_wk_nonstag(p,deriv) result(dp)

    type (derivative_t)         :: deriv
    real(kind=real_kind), intent(in) :: p(np,np)

    real(kind=real_kind)             :: dp(np,np,2)

    ! Local
      
    integer i
    integer j
    integer l
    logical, parameter :: UseUnroll = .TRUE.

    real(kind=real_kind)  sumx00,sumx01
    real(kind=real_kind)  sumy00,sumy01
    real(kind=real_kind)  sumx10,sumx11
    real(kind=real_kind)  sumy10,sumy11

    !JMD ================================
    !JMD 2*np*np*np Flops 
    !JMD ================================

!   print *, "gradient_wk_nonstag"
    if(modulo(np,2) .eq. 0 .and. UseUnroll) then
! this is just loop unrolling - a good compiler should do it for you jpe

       do j=1,np,2
          do l=1,np,2
             sumx00=0.0d0
             sumx01=0.0d0
             sumx10=0.0d0
             sumx11=0.0d0

             sumy00=0.0d0
             sumy01=0.0d0
             sumy10=0.0d0
             sumy11=0.0d0

             do i=1,np
                sumx00 = sumx00 + deriv%Dvv_twt(i,l  )*p(i,j  )
                sumx01 = sumx01 + deriv%Dvv_twt(i,l+1)*p(i,j  )
                sumx10 = sumx10 + deriv%Dvv_twt(i,l  )*p(i,j+1)
                sumx11 = sumx11 + deriv%Dvv_twt(i,l+1)*p(i,j+1)

                sumy00 = sumy00 + deriv%Mvv_twt(i,l  )*p(i,j  )
                sumy01 = sumy01 + deriv%Mvv_twt(i,l+1)*p(i,j  )
                sumy10 = sumy10 + deriv%Mvv_twt(i,l  )*p(i,j+1)
                sumy11 = sumy11 + deriv%Mvv_twt(i,l+1)*p(i,j+1)
             end do

             deriv%vvtempt(j  ,l  ,1) = sumx00
             deriv%vvtempt(j  ,l+1,1) = sumx01
             deriv%vvtempt(j+1,l  ,1) = sumx10
             deriv%vvtempt(j+1,l+1,1) = sumx11

             deriv%vvtempt(j  ,l  ,2) = sumy00
             deriv%vvtempt(j  ,l+1,2) = sumy01
             deriv%vvtempt(j+1,l  ,2) = sumy10
             deriv%vvtempt(j+1,l+1,2) = sumy11

          end do
       end do
       ! vvtempt1 = p'*Dvv_twt
       ! vvtempt2 = p'*Mvv_twt
       ! dp1 = dy*Mvv_twt*vvtempt1' = dy*Mvv_twt*(p'*Dvv_twt)' = dy*Mvv_twt*Dvv_twt'*p
       ! dp2 = dx*Dvv_twt*vvtempt2' = dx*Dvv_twt*(p'*Mvv_twt)' = dx*Dvv_twt*Mvv_twt'*p
       !     New formulation 
       ! dp1 = dy*MvvDvvt*p
       ! dp2 = dx*DvvMvvt*p
       ! MvvDvvt = Mvv_twt*Dvv_twt'
       ! DvvMvvt = Dvv_twt*Mvv_twt'


       !JMD ================================
       !JMD 2*np*np*np Flops 
       !JMD ================================

       do j=1,np,2
          do i=1,np,2
             sumx00=0.0d0
             sumx01=0.0d0
             sumx10=0.0d0
             sumx11=0.0d0
             
             sumy00=0.0d0
             sumy01=0.0d0
             sumy10=0.0d0
             sumy11=0.0d0
             
             do l=1,np
                sumx00 = sumx00 +  deriv%Mvv_twt(l,j  )*deriv%vvtempt(l,i  ,1)
                sumx01 = sumx01 +  deriv%Mvv_twt(l,j+1)*deriv%vvtempt(l,i  ,1)
                sumx10 = sumx10 +  deriv%Mvv_twt(l,j  )*deriv%vvtempt(l,i+1,1)
                sumx11 = sumx11 +  deriv%Mvv_twt(l,j+1)*deriv%vvtempt(l,i+1,1)

                sumy00 = sumy00 +  deriv%Dvv_twt(l,j  )*deriv%vvtempt(l,i  ,2)
                sumy01 = sumy01 +  deriv%Dvv_twt(l,j+1)*deriv%vvtempt(l,i  ,2)
                sumy10 = sumy10 +  deriv%Dvv_twt(l,j  )*deriv%vvtempt(l,i+1,2)
                sumy11 = sumy11 +  deriv%Dvv_twt(l,j+1)*deriv%vvtempt(l,i+1,2)
             end do
             
             dp(i  ,j  ,1) = sumx00
             dp(i  ,j+1,1) = sumx01
             dp(i+1,j  ,1) = sumx10
             dp(i+1,j+1,1) = sumx11
             
             dp(i  ,j  ,2) = sumy00
             dp(i  ,j+1,2) = sumy01
             dp(i+1,j  ,2) = sumy10
             dp(i+1,j+1,2) = sumy11

          end do
       end do
    else

       do j=1,np
          do l=1,np
             sumx00=0.0d0

             sumy00=0.0d0

             do i=1,np
                sumx00 = sumx00 + deriv%Dvv_twt(i,l  )*p(i,j  )

                sumy00 = sumy00 + deriv%Mvv_twt(i,l  )*p(i,j  )
             end do

             deriv%vvtempt(j  ,l  ,1) = sumx00

             deriv%vvtempt(j  ,l  ,2) = sumy00

          end do
       end do

       !JMD ================================
       !JMD 2*np*np*np Flops 
       !JMD ================================

       do j=1,np
          do i=1,np
             sumx00=0.0d0
             
             sumy00=0.0d0
             
             do l=1,np
                sumx00 = sumx00 +  deriv%Mvv_twt(l,j  )*deriv%vvtempt(l,i  ,1)

                sumy00 = sumy00 +  deriv%Dvv_twt(l,j  )*deriv%vvtempt(l,i  ,2)
             end do
             
             dp(i  ,j  ,1) = sumx00
             
             dp(i  ,j  ,2) = sumy00

          end do
       end do
    end if
  end function gradient_wk_nonstag

!  ================================================
!  gradient_str_stag:
! 
!  Compute the *strong* form gradient:
!  maps scalar field on the pressure grid to the
!  velocity grid
!  ================================================

  function gradient_str_stag(p,deriv) result(dp)

    type (derivative_stag_t)         :: deriv
    real(kind=real_kind), intent(in) :: p(np,np)

    real(kind=real_kind)             :: dp(np,np,2)

    ! Local
      
    integer i
    integer j
    integer l

    logical, parameter :: UseUnroll=.TRUE.

    real(kind=real_kind)  sumx00,sumx01
    real(kind=real_kind)  sumy00,sumy01
    real(kind=real_kind)  sumx10,sumx11
    real(kind=real_kind)  sumy10,sumy11

#ifdef DEBUG
    print *, "gradient_str_stag"
#endif
if(MODULO(np,2) == 0 .and. UseUnroll) then 
    do j=1,np,2
       do l=1,np,2
          sumx00=0.0d0
          sumx01=0.0d0
          sumx10=0.0d0
          sumx11=0.0d0

          sumy00=0.0d0
          sumy01=0.0d0
          sumy10=0.0d0
          sumy11=0.0d0

          do i=1,np
             sumx00 = sumx00 + deriv%Dpv(i,l  )*p(i,j  )
             sumx01 = sumx01 + deriv%Dpv(i,l+1)*p(i,j  )
             sumx10 = sumx10 + deriv%Dpv(i,l  )*p(i,j+1)
             sumx11 = sumx11 + deriv%Dpv(i,l+1)*p(i,j+1)

             sumy00 = sumy00 + deriv%M_t(i,l  )*p(i,j  )
             sumy01 = sumy01 + deriv%M_t(i,l+1)*p(i,j  )
             sumy10 = sumy10 + deriv%M_t(i,l  )*p(i,j+1)
             sumy11 = sumy11 + deriv%M_t(i,l+1)*p(i,j+1)
          end do

          deriv%vtempt(j  ,l  ,1) = sumx00
          deriv%vtempt(j  ,l+1,1) = sumx01
          deriv%vtempt(j+1,l  ,1) = sumx10
          deriv%vtempt(j+1,l+1,1) = sumx11

          deriv%vtempt(j  ,l  ,2) = sumy00
          deriv%vtempt(j  ,l+1,2) = sumy01
          deriv%vtempt(j+1,l  ,2) = sumy10
          deriv%vtempt(j+1,l+1,2) = sumy11

       end do
    end do




    do j=1,np,2
       do i=1,np,2
          sumx00=0.0d0
          sumx01=0.0d0
          sumx10=0.0d0
          sumx11=0.0d0

          sumy00=0.0d0
          sumy01=0.0d0
          sumy10=0.0d0
          sumy11=0.0d0

          do l=1,np
             sumx00 = sumx00 +  deriv%M_t(l,j  )*deriv%vtempt(l,i  ,1)
             sumx01 = sumx01 +  deriv%M_t(l,j+1)*deriv%vtempt(l,i  ,1)
             sumx10 = sumx10 +  deriv%M_t(l,j  )*deriv%vtempt(l,i+1,1)
             sumx11 = sumx11 +  deriv%M_t(l,j+1)*deriv%vtempt(l,i+1,1)

             sumy00 = sumy00 +  deriv%Dpv(l,j  )*deriv%vtempt(l,i  ,2)
             sumy01 = sumy01 +  deriv%Dpv(l,j+1)*deriv%vtempt(l,i  ,2)
             sumy10 = sumy10 +  deriv%Dpv(l,j  )*deriv%vtempt(l,i+1,2)
             sumy11 = sumy11 +  deriv%Dpv(l,j+1)*deriv%vtempt(l,i+1,2)
          end do

          dp(i  ,j  ,1) = sumx00
          dp(i  ,j+1,1) = sumx01
          dp(i+1,j  ,1) = sumx10
          dp(i+1,j+1,1) = sumx11

          dp(i  ,j  ,2) = sumy00
          dp(i  ,j+1,2) = sumy01
          dp(i+1,j  ,2) = sumy10
          dp(i+1,j+1,2) = sumy11

       end do
    end do
else
    do j=1,np
       do l=1,np
          sumx00=0.0d0
          sumy00=0.0d0
          do i=1,np
             sumx00 = sumx00 + deriv%Dpv(i,l  )*p(i,j  )
             sumy00 = sumy00 + deriv%M_t(i,l  )*p(i,j  )
   	   enddo
	   deriv%vtempt(j  ,l  ,1) = sumx00
   	   deriv%vtempt(j  ,l  ,2) = sumy00
	enddo
    enddo
    do j=1,np
       do i=1,np
          sumx00=0.0d0
          sumy00=0.0d0
          do l=1,np
             sumx00 = sumx00 +  deriv%M_t(l,j  )*deriv%vtempt(l,i  ,1)
             sumy00 = sumy00 +  deriv%Dpv(l,j  )*deriv%vtempt(l,i  ,2)
	  enddo
          dp(i  ,j  ,1) = sumx00
	  dp(i  ,j  ,2) = sumy00
       enddo
    enddo
endif

  end function gradient_str_stag

!  ================================================
!  gradient_str_nonstag:
!
!  Compute the *strong* gradient on the velocity grid
!  of a scalar field on the velocity grid
!  ================================================

  function gradient_str_nonstag(s,deriv) result(ds)

    type (derivative_t)              :: deriv
    real(kind=real_kind), intent(in) :: s(np,np)

    real(kind=real_kind) :: ds(np,np,2)

    integer i
    integer j
    integer l
    logical, parameter :: UseUnroll = .TRUE.

    real(kind=real_kind) ::  dsdx00,dsdx01
    real(kind=real_kind) ::  dsdx10,dsdx11

    real(kind=real_kind) ::  dsdy00,dsdy01
    real(kind=real_kind) ::  dsdy10,dsdy11
#ifdef DEBUG
    print *, "gradient_str_nonstag"
!   write(17) np,s,deriv
#endif
    if(modulo(np,2) .eq. 0 .and. UseUnroll) then
       do j=1,np,2
          do l=1,np,2
             dsdx00=0.0d0
             dsdx01=0.0d0
             dsdx10=0.0d0
             dsdx11=0.0d0

             dsdy00=0.0d0
             dsdy01=0.0d0
             dsdy10=0.0d0
             dsdy11=0.0d0

             do i=1,np
                dsdx00 = dsdx00 + deriv%Dvv(i,l  )*s(i,j  )
                dsdx01 = dsdx01 + deriv%Dvv(i,l+1)*s(i,j  )
                dsdx10 = dsdx10 + deriv%Dvv(i,l  )*s(i,j+1)
                dsdx11 = dsdx11 + deriv%Dvv(i,l+1)*s(i,j+1)

                dsdy00 = dsdy00 + deriv%Dvv(i,l  )*s(j  ,i)
                dsdy01 = dsdy01 + deriv%Dvv(i,l+1)*s(j  ,i)
                dsdy10 = dsdy10 + deriv%Dvv(i,l  )*s(j+1,i)
                dsdy11 = dsdy11 + deriv%Dvv(i,l+1)*s(j+1,i)
             end do
#ifdef DEBUG
             if(j.eq.3.and.l.eq.1) then
                print *, dsdx00
             endif
#endif
             ds(l  ,j  ,1) = dsdx00
             ds(l+1,j  ,1) = dsdx01
             ds(l  ,j+1,1) = dsdx10
             ds(l+1,j+1,1) = dsdx11

             ds(j  ,l  ,2) = dsdy00
             ds(j  ,l+1,2) = dsdy01
             ds(j+1,l  ,2) = dsdy10
             ds(j+1,l+1,2) = dsdy11

          end do

       end do
    else
       do j=1,np
          do l=1,np
             dsdx00=0.0d0

             dsdy00=0.0d0

             do i=1,np
                dsdx00 = dsdx00 + deriv%Dvv(i,l  )*s(i,j  )

                dsdy00 = dsdy00 + deriv%Dvv(i,l  )*s(j  ,i)
             end do
             ds(l  ,j  ,1) = dsdx00

             ds(j  ,l  ,2) = dsdy00

          end do

       end do
    end if
  end function gradient_str_nonstag

!  ================================================
!  vorticity:
!
!  Compute the vorticity of the velocity field on the
!  velocity grid
!  ================================================

  function vorticity(v,deriv) result(vort)

    type (derivative_t)              :: deriv
    real(kind=real_kind), intent(in) :: v(np,np,2)

    real(kind=real_kind) :: vort(np,np)

    integer i
    integer j
    integer l
    
    logical, parameter :: UseUnroll = .TRUE.

    real(kind=real_kind) ::  dvdx00,dvdx01
    real(kind=real_kind) ::  dvdx10,dvdx11

    real(kind=real_kind) ::  dudy00,dudy01
    real(kind=real_kind) ::  dudy10,dudy11

if(MODULO(np,2) == 0 .and. UseUnroll) then 
    do j=1,np,2
       do l=1,np,2

          dudy00=0.0d0
          dudy01=0.0d0
          dudy10=0.0d0
          dudy11=0.0d0

          dvdx00=0.0d0
          dvdx01=0.0d0
          dvdx10=0.0d0
          dvdx11=0.0d0

          do i=1,np

             dvdx00 = dvdx00 + deriv%Dvv(i,l  )*v(i,j  ,2)
             dvdx01 = dvdx01 + deriv%Dvv(i,l+1)*v(i,j  ,2)
             dvdx10 = dvdx10 + deriv%Dvv(i,l  )*v(i,j+1,2)
             dvdx11 = dvdx11 + deriv%Dvv(i,l+1)*v(i,j+1,2)

             dudy00 = dudy00 + deriv%Dvv(i,l  )*v(j  ,i,1)
             dudy01 = dudy01 + deriv%Dvv(i,l+1)*v(j  ,i,1)
             dudy10 = dudy10 + deriv%Dvv(i,l  )*v(j+1,i,1)
             dudy11 = dudy11 + deriv%Dvv(i,l+1)*v(j+1,i,1)

          end do

          vort(l  ,j  ) = dvdx00
          vort(l+1,j  ) = dvdx01
          vort(l  ,j+1) = dvdx10
          vort(l+1,j+1) = dvdx11

          deriv%vvtemp(j  ,l  ) = dudy00
          deriv%vvtemp(j  ,l+1) = dudy01
          deriv%vvtemp(j+1,l  ) = dudy10
          deriv%vvtemp(j+1,l+1) = dudy11

        end do
    end do
else
    do j=1,np
       do l=1,np

          dudy00=0.0d0
	  dvdx00=0.0d0

          do i=1,np
             dvdx00 = dvdx00 + deriv%Dvv(i,l  )*v(i,j  ,2)
             dudy00 = dudy00 + deriv%Dvv(i,l  )*v(j  ,i,1)
	  enddo
 
	  vort(l  ,j  ) = dvdx00
	  deriv%vvtemp(j  ,l  ) = dudy00
	enddo
     enddo

endif

    do j=1,np
       do i=1,np
          vort(i,j)=vort(i,j)-deriv%vvtemp(i,j)
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
          do i=1,np
             sumx00 = sumx00 + deriv%Cfvm(i,l  )*v(i,j  )
          enddo
          vtemp(j  ,l) = sumx00
        enddo
    enddo
    do j=1,nc
       do i=1,nc
          sumx00=0.0d0
          do l=1,np
             sumx00 = sumx00 + deriv%Cfvm(l,j  )*vtemp(l,i)
          enddo
          p(i  ,j  ) = sumx00
       enddo
    enddo
  end function interpolate_gll2fvm_points
  !  ================================================
  !  interpolate_gll2spelt_points:
  !
  !  shape function interpolation from data on GLL grid the spelt grid
  !  Author: Christoph Erath
  !  ================================================
  function interpolate_gll2spelt_points(v,deriv) result(p)
    real(kind=real_kind), intent(in) :: v(np,np)
    type (derivative_t)         :: deriv
    real(kind=real_kind) :: p(nep,nep)

    ! Local
    integer i,j,l

    real(kind=real_kind)  sumx00,sumx01
    real(kind=real_kind)  sumx10,sumx11
    real(kind=real_kind)  vtemp(np,nep)

    do j=1,np
       do l=1,nep
          sumx00=0.0d0
          do i=1,np
             sumx00 = sumx00 + deriv%Sfvm(i,l  )*v(i,j  )
          enddo
          vtemp(j  ,l) = sumx00
        enddo
    enddo
    do j=1,nep
       do i=1,nep
          sumx00=0.0d0
          do l=1,np
             sumx00 = sumx00 + deriv%Sfvm(l,j  )*vtemp(l,i)
          enddo
          p(i  ,j  ) = sumx00
       enddo
    enddo
  end function interpolate_gll2spelt_points

!  ================================================
!  interpolate_gll2fvm_corners:
!
!  shape funtion interpolation from data on GLL grid to physics grid
!
!  ================================================
  function interpolate_gll2fvm_corners(v,deriv) result(p)

    real(kind=real_kind), intent(in) :: v(np,np)
    type (derivative_t)         :: deriv
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
          do i=1,np
             sumx00 = sumx00 + deriv%Mfvm(i,l  )*v(i,j  )
          enddo
          vtemp(j  ,l) = sumx00
        enddo
    enddo
    do j=1,nc+1
       do i=1,nc+1
          sumx00=0.0d0
          do l=1,np
             sumx00 = sumx00 + deriv%Mfvm(l,j  )*vtemp(l,i)
          enddo
          p(i  ,j  ) = sumx00
       enddo
    enddo
  end function interpolate_gll2fvm_corners



!  ================================================
!  remap_phys2gll:
!
!  interpolate to an equally spaced (in reference element coordinate system)
!  "physics" grid to the GLL grid
!
!  1st order, monotone, conservative
!  ================================================
  function remap_phys2gll(pin,nphys) result(pout)
    real(kind=real_kind), intent(in) :: pin(nphys*nphys)
    integer :: nphys
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
    real(kind=longdouble_kind) :: gll_edges(np+1),phys_edges(nphys+1)
    type(quadrature_t) :: gll_pts
    ! setup (most be done on masterthread only) since all data is static
#if (! defined ELEMENT_OPENMP)
!OMP MASTER
#endif
    if (nphys_init/=nphys) then
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
          if (x2>1+tol) call abortmp('ERROR: did not find next intersection point')
             
          count=count+1
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
       if (count>nintersect) call abortmp('ERROR: nintersect was too small')
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
    endif
#if (! defined ELEMENT_OPENMP)
    !OMP ENDMASTER
    !OMP BARRIER
#endif


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


  function gradient_sphere(s,deriv,Dinv) result(ds)
!
!   input s:  scalar
!   output  ds: spherical gradient of s, lat-lon coordinates
!

    type (derivative_t)              :: deriv
    real(kind=real_kind), intent(in), dimension(2,2,np,np) :: Dinv
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
          ds(i,j,1)=Dinv(1,1,i,j)*v1(i,j) + Dinv(2,1,i,j)*v2(i,j)
          ds(i,j,2)=Dinv(1,2,i,j)*v1(i,j) + Dinv(2,2,i,j)*v2(i,j)
       enddo
    enddo

    end function gradient_sphere




  function gradient_sphere_wk(s,deriv,elem) result(ds)
!
!   input s:  scalar
!   output  ds: weak gradient, lat/lon coordinates
!
!   integral[ phivec dot grad(s) ] 
!        = phivec  dot  gradient_sphere(p) * spheremp(i,j) 
!   if phivec_contra(:,:,1)=phi (cardinal function)
!   if phivec_contra(:,:,2)=0
!        = gradp_covariant(:,:,1) * spheremp(i,j)
!   if phivec_contra(:,:,1)=0
!   if phivec_contra(:,:,2)=phi (cardinal function)
!        = gradp_covariant(:,:,2) * spheremp(i,j)
!
!   weak form is thus defined over all contra cardinal functions:   
!      integral[ div(phivec) s ] = sum  spheremp()* phi_x() * s
!                                  sum  spheremp()* phi_y() * s
!      
!

    type (derivative_t)              :: deriv
    type (element_t)              :: elem
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
          do j=1,np
             ! phi(m)_x  sum over first index, second index fixed at n
             dscov(m,n,1)=dscov(m,n,1)-(elem%mp(j,n)*elem%metdet(m,n)*s(j,n)*deriv%Dvv(m,j) )*rrearth
             ! phi(n)_y  sum over second index, 1st index fixed at m
             dscov(m,n,2)=dscov(m,n,2)-(elem%mp(m,j)*elem%metdet(m,n)*s(m,j)*deriv%Dvv(n,j) )*rrearth
          enddo
       enddo
    enddo

#if 0
    ! slow form:
    do m=1,np
       do n=1,np
          vcontra=0
          vcontra(m,n,1)=1

          ! contra->latlon:
          v(:,:,1)=(elem%D(1,1,:,:)*vcontra(:,:,1) + elem%D(1,2,:,:)*vcontra(:,:,2))
          v(:,:,2)=(elem%D(2,1,:,:)*vcontra(:,:,1) + elem%D(2,2,:,:)*vcontra(:,:,2))


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
          v(:,:,1)=(elem%D(1,1,:,:)*vcontra(:,:,1) + elem%D(1,2,:,:)*vcontra(:,:,2))
          v(:,:,2)=(elem%D(2,1,:,:)*vcontra(:,:,1) + elem%D(2,2,:,:)*vcontra(:,:,2))

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
    ! convert covariant to latlon 
    ds(:,:,1)=elem%Dinv(1,1,:,:)*dscov(:,:,1) + elem%Dinv(2,1,:,:)*dscov(:,:,2)
    ds(:,:,2)=elem%Dinv(1,2,:,:)*dscov(:,:,1) + elem%Dinv(2,2,:,:)*dscov(:,:,2)

    end function gradient_sphere_wk




  function ugradv_sphere(u,v,deriv,elem) result(ugradv)
!
!   input:  vectors u and v  (latlon coordinates)
!   output: vector  [ u dot grad ] v  (latlon coordinates)
!
    type (derivative_t)              :: deriv
    type (element_t)                 :: elem
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
       ugradv(:,:,component)=sum( dum_cart(:,:,:)*elem%vec_sphere2cart(:,:,:,component) ,3)
    end do

  end function ugradv_sphere



  function curl_sphere(s,deriv,elem) result(ds)
!
!   input s:  scalar  (assumed to be  s khat)
!   output  curl(s khat) in lat-lon coordinates
! 
!   This subroutine can be used to compute divergence free velocity fields,
!   since div(ds)=0
!
!    first compute:  
!    curl(s khat) = (1/jacobian) ( ds/dy, -ds/dx ) in contra-variant coordinates
!    then map to lat-lon
!
    type (derivative_t)              :: deriv
    type (element_t)                 :: elem
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
          ds(i,j,1)=(elem%D(1,1,i,j)*v1(i,j) + elem%D(1,2,i,j)*v2(i,j))/elem%metdet(i,j)
          ds(i,j,2)= (elem%D(2,1,i,j)*v1(i,j) + elem%D(2,2,i,j)*v2(i,j))/elem%metdet(i,j)
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
    type (derivative_t)              :: deriv
    type (element_t)                 :: elem
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
          vtemp(i,j,1)=(elem%Dinv(1,1,i,j)*v(i,j,1) + elem%Dinv(1,2,i,j)*v(i,j,2))
          vtemp(i,j,2)=(elem%Dinv(2,1,i,j)*v(i,j,1) + elem%Dinv(2,2,i,j)*v(i,j,2))
       enddo
    enddo

    do n=1,np
       do m=1,np

          div(m,n)=0
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
                gtemp(i,j,1)=(elem%D(1,1,i,j)*ggtemp(i,j,1) + elem%D(2,1,i,j)*ggtemp(i,j,2))
                gtemp(i,j,2)=(elem%D(1,2,i,j)*ggtemp(i,j,1) + elem%D(2,2,i,j)*ggtemp(i,j,2))
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
    type (derivative_t)              :: deriv
    type (element_t)                 :: elem
    real(kind=real_kind) :: result(np,np)

    ! Local
    real(kind=real_kind) :: ucontra(np,np,2)  ! in lat-lon coordinates
    integer i,j

    ! latlon->contra
    do j=1,np
       do i=1,np
          ucontra(i,j,1)=(elem%Dinv(1,1,i,j)*v(i,j,1) + elem%Dinv(1,2,i,j)*v(i,j,2))
          ucontra(i,j,2)=(elem%Dinv(2,1,i,j)*v(i,j,1) + elem%Dinv(2,2,i,j)*v(i,j,2))
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
    type (derivative_t)              :: deriv
    type (element_t)                 :: elem
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
             ucontra(i,j,1)=(elem%Dinv(1,1,i,j)*v(i,j,1) + elem%Dinv(1,2,i,j)*v(i,j,2))
             ucontra(i,j,2)=(elem%Dinv(2,1,i,j)*v(i,j,1) + elem%Dinv(2,2,i,j)*v(i,j,2))
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

    

  function vorticity_sphere(v,deriv,elem) result(vort)
!
!   input:  v = velocity in lat-lon coordinates
!   ouput:  spherical vorticity of v
!

    type (derivative_t)              :: deriv
    type (element_t)                 :: elem
    real(kind=real_kind), intent(in) :: v(np,np,2)

    real(kind=real_kind) :: vort(np,np)

    integer i
    integer j
    integer l
    
    real(kind=real_kind) ::  dvdx00
    real(kind=real_kind) ::  dudy00
    real(kind=real_kind) ::  vco(np,np,2)
    real(kind=real_kind) ::  vtemp(np,np)

    ! convert to covariant form
    do j=1,np
       do i=1,np
          vco(i,j,1)=(elem%D(1,1,i,j)*v(i,j,1) + elem%D(2,1,i,j)*v(i,j,2))
          vco(i,j,2)=(elem%D(1,2,i,j)*v(i,j,1) + elem%D(2,2,i,j)*v(i,j,2))
       enddo
    enddo

    do j=1,np
       do l=1,np

          dudy00=0.0d0
	  dvdx00=0.0d0

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

      type (derivative_t)              :: deriv
      type (element_t)                 :: elem
      real(kind=real_kind), intent(in) :: v(np,np,2)

      real(kind=real_kind) :: vort(np,np)

      integer i
      integer j
      integer l

      real(kind=real_kind) ::  dvdx00
      real(kind=real_kind) ::  dudy00
      real(kind=real_kind) ::  vco(np,np,2)
      real(kind=real_kind) :: vtemp(np,np)
      real(kind=real_kind) :: rdx
      real(kind=real_kind) :: rdy

! dx,dy are no longer initialized - this routine must be updated 
! see vorticity_sphere above
!      rdx=2.0D0/(elem%dx*rrearth) ! strong derivative inverse x length
!      rdy=2.0D0/(elem%dy*rrearth) ! strong derivative inverse y length

      ! convert to covariant form
                                                                    
      do j=1,np
         do i=1,np
            vco(i,j,1)=(elem%D(1,1,i,j)*v(i,j,1) + elem%D(2,1,i,j)*v(i,j,2))
            vco(i,j,2)=(elem%D(1,2,i,j)*v(i,j,1) + elem%D(2,2,i,j)*v(i,j,2))


         enddo
      enddo

                                                                                                               
      do j=1,np
         do l=1,np
          
            dudy00=0.0d0
            dvdx00=0.0d0

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
            vort(i,j)=elem%rmetdet(i,j)*(rdx*vort(i,j)-rdy*vtemp(i,j))
         end do 
      end do 
     
  end function vorticity_sphere_diag



  function divergence_sphere(v,deriv,elem) result(div)
!
!   input:  v = velocity in lat-lon coordinates
!   ouput:  div(v)  spherical divergence of v
!


    real(kind=real_kind), intent(in) :: v(np,np,2)  ! in lat-lon coordinates
    type (derivative_t)              :: deriv
    type (element_t)                 :: elem
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
          gv(i,j,1)=elem%metdet(i,j)*(elem%Dinv(1,1,i,j)*v(i,j,1) + elem%Dinv(1,2,i,j)*v(i,j,2))
          gv(i,j,2)=elem%metdet(i,j)*(elem%Dinv(2,1,i,j)*v(i,j,1) + elem%Dinv(2,2,i,j)*v(i,j,2))
       enddo
    enddo

    ! compute d/dx and d/dy         
    do j=1,np
       do l=1,np
          dudx00=0.0d0
          dvdy00=0.0d0
          do i=1,np
             dudx00 = dudx00 + deriv%Dvv(i,l  )*gv(i,j  ,1)
             dvdy00 = dvdy00 + deriv%Dvv(i,l  )*gv(j  ,i,2)
          end do
          div(l  ,j  ) = dudx00
          vvtemp(j  ,l  ) = dvdy00
       end do
    end do

    do j=1,np
       do i=1,np
          div(i,j)=(div(i,j)+vvtemp(i,j))*(elem%rmetdet(i,j)*rrearth)
       end do
    end do
    
  end function divergence_sphere


!three types of viscosity are supported right now:
! (1) const hv, i.e., the operator nu * (\div \grad)^hypervis_order
! (2) variable-within-element (or just variable) hv, the operator nu * (viscosity \div \grad )^hypervis_order
! (3) tensor hv,  nu * ( \div * tensor * \grad )^hypervis_order

! the switch between (2) and (3) is here, in preproc directive #TENSORHV
! the switch between (1) and (2,3) is in namelist variable:
! hypervis_power =0 for const hv (1),  <>0 for variable hv (2) and tensor hv (3)
! if using (2), it is required to set also fine_ne, max_hypervis_courant

#undef TENSORHV
!#define TENSORHV


  function laplace_sphere_wk(s,deriv,elem,viscosity) result(laplace)
!
!   input:  s = scalar
!   ouput:  -< grad(PHI), grad(s) >   = weak divergence of grad(s)
!     note: for this form of the operator, grad(s) does not need to be made C0
!            
    real(kind=real_kind), intent(in) :: s(np,np) 
    real(kind=real_kind), pointer, dimension(:,:) :: viscosity
    type (derivative_t)              :: deriv
    type (element_t)                 :: elem
    real(kind=real_kind)             :: laplace(np,np)
    real(kind=real_kind)             :: laplace2(np,np)
    integer i,j

    ! Local
    real(kind=real_kind) :: grads(np,np,2), oldgrads(np,np,2)

    grads=gradient_sphere(s,deriv,elem%Dinv)
 
#ifndef TENSORHV
! const or variable viscosity, (1) or (2)
    if (ASSOCIATED(viscosity)) then
        grads(:,:,1) = grads(:,:,1)*viscosity(:,:)
        grads(:,:,2) = grads(:,:,2)*viscosity(:,:)
    end if
#else
! tensor hv, (3)
    if (ASSOCIATED(viscosity)) then
      oldgrads=grads
      do j=1,np
	do i=1,np
	  grads(i,j,1) = sum(oldgrads(i,j,:)*elem%tensorVisc(1,:,i,j))
	  grads(i,j,2) = sum(oldgrads(i,j,:)*elem%tensorVisc(2,:,i,j))
	end do
      end do
    endif
#endif

    ! note: divergnece_sphere and divergence_sphere_wk are identical *after* bndry_exchange
    ! if input is C_0.  Here input is not C_0, so we should use divergence_sphere_wk().  

    laplace=divergence_sphere_wk(grads,deriv,elem)

  end function laplace_sphere_wk


  function vlaplace_sphere_wk(v,deriv,elem,viscosity,nu_ratio) result(laplace)
!
!   input:  v = vector in lat-lon coordinates
!   ouput:  weak laplacian of v, in lat-lon coordinates

    real(kind=real_kind), intent(in) :: v(np,np,2) 
    real(kind=real_kind), pointer, dimension(:,:) :: viscosity
    type (derivative_t)              :: deriv
    type (element_t)                 :: elem
    real(kind=real_kind), optional :: nu_ratio
    real(kind=real_kind) :: laplace(np,np,2)

    if (which_vlaplace .eq. 2) then
      laplace=cartesian_laplace_sphere_wk(v,deriv,elem,viscosity,nu_ratio)
    else
      laplace=vector_identities_laplace_sphere_wk(v,deriv,elem,viscosity,nu_ratio)
    endif

  end function vlaplace_sphere_wk



  function cartesian_laplace_sphere_wk(v,deriv,elem,viscosity,nu_ratio) result(laplace)
!
!   input:  v = vector in lat-lon coordinates
!   ouput:  weak laplacian of v, in lat-lon coordinates

    real(kind=real_kind), intent(in) :: v(np,np,2) 
    real(kind=real_kind), pointer, dimension(:,:) :: viscosity
    type (derivative_t)              :: deriv
    type (element_t)                 :: elem
    real(kind=real_kind) :: laplace(np,np,2)
    real(kind=real_kind), optional :: nu_ratio
    ! Local

    integer component
    real(kind=real_kind) :: dum_cart(np,np,3)


    ! latlon -> cartesian
    do component=1,3
       dum_cart(:,:,component)=sum( elem%vec_sphere2cart(:,:,component,:)*v(:,:,:) ,3)
    end do

    ! Do laplace on cartesian comps
    do component=1,3
       dum_cart(:,:,component) = laplace_sphere_wk(dum_cart(:,:,component),deriv,elem,viscosity)
    enddo

    ! cartesian -> latlon
    do component=1,2
       ! vec_sphere2cart is its own pseudoinverse.
       laplace(:,:,component)=sum( dum_cart(:,:,:)*elem%vec_sphere2cart(:,:,:,component) ,3)
    end do 

  end function cartesian_laplace_sphere_wk



  function vector_identities_laplace_sphere_wk(v,deriv,elem,viscosity,nu_ratio) result(laplace)
!
!   input:  v = vector in lat-lon coordinates
!   ouput:  weak laplacian of v, in lat-lon coordinates
!   note: integrals must be performed in contra variant coordinates,
!         convert to lat-lon after computing integral
!
!   laplace(v) =  grad(div) -  curl(vor) 
!   weak form   < PHI , grad(div) > - < PHI, curl(vor*khat) >    
!   by parts:   -< div(PHI) div >   - < vor curl(PHI) >      
!             
!
! used vector identity: div(F cross G)=G dot curl F - F dot curl G
!                  OR:    < G, curl F> = < F, curl G >
!
!  NOTE: for the equation:   < PHI, LAPLACE > =   -< div(PHI) div >   - < vor curl(PHI) >      
!  if LAPLACE is in covarient, we test with PHI = (1,0) and (0,1) in contra-variant
!  if LAPLACE is in contra-varient, we test with PHI = (1,0) and (0,1) in co-varient
!
!  compute with two different test functions:
!  (then transform back to lat-lon (since we output in lat-lon))
!  test function 1:  contra:  (phi,0)        covarient:   (met11 phi, met21 phi)
!  test function 2:  contra:  (0,phi)        covariant:   (met12 phi, met22 phi)
!  
!  div acts on contra components
!   < div(PHI) div >  =  <  1/g (g phi)_x div >  = <  phi_x div >    (test 1)
!                        <  1/g (g phi)_y div >  = <  phi_y div >    (test 2)
!
!
!  curl  acts on co-variant
!   < curl(PHI) vor >  =   <  1/g [-(met11 phi )_y + (met21 phi)_x  ] vor > 
!                          <  1/g [-(met12 phi )_y + (met22 phi)_x  ] vor >
!                      =   <  1/g [-met11 phi_y + met21 phi_x  ] vor > 
!                          <  1/g [-met12 phi_y + met22 phi_x  ] vor >
!             
!
!  curl acting on co-variant test functions:
!  test function 1:  co:  (phi,0)        curl(PHI) = -phi_y
!  test function 2:  co:  (0,phi)        curl(PHI) = phi_x
!   < curl(PHI) vor >  =  <  1/g -phi_y vor  >  = <  1/g -phy_y vor  >    (test 1)
!                         <  1/g phi_x vor >    = <  1/g phi_x vor >      (test 2)
!
!                         
!   
!
! NOTE:  SEM can compute (g11 phi)_y in two ways:
!        project, than derivative:     met11 phi_y   (we are using this formula)
!        expand first:                 met11 phi_y  + met11_y phi
!
!
! compare to weak divergence:  < grad(PHI), v >
!  if v is in contra,  grad(PHI) in covariant = (phi_x, phi_y)
!  < phi_x v1  + phi_y v2 >   THUS:  < phi_x div > loop should look like d/dx
!  loop in divergence_sphere_wk()
!
! NOTE: dont forget < u v > = integral g*u*v = sum mv()*metdet()*u*v  
!        with g=metdet(), spheremv=mv()*metdet()
!
!  < phy_y div > = sum spheremv * phy_y * div
!  < 1/g F >     = sum mv * F
!
    real(kind=real_kind), intent(in) :: v(np,np,2) 
    real(kind=real_kind), pointer, dimension(:,:) :: viscosity
    type (derivative_t)              :: deriv
    type (element_t)                 :: elem
    real(kind=real_kind) :: laplace(np,np,2)
    real(kind=real_kind), optional :: nu_ratio
    ! Local

    integer i,j,l,m,n
    real(kind=real_kind) :: vor(np,np),div(np,np)
    real(kind=real_kind) :: v1,v2,div1,div2,vor1,vor2,phi_x,phi_y

    div=divergence_sphere(v,deriv,elem)
    vor=vorticity_sphere(v,deriv,elem)


    if (ASSOCIATED(viscosity)) then
       div = div*viscosity
       vor = vor*viscosity
    end if
    if (present(nu_ratio)) div = nu_ratio*div

#undef DIVCONTRA_VORCO

    do n=1,np
       do m=1,np

          div1=0; div2=0;
          vor1=0; vor2=0; 

          do l=1,np
             phi_x=deriv%Dvv(m,l)*rrearth ! (m,n) cardinal function, d/dx  
                                        ! phi_x(i,j) = 0  for j<>n.  so treat this as phi_x(i,n)

             phi_y=deriv%Dvv(n,l)*rrearth ! (m,n) cardinal function, d/dy
                                        ! phi_y(i,j) = 0  for i<>m.  so treat this as phi_x(m,j)

             div1=div1 + elem%spheremp(l,n)*div(l,n)*phi_x
             div2=div2 + elem%spheremp(m,l)*div(m,l)*phi_y
             
#ifdef DIVCONTRA_VORCO
             vor1=vor1 - elem%mp(m,l)*vor(m,l)*phi_y
             vor2=vor2 + elem%mp(l,n)*vor(l,n)*phi_x
#else
             vor1=vor1 - elem%mp(m,l)*vor(m,l)*elem%met(1,1,m,l)*phi_y &
                         + elem%mp(l,n)*vor(l,n)*elem%met(2,1,l,n)*phi_x

             vor2=vor2 - elem%mp(m,l)*vor(m,l)*elem%met(1,2,m,l)*phi_y &
                         + elem%mp(l,n)*vor(l,n)*elem%met(2,2,l,n)*phi_x

#endif

          enddo
#ifdef DIVCONTRA_VORCO
          v1=-div1
          v2=-div2

          !  (v1,v2) = divergence componet tested against contra-variant, so result is CO-variant
          laplace(m,n,1)=elem%Dinv(1,1,m,n)*v1 + elem%Dinv(2,1,m,n)*v2   ! co->latlon
          laplace(m,n,2)=elem%Dinv(1,2,m,n)*v1 + elem%Dinv(2,2,m,n)*v2   ! co->latlon

          v1=-vor1
          v2=-vor2
          !  (v1,v2) = vorticity component tested against co-variant.  so result is CONTRA 
          laplace(m,n,1)=laplace(m,n,1) + elem%D(1,1,m,n)*v1 + elem%D(1,2,m,n)*v2   ! contra->latlon
          laplace(m,n,2)=laplace(m,n,2) + elem%D(2,1,m,n)*v1 + elem%D(2,2,m,n)*v2   ! contra->latlon
#else
          v1=-( div1 + vor1 )
          v2=-( div2 + vor2 )

          !  (v1,v2) = RHS tested agains contra-variant delta functions, so result is CO-varient
          laplace(m,n,1)=elem%Dinv(1,1,m,n)*v1 + elem%Dinv(2,1,m,n)*v2   ! co->latlon
          laplace(m,n,2)=elem%Dinv(1,2,m,n)*v1 + elem%Dinv(2,2,m,n)*v2   ! co->latlon
#endif
          ! add in correction so we dont damp rigid rotation
#define UNDAMPRR
#ifdef UNDAMPRR
          laplace(m,n,1)=laplace(m,n,1) + 2*elem%spheremp(m,n)*v(m,n,1)*(rrearth**2)
          laplace(m,n,2)=laplace(m,n,2) + 2*elem%spheremp(m,n)*v(m,n,2)*(rrearth**2)
#endif
       enddo
    enddo
  end function vector_identities_laplace_sphere_wk



!-----------------------------------------------------------------------------------


  function gll_to_dgmodal(p,deriv) result(phat)
!
!   input:  v = velocity in lat-lon coordinates
!   ouput:  phat = Legendre coefficients
!
!   Computes  < g dot p  > = SUM  g(i,j) p(i,j) w(i) w(j)
!   (the quadrature approximation on the *reference element* of the integral of p against
!    all Legendre polynomials up to degree npdg
!
!   for npdg < np, this routine gives the (exact) modal expansion of p/spheremp()
!
    real(kind=real_kind), intent(in) :: p(np,np) 
    type (derivative_t)              :: deriv
    real(kind=real_kind) :: phat(npdg,npdg)

    ! Local
    integer i,j,m,n
    real(kind=real_kind) :: A(np,npdg)
    A=0
    phat=0

    ! N^3 tensor product formulation:
    do m=1,npdg
    do j=1,np
    do i=1,np
       A(j,m)=A(j,m)+( p(i,j)*deriv%Mvv_twt(i,i)*deriv%Mvv_twt(j,j)  )*deriv%legdg(m,i)
    enddo
    enddo
    enddo

    do n=1,npdg
    do m=1,npdg
    do j=1,np
       phat(m,n)=phat(m,n)+A(j,m)*deriv%legdg(n,j)
    enddo
    enddo
    enddo
    
#if 0
    do m=1,npdg
       do n=1,npdg
          do j=1,np
             do i=1,np
                gmn = deriv%legdg(m,i)*deriv%legdg(n,j) ! basis function
                phat(m,n)=phat(m,n)+gmn*p(i,j)*deriv%Mvv_twt(i,i)*deriv%Mvv_twt(j,j)
             enddo
          enddo
       enddo
    enddo
#endif
  end function

  function dgmodal_to_gll(phat,deriv) result(p)
!
!   input:  phat = coefficients of Legendre expansion
!   ouput:  p    = sum expansion to evaluate phat at GLL points
!
    real(kind=real_kind) :: p(np,np) 
    type (derivative_t)  :: deriv
    real(kind=real_kind) :: phat(npdg,npdg)
    ! Local
    integer i,j,m,n
    real(kind=real_kind) :: A(npdg,np)

    p(:,:)=0
    ! tensor product version
    A=0
    do i=1,np
    do n=1,npdg
    do m=1,npdg
       A(n,i)=A(n,i)+phat(m,n)*deriv%legdg(m,i)
    enddo
    enddo
    enddo
    do j=1,np
    do i=1,np
    do n=1,npdg
       p(i,j) = p(i,j)+A(n,i)*deriv%legdg(n,j)
    enddo
    enddo
    enddo

#if 0
    do j=1,np
       do i=1,np
          do m=1,npdg
             do n=1,npdg
                p(i,j)=p(i,j)+phat(m,n)*deriv%legdg(m,i)*deriv%legdg(n,j) 
             enddo
          enddo
       enddo
    enddo
#endif
  end function

  ! Given a field defined on the unit element, [-1,1]x[-1,1]
  ! sample values, sampled_val, and integration weights, metdet,
  ! at a number, np, of Gauss-Lobatto-Legendre points. Divide
  ! the square up into intervals by intervals sub-squares so that
  ! there are now intervals**2 sub-cells.  Integrate the 
  ! function defined by sampled_val and metdet over each of these
  ! sub-cells and return the integrated values as an 
  ! intervals by intervals matrix.
  !
  ! Efficiency is obtained by computing and caching the appropriate
  ! integration matrix the first time the function is called.
  function subcell_integration(sampled_val, metdet, np, intervals) result(values)

    implicit none

    integer              , intent(in)  :: np
    integer              , intent(in)  :: intervals
    real (kind=real_kind), intent(in)  :: sampled_val(np,np)
    real (kind=real_kind), intent(in)  :: metdet     (np,np)
    real (kind=real_kind)              :: values(intervals,intervals)

    real (kind=real_kind)              :: V          (np,np)
    integer i,j

    V  = sampled_val * metdet


    if (.not.ALLOCATED(integration_matrix)      .or. &
        SIZE(integration_matrix,1).ne.intervals .or. &
        SIZE(integration_matrix,2).ne.np) then
      call allocate_subcell_integration_matrix(np,intervals)
    end if

    ! Multiply the sampled values by the weighted jacobians.  
    ! Symmetry allows us to write this as J^t V J
    ! where J is a vector.  

    values = MATMUL(integration_matrix, &
             MATMUL(V,TRANSPOSE(integration_matrix)))

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

  end subroutine allocate_subcell_integration_matrix



end module derivative_mod








