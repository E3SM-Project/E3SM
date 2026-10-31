#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

! ---------------------------------------------------------------------------------------------------
! Geometry-independent element metric terms.
!
! metric_atomic() computes D/Dinv/met/metinv/metdet, the element length scales, and the tensor
! (hyper)viscosity coefficients tensorVisc and tensorVisc_2.  The math is identical for spherical
! and planar geometry; only two things differ:
!
!   1) the reference-element-to-physical map, supplied by the caller as the dmap_fn argument
!      (cube_mod::Dmap for the sphere, planar_mod::plane_Dmap for the plane), and
!   2) the length scaling, which enters only through scale_factor / scale_factor_inv
!      (rearth / 1/rearth for the sphere, 1 / 1 for the plane -- see namelist_mod::readnl).
!
! This module deliberately sits BELOW cube_mod and planar_mod in the dependency graph: taking the
! map as a procedure argument is what lets both of them use it without a circular dependency.
! ---------------------------------------------------------------------------------------------------
module metric_mod
  use kinds,                  only : real_kind, longdouble_kind
  use dimensions_mod,         only : np
  use element_mod,            only : element_t
  use coordinate_systems_mod, only : cartesian3D_t, cartesian2D_t
  use control_mod,            only : hypervis_scaling, laplace_scaling, cubed_sphere_map, geometry
  use physical_constants,     only : scale_factor, scale_factor_inv
  use parallel_mod,           only : abortmp

  implicit none
  private

  public :: metric_atomic
  public :: dmap_i

  ! -------------------------------------------------------------------
  ! Signature shared by cube_mod::Dmap and planar_mod::plane_Dmap.
  ! plane_Dmap ignores every argument but D (its map is affine), so the
  ! caller can pass the gll/cartp/facenum arguments unconditionally.
  abstract interface
    subroutine dmap_i(D, a, b, corners3D, ref_map, cartp, facenum)
      import :: real_kind, cartesian3D_t, cartesian2D_t, np
      real (kind=real_kind), intent(out) :: D(2,2)
      real (kind=real_kind), intent(in)  :: a,b
      type (cartesian3D_t)               :: corners3D(4)
      integer                            :: ref_map
      type (cartesian2D_t), optional     :: cartp(np,np)
      integer, optional                  :: facenum
    end subroutine dmap_i
  end interface

contains

  ! =========================================
  ! metric_atomic:
  !
  ! Initialize element metric terms (atomic), for either geometry:
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
  ! so if we want to tweak the mapping by a factor alpha (so he weights add up to 4pi, for example)
  ! we take:
  !    NEW       OLD     
  !       D = sqrt(alpha) D  and then rederive all quantities.  
  !    detD = alpha detD
  !    
  ! where alpha = domain_size/SEMarea, SEMarea = global sum elem(ie)%mv(i,j)*elem(ie)%metdet(i,j)
  ! 
  ! =========================================

  subroutine metric_atomic(elem,gll_points,alpha,dmap_fn)

    type (element_t) :: elem
    real(kind=real_kind) :: alpha
    real (kind=longdouble_kind)      :: gll_points(np)
    procedure(dmap_i)                :: dmap_fn   ! Dmap (sphere) or plane_Dmap (plane)
    ! Local variables
    integer :: i,j

    real (kind=real_kind) :: norm
    real (kind=real_kind) :: detD      ! determinant of vector field mapping matrix.  

    real (kind=real_kind) :: x1        ! 1st cube face coordinate
    real (kind=real_kind) :: x2        ! 2nd cube face coordinate
    real (kind=real_kind) :: M(2,2),E(2,2),eig(2),DE(2,2),DEL(2,2),V(2,2), lamStar1, lamStar2
    integer :: imaxM(2)
    real (kind=real_kind) :: l1, l2, min_svd,max_svd,max_normDinv

    
    ! ==============================================
    ! Initialize differential mapping operator
    ! to and from vector fields on the sphere to 
    ! contravariant vector fields on the cube
    ! i.e. dM/dx^i in Sadourney (1972) and it's 
    ! inverse
    ! ==============================================

    max_svd = 0.0d0
    max_normDinv = 0.0d0
    min_svd = 1d99
    do j=1,np
       do i=1,np
          x1=gll_points(i)
          x2=gll_points(j)
          call dmap_fn(elem%D(i,j,:,:),x1,x2,elem%corners3D,cubed_sphere_map,elem%cartp,elem%facenum)


          ! Numerical metric tensor based on analytic D: met = D^T times D
          ! (D maps between sphere and reference element)
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
             call abortmp('Impossible error in metric_mod.F90::metric_atomic(), geometry='//trim(geometry))
          endif

          ! the other eigenvector is orthgonal:
          E(1,2)=-E(2,1)
          E(2,2)= E(1,1)

!normalize columns
	  E(:,1)=E(:,1)/sqrt(sum(E(:,1)*E(:,1))); 
	  E(:,2)=E(:,2)/sqrt(sum(E(:,2)*E(:,2))); 


! OBTAINING TENSOR FOR HV:

! Instead of the traditional scalar Laplace operator \grad \cdot \grad
! we introduce \grad \cdot V \grad
! where V = D E LAM LAM^* E^T D^T. 
! Recall (metric_tensor)^{-1}=(D^T D)^{-1} = E LAM E^T.
! Here, LAM = diag( 4/((np-1)dx)^2 , 4/((np-1)dy)^2 ) = diag(  4/(dx_elem)^2, 4/(dy_elem)^2 )
! Note that metric tensors and LAM correspondingly are quantities on a unit sphere.

! This motivates us to use V = D E LAM LAM^* E^T D^T
! where LAM^* = diag( nu1, nu2 ) where nu1, nu2 are HV coefficients scaled like (dx)^{hv_scaling/2}, (dy)^{hv_scaling/2}.
! (Halves in powers come from the fact that HV consists of two Laplace iterations.)

! Originally, we took LAM^* = diag(
!  1/(eig(1)**(hypervis_scaling/4.0d0))*(rearth**(hypervis_scaling/2.0d0))
!  1/(eig(2)**(hypervis_scaling/4.0d0))*(rearth**(hypervis_scaling/2.0d0)) ) = 
!  = diag( lamStar1, lamStar2)
!  \simeq ((np-1)*dx_sphere / 2 )^hv_scaling/2 = SQRT(OPERATOR_HV)
! because 1/eig(...) \simeq (dx_on_unit_sphere)^2 .
! Introducing the notation OPERATOR = lamStar^2 is useful for conversion formulas.

! This leads to the following conversion formula: nu_const is nu used for traditional HV on uniform grids
! nu_tensor = nu_const * OPERATOR_HV^{-1}, so
! nu_tensor = nu_const *((np-1)*dx_sphere / 2 )^{ - hv_scaling} or
! nu_tensor = nu_const *(2/( (np-1) * dx_sphere) )^{hv_scaling} .
! dx_sphere = 2\pi *rearth/(np-1)/4/NE
! [nu_tensor] = [meter]^{4-hp_scaling}/[sec]

! (1) Later developments:
! Apply tensor V only at the second Laplace iteration. Thus, LAM^* should be scaled as (dx)^{hv_scaling}, (dy)^{hv_scaling},
! see this code below:
!          DEL(1:2,1) = (lamStar1**2) *eig(1)*DE(1:2,1)
!          DEL(1:2,2) = (lamStar2**2) *eig(2)*DE(1:2,2)

! (2) Later developments:
! (The derivation below is written for the sphere, where scale_factor==rearth.
!  For planar geometry scale_factor==1 and the rearth factors simply drop out.)
! Bringing [nu_tensor] to 1/[sec]:
!	  lamStar1=1/(eig(1)**(hypervis_scaling/4.0d0)) *(scale_factor**2.0d0)
!	  lamStar2=1/(eig(2)**(hypervis_scaling/4.0d0)) *(scale_factor**2.0d0)
! OPERATOR_HV = ( (np-1)*dx_unif_sphere / 2 )^{hv_scaling} * rearth^4
! Conversion formula:
! nu_tensor = nu_const * OPERATOR_HV^{-1}, so
! nu_tensor = nu_const *( 2*rearth /((np-1)*dx))^{hv_scaling} * rearth^{-4.0}.

! For the baseline coefficient nu=1e15 for NE30, 
! nu_tensor=7e-8 (BUT RUN TWICE AS SMALL VALUE FOR NOW) for hv_scaling=3.2
! and 
! nu_tensor=1.3e-6 for hv_scaling=4.0.


!matrix D*E
          DE(1,1)=sum(elem%D(i,j,1,:)*E(:,1))
          DE(1,2)=sum(elem%D(i,j,1,:)*E(:,2))
          DE(2,1)=sum(elem%D(i,j,2,:)*E(:,1))
          DE(2,2)=sum(elem%D(i,j,2,:)*E(:,2))

	  lamStar1=1/(eig(1)**(hypervis_scaling/4.0d0)) *(scale_factor**2.0d0)
	  lamStar2=1/(eig(2)**(hypervis_scaling/4.0d0)) *(scale_factor**2.0d0)

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


!
!         Tensor with scalings=2 for regular laplace operator (spnge layer)
	  lamStar1=1/(eig(1)**(laplace_scaling/2.0d0))*(scale_factor**2.0d0)
	  lamStar2=1/(eig(2)**(laplace_scaling/2.0d0))*(scale_factor**2.0d0)
          DEL(1:2,1) = lamStar1 *eig(1)*DE(1:2,1)
          DEL(1:2,2) = lamStar2 *eig(2)*DE(1:2,2)
          V(1,1)=sum(DEL(1,:)*DE(1,:))
          V(1,2)=sum(DEL(1,:)*DE(2,:))
          V(2,1)=sum(DEL(2,:)*DE(1,:))
          V(2,2)=sum(DEL(2,:)*DE(2,:))

	  elem%tensorVisc_2(i,j,:,:)=V(:,:)

       end do
    end do

!    see Paul Ullrich writeup:   
!    max_normDinv might be a tighter bound than max_svd for deformed elements
!    max_svd >= max_normDinv/sqrt(2), with equality holding if |g^x| = |g^y| 
!    elem%normDinv=max_normDinv/sqrt(2)     

    ! this norm is consistent with length scales defined below:
    elem%normDinv=max_svd


    ! compute element length scales, based on SVDs, in km:
    elem%dx_short = 1.0d0/(max_svd*0.5d0*dble(np-1)*scale_factor_inv*1000.0d0)
    elem%dx_long  = 1.0d0/(min_svd*0.5d0*dble(np-1)*scale_factor_inv*1000.0d0)

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

end module metric_mod
