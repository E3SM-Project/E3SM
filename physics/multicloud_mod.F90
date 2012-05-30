#include "timer.h"
module multicloud_mod

  use kinds, only              : real_kind
  use dimensions_mod, only     : nelemd,nlev, nlevp, np, qsize
  use physical_constants, only : rgas, Rwater_vapor, kappa, g, rearth, rrearth, p0
  use derivative_mod, only     : gradient, vorticity, gradient_wk, derivative_t, divergence
  use element_mod, only        : element_t
  use filter_mod, only         : filter_t, filter_P, preq_filter
  use hybvcoord_mod, only      : hvcoord_t
  use time_mod, only           : TimeLevel_t, smooth
  use prim_si_mod, only        : preq_pressure
  use diffusion_mod, only      : scalar_diffusion, diffusion_init
  use control_mod, only        : integration, test_case, filter_freq_advection
#ifndef CAM
  use physics_mod, only        : elem_physics_t, Prim_Condense,getsurfpress,Temp2PotTemp
#endif
  use edge_mod, only           : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack, initedgebuffer
  use hybrid_mod, only         : hybrid_t
  use bndry_mod, only          : bndry_exchangev
  use column_types_mod, only    : ColumnModelMulticloud_t
  use physical_constants, only : kappa
  use aquaplanet, only         : cool_max, presc_cooling_mc

  use interpolate_mod, only : interpolate_t, setup_latlon_interp, interpdata_t, &
       get_interp_parameter, get_interp_lat, get_interp_lon, interpolate_scalar, interpolate_vector


  implicit none

  private

  public :: Multicloud_init
  public :: Prim_Advance_Multicloud
  public :: verticalprojectionvector
  public :: verticalprojection
  public :: verticalprojection1D
  public :: VerticalIntegral1D
  public :: ProjVandDQ
  public :: LambdaSwitch
  public :: ApplyShearDamping
  public :: InitShearDamping

  type (EdgeBuffer_t) :: edgeMc	
  type (EdgeBuffer_t) :: edgeS1,edgeS2
  type (EdgeBuffer_t) :: edge3	
  type (EdgeBuffer_t) :: edge4

  ! Interpolation hidden data

  type(interpdata_t),   allocatable :: interpdata(:)
  real(kind=real_kind), allocatable :: udamp(:),udamp_copy(:),lat(:)

contains

  subroutine Multicloud_Init()

    call Prim_Advance_Multicloud_init()

  end subroutine Multicloud_Init

  function LambdaSwitch(T,Tm,Tp,A,B,Lstar) result(lambda)

    real (kind=real_kind),intent(in) :: T,Tm,Tp,A,B,Lstar
    real (kind=real_kind)            :: lambda

    !lambda=0.0D0

    if(T<Tm)then
       lambda=Lstar
    else if(T>Tp)then
       lambda=1.0D0
    else
       lambda=A*T+B
    end if

  end function LambdaSwitch

  subroutine VerticalProjectionVector(v,proj,w,lnps,hvcoord)

    real (kind=real_kind), dimension(np,np,2,nlev) :: v
    real (kind=real_kind), dimension(nlev)         :: w
    real (kind=real_kind), dimension(np,np)        :: lnps
    type (hvcoord_t)                               :: hvcoord

    ! LOCAL
    integer :: i,j,k
    real (kind=real_kind)                          :: delp
    real (kind=real_kind), dimension(np,np,nlev+1) :: ph
    real (kind=real_kind), dimension(np,np,2)      :: proj
    real (kind=real_kind), dimension(np,np)        :: ps
    real (kind=real_kind), dimension(np,np,nlev)   :: dp

    ps = Getsurfpress(lnps)

    do k=1,nlev+1
       do j=1,np
          do i=1,np
             ph(i,j,k)   = hvcoord%hyai(k)*hvcoord%ps0 + hvcoord%hybi(k)*ps(i,j)
          end do
       end do
    end do

    do k=1,nlev
       do j=1,np
          do i=1,np
             delp = ph(i,j,nlev+1) - ph(i,j,1)
             dp(i,j,k)  = (ph(i,j,k+1) - ph(i,j,k))/delp
          end do
       end do
    enddo

    proj=0.0D0

    do k=1,nlev
       do j=1,np
          do i=1,np
             proj(i,j,1) = proj(i,j,1) + dp(i,j,k)*v(i,j,1,k)*w(k)
             proj(i,j,2) = proj(i,j,2) + dp(i,j,k)*v(i,j,2,k)*w(k)
          enddo
       enddo
    enddo

  end subroutine VerticalProjectionVector

  subroutine VerticalProjection(v,proj,w,lnps,hvcoord)

    real (kind=real_kind), dimension(np,np,nlev)   :: v
    real (kind=real_kind), dimension(nlev)         :: w
    real (kind=real_kind), dimension(np,np)        :: lnps
    type (hvcoord_t)                               :: hvcoord

    ! LOCAL
    integer :: i,j,k
    real (kind=real_kind)                          :: delp
    real (kind=real_kind), dimension(np,np,nlev+1) :: ph
    real (kind=real_kind), dimension(np,np)        :: proj
    real (kind=real_kind), dimension(np,np)        :: ps
    real (kind=real_kind), dimension(np,np,nlev)   :: dp

    ps = Getsurfpress(lnps)

    !print *,"INSIDE PROJECTION  ", sum(abs(v(:,:,:)))

    do k=1,nlev+1
       do j=1,np
          do i=1,np
             ph(i,j,k)   = hvcoord%hyai(k)*hvcoord%ps0 + hvcoord%hybi(k)*ps(i,j)
          end do
       end do
    end do

    do k=1,nlev
       do j=1,np
          do i=1,np
             delp = ph(i,j,nlev+1) - ph(i,j,1)
             dp(i,j,k)  = (ph(i,j,k+1) - ph(i,j,k))/delp
          end do
       end do
    enddo

    proj(:,:) =0.0D0

    do k=1,nlev
       do j=1,np
          do i=1,np
             proj(i,j) = proj(i,j) + dp(i,j,k)*v(i,j,k)*w(k)
          enddo
       enddo
    enddo

  end subroutine VerticalProjection


  subroutine VerticalProjection1D(v,proj,w,lnps,hvcoord)

    real (kind=real_kind), dimension(nlev) :: v
    real (kind=real_kind), dimension(nlev) :: w
    real (kind=real_kind)                  :: lnps
    type (hvcoord_t)                       :: hvcoord

    ! LOCAL
    integer :: i,j,k
    real (kind=real_kind)                    :: delp
    real (kind=real_kind), dimension(nlev+1) :: ph
    real (kind=real_kind)                    :: proj
    real (kind=real_kind)                    :: ps
    real (kind=real_kind), dimension(nlev)   :: dp

    logical                                  :: Debug=.FALSE.

    ps = exp(lnps)

    if(Debug)print *,">>> PS = ", ps

    do k=1,nlev+1
       ph(k)   = hvcoord%hyai(k)*hvcoord%ps0 + hvcoord%hybi(k)*ps
    end do

    if(Debug)print *,">>> PH = ", ph

    delp = 1.0D0/(ph(nlev+1) - ph(1))

    if(Debug)print *,"delp = ", delp

    do k=1,nlev
       dp(k)  = (ph(k+1) - ph(k))*delp
    end do


    proj=0.0D0

    do k=1,nlev
       proj = proj + dp(k)*v(k)*w(k)
    enddo

  end subroutine VerticalProjection1D


  !===================================================
  ! psi_n(p)  = - 1/(p_b-p_t) \Int(p_t,p) phi_n(p')dp'
  !===================================================
  subroutine VerticalIntegral1D(v,proj,w,lnps,hvcoord)

    real (kind=real_kind), dimension(nlev) :: v
    real (kind=real_kind), dimension(nlev) :: w
    real (kind=real_kind)                  :: lnps
    type (hvcoord_t)                       :: hvcoord

    ! LOCAL
    integer :: i,j,k
    real (kind=real_kind)                    :: delp
    real (kind=real_kind), dimension(nlev+1) :: ph
    real (kind=real_kind)                    :: ps
    real (kind=real_kind), dimension(nlev)   :: proj
    real (kind=real_kind), dimension(nlev)   :: dp

    logical                                  :: Debug=.FALSE.

    ps = exp(lnps)

    do k=1,nlev+1
       ph(k)   = hvcoord%hyai(k)*hvcoord%ps0 + hvcoord%hybi(k)*ps
    end do

    delp = -1.0D0/(ph(nlev+1) - ph(1))

    do k=1,nlev
       dp(k)  = (ph(k+1) - ph(k))*delp
    end do


    proj(1)=0.0D0

    do k=2,nlev
       proj(k) = proj(k-1) + dp(k)*v(k)*w(k)
    enddo

  end subroutine VerticalIntegral1D


  ! ================================================
  ! Computes the integral of [v dq/deta]
  ! ================================================
  subroutine ProjVandDQ(v,q,proj,lnps,hvcoord)

    real (kind=real_kind), dimension(np,np,nlev)   :: v
    real (kind=real_kind), dimension(np,np,nlev)   :: q
    real (kind=real_kind), dimension(np,np)        :: lnps
    type (hvcoord_t)                               :: hvcoord

    ! LOCAL
    integer :: i,j,k
    real (kind=real_kind), dimension(np,np,nlev+1) :: ph
    real (kind=real_kind), dimension(np,np)        :: proj
    real (kind=real_kind), dimension(np,np)        :: ps
    real (kind=real_kind), dimension(np,np,nlev)   :: dp
    real (kind=real_kind), dimension(np,np,nlev)   :: dq

    !ps = Getsurfpress(lnps)

    ps =EXP(lnps)

    do k=1,nlev+1
       do j=1,np
          do i=1,np
             ph(i,j,k)   = hvcoord%hyai(k)*hvcoord%ps0 + hvcoord%hybi(k)*ps(i,j)
          end do
       end do
    end do


    ! We suppose here Neumann for q

    k=1
    do j=1,np
       do i=1,np
          dq(i,j,k) = 0.0D0
          dp(i,j,k) = (ph(i,j,k+1) - ph(i,j,k))
       end do
    end do

    k=nlev
    do j=1,np
       do i=1,np
          dq(i,j,k) = 0.0D0
          dp(i,j,k) = (ph(i,j,k+1) - ph(i,j,k))
       end do
    end do

    do k=2,nlev-1
       do j=1,np
          do i=1,np
             dq(i,j,k) = 0.5D0*(q(i,j,k+1) - q(i,j,k-1))
             dp(i,j,k) = (ph(i,j,k+1) - ph(i,j,k))
          end do
       end do
    enddo

    proj(:,:)=0.0D0

    do k=1,nlev
       do j=1,np
          do i=1,np
             proj(i,j) = proj(i,j) + v(i,j,k)*dq(i,j,k)
          enddo
       enddo
    enddo

  end subroutine ProjVandDQ



  ! PRIVATE

  subroutine Prim_Advance_Multicloud_Init()

    call initEdgeBuffer(edgeMc,2)
    call Diffusion_Multicloud_init()

  end subroutine Prim_Advance_Multicloud_Init

  subroutine diffusion_Multicloud_init()

       call initEdgeBuffer(edgeS1,2)
       call initEdgeBuffer(edgeS2,4)
       call initEdgeBuffer(edge3, 6)
       call initEdgeBuffer(edge4, 8)

  end subroutine diffusion_Multicloud_init

  !=============================================================
  !
  ! Advances in time the advection equation of the
  ! Galerkin projected primitive equns and advances the ODE for
  ! the bnd layer model
  !
  !=============================================================
  subroutine Prim_Advance_Multicloud(elem,elem_physics, cm,hvcoord,deriv,flt,hybrid,dt,tl,nets,nete)

    implicit none
    type (element_t), intent(inout)   :: elem(:)
    type (elem_physics_t), intent(inout)   :: elem_physics(:)
    type (ColumnModelMulticloud_t)    :: cm
    type (derivative_t), intent(in)   :: deriv
    type (hvcoord_t)                  :: hvcoord
    type (filter_t)                   :: flt

    type (hybrid_t),     intent(in)   :: hybrid

    real(kind=real_kind) , intent(in) :: dt
    type (TimeLevel_t)                :: tl

    integer,intent(in)                :: nets,nete

    real(kind=real_kind) :: dt2, qtens, tebtens
    real(kind=real_kind) :: rdx, rdy
    real(kind=real_kind) :: v1, v2, vcon1, vcon2
    real(kind=real_kind) :: ubar1,ubar2,q,Hd,Hd0,D0,D00,D,Pbig,deltem,deltem0,tem,tem0,teb
    real(kind=real_kind) :: lws,lws0,Qc,Qc0,ig, Qcongest,Qcongest0,Pbig0
    real(kind=real_kind) :: Acst,Bcst


    real (kind=real_kind), dimension(np,np,2) :: gradqmc, gradteb
    real (kind=real_kind), dimension(np,np,2) :: gv,gubar
    real (kind=real_kind), dimension(np,np)   :: div,divpert,divubar,rmetdet
    real (kind=real_kind), dimension(np,np)   :: t1,t2,t1_0,t2_0,ps,ps0
    real (kind=real_kind), dimension(np,np,nlev)   :: p,pot,pot0,p_0,zfull
    real (kind=real_kind), dimension(np,np,nets:nete) :: Hd2D

    real (kind=real_kind) :: cooling(nlev),an,anm1,tau_damp
    real (kind=real_kind) :: tp1,tp2,up1(2),up2(2)

    integer :: i,j,k,ie
    integer :: n0,nm1,np1,nfilt, n0v, nstep

    logical :: Debug=.FALSE.

    logical :: flag_deltem=.FALSE.

    n0    = tl%n0
    nm1   = tl%nm1
    np1   = tl%np1

    an = 1.5D0
    anm1=1.0D0-an

    nfilt = n0
    nstep = tl%nstep

    n0v = n0

    dt2 = 2.0_real_kind * dt
    ! Use real vel at k=nlev instead
    ! of ubar for teb ... later


#if 1
    if(nstep>0 .and. filter_freq_advection > 0 .and. modulo(nstep,filter_freq_advection)==0) then

       do ie=nets,nete
          call filter_P(elem_physics(ie)%qmc(:,:,nfilt),  flt)
          call filter_P(elem_physics(ie)%teb(:,:,nfilt),flt)
          do j=1,np
             do i=1,np
                elem_physics(ie)%qmc(i,j,nfilt) = elem_physics(ie)%mp(i,j)*elem_physics(ie)%qmc(i,j,nfilt)
                elem_physics(ie)%teb(i,j,nfilt) = elem_physics(ie)%mp(i,j)*elem_physics(ie)%teb(i,j,nfilt)
             end do
          end do
          call edgeVpack(edgeMc,elem_physics(ie)%qmc(:,:,nfilt),1,0,elem(ie)%desc)
          call edgeVpack(edgeMc,elem_physics(ie)%teb(:,:,nfilt),1,1,elem(ie)%desc)
       end do

       call bndry_exchangeV(hybrid,edgeMc)

       do ie=nets,nete
          call edgeVunpack(edgeMc,elem_physics(ie)%qmc(:,:,nfilt),1,0,elem(ie)%desc)
          call edgeVunpack(edgeMc,elem_physics(ie)%teb(:,:,nfilt),1,1,elem(ie)%desc)
          do j=1,np	
             do i=1,np
                elem_physics(ie)%qmc(i,j,nfilt) = elem_physics(ie)%rmp(i,j)*elem_physics(ie)%qmc(i,j,nfilt)
                elem_physics(ie)%teb(i,j,nfilt) = elem_physics(ie)%rmp(i,j)*elem_physics(ie)%teb(i,j,nfilt)
             enddo
          enddo
       end do

    end if
#endif

    do ie=nets,nete

       rdx=2.0_real_kind/(elem(ie)%dx*rearth) ! strong derivative inverse x length
       rdy=2.0_real_kind/(elem(ie)%dy*rearth) ! strong derivative inverse y length

       !
       ! Project the velocity onto the 2 barotropic modes selected and do an average Ubar
       !
       ! V1
       call VerticalProjectionVector(elem(ie)%state%v(:,:,:,:,n0v),elem_physics(ie)%uproj1(:,:,:,n0v),&
            cm%D%phi(:,1),elem(ie)%state%lnps(:,:,n0v),hvcoord)
       ! V2
       call VerticalProjectionVector(elem(ie)%state%v(:,:,:,:,n0v),elem_physics(ie)%uproj2(:,:,:,n0v),&
            cm%D%phi(:,2),elem(ie)%state%lnps(:,:,n0v),hvcoord)
       ! Ubar
       call VerticalProjectionVector(elem(ie)%state%v(:,:,:,:,n0v),elem_physics(ie)%ubar(:,:,:,n0v),&
            cm%D%cstmode(:),elem(ie)%state%lnps(:,:,n0v),hvcoord)

       ! 01/20/2009
       do j=1,np	
          do i=1,np
             up1 = elem_physics(ie)%uproj1(i,j,:,n0v)
             up2 = elem_physics(ie)%uproj2(i,j,:,n0v)
             elem_physics(ie)%uproj1(i,j,1,n0v)=up1(1)*cm%D%Aphi(1,1)+up2(1)*cm%D%Aphi(1,2)
             elem_physics(ie)%uproj1(i,j,2,n0v)=up1(2)*cm%D%Aphi(1,1)+up2(2)*cm%D%Aphi(1,2)

             elem_physics(ie)%uproj2(i,j,1,n0v)=up1(1)*cm%D%Aphi(2,1)+up2(1)*cm%D%Aphi(2,2)
             elem_physics(ie)%uproj2(i,j,2,n0v)=up1(2)*cm%D%Aphi(2,1)+up2(2)*cm%D%Aphi(2,2)
          enddo
       enddo

       do j=1,np	
          do i=1,np
             q  = elem_physics(ie)%qmc(i,j,n0v)
             v1 = (elem_physics(ie)%uproj1(i,j,1,n0v) + cm%alpha*elem_physics(ie)%uproj2(i,j,1,n0v) + 0.0D0*elem_physics(ie)%ubar(i,j,1,n0v))*q
             v2 = (elem_physics(ie)%uproj1(i,j,2,n0v) + cm%alpha*elem_physics(ie)%uproj2(i,j,2,n0v) + 0.0D0*elem_physics(ie)%ubar(i,j,2,n0v))*q

             ! project from sphere to contravariant velocities
             vcon1 = elem(ie)%Dinv(1,1,i,j)*v1 + elem(ie)%Dinv(1,2,i,j)*v2
             vcon2 = elem(ie)%Dinv(2,1,i,j)*v1 + elem(ie)%Dinv(2,2,i,j)*v2

             gv(i,j,1) = elem(ie)%metdet(i,j)*vcon1
             gv(i,j,2) = elem(ie)%metdet(i,j)*vcon2

          enddo
       enddo

       div(:,:) = divergence(gv,deriv)

       do j=1,np	
          do i=1,np
             v1 = elem_physics(ie)%uproj1(i,j,1,n0v)*cm%D%Qt1 + elem_physics(ie)%uproj2(i,j,1,n0v)*cm%D%Qt2
             v2 = elem_physics(ie)%uproj1(i,j,2,n0v)*cm%D%Qt1 + elem_physics(ie)%uproj2(i,j,2,n0v)*cm%D%Qt2

             ! project from sphere to contravariant velocities
             vcon1 = elem(ie)%Dinv(1,1,i,j)*v1 + elem(ie)%Dinv(1,2,i,j)*v2
             vcon2 = elem(ie)%Dinv(2,1,i,j)*v1 + elem(ie)%Dinv(2,2,i,j)*v2

             gv(i,j,1) = elem(ie)%metdet(i,j)*vcon1
             gv(i,j,2) = elem(ie)%metdet(i,j)*vcon2

             ! this was eliminated
             !v1 = elem_physics(ie)%ubar(i,j,1,n0v)*cm%Qt0
             !v2 = elem_physics(ie)%ubar(i,j,2,n0v)*cm%Qt0

             v1 = elem_physics(ie)%ubar(i,j,1,n0v)
             v2 = elem_physics(ie)%ubar(i,j,2,n0v)

             ! project from sphere to contravariant velocities
             vcon1 = elem(ie)%Dinv(1,1,i,j)*v1 + elem(ie)%Dinv(1,2,i,j)*v2
             vcon2 = elem(ie)%Dinv(2,1,i,j)*v1 + elem(ie)%Dinv(2,2,i,j)*v2

             ! 5/21/2009 was missing q
             q  = elem_physics(ie)%qmc(i,j,n0v)
             gubar(i,j,1) = elem(ie)%metdet(i,j)*vcon1*q
             gubar(i,j,2) = elem(ie)%metdet(i,j)*vcon2*q

          enddo
       enddo

       divpert(:,:)   = divergence(gv,deriv)
       divubar(:,:)   = divergence(gubar,deriv)
       gradqmc(:,:,:) = gradient(elem_physics(ie)%qmc(:,:,n0),  deriv)
       gradteb(:,:,:) = gradient(elem_physics(ie)%teb(:,:,n0),  deriv)

       rmetdet(:,:) = 1.0_real_kind/elem(ie)%metdet(:,:)

       ps(:,:)  = Getsurfpress(elem(ie)%state%lnps(:,:,nm1))
       ps0(:,:) = Getsurfpress(elem(ie)%state%lnps(:,:,n0))

       do k=1,nlev
          do j=1,np
             do i=1,np
                p(i,j,k)   = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*ps(i,j)
                p_0(i,j,k) = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*ps0(i,j)
             enddo
          enddo
       enddo

       ! ASC REMOVE BG T(eta) STATE
       !do k=1,nlev
          pot0(:,:,:) = Temp2PotTemp(elem(ie)%state%T(:,:,:,n0),p_0(:,:,:))-elem_physics(ie)%pot0(:,:,:)
          pot(:,:,:)  = Temp2PotTemp(elem(ie)%state%T(:,:,:,nm1),p(:,:,:))-elem_physics(ie)%pot0(:,:,:)
          if(Debug)print *,"Potential temp     at level ",k, " is ", sum(abs(pot(:,:,k)))
          if(Debug)print *,"Potential temp t=0 at level ",k, " is ", sum(abs(elem_physics(ie)%pot0(:,:,k)))

       !enddo

       call VerticalProjection(pot(:,:,:),t1(:,:),cm%D%psi(:,1),elem(ie)%state%lnps(:,:,nm1),hvcoord)
       call VerticalProjection(pot(:,:,:),t2(:,:),cm%D%psi(:,2),elem(ie)%state%lnps(:,:,nm1),hvcoord)
       call VerticalProjection(pot0(:,:,:),t1_0(:,:),cm%D%psi(:,1),elem(ie)%state%lnps(:,:,n0),hvcoord)
       call VerticalProjection(pot0(:,:,:),t2_0(:,:),cm%D%psi(:,2),elem(ie)%state%lnps(:,:,n0),hvcoord)

       ! new 01/20/2009
       do j=1,np	
          do i=1,np
             tp1 = t1(i,j)
             tp2 = t2(i,j)
             t1(i,j) = tp1*cm%D%Apsi(1,1) + tp2*cm%D%Apsi(1,2)
             t2(i,j) = tp1*cm%D%Apsi(2,1) + tp2*cm%D%Apsi(2,2)
             tp1 = t1_0(i,j)
             tp2 = t2_0(i,j)
             t1_0(i,j) = tp1*cm%D%Apsi(1,1) + tp2*cm%D%Apsi(1,2)
             t2_0(i,j) = tp1*cm%D%Apsi(2,1) + tp2*cm%D%Apsi(2,2)
          enddo
       enddo

       if(Debug)print *,"theta1 = ", sum(abs(t1))
       if(Debug)print *,"theta2 = ", sum(abs(t2))
       if(Debug)print *,"QMC = ",sum(abs(elem_physics(ie)%qmc(:,:,:)))

       do j=1,np	
          do i=1,np

             ubar1 = elem(ie)%state%v(i,j,1,nlev,n0v)
             ubar2 = elem(ie)%state%v(i,j,2,nlev,n0v)

             ! project from sphere to contravariant velocities
             vcon1 = elem(ie)%Dinv(1,1,i,j)*ubar1 + elem(ie)%Dinv(1,2,i,j)*ubar2
             vcon2 = elem(ie)%Dinv(2,1,i,j)*ubar1 + elem(ie)%Dinv(2,2,i,j)*ubar2

             !qtens   = vcon1*gradqmc(i,j,1)   + vcon2*gradqmc(i,j,2)
             qtens = 0.0D0

             tebtens = vcon1*gradteb(i,j,1)   + vcon2*gradteb(i,j,2) ! ASC RHS To be defined by bnd layer model
	    ! tebtens=0.d0
             tem     = elem_physics(ie)%qmc(i,j,nm1) + cm%D%psi1avg*(t1(i,j)+ cm%alpha2*t2(i,j)) !*cm%psi2avg)

             deltem = elem_physics(ie)%teb(i,j,nm1)-tem+cm%TebMinTem
             teb    = elem_physics(ie)%teb(i,j,n0)
             if(Debug)print *,"teb       = ", teb
             if(Debug)print *,"tem       = ", tem
             if(Debug)print *,"teb - tem = ", deltem

             !if(.not.flag_deltem)then
             !   if(DSIGN(1.0D0,deltem) < 0.5D0)flag_deltem=.TRUE.
             !end if

             !deltem = max(-deltem,deltem)

             lws =  LambdaSwitch(deltem,cm%Tminus,cm%Tplus,cm%A,cm%B,cm%lambdastar)
             if(Debug)print *,"LWS = ",lws

             ! At "n-1"
             Qc   = max(0.0D0,cm%D%Q0c + (1.0D0/cm%Tau_convec)*(cm%a1*elem_physics(ie)%teb(i,j,nm1) &
                  + cm%a2*elem_physics(ie)%qmc(i,j,nm1) -cm%a0*(t1(i,j)+ cm%gamma2*t2(i,j))))
             Hd   = ((1.0D0 -lws)/(1.0D0 - cm%lambdastar))*Qc
             D0   = (cm%D%m0/cm%Q0R1)*max(0.0D0,(cm%Q0R1 + cm%mu*(elem_physics(ie)%Hs(i,j,nm1) &
                  - elem_physics(ie)%Hc(i,j,nm1))))*deltem

             ! At "n"
             tem0    = elem_physics(ie)%qmc(i,j,n0) + cm%D%psi1avg*(t1_0(i,j)+ cm%alpha2*t2_0(i,j)) !*cm%psi2avg)
             deltem0 = elem_physics(ie)%teb(i,j,n0)-tem0+cm%TebMinTem

            ! if(.not.flag_deltem)then
            !    if(DSIGN(1.0D0,deltem0) < 0.5D0)flag_deltem=.TRUE.
            ! end if

             !deltem0 = max(-deltem0,deltem0)


             lws0   =  LambdaSwitch(deltem0,cm%Tminus,cm%Tplus,cm%A,cm%B,cm%lambdastar)
             Qc0    = max(0.0D0,cm%D%Q0c + (1.0D0/cm%Tau_convec)*(cm%a1*elem_physics(ie)%teb(i,j,n0) &
                  + cm%a2*elem_physics(ie)%qmc(i,j,n0) -cm%a0*(t1_0(i,j)+ cm%gamma2*t2_0(i,j)))) 
             Hd0    = ((1.0D0 -lws0)/(1.0D0 - cm%lambdastar))*Qc0
             D00    = (cm%D%m0/cm%Q0R1)*max(0.0D0,(cm%Q0R1 + cm%mu*(elem_physics(ie)%Hs(i,j,n0) &
                  - elem_physics(ie)%Hc(i,j,n0))))*deltem0

             if(Debug)print *,"D0  = ",D0

             if(cm%closure == 1)then
               D        = lws*D0
               Qcongest = D/cm%Htall
               Qcongest0 = lws0*D00/cm%Htall
               D0 =lws0*D00
             else
                Qcongest = max(0.0D0,cm%D%Q0c + (1.0D0/cm%Tau_convec)*(elem_physics(ie)%teb(i,j,nm1) &
                   -cm%a0prime*(t1(i,j)+ cm%gamma2prime*t2(i,j))))
                Qcongest0 = max(0.0D0,cm%D%Q0c + (1.0D0/cm%Tau_convec)*(elem_physics(ie)%teb(i,j,n0) &
                   -cm%a0prime*(t1_0(i,j)+ cm%gamma2prime*t2_0(i,j))))
                D    = D0
                D0=D00
             end if

             ! 7/1/09 we should remove max
             !Pbig = max(0.0D0,Hd*cm%psi1avgtrunc + (elem_physics(ie)%Hc(i,j,nm1) - elem_physics(ie)%Hs(i,j,nm1))*cm%psi2avgtrunc) ! Precipitation
             !Pbig0 = max(0.0D0,Hd0*cm%psi1avgtrunc + (elem_physics(ie)%Hc(i,j,n0) - elem_physics(ie)%Hs(i,j,n0))*cm%psi2avgtrunc) ! Precipitation
             Pbig  = Hd*cm%D%psi1avgtrunc  + (1.0D0-cm%csr)*(elem_physics(ie)%Hc(i,j,nm1) - elem_physics(ie)%Hs(i,j,nm1))*cm%D%psi2avgtrunc ! Precipitation
             Pbig0 = Hd0*cm%D%psi1avgtrunc + (1.0D0-cm%csr)*(elem_physics(ie)%Hc(i,j,n0) - elem_physics(ie)%Hs(i,j,n0))*cm%D%psi2avgtrunc ! Precipitation


             ! just to store the info
             ! (Pmc is one dt back)
             elem_physics(ie)%Pmc(i,j) = Pbig0
             elem_physics(ie)%Hd(i,j)  = Hd0

             if(Hd0 < 0.0D0)then
                print *,"Hd0                = ", Hd0
                print *,"Lambdaswitch       = ", lws0
                print *,"Qc0                = ", Qc0
                print *,"Facteur de Hd0     = ",((1.0D0 -lws0)/(1.0D0 - cm%lambdastar))
                print *,"Tminus             = ", cm%Tminus
                print *,"Tplus              = ", cm%Tplus
                print *,"A                  = ", cm%A
                print *,"B                  = ", cm%B
                print *,"lambdastar         = ", cm%lambdastar
             end if

             if(cm%closure == 1)then
                elem_physics(ie)%D(i,j)   = lws*D0/cm%Htall
             else
                elem_physics(ie)%D(i,j)   = D0/cm%Htall
             end if

             Hd2D(i,j,ie)            = Hd0

             ! Avg moisture
             ! 5/21/2009
             ! 1) if ubar contrib in div set to zero THEN divubar  needs to be there
             ! 2) if ubar contrib in div not zero    THEN divbar = 0
             ! we try option 1 to see if the perturbations were swamped by the bar quantities

             elem_physics(ie)%qmc(i,j,np1) = elem_physics(ie)%mp(i,j)*(elem_physics(ie)%qmc(i,j,nm1) &
                  - dt2*(div(i,j)+divubar(i,j)+divpert(i,j))*rmetdet(i,j) &
                  + dt2*(-Pbig0 + D0/cm%Htall)*elem_physics(ie)%mask(i,j))             
             
             ! Boundary layer potential temperature
             
             !elem_physics(ie)%teb(i,j,np1) = elem_physics(ie)%mp(i,j)*(elem_physics(ie)%teb(i,j,nm1) - dt2*(tebtens) &
             !     + dt2*((cm%TstarMinTeb-elem_physics(ie)%teb(i,j,n0))*cm%D%Tau_evap_Inv -D0/cm%h)*elem(ie)%mask(i,j))
             
             ! added the varying Sea surface temperature:
             !elem_physics(ie)%teb(i,j,np1) = elem_physics(ie)%mp(i,j)*(elem_physics(ie)%teb(i,j,nm1) - dt2*(tebtens) &
              !    + dt2*((elem_physics(ie)%delthetasurf(i,j)-elem_physics(ie)%teb(i,j,n0))*cm%D%Tau_evap_Inv -D0/cm%h)*elem(ie)%mask(i,j))
          

             !semi-analytic integration of theta_eb equation, to prevent stiffness: Oct 12, 2009
            
              Bcst    = (cm%D%m0/cm%Q0R1)*max(0.0D0,(cm%Q0R1 + cm%mu*(elem_physics(ie)%Hs(i,j,n0) &
                  - elem_physics(ie)%Hc(i,j,n0))))/cm%h
            
             Acst = cm%D%Tau_evap_Inv + Bcst  
             tem0 = tem0 - cm%TebMinTem
 
             elem_physics(ie)%teb(i,j,np1) = elem_physics(ie)%mp(i,j)*((elem_physics(ie)%teb(i,j,nm1)*dexp(-Acst*dt2) +&
                                           (elem_physics(ie)%delthetasurf(i,j)*cm%D%Tau_evap_Inv + Bcst*tem0)*(1.d0 - dexp(-Acst*dt2))/Acst)&
                    *elem_physics(ie)%mask(i,j) - dt2*tebtens)

                 
             
             ! ODE 1
             elem_physics(ie)%Hs(i,j,np1)  = (elem_physics(ie)%Hs(i,j,n0) &
                  +   an*dt*(cm%alpha_strat*Hd0 -elem_physics(ie)%Hs(i,j,n0))*elem_physics(ie)%mask(i,j)/cm%Tau_strat &
                  + anm1*dt*(cm%alpha_strat*Hd -elem_physics(ie)%Hs(i,j,nm1))*elem_physics(ie)%mask(i,j)/cm%Tau_strat)

             ! ODE 2
             elem_physics(ie)%Hc(i,j,np1)  = (elem_physics(ie)%Hc(i,j,n0) &
                  +  an*dt*(Qcongest0*cm%alpha_congest*(lws0-cm%lambdastar)/(1.0D0-cm%lambdastar)&
                  -elem_physics(ie)%Hc(i,j,n0))*elem_physics(ie)%mask(i,j)/cm%Tau_congest &
                  +anm1*dt*(Qcongest*cm%alpha_congest*(lws-cm%lambdastar)/(1.0D0-cm%lambdastar)&
                  -elem_physics(ie)%Hc(i,j,nm1))*elem_physics(ie)%mask(i,j)/cm%Tau_congest)
          enddo
       enddo

       if(Debug)print *,"qmc = ",sum(abs(elem_physics(ie)%qmc(:,:,:)))
       if(Debug)print *,"teb = ",sum(abs(elem_physics(ie)%teb(:,:,:)))
       if(Debug)print *,"Hs  = ",sum(abs(elem_physics(ie)%Hs(:,:,:)))
       if(Debug)print *,"Hc  = ",sum(abs(elem_physics(ie)%Hc(:,:,:)))
       if(Debug)print *,"t1  = ",sum(abs(t1(:,:)))
       if(Debug)print *,"t2  = ",sum(abs(t2(:,:)))

       call edgeVpack(edgeMc,elem_physics(ie)%qmc(:,:,np1),1,0,elem(ie)%desc)
       call edgeVpack(edgeMc,elem_physics(ie)%teb(:,:,np1),1,1,elem(ie)%desc)

    end do

    call bndry_exchangeV(hybrid,edgeMc)

    do ie=nets,nete
       call edgeVunpack(edgeMc,elem_physics(ie)%qmc(:,:,np1),1,0,elem(ie)%desc)
       call edgeVunpack(edgeMc,elem_physics(ie)%teb(:,:,np1),1,1,elem(ie)%desc)

       do j=1,np	
          do i=1,np
             elem_physics(ie)%qmc(i,j,np1) = elem_physics(ie)%rmp(i,j)*elem_physics(ie)%qmc(i,j,np1)
             elem_physics(ie)%teb(i,j,np1) = elem_physics(ie)%rmp(i,j)*elem_physics(ie)%teb(i,j,np1)
          enddo
       enddo
    enddo

    ! diffuse
    !call multicloud_diffusion(elem, nets,nete,np1,deriv,dt2,hybrid)

    cooling=presc_cooling_mc(hvcoord,cm,np,nlev,0) ! cooling can change in time

    do ie=nets,nete

       ps(:,:) = Getsurfpress(elem(ie)%state%lnps(:,:,n0))

       do j=1,np	
          do i=1,np
             ! clip
             !elem_physics(ie)%qmc(i,j,np1)  = max(1.e-12_real_kind, elem_physics(ie)%qmc(i,j,np1))

             ! Robert-Asselin filter
             elem_physics(ie)%qmc(i,j,n0)  = elem_physics(ie)%qmc(i,j,n0) + &
                  smooth*(elem_physics(ie)%qmc(i,j,nm1) - &
                  2.0D0*elem_physics(ie)%qmc(i,j,n0) + elem_physics(ie)%qmc(i,j,np1))

             elem_physics(ie)%teb(i,j,n0)  = elem_physics(ie)%teb(i,j,n0) + &
                  smooth*(elem_physics(ie)%teb(i,j,nm1) - &
                  2.0D0*elem_physics(ie)%teb(i,j,n0) + elem_physics(ie)%teb(i,j,np1))
          enddo
       enddo

       ! Store Qheating at n for the forcing to take place
       ! at n-1 on the next Dycore step

       ig = 1.0D0/g

       do k=1,nlev
          do j=1,np
             do i=1,np
                p(i,j,k) = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*ps(i,j)
             enddo
          enddo
       enddo

       do k=1,nlev
          do j=1,np	
             do i=1,np
                elem_physics(ie)%QHeating(i,j,k) =    (1.4D0*((Hd2D(i,j,ie)*cm%D%psitrunc(k,1) &
                     + (elem_physics(ie)%Hc(i,j,n0) - elem_physics(ie)%Hs(i,j,n0))*(cm%D%psitrunc(k,2) - cm%csr*cm%D%psi2avgtrunc)) - cooling(k))*elem_physics(ie)%mask(i,j) &
                     +  cm%Tau_r*elem_physics(ie)%pot0(i,j,k))*(p(i,j,k)/p0)**(kappa) &
                     -  cm%Tau_r*elem(ie)%state%T(i,j,k,n0)
             enddo
          enddo
       enddo

       do k=1,nlev
          if(Debug) print *," LEVEL = ",k
          if(Debug) print *,"Tau_r*pot0*conversion  = ",sum(abs(cm%Tau_r*elem_physics(ie)%pot0(:,:,k)*(p(:,:,k)/ps(:,:))**(kappa)))
          if(Debug) print *,"Tau_r*T                = ",sum(abs(cm%Tau_r*elem(ie)%state%T(:,:,k,n0)))
          if(Debug) print *,"COOLING                = ", cooling(k)
          if(Debug) print *,"Convective heating     = ",sum(abs((Hd2D(:,:,ie)*cm%D%psitrunc(k,1) &
               + (elem_physics(ie)%Hc(:,:,n0) - elem_physics(ie)%Hs(:,:,n0))*cm%D%psitrunc(k,2))))/REAL(np*np,kind=real_kind)
       end do

      if(Debug)print *,"QHeating  = ",sum(abs(elem_physics(ie)%QHeating(:,:,:)))

    end do

    ! ERROR DIAGNOSTIC

  !  if(flag_deltem)then
  !     print *,"DELTA THETA_em was negative at at least one grid point"
  !  end if

	
  end subroutine Prim_Advance_Multicloud



  subroutine multicloud_diffusion(elem, elem_physics, nets,nete,np1,deriv,dt2,hybrid)

    use hybrid_mod, only : hybrid_t
    use physical_constants, only : rearth
    use element_mod, only : element_t
    use bndry_mod, only : bndry_exchangev
    use control_mod, only : nu_mc
    use derivative_mod, only : gradient
    implicit none
    type(element_t), intent(inout), target :: elem(:)
    type(elem_physics_t), intent(inout), target :: elem_physics(:)
    integer, intent(in)                :: nets,nete
    integer, intent(in)                :: np1
    type(derivative_t), intent(in)     :: deriv
    real(kind=real_kind),intent(in)    :: dt2
    type(hybrid_t),intent(in)          :: hybrid

    ! ======================
    ! Local Variables
    ! ======================
    real(kind=real_kind) :: temp(np,np)

    real(kind=real_kind) :: grad_qmc_np1(np,np,2,nets:nete)
    real(kind=real_kind) :: lap_qmc_np1(np,np,nets:nete)

    real(kind=real_kind) :: gradqmc_tmp(np,np,2)
    real(kind=real_kind) :: divqmc_tmp(np,np)

    real(kind=real_kind) :: grad_teb_np1(np,np,2,nets:nete)
    real(kind=real_kind) :: lap_teb_np1(np,np,nets:nete)

    real(kind=real_kind) :: gradteb_tmp(np,np,2)
    real(kind=real_kind) :: divteb_tmp(np,np)

    real(kind=real_kind) :: v1,v2
    real(kind=real_kind) :: rdx,rdy

    real (kind=real_kind), dimension(np,np)      :: rmetdetv

    real(kind=real_kind), dimension(:,:), pointer :: mv
    real(kind=real_kind), dimension(:,:), pointer :: rmp
    real(kind=real_kind), dimension(:,:), pointer :: metdet
    real(kind=real_kind), dimension(:,:,:,:), pointer :: met
    real(kind=real_kind), dimension(:,:,:,:), pointer :: metinv
    real(kind=real_kind), dimension(:,:,:,:), pointer :: D
    real(kind=real_kind), dimension(:,:,:,:), pointer :: Dinv

    integer :: i,j,ie

    do ie=nets,nete

       rdx=2.0_real_kind/(elem(ie)%dx*rearth) ! strong derivative inverse x length
       rdy=2.0_real_kind/(elem(ie)%dy*rearth) ! strong derivative inverse y length

       mv       => elem_physics(ie)%mp
       metinv   => elem(ie)%metinv
       metdet   => elem(ie)%metdet

       gradqmc_tmp(:,:,:) = gradient(elem_physics(ie)%qmc(:,:,np1),deriv)
       gradteb_tmp(:,:,:) = gradient(elem_physics(ie)%teb(:,:,np1),deriv)

       do j=1,np
          do i=1,np
             v1 = mv(i,j)*gradqmc_tmp(i,j,1)
             v2 = mv(i,j)*gradqmc_tmp(i,j,2)
             grad_qmc_np1(i,j,1,ie) = metdet(i,j)*(metinv(1,1,i,j)*v1 + metinv(1,2,i,j)*v2)
             grad_qmc_np1(i,j,2,ie) = metdet(i,j)*(metinv(2,1,i,j)*v1 + metinv(2,2,i,j)*v2)

             v1 = mv(i,j)*gradteb_tmp(i,j,1)
             v2 = mv(i,j)*gradteb_tmp(i,j,2)
             grad_teb_np1(i,j,1,ie) = metdet(i,j)*(metinv(1,1,i,j)*v1 + metinv(1,2,i,j)*v2)
             grad_teb_np1(i,j,2,ie) = metdet(i,j)*(metinv(2,1,i,j)*v1 + metinv(2,2,i,j)*v2)
          end do
       end do

       call edgeVpack(edgeS2,grad_qmc_np1(:,:,:,ie),2,0,elem(ie)%desc)
       call edgeVpack(edgeS2,grad_teb_np1(:,:,:,ie),2,2,elem(ie)%desc)
       call edgerotate(edgeS2,4,0,elem(ie)%desc)

    end do

    call bndry_exchangeV(hybrid,edgeS2)

    do ie=nets,nete
       rdx=2.0_real_kind/(elem(ie)%dx*rearth) ! strong derivative inverse x length
       rdy=2.0_real_kind/(elem(ie)%dy*rearth) ! strong derivative inverse y length

       mv       => elem_physics(ie)%mp
       rmp      => elem_physics(ie)%rmp
       metdet   => elem(ie)%metdet
       metinv   => elem(ie)%metinv
       rmetdetv(:,:)=1.0_real_kind/elem(ie)%metdet(:,:)
       D        => elem(ie)%D
       Dinv     => elem(ie)%Dinv

       call edgeVunpack(edgeS2, grad_qmc_np1(:,:,:,ie), 2, 0, elem(ie)%desc)
       call edgeVunpack(edgeS2, grad_teb_np1(:,:,:,ie), 2, 2, elem(ie)%desc)

       do j=1,np
          do i=1,np
             grad_qmc_np1(i,j,1,ie) = rmp(i,j)*grad_qmc_np1(i,j,1,ie)
             grad_qmc_np1(i,j,2,ie) = rmp(i,j)*grad_qmc_np1(i,j,2,ie)
             grad_teb_np1(i,j,1,ie) = rmp(i,j)*grad_teb_np1(i,j,1,ie)
             grad_teb_np1(i,j,2,ie) = rmp(i,j)*grad_teb_np1(i,j,2,ie)
          end do
       end do

       divqmc_tmp(:,:) = divergence(grad_qmc_np1(:,:,:,ie),deriv)
       divteb_tmp(:,:) = divergence(grad_teb_np1(:,:,:,ie),deriv)

       do j=1,np
          do i=1,np
             lap_qmc_np1(i,j,ie)   = mv(i,j)*rmetdetv(i,j)*divqmc_tmp(i,j)
             lap_teb_np1(i,j,ie)   = mv(i,j)*rmetdetv(i,j)*divteb_tmp(i,j)
          end do
       end do

       call edgeVpack(edgeS1, lap_qmc_np1(1,1,ie),1, 0,elem(ie)%desc)
       call edgeVpack(edgeS1, lap_teb_np1(1,1,ie),1, 1,elem(ie)%desc)

    end do

    call bndry_exchangeV(hybrid,edgeS1)

    do ie=nets,nete

       rmp           => elem_physics(ie)%rmp
       metdet        => elem(ie)%metdet
       rmetdetv(:,:) = 1.0_real_kind/elem(ie)%metdet(:,:)


       call edgeVunpack(edgeS1, lap_qmc_np1(1,1,ie), 1, 0, elem(ie)%desc)
       call edgeVunpack(edgeS1, lap_teb_np1(1,1,ie), 1, 0, elem(ie)%desc)

       do j=1,np
          do i=1,np
             lap_qmc_np1(i,j,ie)   = rmp(i,j)*lap_qmc_np1(i,j,ie)
             lap_teb_np1(i,j,ie)   = rmp(i,j)*lap_teb_np1(i,j,ie)
             elem_physics(ie)%qmc(i,j,np1)   = elem_physics(ie)%qmc(i,j,np1)   + nu_mc*dt2*lap_qmc_np1(i,j,ie)
             elem_physics(ie)%teb(i,j,np1)   = elem_physics(ie)%teb(i,j,np1)   + nu_mc*dt2*lap_teb_np1(i,j,ie)
          end do
       end do

    end do

  end subroutine multicloud_diffusion

  subroutine InitShearDamping(elem,hybrid,nets,nete)
    use parallel_mod, only : parallel_t
    type (parallel_t)                    :: par
    type (element_t), intent(inout)   :: elem(:)
    type (hybrid_t),     intent(in)   :: hybrid
    integer,intent(in)                :: nets,nete

    integer :: nlat,nlon

    allocate(interpdata(nelemd))

    call setup_latlon_interp(elem,interpdata, hybrid%par)

    nlat = get_interp_parameter('nlat')
    nlon = get_interp_parameter('nlon')

    allocate(lat(nlat))

    lat = get_interp_lat()

    ! space to compute the interpolation
    allocate(udamp(nlat))
    allocate(udamp_copy(nlat))

  end subroutine InitShearDamping

  !
  ! Computes a damping for the shear component of the flow
  !
  subroutine ApplyShearDamping(elem,elem_physics, tau_damp,hybrid,tl,nets,nete)

#ifdef _MPI
    use parallel_mod, only: mpi_sum, mpireal_t
#endif

    type (element_t),    intent(inout):: elem(:)
    type (elem_physics_t),    intent(inout):: elem_physics(:)
    real(kind=real_kind),intent(in)   :: tau_damp
    type (hybrid_t),     intent(in)   :: hybrid
    type (TimeLevel_t)                :: tl
    integer,             intent(in)   :: nets,nete

    ! local

    integer :: ie,i,j,k,ilat,ierr,nlat
    integer :: na,nb,nm,ii
    integer :: ncnt
    integer :: st,en ! start-end in array
    real(kind=real_kind), allocatable :: var(:)
    real(kind=real_kind) :: xp,dlat

    ncnt = sum(interpdata(nets:nete)%n_interp)

    allocate(var(ncnt))

    udamp(:)=0.0D0

    ! Perform interpolation to lat-lon grid

    nlat = get_interp_parameter('nlat')

    dlat = 1.0D0/real(nlat,kind=real_kind)

    st=1
    do ie=nets,nete
       en=st+interpdata(ie)%n_interp
       call interpolate_scalar(interpdata(ie),elem_physics(ie)%ubar(:,:,1,tl%nm1),np,var(st:en))
       do i=1,interpdata(ie)%n_interp
          ilat = interpdata(ie)%ilat(i)
          udamp(ilat) = udamp(ilat) + var(st-1+i)*dlat
       end do
       st=st+interpdata(ie)%n_interp
    end do

    ! share 1D function

#ifdef _MPI
    do i=1,nlat
       udamp_copy(i) = udamp(i)
       udamp(i)      = 0.0D0
    end do
    call MPI_Allreduce(udamp_copy(1),udamp(1),nlat,MPIreal_t,MPI_SUM,hybrid%par%comm,ierr)
#endif

    ! linear interpolation back to cubed sphere

    do ie=nets,nete
       do j=1,np
          do i=1,np

             xp = elem(ie)%spherep(i,j)%lat
             na = 1
             nb = nlat

             do
                if  ((nb-na) <=  1)  exit
                nm = (nb + na)/2
                if (xp  >  lat(nm)) then
                   na = nm
                else
                   nb = nm
                endif
             enddo

             ii   = na
             dlat = 1.0D0/(lat(ii+1)-lat(ii))

             do k=1,nlev

                elem(ie)%derived%FM(i,j,1,k,tl%nm1) = elem(ie)%derived%FM(i,j,1,k,tl%nm1)  &
                     - tau_damp*(udamp(ii+1)*(xp - lat(ii)) - udamp(ii)*(xp - lat(ii+1)))*dlat 
             end do
          end do
       end do

       !print *,elem(ie)%state%FM(:,:,1,1,tl%nm1)

    end do
    !print *,"TAU =", tau_damp
    !print *,"DAMPING = ",udamp



   deallocate(var)

  end subroutine ApplyShearDamping

end module multicloud_mod
