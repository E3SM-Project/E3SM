module gcrk_mod                                                                        
  use grid_init_mod, only : dxi, dzi, dt, dti, gac

  use prof_init_mod, only : rho0

  implicit none
  real :: eer, eem
  integer :: niter, nitsm, icount

  integer :: itmn, iprc
  
  real, allocatable :: r(:), qr(:), ar(:) 
  real, private :: dxil, dzil
contains                                                         
  subroutine gcrk_1(p,pfx,pfz,u,w,itr,eps0, nx, nz, lord, nplt)
    implicit none
    integer, intent(in) :: nx, nz, lord, nplt
    integer, intent(inout) :: itr
    REAL, intent(in) :: u(nx*nz), w(nx*nz) 
    real, intent(out) :: pfx(nx*nz), pfz(nx*nz)
    real , intent(inout) :: p(nx*nz)


    real::  x (nx*nz, lord), ax (nx*nz, lord), ax2 (lord), axar (lord),   &
         del (lord), rho2d (nx, nz) , eps0
    real :: beta, dvmx, rl2,  rax, plmx,  eps, eer0, rl20, eem0
    real :: snorm, epa
    integer :: it, ii,  i1, ier, i, k, ner, nlc, nxz

!convergence test modes **************************************************
    LOGICAL, parameter ::  ctest=.false.
    real :: err(0:nplt),xitr(0:nplt)       

    if(ctest) then                       
       itr=6000/lord                        
       ner=60                               
       snorm=1./real(nx*nz)                   
       eps0=1.e-15                          
    endif
!convergence test modes **************************************************
#ifdef TESTMODE
      call tbeg('gcrk')
#endif
      nxz = nx*nz
      do k=1,nz
      do i=1,nx
      rho2d(i,k)=rho0(k)*gac(k)
      enddo
      enddo

      eps=eps0*dti
      epa=1.e-30
      nlc=0

      call precon_1(r,qr,ar,pfx,pfz,rho2d,iprc,0,nx,nz)

      do k=1,nxz
        r(k)=0.
       ar(k)=0.
       qr(k)=0.
      enddo
      do i=1,lord
       do k=1,nxz
         x(k,i)=0.
        ax(k,i)=0.
       enddo
      enddo
      call prforc_1(p,pfx,pfz,u,w,nx,nz)
      call rhsdiv_1(pfx,pfz,rho2d,r,-1, nx, nz)
      call precon_1(r,qr,ar,pfx,pfz,rho2d,iprc,1,nx,nz)
      eer0=0.
      eem0=-1.e15
      rl20=0.
      do k=1,nxz
      eer0=eer0+qr(k)**2
      eem0=max(eem0,abs(qr(k)))
      rl20=rl20+r(k)**2
      enddo
      eer0=max(eer0,epa)
      eem0=max(eem0,epa)
      rl20=max(rl20,epa)
!convergence test modes **************************************************
      if(ctest) then                                                    
      do ier=0,nplt                                                     
      err(ier)=eps                                                      
      enddo                                                             
      eer=-1.e15                                                        
      do k=1,nxz                                                         
        eer=max(eer,abs(r(k)))                                          
      end do                                                            
      err(0)=eer                                                        
      print 300,  err(0)                                                
  300 format(4x,e11.4,' residual error at it=1')                        
      endif                                                             
!convergence test modes **************************************************
       do k=1,nxz
        x(k,1)=qr(k)
       enddo
      call laplc_1(x(1,1),ax(1,1),pfx,pfz,nx,nz)

      do 100 it=1,itr
       do i=1,lord
        ax2(i)=0.
        rax=0.
         do k=1,nxz
          rax=rax+r(k)*ax(k,i)
          ax2(i)=ax2(i)+ax(k,i)*ax(k,i)
         enddo
        ax2(i)=max(epa,ax2(i))
        beta=-rax/ax2(i)
        dvmx=-1.e15
        rl2=0.
         do k=1,nxz
          p(k)=p(k)+beta* x(k,i)
          r(k)=r(k)+beta*ax(k,i)
          dvmx=max(dvmx,abs(r(k)))
          rl2=rl2+r(k)*r(k)
         enddo
       if(dvmx.le.eps.and.it.ge.itmn) go to 200
       if(rl2.ge.rl20.and..not.ctest) go to 200
          rl20=max(rl2,epa)
       call precon_1(r,qr,ar,pfx,pfz,rho2d,iprc,1,nx,nz)
       call laplc_1(qr,ar,pfx,pfz,nx,nz)
        nlc=nlc+1
         do ii=1,i
          axar(ii)=0.
           do k=1,nxz
            axar(ii)=axar(ii)+ax(k,ii)*ar(k)
           enddo
          del(ii)=-axar(ii)/ax2(ii)
         enddo
        if(i.lt.lord) then
          do k=1,nxz
            x(k,i+1)=qr(k)
           ax(k,i+1)=ar(k)
          enddo
           do ii=1,i
            do k=1,nxz
              x(k,i+1)= x(k,i+1)+del(ii)* x(k,ii)
             ax(k,i+1)=ax(k,i+1)+del(ii)*ax(k,ii)
            enddo
           enddo
        else
          do k=1,nxz
            x(k,1)=qr(k)+del(1)* x(k,1)
           ax(k,1)=ar(k)+del(1)*ax(k,1)
          enddo
           do ii=2,i
            do k=1,nxz
              x(k,1 )= x(k,1)+del(ii)* x(k,ii)
              x(k,ii)=0.
             ax(k,1 )=ax(k,1)+del(ii)*ax(k,ii)
             ax(k,ii)=0.
            enddo
           enddo
        endif
!convergence test modes **************************************************
      if(ctest) then                                                    
      if(nlc/ner*ner.eq.nlc) then                                       
      ier=nlc/ner                                                       
      eer=-1.e15                                                        
      do 50 k=1,nxz                                                      
   50 eer=max(eer,abs(r(k)))                                            
      err(ier)=eer                                                      
      endif                                                             
      endif                                                             
!convergence test modes **************************************************
       enddo
  100 continue
  200 continue
      eer=0.
      eem=-1.e15
      do k=1,nxz
      eer=eer+qr(k)**2
      eem=max(eem,abs(qr(k)))
      enddo
      eer=eer/eer0
      eem=eem/eem0
      niter=nlc
      nitsm=nitsm+niter
      icount=icount+1

!convergence test modes **************************************************
!     if(ctest) then                                                    * 
!     print 301, (err(ier),ier=1,nplt,1)                                *
! 301 format(4x,5e11.4)                                                 *
!     do 400 ier=0,nplt                                                 *
!     xitr(ier)=ier*ner                                                 *
! 400 err(ier)=alog10(err(ier)*dt )                                     *
!     plmx=real(itr*lord)                                               *
!     call set(.1,.9,.1,.9,0.,plmx,-10.,0.,1)                           *
!     call labmod('(i4)','(f5.0)',4,4,2,2,20,20,0)                      *
!     call periml(4,10,5,2)                                             *
!     call dashdc('$$$$$$$$$$$$$$$$$$$$',10,12)                         *
!     call curved(xitr,err,nplt+1)                                      *
!     i1=int(102.4+409.6)                                               *
!     call wtstr(cpux(i1),cpuy(50),'niter',3,0,0)                       *
!     call wtstr(cpux(17),cpuy(i1),'log e',3,90,0)                      *
!     call frame                                                        *
!     endif                                                             *
!convergence test modes **************************************************
#ifdef TESTMODE
      call tend('gcrk')
#endif
      return
    end subroutine gcrk_1


    subroutine precon_1(rhs,p,r,c11,c33,d,iflg,jfl, nx, nz)
      implicit none
      integer, intent(in) :: nx, nz

      real, intent(in) :: rhs(nx,nz), d(nx,nz)
      real, intent(out) ::  p(nx,nz), c11(nx,nz), c33(nx,nz), r(nx,nz)
      integer, intent(in) :: iflg, jfl
      real:: e (nx, 0:nz-1), f (nx, 0:nz-1), px(nx,nz), dgh(nx,nz),  &
           po(nx,nz)                                                         
      real :: beta=-1.e15, beti, dn
      integer, parameter :: itr=2, line=1
      integer :: k, it, i
    !    DATA beta / - 1.e15 / 
    !    DATA itr, line / 2, 1 / 
      real, parameter :: omg=0.7, oms=1.-omg
      real :: dxi2, dzi2

!BGL  arrays for vector routines
      real rec_dn(nx), rec_1me(nx)
      integer :: nxz

      dxi2=0.25*dxi*dxi
      dzi2=0.25*dzi*dzi 

      nxz=nx*nz
      if(iflg.eq.0) then
       do i=1,nxz
        p(i,1)=rhs(i,1)
       enddo
      return
      endif
#ifdef TESTMODE
      call tbeg('  precon')
#endif
!      omg=.7
!      oms=1.-omg
!      dxi2=0.25*dxi*dxi
!      dzi2=0.25*dzi*dzi
      do k=1,nz
        do i=1,nx
         c33(i,k)=d(i,k)*dzi2
         c11(i,k)=d(i,k)*dxi2
         dgh(i,k)=0.
          po(i,k)=0.
           p(i,k)=0.
           r(i,k)=0.
        enddo

        if(line.eq.1) then
          dgh(1,k)=c11(2,k)+c11(nx-1,k)
          do i=2,nx-1
            dgh(i,k)=c11(i+1,k)+c11(i-1,k)
          enddo
          dgh(nx,k)=c11(2,k)+c11(nx-1,k)
        endif
      enddo

      if(jfl.eq.0) then
      if(line.eq.0) then
      beta=-1.e15
      do i=1,nxz
      beta=amax1(beta,abs(c11(i,1))/d(i,1))
      enddo
      beta=0.5/beta
      else
      beta=1.
      endif
#ifdef TESTMODE
      call tend('  precon')
#endif
      return
      endif
      beti=1./beta*(1-line)

      do 100 it=1,itr
      do i=1,nx
        rec_dn(i) = d(i,1)*beti+2.*c33(i,2)+dgh(i,1)
      end do
      call vrec(rec_dn, rec_dn, nx)
      do i=1,nx
        r(i,1)=r(i,1)+d(i,1)*(beti*p(i,1)-rhs(i,1)) &
             +dgh(i,1)*p(i,1)
      enddo
      do i=1,nx
       e(i,0)=1.
       f(i,0)=0.
       e(i,1)=2.*c33(i,2)*rec_dn(i)
       f(i,1)=     r(i,1)*rec_dn(i)
      enddo
      do k=2,nz-1
        do i=1,nx
          r(i,k)=r(i,k)+d(i,k)*(beti*p(i,k)-rhs(i,k)) &
                      +dgh(i,k)*p(i,k)
        enddo
        do i=1,nx
          rec_dn(i) = c33(i,k+1)+c33(i,k-1)*(1.-e(i,k-2))+d(i,k)*beti &
                   + dgh(i,k)
        end do
        call vrec(rec_dn, rec_dn, nx)
        do i=1,nx
          e(i,k)= c33(i,k+1)*rec_dn(i)
          f(i,k)=(c33(i,k-1)*f(i,k-2)+r(i,k))*rec_dn(i)
        enddo
      enddo
       do i=1,nx
         rec_dn(i) = d(i,nz)*beti+2.*(1.-e(i,nz-2))*c33(i,nz-1) &
                  + dgh(i,nz)
       end do
       call vrec(rec_dn, rec_dn, nx)
       rec_1me(1:nx) = 1.-e(1:nx,nz-1)
       call vrec(rec_1me, rec_1me, nx)
       do i=1,nx
         r(i,nz)=r(i,nz)+d(i,nz)*(beti*p(i,nz)-rhs(i,nz)) &
                     +dgh(i,nz)*p(i,nz)
       end do
       do i=1,nx
        p(i,nz)=(r(i,nz)+2.*f(i,nz-2)*c33(i,nz-1))*rec_dn(i)
        p(i,nz-1)=f(i,nz-1)*rec_1me(i)
       enddo

      do k=nz-2,1,-1
       do i=1,nx
        p(i,k)=e(i,k)*p(i,k+2)+f(i,k)
       enddo
      enddo


      if(line.eq.1) then
        do k=1,nz
          do i=1,nx
            p(i,k)=oms*po(i,k)+omg*p(i,k)
            po(i,k)=     p(i,k)
          enddo
        enddo
      endif

      if(it.eq.itr) go to 101
      do k=1,nz
        px(1,k)=c11(1,k)*(p(2,k)-p(nx-1,k))
        do i=2,nx-1
          px(i,k)=c11(i,k)*(p(i+1,k)-p(i-1,k))
        enddo
        px(nx,k)=c11(nx,k)*(p(2,k)-p(nx-1,k))

        r(1,k)=(px(2,k)-px(nx-1,k))
        do i=2,nx-1
          r(i,k)=px(i+1,k)-px(i-1,k)
        end do
        r(nx,k)=(px(2,k)-px(nx-1,k))
      enddo

  100 continue
  101 continue
#ifdef TESTMODE
      call tend('  precon')
#endif
      return

      end subroutine precon_1


      subroutine prforc_1(p,pfx,pfz,u,w,nx,nz)    
        implicit none
        integer, intent(in) :: nx, nz
        real, intent(in) :: p (nx,nz), u (nx,nz), w (nx,nz)                                                      
        real, intent(out) :: pfx (nx,nz), pfz (nx,nz)                          
        real:: px (nx, nz), pz (nx, nz) 
        integer :: i, k, nxz

        nxz = nx*nz
      dxil=.5*dxi
      dzil=.5*dzi

!compute pressure derivatives everywhere
      do 18 k=1,nz
      do 1 i=2,nx-1
    1 px(i,k)=     dxil*(p(i+1,k)-p(i-1,k))
      px(1,k)=dxil*(p(2,k)-p(nx-1,k))
      px(nx,k)=dxil*(p(2,k)-p(nx-1,k))
   18 continue
      do 38 i=1,nx
      do 3 k=2,nz-1
    3 pz(i,k)=dzil*(p(i,k+1)-p(i,k-1))
      pz(i,1)= dzi*(p(i,2)-p(i,1))
   38 pz(i,nz)= dzi*(p(i,nz)-p(i,nz-1))

!compute interior pressure forces
      do 21 i=1,nx
      do 10 k=2,nz-1
      pfx(i,k)=u(i,k)-px(i,k)
   10 pfz(i,k)=w(i,k)-pz(i,k)
      pfz(i,1)=0.
      pfz(i,nz)=0.
      pfx(i,1)=u(i,1)-px(i,1)
   21 pfx(i,nz)=u(i,nz)-px(i,nz)

      return
      end subroutine prforc_1


      subroutine laplc_1(p,r,u,w,nx,nz)
        implicit none
        integer, intent(in) :: nx, nz
        real, intent(in) ::  p (nx,nz)
        real, intent(out) :: r (nx,nz), u (nx,nz), w (nx,nz) 
        real :: px (nx, nz), pz (nx, nz) , coef
        integer :: i, k
!BGL  additional variables for optimization
      real rscale
#ifdef TESTMODE
      call tbeg('  lapl')
#endif
      dxil=.5*dxi
      dzil=.5*dzi

!compute pressure derivatives everywhere

      do k=1,nz
         px(1,k) = dxil*(p(2,k)-p(nx-1,k))
         do i=2,nx-1
            px(i,k) = dxil*(p(i+1,k)-p(i-1,k))
         end do
         px(nx,k) = dxil*(p(2,k)-p(nx-1,k))
      end do
      do k=2,nz-1
         do i=1,nx
            pz(i,k) = dzil*(p(i,k+1)-p(i,k-1))
            u(i,k) = -px(i,k)
            w(i,k) = -pz(i,k)
         end do
      enddo
      do i=1,nx
         pz(i,1) = dzi*(p(i,2)-p(i,1))
         w(i,1) = 0.
         u(i,1) = -px(i,1)
      end do
      do i=1,nx
         pz(i,nz) = dzi*(p(i,nz)-p(i,nz-1))
         w(i,nz) = 0.
         u(i,nz) = -px(i,nz)
      end do
      do k=1,nz
         coef=rho0(k)*gac(k)
         do i=1,nx
            u(i,k)=coef*u(i,k)
            w(i,k)=coef*w(i,k)
         end do
      end do

!compute interior pressure forces


!compute laplacian
      do 99 k=1,nz
        r(1,k)=dxil*(u(2,k)-u(nx-1,k))
        do i=2,nx-1
          r(i,k)=dxil*(u(i+1,k)-u(i-1,k))
        end do
        r(nx,k)=dxil*(u(2,k)-u(nx-1,k))

        if (k .eq. 1) then
          do i=1,nx
            r(i,1)=r(i,1)+dzi *(w(i,2)-w(i,1)) 
          end do
        else if (k .eq. nz) then
          do i=1,nx
            r(i,nz)=r(i,nz)+dzi *(w(i,nz)-w(i,nz-1))
          end do
        else
          do i=1,nx
            r(i,k)=r(i,k)+dzil*(w(i,k+1)-w(i,k-1))
          end do
        endif

        rscale = 1.0/(rho0(k)*gac(k))
        do i=1,nx
          r(i,k)=-r(i,k)*rscale
        end do
   99 end do
#ifdef TESTMODE
      call tend('  lapl')
#endif
      return
      end subroutine laplc_1

      SUBROUTINE rhsdiv_1 (u, w, d, r, iflg, nx, nz) 
        implicit none
        integer, intent(in) :: iflg, nx, nz
        real :: u (nx, nz), w (nx, nz), d (nx, nz), r (nx, nz) 
        integer :: i, j, k, nxz


!BGL  variables for optimization
      real rflg, recd(nx)

      nxz = nx*nz
!      n=n1
!      l=l1
!      nl=n*l

      do 200 i=1,nxz
  200 r(i,1)=0.

      dxil=.5*dxi
      dzil=.5*dzi

      do 1 j=1,nz
      do 1 i=2,nx-1
    1 r(i,j)=dxil*(u(i+1,j)*d(i+1,j)-u(i-1,j)*d(i-1,j))
      do 11 j=1,nz
      r(1,j)=dxil*(u(2,j)*d(2,j)-u(nx-1,j)*d(nx-1,j))
      r(nx,j)=dxil*(u(2,j)*d(2,j)-u(nx-1,j)*d(nx-1,j))
   11 continue
      do 3 k=2,nz-1
      do 3 i=1,nx
    3 r(i,k)=r(i,k) &
             +dzil*(w(i,k+1)*d(i,k+1)-w(i,k-1)*d(i,k-1))
      do 13 i=1,nx
      r(i,1)=r(i,1)+dzi*(w(i,2)*d(i,2)-w(i,1)*d(i,1))
   13 r(i,nz)=r(i,nz)+dzi*(w(i,nz)*d(i,nz)-w(i,nz-1)*d(i,nz-1))

      if(iflg.ne.0) then
      rflg = real(iflg)
      do k=1,nz
        call vrec(recd, d(1,k), nx)
        do i=1,nx
          r(i,k)=rflg*r(i,k)*recd(i)
        end do
    4 end do
      endif

      return
      end subroutine rhsdiv_1

  end module gcrk_mod
