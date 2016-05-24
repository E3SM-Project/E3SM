module mpdat_mod
!  integer, private, PARAMETER ::n1 = nx + 1, n2 = nz + 1

!  integer, private, PARAMETER :: n1m = n1 - 1, n2m = n2 - 1

!  real, private:: v1 (n1, n2m), v2 (n1m, n2), f1 (n1, n2m), f2 (n1m, n2),   &
!       cp (n1m, n2m), cn (n1m, n2m), mx (n1m, n2m), mn (n1m, n2m)
!  real, public :: sprecip(nx)

  real, allocatable:: v1(:,:), v2(:,:), f1(:,:), f2(:,:),   &
       cp(:,:), cn(:,:), mx(:,:), mn(:,:)
contains    
  subroutine mpdat_2d(u1,u2,x,h,iflg,liner, sprecip, nx, nz)
    implicit none
    integer, intent(in) :: nx, nz, iflg, liner
    real, intent(out) :: sprecip(nx)
    integer ::n1, n2 

    real:: u1(nx+1, nz), u2(nx, nz+1), x(nx, nz), h(nx, nz) 
    
    integer, PARAMETER :: iord0 = 2, isor = 1, nonos = 1, idiv = 0
    real, parameter:: ep = 1.e-12 
    integer :: ip, im, i0, k, jp, j, iord, jm, i
    real :: donor, a, y1, y2
    real ::  x1, x2
    real :: vcorr, b
    real :: vcor31, x0, x3
    real :: vcor32
    real :: vdiv1, a1, a2, a3, r    
    real :: vdiv2, b1, b2, b3, b4
    real :: pn, y
    real :: pp

!BGL  additional variables for optimization
      real rech(nx+1), recy(nx+1), tmpp(nx+1), tmpm(nx+1)
      real recp(nx+1), recn(nx+1)
      real vdyf(nx+1), recx(nx+1)
      integer im_list(nx+1), ip_list(nx+1)

!
#ifdef USE_FSEL
      donor(y1,y2,a) = a*fsel(a,y1,y2)
#else
      donor(y1,y2,a)=cvmgm(y2,y1,a)*a
#endif
!BGL  vdyf(x1,x2,a,r)=(abs(a)-a**2*r)*(abs(x2)-abs(x1))
!BGL 1                               /(abs(x2)+abs(x1)+ep)
      vcorr(a,b,y1,y2)=-0.125*a*b*y1*y2
      vcor31(a,x0,x1,x2,x3,r)= -(a -3.*abs(a)*a/r+2.*a**3/r**2)/3. &
           *(abs(x0)+abs(x3)-abs(x1)-abs(x2))  &
           /(abs(x0)+abs(x3)+abs(x1)+abs(x2)+ep)
      vcor32(a,b,y1,y2,r)=0.25*b/r*(abs(a)-2.*a**2/r)*y1/y2
      vdiv1(a1,a2,a3,r)=0.25*a2*(a3-a1)/r
      vdiv2(a,b1,b2,b3,b4,r)=0.25*a*(b1+b2-b3-b4)/r
      pp(y)= max(0.,y)
      pn(y)=-min(0.,y)
#ifdef TESTMODE
      call tbeg('mpdat')
#endif
      n1=nx+1
      n2=nz+1

      iord=iord0
      if(isor.eq.3) iord=max(iord,3)
      if(liner.eq.1) iord=1


      do j=1,n2-1
        do i=1,n1
          v1(i,j) = u1(i,j)
        end do 
      end do 

      do j=1,n2
        do i=1,n1-1
          v2(i,j) = u2(i,j)
        end do 
      enddo
      
      if(nonos.eq.1) then
        do i = 1, nx
          im_list(i) = i - 1 + (n1 - i)/nx*(n1 - 2)
          ip_list(i) = i + 1 - i/nx*(n1 - 2)
        end do

        do j=1,nz
          jm=max(j-1,1  )
          jp=min0(j+1,nz)
          do i=1,nx
            im = im_list(i)
            ip = ip_list(i)
            mx(i,j)=max(x(im,j),x(i,j),x(ip,j),x(i,jm),x(i,jp))
            mn(i,j)=min(x(im,j),x(i,j),x(ip,j),x(i,jm),x(i,jp))
          end do
        end do
      endif
 
      do 3 k=1,iord
 
      do 331 j=1,n2-1
      do 331 i=2,n1-1
  331 f1(i,j)=donor(x(i-1,j),x(i,j),v1(i,j))
      do j=1,n2-1
      f1(1 ,j)=f1(n1-1,j)
      f1(n1,j)=f1(2,j)
      enddo
      do 332 j=2,n2-1
      do 332 i=1,n1-1
  332 f2(i,j)=donor(x(i,j-1),x(i,j),v2(i,j))
      if (iflg.eq.6 .and. k.eq.1) then
        do i=1,nx
          f2(i, 1)=donor(x(i,  1),x(i,  1),v2(i, 1))
          sprecip(i)=-f2(i,1)
          f2(i,n2)=donor(x(i,nz),x(i,nz),v2(i,n2))
        end do
      else
        do i=1,nx
          f2(i, 1)=-f2(i,  2)
          f2(i,n2)=-f2(i,nz)
        end do
      end if
  
      do 333 j=1,n2-1
        call vrec(rech, h(1,j), n1-1)
        do i=1,n1-1
          x(i,j)=x(i,j)-(f1(i+1,j)-f1(i,j)+f2(i,j+1)-f2(i,j))*rech(i)
        end do
  333 end do
      if(k.eq.iord) go to 6

      do 49 j=1,n2-1
      do 49 i=1,n1
      f1(i,j)=v1(i,j)
   49 v1(i,j)=0.
      do 50 j=1,n2
      do 50 i=1,n1-1
      f2(i,j)=v2(i,j)
   50 v2(i,j)=0.
      do 51 j=2,n2-2
        do i=2,n1-1
          rech(i) = 0.5*(h(i-1,j) + h(i,j))
          tmpp(i) = abs(x(i-1,j+1)) + abs(x(i,j+1))
          tmpm(i) = abs(x(i-1,j-1)) + abs(x(i,j-1))
          recy(i) = (tmpp(i) + tmpm(i) + ep)*rech(i)
        end do

        call vrec(rech(2), rech(2), n1-2)
        call vrec(recy(2), recy(2), n1-2)

        do i=2,n1-1
          x1 = x(i-1,j)
          x2 = x(i,j)
          a  = f1(i,j)
          r  = rech(i)
          vdyf(i) = (abs(a)-a**2*r)*(abs(x2)-abs(x1))
          recx(i) = abs(x2)+abs(x1)+ep
        end do

        call vrec(recx(2), recx(2), n1-2)

        do i=2,n1-1
          v1(i,j)=vdyf(i)*recx(i) &
                +vcorr(f1(i,j),  &
                       f2(i-1,j) + f2(i-1,j+1) + f2(i,j+1) + f2(i,j), &
                       tmpp(i) - tmpm(i), &
                       recy(i))
        end do
   51 end do
      if(idiv.eq.1) then
      do 511 j=2,n2-2
      do 511 i=2,n1-1
  511 v1(i,j)=v1(i,j) &
         -vdiv1(f1(i-1,j),f1(i,j),f1(i+1,j),.5*(h(i-1,j)+h(i,j))) &
         -vdiv2(f1(i,j),f2(i-1,j+1),f2(i,j+1),f2(i-1,j),f2(i,j),  &
                      .5*(h(i-1,j)+h(i,j)))
      endif
      do 52 j=2,n2-1
        do i=2,n1-2
          rech(i) = 0.5*(h(i,j-1) + h(i,j))
          tmpp(i) = abs(x(i+1,j-1)) + abs(x(i+1,j))
          tmpm(i) = abs(x(i-1,j-1)) + abs(x(i-1,j))
          recy(i) = (tmpp(i) + tmpm(i) + ep)*rech(i)
        end do

        call vrec(rech(2), rech(2), n1-3)
        call vrec(recy(2), recy(2), n1-3)

        do i=2,n1-2
          x1 = x(i,j-1)
          x2 = x(i,j)
          a  = f2(i,j)
          r  = rech(i)
          vdyf(i) = (abs(a)-a**2*r)*(abs(x2)-abs(x1))
          recx(i) = abs(x2)+abs(x1)+ep
        end do

        call vrec(recx(2), recx(2), n1-2)

        do i=2,n1-2
          v2(i,j)=vdyf(i)*recx(i) &
                +vcorr(f2(i,j),  &
                       f1(i,j-1) + f1(i,j) + f1(i+1,j) + f1(i+1,j-1), &
                       tmpp(i) - tmpm(i), &
                       recy(i))
        end do
   52 end do

      i0=n1-2
      do j=2,n2-1
      rech(1) = 2.0/(h(1,j-1) + h(1,j))
      tmpp(1) = abs(x(2,j-1))  + abs(x(2,j))
      tmpm(1) = abs(x(i0,j-1)) + abs(x(i0,j))
      recy(1) = rech(1)/(tmpp(1) + tmpm(1) + ep)
      x1 = x(1,j-1)
      x2 = x(1,j)
      a  = f2(1,j)
      r = rech(1)
      vdyf(1) = (abs(a)-a**2*r)*(abs(x2)-abs(x1))
      recx(1) = abs(x2)+abs(x1)+ep
      v2(1,j)=vdyf(1)/recx(1) &
            +vcorr(f2(1,j),  &
                   f1(1,j-1) + f1(1,j) + f1(2,j) + f1(2,j-1), &
                   tmpp(1) - tmpm(1), &
                   recy(1))
      v2(n1-1,j)=v2(1,j)
      enddo

      if(idiv.eq.1) then
      do 521 j=2,n2-1
      do 521 i=1,n1-1
  521 v2(i,j)=v2(i,j) &
         -vdiv1(f2(i,j-1),f2(i,j),f2(i,j+1),.5*(h(i,j-1)+h(i,j))) &
         -vdiv2(f2(i,j),f1(i+1,j),f1(i+1,j-1),f1(i,j-1),f1(i,j),  &
                      .5*(h(i,j-1)+h(i,j)))
      endif
      if(isor.eq.3) then
      do 61 j=2,n2-2
      do 61 i=3,n1-2
   61 v1(i,j)=v1(i,j)     +vcor31(f1(i,j), &
             x(i-2,j),x(i-1,j),x(i,j),x(i+1,j),.5*(h(i-1,j)+h(i,j)))
      do j=2,n2-2
      v1(2,j)=v1(2,j)     +vcor31(f1(2,j), &
             x(n1-2,j),x(1,j),x(2,j),x(3,j),.5*(h(1,j)+h(2,j)))
      v1(n1-1,j)=v1(n1-1,j) +vcor31(f1(n1-1,j),x(n1-3,j),x(n1-2,j), &
                       x(n1-1,j),x(2,j),.5*(h(n1-2,j)+h(n1-1,j)))
      enddo
      do 62 j=2,n2-2
      do 62 i=2,n1-1
   62 v1(i,j)=v1(i,j) &
           +vcor32(f1(i,j),f2(i-1,j)+f2(i-1,j+1)+f2(i,j+1)+f2(i,j), &
           abs(x(i,j+1))-abs(x(i,j-1))-abs(x(i-1,j+1))+abs(x(i-1,j-1)), &
           abs(x(i,j+1))+abs(x(i,j-1))+abs(x(i-1,j+1))+abs(x(i-1,j-1))+ep, &
           .5*(h(i-1,j)+h(i,j)))
      do 63 j=3,n2-2
      do 63 i=1,n1-1
   63 v2(i,j)=v2(i,j)     +vcor31(f2(i,j), &
             x(i,j-2),x(i,j-1),x(i,j),x(i,j+1),.5*(h(i,j-1)+h(i,j)))
      do 64 j=3,n2-2
      do 64 i=2,n1-2
64       v2(i,j)=v2(i,j) &
              +vcor32(f2(i,j),f1(i,j-1)+f1(i+1,j-1)+f1(i+1,j)+f1(i,j), &
              abs(x(i+1,j))-abs(x(i-1,j))-abs(x(i+1,j-1))+abs(x(i-1,j-1)), &
              abs(x(i+1,j))+abs(x(i-1,j))+abs(x(i+1,j-1))+abs(x(i-1,j-1))+ep, &
              .5*(h(i,j-1)+h(i,j)))
      do 641 j=3,n2-2
      v2(1,j)=v2(1,j) &
      +vcor32(f2(1,j),f1(1,j-1)+f1(2,j-1)+f1(2,j)+f1(1,j), &
        abs(x(2,j))-abs(x(n1-2,j))-abs(x(2,j-1))+abs(x(n1-2,j-1)), &
        abs(x(2,j))+abs(x(n1-2,j))+abs(x(2,j-1))+abs(x(n1-2,j-1))+ep, &
                        .5*(h(1,j-1)+h(1,j)))
  641 v2(n1-1,j)=v2(1,j)
      endif
 
        do j=1,nz
          v1( 1,j)=v1(nx,j)
          v1(n1,j)=v1(  2,j)
        end do

        do i=1,nx
          v2(i, 1)=-v2(i,  2)
          v2(i,n2)=-v2(i,nz)
        enddo

                  if(nonos.eq.1) then
!                 non-osscilatory option
      do 402 j=1,nz 
      do 4021 i=2,n1-1
 4021 f1(i,j)=donor(x(i-1,j),x(i,j),v1(i,j))
      f1(1 ,j)=f1(nx,j)
      f1(n1,j)=f1(2  ,j)
  402 continue
     
      do j=2,nz
        do i=1,nx
          f2(i,j)=donor(x(i,j-1),x(i,j),v2(i,j))
        end do
      end do
      if(iflg.ne.6) then
        do i=1,nx
          f2(i, 1)=-f2(i,  2)
          f2(i,n2)=-f2(i,nz)
        end do
      else
        do i=1,nx
          f2(i, 1)=0.
          f2(i,n2)=0.
        end do
      endif

      do 405 j=1,nz   
        do i=1,nx
          recp(i) =  &
         (pn(f1(i+1,j))+pp(f1(i,j))+pn(f2(i,j+1))+pp(f2(i,j))+ep)
          recn(i) = &
         (pp(f1(i+1,j))+pn(f1(i,j))+pp(f2(i,j+1))+pn(f2(i,j))+ep)
        end do

        call vrec(recp, recp, nx)
        call vrec(recn, recn, nx)

        jm=max(j-1,1  )
        jp=min(j+1,nz)
        do i=1,nx
          im = im_list(i)
          ip = ip_list(i)
          mx(i,j)=max(x(im,j),x(i,j),x(ip,j),x(i,jm),x(i,jp),mx(i,j))
          mn(i,j)=min(x(im,j),x(i,j),x(ip,j),x(i,jm),x(i,jp),mn(i,j))
          cp(i,j)=(mx(i,j)-x(i,j))*h(i,j)*recp(i)
          cn(i,j)=(x(i,j)-mn(i,j))*h(i,j)*recn(i)
        end do

        do i=2,nx 
           v1(i,j)=pp(v1(i,j))* &
                ( min(1.,cp(i,j),cn(i-1,j))*pp(sign(1., x(i-1,j))) &
                +min(1.,cp(i-1,j),cn(i,j))*pp(sign(1.,-x(i-1,j))) ) &
                -pn(v1(i,j))* &
                ( min(1.,cp(i-1,j),cn(i,j))*pp(sign(1., x(i ,j ))) &
                +min(1.,cp(i,j),cn(i-1,j))*pp(sign(1.,-x(i ,j ))) )
        end do
        v1( 1,j)=v1(nx,j)
        v1(n1,j)=v1( 2 ,j)
  405 end do

      do 406 j=2,nz 
        do i=1,nx 
          v2(i,j)=pp(v2(i,j))* &
         ( min(1.,cp(i,j),cn(i,j-1))*pp(sign(1., x(i,j-1))) &
          +min(1.,cp(i,j-1),cn(i,j))*pp(sign(1.,-x(i,j-1))) ) &
              -pn(v2(i,j))* &
         ( min(1.,cp(i,j-1),cn(i,j))*pp(sign(1., x(i ,j ))) &
          +min(1.,cp(i,j),cn(i,j-1))*pp(sign(1.,-x(i ,j ))) )
        end do
  406 end do
                  endif
    3                      continue
    6 continue
#ifdef TESTMODE
      call tend('mpdat')
#endif
      return
    end subroutine mpdat_2d
#if !defined(_AIX) && !defined(_BGL)
  function cvmgm(a,ab,ac)
       if(ac.lt.0) then
        cvmgm=a
       else
        cvmgm=ab
       endif
       return
       end
#endif


 end module mpdat_mod
