
subroutine shear_prod3D(def2)
	
use vars
implicit none
	
real def2(nx,ny,nzm)
	
real rdx0,rdx,rdx_up,rdx_dn
real rdy0,rdy,rdy_up,rdy_dn
real rdz,rdzw_up,rdzw_dn
integer i,j,k,ib,ic,jb,jc,kb,kc

rdx0=1./dx 
rdy0=1./dy

do k=2,nzm-1  

 kb=k-1
 kc=k+1
 rdz = 1./(dz*adz(k))
 rdzw_up = 1./(dz*adzw(kc))
 rdzw_dn = 1./(dz*adzw(k))
 rdx=rdx0 * sqrt(dx*rdz) ! take into account grid anisotropy
 rdy=rdy0 * sqrt(dy*rdz) 
 rdx_up=rdx0 * sqrt(dx*rdzw_up) 
 rdy_up=rdy0 * sqrt(dy*rdzw_up) 
 rdx_dn=rdx0 * sqrt(dx*rdzw_dn) 
 rdy_dn=rdy0 * sqrt(dy*rdzw_dn) 

 do j=1,ny
   jb=j-YES3D
   jc=j+YES3D
   do i=1,nx
     ib=i-1
     ic=i+1
	 
      def2(i,j,k)=2.* ( &
          ( (u(ic,j,k)-u(i,j,k))*rdx)**2+ &
          ( (v(i,jc,k)-v(i,j,k))*rdy)**2+ &
          ( (w(i,j,kc)-w(i,j,k))*rdz)**2 ) &
        + 0.25 * ( &
          ( (u(ic,jc,k)-u(ic,j ,k))*rdy+(v(ic,jc,k)-v(i ,jc,k))*rdx )**2 +  &
          ( (u(i ,jc,k)-u(i ,j ,k))*rdy+(v(i ,jc,k)-v(ib,jc,k))*rdx )**2 +  &
          ( (u(ic,j ,k)-u(ic,jb,k))*rdy+(v(ic,j ,k)-v(i ,j ,k))*rdx )**2 +  &
          ( (u(i ,j ,k)-u(i ,jb,k))*rdy+(v(i ,j ,k)-v(ib,j ,k))*rdx )**2 )   
      def2(i,j,k)=def2(i,j,k) &
        + 0.25 * ( &
          ( (u(ic,j,kc)-u0(kc)-u(ic,j, k)+u0(k))*rdzw_up+ &
            (w(ic,j,kc)-w(i ,j,kc))*rdx_up )**2 + &
          ( (u(i ,j,kc)-u0(kc)-u(i ,j, k)+u0(k))*rdzw_up+ &
            (w(i ,j,kc)-w(ib,j,kc))*rdx_up )**2 + &
          ( (u(ic,j,k )-u0(k)-u(ic,j,kb)+u0(kb))*rdzw_dn+ &
            (w(ic,j,k )-w(i ,j,k ))*rdx_dn )**2 + &
          ( (u(i ,j,k )-u0(k)-u(i ,j,kb)+u0(kb))*rdzw_dn+ &
            (w(i ,j,k )-w(ib,j,k ))*rdx_dn )**2 )
      def2(i,j,k)=def2(i,j,k) &	
        + 0.25 * ( & 
          ( (v(i,jc,kc)-v0(kc)-v(i,jc, k)+v0(k))*rdzw_up+ &
            (w(i,jc,kc)-w(i,j ,kc))*rdy_up )**2 + &
          ( (v(i,j ,kc)-v0(kc)-v(i,j , k)+v0(k))*rdzw_up+ &
            (w(i,j ,kc)-w(i,jb,kc))*rdy_up )**2 + &
          ( (v(i,jc,k )-v0(k)-v(i,jc,kb)+v0(kb))*rdzw_dn+ &
            (w(i,jc,k )-w(i,j ,k ))*rdy_dn )**2 + &
          ( (v(i,j ,k )-v0(k)-v(i,j ,kb)+v0(kb))*rdzw_dn+ &
            (w(i,j ,k )-w(i,jb,k ))*rdy_dn )**2 )

    end do
 end do
end do ! k


k=1
kc=k+1

rdz = 1./(dz*adz(k))
rdzw_up = 1./(dz*adzw(kc))
rdx=rdx0 * sqrt(dx*rdz) ! take into account grid anisotropy
rdy=rdy0 * sqrt(dy*rdz) 
rdx_up=rdx0 * sqrt(dx*rdzw_up) 
rdy_up=rdy0 * sqrt(dy*rdzw_up) 
	 
do j=1,ny
  jb=j-YES3D
  jc=j+YES3D
  do i=1,nx
     ib=i-1
     ic=i+1
	 	
      def2(i,j,k)=2.* ( &
          ( (u(ic,j,k)-u(i,j,k))*rdx)**2+ &
          ( (v(i,jc,k)-v(i,j,k))*rdy)**2+ &
          ( (w(i,j,kc)-w(i,j,k))*rdz)**2 ) &
        + 0.25 * ( &
          ( (u(ic,jc,k)-u(ic,j ,k))*rdy+(v(ic,jc,k)-v(i ,jc,k))*rdx )**2 +  &
          ( (u(i ,jc,k)-u(i ,j ,k))*rdy+(v(i ,jc,k)-v(ib,jc,k))*rdx )**2 +  &
          ( (u(ic,j ,k)-u(ic,jb,k))*rdy+(v(ic,j ,k)-v(i ,j ,k))*rdx )**2 +  &
          ( (u(i ,j ,k)-u(i ,jb,k))*rdy+(v(i ,j ,k)-v(ib,j ,k))*rdx )**2 )   &
	 + 0.5 * ( &
          ( (v(i,jc,kc)-v0(kc)-v(i,jc, k)+v0(k))*rdzw_up+ &
            (w(i,jc,kc)-w(i,j ,kc))*rdy_up )**2 + &
          ( (v(i,j ,kc)-v0(kc)-v(i,j , k)+v0(k))*rdzw_up+ &
          (w(i,j ,kc)-w(i,jb,kc))*rdy_up )**2 ) &
	 + 0.5 * ( &
          ( (u(ic,j,kc)-u0(kc)-u(ic,j, k)+u0(k))*rdzw_up+ &
          (w(ic,j,kc)-w(i ,j,kc))*rdx_up )**2 + &
          ( (u(i ,j,kc)-u0(kc)-u(i ,j, k)+u0(k))*rdzw_up+ &
          (w(i ,j,kc)-w(ib,j,kc))*rdx_up )**2 )


   end do 
end do
	 
	
k=nzm
kc=k+1
kb=k-1

rdz = 1./(dz*adz(k))
rdzw_dn = 1./(dz*adzw(k))
rdx=rdx0 * sqrt(dx*rdz) ! take into account grid anisotropy
rdy=rdy0 * sqrt(dy*rdz) 
rdx_dn=rdx0 * sqrt(dx*rdzw_dn) 
rdy_dn=rdy0 * sqrt(dy*rdzw_dn) 

do j=1,ny
  jb=j-1*YES3D
  jc=j+1*YES3D
  do i=1,nx
      ib=i-1
      ic=i+1
      def2(i,j,k)=2.* ( &
           ( (u(ic,j,k)-u(i,j,k))*rdx)**2+ &
           ( (v(i,jc,k)-v(i,j,k))*rdy)**2+ &
           ( (w(i,j,kc)-w(i,j,k))*rdz)**2 ) &
       + 0.25 * ( &
           ( (u(ic,jc,k)-u(ic,j ,k))*rdy+(v(ic,jc,k)-v(i ,jc,k))*rdx )**2 +  &
           ( (u(i ,jc,k)-u(i ,j ,k))*rdy+(v(i ,jc,k)-v(ib,jc,k))*rdx )**2 +  &
           ( (u(ic,j ,k)-u(ic,jb,k))*rdy+(v(ic,j ,k)-v(i ,j ,k))*rdx )**2 +  &
           ( (u(i ,j ,k)-u(i ,jb,k))*rdy+(v(i ,j ,k)-v(ib,j ,k))*rdx )**2 )   &
       + 0.5 * ( &
           ( (v(i,jc,k )-v0(k)-v(i,jc,kb)+v0(kb))*rdzw_dn+ &
             (w(i,jc,k )-w(i,j ,k ))*rdy_dn )**2 + &
           ( (v(i,j ,k )-v0(k)-v(i,j ,kb)+v0(kb))*rdzw_dn+ &
             (w(i,j ,k )-w(i,jb,k ))*rdy_dn )**2 ) &
 	+ 0.5 * ( &
           ( (u(ic,j,k )-u0(k)-u(ic,j,kb)+u0(kb))*rdzw_dn+ &
             (w(ic,j,k )-w(i ,j,k ))*rdx_dn )**2 + &
           ( (u(i ,j,k )-u0(k)-u(i ,j,kb)+u0(kb))*rdzw_dn+ &
             (w(i ,j,k )-w(ib,j,k ))*rdx_dn )**2 )
  end do 
end do
	
end

