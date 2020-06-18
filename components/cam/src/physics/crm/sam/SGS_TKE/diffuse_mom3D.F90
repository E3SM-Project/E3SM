module diffuse_mom3D_mod
  use params, only: asyncid
  implicit none

contains

  subroutine diffuse_mom3D(ncrms,grdf_x, grdf_y, grdf_z, dimx1_d, dimx2_d, dimy1_d, dimy2_d, tk)

    !        momentum tendency due to SGS diffusion

    use vars
    use params, only: crm_rknd
    use openacc_utils
    implicit none
    integer, intent(in) :: ncrms
    integer :: dimx1_d, dimx2_d, dimy1_d, dimy2_d
    real(crm_rknd) tk(ncrms,dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm) ! SGS eddy viscosity
    real(crm_rknd) grdf_x(ncrms,nzm)! grid factor for eddy diffusion in x
    real(crm_rknd) grdf_y(ncrms,nzm)! grid factor for eddy diffusion in y
    real(crm_rknd) grdf_z(ncrms,nzm)! grid factor for eddy diffusion in z
    real(crm_rknd) rdx2,rdy2,rdz2,rdz,rdx25,rdy25
    real(crm_rknd) rdx21,rdy21,rdx251,rdy251,rdz25
    real(crm_rknd) dxy,dxz,dyx,dyz,dzx,dzy
    integer i,j,k,ic,ib,jb,jc,kc,kcu,icrm
    real(crm_rknd) tkx, tky, tkz, rhoi, iadzw, iadz
    real(crm_rknd), allocatable :: fu(:,:,:,:)
    real(crm_rknd), allocatable :: fv(:,:,:,:)
    real(crm_rknd), allocatable :: fw(:,:,:,:)

    allocate( fu(ncrms,0:nx,0:ny,nz) )
    allocate( fv(ncrms,0:nx,0:ny,nz) )
    allocate( fw(ncrms,0:nx,0:ny,nz) )
    call prefetch( fu )
    call prefetch( fv )
    call prefetch( fw )

    rdx2=1./(dx*dx)
    rdy2=1./(dy*dy)
    rdx25=0.25*rdx2
    rdy25=0.25*rdy2
    dxy=dx/dy
    dyx=dy/dx

    !$acc parallel loop collapse(4) async(asyncid)
    do k=1,nzm
      do j=1,ny
        do i=0,nx
          do icrm = 1 , ncrms
            jb=j-1
            kc=k+1
            kcu=min(kc,nzm)
            ic=i+1
            dxz=dx/(dz(icrm)*adzw(icrm,kc))
            rdx21=rdx2    * grdf_x(icrm,k)
            rdx251=rdx25  * grdf_x(icrm,k)
            tkx=rdx21*tk(icrm,i,j,k)
            fu(icrm,i,j,k)=-2.*tkx*(u(icrm,ic,j,k)-u(icrm,i,j,k))
            tkx=rdx251*(tk(icrm,i,j,k)+tk(icrm,i,jb,k)+tk(icrm,ic,j,k)+tk(icrm,ic,jb,k))
            fv(icrm,i,j,k)=-tkx*(v(icrm,ic,j,k)-v(icrm,i,j,k)+(u(icrm,ic,j,k)-u(icrm,ic,jb,k))*dxy)
            tkx=rdx251*(tk(icrm,i,j,k)+tk(icrm,ic,j,k)+tk(icrm,i,j,kcu)+tk(icrm,ic,j,kcu))
            fw(icrm,i,j,k)=-tkx*(w(icrm,ic,j,kc)-w(icrm,i,j,kc)+(u(icrm,ic,j,kcu)-u(icrm,ic,j,k))*dxz)
          enddo
        enddo
      enddo
    enddo
    !$acc parallel loop collapse(4) async(asyncid)
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            kc=k+1
            ib=i-1
            dudt(icrm,i,j,k,na)=dudt(icrm,i,j,k,na)-(fu(icrm,i,j,k)-fu(icrm,ib,j,k))
            dvdt(icrm,i,j,k,na)=dvdt(icrm,i,j,k,na)-(fv(icrm,i,j,k)-fv(icrm,ib,j,k))
            dwdt(icrm,i,j,kc,na)=dwdt(icrm,i,j,kc,na)-(fw(icrm,i,j,k)-fw(icrm,ib,j,k))
          enddo
        enddo
      enddo
    enddo

    !$acc parallel loop collapse(4) async(asyncid)
    do k=1,nzm
      do j=0,ny
        do i=1,nx
          do icrm = 1 , ncrms
            jc=j+1
            kc=k+1
            kcu=min(kc,nzm)
            ib=i-1
            dyz=dy/(dz(icrm)*adzw(icrm,kc))
            rdy21=rdy2    * grdf_y(icrm,k)
            rdy251=rdy25  * grdf_y(icrm,k)
            tky=rdy21*tk(icrm,i,j,k)
            fv(icrm,i,j,k)=-2.*tky*(v(icrm,i,jc,k)-v(icrm,i,j,k))
            tky=rdy251*(tk(icrm,i,j,k)+tk(icrm,ib,j,k)+tk(icrm,i,jc,k)+tk(icrm,ib,jc,k))
            fu(icrm,i,j,k)=-tky*(u(icrm,i,jc,k)-u(icrm,i,j,k)+(v(icrm,i,jc,k)-v(icrm,ib,jc,k))*dyx)
            tky=rdy251*(tk(icrm,i,j,k)+tk(icrm,i,jc,k)+tk(icrm,i,j,kcu)+tk(icrm,i,jc,kcu))
            fw(icrm,i,j,k)=-tky*(w(icrm,i,jc,kc)-w(icrm,i,j,kc)+(v(icrm,i,jc,kcu)-v(icrm,i,jc,k))*dyz)
          enddo
        enddo
      enddo
    enddo
    !$acc parallel loop collapse(4) async(asyncid)
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            jb=j-1
            kc=k+1
            dudt(icrm,i,j,k,na)=dudt(icrm,i,j,k,na)-(fu(icrm,i,j,k)-fu(icrm,i,jb,k))
            dvdt(icrm,i,j,k,na)=dvdt(icrm,i,j,k,na)-(fv(icrm,i,j,k)-fv(icrm,i,jb,k))
            dwdt(icrm,i,j,kc,na)=dwdt(icrm,i,j,kc,na)-(fw(icrm,i,j,k)-fw(icrm,i,jb,k))
          enddo
        enddo
      enddo
    enddo

    !$acc parallel loop collapse(2) async(asyncid)
    do k = 1 , nzm
      do icrm = 1 , ncrms
        uwsb(icrm,k)=0.
        vwsb(icrm,k)=0.
      enddo
    enddo

    !-------------------------
    !$acc parallel loop collapse(4) async(asyncid)
    do k=1,nzm-1
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            jb=j-1
            kc=k+1
            ib=i-1
            rdz=1./dz(icrm)
            rdz2 = rdz*rdz * grdf_z(icrm,k)
            rdz25 = 0.25*rdz2
            iadz = 1./adz(icrm,k)
            iadzw= 1./adzw(icrm,kc)
            dzx=dz(icrm)/dx
            dzy=dz(icrm)/dy
            tkz=rdz2*tk(icrm,i,j,k)
            fw(icrm,i,j,kc)=-2.*tkz*(w(icrm,i,j,kc)-w(icrm,i,j,k))*rho(icrm,k)*iadz
            tkz=rdz25*(tk(icrm,i,j,k)+tk(icrm,ib,j,k)+tk(icrm,i,j,kc)+tk(icrm,ib,j,kc))
            fu(icrm,i,j,kc)=-tkz*( (u(icrm,i,j,kc)-u(icrm,i,j,k))*iadzw + (w(icrm,i,j,kc)-w(icrm,ib,j,kc))*dzx)*rhow(icrm,kc)
            tkz=rdz25*(tk(icrm,i,j,k)+tk(icrm,i,jb,k)+tk(icrm,i,j,kc)+tk(icrm,i,jb,kc))
            fv(icrm,i,j,kc)=-tkz*( (v(icrm,i,j,kc)-v(icrm,i,j,k))*iadzw + (w(icrm,i,j,kc)-w(icrm,i,jb,kc))*dzy)*rhow(icrm,kc)
            !$acc atomic update
            uwsb(icrm,kc)=uwsb(icrm,kc)+fu(icrm,i,j,kc)
            !$acc atomic update
            vwsb(icrm,kc)=vwsb(icrm,kc)+fv(icrm,i,j,kc)
          enddo
        enddo
      enddo
    enddo

    !$acc parallel loop collapse(3) async(asyncid)
    do j=1,ny
      do i=1,nx
        do icrm = 1 , ncrms
          rdz=1./dz(icrm)
          rdz2 = rdz*rdz * grdf_z(icrm,nzm-1)
          tkz=rdz2*grdf_z(icrm,nzm)*tk(icrm,i,j,nzm)
          fw(icrm,i,j,nz)=-2.*tkz*(w(icrm,i,j,nz)-w(icrm,i,j,nzm))/adz(icrm,nzm)*rho(icrm,nzm)
          fu(icrm,i,j,1)=fluxbu(icrm,i,j) * rdz * rhow(icrm,1)
          fv(icrm,i,j,1)=fluxbv(icrm,i,j) * rdz * rhow(icrm,1)
          fu(icrm,i,j,nz)=fluxtu(icrm,i,j) * rdz * rhow(icrm,nz)
          fv(icrm,i,j,nz)=fluxtv(icrm,i,j) * rdz * rhow(icrm,nz)
          !$acc atomic update
          uwsb(icrm,1) = uwsb(icrm,1) + fu(icrm,i,j,1)
          !$acc atomic update
          vwsb(icrm,1) = vwsb(icrm,1) + fv(icrm,i,j,1)
        enddo
      enddo
    enddo

    !$acc parallel loop collapse(4) async(asyncid)
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            kc=k+1
            rhoi = 1./(rho(icrm,k)*adz(icrm,k))
            dudt(icrm,i,j,k,na)=dudt(icrm,i,j,k,na)-(fu(icrm,i,j,kc)-fu(icrm,i,j,k))*rhoi
            dvdt(icrm,i,j,k,na)=dvdt(icrm,i,j,k,na)-(fv(icrm,i,j,kc)-fv(icrm,i,j,k))*rhoi
          enddo
        enddo
      enddo ! k
    enddo

    !$acc parallel loop collapse(4) async(asyncid)
    do k=2,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            rhoi = 1./(rhow(icrm,k)*adzw(icrm,k))
            dwdt(icrm,i,j,k,na)=dwdt(icrm,i,j,k,na)-(fw(icrm,i,j,k+1)-fw(icrm,i,j,k))*rhoi
          enddo
        enddo
      enddo ! k
    enddo

    deallocate( fu )
    deallocate( fv )
    deallocate( fw )

  end subroutine diffuse_mom3D

end module diffuse_mom3D_mod
