module diffuse_mom2D_mod
  use params, only: asyncid
  implicit none

contains

  subroutine diffuse_mom2D(ncrms,grdf_x, grdf_z, dimx1_d, dimx2_d, dimy1_d, dimy2_d, tk)

    !        momentum tendency due to SGS diffusion

    use vars
    use params, only: crm_rknd
    use openacc_utils
    implicit none
    integer, intent(in) :: ncrms
    integer :: dimx1_d, dimx2_d, dimy1_d, dimy2_d
    real(crm_rknd) tk(ncrms,dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm) ! SGS eddy viscosity
    real(crm_rknd) grdf_x(ncrms,nzm)! grid factor for eddy diffusion in x
    real(crm_rknd) grdf_z(ncrms,nzm)! grid factor for eddy diffusion in z

    real(crm_rknd) rdx2,rdz2,rdz,rdx25,rdz25,rdx21,rdx251
    real(crm_rknd) dxz,dzx

    integer i,j,k,ic,ib,kc,kcu,icrm
    real(crm_rknd) tkx, tkz, rhoi, iadzw, iadz
    real(crm_rknd), allocatable :: fu(:,:,:,:)
    real(crm_rknd), allocatable :: fv(:,:,:,:)
    real(crm_rknd), allocatable :: fw(:,:,:,:)
    integer :: numgangs  !For working around PGI bugs where PGI did not allocate enough gangs

    allocate( fu(ncrms,0:nx,1,nz) )
    allocate( fv(ncrms,0:nx,1,nz) )
    allocate( fw(ncrms,0:nx,1,nz) )
    call prefetch( fu )
    call prefetch( fv )
    call prefetch( fw )

    rdx2=1./dx/dx
    rdx25=0.25*rdx2

    j=1

    !For working around PGI bugs where PGI did not allocate enough gangs
    numgangs = ceiling( ncrms*nzm*nx/128. )
    !$acc parallel loop gang vector collapse(3) vector_length(128) num_gangs(numgangs) async(asyncid)
    do k=1,nzm
      do i=0,nx
        do icrm = 1 , ncrms
          kc=k+1
          kcu=min(kc,nzm)
          dxz=dx/(dz(icrm)*adzw(icrm,kc))
          rdx21=rdx2 * grdf_x(icrm,k)
          rdx251=rdx25 * grdf_x(icrm,k)
          ic=i+1
          tkx=rdx21*tk(icrm,i,j,k)
          fu(icrm,i,j,k)=-2.*tkx*(u(icrm,ic,j,k)-u(icrm,i,j,k))
          fv(icrm,i,j,k)=-tkx*(v(icrm,ic,j,k)-v(icrm,i,j,k))
          tkx=rdx251*(tk(icrm,i,j,k)+tk(icrm,ic,j,k)+tk(icrm,i,j,kcu)+tk(icrm,ic,j,kcu))
          fw(icrm,i,j,k)=-tkx*(w(icrm,ic,j,kc)-w(icrm,i,j,kc)+(u(icrm,ic,j,kcu)-u(icrm,ic,j,k))*dxz)
        end do
      end do
    end do
    !For working around PGI bugs where PGI did not allocate enough gangs
    numgangs = ceiling( ncrms*nzm*nx/128. )
    !$acc parallel loop gang vector collapse(3) vector_length(128) num_gangs(numgangs) async(asyncid)
    do k=1,nzm
      do i=1,nx
        do icrm = 1 , ncrms
          kc=k+1
          ib=i-1
          dudt(icrm,i,j,k,na)=dudt(icrm,i,j,k,na)-(fu(icrm,i,j,k)-fu(icrm,ib,j,k))
          dvdt(icrm,i,j,k,na)=dvdt(icrm,i,j,k,na)-(fv(icrm,i,j,k)-fv(icrm,ib,j,k))
          dwdt(icrm,i,j,kc,na)=dwdt(icrm,i,j,kc,na)-(fw(icrm,i,j,k)-fw(icrm,ib,j,k))
        end do
      end do
    end do

    !-------------------------

    !$acc parallel loop collapse(2) async(asyncid)
    do k = 1 , nzm
      do icrm = 1 , ncrms
        uwsb(icrm,k)=0.
        vwsb(icrm,k)=0.
      enddo
    enddo

    !For working around PGI bugs where PGI did not allocate enough gangs
    numgangs = ceiling( ncrms*(nzm-1)*nx/128. )
    !$acc parallel loop gang vector collapse(3) vector_length(128) num_gangs(numgangs) async(asyncid)
    do k=1,nzm-1
      do i=1,nx
        do icrm = 1 , ncrms
          kc=k+1
          rdz=1./dz(icrm)
          rdz2=rdz*rdz *grdf_z(icrm,k)
          rdz25=0.25*rdz2
          dzx=dz(icrm)/dx
          iadz = 1./adz(icrm,k)
          iadzw= 1./adzw(icrm,kc)
          ib=i-1
          tkz=rdz2*tk(icrm,i,j,k)
          fw(icrm,i,j,kc)=-2.*tkz*(w(icrm,i,j,kc)-w(icrm,i,j,k))*rho(icrm,k)*iadz
          tkz=rdz25*(tk(icrm,i,j,k)+tk(icrm,ib,j,k)+tk(icrm,i,j,kc)+tk(icrm,ib,j,kc))
          fu(icrm,i,j,kc)=-tkz*( (u(icrm,i,j,kc)-u(icrm,i,j,k))*iadzw + (w(icrm,i,j,kc)-w(icrm,ib,j,kc))*dzx)*rhow(icrm,kc)
          fv(icrm,i,j,kc)=-tkz*(v(icrm,i,j,kc)-v(icrm,i,j,k))*iadzw*rhow(icrm,kc)
          !$acc atomic update
          uwsb(icrm,kc)=uwsb(icrm,kc)+fu(icrm,i,j,kc)
          !$acc atomic update
          vwsb(icrm,kc)=vwsb(icrm,kc)+fv(icrm,i,j,kc)
        end do
      end do
    end do

    !$acc parallel loop collapse(2) async(asyncid)
    do i=1,nx
      do icrm = 1 , ncrms
        rdz=1./dz(icrm)
        rdz2=rdz*rdz *grdf_z(icrm,nzm-1)
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
      end do
    end do

    !$acc parallel loop collapse(3) async(asyncid)
    do k=1,nzm
      do i=1,nx
        do icrm = 1 , ncrms
          kc=k+1
          rhoi = 1./(rho(icrm,k)*adz(icrm,k))
          dudt(icrm,i,j,k,na)=dudt(icrm,i,j,k,na)-(fu(icrm,i,j,kc)-fu(icrm,i,j,k))*rhoi
          dvdt(icrm,i,j,k,na)=dvdt(icrm,i,j,k,na)-(fv(icrm,i,j,kc)-fv(icrm,i,j,k))*rhoi
        end do
      end do ! k
    end do ! k

    !$acc parallel loop collapse(3) async(asyncid)
    do k=2,nzm
      do i=1,nx
        do icrm = 1 , ncrms
          rhoi = 1./(rhow(icrm,k)*adzw(icrm,k))
          dwdt(icrm,i,j,k,na)=dwdt(icrm,i,j,k,na)-(fw(icrm,i,j,k+1)-fw(icrm,i,j,k))*rhoi
        end do
      end do ! k
    end do ! k

    deallocate( fu )
    deallocate( fv )
    deallocate( fw )

  end subroutine diffuse_mom2D


end module diffuse_mom2D_mod
