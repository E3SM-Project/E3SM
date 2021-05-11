module diffuse_scalar3D_mod
  implicit none

contains

  subroutine diffuse_scalar3D (ncrms,dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,grdf_y,grdf_z,field,fluxb,fluxt,tkh,rho,rhow,flux)

    use grid
    use params
    use task_util_mod, only: task_rank_to_index
    use openacc_utils
    implicit none
    integer, intent(in) :: ncrms
    ! input
    integer :: dimx1_d,dimx2_d,dimy1_d,dimy2_d
    real(crm_rknd) grdf_x(ncrms,nzm)! grid factor for eddy diffusion in x
    real(crm_rknd) grdf_y(ncrms,nzm)! grid factor for eddy diffusion in y
    real(crm_rknd) grdf_z(ncrms,nzm)! grid factor for eddy diffusion in z
    real(crm_rknd) field(ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) ! scalar
    real(crm_rknd) tkh(ncrms,0:nxp1,1-YES3D:nyp1,nzm) ! eddy conductivity
    real(crm_rknd) fluxb(ncrms,nx,ny)   ! bottom flux
    real(crm_rknd) fluxt(ncrms,nx,ny)   ! top flux
    real(crm_rknd) rho(ncrms,nzm)
    real(crm_rknd) rhow(ncrms,nz)
    real(crm_rknd) flux(ncrms,nz)
    ! local
    real(crm_rknd), allocatable :: flx_x(:,:,:,:)
    real(crm_rknd), allocatable :: flx_y(:,:,:,:)
    real(crm_rknd), allocatable :: flx_z(:,:,:,:)
    real(crm_rknd), allocatable :: dfdt (:,:,:,:)
    real(crm_rknd) rdx2,rdy2,rdz2,rdz,rdx5,rdy5,rdz5,tmp
    real(crm_rknd) dxy,dyx,tkx,tky,tkz,rhoi
    integer i,j,k,ib,ic,jb,jc,kc,kb,icrm

    if(.not.dosgs) return

    allocate( flx_x(ncrms,0:nx,0:ny,0:nzm) )
    allocate( flx_y(ncrms,0:nx,0:ny,0:nzm) )
    allocate( flx_z(ncrms,0:nx,0:ny,0:nzm) )
    allocate( dfdt (ncrms,nx,ny,nz) )
    call prefetch( flx_x )
    call prefetch( flx_y )
    call prefetch( flx_z )
    call prefetch( dfdt  )

    rdx2=1.D0/(dx*dx)
    rdy2=1.D0/(dy*dy)
    dxy=dx/dy
    dyx=dy/dx

    !$acc parallel loop collapse(4) async(asyncid)
    do k = 1 , nzm
      do j = 1 , ny
        do i = 1 , nx
          do icrm = 1 , ncrms
            dfdt(icrm,i,j,k)=0.
          enddo
        enddo
      enddo
    enddo

    !-----------------------------------------
    if(dowallx) then
      if(mod(rank,nsubdomains_x).eq.0) then
        !$acc parallel loop collapse(3) async(asyncid)
        do k=1,nzm
          do j=1,ny
            do icrm = 1 , ncrms
              field(icrm,0,j,k) = field(icrm,1,j,k)
            enddo
          enddo
        enddo
      endif
      if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
        !$acc parallel loop collapse(3) async(asyncid)
        do k=1,nzm
          do j=1,ny
            do icrm = 1 , ncrms
              field(icrm,nx+1,j,k) = field(icrm,nx,j,k)
            enddo
          enddo
        enddo
      endif
    endif

    if(dowally) then
      if(rank.lt.nsubdomains_x) then
        !$acc parallel loop collapse(3) async(asyncid)
        do k=1,nzm
          do i=1,nx
            do icrm = 1 , ncrms
              field(icrm,i,1-YES3D,k) = field(icrm,i,1,k)
            enddo
          enddo
        enddo
      endif
      if(rank.gt.nsubdomains-nsubdomains_x-1) then
        !$acc parallel loop collapse(3) async(asyncid)
        do k=1,nzm
          do i=1,nx
            do icrm = 1 , ncrms
              field(icrm,i,ny+YES3D,k) = field(icrm,i,ny,k)
            enddo
          enddo
        enddo
      endif
    endif

    if(dowally) then
      !$acc parallel loop collapse(3) async(asyncid)
      do k=1,nzm
        do i=1,nx
          do icrm = 1 , ncrms
            field(icrm,i,1-YES3D,k) = field(icrm,i,1,k)
          enddo
        enddo
      enddo
      !$acc parallel loop collapse(3) async(asyncid)
      do k=1,nzm
        do i=1,nx
          do icrm = 1 , ncrms
            field(icrm,i,ny+YES3D,k) = field(icrm,i,ny,k)
          enddo
        enddo
      enddo
    endif

    !  Horizontal diffusion:
    !$acc parallel loop collapse(4) async(asyncid)
    do k=1,nzm
      do j=0,ny
        do i=0,nx
          do icrm = 1 , ncrms
            if (j >= 1) then
              ic=i+1
              rdx5=0.5D0*rdx2  * grdf_x(icrm,k)
              tkx=rdx5*(tkh(icrm,i,j,k)+tkh(icrm,ic,j,k))
              flx_x(icrm,i,j,k)=-tkx*(field(icrm,ic,j,k)-field(icrm,i,j,k))
            endif
            if (i >= 1) then
              jc=j+1
              rdy5=0.5D0*rdy2  * grdf_y(icrm,k)
              tky=rdy5*(tkh(icrm,i,j,k)+tkh(icrm,i,jc,k))
              flx_y(icrm,i,j,k)=-tky*(field(icrm,i,jc,k)-field(icrm,i,j,k))
            endif
          enddo
        enddo
      enddo
    enddo
    !$acc parallel loop collapse(4) async(asyncid)
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            ib=i-1
            dfdt(icrm,i,j,k)=dfdt(icrm,i,j,k)-(flx_x(icrm,i,j,k)-flx_x(icrm,ib,j,k))
            jb=j-1
            dfdt(icrm,i,j,k)=dfdt(icrm,i,j,k)-(flx_y(icrm,i,j,k)-flx_y(icrm,i,jb,k))
          enddo
        enddo
      enddo ! k
    enddo

    !  Vertical diffusion:
    !$acc parallel loop collapse(2) async(asyncid)
    do k = 1 , nzm
      do icrm = 1 , ncrms
        flux(icrm,k) = 0.
      enddo
    enddo

    !$acc parallel loop collapse(4) async(asyncid)
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            if (k <= nzm-1) then
              kc=k+1
              rhoi = rhow(icrm,kc)/adzw(icrm,kc)
              rdz2=1.D0/(dz(icrm)*dz(icrm))
              rdz5=0.5D0*rdz2 * grdf_z(icrm,k)
              tkz=rdz5*(tkh(icrm,i,j,k)+tkh(icrm,i,j,kc))
              flx_z(icrm,i,j,k)=-tkz*(field(icrm,i,j,kc)-field(icrm,i,j,k))*rhoi
              !$acc atomic update
              flux(icrm,kc) = flux(icrm,kc) + flx_z(icrm,i,j,k)
            elseif (k == nzm) then
              tmp=1.D0/adzw(icrm,nz)
              rdz=1.D0/dz(icrm)
              flx_z(icrm,i,j,0)=fluxb(icrm,i,j)*rdz*rhow(icrm,1)
              flx_z(icrm,i,j,nzm)=fluxt(icrm,i,j)*rdz*tmp*rhow(icrm,nz)
              !$acc atomic update
              flux(icrm,1) = flux(icrm,1) + flx_z(icrm,i,j,0)
            endif
          enddo
        enddo
      enddo
    enddo

    !$acc parallel loop collapse(4) async(asyncid)
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            kb=k-1
            rhoi = 1.D0/(adz(icrm,k)*rho(icrm,k))
            dfdt(icrm,i,j,k)=dtn*(dfdt(icrm,i,j,k)-(flx_z(icrm,i,j,k)-flx_z(icrm,i,j,kb))*rhoi)
            field(icrm,i,j,k)=field(icrm,i,j,k)+dfdt(icrm,i,j,k)
          enddo
        enddo
      enddo
    enddo

    deallocate( flx_x )
    deallocate( flx_y )
    deallocate( flx_z )
    deallocate( dfdt  )

  end subroutine diffuse_scalar3D

end module diffuse_scalar3D_mod
