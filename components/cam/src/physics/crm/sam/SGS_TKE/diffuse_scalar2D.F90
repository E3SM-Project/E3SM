module diffuse_scalar2D_mod
  implicit none

contains
  subroutine diffuse_scalar2D (ncrms,dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,grdf_z,field,fluxb,fluxt,tkh,rho,rhow,flux)

    use grid
    use params
    use openacc_utils
    implicit none
    integer, intent(in) :: ncrms
    ! input
    integer :: dimx1_d,dimx2_d,dimy1_d,dimy2_d
    real(crm_rknd) grdf_x(ncrms,nzm)! grid factor for eddy diffusion in x
    real(crm_rknd) grdf_z(ncrms,nzm)! grid factor for eddy diffusion in z
    real(crm_rknd) field(ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) ! scalar
    real(crm_rknd) tkh(ncrms,dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm) ! SGS eddy conductivity
    real(crm_rknd) fluxb(ncrms,nx,ny)   ! bottom flux
    real(crm_rknd) fluxt(ncrms,nx,ny)   ! top flux
    real(crm_rknd) rho(ncrms,nzm)
    real(crm_rknd) rhow(ncrms,nz)
    real(crm_rknd) flux(ncrms,nz)
    ! local
    real(crm_rknd), allocatable :: flx(:,:,:,:)
    real(crm_rknd), allocatable :: dfdt(:,:,:,:)
    real(crm_rknd) rdx2,rdz2,rdz,rdx5,rdz5,tmp
    real(crm_rknd) tkx,tkz,rhoi
    integer i,j,k,ib,ic,kc,kb,icrm
    integer :: numgangs  !For working around PGI bug where it didn't create enough OpenACC gangs

    if(.not.dosgs.and..not.docolumn) return

    rdx2=1./(dx*dx)
    j=1

    allocate( flx(ncrms,0:nx,1,0:nzm) )
    allocate( dfdt(ncrms,nx,ny,nzm) )
    call prefetch( flx  )
    call prefetch( dfdt )

    !For working around PGI bug where it didn't create enough OpenACC gangs
    numgangs = ceiling(ncrms*nzm*ny*nx/128.)
    !$acc parallel loop vector_length(128) num_gangs(numgangs) collapse(3) async(asyncid)
    do k = 1 , nzm
      do i = 1 , nx
        do icrm = 1 , ncrms
          dfdt(icrm,i,j,k)=0.
        enddo
      enddo
    enddo

    if(dowallx) then
      if(mod(rank,nsubdomains_x).eq.0) then
        !$acc parallel loop collapse(2) async(asyncid)
        do k=1,nzm
          do icrm = 1 , ncrms
            field(icrm,0,j,k) = field(icrm,1,j,k)
          enddo
        enddo
      endif
      if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
        !$acc parallel loop collapse(2) async(asyncid)
        do k=1,nzm
          do icrm = 1 , ncrms
            field(icrm,nx+1,j,k) = field(icrm,nx,j,k)
          enddo
        enddo
      endif
    endif

    if(.not.docolumn) then
      !$acc parallel loop collapse(3) async(asyncid)
      do k=1,nzm
        do i=0,nx
          do icrm = 1 , ncrms
            rdx5=0.5*rdx2  *grdf_x(icrm,k)
            ic=i+1
            tkx=rdx5*(tkh(icrm,i,j,k)+tkh(icrm,ic,j,k))
            flx(icrm,i,j,k)=-tkx*(field(icrm,ic,j,k)-field(icrm,i,j,k))
          enddo
        enddo
      enddo
      !$acc parallel loop collapse(3) async(asyncid)
      do k=1,nzm
        do i=1,nx
          do icrm = 1 , ncrms
            ib=i-1
            dfdt(icrm,i,j,k)=dfdt(icrm,i,j,k)-(flx(icrm,i,j,k)-flx(icrm,ib,j,k))
          enddo
        enddo
      enddo
    endif

    !$acc parallel loop collapse(2) async(asyncid)
    do k = 1 , nzm
      do icrm = 1 , ncrms
        flux(icrm,k) = 0.
      enddo
    enddo

    !$acc parallel loop collapse(3) async(asyncid)
    do k=1,nzm
      do i=1,nx
        do icrm = 1 , ncrms
          if (k <= nzm-1) then
            kc=k+1
            rhoi = rhow(icrm,kc)/adzw(icrm,kc)
            rdz2=1./(dz(icrm)*dz(icrm))
            rdz5=0.5*rdz2 * grdf_z(icrm,k)
            tkz=rdz5*(tkh(icrm,i,j,k)+tkh(icrm,i,j,kc))
            flx(icrm,i,j,k)=-tkz*(field(icrm,i,j,kc)-field(icrm,i,j,k))*rhoi
            !$acc atomic update
            flux(icrm,kc) = flux(icrm,kc) + flx(icrm,i,j,k)
          elseif (k == nzm) then
            tmp=1./adzw(icrm,nz)
            rdz=1./dz(icrm)
            flx(icrm,i,j,0)=fluxb(icrm,i,j)*rdz*rhow(icrm,1)
            flx(icrm,i,j,nzm)=fluxt(icrm,i,j)*rdz*tmp*rhow(icrm,nz)
            !$acc atomic update
            flux(icrm,1) = flux(icrm,1) + flx(icrm,i,j,0)
          endif
        enddo
      enddo
    enddo

    !$acc parallel loop collapse(3) async(asyncid)
    do k=1,nzm
      do i=1,nx
        do icrm = 1 , ncrms
          kb=k-1
          rhoi = 1./(adz(icrm,k)*rho(icrm,k))
          dfdt(icrm,i,j,k)=dtn*(dfdt(icrm,i,j,k)-(flx(icrm,i,j,k)-flx(icrm,i,j,kb))*rhoi)
          field(icrm,i,j,k)=field(icrm,i,j,k) + dfdt(icrm,i,j,k)
        enddo
      enddo
    enddo

    deallocate( flx )
    deallocate( dfdt )

  end subroutine diffuse_scalar2D

end module diffuse_scalar2D_mod
