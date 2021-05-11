module advect_scalar3D_mod
  use params, only: asyncid
  use openacc_utils
  implicit none

contains

  subroutine advect_scalar3D (ncrms, f, u, v, w, rho, rhow, flux)
    !     positively definite monotonic advection with non-oscillatory option
    use grid
    use params, only: dowallx, dowally, crm_rknd
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) f(ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real(crm_rknd) u(ncrms,dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm)
    real(crm_rknd) v(ncrms,dimx1_v:dimx2_v, dimy1_v:dimy2_v, nzm)
    real(crm_rknd) w(ncrms,dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz )
    real(crm_rknd) rho(ncrms,nzm)
    real(crm_rknd) rhow(ncrms,nz)
    real(crm_rknd) flux(ncrms,nz)
    real(crm_rknd), allocatable :: mx (:,:,:,:)
    real(crm_rknd), allocatable :: mn (:,:,:,:)
    real(crm_rknd), allocatable :: uuu(:,:,:,:)
    real(crm_rknd), allocatable :: vvv(:,:,:,:)
    real(crm_rknd), allocatable :: www(:,:,:,:)
    real(crm_rknd), allocatable :: iadz (:,:)
    real(crm_rknd), allocatable :: irho (:,:)
    real(crm_rknd), allocatable :: irhow(:,:)
    real(crm_rknd) eps, dd
    integer i,j,k,ic,ib,jc,jb,kc,kb, icrm
    logical nonos
    real(crm_rknd) x1, x2, a, b, a1, a2, y
    real(crm_rknd) andiff,across,pp,pn

    !Statement functions
    andiff(x1,x2,a,b)=(abs(a)-a*a*b)*0.5D0*(x2-x1)
    across(x1,a1,a2)=0.03125D0*a1*a2*x1
    pp(y)= max(real(0.,crm_rknd),y)
    pn(y)=-min(real(0.,crm_rknd),y)

    nonos = .true.
    eps = 1.D-10

    allocate( mx (ncrms,0:nxp1 ,0:nyp1 ,nzm) )
    allocate( mn (ncrms,0:nxp1 ,0:nyp1 ,nzm) )
    allocate( uuu(ncrms,-1:nxp3,-1:nyp2,nzm) )
    allocate( vvv(ncrms,-1:nxp2,-1:nyp3,nzm) )
    allocate( www(ncrms,-1:nxp2,-1:nyp2,nz ) )
    allocate( iadz (ncrms,nzm) )
    allocate( irho (ncrms,nzm) )
    allocate( irhow(ncrms,nzm) )

    call prefetch( mx  )
    call prefetch( mn  )
    call prefetch( uuu )
    call prefetch( vvv )
    call prefetch( www )
    call prefetch( iadz  )
    call prefetch( irho  )
    call prefetch( irhow )

    !$acc parallel loop collapse(3) async(asyncid)
    do j = -1 , nyp2
      do i = -1 , nxp2
        do icrm = 1 , ncrms
          www(icrm,i,j,nz)=0.
        enddo
      enddo
    enddo

    if (dowallx) then
      if (mod(rank,nsubdomains_x).eq.0) then
        !$acc parallel loop collapse(4) async(asyncid)
        do k=1,nzm
          do j=dimy1_u,dimy2_u
            do i=dimx1_u,1
              do icrm = 1 , ncrms
                u(icrm,i,j,k) = 0.
              enddo
            enddo
          enddo
        enddo
      endif
      if (mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
        !$acc parallel loop collapse(4) async(asyncid)
        do k=1,nzm
          do j=dimy1_u,dimy2_u
            do i=nx+1,dimx2_u
              do icrm = 1 , ncrms
                u(icrm,i,j,k) = 0.
              enddo
            enddo
          enddo
        enddo
      endif
    endif

    if (dowally) then
      if (rank.lt.nsubdomains_x) then
        !$acc parallel loop collapse(4) async(asyncid)
        do k=1,nzm
          do j=dimy1_v,1
            do i=dimx1_v,dimx2_v
              do icrm = 1 , ncrms
                v(icrm,i,j,k) = 0.
              enddo
            enddo
          enddo
        enddo
      endif
      if (rank.gt.nsubdomains-nsubdomains_x-1) then
        !$acc parallel loop collapse(4) async(asyncid)
        do k=1,nzm
          do j=ny+1,dimy2_v
            do i=dimx1_v,dimx2_v
              do icrm = 1 , ncrms
                v(icrm,i,j,k) = 0.
              enddo
            enddo
          enddo
        enddo
      endif
    endif

    !-----------------------------------------

    if (nonos) then
      !$acc parallel loop collapse(4) async(asyncid)
      do k=1,nzm
        do j=0,nyp1
          do i=0,nxp1
            do icrm = 1 , ncrms
              kc=min(nzm,k+1)
              kb=max(1,k-1)
              jb=j-1
              jc=j+1
              ib=i-1
              ic=i+1
              mx(icrm,i,j,k)=max(f(icrm,ib,j,k),f(icrm,ic,j,k),f(icrm,i,jb,k),f(icrm,i,jc,k),&
                                 f(icrm,i,j,kb),f(icrm,i,j,kc),f(icrm,i,j,k))
              mn(icrm,i,j,k)=min(f(icrm,ib,j,k),f(icrm,ic,j,k),f(icrm,i,jb,k),f(icrm,i,jc,k),&
                                 f(icrm,i,j,kb),f(icrm,i,j,kc),f(icrm,i,j,k))
            enddo
          enddo
        enddo
      enddo
    endif  ! nonos

    !$acc parallel loop collapse(4) async(asyncid)
    do k=1,nzm
      do j=-1,nyp3
        do i=-1,nxp3
          do icrm = 1 , ncrms
            kb=max(1,k-1)
            if (j <= nyp2                ) uuu(icrm,i,j,k)=max(real(0.,crm_rknd),u(icrm,i,j,k))*f(icrm,i-1,j,k)+&
                                                           min(real(0.,crm_rknd),u(icrm,i,j,k))*f(icrm,i,j,k)
            if (i <= nxp2                ) vvv(icrm,i,j,k)=max(real(0.,crm_rknd),v(icrm,i,j,k))*f(icrm,i,j-1,k)+&
                                                           min(real(0.,crm_rknd),v(icrm,i,j,k))*f(icrm,i,j,k)
            if (i <= nxp2 .and. j <= nyp2) www(icrm,i,j,k)=max(real(0.,crm_rknd),w(icrm,i,j,k))*f(icrm,i,j,kb )+&
                                                           min(real(0.,crm_rknd),w(icrm,i,j,k))*f(icrm,i,j,k)
            if (i == -1 .and. j == -1) then
              flux(icrm,k) = 0.
            endif
          enddo
        enddo
      enddo
    enddo

    !$acc parallel loop collapse(2) async(asyncid)
    do k=1,nzm
      do icrm = 1 , ncrms
        irho(icrm,k) = 1./rho(icrm,k)
        iadz(icrm,k) = 1./adz(icrm,k)
        irhow(icrm,k)=1./(rhow(icrm,k)*adz(icrm,k))
      enddo
    enddo

    !$acc parallel loop collapse(4) async(asyncid)
    do k=1,nzm
      do j=-1,nyp2
        do i=-1,nxp2
          do icrm = 1 , ncrms
            if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny) then
              !$acc atomic update
              flux(icrm,k) = flux(icrm,k) + www(icrm,i,j,k)
            endif
            f(icrm,i,j,k)=f(icrm,i,j,k)-( uuu(icrm,i+1,j,k)-uuu(icrm,i,j,k)  & 
                              + vvv(icrm,i,j+1,k)-vvv(icrm,i,j,k)  &
                              +(www(icrm,i,j,k+1)-www(icrm,i,j,k) )*iadz(icrm,k))*irho(icrm,k)
          enddo
        enddo
      enddo
    enddo

    !$acc parallel loop collapse(4) async(asyncid)
    do k=1,nzm
      do j=0,nyp2
        do i=0,nxp2
          do icrm = 1 , ncrms
            if (j <= nyp1) then
              kc=min(nzm,k+1)
              kb=max(1,k-1)
              dd=2./(kc-kb)/adz(icrm,k)
              jb=j-1
              jc=j+1
              ib=i-1
              uuu(icrm,i,j,k)=andiff(f(icrm,ib,j,k),f(icrm,i,j,k),u(icrm,i,j,k),irho(icrm,k)) &
              -(across(f(icrm,ib,jc,k)+f(icrm,i,jc,k)-f(icrm,ib,jb,k)-f(icrm,i,jb,k), &
              u(icrm,i,j,k), v(icrm,ib,j,k)+v(icrm,ib,jc,k)+v(icrm,i,jc,k)+v(icrm,i,j,k)) &
              +across(dd*(f(icrm,ib,j,kc)+f(icrm,i,j,kc)-f(icrm,ib,j,kb)-f(icrm,i,j,kb)), &
              u(icrm,i,j,k), w(icrm,ib,j,k)+w(icrm,ib,j,kc)+w(icrm,i,j,k)+w(icrm,i,j,kc))) *irho(icrm,k)
            endif
            if (i <= nxp1) then
              kc=min(nzm,k+1)
              kb=max(1,k-1)
              dd=2./(kc-kb)/adz(icrm,k)
              jb=j-1
              ib=i-1
              ic=i+1
              vvv(icrm,i,j,k)=andiff(f(icrm,i,jb,k),f(icrm,i,j,k),v(icrm,i,j,k),irho(icrm,k)) &
              -(across(f(icrm,ic,jb,k)+f(icrm,ic,j,k)-f(icrm,ib,jb,k)-f(icrm,ib,j,k), &
              v(icrm,i,j,k), u(icrm,i,jb,k)+u(icrm,i,j,k)+u(icrm,ic,j,k)+u(icrm,ic,jb,k)) &
              +across(dd*(f(icrm,i,jb,kc)+f(icrm,i,j,kc)-f(icrm,i,jb,kb)-f(icrm,i,j,kb)), &
              v(icrm,i,j,k), w(icrm,i,jb,k)+w(icrm,i,j,k)+w(icrm,i,j,kc)+w(icrm,i,jb,kc))) *irho(icrm,k)
            endif
            if (i <= nxp1 .and. j <= nyp1) then
              kb=max(1,k-1)
              jb=j-1
              jc=j+1
              ib=i-1
              ic=i+1
              www(icrm,i,j,k)=andiff(f(icrm,i,j,kb),f(icrm,i,j,k),w(icrm,i,j,k),irhow(icrm,k)) &
              -(across(f(icrm,ic,j,kb)+f(icrm,ic,j,k)-f(icrm,ib,j,kb)-f(icrm,ib,j,k), &
              w(icrm,i,j,k), u(icrm,i,j,kb)+u(icrm,i,j,k)+u(icrm,ic,j,k)+u(icrm,ic,j,kb)) &
              +across(f(icrm,i,jc,k)+f(icrm,i,jc,kb)-f(icrm,i,jb,k)-f(icrm,i,jb,kb), &
              w(icrm,i,j,k), v(icrm,i,j,kb)+v(icrm,i,jc,kb)+v(icrm,i,jc,k)+v(icrm,i,j,k))) *irho(icrm,k)
            endif
          enddo
        enddo
      enddo
    enddo

    !$acc parallel loop collapse(3) async(asyncid)
    do j = -1 , nyp2
      do i = -1 , nxp2
        do icrm = 1 , ncrms
          www(icrm,i,j,1) = 0.
        enddo
      enddo
    enddo

    !---------- non-osscilatory option ---------------
    if (nonos) then
      !$acc parallel loop collapse(4) async(asyncid)
      do k=1,nzm
        do j=0,nyp1
          do i=0,nxp1
            do icrm = 1 , ncrms
              kc=min(nzm,k+1)
              kb=max(1,k-1)
              jb=j-1
              jc=j+1
              ib=i-1
              ic=i+1
              mx(icrm,i,j,k)=max(f(icrm,ib,j,k),f(icrm,ic,j,k),f(icrm,i,jb,k),f(icrm,i,jc,k),&
                                 f(icrm,i,j,kb),f(icrm,i,j,kc),f(icrm,i,j,k),mx(icrm,i,j,k))
              mn(icrm,i,j,k)=min(f(icrm,ib,j,k),f(icrm,ic,j,k),f(icrm,i,jb,k),f(icrm,i,jc,k),&
                                 f(icrm,i,j,kb),f(icrm,i,j,kc),f(icrm,i,j,k),mn(icrm,i,j,k))
            enddo
          enddo
        enddo
      enddo

      !$acc parallel loop collapse(4) async(asyncid)
      do k=1,nzm
        do j=0,nyp1
          do i=0,nxp1
            do icrm = 1 , ncrms
              kc=min(nzm,k+1)
              jc=j+1
              ic=i+1
              mx(icrm,i,j,k)=rho(icrm,k)*(mx(icrm,i,j,k)-f(icrm,i,j,k))/ &
                        ( pn(uuu(icrm,ic,j,k)) + pp(uuu(icrm,i,j,k))+ &
                          pn(vvv(icrm,i,jc,k)) + pp(vvv(icrm,i,j,k))+ &
                         (pn(www(icrm,i,j,kc)) + pp(www(icrm,i,j,k)))*iadz(icrm,k)+eps)
              mn(icrm,i,j,k)=rho(icrm,k)*(f(icrm,i,j,k)-mn(icrm,i,j,k))/ &
                        ( pp(uuu(icrm,ic,j,k)) + pn(uuu(icrm,i,j,k))+ &
                          pp(vvv(icrm,i,jc,k)) + pn(vvv(icrm,i,j,k))+ &
                         (pp(www(icrm,i,j,kc)) + pn(www(icrm,i,j,k)))*iadz(icrm,k)+eps)
            enddo
          enddo
        enddo
      enddo

      !$acc parallel loop collapse(4) async(asyncid)
      do k=1,nzm
        do j=1,nyp1
          do i=1,nxp1
            do icrm = 1 , ncrms
              if (j <= ny) then
                ib=i-1
                uuu(icrm,i,j,k) = pp(uuu(icrm,i,j,k))*min(real(1.,crm_rknd),mx(icrm,i,j,k), mn(icrm,ib,j,k)) &
                                 -pn(uuu(icrm,i,j,k))*min(real(1.,crm_rknd),mx(icrm,ib,j,k),mn(icrm,i,j,k))
              endif
              if (i <= nx) then
                jb=j-1
                vvv(icrm,i,j,k) = pp(vvv(icrm,i,j,k))*min(real(1.,crm_rknd),mx(icrm,i,j,k), mn(icrm,i,jb,k)) &
                                 -pn(vvv(icrm,i,j,k))*min(real(1.,crm_rknd),mx(icrm,i,jb,k),mn(icrm,i,j,k))
              endif
              if (i <= nx .and. j <= ny) then
                kb=max(1,k-1)
                www(icrm,i,j,k) = pp(www(icrm,i,j,k))*min(real(1.,crm_rknd),mx(icrm,i,j,k), mn(icrm,i,j,kb)) &
                                 -pn(www(icrm,i,j,k))*min(real(1.,crm_rknd),mx(icrm,i,j,kb),mn(icrm,i,j,k))
                !$acc atomic update
                flux(icrm,k) = flux(icrm,k) + www(icrm,i,j,k)
              endif
            enddo
          enddo
        enddo
      enddo
    endif ! nonos

    !$acc parallel loop collapse(4) async(asyncid)
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            ! MK: added fix for very small negative values (relative to positive values)
            !     especially  when such large numbers as
            !     hydrometeor concentrations are advected. The reason for negative values is
            !     most likely truncation error.
            kc=k+1
            f(icrm,i,j,k)=max(real(0.,crm_rknd),f(icrm,i,j,k) -(uuu(icrm,i+1,j,k)-uuu(icrm,i,j,k)+&
                             vvv(icrm,i,j+1,k)-vvv(icrm,i,j,k)+(www(icrm,i,j,k+1)-www(icrm,i,j,k))*iadz(icrm,k))*irho(icrm,k))
          enddo
        enddo
      enddo
    enddo

    deallocate( mx  )
    deallocate( mn  )
    deallocate( uuu )
    deallocate( vvv )
    deallocate( www )
    deallocate( iadz  )
    deallocate( irho  )
    deallocate( irhow )

  end subroutine advect_scalar3D

end module advect_scalar3D_mod
