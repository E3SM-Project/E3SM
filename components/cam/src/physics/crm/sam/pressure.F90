module pressure_mod
  use task_util_mod
  implicit none

contains

  ! Non-blocking receives before blocking sends
  subroutine pressure(ncrms)
    !       Original pressure solver based on horizontal slabs
    !       (C) 1998, 2002 Marat Khairoutdinov
    !       Works only when the number of slabs is equal to the number of processors.
    !       Therefore, the number of processors shouldn't exceed the number of levels nzm
    !       Also, used for a 2D version
    !       For more processors for the given number of levels and 3D, use pressure_big
    use vars
    use params, only: dowallx, dowally, docolumn, crm_rknd
    use press_rhs_mod
    use press_grad_mod
    use fft_mod
    use openacc_utils
    implicit none
    integer, intent(in) :: ncrms
    integer, parameter :: npressureslabs = nsubdomains
    integer, parameter :: nzslab = max(1,nzm / npressureslabs)
    integer, parameter :: nx2=nx_gl+2, ny2=ny_gl+2*YES3D
    integer, parameter :: n3i=3*nx_gl/2+1,n3j=3*ny_gl/2+1
    real(crm_rknd) work(nx2,ny2)
    real(crm_rknd) ftmp(nx2,ny2)
    real(crm_rknd) ftmp_x(nx2)
    real(crm_rknd) ftmp_y(ny2)
    real(8) b,e
    real(8) xi,xj,xnx,xny,ddx2,ddy2,pii,factx,facty
    real(8) alfa(nzm-1),beta(nzm-1)
    integer i, j, k, id, jd, m, n, it, jt, ii, jj, icrm
    integer nyp22
    real(crm_rknd), allocatable :: f (:,:,:,:)       ! global rhs and array for FTP coefficeients
    real(crm_rknd), allocatable :: ff(:,:,:,:)  ! local (subdomain's) version of f
    integer       , allocatable :: iii(:)
    integer       , allocatable :: jjj(:)
    integer       , allocatable :: ifaxi(:)
    integer       , allocatable :: ifaxj(:)
    real(crm_rknd), allocatable :: trigxi(:)
    real(crm_rknd), allocatable :: trigxj(:)
    real(8)       , allocatable :: a(:,:)
    real(8)       , allocatable :: c(:,:)
    integer iwall,jwall
    integer :: numgangs  !For working aroung PGI OpenACC bug where it didn't create enough gangs
    real(8), allocatable :: eign(:,:)

    allocate( f (ncrms,nx2,ny2,nzslab)      )
    allocate( ff(ncrms,nx+1,ny+2*YES3D,nzm) )
    allocate( iii(0:nx_gl) )
    allocate( jjj(0:ny_gl) )
    allocate( ifaxi(100) )
    allocate( ifaxj(100) )
    allocate( trigxi(n3i) )
    allocate( trigxj(n3j) )
    allocate( a(ncrms,nzm) )
    allocate( c(ncrms,nzm) )

    call prefetch( f      )
    call prefetch( ff     )
    call prefetch( iii    )
    call prefetch( jjj    )
    call prefetch( ifaxi  )
    call prefetch( ifaxj  )
    call prefetch( trigxi )
    call prefetch( trigxj )
    call prefetch( a      )
    call prefetch( c      )

    it = 0
    jt = 0

    !-----------------------------------------------------------------

    if(docolumn) return

    if(dowallx) then
      iwall=1
    else
      iwall=0
    endif
    if(RUN2D) then
      nyp22=1
      jwall=0
    else
      nyp22=nyp2
      if(dowally) then
        jwall=2
      else
        jwall=0
      endif
    endif

    allocate(eign(nxp1-iwall,nyp22-jwall))
    call prefetch(eign)
    !-----------------------------------------------------------------
    !  Compute the r.h.s. of the Poisson equation for pressure
    call press_rhs(ncrms)

    !-----------------------------------------------------------------
    !   Form the horizontal slabs of right-hand-sides of Poisson equation
    n = 0
    !$acc parallel loop collapse(4) async(asyncid)
    do k = 1,nzslab
      do j = 1,ny
        do i = 1,nx
          do icrm = 1 , ncrms
            f(icrm,i,j,k) = p(icrm,i,j,k)
          enddo
        enddo
      enddo
    enddo

    !-------------------------------------------------
    ! Perform Fourier transformation for a slab:
    !$acc parallel loop async(asyncid)
    do icrm = 1 , 1
      call fftfax_crm(nx_gl,ifaxi,trigxi)
      if(RUN3D) call fftfax_crm(ny_gl,ifaxj,trigxj)
    enddo
    !$acc parallel loop gang vector collapse(3) private(work,ftmp_x) async(asyncid)
    do k=1,nzslab
      do j = 1 , ny_gl
        do icrm = 1 , ncrms
          !$acc cache(ftmp_x,work)
          ftmp_x = f(icrm,:,j,k)
          call fft991_crm(ftmp_x,work,trigxi,ifaxi,1,nx2,nx_gl,1,-1)
          f(icrm,:,j,k) = ftmp_x
        enddo
      enddo
    enddo
    if(RUN3D) then
      !$acc parallel loop gang vector collapse(3) private(work,ftmp_y) async(asyncid)
      do k=1,nzslab
        do i = 1 , nx_gl+1
          do icrm = 1 , ncrms
            !$acc cache(ftmp_y,work)
            ftmp_y = f(icrm,i,:,k)
            call fft991_crm(ftmp_y,work,trigxj,ifaxj,1,nx2,ny_gl,1,-1)
            f(icrm,i,:,k) = ftmp_y
          enddo
        enddo
      enddo
    endif

    !-------------------------------------------------
    !   Send Fourier coeffiecients back to subdomains:
    !$acc parallel loop collapse(4) async(asyncid)
    do k = 1,nzslab
      do j = 1,nyp22-jwall
        do i = 1,nxp1-iwall
          do icrm = 1 , ncrms
            ff(icrm,i,j,k) = f(icrm,i,j,k)
          enddo
        enddo
      enddo
    enddo

    !-------------------------------------------------
    !   Solve the tri-diagonal system for Fourier coeffiecients
    !   in the vertical for each subdomain:
    !$acc parallel loop collapse(2) async(asyncid)
    do k=1,nzm
      do icrm = 1 , ncrms
        a(icrm,k)=rhow(icrm,k)/(adz(icrm,k)*adzw(icrm,k)*dz(icrm)*dz(icrm))
        c(icrm,k)=rhow(icrm,k+1)/(adz(icrm,k)*adzw(icrm,k+1)*dz(icrm)*dz(icrm))
      enddo
    enddo

    !$acc parallel loop collapse(2) async(asyncid)
    do j=1,nyp22-jwall
      do i=1,nxp1-iwall
        ddx2=1._8/(dx*dx)
        ddy2=1._8/(dy*dy)
        pii = 3.14159265358979323846D0
        xnx=pii/nx_gl
        xny=pii/ny_gl
        if(dowally) then
          jd=j+jt-1
          facty = 1.d0
        else
          jd=(j+jt-0.1)/2.
          facty = 2.d0
        endif
        xj=jd
        if(dowallx) then
          id=i+it-1
          factx = 1.d0
        else
          id=(i+it-0.1)/2.
          factx = 2.d0
        endif
        xi=id
        eign(i,j)=(2._8*cos(factx*xnx*xi)-2._8)*ddx2+(2._8*cos(facty*xny*xj)-2._8)*ddy2
      enddo
    enddo

    !For working aroung PGI OpenACC bug where it didn't create enough gangs
    numgangs = ceiling(ncrms*(nyp22-jwall)*(nxp2-iwall)/128.)
    !$acc parallel loop gang vector collapse(3) vector_length(128) num_gangs(numgangs) private(alfa,beta) async(asyncid)
    do j=1,nyp22-jwall
      do i=1,nxp1-iwall
        do icrm = 1 , ncrms
          !$acc cache(alfa,beta)
          if(dowally) then
            jd=j+jt-1
          else
            jd=(j+jt-0.1)/2.
          endif
          if(dowallx) then
            id=i+it-1
          else
            id=(i+it-0.1)/2.
          endif
          if(id+jd.eq.0) then
            b=1._8/(eign(i,j)*rho(icrm,1)-a(icrm,1)-c(icrm,1))
            alfa(1)=-c(icrm,1)*b
            beta(1)=ff(icrm,i,j,1)*b
          else
            b=1._8/(eign(i,j)*rho(icrm,1)-c(icrm,1))
            alfa(1)=-c(icrm,1)*b
            beta(1)=ff(icrm,i,j,1)*b
          endif
          do k=2,nzm-1
            e=1._8/(eign(i,j)*rho(icrm,k)-a(icrm,k)-c(icrm,k)+a(icrm,k)*alfa(k-1))
            alfa(k)=-c(icrm,k)*e
            beta(k)=(ff(icrm,i,j,k)-a(icrm,k)*beta(k-1))*e
          enddo
          ff(icrm,i,j,nzm)=(ff(icrm,i,j,nzm)-a(icrm,nzm)*beta(nzm-1))/(eign(i,j)*rho(icrm,nzm)-a(icrm,nzm)+a(icrm,nzm)*alfa(nzm-1))
          do k=nzm-1,1,-1
            ff(icrm,i,j,k)=alfa(k)*ff(icrm,i,j,k+1)+beta(k)
          enddo
        enddo
      enddo
    enddo

    !-----------------------------------------------------------------
    n = 0
    !$acc parallel loop collapse(4) async(asyncid)
    do k = 1,nzslab
      do j = 1,nyp22-jwall
        do i = 1,nxp1-iwall
          do icrm = 1 , ncrms
            f(icrm,i,j,k) = ff(icrm,i,j,k)
          enddo
        enddo
      enddo
    enddo

    !-------------------------------------------------
    !   Perform inverse Fourier transformation:
    if(RUN3D) then
      !$acc parallel loop gang vector collapse(3) private(ftmp_y,work) async(asyncid)
      do k=1,nzslab
        do i = 1 , nx_gl+1
          do icrm = 1 , ncrms
            !$acc cache(ftmp_y,work)
            ftmp_y = f(icrm,i,:,k)
            call fft991_crm(ftmp_y,work,trigxj,ifaxj,1,nx2,ny_gl,1,+1)
            f(icrm,i,:,k) = ftmp_y
          enddo
        enddo
      enddo
    endif
    !$acc parallel loop gang vector collapse(3) private(ftmp_x,work) async(asyncid)
    do k=1,nzslab
      do j = 1 , ny_gl
        do icrm = 1 , ncrms
          !$acc cache(ftmp_x,work)
          ftmp_x = f(icrm,:,j,k)
          call fft991_crm(ftmp_x,work,trigxi,ifaxi,1,nx2,nx_gl,1,+1)
          f(icrm,:,j,k) = ftmp_x
        enddo
      enddo
    enddo

    !-----------------------------------------------------------------
    !   Fill the pressure field for each subdomain:
    !$acc parallel loop async(asyncid)
    do icrm = 1,1
      do i=1,nx_gl
        iii(i)=i
      enddo
      iii(0)=nx_gl
      do j=1,ny_gl
        jjj(j)=j
      enddo
      jjj(0)=ny_gl
    enddo

    n = 0
    !$acc parallel loop collapse(4) async(asyncid)
    do k = 1,nzslab
      do j = 1-YES3D,ny
        do i = 0,nx
          do icrm = 1 , ncrms
            jj=jjj(j)
            ii=iii(i)
            p(icrm,i,j,k) = f(icrm,ii,jj,k)
          enddo
        enddo
      enddo
    enddo

    !  Add pressure gradient term to the rhs of the momentum equation:
    call press_grad(ncrms)

    deallocate(eign)

    deallocate( f  )
    deallocate( ff )
    deallocate( iii )
    deallocate( jjj )
    deallocate( ifaxi )
    deallocate( ifaxj )
    deallocate( trigxi )
    deallocate( trigxj )
    deallocate( a )
    deallocate( c )

  end subroutine pressure

end module pressure_mod
