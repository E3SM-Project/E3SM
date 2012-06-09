#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module checksum_mod
  use edge_mod, only : ghostbuffer_t, edgebuffer_t, ghostbuffer3d_t
  implicit none

  private

  ! test_ghost is called from within a openMP parrallel region,
  ! and these buffers have to be thread-shared
  type (ghostBuffer_t)   :: ghostbuf,ghostbuf_cv
  type (ghostBuffer3d_t) :: ghostbuf3d
  type (edgeBuffer_t)    :: edge1

  public  :: testchecksum, test_ghost
  private :: genchecksum


contains
  !===============================================================================
  subroutine testchecksum(elem, par,GridEdge)
    use element_mod, only : element_t
    use parallel_mod, only : parallel_t, iam, syncmp
    use gridgraph_mod, only : gridedge_t,printchecksum
    use edge_mod, only : initedgebuffer, edgevpack, edgevunpack, edgedgvpack, &
         edgedgvunpack, edgebuffer_t
    use bndry_mod, only : bndry_exchangev
    use kinds, only : real_kind
    use schedule_mod, only : schedule_t, schedule, checkschedule
    use dimensions_mod, only : np, nlev, nelem, nelemd

    implicit none
!#include <stats.h>
    type (element_t), intent(in) :: elem(:)
    !================Arguments===============================
    type (parallel_t),   intent(in)           :: par
    type (GridEdge_t),   intent(in),target    :: GridEdge(:)


    !==============Local Allocatables==========================================
    real (kind=real_kind),allocatable      :: TestPattern_g(:,:,:,:), &
         Checksum_g(:,:,:,:),    &
         TestPattern_l(:,:,:,:), &
         Checksum_l(:,:,:,:), &
         Checksum_dg_l(:,:,:,:)
    type (EdgeBuffer_t)                         :: buffer

    !==============Local Temporaries===========================================
    type (Schedule_t),   pointer                :: pSchedule
    logical                                     :: ChecksumError
    integer                                     :: ix,iy
    integer                                     :: numlev
    integer                                     :: il,ig,ie,k
    integer                                     :: i
    integer                                     :: kptr,ielem
    integer                                     :: nSendCycles,nRecvCycles
    !===================Local Parameters=======================================
    logical, parameter                          :: PrChecksum = .FALSE.
    logical, parameter                          :: VerboseTiming=.FALSE.
    !==========================================================================


    allocate(TestPattern_g(np,np,nlev,nelem))
    allocate(Checksum_g(np,np,nlev,nelem))

    call genchecksum(TestPattern_g,Checksum_g,GridEdge)


    ! Setup the pointer to proper Schedule
#ifdef _PREDICT
    pSchedule => Schedule(iam)
#else
    pSchedule => Schedule(1)
#endif


    allocate(TestPattern_l(np,np,nlev,nelemd))
    allocate(   Checksum_l(np,np,nlev,nelemd))
    allocate(Checksum_dg_l(0:np+1,0:np+1,nlev,nelemd))

    do il=1,nelemd
       Checksum_l(:,:,:,il) = 0.d0
       Checksum_dg_l(:,:,:,il) = 0.d0
       ig = pSchedule%Local2Global(il)
       TestPattern_l(:,:,:,il) = TestPattern_g(:,:,:,ig)
       Checksum_l(:,:,:,il)    = TestPattern_l(:,:,:,il)
    enddo


    if(PrChecksum) then
       print *,'testchecksum:  The GLOBAL pattern is: '
       call PrintCheckSum(TestPattern_g,Checksum_g)
    endif

    if(PrChecksum) then
       print *,'testchecksum:  The LOCAL pattern is: '
       call PrintCheckSum(TestPattern_l,Checksum_l)
    endif

    nSendCycles = pSchedule%nSendCycles
    nRecvCycles = pSchedule%nRecvCycles

    !=======================================
    !  Allocate the communication Buffers
    !=======================================

    call initEdgeBuffer(buffer,nlev)

    !=======================================
    !  Synch everybody up
    !=======================================

    call syncmp(par)

    !==================================================
    !   Pack up the communication buffer
    !==================================================

    do ie=1,nelemd
       kptr=0
       numlev=nlev
       call edgeVpack(buffer,TestPattern_l(1,1,1,ie),numlev,kptr,elem(ie)%desc)
    enddo
    print *,'testchecksum: after call to edgeVpack'

    !==================================================
    !  Perform the boundary exchange
    !==================================================

    call bndry_exchangeV(par,buffer)

    !==================================================
    !   UnPack and accumulate the communication buffer
    !==================================================

    do ie=1,nelemd
       kptr   = 0
       numlev = nlev
       call edgeVunpack(buffer,Checksum_l(1,1,1,ie),numlev,kptr,elem(ie)%desc)
    enddo
    print *,'testchecksum: after call to edgeVunpack'

    !=====================================
    !  Printout the distributed checksum
    !=====================================

    if(PrChecksum) then
       call PrintCheckSum(TestPattern_l,Checksum_l)
    endif

    !============================================================
    !  Verify the distributed checksum against the serial version
    !============================================================
    do i=1,nelemd
       ig       = pSchedule%Local2Global(i)
       ChecksumError=.FALSE.
       do k=1,nlev
          do iy=1,np
             do ix=1,np
                if(Checksum_l(ix,iy,k,i) .ne. Checksum_g(ix,iy,k,ig)) then
                   write(*,100) iam , INT(TestPattern_l(ix,iy,k,i)),   &
                        INT(Checksum_g(ix,iy,k,ig)), INT(Checksum_l(ix,iy,k,i))
                   ChecksumError=.TRUE.
                endif
             enddo
          enddo
       enddo
       if(PrChecksum) then
          if(.NOT. ChecksumError) print *, 'IAM: ',iam,'testchecksum: Element', &
               pSchedule%Local2Global(i),'Verified'
       endif
    enddo
    ! =========================================
    ! Perform checksum for DG boundary exchange
    ! =========================================

    !==================================================
    !   Pack up the communication buffer
    !==================================================

    print *,'testchecksum: before call to edgeVpack'
    do ielem=1,nelemd
       kptr=0
       numlev=nlev
       call edgeDGVpack(buffer,Checksum_l(1,1,1,ielem),numlev,kptr,elem(ielem)%desc)
    enddo
    print *,'testchecksum: after call to edgeVpack'

    !==================================================
    !  Perform the boundary exchange
    !==================================================

    call bndry_exchangeV(par,buffer)

    !==================================================
    !   UnPack and accumulate the communication buffer
    !==================================================

    print *,'testchecksum: before call to edgeVunpack'
    do ielem=1,nelemd
       kptr   = 0
       numlev = nlev
       call edgeDGVunpack(buffer,Checksum_dg_l(0,0,1,ielem),numlev,kptr,elem(ielem)%desc)
    enddo
    print *,'testchecksum: after call to edgeDGVunpack'

    !==================================================
    !   Check Correctness for DG communication
    !==================================================
    do ie=1,nelemd
       ig       = pSchedule%Local2Global(ie)
       ChecksumError=.FALSE.
       do k=1,nlev
          do i=1,np

             ! =========================
             ! Check South Flux
             ! =========================
             if(Checksum_l(i,1,k,ie) .ne. Checksum_dg_l(i,0,k,ie)) then
                write(*,100) iam , INT(TestPattern_l(i,1,k,ie)), &
                     INT(Checksum_g(i,1,k,ig)), INT(Checksum_l(i,1,k,ie))
                ChecksumError=.TRUE.
             endif

             ! =========================
             ! Check East Flux
             ! =========================
             if(Checksum_l(np,i,k,ie) .ne. Checksum_dg_l(np+1,i,k,ie)) then
                write(*,100) iam , INT(TestPattern_l(np,i,k,ie)), &
                     INT(Checksum_g(np,i,k,ig)), INT(Checksum_l(np,i,k,ie))
                ChecksumError=.TRUE.
             endif

             ! =========================
             ! Check North Flux
             ! =========================
             if(Checksum_l(i,np,k,ie) .ne. Checksum_dg_l(i,np+1,k,ie)) then
                write(*,100) iam , INT(TestPattern_l(i,np,k,ie)), &
                     INT(Checksum_g(i,np,k,ig)), INT(Checksum_l(i,np,k,ie))
                ChecksumError=.TRUE.
             endif

             ! =========================
             ! Check West Flux
             ! =========================
             if(Checksum_l(1,i,k,ie) .ne. Checksum_dg_l(0,i,k,ie)) then
                write(*,100) iam , INT(TestPattern_l(1,i,k,ie)), &
                     INT(Checksum_g(1,i,k,ig)), INT(Checksum_l(1,i,k,ie))
                ChecksumError=.TRUE.
             endif

          enddo
       enddo
       if(PrChecksum) then
          if(.NOT. ChecksumError) print *, 'IAM: ',iam,'testchecksum: Element', &
               pSchedule%Local2Global(ie),'Verified'
       endif
    enddo




    !mem  print *,'testchecksum: before call to deallocate'
    deallocate(buffer%buf)
    deallocate(TestPattern_g)
    deallocate(Checksum_g)
    deallocate(TestPattern_l)
    deallocate(Checksum_l)
    deallocate(Checksum_dg_l)

    call CheckSchedule()


100 format('IAM:',i3,'testchecksum: Error with checksum:',I10,'--->',I10,' !=',I10)

  end subroutine testchecksum

  subroutine genchecksum(TestField,Checksum,GridEdge)
    use kinds, only : real_kind
    use gridgraph_mod, only : gridedge_t, edgeindex_t

    implicit none

    type (GridEdge_t), intent(in),target   :: GridEdge(:)
    real(kind=real_kind), intent(inout),target  :: TestField(:,:,:,:)
    real(kind=real_kind), intent(inout),target  :: Checksum(:,:,:,:)

    type (EdgeIndex_t), pointer :: sIndex,dIndex
    type (GridEdge_t), pointer  :: gEdge

    integer                     :: i,ix,iy,k,iedge
    integer                     :: nelem_edge,nwords
    integer                     :: idest,isrc
    logical,parameter           :: PrChecksum =.FALSE.
#ifdef TESTGRID
    nelem_edge = SIZE(GridEdge)

    !  Setup the test pattern
    do i=1,nelem
       do k=1,nlev
          do iy=1,np
             do ix=1,np
                TestField(ix,iy,k,i)= ix + (iy-1)*np + 100*i + 1000000*(k-1)
             enddo
          enddo
       enddo
    enddo

    !  Initalize the checksum array to be the test pattern
    Checksum=TestField

    !  Calculate the checksum
    do iedge = 1,nelem_edge
       gEdge  => GridEdge(iedge)
       isrc   =  gEdge%tail%number
       idest  =  gEdge%head%number
       nwords =  gEdge%wgtV
       sIndex => gEdge%TailIndex
       dIndex => gEdge%HeadIndex


       do k=1,nlev
          do i=1,nwords
             Checksum(GridEdge(iedge)%HeadIndex%ixV(i),GridEdge(iedge)%HeadIndex%iyV(i),k,idest) = &
                  Checksum(GridEdge(iedge)%HeadIndex%ixV(i),GridEdge(iedge)%HeadIndex%iyV(i),k,idest) +  &
                  TestField(GridEdge(iedge)%TailIndex%ixV(i),GridEdge(iedge)%TailIndex%iyV(i),k,isrc)
          enddo
       enddo

    enddo

#if 0
    if(PrChecksum .AND. par%masterproc) then
       call PrintChecksum(TestField,Checksum)
    endif
#endif

100 format(1x,'element=',I2)
110 format(1x,I4,2x,'checksum = ',I5)
#endif
  end subroutine genchecksum






  subroutine test_ghost(hybrid,elem,nets,nete)
!
! MT 7/2011  test new ghost exchange to populate ghost cells
!
  use kinds, only : real_kind
  use parallel_mod, only : syncmp
  use dimensions_mod, only : np, nlev, max_corner_elem
  use hybrid_mod, only : hybrid_t
  use element_mod, only : element_t
  use bndry_mod, only : ghost_exchangevfull, bndry_exchangev, ghost_exchangev3d
  use edge_mod, only : ghostbuffer_t, ghostvpackfull, ghostvunpackfull, initghostbufferfull,&
       freeghostbuffer, edgebuffer_t, freeedgebuffer, initedgebuffer,&
       edgevpack,edgevunpack, &
       initghostbuffer3d, ghostvpack3d, ghostvunpack3d

  use control_mod, only : north,south,east,west,neast, nwest, seast, swest

  implicit none

  integer, parameter :: NHC=1

  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  integer :: nets,nete

  real (kind=real_kind) :: pin (        np,          np,  nlev,  nets:nete)  !CE: fvm tracer
  real (kind=real_kind) :: pout(1-np:np+np,  1-np:np+np,  nlev,  nets:nete)  !CE: fvm tracer


  integer :: mult(5:8)
  real (kind=real_kind) :: vout(1-NHC : np+NHC,   1-NHC : np+NHC, nlev, nets:nete)  !CE: fvm tracer

  real (kind=real_kind) :: sw  (1-NHC : 1     ,   1-NHC : 1       , nlev, max_corner_elem-1, nets:nete)
  real (kind=real_kind) :: se  (   np : np+NHC,   1-NHC : 1       , nlev, max_corner_elem-1, nets:nete)
  real (kind=real_kind) :: ne  (   np : np+NHC,      np : np+NHC  , nlev, max_corner_elem-1, nets:nete)
  real (kind=real_kind) :: nw  (1-NHC : 1     ,      np : np+NHC  , nlev, max_corner_elem-1, nets:nete)

  real (kind=real_kind) :: cin(2,2,nlev,nets:nete)  !CE: fvm tracer
  real (kind=real_kind) :: cout(-1:4,-1:4,nlev,nets:nete)  !CE: fvm tracer
  integer :: i,j,k,m,n,ie,kptr,np1,np2,nc,nc1,nc2
  logical :: fail,fail1,fail2
  real (kind=real_kind) :: tol=.1
  call syncmp(hybrid%par)
  if (hybrid%par%masterproc) print *,'testing ghost exchange'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! first test on the Gauss Grid with same number of ghost cells:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nc=2   ! test using GLL interior points
  nc1=-1
  nc2=4

  np1=1-np
  np2=np+np

  call initghostbufferfull(ghostbuf,   nlev,np)
  call initghostbufferfull(ghostbuf_cv,nlev,nc)
  call initghostbuffer3d(ghostbuf3d,   nlev,np, NHC)

  do k=1,nlev
  do ie=nets,nete
     !pin(:,:,1,ie)=elem(ie)%spherep(:,:)%lat + elem(ie)%spherep(:,:)%lon
     pin(:,:,k,ie)=elem(ie)%gdofp(:,:)
     cin(1,1,k,ie)=pin(1,1,k,ie)
     cin(nc,nc,k,ie)=pin(np,np,k,ie)
     cin(1,nc,k,ie)=pin(1,np,k,ie)
     cin(nc,1,k,ie)=pin(np,1,k,ie)
  enddo
  enddo
  pout=0
  cout=0

#if 0
  ! DSS pin, to make all edge points agree exactly:
  call initedgebuffer(edge1,nlev)
  do ie=nets,nete
     kptr=0
     do k=1,nlev
        pin(:,:,k,ie)=pin(:,:,k,ie)*elem(ie)%spheremp(:,:)
     enddo
     call edgeVpack(edge1, pin(:,:,:,ie),nlev,kptr,elem(ie)%desc)
  end do
  call bndry_exchangeV(hybrid,edge1)
  do ie=nets,nete
     kptr=0
     call edgeVunpack(edge1, pin(:,:,:,ie), nlev, kptr, elem(ie)%desc)
     do k=1,nlev
        pin(:,:,k,ie)=pin(:,:,k,ie)*elem(ie)%rspheremp(:,:)
     enddo
  enddo
#endif



  cout=0
  do ie=nets,nete
     kptr=0
     call ghostVpackfull(ghostbuf, pin(:,:,:,ie),1,np,np,nlev,kptr,elem(ie)%desc)
     call ghostVpackfull(ghostbuf_cv, cin(:,:,:,ie),1,nc,nc,nlev,kptr,elem(ie)%desc)
  end do

  call ghost_exchangeVfull(hybrid,ghostbuf,np)
  call ghost_exchangeVfull(hybrid,ghostbuf_cv,nc)

  !call syncmp(hybrid%par)
  do ie=nets,nete
     kptr=0
     call ghostVunpackfull(ghostbuf, pout(:,:,:,ie), np1,np2,np,nlev, kptr, elem(ie)%desc)
     call ghostVunpackfull(ghostbuf_cv, cout(:,:,:,ie), nc1,nc2,nc,nlev, kptr, elem(ie)%desc)
  enddo


  do ie=nets,nete
     ! check interior was not changed:
     fail=.false.
     do i=1,np
        do j=1,np
           if (pout(i,j,1,ie) /= 0) fail=.true.
        enddo
     enddo
     ! check i=1,nc edges
     do j=1,np
        if (abs(pin(1,j,1,ie)-pout(0,j,1,ie)) .gt. tol )  fail=.true.
        if (abs(pin(np,j,1,ie) - pout(np+1,j,1,ie)).gt.tol) fail=.true.
     enddo
     ! check j=1,np edges
     do i=1,np
        if (abs(pin(i,1,1,ie)-pout(i,0,1,ie)) .gt. tol )  fail=.true.
        if (abs(pin(i,np,1,ie) - pout(i,np+1,1,ie)).gt.tol) fail=.true.
     enddo

!       nc +--------+
!        ^ | nw  ne |
!     j  | |        |
!        1 | sw  se |
!          +--------+
!           1 --> nc
!              i

     ! check SW corner ghost cells
     if ( elem(ie)%desc%putmapP_ghost(swest) /= -1) then
        ! i=1   vs i=0    j=np1,0
        ! j=1   vs j=0    i=np1,0
        do i=np1,0
           if (abs(pout(i,1,1,ie)-pout(i,0,1,ie)) .gt. tol )  fail=.true.
           if (abs(pout(1,i,1,ie)-pout(0,i,1,ie)).gt.tol) fail=.true.
        enddo

        ! check that cout maches pout corners
        if (abs(cout(nc1,nc1,1,ie)-pout(np1,np1,1,ie)).gt.tol) fail=.true.
        if (abs(cout(nc1+1,nc1,1,ie)-pout(np1+np-1,np1,1,ie)).gt.tol) fail=.true.
        if (abs(cout(nc1,nc1+1,1,ie)-pout(np1,np1+np-1,1,ie)).gt.tol) fail=.true.
        if (abs(cout(nc1+1,nc1+1,1,ie)-pout(np1+np-1,np1+np-1,1,ie)).gt.tol) fail=.true.

     endif

     !check SE corner
     if ( elem(ie)%desc%putmapP_ghost(seast) /= -1) then
        do i=np+1,np2
           if (abs(pout(i,1,1,ie)-pout(i,0,1,ie)) .gt. tol )  fail=.true.
        enddo
        do j=np1,0
           if (abs(pout(np+1,j,1,ie)-pout(np,j,1,ie)).gt.tol) fail=.true.
        enddo

        ! check that cout maches pout corners
        if (abs(cout(nc2,nc1,1,ie)-pout(np2,np1,1,ie)).gt.tol) fail=.true.
        if (abs(cout(nc2-1,nc1,1,ie)-pout(np2-np+1,np1,1,ie)).gt.tol) fail=.true.
        if (abs(cout(nc2,nc1+1,1,ie)-pout(np2,np1+np-1,1,ie)).gt.tol) fail=.true.
        if (abs(cout(nc2-1,nc1+1,1,ie)-pout(np2-np+1,np1+np-1,1,ie)).gt.tol) fail=.true.

     endif


     !check NE corner
     if ( elem(ie)%desc%putmapP_ghost(neast) /= -1) then
        do i=np+1,np2
           if (abs(pout(i,np,1,ie)-pout(i,np+1,1,ie)) .gt. tol )  fail=.true.
        enddo
        do j=np+1,np2
           if (abs(pout(np,j,1,ie)-pout(np+1,j,1,ie)).gt.tol) fail=.true.
        enddo

        ! check that cout maches pout corners
        if (abs(cout(nc2,nc2,1,ie)-pout(np2,np2,1,ie)).gt.tol) fail=.true.
        if (abs(cout(nc2-1,nc2,1,ie)-pout(np2-np+1,np2,1,ie)).gt.tol) fail=.true.
        if (abs(cout(nc2,nc2-1,1,ie)-pout(np2,np2-np+1,1,ie)).gt.tol) fail=.true.
        if (abs(cout(nc2-1,nc2-1,1,ie)-pout(np2-np+1,np2-np+1,1,ie)).gt.tol) fail=.true.

     endif

     !check NW corner
     if ( elem(ie)%desc%putmapP_ghost(nwest) /= -1) then
        do i=np1,0
           if (abs(pout(i,np,1,ie)-pout(i,np+1,1,ie)) .gt. tol )  fail=.true.
        enddo
        do j=np+1,np2
           if (abs(pout(1,j,1,ie)-pout(0,j,1,ie)).gt.tol) fail=.true.
        enddo

        ! check that cout maches pout corners
        if (abs(cout(nc1,nc2,1,ie)-pout(np1,np2,1,ie)).gt.tol) fail=.true.
        if (abs(cout(nc1+1,nc2,1,ie)-pout(np1+np-1,np2,1,ie)).gt.tol) fail=.true.
        if (abs(cout(nc1,nc2-1,1,ie)-pout(np1,np2-np+1,1,ie)).gt.tol) fail=.true.
        if (abs(cout(nc1+1,nc2-1,1,ie)-pout(np1+np-1,np2-np+1,1,ie)).gt.tol) fail=.true.

     endif





     !if (hybrid%par%masterproc .and. ie==nets) fail=.true.
     ! only print masterproc:
     !if ( .not. hybrid%par%masterproc ) fail=.false.

     if (fail) then
           print *,'ie=',ie
           do j=np,1,-1
           !   write(*,'(99f9.4)') (pin(i,j,1,ie),i=1,np),(pout(i,j,1,ie),i=1,np)
           enddo
           print *
           print *,'neast,nwest,seast,swest',&
                elem(ie)%desc%putmapP_ghost(neast), &
                elem(ie)%desc%putmapP_ghost(nwest), &
                elem(ie)%desc%putmapP_ghost(seast), &
                elem(ie)%desc%putmapP_ghost(swest)
           print *,'reverse:'
           print *,elem(ie)%desc%reverse(nwest),elem(ie)%desc%reverse(north),elem(ie)%desc%reverse(neast)
           print *,elem(ie)%desc%reverse(west),' ',elem(ie)%desc%reverse(east)
           print *,elem(ie)%desc%reverse(swest),elem(ie)%desc%reverse(south),elem(ie)%desc%reverse(seast)

#if 1
           do j=np+np,np+1,-1
              write(*,'(99f9.0)') (pout(i,j,1,ie),i=np1,np2)
           enddo
           do j=np,1,-1
              write(*,'(99f9.0)') (pout(i,j,1,ie),i=1-np,0),(pin(i,j,1,ie),i=1,np),(pout(i,j,1,ie),i=np+1,np+np)
           enddo
           do j=0,1-np,-1
              write(*,'(99f9.0)') (pout(i,j,1,ie),i=np1,np2)
           enddo
#endif
           print *,'CV grid:'
           do j=nc+nc,nc+1,-1
              write(*,'(99f9.0)') (cout(i,j,1,ie),i=nc1,nc2)
           enddo
           do j=nc,1,-1
              write(*,'(99f9.0)') (cout(i,j,1,ie),i=1-nc,0),(cin(i,j,1,ie),i=1,nc),(cout(i,j,1,ie),i=nc+1,nc+nc)
           enddo
           do j=0,1-nc,-1
              write(*,'(99f9.0)') (cout(i,j,1,ie),i=nc1,nc2)
           enddo
           stop 'ERROR:  ghost exchange failed consistency check 0'
     endif
  enddo


#if 1
  cout=0
  vout=0
  do ie=nets,nete
     kptr=0
     call ghostVpack3d(ghostbuf3d, pin(:,:,:,ie), nlev, kptr, elem(ie)%desc)
  end do

  call ghost_exchangeV3d(hybrid, ghostbuf3d)

  do ie=nets,nete
     kptr=0
     call ghostVunpack3d(ghostbuf3d, vout(:,:,:,ie), nlev, kptr, elem(ie)%desc, sw,se,nw,ne,mult)
  enddo

  if (hybrid%masterthread) print *, "Checking results of ghostVunpack3d against ghostVunpackfull."

! real (kind=real_kind) :: pout(1-np  : np+np,  1-np  : np+np,  nlev, nets:nete)
! real (kind=real_kind) :: vout(1-NHC : np+NHC, 1-NHC : np+NHC, nlev, nets:nete)

  fail=.false.
  do ie=nets,nete
    do i=1,np
      do j=1,np
        do k=1,nlev
          if (pout(i,j,k,ie) /= vout(i,j,k,ie)) fail=.true.
        enddo
      enddo
    enddo
  enddo
  if (fail) then
    print *, ' Comparison of interrior ghostVunpackfull and ghostVunpack3d fialed.'
    stop 'ERROR:  ghost exchange failed consistency check 1'
  endif
  fail=.false.
  do ie=nets,nete
    do i=1-NHC, np+NHC
      if (i.lt.1 .or. np.lt.i) then
        do j=1-NHC, np+NHC
          if (j.lt.1 .or. np.lt.j) then
            do k=1,nlev
              m=i+1
              n=j+1
              if (i.lt.1) m=i-1
              if (j.lt.1) n=j-1
              if (pout(m,n,k,ie) /= vout(i,j,k,ie) .and. pout(m,n,k,ie) /= 0.0) fail=.true.
              if (pout(m,n,k,ie) /= vout(i,j,k,ie) .and. pout(m,n,k,ie) /= 0.0) then
                 print *," Error in Checksum:",i,j,m,n,pout(m,n,k,ie),vout(i,j,k,ie), fail, "  ************ "
              endif
            enddo
          endif
        enddo
      endif
    enddo
  enddo
  if (fail) then
    print *, ' Comparison of exterior ghostVunpackfull and ghostVunpack3d fialed.'
    stop 'ERROR:  ghost exchange failed consistency check 2'
  endif
  if (hybrid%masterthread) print *, "Check passed: Results of ghostVunpack3d and ghostVunpackfull match."

#endif
  call freeghostbuffer(ghostbuf)
  call freeghostbuffer(ghostbuf_cv)
  call syncmp(hybrid%par)
  if (hybrid%par%masterproc) print *,'done'
end subroutine

end module checksum_mod
