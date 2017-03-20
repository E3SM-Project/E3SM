!===============================================================================
! SVN $Id: smooth_mod.F90 56089 2013-12-18 00:50:07Z mlevy@ucar.edu $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/tools/mapping/trunk_tags/mapping_141106/gen_mapping_files/runoff_to_ocn/src/smooth_mod.F90 $
!===============================================================================

MODULE smooth_mod

#define _NEW 1
#define _BREADTH 1
   use map_mod

   implicit none

!-------------------------------------------------------------------------------
!     module variables
!-------------------------------------------------------------------------------

! from map_mod:
!   real,parameter   :: rEarth   =  6.37122e6    ! radius of earth ~ m
!   real,parameter   :: DtoR     =  3.14159265358979323846/180.0
   real(r8),parameter   :: DtoR = DEGtoRAD
   real(r8),allocatable :: garr(:,:,:)   ! dummy runoff (nx,ny,nbasin) in kg/s/m^2
   integer, parameter :: maxLinear = 120000

!===============================================================================
CONTAINS
!===============================================================================

SUBROUTINE smooth_init(map_in, map_out)

   implicit none

   !--- arguments ---
   type(sMatrix),intent( in)   :: map_in   ! original unsmoothed, matrix
   type(sMatrix),intent(inout) :: map_out  ! smoothing matrix

   !--- local ---
   integer :: i,j,n ! indicies: row, col, sparse matrix
   integer :: rcode ! return code
   integer :: ns    ! number of links
   integer :: nactive
   integer :: jmd_count

   integer,allocatable :: i_count(:)

   character(*),parameter :: subName = "(smooth_init) "

!-------------------------------------------------------------------------------
! PURPOSE:
! o Given map_in create map_out which maps from domain_b of
!   map_in to itself.  This is a template for the square smoothing
!   map, which we assume is needed for domain_b of map_in.
! Notes:
! o Deallocation below copied from map_dup routine.  Depending on the context,
!   this routine might be given a map_out that has already been allocated, in
!   which case it has to be deallocated before being (re) allocated
! o the number of non-zero rows per column depends on the resolution of the
!   grid and the radius of smoothing (a nml parameter), so it's hard to guess
!   how big this smoothing matrix might be.  To date (March 2012), the goal at
!   NCAR has been to smooth each src cell over, say, 10-40 cells. A variety of
!   So a good guess might be, say, 40 x (the number of src cells).
!   But typically some cells, where the grid is finer, may get smooth over many
!   more, and the nml radius might be set much larger than typically desired.
!   Thus the current guess is close (2 ^ 31 ) - 1
!   which is the maximum size that can be indexed with a 4-byte integer.
!-------------------------------------------------------------------------------

   !------------------------------------------------
   ! de-allocate space (if already allocated)
   !------------------------------------------------
   deallocate(map_out%  xc_a,STAT=rcode)
   deallocate(map_out%  yc_a,STAT=rcode)
   deallocate(map_out%  xv_a,STAT=rcode)
   deallocate(map_out%  yv_a,STAT=rcode)
   deallocate(map_out%mask_a,STAT=rcode)
   deallocate(map_out%area_a,STAT=rcode)

   deallocate(map_out%  xc_b,STAT=rcode)
   deallocate(map_out%  yc_b,STAT=rcode)
   deallocate(map_out%  xv_b,STAT=rcode)
   deallocate(map_out%  yv_b,STAT=rcode)
   deallocate(map_out%mask_b,STAT=rcode)
   deallocate(map_out%area_b,STAT=rcode)

   deallocate(map_out%frac_a,STAT=rcode)
   deallocate(map_out%frac_b,STAT=rcode)

   deallocate(map_out%s     ,STAT=rcode)
   deallocate(map_out%row   ,STAT=rcode)
   deallocate(map_out%col   ,STAT=rcode)
   deallocate(map_out%sn1   ,STAT=rcode)
   deallocate(map_out%sn2   ,STAT=rcode)

   !------------------------------------------------
   ! allocate space for smoothing matrix (must guess its size)
   !------------------------------------------------
   allocate(i_count(map_in%n_b))
   i_count = 0
   do n=1,map_in%n_s
      i = map_in%row(n)
      if (map_in%s(n) > 0.0) i_count(i) = i_count(i) + 1
   end do
   nactive = count(i_count > 0) ! number of ocn cells that might get mapped to
   deallocate(i_count)
!  nactive = count(map_in%mask_b > 0) ! number of active ocn cells
!  nactive = map_in%n_s               ! number of ocn cells receiving runoff (an upper bound)
!  ns = nactive*(map_in%n_b/(20*20)) ! (about 1000 miles/1000 miles smoothing) guess a sufficient size for smoothing matrix
!  ns = nactive*10000
!  ns = nactive*1000
   ns = 2147000000  ! max 4-byte int
   write(*,*) subName,' orig matrix a-domain n active: ',count(map_in%mask_a > 0)
   write(*,*) subName,' orig matrix b-domain n active: ',count(map_in%mask_b > 0)
   write(*,*) subName,' orig matrix a-domain size    : ',map_in%n_a
   write(*,*) subName,' orig matrix b-domain size    : ',map_in%n_b
   write(*,*) subName,' orig matrix            ns    : ',map_in%n_s
   write(*,*) subName,' smoothing matrix guess ns    : ',ns
   write(*,*) subName,' nactive                      = ',nactive
!  write(*,*) subName,' map_in%n_b                   = ',map_in%n_b

   allocate(map_out%  xc_a(            map_in%n_b) )
   allocate(map_out%  yc_a(            map_in%n_b) )
   allocate(map_out%  xv_a(map_in%nv_b,map_in%n_b) )
   allocate(map_out%  yv_a(map_in%nv_b,map_in%n_b) )
   allocate(map_out%mask_a(            map_in%n_b) )
   allocate(map_out%area_a(            map_in%n_b) )

   allocate(map_out%  xc_b(            map_in%n_b) )
   allocate(map_out%  yc_b(            map_in%n_b) )
   allocate(map_out%  xv_b(map_in%nv_b,map_in%n_b) )
   allocate(map_out%  yv_b(map_in%nv_b,map_in%n_b) )
   allocate(map_out%mask_b(            map_in%n_b) )
   allocate(map_out%area_b(            map_in%n_b) )

   allocate(map_out%frac_a(map_in%n_b) )
   allocate(map_out%frac_b(map_in%n_b) )

   allocate(map_out%s  (ns))
   allocate(map_out%row(ns))
   allocate(map_out%col(ns))
   allocate(map_out%sn1(map_in%n_b) )
   allocate(map_out%sn2(map_in%n_b) )

   !------------------------------------------------
   ! set values
   !------------------------------------------------
   map_out%   n_a = map_in%   n_b
   map_out%dims_a = map_in%dims_b
   map_out%  ni_a = map_in%  ni_b
   map_out%  nj_a = map_in%  nj_b
   map_out%  nv_a = map_in%  nv_b
   map_out%  xc_a = map_in%  xc_b
   map_out%  yc_a = map_in%  yc_b
   map_out%  xv_a = map_in%  xv_b
   map_out%  yv_a = map_in%  yv_b
!  map_out%mask_a = map_in%mask_b ! all active ocn cells
   map_out%area_a = map_in%area_b

   !--- compute minimal src domain mask for smoothing matrix ---
   jmd_count = 0
   map_out%mask_a = 0
   do n=1,map_in%n_s
      i = map_in%row(n)     ! this ocn cell could get runoff
      map_out%mask_a(i) = 1 ! this ocn cell's runoff get's smoothed
      if(map_in%s(n) > 0.0) then
        jmd_count = jmd_count+1
      endif
   end do
   write(*,*) subName,'number of source points is =  ',jmd_count
   write(*,*) subName,'map_in%ns                  =  ',map_in%n_s

   map_out%   n_b = map_in%   n_b
   map_out%dims_b = map_in%dims_b
   map_out%  ni_b = map_in%  ni_b
   map_out%  nj_b = map_in%  nj_b
   map_out%  nv_b = map_in%  nv_b
   map_out%  xc_b = map_in%  xc_b
   map_out%  yc_b = map_in%  yc_b
   map_out%  xv_b = map_in%  xv_b
   map_out%  yv_b = map_in%  yv_b
   map_out%mask_b = map_in%mask_b
   map_out%area_b = map_in%area_b

!  map_out%frac_a = map_in%frac_b
!  map_out%frac_b = map_in%frac_b
   map_out%frac_a = 1.0
   map_out%frac_b = 1.0

   map_out%n_s    = ns
   map_out%s      = 1.0
   map_out%row    = 1
   map_out%col    = 1
   map_out%sn1    = map_in%sn1
   map_out%sn2    = map_in%sn2

   map_out%title      = "CCSM conservative smoothing map"
   map_out%normal     = map_in%normal
   map_out%method     = "created using SVN $Id: smooth_mod.F90 56089 2013-12-18 00:50:07Z mlevy@ucar.edu $"
   map_out%history    = map_in%history
   map_out%convention = map_in%convention
   map_out%domain_a   = map_in%domain_b
   map_out%domain_b   = map_in%domain_b

END SUBROUTINE smooth_init

!===============================================================================

SUBROUTINE smooth(map,efold,rmax)

   implicit none

   !--- arguments ---
   type(sMatrix), intent(inout) :: map
   real(kind=r8), intent(in) :: efold                    ! efold scale (m)
   real(kind=r8), intent(in) :: rmax                     ! max smoothing distance (m)

   !--- local ---
   integer         :: n ! loop over matrix elements
   integer         :: i,j,i0,i1,j0,j1,ic,jc,itmp,jtmp,k
   integer         :: ii,jj,ip,np,ni,nj,n_s,ngood,nbox
   integer         :: itest,jtest,k1,k2,k3,k4,gdcnt,cnt
   integer         :: mingd,maxgd,avggd,mingdi,mingdj,maxgdi,maxgdj
   integer         :: minbx,maxbx,minbxi,minbxj,maxbxi,maxbxj
   integer         :: onegd,nactivea,nactiveb,ntot
   integer, allocatable :: iind(:),jind(:),imask(:,:)
   integer, allocatable :: imask_jmd(:,:)
   real(r8), allocatable :: rdist(:,:),areaa(:,:),wgt(:,:),wgt2(:,:)
   real(r8), allocatable :: diff(:,:)
   real(r8), allocatable :: s(:),row(:),col(:)
   real(r8) :: wgtsum
   real(r8), parameter :: rmiss = 1.e30

   integer :: length
   integer, allocatable :: i2ind(:),j2ind(:)
   integer, allocatable :: indxLinear(:,:)
   integer :: t00
   integer :: i2,j2

   integer,allocatable :: nDest(:)    ! diag: # of cells (size of smoothed footprint)
   integer :: minDest,maxDest,avgDest ! diag: min,max,avg # of cells (in smoothed footprint)
   integer :: nDest10, nDest100, nDest1000, nDest10000 ! for diag of size of footprint

   integer         :: nMod  ! for time-stamp progress report
   character( 8)   :: dstr  ! wall clock date
   character(10)   :: tstr  ! wall clock time

   !--- formats ---
   character(*),parameter :: F00 = "('(smooth) ',4a)"
   character(*),parameter :: F1  = "('(smooth) ',a,2i11)"
   character(*),parameter :: F2  = "('(smooth) ',a,2F11.6)"
   character(*),parameter :: F3  = "('(smooth) ',a,es18.7)"
   character(*),parameter :: F6  = "('(smooth) ',a,3i6,es18.7,i6)"
   character(*),parameter :: F9  = "('(smooth) ',a,i6,' @  (',i8,',',i8,')')"
   character(*),parameter :: F12 = "('(smooth) ',1x,a4,2('-',a2),2x,a2,2(':',a2),' j=',i9,', ',f7.3,'% complete')"

   character(*),parameter :: subName = "(smooth) "

   integer :: level
   integer :: strPtr
   integer :: kStart
!-------------------------------------------------------------------------------
! PURPOSE:
! o Given a square sMatrix, map, generate the links for local smoothing
!   of input fluxes with following properties:
!       -must be conservative
!       -must redistribute onto active points only
!       -uses a 2D Gaussian defined by efold & rmax for weights
!       -weights are a function of shortest active distance, so that
!        weights are zero across land isthmuses
!
!             -------------------------------------------------
!             |               |               |               |
!             |               |               |               |
!             |               |               |               |
!             |   ip=2        |    ip=3       |    ip=4       |
!    j1       |               |               |       ^       |
!             |               |               |       |       |
!             |               |               |       |d2     |
!             |               |               |       |       |
!             ----------------------------------------|--------
!             |               |               |       |       |
!             |               |               |       |       |
!             |   ip=1        |   (ic,jc)     |    (ii,jj)    |
!             |               |               |d1     |       |  d3
!             |               |    source<------------x-------|------->
!             |               |    cell       |       |       |
!             |               |               |    ip=5       |
!             |               |               |       |       |
!             ----------------------------------------|--------
!             |               |               |       |d4     |
!             |               |               |       |       |
!             |   ip=8        |    ip=7       |    ip=6       |
!    j0       |               |               |       |       |
!             |               |               |       v       |
!             |               |               |               |
!             |               |               |               |
!             |               |               |               |
!             -------------------------------------------------
!
!                    i0                              i1
!
!-------------------------------------------------------------------------------

   ni = map%ni_a
   nj = map%nj_a
   allocate(iind(map%n_a))
   allocate(jind(map%n_a))
   allocate(wgt  (ni,nj))
   allocate(wgt2 (ni,nj))
   allocate(diff (ni,nj))
   allocate(rdist(ni,nj))
   allocate(areaa(ni,nj))
   allocate(imask(ni,nj))
   allocate(imask_JMD(ni,nj))
   do j=1,map%n_a
      iind(j) = mod(j-1,ni)+1
      jind(j) = int((j-1)/ni)+1
      imask(iind(j),jind(j)) = map%mask_b(j)
      areaa(iind(j),jind(j)) = map%area_a(j)
   enddo



   length = 0
   allocate(i2ind(maxLinear))
   allocate(j2ind(maxLinear))
   allocate(indxLinear(ni,nj))
   i2ind=0;j2ind=0
#ifdef _NEW
   do i=1,map%n_b
       ii = iind(i)
       jj = jind(i)
       indxLinear(ii,jj) = i
   end do
#endif


   itest = 70
   jtest = 87

   n_s = 0
   minbx = 1000
   maxbx = 0
   mingd = 1000
   maxgd = 0
   avggd = 0
   onegd = 0
   write(6,F1) 'ni = ',ni
   write(6,F1) 'nj = ',nj
   write(6,F1) 'n_a = ',map%n_a
   write(6,F1) 'n_b = ',map%n_b
   write(6,F1) 'n_s = ',map%n_s
   write(6,F2) 'range(xc_a) = ',minval(map%xc_a(1:map%n_a)), &
  &     maxval(map%xc_a(1:map%n_a))
   write(6,F2) 'range(yc_a) = ',minval(map%yc_a(1:map%n_a)), &
  &     maxval(map%yc_a(1:map%n_a))
   write(6,F2) 'using efold (km) of = ',efold/1000.
   write(6,F2) 'using max radius (km) of = ',rmax/1000.

   !----------------------------------------------------------------------------
   ! compute smoothing matrix
   !----------------------------------------------------------------------------
   call shr_timer_get  (t00,subName//"compute smoothing matrix")
   call shr_timer_start(t00)

   !--- set up progress report print interval ---
   nMod = 1
   do while ( map%n_a > 500*nMod) ! want between 50 & 500 print statements
      nMod = 10*nMod
   end do
!  write(*,*) subName,"nMod = ",nMod
   nMod = 1000

   do j=1,map%n_a ! loop over all source points

     !--- progress report ---
     if (mod(j-1,nMod)==0) then
        call date_and_time(dstr,tstr)
        write(6,F12) dstr(1:4),dstr(5:6),dstr(7:8),tstr(1:2),tstr(3:4),tstr(5:6),j,100.0*float(j)/float(map%n_a)
     endif

     if (map%mask_a(j) /= 0) then ! only consider active source points
        !-------------------------------------------------------------
        ! find all cells within max radius of a given src cell j
        !-------------------------------------------------------------

        ic = iind(j)
        jc = jind(j)
        length = 1
        j2ind=0;i2ind=0

        rdist = 0.0
        wgt = 0.0
        wgt2 = 0.0
        wgtsum = 0.0
        nbox = 0
        gdcnt = 0
        kStart = 1

#ifdef _BREADTH
        i2ind(1) = ic
        j2ind(1) = jc
        imask_JMD = -1000
        where(imask ==  1) imask_JMD = -1
        level = 1
        strPtr=1
        imask_JMD(ic,jc) = 0
        call breadth_setDist(ni,nj,level,map%xc_a,map%yc_a,imask_JMD, &
                rdist,rmax,i2ind,j2ind, strPtr, length)

#else
        !--- recursive function to find dest cells ---

        call depth_setDist(ic,jc,ni,nj,map%xc_a,map%yc_a,0.0,imask,rdist,rmax,i2ind,j2ind,length)
#endif

!        print *,'smooth: after _setDist: length is: ',length-kStart+1
        if(length > maxLinear) then
           print *,'Error need to increase maxLinear length is:', length-kStart+1
        endif

#if _DBG
        print *,'source point: ',ic,jc
        do k=kStart,length
           i2 = i2ind(k)
           j2 = j2ind(k)
           print *,'neighbor points: (i,j): ',i2,j2,rdist(i2,j2)
        enddo
        stop 'smooth: After first call to _setDist'
#endif

#ifdef _NEW
       !-------------------------------------------------------------
       ! 1st pass computing smoothing matrix weights for column j
       !-------------------------------------------------------------
       wgtsum = 0.d0
       do k = kStart,length ! loop over destination points for input cell j
          ii = i2ind(k)
          jj = j2ind(k)
          wgt(ii,jj) = exp(-rdist(ii,jj)/efold)*areaa(ii,jj)
          wgtsum = wgtsum + wgt(ii,jj)
       enddo

       !-------------------------------------------------------------
       ! final pass computing smoothing matrix weights for column j
       !-------------------------------------------------------------
       do k=kStart,length ! loop over destination points for input cell j
          ii = i2ind(k)
          jj = j2ind(k)
          i2 = indxLinear(ii,jj)
          n_s = n_s + 1
          if (n_s > map%n_s) then
              write(6,F1)  'ERROR: smoothing matrix n_s > ',map%n_s
              write(6,F00) 'ERROR: initial guess of matrix size is too small'
              stop 'subName'
          end if
!         write(6,F1) 'n_s = ',n_s
          map%s(n_s) = (map%area_a(j)/map%area_b(i2))*(wgt(ii,jj)/wgtsum)
          map%col(n_s) = j
          map%row(n_s) = i2
       end do
#else
       !-------------------------------------------------------------
       ! 1st pass computing smoothing matrix weights for column j
       !-------------------------------------------------------------
       where(rdist < 1.e10)
          wgt = exp(-rdist/efold)*areaa
       end where
       wgtsum = sum(wgt)

       !-------------------------------------------------------------
       ! final pass computing smoothing matrix weights for column j
       !-------------------------------------------------------------
       do i=1,map%n_b ! loop over all output domain cells
          ii = iind(i)
          jj = jind(i)
          if (wgt(ii,jj) > 0.0) then
             n_s = n_s + 1
!            write(6,F1) 'n_s = ',n_s
             map%s(n_s) = (map%area_a(j)/map%area_b(i))*(wgt(ii,jj)/wgtsum)
             map%col(n_s) = j
             map%row(n_s) = i
          endif
       end do
#endif
!
!     else
!        !---------------------------------------------------
!        ! if mask_a == 0 (is not considered for smoothing)
!        ! but mask_b != 0 (is an active ocn point)
!        ! then fill in 1.0 on diagonal
!        !      just in case runoff appears where mask_a == 0
!        !---------------------------------------------------
!        ! BK: comment out this else-option...
!        ! don't put 1's on the diagonal -- rather assume that
!        ! if mask_a == 0 (is not considered for smoothing)
!        ! then no runoff will appear in this cell
!        !---------------------------------------------------
!        if(map%mask_b(j) /= 0) then
!           n_s = n_s + 1
!           ii = iind(j)
!           jj = jind(j)
!           map%s(n_s) = 1.0
!           map%col(n_s) = j
!           map%row(n_s) = j
!        endif
      endif                !  if mask_a /= 0
   enddo                !  loop over j=1,map%n_a
   call shr_timer_stop (t00)


   write(6,F3) "done.       "

   !----------------------------------------------------------------------------
   ! document the basic properties of the smoothing matrix
   !----------------------------------------------------------------------------
   allocate(nDest(map%n_a))
   nDest(:) = 0
   do n=1,n_s  ! loop over all non-zero matrix elements
      j = map%col(n)          ! corresponding src cell
      nDest(j) = nDest(j) + 1 ! count # dest cells for each src cell
   end do

   n = 0
   minDest = n_s
   maxDest = 0
   avgDest = 0
   nDest10    = 0
   nDest100   = 0
   nDest1000  = 0
   nDest10000 = 0
   do j=1,map%n_a  ! loop over all input grid cells
      if (nDest(j) > 0) then
         n = n + 1
         minDest = min(minDest,nDest(j))
         maxDest = max(maxDest,nDest(j))
         avgDest = avgDest +   nDest(j)
         if      (nDest(j) < 10 ) then
            nDest10    = nDest10    + 1
         else if (nDest(j) < 100 ) then
            nDest100   = nDest100   + 1
         else if (nDest(j) < 1000) then
            nDest1000  = nDest1000  + 1
         else
            nDest10000 = nDest10000 + 1
         end if
      end if
   end do
   deallocate(nDest)
   avgDest = nint(avgDest/float(n))
 ! avggd = int(avggd/count(map%mask_a == 1))

   write(6,F1) "map n_s   = ",n_s
   write(6,F3) "min wgt     ",minval(wgt)
   write(6,F3) "max wgt     ",maxval(wgt)
   write(6,F3) "wgtsum      ",wgtsum
   write(6,F3) "min S       ",minval(map%s  (1:n_s))
   write(6,F3) "max S       ",maxval(map%s  (1:n_s))
   write(6,F1) "min row     ",minval(map%row(1:n_s))
   write(6,F1) "max row     ",maxval(map%row(1:n_s))
   write(6,F1) "min col     ",minval(map%col(1:n_s))
   write(6,F1) "max col     ",maxval(map%col(1:n_s))
   write(6,*) subName,"min # dst cells per src cell = ",minDest
   write(6,*) subName,"max # dst cells per src cell = ",maxDest
   write(6,*) subName,"avg # dst cells per src cell = ",avgDest
   write(6,*) subName,"# src cells smoothed to   0-  9 cells = ",nDest10
   write(6,*) subName,"# src cells smoothed to  10- 99 cells = ",nDest100
   write(6,*) subName,"# src cells smoothed to 100-999 cells = ",nDest1000
   write(6,*) subName,"# src cells smoothed to   > 999 cells = ",nDest10000

!  write(6,F9) "min # dest pts : ",
!  write(6,F9) "min # dest pts : ",mingd,mingdi,mingdj
!  write(6,F9) "max # dest pts : ",maxgd,maxgdi,maxgdj
!  write(6,F6) "avg # dest pts : ",avggd
!  write(6,F6) "src # having only 1 dest  : ",onegd
!  write(6,F9) "min # boxes    : ",minbx,minbxi,minbxj
!  write(6,F9) "max # boxes    : ",maxbx,maxbxi,maxbxj
!  write(6,F1) 'now, cnt = ',cnt

   !----------------------------------------------------------------------------
   ! dealloc work arrays, resize the mapping, now that its size is known
   !----------------------------------------------------------------------------
   deallocate(wgt)
   deallocate(rdist)
   deallocate(iind)
   deallocate(jind)
   deallocate(i2ind,j2ind)

   allocate(s  (n_s))
   allocate(row(n_s))
   allocate(col(n_s))
   s   = map%s  (1:n_s)
   row = map%row(1:n_s)
   col = map%col(1:n_s)

   deallocate(map%s)
   deallocate(map%row)
   deallocate(map%col)
   map%n_s = n_s
   allocate(map%s  (n_s))
   allocate(map%row(n_s))
   allocate(map%col(n_s))

   map%s   = s
   map%row = row
   map%col = col
   deallocate(s)
   deallocate(row)
   deallocate(col)

   write(6,F00) 'exit subroutine.'

END SUBROUTINE smooth

!===============================================================================

integer FUNCTION iadd(i,di,ni)

   implicit none

   !--- arguments ---
   integer :: i, di, ni

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (di < 0) then
     iadd = mod(i+di-1+ni,ni)+1
   else
!    why was this coded like this??
!    iadd = mod(i,ni)+di
     iadd = mod(i+di,ni)
   endif

END FUNCTION iadd

!===============================================================================

integer FUNCTION jadd(j,dj,nj)

   implicit none

   !--- arguments ---
   integer :: j, dj, nj

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (dj < 0) then
     jadd = max(j+dj,1)
   else
     jadd = min(j+dj,nj)
   endif

END FUNCTION jadd

!===============================================================================

real(r8) FUNCTION distance(x0,y0,x1,y1)

   implicit none

   !--- arguments ---
   real(r8), intent(in) :: x0,y0,x1,y1

   !--- local ---
   real(r8) :: dph,dth,myx0,myx1

!-------------------------------------------------------------------------------
! PURPOSE:
!  o Return the distance in meters between two points on the
!    Earth with coordinates (x0,y0) and (x1,y1) in degrees.
!  o NOTE: does not take into account curvature of the grid...
!-------------------------------------------------------------------------------

   myx0 = x0
   myx1 = x1
   if (x1-x0>180.0) myx0 = x0 + 360.0
   if (x0-x1>180.0) myx1 = x1 + 360.0
   dph  = (myx1-myx0)*cos(DtoR*(y0+y1)/2.0)
   dth  = (y1-y0)
   distance = sqrt(dth**2 + dph**2)*DtoR*rEarth

END FUNCTION distance

!===============================================================================

real(r8) FUNCTION falloff(x1,y1,x2,y2,s)

   implicit none

   !--- arguments ---
   real(r8), intent(in) :: x1,y1,x2,y2,s
   real(r8) :: dr

!-------------------------------------------------------------------------------
! PURPOSE:
!  o Compute exponential fall-off factor given two lon/lat
!    oordinates (x1,y1) and (x2,y2) in degrees.
!  o NOTE: does not take into account curvature of the grid...
!  o Function is g(d) = exp(-d/s)
!                       for distance d, efold scale length s (m)
!-------------------------------------------------------------------------------

    dr = distance(x1,y1,x2,y2)
    falloff = exp(-dr/s)

end FUNCTION falloff

!===============================================================================

SUBROUTINE test_smooth(map,field,value,source)

   implicit none

   !--- arguments ---
   type(sMatrix), intent(in) :: map     ! smoothing map
   real(r8), intent(inout) :: field(:,:)    ! output flux field
   real(r8), intent(in) :: value            ! source flux value
   integer, intent(in) :: source(:,:)   ! source points

   !--- local ---
   integer :: i,n,isrc,jsrc,idst,jdst
   integer :: nsrc
   integer, allocatable :: mask(:,:)

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   nsrc = size(source,2)
   field = 0.0
   allocate(mask(map%ni_b,map%nj_b))
   mask = reshape(map%mask_b,(/map%ni_b,map%nj_b/))
   do n=1,map%n_s
       isrc = mod(map%col(n)-1,map%ni_a)+1
       jsrc = int((map%col(n)-1)/map%ni_a)+1
       idst = mod(map%row(n)-1,map%ni_b)+1
       jdst = int((map%row(n)-1)/map%ni_b)+1
       do i=1,nsrc
         if (isrc == source(1,i) .and. jsrc == source(2,i)) then
            field(idst,jdst) = field(idst,jdst)+map%s(n)*value
         endif
       enddo
   end do

   where(mask == 0) field = 1.e30

END SUBROUTINE test_smooth

!===============================================================================

SUBROUTINE smooth_field_write(map,field,filename)

   implicit none

   !--- arguments ---
   type(sMatrix), intent(in) :: map       ! smoothing map
   real(r8), intent(in) :: field(:,:)         ! output flux field
   character(*) , intent(in) :: filename  ! name of data file

   !--- local ---
   integer :: i,m,n,isrc,jsrc,idst,jdst
   integer :: nsrc
   real(r8)    :: spv                         ! field missing value
   character(len= 8)     :: cdate   ! wall clock date
   character(len=10)     :: ctime   ! wall clock time
   character(len=240)     :: str     ! variable length char string
   character(len=240)     :: attstr  ! netCDF attribute name string
   integer                :: rcode   ! netCDF routine return code
   integer                :: fid     ! netCDF file      ID
   integer                :: vid     ! netCDF variable  ID
   integer                :: did     ! netCDF dimension ID
   integer                :: vdid(2) ! netCDF dimension ID

   !--- formats ---
   character(len=*),parameter :: F00 = "('(map_write) ',3a)"

!-------------------------------------------------------------------------------
! PURPOSE:
! o writes a field computed from smoothing "map" for viewing
! o assumes field grid is equivalent to b grid of "map"
!-------------------------------------------------------------------------------

   spv = maxval(field)

   !-----------------------------------------------------------------
   ! create a new nc file
   !-----------------------------------------------------------------
!  rcode = nf_create(trim(filename),NF_CLOBBER,fid)
   rcode = nf_create(trim(filename),IOR(NF_CLOBBER,NF_64BIT_OFFSET),fid)
   if (rcode.ne.NF_NOERR) write(*,F00) nf_strerror(rcode)

   !-----------------------------------------------------------------
   ! global attributes
   !-----------------------------------------------------------------
    str  = map%title
   rcode = nf_put_att_text(fid,NF_GLOBAL,'title'      ,len_trim(str),str)
    str  = map%normal
   rcode = nf_put_att_text(fid,NF_GLOBAL,'normalization',len_trim(str),str)
    str  = map%method
   rcode = nf_put_att_text(fid,NF_GLOBAL,'map_method' ,len_trim(str),str)
    str  = map%history
   call date_and_time(cdate,ctime) ! f90 intrinsic
    str = 'File created: '//cdate(1:4)//'-'//cdate(5:6)//'-'//cdate(7:8) &
   &                 //' '//ctime(1:2)//':'//ctime(3:4)//':'//ctime(5:6)
   rcode = nf_put_att_text(fid,NF_GLOBAL,'history'    ,len_trim(str),str)
    str  = map%convention
   rcode = nf_put_att_text(fid,NF_GLOBAL,'conventions',len_trim(str),str)
    str  = map%domain_a
   rcode = nf_put_att_text(fid,NF_GLOBAL,'domain_a'   ,len_trim(str),str)
    str  = map%domain_b
   rcode = nf_put_att_text(fid,NF_GLOBAL,'domain_b'   ,len_trim(str),str)

   !-----------------------------------------------------------------
   ! dimension data
   !-----------------------------------------------------------------
   rcode = nf_def_dim(fid, 'ni_b', map%ni_b, did) ! # of points wrt i
   rcode = nf_def_dim(fid, 'nj_b', map%nj_b, did) ! # of points wrt j

   !-----------------------------------------------------------------
   ! define data
   !-----------------------------------------------------------------

   rcode = nf_inq_dimid(fid,'ni_b',vdid(1))
   rcode = nf_inq_dimid(fid,'nj_b',vdid(2))
   rcode = nf_def_var  (fid,'field',NF_DOUBLE,2,vdid,vid)
   str   = 'output of call to test_smooth'
   rcode = nf_put_att_text(fid,vid,"description",len_trim(str),str)
   rcode = nf_put_att_double(fid,vid,"missing_value",NF_FLOAT,1,spv)

   !-----------------------------------------------------------------
   ! put data
   !-----------------------------------------------------------------
   rcode = nf_enddef(fid)

   rcode = nf_inq_varid     (fid,  'field',vid)
   rcode = nf_put_var_double(fid, vid, field)

   rcode = nf_close(fid)

   if (rcode.ne.NF_NOERR) write(*,F00) nf_strerror(rcode)

END SUBROUTINE smooth_field_write

!===============================================================================

SUBROUTINE dummyflux(imt,jmt,mask,xt,yt,at,runfile,efold)

!-------------------------------------------------------------------------------
!    Based on gconst of /home/dataproc/yeager/RUNOFF/fortran/garray.f
!    which was invoked as a shared object by the NCL routine
!    /home/dataproc/yeager/RUNOFF/construct_runoff_netcdf_dummy.ncl
!
! PURPOSE:
!  o Construct runoff mapping based on a 'runfile' data file which
!    partitions known river runoffs into a 19-basin scheme by Large.
!  o runfile specifies rivers/coastal discharges (kg/s) for each basin along
!    with lon/lat boxes within which this discharge is to be distributed.
!    The distribution array is g(:,:,:) with units (kg/s/m^2).
!    This routine distributes 'coastal' runoff as:
!
!        g(i,j,thisbasin) = [total coastal discharge (kg/s)]/[total
!                               distribution area (m^2)]
!
!    This routine distributes 'river' runoff as:
!
!        g(i,j,thisriver) = [exponential falloff factor]*[total river
!                               discharge (kg/s)]/[total exponential-weighted
!                               distribution area (m^2)]
!
!    Areas (from array at) are based on ocean grid TAREA.
!
!-------------------------------------------------------------------------------

   !--- local ---
   integer :: imt,jmt,jm,jp,im,ip,nr,n,i,j
   real(r8)    :: rlat,rlon,xf
   integer, dimension(imt,jmt) :: mask  ! mask of grid (0 for land)
   real(r8), dimension(imt,jmt) :: xt,yt    ! lon,lat of grid
   real(r8), dimension(imt,jmt) :: at       ! area of grid (m^2)
   real(r8), dimension(imt,jmt) :: A        ! masked out area array
   real(r8), dimension(imt,jmt) :: fi       ! exponential factor
   real(r8), allocatable :: total(:)        ! total runoff by basin (kg/s)
   character*47 :: runfile
   character*27 :: rname
   integer, parameter :: nrmax=70       ! max number of rivers
   integer :: nb,nbo,ibasin,nrivers,msea
   real(r8) :: efold                        ! efold scale (m) for falloff
   real(r8) :: rarea                        ! sum of destination cell areas
   real(r8) :: firarea                      ! exponential-weighted sum of
                                        ! destination cell areas
   real(r8) :: sumb,sumd,gmax
   real(r8), dimension(2,nrmax) :: blon,blat
   real(r8), dimension(nrmax) :: clon,clat,nref,discharge
   character*40, allocatable :: bname(:)

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   ! Check that received grid lon/lat arrays are reasonable
   if (minval(xt).lt.0.0 .or. maxval(xt).gt.360.0) then
      write(6,*) 'xt out of range'
      stop
   endif
   if (minval(yt).lt.-90.0 .or. maxval(yt).gt.90.0) then
      write(6,*) 'yt out of range'
      stop
   endif

   open(13,file=runfile,status='old',form='formatted')
   read(13,*)
   read(13,*) nbo
   nb = nbo+1           ! nbo seperate basins, plus global
   allocate(garr(imt,jmt,nb))
   allocate(total(nb))
   allocate(bname(nb))
   total = 0.0
   garr = 0.0

   bname(1) = '    Global'

   ! Process each of the basins
   do n=2,nb
     read(13,99) ibasin,nrivers,msea,clon(1),clat(1),bname(n)
     sumd = 0.0
     do nr=1,nrivers
       read (13,100) nref(nr),rname,discharge(nr),clon(nr),clat(nr), &
     &        blon(1,nr),blon(2,nr),blat(1,nr),blat(2,nr)
       sumd = sumd + discharge(nr)
     enddo
     total(n) = sumd    ! includes coastal runoff

     do nr=1,nrivers
       xf = discharge(nr) / sumd
       rarea = 0.0
       firarea = 0.0
       do j=1,jmt
         do i=1,imt
!             m=(j-1) * (imt-2) + i - 1
           rlon = xt(i,j)
           rlat = yt(i,j)
           if (blon(1,nr).lt.0.0.and.xt(i,j).gt.180.) rlon = xt(i,j)-360.0

           ! If the grid point lies in the specified lon/lat box
           ! where the discarge is to be distributed...
           if ((blat(1,nr).le.rlat).and.(rlat.le.blat(2,nr)).and. &
     &         (mask(i,j).ne.0).and.(blon(1,nr).le.rlon).and. &
     &         (rlon.le.blon(2,nr))) then

              if (nref(nr).eq.0) then   ! for Coastal runoff
                                        ! check that destination cell
                                        ! is adjacent to land
                jm = jadd(j,-1,jmt)
                jp = jadd(j,1,jmt)
                im = iadd(i,-1,imt)
                ip = iadd(i,1,imt)
!       exception for (OLD) grids with wraparound points
        if ((imt.eq.102.or.imt.eq.152).and.i.eq.1) im = imt-1
        if ((imt.eq.102.or.imt.eq.152).and.i.eq.imt) ip = 2

                if ((mask(im,j).eq.0).or.(mask(im,jp).eq.0) &
     & .or.(mask(i,jp).eq.0).or.(mask(ip,jp).eq.0).or. &
     & (mask(ip,j).eq.0).or.(mask(ip,jm).eq.0).or. &
     & (mask(i,jm).eq.0).or.(mask(im,jm).eq.0)) then
                  A(i,j) = at(i,j)
                  rarea = rarea + A(i,j)
                else
                  A(i,j) = 0.0
                endif
              else                      ! for River runoff
                 A(i,j)  = at(i,j)
                 rarea = rarea + A(i,j)
                 fi(i,j) = falloff(rlon,rlat,clon(nr),clat(nr),efold)
                 firarea = firarea + fi(i,j)*A(i,j)
              endif
           else
              A(i,j)  = 0.0
              fi(i,j)  = 0.0
           endif
         enddo
       enddo

       do j=1,jmt
         do i=1,imt
           !  for unitless
           !  garr(i,j,n) = garr(i,j,n) + xf * A(i,j) / rarea
           !  for [kg/m^2/s]
           if (A(i,j) .gt. 0.) then
             if (nref(nr).eq.0) then
               garr(i,j,n) = garr(i,j,n) + discharge(nr)/rarea
             else
               garr(i,j,n) = garr(i,j,n) + fi(i,j)*discharge(nr) / firarea
             endif
           endif
         enddo
       enddo

     enddo     ! do nr=1,nrivers
   enddo       ! do n=2,nb

! Compute global total runoff
   do n=2,nb
     total(1) = total(1) + total(n)
   enddo
   do j=1,jmt
     do i=1,imt
       do n=2,nb
          garr(i,j,1) = garr(i,j,1) + garr(i,j,n)
       enddo
     enddo
   enddo

! Test runoff conservation and write summary table to runoff.doc
   write(6,*) 'Results of Data Runoff Distribution:'
   write(6,101) 'Basin','Global g*area (kg/s)','Data (kg/s)','Difference'
   do n=1,nb
     gmax = 0.0
     sumb = 0.0
     do j=1,jmt
       do i=1,imt
         sumb = sumb + garr(i,j,n)*at(i,j)
         if (garr(i,j,n).gt.gmax) gmax = garr(i,j,n)
       enddo
     enddo
     write(6,102) n,bname(n),sumb,total(n),(sumb - total(n))
   enddo

99     format(3i4,2f8.1,a30)
100    format(i4,2x,a27,f12.1,6f7.1)
101    format(a10,11x,a25,4x,a11,4x,a17)
102    format(i2,2x,a26,2x,f14.3,3x,f14.3,3x,f14.3)

END SUBROUTINE dummyflux

!===============================================================================

SUBROUTINE dummyweights(map,at)

   implicit none

   !--- arguments ---
   type(sMatrix), intent(inout) :: map  ! scrip map template to be filled
   real(r8)         :: at(:,:)              ! ocn TAREA in m^2

   !--- local ---
   integer :: ni,nj,imt,jmt,i,j,m,n_s,nbasin
   real(r8)    :: sumb
   real(r8), allocatable  :: col(:,:),row(:,:),s(:,:),garr2d(:,:)
   real(r8), allocatable  :: col1d(:),row1d(:),s1d(:),cplratio(:),at1d(:)
   logical, allocatable  :: smask(:,:)
   real(r8), allocatable  :: total(:)

!-------------------------------------------------------------------------------
! PURPOSE:
! o Based on M. Hecht's /fs/cgd/home0/hecht/csm/runoff/transfer_matrix.ncl
! o Compute correct weights to map 19-element 'normalized' runoff flux vector
!   (kg/s/m^2) onto the POP ocean grid.
!
!       R(j)  = net runoff (kg/s) into basin j
!       dA(j) = fictitious area (m^2) associated with basin j (1/19 of globe)
!               from scrip
!       F(j)  = "flux" sent to cpl by dlnd = R(j)/dA(j) (kg/s/m^2)
!               (This flux is stored in the data_aroff netcdf file)
!       garr(i)  = flux which should go into ocean cell i, as determined by
!                routine dummyflux (sum over i is known to preserve R when
!                multiplied by ocean areas!)
!       dA_o(i) = area (m^2) associated with ocean cell i, from POP
!       dA_s(i) = area (m^2) associated with ocean cell i, from scrip
!
!   Then in order to preserve R(i), the mapping weights should be
!       s(i,j) = dA(j)*[garr(i)/R(j)]*[dA_o(i)/dA_s(i)]
!
!   The "flux" mapping, then will be:
!       Fo(i) = sum_over_j{s(i,j)*F(j)}
!             = sum_over_j{dA(j)   garr(i)    dA_o(i)   R(j)  }
!                         {     * --------- * ------- * ----  }
!                         {        R(j)       dA_s(i)   dA(j) }
!             = garr(i) * dA_o(i)/dA_s(i)
!
!   Before sending to the ocean, the couple will multiply Fo(i) by
!   an area correction factor, so that
!       F*o(i) = runoff flux received by ocean model (kg/s/m^2)
!              = Fo(i) * dA_s(i)/dA_o(i)
!              = garr(i), by definition the correct flux into ocean cell i
!
!   Therefore, sum_over_i{F*o(i)*dAo(i)} = sum_over_j{R(j)}
!
!-------------------------------------------------------------------------------

   imt = map%ni_b
   jmt = map%nj_b
   ni = map%ni_a
   nj = map%nj_a

   allocate(col(map%n_b,map%n_a))
   allocate(row(map%n_b,map%n_a))
   allocate(s(map%n_b,map%n_a))
   allocate(total(map%n_a))
   allocate(smask(map%n_b,map%n_a))
   nbasin = size(garr,3)
   allocate(garr2d(map%n_b,nbasin-1))
   allocate(cplratio(map%n_b))
   allocate(at1d(map%n_b))

   ! This is the areafact ratio that cpl will multiply
   ! runoff flux with before sending to ocean
   cplratio = (map%area_b*(rEarth**2))/reshape(at,(/map%n_b/))

   garr2d = reshape(garr(:,:,2:nbasin),(/map%n_b,nbasin-1/))
   at1d = reshape(at,(/map%n_b/))

   !-----------------------------------------------------------------
   ! compute the 19-basin total runoffs (kg/s)
   !-----------------------------------------------------------------
   do j=1,map%n_a
     sumb = 0.0
     do i=1,map%n_b
       sumb = sumb + garr2d(i,j)*at1d(i)
     enddo
     total(j) = sumb
   enddo

   ! Using garr, fill in the mapping matrix
   do i=1,map%n_b
     do j=1,map%n_a
       col(i,j) = j
       row(i,j) = i
       s(i,j) = map%area_a(j)*(rEarth**2)*garr2d(i,j)/(total(j)*cplratio(i))
    !  s(i,j) = map%area_a(j)*(rEarth**2)*garr2d(i,j)/(total(j))
     enddo
   enddo

   !-----------------------------------------------
   ! resize the mapping, now that its size is known
   !-----------------------------------------------
   smask = .false.
   where(s.gt.0.0) smask=.true.
   n_s = count(smask)
   allocate(s1d(n_s))
   allocate(row1d(n_s))
   allocate(col1d(n_s))
   s1d = pack(s,smask)
   col1d = pack(col,smask)
   row1d = pack(row,smask)

   deallocate(map%s)
   deallocate(map%row)
   deallocate(map%col)
   map%n_s = n_s
   allocate(map%s(n_s))
   allocate(map%row(n_s))
   allocate(map%col(n_s))
   map%s = s1d
   map%row = row1d
   map%col = col1d

END SUBROUTINE dummyweights

!===============================================================================

SUBROUTINE dummy_aroff_write(map,srcfile,datafile,at)

   implicit none

   !--- arguments ---
   type(sMatrix), intent(in) :: map       ! scrip map template
   character(*) , intent(in) :: srcfile   ! name of ascii river runoff file
   character(*) , intent(in) :: datafile  ! name of runoff data file
   real(r8) , intent(in) :: at(:,:)           ! ocean TAREA (m^2)

   !--- local ---
   integer :: ni,nj,imt,jmt,i,j,m,n
   real(r8)    :: sumb
   real(r8), allocatable     :: work(:,:),work2(:,:,:),basinroff(:,:)
   integer, allocatable  :: iwork(:,:)
   character(len= 8)     :: cdate   ! wall clock date
   character(len=10)     :: ctime   ! wall clock time
   character(len=240)     :: str     ! variable length char string
   character(len=240)     :: attstr  ! netCDF attribute name string
   integer                :: rcode   ! netCDF routine return code
   integer                :: fid     ! netCDF file      ID
   integer                :: vid     ! netCDF variable  ID
   integer                :: did     ! netCDF dimension ID
   integer                :: vdid(3) ! netCDF dimension ID

   !--- formats ---
   character(len=*),parameter :: F00 = "('(dummy_aroff_write) ',3a)"

!-------------------------------------------------------------------------------
! PURPOSE:
! o Once garr has been defined after a call to dummyflux, this routine
!   will output a data.runoff.nc file for 19-basin data runoff.  This
!   file is needed by dlnd6, and both "domain.runoff.nc" and "data.runoff.nc"
!   get softlinked to this file in the dlnd buildnml script.  The dlnd
!   namelist parameter data_aroff will get set to this file.
!
! o The contents of the netcdf are fictitious 19-basin grid information (derived
!   from the r19.nc file through the sMatrix "map") and runoff values for each
!   basin in kg/s/m^2.  These values are computed as the total runoff for each
!   basin as specified in srcfile (kg/s) normalized by a fictitious area
!   (1/19 of the globe) in m^2.
!-------------------------------------------------------------------------------

   imt = map%ni_b
   jmt = map%nj_b
   ni = map%ni_a
   nj = map%nj_a
   allocate(iwork(ni,nj))
   allocate(work(ni,nj))
   allocate(work2(ni,nj,4))
   allocate(basinroff(ni,nj))

   !-----------------------------------------------------------------
   ! compute the normalized 19-basin totals (kg/s/m^2)
   ! (This is total of each of 19 basins (kg/s) divided by bogus
   !  'basin area' = (area of Earth)/19 ).
   !-----------------------------------------------------------------
   do n=1,ni
     sumb = 0.0
     do j=1,jmt
       do i=1,imt
         sumb = sumb + garr(i,j,n+1)*at(i,j)
       enddo
     enddo
     basinroff(n,1) = sumb/(map%area_a(n)*rEarth*rEarth)
   enddo

   !-----------------------------------------------------------------
   ! create the 19-basin data_aroff netcdf used by dlnd6
   !-----------------------------------------------------------------
!  rcode = nf_create(trim(datafile),NF_CLOBBER,fid)
   rcode = nf_create(trim(datafile),IOR(NF_CLOBBER,NF_64BIT_OFFSET),fid)
   if (rcode.ne.NF_NOERR) write(*,F00) nf_strerror(rcode)

   !-----------------------------------------------------------------
   ! global attributes
   !-----------------------------------------------------------------
    str  = 'runoff data, fictitious 19 basin domain'
   rcode = nf_put_att_text(fid,NF_GLOBAL,'title'      ,len_trim(str),str)
   call date_and_time(cdate,ctime) ! f90 intrinsic
    str = 'File created: '//cdate(1:4)//'-'//cdate(5:6)//'-'//cdate(7:8) &
   &                 //' '//ctime(1:2)//':'//ctime(3:4)//':'//ctime(5:6)
   rcode = nf_put_att_text(fid,NF_GLOBAL,'history'    ,len_trim(str),str)
    str  = 'created using dummy runoff tools in /cgd/oce/yeager/POP_tools/'// &
   &    'runoff/runoff_smooth_mod.F90'
   rcode = nf_put_att_text(fid,NF_GLOBAL,'creation',len_trim(str),str)
    str  = 'source data file was '//trim(srcfile)
   rcode = nf_put_att_text(fid,NF_GLOBAL,'data'   ,len_trim(str),str)

   !-----------------------------------------------------------------
   ! dimensions
   !-----------------------------------------------------------------
   rcode = nf_def_dim(fid, 'ni', ni, did) ! # of points wrt i
   rcode = nf_def_dim(fid, 'nj', nj, did)  ! # of points wrt j
   rcode = nf_def_dim(fid, 'nv', 4, did)  ! # of points wrt j

   !-----------------------------------------------------------------
   ! variables
   !-----------------------------------------------------------------
   rcode = nf_inq_dimid(fid,'nv',vdid(1))
   rcode = nf_inq_dimid(fid,'ni',vdid(2))
   rcode = nf_inq_dimid(fid,'nj',vdid(3))
   rcode = nf_def_var(fid,'xc',NF_DOUBLE,2,vdid(2:3),vid)
   str   = 'degrees'
   rcode = nf_put_att_text(fid,vid,'units',len_trim(str),str)
   str   = 'longitude of runoff grid cell center'
   rcode = nf_put_att_text(fid,vid,'long_name',len_trim(str),str)
   rcode = nf_def_var(fid,'yc',NF_DOUBLE,2,vdid(2:3),vid)
   str   = 'degrees'
   rcode = nf_put_att_text(fid,vid,'units',len_trim(str),str)
   str   = 'latitude of runoff grid cell center'
   rcode = nf_put_att_text(fid,vid,'long_name',len_trim(str),str)
   rcode = nf_def_var(fid,'xv',NF_DOUBLE,3,vdid,vid)
   str   = 'degrees'
   rcode = nf_put_att_text(fid,vid,'units',len_trim(str),str)
   str   = 'longitudes of runoff grid cell vertices'
   rcode = nf_put_att_text(fid,vid,'long_name',len_trim(str),str)
   rcode = nf_def_var(fid,'yv',NF_DOUBLE,3,vdid,vid)
   str   = 'degrees'
   rcode = nf_put_att_text(fid,vid,'units',len_trim(str),str)
   str   = 'latitudes of runoff grid cell vertices'
   rcode = nf_put_att_text(fid,vid,'long_name',len_trim(str),str)
   rcode = nf_def_var(fid,'mask',NF_SHORT,2,vdid(2:3),vid)
   str   = 'unitless'
   rcode = nf_put_att_text(fid,vid,'units',len_trim(str),str)
   str   = 'runoff grid domain mask'
   rcode = nf_put_att_text(fid,vid,'long_name',len_trim(str),str)
   rcode = nf_def_var(fid,'area',NF_DOUBLE,2,vdid(2:3),vid)
   str   = 'square radians'
   rcode = nf_put_att_text(fid,vid,'units',len_trim(str),str)
   str   = 'area of runoff grid cell'
   rcode = nf_put_att_text(fid,vid,'long_name',len_trim(str),str)
   rcode = nf_def_var(fid,'runoff',NF_DOUBLE,2,vdid(2:3),vid)
   str   = 'kg/s/m^2'
   rcode = nf_put_att_text(fid,vid,'units',len_trim(str),str)
   str   = 'Basin total runoff normalized by 1/19 of Globe'
   rcode = nf_put_att_text(fid,vid,'long_name',len_trim(str),str)
   rcode = nf_put_att_double(fid,vid,"_FillValue",NF_DOUBLE,1,-9999.)
   rcode = nf_def_var(fid,'rEarth',NF_DOUBLE,0,vdid,vid)
   str   = 'meters'
   rcode = nf_put_att_text(fid,vid,'units',len_trim(str),str)
   str   = 'Radius of the Earth'
   rcode = nf_put_att_text(fid,vid,'long_name',len_trim(str),str)

   !-----------------------------------------------------------------
   ! finish writing it
   !-----------------------------------------------------------------
   rcode = nf_enddef(fid)
   rcode = nf_inq_varid(fid,'xc',vid)
   work = reshape(map%xc_a,(/ni,nj/))
   rcode = nf_put_var_double(fid, vid, work)
   rcode = nf_inq_varid(fid,'yc',vid)
   work = reshape(map%yc_a,(/ni,nj/))
   rcode = nf_put_var_double(fid, vid, work)
   rcode = nf_inq_varid(fid,'xv',vid)
   work2 = reshape(map%xv_a,(/ni,nj,4/))
   rcode = nf_put_var_double(fid, vid, work2)
   rcode = nf_inq_varid(fid,'yv',vid)
   work2 = reshape(map%yv_a,(/ni,nj,4/))
   rcode = nf_put_var_double(fid, vid, work2)
   rcode = nf_inq_varid(fid,'mask',vid)
   iwork = reshape(map%mask_a,(/ni,nj/))
   rcode = nf_put_var_int(fid, vid, iwork)
   rcode = nf_inq_varid(fid,'area',vid)
   work = reshape(map%area_a,(/ni,nj/))
   rcode = nf_put_var_double(fid, vid, work)
   rcode = nf_inq_varid(fid,'rEarth',vid)
   rcode = nf_put_var_double(fid, vid, rEarth)
   rcode = nf_inq_varid(fid,'runoff',vid)
   rcode = nf_put_var_double(fid, vid, basinroff)
   rcode = nf_close(fid)
   if (rcode.ne.NF_NOERR) write(*,F00) nf_strerror(rcode)

END SUBROUTINE dummy_aroff_write

!===============================================================================

SUBROUTINE POP_TAREA_compute(filename,nx,ny,area)

   implicit none

   !--- arguments ---
   character(*) , intent(in)    :: filename  ! name of ocn grid binary
   integer      :: nx,ny
   real(r8) , allocatable, intent(out)   :: area(:,:)

   !--- local ---
   real(r8) , allocatable   :: work(:,:),htn(:,:),hte(:,:),dxt(:,:),dyt(:,:)

   !--- formats ---
   character(len=*),parameter :: F00 = "('(POP_TAREA_compute) ',3a)"

!-------------------------------------------------------------------------------
! PURPOSE:
! o Compute POP TAREA array from grid.ieeer8 binary, as it is done
!   in the POP code.
!-------------------------------------------------------------------------------

   allocate(work(nx,ny))
   allocate(htn(nx,ny))
   allocate(hte(nx,ny))
   allocate(dxt(nx,ny))
   allocate(dyt(nx,ny))
   allocate(area(nx,ny))

   open(10,file=filename,form='unformatted',access='direct', &
  &     recl=nx*ny*8)
   read(10,rec=1) work
   read(10,rec=2) work
   read(10,rec=3) htn
   read(10,rec=4) hte
   close(10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Below code fragment is from POP grid.F
!  slightly modified
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !***
      !*** construct T-grid cell widths
      !***
      work = cshift(htn,-1,2)
      dxt = 0.5*(htn + work)
      dxt(:,1) = dxt(:,2)       ! reasonable kluge

      work = cshift(hte,-1,1)
      dyt = 0.5*(hte + work)

      where (dxt == 0.0) dxt=1.0
      where (dyt == 0.0) dyt=1.0

      area = dxt*dyt

END SUBROUTINE POP_TAREA_compute

!===============================================================================

SUBROUTINE POP_2D_read(filename,varname,data)

   !--- modules ---

   implicit none

   !--- includes ---
#include <netcdf.inc>

   !--- arguments ---
   character(*) , intent(in)    :: filename  ! name of data file
   character(*) , intent(in)    :: varname  ! name of data file
   real(r8) , intent(inout)         :: data(:,:)

   !--- local ---
   integer         :: nxin,nyin,nx,ny,ndims,itype     ! generic indicies
   integer         :: fid,vid,did,rcode
   integer, dimension(2)  :: dimid

   !--- formats ---
   character(len=*),parameter :: F00 = "('(POP_2D_read) ',3a)"

!-------------------------------------------------------------------------------
! PURPOSE:
! o reads map matrix information from netCDF data file
!
! NOTE:
!-------------------------------------------------------------------------------

   nxin = size(data,1)
   nyin = size(data,2)

   !-----------------------------------------------
   ! read in the variable
   !-----------------------------------------------
   rcode = nf_open(filename,NF_NOWRITE,fid)
   rcode = nf_inq_varid(fid,varname,vid)
   if (rcode.ne.NF_NOERR) write(*,F00) nf_strerror(rcode)
   rcode = nf_inq_varndims(fid,vid,ndims)
   if (ndims.ne.2) then
        write(*,F00) 'variable did not have dim 2'
        stop
   endif
   rcode = nf_inq_vardimid(fid,vid,dimid)
   rcode = nf_inq_dimlen(fid, dimid(1)   , nx  )
   rcode = nf_inq_dimlen(fid, dimid(2)   , ny  )
   if (nx.ne.nxin.or.ny.ne.nyin) then
        write(*,F00) 'variable does not have correct size',nx,ny
        stop
   endif
   rcode = nf_inq_vartype(fid,vid,itype)

   select case (itype)
   case (NF_FLOAT)
      rcode = nf_get_var_real(fid,vid,data)
   case default
      write(*,F00) 'can only handle float'
      stop
   end select

   rcode = nf_close(fid)

END SUBROUTINE POP_2D_read

!===============================================================================
recursive SUBROUTINE breadth_setDist(ni,nj,level,x,y,mask, rdist,rmax,i2ind,j2ind,strptr,length)

!-----------------------------------------------------------------------
!  Recursively sets a distance from a seed point (i0,j0)
!-----------------------------------------------------------------------

   integer :: i0,j0
   integer :: ni, nj
   integer :: kk
   integer :: ie,je,iw,jw,in,jn,is,js,k,j
   integer :: ine,jne,inw,jnw,ise,jse,isw,jsw
   integer :: level
   integer, dimension(ni,nj) :: mask
   real(r8), intent(inout),dimension(ni,nj) :: rdist
   real(kind=r8) :: rmax
   real(r8) :: d
   real(r8) :: ds,dn,de,dw
   real(r8) :: dne, dnw, dse, dsw
   real(r8), dimension(ni*nj) :: x,y
   integer, dimension(:) :: i2ind,j2ind
   integer :: strptr,strPtr2
   integer :: length,length2
   logical, parameter :: Debug = .false.

   if(Debug) print *,'breadth search: level is: ',level,'length is: ',length
   if(Debug) print *,'breadth_setDist:', rmax

   strPtr2=length
   length2=length
   do kk=strPtr,length2
      i0 = i2ind(kk)
      j0 = j2ind(kk)
      ! ------------------------------------------
      ! only perform the search if patch length
      ! is below cutoff threshold
      ! ------------------------------------------
      if(rdist(i0,j0) < rmax) then
      if(Debug) print *,'rdist(i0,j0) rmax: ',rdist(i0,j0),rmax

      j = (j0-1)*ni+i0

      !-----------------
      ! West neighbor
      !-----------------
      iw = iadd(i0,-1,ni)
      jw = j0
      if(mask(iw,jw) == -1) then
          k = (j0-1)*ni+iw
          dw = distance(x(k),y(k),x(j),y(j))+rdist(i0,j0)
          mask(iw,jw) = level
          rdist(iw,jw) = dw
          length = length + 1
          i2ind(length) = iw
          j2ind(length) = jw
          if(Debug) print *,'distance to west point: ',dw
      endif

      !-----------------
      ! East neighbor
      !-----------------
      ie = iadd(i0,1,ni)
      je = j0
      if(mask(ie,je) == -1) then
         k = (j0-1)*ni+ie
          de = distance(x(k),y(k),x(j),y(j))+rdist(i0,j0)
          mask(ie,je) = level
          rdist(ie,je) = de
          length = length + 1
          i2ind(length) = ie
          j2ind(length) = je
          if(Debug) print *,'distance to east point: ',de
      endif

      !----------------
      ! North Neighbor
      !----------------
      jn = jadd(j0,1,nj)
      in = i0
      if(mask(in,jn) == -1) then
          k = (jn-1)*ni+i0
          dn = distance(x(k),y(k),x(j),y(j))+rdist(i0,j0)
          mask(in,jn) = level
          rdist(in,jn) = dn
          length = length + 1
          i2ind(length) = in
          j2ind(length) = jn
          if(Debug) print *,'distance to north point: ',dn
      endif

      !----------------
      ! South Neighbor
      !----------------
      js = jadd(j0,-1,nj)
      is = i0
      if(mask(is,js) == -1) then
         k = (js-1)*ni+i0
         ds = distance(x(k),y(k),x(j),y(j))+rdist(i0,j0)
         mask(is,js) = level
         rdist(is,js) = ds
         length = length + 1
         i2ind(length) = is
         j2ind(length) = js
         if(Debug) print *,'distance to south point: ',ds
      endif
    endif

#if 0
! Note: this was meant to add diagonal search to smoother, but
!       it doesn't work so it is not run at this time
      !-----------------
      ! North West neighbor
      !-----------------
      inw = iadd(i0,-1,ni)
      jnw = jadd(j0,1,nj)
      if(mask(inw,jnw) == -1) then
          k = (jnw-1)*ni+inw
          dnw = distance(x(k),y(k),x(j),y(j))+rdist(i0,j0)
          mask(inw,jnw) = level
          rdist(inw,jnw) = dnw
          length = length + 1
          i2ind(length) = inw
          j2ind(length) = jnw
          if(Debug) print *,'distance to northwest point: ',dnw
      endif

      !-----------------
      ! North East neighbor
      !-----------------
      ine = iadd(i0,1,ni)
      jne = jadd(j0,1,nj)
      if(mask(ine,jne) == -1) then
          k = (jne-1)*ni+ine
          dne = distance(x(k),y(k),x(j),y(j))+rdist(i0,j0)
          mask(ine,jne) = level
          rdist(ine,jne) = dne
          length = length + 1
          i2ind(length) = ine
          j2ind(length) = jne
          if(Debug) print *,'distance to northwest point: ',dne
      endif

      !-----------------
      ! South East neighbor
      !-----------------
      ise = iadd(i0,1,ni)
      jse = jadd(j0,-1,nj)
      if(mask(ise,jse) == -1) then
          k = (jse-1)*ni+ise
          dse = distance(x(k),y(k),x(j),y(j))+rdist(i0,j0)
          mask(ise,jse) = level
          rdist(ise,jse) = dse
          length = length + 1
          i2ind(length) = ise
          j2ind(length) = jse
          if(Debug) print *,'distance to northwest point: ',dse
      endif

      !-----------------
      ! South West neighbor
      !-----------------
      isw = iadd(i0,-1,ni)
      jsw = jadd(j0,-1,nj)
      if(mask(isw,jsw) == -1) then
          k = (jsw-1)*ni+isw
          dsw = distance(x(k),y(k),x(j),y(j))+rdist(i0,j0)
          mask(isw,jsw) = level
          rdist(isw,jsw) = dsw
          length = length + 1
          i2ind(length) = isw
          j2ind(length) = jsw
          if(Debug) print *,'distance to northwest point: ',dsw
      endif
#endif

  enddo
  level = level + 1
  if(length > strptr2) then
      call breadth_setDist(ni,nj,level,x,y,mask, rdist,rmax,i2ind,j2ind,strptr2,length)
  endif

END SUBROUTINE breadth_setDist

recursive SUBROUTINE depth_setDist(i0,j0,ni,nj,x,y,value,mask,rdist,rmax,i2ind,j2ind,length)

!-----------------------------------------------------------------------
!  Recursively sets a distance from a seed point (i0,j0)
!-----------------------------------------------------------------------

   integer :: i0, j0, ni, nj
   integer :: ie,iw,jn,js,k,j
   integer, dimension(ni,nj) :: mask
   real(r8), dimension(ni,nj) :: rdist
   real(r8) :: rmax,d,value
   real(r8), dimension(ni*nj) :: x,y
   integer, dimension(:) :: i2ind,j2ind
   integer :: length

   if (mask(i0,j0) /= 0 .and. value < rmax .and. value < rdist(i0,j0)) then
      j = (j0-1)*ni+i0
      if(rdist(i0,j0) > 1.e10) then
         ! first visit to point
         length = length + 1
         i2ind(length) = i0
         j2ind(length) = j0
      endif
      rdist(i0,j0)=value
      mask(i0,j0)=1

      iw = iadd(i0,-1,ni)
      k = (j0-1)*ni+iw
      d = distance(x(k),y(k),x(j),y(j))+rdist(i0,j0)
      call depth_setDist(iw,j0,ni,nj,x,y,d,mask,rdist,rmax,i2ind,j2ind,length)

      ie = iadd(i0,1,ni)
      k = (j0-1)*ni+ie
      d = distance(x(k),y(k),x(j),y(j))+rdist(i0,j0)
      call depth_setDist(ie,j0,ni,nj,x,y,d,mask,rdist,rmax,i2ind,j2ind,length)

      jn = jadd(j0,1,nj)
      k = (jn-1)*ni+i0
      d = distance(x(k),y(k),x(j),y(j))+rdist(i0,j0)
      call depth_setDist(i0,jn,ni,nj,x,y,d,mask,rdist,rmax,i2ind,j2ind,length)

      js = jadd(j0,-1,nj)
      k = (js-1)*ni+i0
      d = distance(x(k),y(k),x(j),y(j))+rdist(i0,j0)
      call depth_setDist(i0,js,ni,nj,x,y,d,mask,rdist,rmax,i2ind,j2ind,length)

   endif

END SUBROUTINE depth_setDist

!===============================================================================
END MODULE smooth_mod
!===============================================================================
