subroutine initcom

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialize Model commons, including COMCON, COMHYB, COMMAP, COMSPE,
! and COMTRCNM
! 
! Method: 
! 
! Author: 
! Original version:  CCM1
! Standardized:      L. Bath, Jun 1992
!                    L. Buja, Feb 1996
!
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------
   use shr_kind_mod,    only: r8 => shr_kind_r8
   use pmgrid,          only: plat, plev, plon, plevp
   use scanslt,         only: j1, platd, nlonex
   use grmult,          only: grmult_init
   use pspect
   use comspe
   use rgrid,           only: nlon, nmmax, wnummax, beglatpair, fullgrid
   use gauaw_mod,       only: gauaw
   use commap
   use physconst,       only: rair, rearth, ra
   use time_manager,    only: get_step_size
   use abortutils,      only: endrun
   use spmd_utils,      only: masterproc
   use hycoef,          only: hypd, hypm
   use sld_control_mod, only: ifax, trig
   use cam_logfile,     only: iulog
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
!
! Local workspace
!
   real(r8) zsi(plat)      ! sine of latitudes
   real(r8) zw(plat)       ! Gaussian weights
   real(r8) zra2           ! ra squared
   real(r8) zalp(2*pspt)   ! Legendre function array
   real(r8) zdalp(2*pspt)  ! Derivative array
   real(r8) zslat          ! sin of lat  and cosine of colatitude

   integer i           ! longitude index
   integer j           ! Latitude index
   integer k           ! Level index
   integer kk          ! Level index
   integer kkk         ! Level index
   integer m           ! Index for legendre array
   integer nkk         ! Print control variables
   integer ik1         ! Print index temporary variable
   integer ik2         ! Print index temporary variable
   integer itmp        ! Dimension of polynomial arrays temporary.
   integer iter        ! Iteration index
   real(r8)    zdt         ! Time step for settau

   logical lprint      ! Debug print flag
   integer irow        ! Latitude pair index
   integer lat         ! Latitude index

   real(r8) xlat           ! Latitude (radians)
   real(r8) pi             ! Mathematical pi (3.14...)
   real(r8) dtime          ! timestep size [seconds]
!
!-----------------------------------------------------------------------
!
   lprint = masterproc .and. .FALSE.

   dtime = get_step_size()

   call hdinti  (rearth  ,dtime   )

!
! NMAX dependent arrays
!
   if (pmmax.gt.plon/2) then
      call endrun ('INITCOM:mmax=ptrm+1 .gt. plon/2')
   end if
   zra2 = ra*ra
   do j=2,pnmax
      sq(j)  = j*(j-1)*zra2
      rsq(j) = 1._r8/sq(j)
   end do
   sq(1)  = 0._r8
   rsq(1) = 0._r8
!
! MMAX dependent arrays
!
   do j=1,pmmax
      xm(j) = j-1
   end do
!
! Gaussian latitude dependent arrays
!
   call gauaw(zsi     ,zw      ,plat    )
   do irow=1,plat/2
      slat(irow) = zsi(irow)
      w(irow)              = zw(irow)
      w(plat - irow + 1)   = zw(irow)
      cs(irow)  = 1._r8 - zsi(irow)*zsi(irow)
      xlat = asin(slat(irow))
      clat(irow) = -xlat
      clat(plat - irow + 1) = xlat
   end do

   do lat=1,plat
      latdeg(lat) = clat(lat)*45._r8/atan(1._r8)
   end do
!
! Integration matrices of hydrostatic equation(href) and conversion
! term(ecref).  ecref calculated to be consistent with continuity Eq.;
! href calculated to conserve energy.
!
   do k=1,plev
      do kk=1,plev
         href(kk,k) = 0._r8
         ecref(kk,k) = 0._r8
      end do
   end do
!
! Mean atmosphere energy conversion term is consistent with continiuty
! Eq.  In ecref, 1st index = column; 2nd index = row of matrix.
! Mean atmosphere energy conversion term is energy conserving
!
   do k=1,plev
      ecref(k,k) = 0.5_r8/hypm(k) * hypd(k)
      do kk=1,k-1
         ecref(kk,k) = 1._r8/hypm(k) * hypd(kk)
      end do
   end do
!
! Reference hydrostatic integration matrix consistent with conversion
! term for energy conservation.  In href, 1st index = column; 
! 2nd index = row of matrix.
!
   do k = 1,plev
      do kk = k,plev
         href(kk,k) = ecref(k,kk)*hypd(kk)/hypd(k)
      end do
   end do
!
! Print statements
!
   if (lprint) then
      nkk = plev/13
      if (mod(plev,13).ne.0) nkk = nkk + 1
      write(iulog,*)' '
      write(iulog,*)'INITCOM: Hydrostatic matrix href'
      do kk=1,nkk
         ik1 = 1 + (kk-1)*13
         ik2 = min0( ik1+12, plev )
         write(iulog,9920) (k,k=ik1,ik2)
         do kkk=1,plev
            write(iulog,9910) kkk,(href(kkk,k),k=ik1,ik2)
         end do
      end do
      write(iulog,*)' '
      write(iulog,*)'INITCOM: Thermodynamic matrix ecref'
      do kk=1,nkk
         ik1 = 1 + (kk-1)*13
         ik2 = min0( ik1+12, plev )
         write(iulog,9920) (k,k=ik1,ik2)
         do kkk=1,plev
            write(iulog,9910) kkk,(ecref(kkk,k),k=ik1,ik2)
         end do
      end do
   end if
!
! Multiply href by r
!
   do k=1,plev
      do kk=1,plev
         href(kk,k) = href(kk,k)*rair
      end do
   end do
!
! Compute truncation parameters
!
   if (masterproc) then
      write(iulog,9950) ptrm,ptrn,ptrk
   end if
!
! Compute semi-implicit timestep constants (COMSPE)
!
   zdt = 0.5_r8*dtime
!
   call settau(zdt)
!
! Compute constants related to Legendre transforms
! Compute and reorder ALP and DALP
!
   do j=1,plat/2
      zslat = slat(j)
      itmp = 2*pspt - 1
      call phcs  (zalp    ,zdalp   ,itmp    ,zslat    )
      call reordp(j       ,itmp    ,zalp    ,zdalp   )
   end do

   call grmult_init()
!
! Determine whether full or reduced grid
!
   fullgrid = .true.
   do j=1,plat
      if (masterproc) then
         write(iulog,*)'nlon(',j,')=',nlon(j),' wnummax(',j,')=',wnummax(j)
      end if
      if (nlon(j).lt.plon) fullgrid = .false.
   end do
!
! Mirror latitudes south of south pole
!
   lat = 1
   do j=j1-2,1,-1
      nlonex(j) = nlon(lat)
      lat = lat + 1
   end do
   nlonex(j1-1) = nlon(1)     ! south pole
!
! Real latitudes
!
   j = j1
   do lat=1,plat
      nlonex(j) = nlon(lat)
      j = j + 1
   end do
   nlonex(j1+plat) = nlon(plat)  ! north pole
!
! Mirror latitudes north of north pole
!
   lat = plat
   do j=j1+plat+1,platd
      nlonex(j) = nlon(lat)
      lat = lat - 1
   end do
!
! Longitude array
!
   pi = 4.0_r8*atan(1.0_r8)
   do lat=1,plat
      do i=1,nlon(lat)
         londeg(i,lat) = (i-1)*360._r8/nlon(lat)
         clon(i,lat)   = (i-1)*2.0_r8*pi/nlon(lat)
      end do
   end do

   do j=1,plat/2
      nmmax(j) = wnummax(j) + 1
   end do
   do m=1,pmmax
      do irow=1,plat/2
         if (nmmax(irow) .ge. m) then
            beglatpair(m) = irow
            goto 10
         end if
      end do
      call endrun ('INITCOM: Should not ever get here')
10    continue
   end do
!
! Set up trigonometric tables for fft
!
   do j=1,plat
      call set99(trig(1,j),ifax(1,j),nlon(j))
   end do
!
! Set flag indicating dynamics grid is now defined.
! NOTE: this ASSUMES initcom is called after spmdinit.  The setting of nlon done here completes
! the definition of the dynamics grid.
!

   return

9910 format( 1x,i3,13f9.5)
9920 format(/,      13i9)
9950 format(/,'     Truncation Parameters',/,'     NTRM = ',i4,/, &
      '     NTRN = ',i4,/,'     NTRK = ',i4,/)

end subroutine initcom
