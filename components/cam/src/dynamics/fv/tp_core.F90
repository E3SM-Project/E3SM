#if defined( UNICOSMP ) || defined ( NEC_SX )
#define VECTORIZE
#endif
module tp_core
!BOP
!
! !MODULE: tp_core --- Utilities for the transport core
!
! !USES:
   use shr_kind_mod, only : r8 => shr_kind_r8

!
! !PUBLIC MEMBER FUNCTIONS:
      public tp2c, tp2d, xtp, xtpv, fxppm, xmist, steepx, lmppm
      public huynh, ytp, ymist, fyppm, tpcc, ycc
!
! !DESCRIPTION:
!
!      This module provides 
!
!      \begin{tabular}{|l|l|} \hline \hline
!       tp2c  &   \\ \hline
!       tp2d  &   \\ \hline 
!       xtp  &   \\ \hline 
!       fxppm  &   \\ \hline 
!       xmist  &   \\ \hline 
!       steepx  &   \\ \hline 
!       lmppm  &   \\ \hline 
!       huynh  &   \\ \hline 
!       ytp  &   \\ \hline 
!       ymist  &   \\ \hline 
!       fyppm  &   \\ \hline 
!       tpcc  &   \\ \hline 
!       ycc  &   \\ \hline
!                                \hline
!      \end{tabular}
!
! !REVISION HISTORY:
!   01.01.15   Lin        Routines coalesced into this module
!   01.03.26   Sawyer     Additional ProTeX documentation
!   03.11.19   Sawyer     Merged in CAM changes by Mirin
!   04.10.07   Sawyer     ompinner now from dynamics_vars
!   05.03.25   Todling    shr_kind_r8 can only be referenced once (MIPSpro-7.4.2)
!   05.05.25   Sawyer     Merged CAM and GEOS5 versions (mostly CAM)
!   06.09.06   Sawyer     Turned "magic numbers" into F90 parameters
!
!EOP
!-----------------------------------------------------------------------

! Magic numbers used in this module

   private
   real(r8), parameter ::  D0_0                    =  0.0_r8
   real(r8), parameter ::  D0_05                   =  0.05_r8
   real(r8), parameter ::  D0_25                   =  0.25_r8
   real(r8), parameter ::  D0_5                    =  0.5_r8
   real(r8), parameter ::  D1_0                    =  1.0_r8
   real(r8), parameter ::  D2_0                    =  2.0_r8
   real(r8), parameter ::  D3_0                    =  3.0_r8
   real(r8), parameter ::  D4_0                    =  4.0_r8
   real(r8), parameter ::  D8_0                    =  8.0_r8
   real(r8), parameter ::  D12_0                   = 12.0_r8
   real(r8), parameter ::  D24_0                   = 24.0_r8

CONTAINS

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: tp2c --- Perform transport on a C grid
!
! !INTERFACE: 
 subroutine tp2c(dh, va, h, crx, cry, im, jm,                            &
                 iord, jord, ng, fx, fy, ffsl,                           &
                 rcap, acosp, xfx, yfx, cosp, id, jfirst, jlast)
!-----------------------------------------------------------------------

 implicit none

! !INPUT PARAMETERS:
   integer im, jm                  ! Dimensions
   integer jfirst, jlast           ! Latitude strip
   integer iord, jord              ! Interpolation order in x,y
   integer ng                      ! Max. NS dependencies
   integer id                      ! density (0)  (mfx = C)
   real (r8) rcap                  ! Ask S.-J. (polar constant?)
   real (r8) acosp(jm)             ! Ask S.-J. (difference to cosp??)
   logical ffsl(jm)                ! Use flux-form semi-Lagrangian trans.?
                                      ! (N*NG S*NG)
   real (r8) cosp(jm)                   ! Critical angle
   real (r8) va(im,jfirst:jlast)        ! Courant  (unghosted)
   real (r8) h(im,jfirst-ng:jlast+ng)   ! Pressure ( N*NG S*NG )
   real (r8) crx(im,jfirst-ng:jlast+ng) ! Ask S.-J. ( N*NG S*NG )
   real (r8) cry(im,jfirst:jlast+1)     ! Ask S.-J. ( N like FY )
   real (r8) xfx(im,jfirst:jlast)       ! Ask S.-J. ( unghosted like FX )
   real (r8) yfx(im,jfirst:jlast+1)     ! Ask S.-J. ( N like FY )

! !OUTPUT PARAMETERS:
   real (r8) dh(im,jfirst:jlast)        ! Ask S.-J. ( unghosted )
   real (r8) fx(im,jfirst:jlast)        ! Flux in x ( unghosted )
   real (r8) fy(im,jfirst:jlast+1)      ! Flux in y ( N, see tp2c )

! !DESCRIPTION:
!     Perform transport on a C grid.   The number of ghost
!     latitudes (NG) depends on what method (JORD) will be used
!     subsequentally.    NG is equal to MIN(ABS(JORD),3).
!     Ask S.-J. how exactly this differs from TP2C.
!
! !REVISION HISTORY:
!
!EOP
!-----------------------------------------------------------------------
!BOC
   integer i, j, js2g0, jn2g0
   real (r8)  sum1

   js2g0  = max(2,jfirst)          !  No ghosting
   jn2g0  = min(jm-1,jlast)        !  No ghosting

   call tp2d(va, h, crx, cry, im, jm, iord, jord, ng,fx, fy, ffsl,    &
             xfx, yfx, cosp, id, jfirst, jlast)

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
   do j=js2g0,jn2g0
      do i=1,im-1
         dh(i,j) = fx(i,j) - fx(i+1,j) + (fy(i,j)-fy(i,j+1))*acosp(j)
      enddo
      dh(im,j) = fx(im,j) - fx(1,j) + (fy(im,j)-fy(im,j+1))*acosp(j)
   enddo

! Poles
   if ( jfirst ==  1 ) then
!       sum1 = - SUM( fy(1:im, 2) ) * rcap
        sum1 = D0_0
        do i=1,im
          sum1 = sum1 + fy(i,2)
        enddo
          sum1 = -sum1*rcap
        do i=1,im
          dh(i, 1) = sum1
        enddo
   endif
   
   if ( jlast == jm ) then
!       sum1 = SUM( fy(1:im,jm) ) * rcap
        sum1 = D0_0
        do i=1,im
          sum1 = sum1 + fy(i,jm)
        enddo
          sum1 = sum1*rcap
        do i=1,im
          dh(i,jm) = sum1
        enddo
   endif
   return
!EOC
 end subroutine tp2c
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: tp2d --- Perform transport on a D grid
!
! !INTERFACE: 
 subroutine tp2d(va, q, crx, cry, im, jm, iord, jord, ng, fx, fy,        &
                 ffsl, xfx, yfx, cosp, id, jfirst, jlast)
!-----------------------------------------------------------------------
! !USES:

 implicit none

! !INPUT PARAMETERS:
   integer im, jm                    ! Dimensions
   integer jfirst, jlast             ! Latitude strip
   integer iord, jord                ! Interpolation order in x,y
   integer ng                        ! Max. NS dependencies
   integer id                        ! density (0)  (mfx = C)
                                     ! mixing ratio (1) (mfx = mass flux)
   logical ffsl(jm)                  ! Use flux-form semi-Lagrangian trans.?
                                     ! ghosted N*ng S*ng
   real (r8) cosp(jm)                     ! Critical angle
   real (r8) va(im,jfirst:jlast)          ! Courant  (unghosted)
   real (r8) q(im,jfirst-ng:jlast+ng)     ! transported scalar ( N*NG S*NG )
   real (r8) crx(im,jfirst-ng:jlast+ng)   ! Ask S.-J. ( N*NG S*NG )
   real (r8) cry(im,jfirst:jlast+1)       ! Ask S.-J. ( N like FY )
   real (r8) xfx(im,jfirst:jlast)         ! Ask S.-J. ( unghosted like FX )
   real (r8) yfx(im,jfirst:jlast+1)       ! Ask S.-J. ( N like FY )

! !OUTPUT PARAMETERS:
   real (r8) fx(im,jfirst:jlast)          ! Flux in x ( unghosted )
   real (r8) fy(im,jfirst:jlast+1)        ! Flux in y ( N, see tp2c )

! !DESCRIPTION:
!     Perform transport on a D grid.   The number of ghost
!     latitudes (NG) depends on what method (JORD) will be used
!     subsequentally.    NG is equal to MIN(ABS(JORD),3).
!
!
! !REVISION HISTORY:
!   WS  99.04.13:  Added jfirst:jlast concept
!       99.04.21:  Removed j1 and j2 (j1=2, j2=jm-1 consistently)
!       99.04.27:  Removed dc, wk2 as arguments (local to YTP)
!       99.04.27:  Removed adx as arguments (local here)
!   SJL 99.07.26:  ffsl flag added
!   WS  99.09.07:  Restructuring, cleaning, documentation
!   WS  99.10.22:  NG now argument; arrays pruned
!   WS  00.05.14:  Renamed ghost indices as per Kevin's definitions
!
!EOP
!---------------------------------------------------------------------
!BOC

! Local:
   integer i, j, iad, jp, js2g0, js2gng, jn2g0, jn2gng
   real (r8) adx(im,jfirst-ng:jlast+ng)
   real (r8) wk1v(im,jfirst-ng:jlast+ng)
   real (r8)   dm(-im/3:im+im/3)
   real (r8) qtmpv(-im/3:im+im/3,jfirst-ng:jlast+ng)
   real (r8)   al(-im/3:im+im/3)
   real (r8)   ar(-im/3:im+im/3)
   real (r8)   a6(-im/3:im+im/3)

! Number of ghost latitudes
    js2g0  = max(2,jfirst)          !  No ghosting
    js2gng = max(2,jfirst-ng)       !  Number needed on S
    jn2g0  = min(jm-1,jlast)        !  No ghosting
    jn2gng = min(jm-1,jlast+ng)     !  Number needed on N
    iad = 1

    call xtpv(im,  ffsl, wk1v, q, crx, iad, crx,        &
             cosp, 0, dm, qtmpv, al, ar, a6,            &
             jfirst, jlast, js2gng, jn2gng, jm,         &
             1, jm, jfirst-ng, jlast+ng,                &
             jfirst-ng, jlast+ng, jfirst-ng, jlast+ng,  &
             jfirst-ng, jlast+ng, jfirst-ng, jlast+ng)

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
    do j=js2gng,jn2gng               !  adx needed on N*ng S*ng

       do i=1,im-1
          adx(i,j) = q(i,j) + D0_5 *                       &
                     (wk1v(i,j)-wk1v(i+1,j) + q(i,j)*(crx(i+1,j)-crx(i,j)))
       enddo
          adx(im,j) = q(im,j) + D0_5 *                     &
                      (wk1v(im,j)-wk1v(1,j) + q(im,j)*(crx(1,j)-crx(im,j)))
    enddo

! WS 99.09.07 : Split up north and south pole

     if ( jfirst-ng <= 1 ) then
        do i=1,im 
          adx(i, 1) = q(i,1)
        enddo
     endif 
     if ( jlast+ng >= jm ) then
        do i=1,im 
          adx(i,jm) = q(i,jm)
        enddo
     endif

     call ytp(im,jm,fy, adx,cry,yfx,ng,jord,0,jfirst,jlast)

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i,jp)
#endif
      do j=js2g0,jn2g0
        do i=1,im
           jp = j-va(i,j)
           wk1v(i,j) = q(i,j) +D0_5*va(i,j)*(q(i,jp)-q(i,jp+1))
        enddo
      enddo

      call xtpv(im,  ffsl, fx, wk1v, crx, iord, xfx,        &
               cosp, id, dm, qtmpv, al, ar, a6,             &
               jfirst, jlast, js2g0, jn2g0, jm,             &
               1, jm, jfirst, jlast,                        &
               jfirst-ng, jlast+ng, jfirst-ng, jlast+ng,    &
               jfirst, jlast, jfirst-ng, jlast+ng)

    return
!EOC
 end subroutine tp2d
!-----------------------------------------------------------------------

#ifndef VECTORIZE
!-----------------------------------------------------------------------
!BOP
! !IROUTINE: xtpv
!
! !INTERFACE: 
 subroutine xtpv(im, ffslv,  fxv,  qv,  cv,  iord,  mfxv,        &
                cosav, id, dmw, qtmpv, alw, arw, a6w,            &
                jfirst, jlast, jlow, jhigh, jm,                  &
                jl2, jh2, jl3, jh3,                              &
                jl4, jh4, jl5, jh5,                              &
                jl7, jh7, jl11, jh11)
!-----------------------------------------------------------------------

 implicit none
 
! !INPUT PARAMETERS:
   integer id               ! ID = 0: density (mfx = C)
                            ! ID = 1: mixing ratio (mfx is mass flux)

   integer im               ! Total longitudes
   integer iord
   integer jfirst, jlast, jlow, jhigh, jm
   integer jl2, jh2, jl3, jh3, jl4, jh4, jl5, jh5
   integer jl7, jh7, jl11, jh11 
   real (r8) cv(im,jl5:jh5)          ! Courant numbers
   real (r8) qv(im,jl4:jh4)
   real (r8) mfxv(im,jl7:jh7)
   logical ffslv(jl2:jh2)
   real (r8) cosav(jm)

! !INPUT/OUTPUT PARAMETERS:
   real (r8) qtmpv(-im/3:im+im/3,jl11:jh11)   ! Input work arrays:
   real (r8)   dmw(-im/3:im+im/3)
   real (r8)   alw(-im/3:im+im/3)
   real (r8)   arw(-im/3:im+im/3)
   real (r8)   a6w(-im/3:im+im/3)

! !OUTPUT PARAMETERS:
   real (r8) fxv(im,jl3:jh3)

! !DESCRIPTION:
!   
!
! !REVISION HISTORY:
!   99.01.01 Lin      Creation
!   01.03.27 Sawyer   Additional ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC

! Local:
   real (r8)       cos_upw               !critical cosine for upwind
   real (r8)       cos_van               !critical cosine for van Leer
   real (r8)       cos_ppm               !critical cosine for ppm

   parameter (cos_upw = D0_05)       !roughly at 87 deg.
   parameter (cos_van = D0_25)       !roughly at 75 deg.
   parameter (cos_ppm = D0_25)

   integer i, imp, j
   real (r8) qmax, qmin
   real (r8) rut, tmp
   integer iu, itmp, ist
   integer isave(im)
   integer iuw, iue
   real (r8) dm(-im/3:im+im/3)
   real (r8) al(-im/3:im+im/3)
   real (r8) ar(-im/3:im+im/3)
   real (r8) a6(-im/3:im+im/3)

   imp = im + 1

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i,iuw,iue,iu,itmp,isave,tmp,qmax,qmin,dm,rut,ist,al,ar,a6)
#endif
  do j = jlow, jhigh

   do i=1,im
      qtmpv(i,j) = qv(i,j)
   enddo

   if( ffslv(j) ) then
! Flux-Form Semi-Lagrangian transport

! Figure out ghost zone for the western edge:
      iuw =  -cv(1,j)
      iuw = min(0, iuw)
 
      do i=iuw, 0
         qtmpv(i,j) = qv(im+i,j)
      enddo 

! Figure out ghost zone for the eastern edge:
      iue = im - cv(im,j)
      iue = max(imp, iue)

      do i=imp, iue
         qtmpv(i,j) = qv(i-im,j)
      enddo

      if( iord == 1 .or. cosav(j) < cos_upw) then
      do i=1,im
        iu = cv(i,j)
      if(cv(i,j) .le. D0_0) then
        itmp = i - iu
        isave(i) = itmp - 1
      else
        itmp = i - iu - 1
        isave(i) = itmp + 1
      endif
        fxv(i,j) = (cv(i,j)-iu) * qtmpv(itmp,j)
      enddo
      else

      do i=1,im
! 2nd order slope
         tmp = D0_25*(qtmpv(i+1,j) - qtmpv(i-1,j))
         qmax = max(qtmpv(i-1,j), qtmpv(i,j), qtmpv(i+1,j)) - qtmpv(i,j)
         qmin = qtmpv(i,j) - min(qtmpv(i-1,j), qtmpv(i,j), qtmpv(i+1,j))
         dm(i) = sign(min(abs(tmp),qmax,qmin), tmp)
      enddo

 
      do i=iuw, 0
         dm(i) = dm(im+i)
      enddo 

      do i=imp, iue
         dm(i) = dm(i-im)
      enddo

      if(iord .ge. 3 .and. cosav(j) .gt. cos_ppm) then
         call fxppm(im, cv(:,j), mfxv(:,j), qtmpv(:,j), dm, fxv(:,j), iord, al, ar, a6,         &
                    iuw, iue, ffslv(j), isave)
      else
      do i=1,im
            iu  = cv(i,j)
            rut = cv(i,j) - iu
         if(cv(i,j) .le. D0_0) then
            itmp = i - iu
            isave(i) = itmp - 1
            fxv(i,j) = rut*(qtmpv(itmp,j)-dm(itmp)*(D1_0+rut))
         else
            itmp = i - iu - 1
            isave(i) = itmp + 1
            fxv(i,j) = rut*(qtmpv(itmp,j)+dm(itmp)*(D1_0-rut))
         endif
      enddo
      endif

      endif

      do i=1,im
      if(cv(i,j) .ge. D1_0) then
        do ist = isave(i),i-1
           fxv(i,j) = fxv(i,j) + qtmpv(ist,j)
        enddo
      elseif(cv(i,j) .le. -D1_0) then
        do ist = i,isave(i)
           fxv(i,j) = fxv(i,j) - qtmpv(ist,j)
        enddo
      endif
      enddo

      if(id .ne. 0) then
         do i=1,im
            fxv(i,j) =  fxv(i,j)*mfxv(i,j)
         enddo
      endif

   else
! Regular PPM (Eulerian without FFSL extension)

      qtmpv(imp,j) = qv(1,j)
      qtmpv(  0,j) = qv(im,j)

      if(iord == 1 .or. cosav(j) < cos_upw) then
         do i=1,im
            iu = real(i,r8) - cv(i,j)
            fxv(i,j) = mfxv(i,j)*qtmpv(iu,j)
         enddo
      else

         qtmpv(-1,j)    = qv(im-1,j)
         qtmpv(imp+1,j) = qv(2,j)

         if(iord > 0 .or. cosav(j) < cos_van) then
            call xmist(im, qtmpv(:,j), dm, 2)
         else
            call xmist(im, qtmpv(:,j), dm, iord)
         endif

         dm(0) = dm(im)

         if( abs(iord).eq.2 .or. cosav(j) .lt. cos_van ) then
            do i=1,im
               iu = real(i,r8) - cv(i,j)
               fxv(i,j) =  mfxv(i,j)*(qtmpv(iu,j)+dm(iu)*(sign(D1_0,cv(i,j))-cv(i,j)))

!              if(cv(i,j) .le. 0.) then
!                 fxv(i,j) = qtmpv(i,j) - dm(i)*(1.+cv(i,j))
!              else
!                 fxv(i,j) = qtmpv(i-1,j) + dm(i-1)*(1.-cv(i,j))
!              endif
!                 fxv(i,j) = fxv(i,j)*mfxv(i,j)

            enddo
         else
            call fxppm(im, cv(:,j), mfxv(:,j), qtmpv(:,j), dm, fxv(:,j), iord, al, ar, a6,       &
                       iuw, iue, ffslv(j), isave)
         endif
      endif

   endif

  enddo

   return
!EOC
 end subroutine xtpv
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: xmist
!
! !INTERFACE: 
 subroutine xmist(im,  q,  dm,  id)
!-----------------------------------------------------------------------

 implicit none

! !INPUT PARAMETERS:
 integer im                   ! Total number of longitudes
 integer id                   ! ID = 0: density (mfx = C)
                              ! ID = 1: mixing ratio (mfx is mass flux)
 real(r8)  q(-im/3:im+im/3)   ! Input latitude 

! !OUTPUT PARAMETERS:
 real(r8) dm(-im/3:im+im/3)   ! 

! !DESCRIPTION:
!   
!
! !REVISION HISTORY:
!   99.01.01 Lin      Creation
!   01.03.27 Sawyer   Additional ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC

 real(r8) r24
 parameter( r24 = D1_0/D24_0)

 integer i
 real(r8) qmin, qmax

    if(id .le. 2) then
       do i=1,im
          dm(i) = r24*(D8_0*(q(i+1) - q(i-1)) + q(i-2) - q(i+2))
       enddo
    else
       do i=1,im
          dm(i) = D0_25*(q(i+1) - q(i-1))
       enddo
    endif

    if( id < 0 ) return

! Apply monotonicity constraint (Lin et al. 1994, MWR)
      do i=1,im
         qmax = max( q(i-1), q(i), q(i+1) ) - q(i)
         qmin = q(i) - min( q(i-1), q(i), q(i+1) )
         dm(i) = sign( min(abs(dm(i)), qmax, qmin), dm(i) )
      enddo
  return
!EOC
 end subroutine xmist
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: fxppm
!
! !INTERFACE: 
 subroutine fxppm(im, c, mfx,  p, dm, fx, iord, al, ar, a6,        &
                  iuw, iue, ffsl, isave)
!-----------------------------------------------------------------------
!
! !USES:
 implicit none

! !INPUT PARAMETERS:
 integer im, iord
 real (r8)  c(im)
 real (r8) p(-im/3:im+im/3)
 real (r8) dm(-im/3:im+im/3)
 real (r8) mfx(im)
 integer iuw, iue
 logical ffsl
 integer isave(im)

! !INPUT/OUTPUT PARAMETERS:
 real (r8) al(-im/3:im+im/3)
 real (r8) ar(-im/3:im+im/3)
 real (r8) a6(-im/3:im+im/3)

! !OUTPUT PARAMETERS:
 real (r8) fx(im)

! !DESCRIPTION:
!   
!
! !REVISION HISTORY:
!   99.01.01 Lin      Creation
!   01.03.27 Sawyer   Additional ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
 real (r8) r3, r23
 parameter ( r3 = D1_0/D3_0, r23 = D2_0/D3_0 )

 integer i, lmt
 integer iu, itmp
 real (r8) ru
 logical steep

  if( iord == 6 ) then
      steep = .true.
  else
      steep = .false.
  endif

  do i=1,im
     al(i) = D0_5*(p(i-1)+p(i)) + (dm(i-1) - dm(i))*r3
  enddo

  if( steep ) call steepx( im, p, al(1), dm )

     do i=1,im-1
        ar(i) = al(i+1)
     enddo
        ar(im) = al(1)

  if(iord == 7) then
     call huynh(im, ar(1), al(1), p(1), a6(1), dm(1))
  else
     if(iord .eq. 3 .or. iord .eq. 5) then
         do i=1,im
            a6(i) = D3_0*(p(i)+p(i)  - (al(i)+ar(i)))
         enddo
     endif
     lmt = iord - 3
     call lmppm( dm(1), a6(1), ar(1), al(1), p(1), im, lmt )
  endif

  if( ffsl ) then

      do i=iuw, 0
         al(i) = al(im+i)
         ar(i) = ar(im+i)
         a6(i) = a6(im+i)
      enddo

      do i=im+1, iue
         al(i) = al(i-im)
         ar(i) = ar(i-im)
         a6(i) = a6(i-im)
      enddo

      do i=1,im
            iu = c(i)
            ru = c(i) - iu
         if(c(i) .gt. D0_0) then
            itmp = i - iu - 1
            isave(i) = itmp + 1
            fx(i) = ru*(ar(itmp)+D0_5*ru*(al(itmp)-ar(itmp) +     &
                        a6(itmp)*(D1_0-r23*ru)) )
         else
            itmp = i - iu
            isave(i) = itmp - 1
            fx(i) = ru*(al(itmp)-D0_5*ru*(ar(itmp)-al(itmp) +     &
                        a6(itmp)*(D1_0+r23*ru)) )
         endif
      enddo

  else
         al(0) = al(im)
         ar(0) = ar(im)
         a6(0) = a6(im)
      do i=1,im
         if(c(i) .gt. D0_0) then
            fx(i) = ar(i-1) + D0_5*c(i)*(al(i-1) - ar(i-1) +   &
                    a6(i-1)*(D1_0-r23*c(i)) )
      else
            fx(i) = al(i) - D0_5*c(i)*(ar(i) - al(i) +         &
                    a6(i)*(D1_0+r23*c(i)))
      endif
            fx(i) = mfx(i) * fx(i)
      enddo
  endif
  return
!EOC
 end subroutine fxppm
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: steepx
!
! !INTERFACE: 
 subroutine  steepx(im, p, al, dm)
!-----------------------------------------------------------------------

! !USES:
 implicit none

! !INPUT PARAMETERS:
 integer im
 real (r8)  p(-im/3:im+im/3)
 real (r8) dm(-im/3:im+im/3)

! !INPUT/OUTPUT PARAMETERS:
 real (r8) al(im)

! !DESCRIPTION:
!   
!
! !REVISION HISTORY:
!   99.01.01 Lin      Creation
!   01.03.27 Sawyer   Additional ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
 integer i
 real (r8) r3
 parameter ( r3 = D1_0/D3_0 )

 real (r8) dh(0:im)
 real (r8) d2(0:im+1)
 real (r8) eta(0:im)
 real (r8) xxx, bbb, ccc

   do i=0,im
      dh(i) = p(i+1) - p(i)
   enddo

! Needs dh(0:im)
   do i=1,im
      d2(i) = dh(i) - dh(i-1)
   enddo
   d2(0) = d2(im)
   d2(im+1) = d2(1)

! needs p(-1:im+2), d2(0:im+1)
   do i=1,im
      if( d2(i+1)*d2(i-1).lt.D0_0 .and. p(i+1).ne.p(i-1) ) then
          xxx    = D1_0 - D0_5 * ( p(i+2) - p(i-2) ) / ( p(i+1) - p(i-1) )
          eta(i) = max(D0_0, min(xxx, D0_5) )
      else
          eta(i) = D0_0
      endif
    enddo

    eta(0) = eta(im)

! needs eta(0:im), dh(0:im-1), dm(0:im)
   do i=1,im
      bbb = ( D2_0*eta(i  ) - eta(i-1) ) * dm(i-1) 
      ccc = ( D2_0*eta(i-1) - eta(i  ) ) * dm(i  ) 
      al(i) = al(i) + D0_5*( eta(i-1) - eta(i)) * dh(i-1) + (bbb - ccc) * r3
   enddo
   return
!EOC
 end subroutine steepx
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: lmppm
!
! !INTERFACE: 
 subroutine lmppm(dm, a6, ar, al, p, im, lmt)
!-----------------------------------------------------------------------

 implicit none

! !INPUT PARAMETERS:
 integer im   ! Total longitudes
 integer lmt  ! LMT = 0: full monotonicity
              ! LMT = 1: Improved and simplified full monotonic constraint
              ! LMT = 2: positive-definite constraint
              ! LMT = 3: Quasi-monotone constraint
 real(r8) p(im)
 real(r8) dm(im)

! !OUTPUT PARAMETERS:
 real(r8) a6(im)
 real(r8) ar(im)
 real(r8) al(im)

! !DESCRIPTION:
!
!
! !REVISION HISTORY:
!   99.01.01 Lin      Creation
!   01.03.27 Sawyer   Additional ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
 real (r8) r12
 parameter ( r12 = D1_0/D12_0 )

 real (r8) da1, da2, fmin, a6da
 real (r8) dr, dl

 integer i

! LMT = 0: full monotonicity
! LMT = 1: Improved and simplified full monotonic constraint
! LMT = 2: positive-definite constraint
! LMT = 3: Quasi-monotone constraint

  if( lmt == 0 ) then

! Full constraint
  do i=1,im
     if(dm(i) .eq. D0_0) then
         ar(i) = p(i)
         al(i) = p(i)
         a6(i) = D0_0
     else
         da1  = ar(i) - al(i)
         da2  = da1**2
         a6da = a6(i)*da1
         if(a6da .lt. -da2) then
            a6(i) = D3_0*(al(i)-p(i))
            ar(i) = al(i) - a6(i)
         elseif(a6da .gt. da2) then
            a6(i) = D3_0*(ar(i)-p(i))
            al(i) = ar(i) - a6(i)
         endif
     endif
  enddo

  elseif( lmt == 1 ) then

! Improved (Lin 2001?) full constraint
      do i=1,im
           da1 = dm(i) + dm(i)
            dl = sign(min(abs(da1),abs(al(i)-p(i))), da1)
            dr = sign(min(abs(da1),abs(ar(i)-p(i))), da1)
         ar(i) = p(i) + dr
         al(i) = p(i) - dl
         a6(i) = D3_0*(dl-dr)
      enddo

  elseif( lmt == 2 ) then
! Positive definite constraint
      do 250 i=1,im
      if(abs(ar(i)-al(i)) .ge. -a6(i)) go to 250
      fmin = p(i) + D0_25*(ar(i)-al(i))**2/a6(i) + a6(i)*r12
      if(fmin.ge.D0_0) go to 250
      if(p(i).lt.ar(i) .and. p(i).lt.al(i)) then
            ar(i) = p(i)
            al(i) = p(i)
            a6(i) = D0_0
      elseif(ar(i) .gt. al(i)) then
            a6(i) = D3_0*(al(i)-p(i))
            ar(i) = al(i) - a6(i)
      else
            a6(i) = D3_0*(ar(i)-p(i))
            al(i) = ar(i) - a6(i)
      endif
250   continue

  elseif(lmt .eq. 3) then
! Quasi-monotone constraint
      do i=1,im
         da1 = D4_0*dm(i)
          dl = sign(min(abs(da1),abs(al(i)-p(i))), da1)
          dr = sign(min(abs(da1),abs(ar(i)-p(i))), da1)
         ar(i) = p(i) + dr
         al(i) = p(i) - dl
         a6(i) = D3_0*(dl-dr)
      enddo
  endif
  return
!EOC
 end subroutine lmppm
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: huynh --- Enforce Huynh's 2nd constraint in 1D periodic domain
!
! !INTERFACE: 
 subroutine huynh(im, ar, al, p, d2, d1)
!-----------------------------------------------------------------------

! !USES:

 implicit none

! !INPUT PARAMETERS:
 integer im
 real(r8)  p(im)

! !OUTPUT PARAMETERS:
 real(r8) ar(im)
 real(r8) al(im)
 real(r8) d2(im)
 real(r8) d1(im)

! !DESCRIPTION:
!
!   Enforce Huynh's 2nd constraint in 1D periodic domain
!
! !REVISION HISTORY:
!   99.01.01 Lin      Creation
!   01.03.27 Sawyer   Additional ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
 integer  i
 real(r8) pmp
 real(r8) lac
 real(r8) pmin
 real(r8) pmax

! Compute d1 and d2
      d1(1) = p(1) - p(im)
      do i=2,im
         d1(i) = p(i) - p(i-1)
      enddo

      do i=1,im-1
         d2(i) = d1(i+1) - d1(i)
      enddo
      d2(im) = d1(1) - d1(im)

! Constraint for AR
!            i = 1
         pmp   = p(1) + D2_0 * d1(1)
         lac   = p(1) + D0_5 * (d1(1)+d2(im)) + d2(im) 
         pmin  = min(p(1), pmp, lac)
         pmax  = max(p(1), pmp, lac)
         ar(1) = min(pmax, max(ar(1), pmin))

      do i=2, im
         pmp   = p(i) + D2_0*d1(i)
         lac   = p(i) + D0_5*(d1(i)+d2(i-1)) + d2(i-1)
         pmin  = min(p(i), pmp, lac)
         pmax  = max(p(i), pmp, lac)
         ar(i) = min(pmax, max(ar(i), pmin))
      enddo

! Constraint for AL
      do i=1, im-1
         pmp   = p(i) - D2_0*d1(i+1)
         lac   = p(i) + D0_5*(d2(i+1)-d1(i+1)) + d2(i+1)
         pmin  = min(p(i), pmp, lac)
         pmax  = max(p(i), pmp, lac)
         al(i) = min(pmax, max(al(i), pmin))
      enddo

! i=im
         i = im
         pmp    = p(im) - D2_0*d1(1)
         lac    = p(im) + D0_5*(d2(1)-d1(1)) + d2(1)
         pmin   = min(p(im), pmp, lac)
         pmax   = max(p(im), pmp, lac)
         al(im) = min(pmax, max(al(im), pmin))

! compute A6 (d2)
      do i=1, im
         d2(i) = D3_0*(p(i)+p(i)  - (al(i)+ar(i)))
      enddo
    return
!EOC
 end subroutine huynh
!-----------------------------------------------------------------------
#endif

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ytp
!
! !INTERFACE: 
 subroutine ytp(im, jm, fy, q, c, yfx, ng, jord, iv, jfirst, jlast)
!-----------------------------------------------------------------------

! !USES:
 implicit none

! !INPUT PARAMETERS:
 integer im, jm                      !  Dimensions
 integer jfirst, jlast               !  Latitude strip
 integer ng                          !  Max. NS dependencies
 integer jord                        !  order of subgrid dist
 integer iv                          !  Scalar=0, Vector=1
 real (r8) q(im,jfirst-ng:jlast+ng)       !  advected scalar N*jord S*jord
 real (r8) c(im,jfirst:jlast+1)           !  Courant   N (like FY)
 real (r8) yfx(im,jfirst:jlast+1)         !  Backgrond mass flux

! !OUTPUT PARAMETERS:
 real (r8) fy(im,jfirst:jlast+1)          !  Flux      N (see tp2c)

! !DESCRIPTION:
!     This routine calculates the flux FX.  The method chosen
!     depends on the order of the calculation JORD (currently
!     1, 2 or 3).  
!
! !CALLED FROM:
!     cd_core
!     tp2d
!
! !REVISION HISTORY:
!
!  SJL 99.04.13:  Delivery
!  WS  99.04.13:  Added jfirst:jlast concept
!  WS  99.04.21:  Removed j1 and j2 (j1=2, j2=jm-1 consistently)
!                 removed a6,ar,al from argument list
!  WS  99.04.27:  DM made local to this routine
!  WS  99.09.09:  Documentation; indentation; cleaning
!  WS  99.10.22:  Added NG as argument; pruned arrays
!  SJL 99.12.24:  Revised documentation; optimized for better cache usage
!  WS  00.05.14:  Renamed ghost indices as per Kevin's definitions
!
!EOP
!---------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
 integer i, j, jt
 integer js2g0, jn1g1

! work arrays (should pass in eventually for performance enhancement):
 real (r8) dm(im,jfirst-ng:jlast+ng)

!     real (r8) ar(im,jfirst-1:jlast+1)  ! AR needs to be ghosted on NS
!     real (r8) al(im,jfirst-1:jlast+2)  ! AL needs to be ghosted on N2S
!     real (r8) a6(im,jfirst-1:jlast+1)  ! A6 needs to be ghosted on NS

   js2g0  = max(2,jfirst)       ! No ghosting
   jn1g1  = min(jm,jlast+1)     ! Ghost N*1
     
   if(jord == 1) then
#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i,jt)
#endif
        do j=js2g0,jn1g1
          do i=1,im
            jt = real(j,r8) - c(i,j)
            fy(i,j) = q(i,jt)
          enddo
        enddo
   else

!
! YMIST requires q on NS;  Only call to YMIST here
!
        call ymist(im, jm, q, dm, ng, jord, iv, jfirst, jlast)

        if( abs(jord) .ge. 3 ) then
 
          call fyppm(c,q,dm,fy,im,jm,ng,jord,iv,jfirst,jlast)

        else
!
! JORD can either have the value 2 or -2 at this point
!
#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i,jt)
#endif
          do j=js2g0,jn1g1
            do i=1,im
              jt = real(j,r8) - c(i,j)
              fy(i,j) = q(i,jt) + (sign(D1_0,c(i,j))-c(i,j))*dm(i,jt)
            enddo
          enddo
        endif
   endif

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
      do j=js2g0,jn1g1
        do i=1,im
          fy(i,j) = fy(i,j)*yfx(i,j)
        enddo
      enddo
    return
!EOC
 end subroutine ytp
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ymist
!
! !INTERFACE: 
 subroutine ymist(im, jm, q, dm, ng, jord, iv, jfirst, jlast)
!-----------------------------------------------------------------------

! !USES:
 implicit none

! !INPUT PARAMETERS:
 integer im, jm                      !  Dimensions
 integer jfirst, jlast               !  Latitude strip
 integer ng                          !  NS dependencies
 integer jord                        !  order of subgrid distribution
 integer iv                          !  Scalar (==0) Vector (==1)
 real (r8) q(im,jfirst-ng:jlast+ng)  !  transported scalar  N*ng S*ng

! !OUTPUT PARAMETERS:
 real (r8) dm(im,jfirst-ng:jlast+ng)      !  Slope only N*(ng-1) S*(ng-1) used

! !DESCRIPTION:
!     Calculate the slope of the pressure.  The number of ghost
!     latitudes (NG) depends on what method (JORD) will be used
!     subsequentally.    NG is equal to MIN(ABS(JORD),3).
!
! !CALLED FROM:
!     ytp
!
! !REVISION HISTORY:
!
!  SJL 99.04.13:  Delivery
!  WS  99.04.13:  Added jfirst:jlast concept
!  WS  99.09.09:  Documentation; indentation; cleaning
!  SJL 00.01.06:  Documentation
!  WS  00.05.14:  Renamed ghost indices as per Kevin's definitions
!
!EOP
!---------------------------------------------------------------------
!BOC

! Local variables

 integer i, j, jm1, im2, js2gng1, jn2gng1
 real (r8) qmax, qmin, tmp

    js2gng1 = max(2,   jfirst-ng+1)     !  Number needed on S
    jn2gng1 = min(jm-1,jlast+ng-1)      !  Number needed on N

    jm1 = jm - 1
    im2 = im / 2

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
      do j=js2gng1,jn2gng1
        do i=1,im
           dm(i,j) = D0_25*(q(i,j+1) - q(i,j-1))
        enddo
      enddo

   if( iv == 0 ) then

        if ( jfirst-ng <= 1 ) then
! S pole
          do i=1,im2
            tmp = D0_25*(q(i,2)-q(i+im2,2))
            qmax = max(q(i,2),q(i,1), q(i+im2,2)) - q(i,1)
            qmin = q(i,1) - min(q(i,2),q(i,1), q(i+im2,2))
            dm(i,1) = sign(min(abs(tmp),qmax,qmin),tmp)
          enddo

          do i=im2+1,im
            dm(i, 1) =  - dm(i-im2, 1)
          enddo
        endif

        if ( jlast+ng >= jm ) then
! N pole
          do i=1,im2
            tmp = D0_25*(q(i+im2,jm1)-q(i,jm1))
            qmax = max(q(i+im2,jm1),q(i,jm), q(i,jm1)) - q(i,jm)
            qmin = q(i,jm) - min(q(i+im2,jm1),q(i,jm), q(i,jm1))
            dm(i,jm) = sign(min(abs(tmp),qmax,qmin),tmp)
          enddo

          do i=im2+1,im
            dm(i,jm) =  - dm(i-im2,jm)
          enddo
        endif

   else

        if ( jfirst-ng <= 1 ) then
! South
          do i=1,im2
            tmp  = D0_25*(q(i,2)+q(i+im2,2))
            qmax = max(q(i,2),q(i,1), -q(i+im2,2)) - q(i,1)
            qmin = q(i,1) - min(q(i,2),q(i,1),-q(i+im2,2))
            dm(i,1) = sign(min(abs(tmp),qmax,qmin),tmp)
          enddo

          do i=im2+1,im
            dm(i, 1) = dm(i-im2, 1)
          enddo
        endif

        if ( jlast+ng >= jm ) then
! North
          do i=1,im2
            tmp  = -D0_25*(q(i+im2,jm1)+q(i,jm1))
            qmax = max(-q(i+im2,jm1),q(i,jm), q(i,jm1)) - q(i,jm)
            qmin = q(i,jm) - min(-q(i+im2,jm1),q(i,jm), q(i,jm1))
            dm(i,jm) = sign(min(abs(tmp),qmax,qmin),tmp)
          enddo

          do i=im2+1,im
            dm(i,jm) = dm(i-im2,jm)
          enddo
        endif

   endif

   if( jord > 0 ) then
!
! Applies monotonic slope constraint (off if jord less than zero)
!
#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i,qmax,qmin)
#endif
        do j=js2gng1,jn2gng1
          do i=1,im
            qmax = max(q(i,j-1),q(i,j),q(i,j+1)) - q(i,j)
            qmin = q(i,j) - min(q(i,j-1),q(i,j),q(i,j+1))
            dm(i,j) = sign(min(abs(dm(i,j)),qmin,qmax),dm(i,j))
          enddo
        enddo
   endif
    return
!EOC
 end subroutine ymist
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: fyppm
!
! !INTERFACE: 
 subroutine fyppm(c,  q,  dm, flux, im, jm, ng, jord, iv, jfirst, jlast)
!-----------------------------------------------------------------------

! !USES:
 implicit none

! !INPUT PARAMETERS:
 integer im, jm                      !  Dimensions
 integer jfirst, jlast               !  Latitude strip
 integer ng                          !  Max. NS dependencies
 integer jord                        !  Approximation order
 integer iv                          !  Scalar=0, Vector=1
 real (r8)  q(im,jfirst-ng:jlast+ng) !  mean value needed only N*2 S*2
 real (r8) dm(im,jfirst-ng:jlast+ng) !  Slope     needed only N*2 S*2
 real (r8)  c(im,jfirst:jlast+1)     !  Courant   N (like FLUX)

! !INPUT/OUTPUT PARAMETERS:
 real (r8) ar(im,jfirst-1:jlast+1)   ! AR needs to be ghosted on NS
 real (r8) al(im,jfirst-1:jlast+2)   ! AL needs to be ghosted on N2S
 real (r8) a6(im,jfirst-1:jlast+1)   ! A6 needs to be ghosted on NS

! !OUTPUT PARAMETERS:
 real (r8) flux(im,jfirst:jlast+1)   !  Flux      N (see tp2c)

! !DESCRIPTION:
!
!   NG is passed from YTP for convenience -- it is actually 1 more in NS
!   than the actual number of latitudes needed here.  But in the shared-memory 
!   case it becomes 0, which is much cleaner.
!
! !CALLED FROM:
!      ytp
!
! !REVISION HISTORY:
!
!  SJL 99.04.13:  Delivery
!  WS  99.04.19:  Added jfirst:jlast concept; FYPPM only called from YTP
!  WS  99.04.21:  Removed j1, j2  (j1=2, j2=jm-1 consistently)
!                 removed a6,ar,al from argument list
!  WS  99.09.09:  Documentation; indentation; cleaning
!  WS  99.10.22:  Added ng as argument; Pruned arrays
!  WS  00.05.14:  Renamed ghost indices as per Kevin's definitions
!
!EOP
!---------------------------------------------------------------------
!BOC
 real (r8)   r3, r23
 parameter ( r3 = D1_0/D3_0, r23 = D2_0/D3_0 )
 integer i, j, imh, jm1, lmt
 integer js1g1, js2g0, js2g1, jn1g2, jn1g1, jn2g1
 integer jan, jlow, jhigh, ilow, ihigh
 integer ja(jlast-jfirst+3)
!     logical steep

!     if(jord .eq. 6) then
!        steep = .true.
!     else
!        steep = .false.
!     endif

      imh = im / 2
      jm1 = jm - 1

      js1g1  = max(1,jfirst-1)         ! Ghost S*1
      js2g0  = max(2,jfirst)           ! No ghosting
      js2g1  = max(2,jfirst-1)         ! Ghost S*1
      jn1g1  = min(jm,jlast+1)         ! Ghost N*1
      jn1g2  = min(jm,jlast+2)         ! Ghost N*2
      jn2g1  = min(jm-1,jlast+1)       ! Ghost N*1

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
      do j=js2g1,jn1g2                 ! AL needed N2S
        do i=1,im                      ! P, dm ghosted N2S2 (at least)
          al(i,j) = D0_5*(q(i,j-1)+q(i,j)) + r3*(dm(i,j-1) - dm(i,j))
        enddo
      enddo

! Yeh's steepening procedure; to be implemented
!     if(steep) call steepy(im,   jm,   jfirst,   jlast,       &
!                           ng,    q,       al,   dm )

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
      do j=js1g1,jn2g1                 ! AR needed NS
        do i=1,im
          ar(i,j) = al(i,j+1)          ! AL ghosted N2S
        enddo
      enddo

! WS 990726 :  Added condition to decide if poles are on this processor

! Poles:

   if( iv == 0 ) then

        if ( jfirst == 1 ) then
          do i=1,imh
            al(i,    1) = al(i+imh,2)
            al(i+imh,1) = al(i,    2)
          enddo
        endif

        if ( jlast == jm ) then
          do i=1,imh
            ar(i,    jm) = ar(i+imh,jm1)
            ar(i+imh,jm) = ar(i,    jm1)
          enddo
        endif

   else

        if ( jfirst == 1 ) then
          do i=1,imh
            al(i,    1) = -al(i+imh,2)
            al(i+imh,1) = -al(i,    2)
          enddo
        endif

        if ( jlast == jm ) then
          do i=1,imh
            ar(i,    jm) = -ar(i+imh,jm1)
            ar(i+imh,jm) = -ar(i,    jm1)
          enddo
        endif

   endif

   if( jord == 3 .or. jord == 5 ) then
#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
      do j=js1g1,jn1g1               ! A6 needed NS
        do i=1,im
          a6(i,j) = D3_0*(q(i,j)+q(i,j) - (al(i,j)+ar(i,j)))
        enddo
      enddo
   endif

      lmt = jord - 3

!       do j=js1g1,jn1g1             !  A6, AR, AL needed NS
!         call lmppm(dm(1,j),a6(1,j),ar(1,j),al(1,j),q(1,j),im,lmt)
!       enddo

#ifdef VECTORIZE
        jan = 1
        ja(1) = 1
        ilow = 1
        ihigh = im*(jn1g1-js1g1+1)
        jlow = 1
        jhigh = 1
        call lmppmv(dm(1,js1g1), a6(1,js1g1), ar(1,js1g1),               &
                   al(1,js1g1),  q(1,js1g1), im*(jn1g1-js1g1+1), lmt,    &
                   jan, ja, ilow, ihigh, jlow, jhigh, jlow, jhigh)
#else
        call lmppm(dm(1,js1g1), a6(1,js1g1), ar(1,js1g1),               &
                   al(1,js1g1),  q(1,js1g1), im*(jn1g1-js1g1+1), lmt)
#endif

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
      do j=js2g0,jn1g1                 ! flux needed N
        do i=1,im
          if(c(i,j).gt.D0_0) then
            flux(i,j) = ar(i,j-1) + D0_5*c(i,j)*(al(i,j-1) - ar(i,j-1) +  &
                        a6(i,j-1)*(D1_0-r23*c(i,j)) )
          else
            flux(i,j) = al(i,j) - D0_5*c(i,j)*(ar(i,j) - al(i,j) +        &
                        a6(i,j)*(D1_0+r23*c(i,j)))
          endif
        enddo
      enddo
    return
!EOC
 end subroutine fyppm 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: tpcc
!
! !INTERFACE: 
 subroutine tpcc(va,   ymass,  q,   crx,  cry,  im,   jm,  ng_c, ng_d,   &
                 iord, jord,   fx,  fy,   ffsl, cose, jfirst, jlast,     &
                 dm,   qtmp,   al,  ar,   a6 )       
!-----------------------------------------------------------------------

! !USES:
 implicit none

! !INPUT PARAMETERS:
  integer im, jm                    ! Dimensions
  integer ng_c                      ! 
  integer ng_d                      ! 
  integer jfirst, jlast             ! Latitude strip
  integer iord, jord                ! Interpolation order in x,y
  logical ffsl(jm)                  ! Flux-form semi-Lagrangian transport?
  real (r8) cose(jm)                ! Critical cosine  (replicated)
  real (r8) va(im,jfirst:jlast)     ! Courant (unghosted like FX)
  real (r8) q(im,jfirst-ng_d:jlast+ng_d) !
  real (r8) crx(im,jfirst-ng_c:jlast+ng_c)
  real (r8) cry(im,jfirst:jlast)    ! Courant # (ghosted like FY)
  real (r8) ymass(im,jfirst:jlast)  ! Background y-mass-flux (ghosted like FY)

! Input 1D work arrays:
  real (r8)   dm(-im/3:im+im/3)
  real (r8) qtmp(-im/3:im+im/3)
  real (r8)   al(-im/3:im+im/3)
  real (r8)   ar(-im/3:im+im/3)
  real (r8)   a6(-im/3:im+im/3)

! !OUTPUT PARAMETERS:
  real (r8) fx(im,jfirst:jlast)     ! Flux in x (unghosted)
  real (r8) fy(im,jfirst:jlast)     ! Flux in y (unghosted since iv==0)

! !DESCRIPTION:
!     In this routine the number 
!     of north ghosted latitude min(abs(jord),2), and south ghosted
!     latitudes is XXXX
!
! !CALLED FROM:
!     cd_core
!
! !REVISION HISTORY:
!   SJL 99.04.13:  Delivery
!   WS  99.04.13:  Added jfirst:jlast concept
!   WS  99.05.10:  Replaced JNP with JM, JMR with JM-1, IMR with IM
!   WS  99.05.10:  Removed fvcore.h and JNP, IMH, IML definitions
!   WS  99.10.20:  Pruned arrays
!   WS  00.05.14:  Renamed ghost indices as per Kevin's definitions
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

  real (r8) adx(im,jfirst-1:jlast+2)
  integer north, south
  integer i, j, jp, im2, js2g0, js2gs, jn2g0, jn1g0, jn1gn
  real (r8) wk1v(im,jfirst-1:jlast+2)
  real (r8) fx1(im)
  real (r8) qtmpv(-im/3:im+im/3,jfirst-1:jlast+2)

    im2 = im/2
    north = min(2,abs(jord))         ! north == 1 or 2
    south = north-1                  ! south == 0 or 1
    js2g0 = max(2,jfirst)
    js2gs = max(2,jfirst-south)
    jn2g0 = min(jm-1,jlast)
    jn1gn = min(jm,jlast+north)
    jn1g0 = min(jm,jlast)

! This loop must be ghosted N*NG, S*NG

    call xtpv( im, ffsl, wk1v, q, crx, 1, crx,                  &
             cose, 0, dm, qtmpv, al, ar, a6,                    &
             jfirst, jlast, js2gs, jn1gn, jm,                   &
             1, jm, jfirst-1, jlast+2,                          &
             jfirst-ng_d, jlast+ng_d, jfirst-ng_c, jlast+ng_c,  &
             jfirst-ng_c, jlast+ng_c, jfirst-1, jlast+2)

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
      do j=js2gs,jn1gn

        do i=1,im-1
          adx(i,j) = q(i,j) + D0_5 *                   &
                     (wk1v(i,j)-wk1v(i+1,j) + q(i,j)*(crx(i+1,j)-crx(i,j)))
        enddo

        adx(im,j) = q(im,j) + D0_5 *                   &
                    (wk1v(im,j)-wk1v(1,j) + q(im,j)*(crx(1,j)-crx(im,j)))
      enddo

      call ycc(im, jm, fy, adx, cry, ymass, jord, 0,jfirst,jlast)

! For Scalar only!!!
      if ( jfirst == 1 ) then   ! ( jfirst -ng_d <= 1 ) fails when 
                                ! ng_d=3, ng_c=2, jlast-jfirst+1 = 3
        do i=1,im2
          q(i,1) = q(i+im2,  2)
        enddo
        do i=im2+1,im
           q(i,1) = q(i-im2,  2)
        enddo
      endif

      if ( jlast == jm ) then
        do i=1,im2
          fx1(i) = q(i+im2,jm)
        enddo
        do i=im2+1,im
           fx1(i) = q(i-im2,jm)
        enddo

        do i=1,im
          if(va(i,jm) .gt. D0_0) then
            adx(i,jm) = q(i,jm) + D0_5*va(i,jm)*(q(i,jm-1)-q(i,jm))
          else
            adx(i,jm) = q(i,jm) + D0_5*va(i,jm)*(q(i,jm)-fx1(i))
          endif
        enddo
      endif

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i,jp)
#endif
      do j=js2g0,jn2g0
        do i=1,im
          jp = j-va(i,j)
! jp = j     if va < 0
! jp = j -1  if va < 0
! q needed max(1, jfirst-1)
          adx(i,j) = q(i,j) + D0_5*va(i,j)*(q(i,jp)-q(i,jp+1))
        enddo
      enddo

      call xtpv( im, ffsl, fx, adx, crx, iord, crx,         &
               cose, 0, dm, qtmpv, al, ar, a6,              &
               jfirst, jlast, js2g0, jn1g0, jm,             &
               1, jm, jfirst, jlast,                        &
               jfirst-1, jlast+2,jfirst-ng_c, jlast+ng_c,   &
               jfirst-ng_c, jlast+ng_c, jfirst-1, jlast+2)

    return
!EOC
 end subroutine tpcc
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ycc
!
! !INTERFACE: 
 subroutine ycc(im, jm, fy, q, vc, ymass, jord, iv, jfirst, jlast)
!-----------------------------------------------------------------------

! !USES:
 implicit none

! !INPUT PARAMETERS:
 integer im, jm                      !  Dimensions
 integer jfirst, jlast               !  Latitude strip
 integer jord                        !  Approximation order
 integer iv                          !  Scalar=0, Vector=1
 real (r8) q(im,jfirst-1-iv:jlast+2)      !  Field (N*2 S*(iv+1))
 real (r8) vc(im,jfirst-iv:jlast)         !  Courant  (like FY)
 real (r8) ymass(im,jfirst-iv:jlast)      !  background mass flux

! !OUTPUT PARAMETERS:
 real (r8) fy(im,jfirst-iv:jlast)         !  Flux (S if iv=1)

! !DESCRIPTION:
!     Will Sawyer's note: In this routine the number 
!     of ghosted latitudes NG is min(abs(jord),2).  The scalar/vector
!     flag determines whether the flux FY needs to be ghosted on the
!     south.  If called from CD\_CORE (iv==1) then it does, if called
!     from TPCC (iv==0) it does not.  
!
! !CALLED FROM:
!     cd_core
!     tpcc
!
! !REVISION HISTORY:
!
!   SJL 99.04.13:  Delivery
!   WS  99.04.19:  Added jfirst:jlast concept
!   WS  99.04.27:  DC removed as argument (local to this routine); DC on N
!   WS  99.05.10:  Replaced JNP with JM, JMR with JM-1, IMR with IM
!   WS  99.05.10:  Removed fvcore.h
!   WS  99.07.27:  Built in tests for SP or NP
!   WS  99.09.09:  Documentation; indentation; cleaning; pole treatment
!   WS  99.09.14:  Loop limits
!   WS  00.05.14:  Renamed ghost indices as per Kevin's definitions
!
!EOP
!---------------------------------------------------------------------
!BOC

! !LOCAL VARIABLES:
  real (r8) dc(im,jfirst-iv:jlast+1)
  real (r8) qmax, qmin
  integer i, j, jt, im2, js2giv, js3giv, jn2g1, jn2g0


   im2 = im/2

   js2giv = max(2,jfirst-iv)
   js3giv = max(3,jfirst-iv)
   jn2g1  = min(jm-1,jlast+1)
   jn2g0  = min(jm-1,jlast)
      
   if(jord == 1) then
#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i,jt)
#endif
        do j=js2giv,jn2g0                      ! FY needed on S*iv
          do i=1,im
! jt=j if vc > 0; jt=j+1 if vc <=0
            jt = real(j+1,r8)  - vc(i,j)         ! VC ghosted like fy
            fy(i,j) = q(i,jt)*ymass(i,j)       ! ymass ghosted like fy
          enddo                                ! q ghosted N*1, S*iv
        enddo

   else

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
        do j=js3giv,jn2g1                      ! dc needed N*1, S*iv
          do i=1,im
            dc(i,j) = D0_25*(q(i,j+1)-q(i,j-1)) ! q ghosted N*2, S*(iv+1)
          enddo
        enddo

        if(iv.eq.0) then
! Scalar.

! WS 99.07.27 : Split loops in SP and NP regions, added SP/NP tests

          if ( jfirst-iv <= 2 ) then
            do i=1,im2
              dc(i, 2) = D0_25 * ( q(i,3) - q(i+im2,2) )
            enddo

            do i=im2+1,im
              dc(i, 2) = D0_25 * ( q(i,3) - q(i-im2,2) )
            enddo
          endif

          if ( jlast == jm ) then
            do i=1,im2
              dc(i,jm) = D0_25 * ( q(i+im2,jm) - q(i,jm-1) )
            enddo

            do i=im2+1,im
              dc(i,jm) = D0_25 * ( q(i-im2,jm) - q(i,jm-1) )
            enddo
          endif

        else
! Vector winds

! WS 99.07.27 : Split loops in SP and NP regions, added SP/NP tests

          if ( jfirst-iv <= 2 ) then
            do i=1,im2
              dc(i, 2) =  D0_25 * ( q(i,3) + q(i+im2,2) )
            enddo

            do i=im2+1,im
              dc(i, 2) =  D0_25 * ( q(i,3) + q(i-im2,2) )
            enddo
          endif

          if ( jlast == jm ) then
            do i=1,im2
              dc(i,jm) = -D0_25 * ( q(i,jm-1) + q(i+im2,jm) )
            enddo

            do i=im2+1,im
              dc(i,jm) = -D0_25 * ( q(i,jm-1) + q(i-im2,jm) )
            enddo
          endif

        endif

        if( jord > 0 ) then
! Monotonic constraint
#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i,qmax,qmin)
#endif
          do j=js3giv,jn2g1            ! DC needed N*1, S*iv
            do i=1,im                  ! P ghosted N*2, S*(iv+1)
              qmax = max(q(i,j-1),q(i,j),q(i,j+1)) - q(i,j)
              qmin = q(i,j) - min(q(i,j-1),q(i,j),q(i,j+1))
              dc(i,j) = sign(min(abs(dc(i,j)),qmin,qmax),dc(i,j))
            enddo
          enddo
!
! WS 99.08.03 : Following loop split into SP and NP part
!
          if ( jfirst-iv <= 2 ) then
            do i=1,im
              dc(i, 2) = D0_0
            enddo
          endif
          if ( jlast == jm ) then
            do i=1,im
              dc(i,jm) = D0_0
            enddo
          endif
        endif

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i,jt)
#endif
       do j=js2giv,jn2g0                   ! fy needed S*iv
         do i=1,im                       
           jt = real(j+1,r8)  - vc(i,j)      ! vc, ymass ghosted like fy
           fy(i,j) = (q(i,jt)+(sign(D1_0,vc(i,j))-vc(i,j))*dc(i,jt))*ymass(i,j)
         enddo
       enddo
    endif
    return
!EOC
 end subroutine ycc
!-----------------------------------------------------------------------

#ifdef VECTORIZE
!-----------------------------------------------------------------------
!BOP
! !IROUTINE: xtpv
!
! !INTERFACE: 
 subroutine xtpv(im, ffslv,  fxv,  qv,  cv,  iord,  mfxv,        &
                cosav, id, dm, qtmpv, al, ar, a6,                &
                jfirst, jlast, jlow, jhigh, jm,                  &
                jl2, jh2, jl3, jh3,                              &
                jl4, jh4, jl5, jh5,                              &
                jl7, jh7, jl11, jh11)
!-----------------------------------------------------------------------

 implicit none
 
! !INPUT PARAMETERS:
   integer id               ! ID = 0: density (mfx = C)
                            ! ID = 1: mixing ratio (mfx is mass flux)

   integer im               ! Total longitudes
   real (r8) cv(im,jl5:jh5)          ! Courant numbers
   real (r8) qv(im,jl4:jh4)
   real (r8) mfxv(im,jl7:jh7)
   logical ffslv(jl2:jh2)
   integer iord
   integer jfirst, jlast, jlow, jhigh, jm
   integer jl2, jh2, jl3, jh3, jl4, jh4, jl5, jh5
   integer jl7, jh7, jl11, jh11 
   real (r8) cosav(jm)

! !INPUT/OUTPUT PARAMETERS:
   real (r8) qtmpv(-im/3:im+im/3,jl11:jh11)   ! Input work arrays:
   real (r8)   dm(-im/3:im+im/3)
   real (r8)   al(-im/3:im+im/3)
   real (r8)   ar(-im/3:im+im/3)
   real (r8)   a6(-im/3:im+im/3)

! !OUTPUT PARAMETERS:
   real (r8) fxv(im,jl3:jh3)

! !DESCRIPTION:
!   
!
! !REVISION HISTORY:
!   99.01.01 Lin      Creation
!   01.03.27 Sawyer   Additional ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC

! Local:
   real (r8)       cos_upw               !critical cosine for upwind
   real (r8)       cos_van               !critical cosine for van Leer
   real (r8)       cos_ppm               !critical cosine for ppm

   parameter (cos_upw = D0_05)       !roughly at 87 deg.
   parameter (cos_van = D0_25)       !roughly at 75 deg.
   parameter (cos_ppm = D0_25)

   real (r8) r24
   parameter (r24 = D1_0/D24_0)

   integer i, imp, j
   real (r8) qmax, qmin
   real (r8) rut, tmp
   real (r8) dmv(-im/3:im+im/3,jlow:jhigh)
   integer iu, itmp, ist
   integer isave(im,jlow:jhigh)
   integer iuwv(jlow:jhigh), iuev(jlow:jhigh)

   integer jatn, jafn, ja
   integer jat(jhigh-jlow+1), jaf(jhigh-jlow+1)
   integer jattn, jatfn, jaftn, jaffn
   integer jatt(jhigh-jlow+1), jatf(jhigh-jlow+1)
   integer jaft(jhigh-jlow+1), jaff(jhigh-jlow+1)
   integer jatftn, jatffn
   integer jatft(jhigh-jlow+1), jatff(jhigh-jlow+1)
   integer jafftn1, jafffn1
   integer jafft1(jhigh-jlow+1), jafff1(jhigh-jlow+1)
   integer jafftn2, jafffn2
   integer jafft2(jhigh-jlow+1), jafff2(jhigh-jlow+1)
   real (r8) qsum((-im/3)-1:im+im/3,jlow:jhigh)   ! work arrays


   jatn = 0
   jafn = 0
   jattn = 0
   jatfn = 0
   jaftn = 0
   jaffn = 0
   jatftn = 0
   jatffn = 0
   jafftn1 = 0
   jafffn1 = 0
   jafftn2 = 0
   jafffn2 = 0
!call ftrace_region_begin("xtpv_1")
   do j = jlow, jhigh
     if (ffslv(j)) then
        jatn = jatn + 1
        jat(jatn) = j
        if( iord == 1 .or. cosav(j) < cos_upw) then
          jattn = jattn + 1
          jatt(jattn) = j
        else
          jatfn = jatfn + 1
          jatf(jatfn) = j
          if(iord .ge. 3 .and. cosav(j) .gt. cos_ppm) then
            jatftn = jatftn + 1
            jatft(jatftn) = j
          else
            jatffn = jatffn + 1
            jatff(jatffn) = j
          endif
        endif
     else
        jafn = jafn + 1
        jaf(jafn) = j
        if( iord == 1 .or. cosav(j) < cos_upw) then
          jaftn = jaftn + 1
          jaft(jaftn) = j
        else
          jaffn = jaffn + 1
          jaff(jaffn) = j
          if(iord > 0 .or. cosav(j) < cos_van) then
            jafftn1 = jafftn1 + 1
            jafft1(jafftn1) = j
          else
            jafffn1 = jafffn1 + 1
            jafff1(jafffn1) = j
          endif
          if( abs(iord).eq.2 .or. cosav(j) .lt. cos_van ) then
            jafftn2 = jafftn2 + 1
            jafft2(jafftn2) = j
          else
            jafffn2 = jafffn2 + 1
            jafff2(jafffn2) = j
          endif
        endif
     endif
   enddo
!call ftrace_region_end("xtpv_1")

   imp = im + 1

   do j = jlow, jhigh
     do i=1,im
       qtmpv(i,j) = qv(i,j)
     enddo
   enddo

! Flux-Form Semi-Lagrangian transport

!call ftrace_region_begin("xtpv_2")
!dir$ concurrent
   do ja = 1, jatn
      j = jat(ja)

! Figure out ghost zone for the western edge:
      iuwv(j) =  -cv(1,j)
      iuwv(j) = min(0, iuwv(j))
 
      do i=iuwv(j), 0
         qtmpv(i,j) = qv(im+i,j)
      enddo 

! Figure out ghost zone for the eastern edge:
      iuev(j) = im - cv(im,j)
      iuev(j) = max(imp, iuev(j))

      do i=imp, iuev(j)
         qtmpv(i,j) = qv(i-im,j)
      enddo

   enddo
!call ftrace_region_end("xtpv_2")

!dir$ concurrent
!call ftrace_region_begin("xtpv_3")
   do ja = 1, jattn
      j = jatt(ja)

      do i=1,im
        iu = cv(i,j)
        if(cv(i,j) .le. D0_0) then
          itmp = i - iu
          isave(i,j) = itmp - 1
        else
          itmp = i - iu - 1
          isave(i,j) = itmp + 1
        endif
        fxv(i,j) = (cv(i,j)-iu) * qtmpv(itmp,j)
      enddo

   enddo
!call ftrace_region_end("xtpv_3")

!dir$ concurrent
!call ftrace_region_begin("xtpv_4")
   do ja = 1, jatfn
      j = jatf(ja)

      do i=1,im
! 2nd order slope
        tmp = D0_25*(qtmpv(i+1,j) - qtmpv(i-1,j))
        qmax = max(qtmpv(i-1,j), qtmpv(i,j), qtmpv(i+1,j)) - qtmpv(i,j)
        qmin = qtmpv(i,j) - min(qtmpv(i-1,j), qtmpv(i,j), qtmpv(i+1,j))
        dmv(i,j) = sign(min(abs(tmp),qmax,qmin), tmp)
      enddo

      do i=iuwv(j), 0
        dmv(i,j) = dmv(im+i,j)
      enddo 

      do i=imp, iuev(j)
        dmv(i,j) = dmv(i-im,j)
      enddo

   enddo
!call ftrace_region_end("xtpv_4")

   call fxppmv(im, cv, mfxv, qtmpv, dmv, fxv, iord,                     &
              iuwv, iuev, ffslv, isave, jatftn, jatft, jlow, jhigh,     &
              jl2, jh2, jl3, jh3, jl5, jh5, jl7, jh7, jl11, jh11)

!dir$ concurrent
!call ftrace_region_begin("xtpv_5")
   do ja = 1, jatffn
      j = jatff(ja)

      do i=1,im
        iu  = cv(i,j)
        rut = cv(i,j) - iu
        if(cv(i,j) .le. D0_0) then
          itmp = i - iu
          isave(i,j) = itmp - 1
          fxv(i,j) = rut*(qtmpv(itmp,j)-dmv(itmp,j)*(D1_0+rut))
        else
          itmp = i - iu - 1
          isave(i,j) = itmp + 1
          fxv(i,j) = rut*(qtmpv(itmp,j)+dmv(itmp,j)*(D1_0-rut))
        endif
      enddo

   enddo
!call ftrace_region_end("xtpv_5")

!dir$ concurrent
!call ftrace_region_begin("xtpv_6")
   do ja = 1, jatn
      j = jat(ja)
      qsum(iuwv(j)-1,j) = D0_0
      do i = iuwv(j), iuev(j)
        qsum(i,j) = qsum(i-1,j) + qtmpv(i,j)
      end do

!
! The boolean terms:
! a)    .and. (isave(i,j) < i)
! b)    .and. (i <= isave(i,j))
! are needed in the IF statements below because I cannot prove to myself
! that the relationship between i and isave are such to guarantee that
! there is always at least one term from qsum (qtmpv,j) contributed to fxv.
!

      do i=1,im
        if(cv(i,j) >= D1_0 .and. (isave(i,j) < i) ) then
            fxv(i,j) = fxv(i,j) + (qsum(i-1,j) - qsum(isave(i,j) - 1,j))
        else if (cv(i,j) <= -D1_0 .and. (i <= isave(i,j)) ) then
            fxv(i,j) = fxv(i,j) - (qsum(isave(i,j),j) - qsum(i-1,j))
        end if
      end do

      if(id .ne. 0) then
         do i=1,im
            fxv(i,j) =  fxv(i,j)*mfxv(i,j)
         enddo
      endif

   enddo
!call ftrace_region_end("xtpv_6")

! Regular PPM (Eulerian without FFSL extension)

!call ftrace_region_begin("xtpv_7")
!dir$ concurrent
!cdir nodep
   do ja = 1, jafn
      j = jaf(ja)

      qtmpv(imp,j) = qv(1,j)
      qtmpv(  0,j) = qv(im,j)
   enddo

!dir$ concurrent
   do ja = 1, jaftn
      j = jaft(ja)

      do i=1,im
         iu = real(i,r8) - cv(i,j)
         fxv(i,j) = mfxv(i,j)*qtmpv(iu,j)
      enddo
   enddo

!dir$ concurrent
!cdir nodep
   do ja = 1, jaffn
      j = jaff(ja)

      qtmpv(-1,j)    = qv(im-1,j)
      qtmpv(imp+1,j) = qv(2,j)

   enddo
!call ftrace_region_end("xtpv_7")

!dir$ concurrent
!call ftrace_region_begin("xtpv_8")
   do ja = 1, jafftn1
      j = jafft1(ja)

! In-line xmist

      do i=1,im
         dmv(i,j) = r24*(D8_0*(qtmpv(i+1,j) - qtmpv(i-1,j)) + qtmpv(i-2,j) - qtmpv(i+2,j))
      enddo

! Apply monotonicity constraint (Lin et al. 1994, MWR)
      do i=1,im
         qmax = max( qtmpv(i-1,j), qtmpv(i,j), qtmpv(i+1,j) ) - qtmpv(i,j)
         qmin = qtmpv(i,j) - min( qtmpv(i-1,j), qtmpv(i,j), qtmpv(i+1,j) )
         dmv(i,j) = sign( min(abs(dmv(i,j)), qmax, qmin), dmv(i,j) )
      enddo

   enddo
!call ftrace_region_end("xtpv_8")

!dir$ concurrent
!call ftrace_region_begin("xtpv_9")
   do ja = 1, jafffn1
      j = jafff1(ja)

! In-line xmist

      if(iord .le. 2) then
         do i=1,im
            dmv(i,j) = r24*(D8_0*(qtmpv(i+1,j) - qtmpv(i-1,j)) + qtmpv(i-2,j) - qtmpv(i+2,j))
         enddo
      else
         do i=1,im
            dmv(i,j) = D0_25*(qtmpv(i+1,j) - qtmpv(i-1,j))
         enddo
      endif

      if( iord >= 0 ) then

! Apply monotonicity constraint (Lin et al. 1994, MWR)
         do i=1,im
            qmax = max( qtmpv(i-1,j), qtmpv(i,j), qtmpv(i+1,j) ) - qtmpv(i,j)
            qmin = qtmpv(i,j) - min( qtmpv(i-1,j), qtmpv(i,j), qtmpv(i+1,j) )
            dmv(i,j) = sign( min(abs(dmv(i,j)), qmax, qmin), dmv(i,j) )
         enddo
      endif

   enddo
!call ftrace_region_end("xtpv_9")

!call ftrace_region_begin("xtpv_10")
!dir$ concurrent
!cdir nodep
   do ja = 1, jaffn
      j = jaff(ja)

      dmv(0,j) = dmv(im,j)

   enddo
!call ftrace_region_end("xtpv_10")

!dir$ concurrent
!call ftrace_region_begin("xtpv_11")
   do ja = 1, jafftn2
      j = jafft2(ja)

      do i=1,im
         iu = real(i,r8) - cv(i,j)
         fxv(i,j) =  mfxv(i,j)*(qtmpv(iu,j)+dmv(iu,j)*(sign(D1_0,cv(i,j))-cv(i,j)))
      enddo

   enddo
!call ftrace_region_end("xtpv_11")

   call fxppmv(im, cv, mfxv, qtmpv, dmv, fxv, iord,                     &
              iuwv, iuev, ffslv, isave, jafffn2, jafff2, jlow, jhigh,   &
              jl2, jh2, jl3, jh3, jl5, jh5, jl7, jh7, jl11, jh11)

   return
!EOC
 end subroutine xtpv
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: fxppmv
!
! !INTERFACE: 

 subroutine fxppmv(im, c, mfx, p, dm, fx, iord,                         &
                  iuw, iue, ffsl, isave, jan, ja, jlow, jhigh,          &
                  jl2, jh2, jl3, jh3, jl5, jh5, jl7, jh7, jl11, jh11)
!-----------------------------------------------------------------------
!
! !USES:
 implicit none

! !INPUT PARAMETERS:
 integer jan, ja(jan), jlow, jhigh, jj, j
 integer jl2, jh2, jl3, jh3, jl5, jh5, jl7, jh7, jl11, jh11
 integer im, iord
 real (r8)  c(im,jl5:jh5)
 real (r8) p(-im/3:im+im/3,jl11:jh11)
 real (r8) dm(-im/3:im+im/3,jlow:jhigh)
 real (r8) mfx(im,jl7:jh7)
 integer iuw(jlow:jhigh), iue(jlow:jhigh)
 logical ffsl(jl2:jh2)
 integer isave(im,jlow:jhigh)

! !OUTPUT PARAMETERS:
 real (r8) fx(im,jl3:jh3)

! !DESCRIPTION:
!   
!
! !REVISION HISTORY:
!   99.01.01 Lin      Creation
!   01.03.27 Sawyer   Additional ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
 real (r8) r3, r23
 parameter ( r3 = D1_0/D3_0, r23 = D2_0/D3_0 )

 integer i, lmt
 integer iu, itmp
 real (r8) ru
 logical steep
 real (r8) al(-im/3:im+im/3,jlow:jhigh)
 real (r8) ar(-im/3:im+im/3,jlow:jhigh)
 real (r8) a6(-im/3:im+im/3,jlow:jhigh)

 integer jbtn, jbfn
 integer jbt(jan), jbf(jan)
 integer ilow, ihigh

 ilow = -im/3
 ihigh = im + im/3

 if( iord == 6 ) then
   steep = .true.
 else
   steep = .false.
 endif

!dir$ concurrent
 do jj = 1, jan
   j = ja(jj)

   do i=1,im
     al(i,j) = D0_5*(p(i-1,j)+p(i,j)) + (dm(i-1,j) - dm(i,j))*r3
   enddo

 enddo

 if (steep) then

   call steepxv( im, p, al, dm, jan, ja, jlow, jhigh, jl11, jh11 )

 endif

!dir$ concurrent
 do jj = 1, jan
   j = ja(jj)

   do i=1,im-1
     ar(i,j) = al(i+1,j)
   enddo
   ar(im,j) = al(1,j)

 enddo

 if(iord == 7) then

   call huynhv(im, ar, al, p, a6, dm, jan, ja, jlow, jhigh, jl11, jh11 )

 else

   if(iord .eq. 3 .or. iord .eq. 5) then

!dir$ concurrent
     do jj = 1, jan
       j = ja(jj)

       do i=1,im
         a6(i,j) = D3_0*(p(i,j)+p(i,j)  - (al(i,j)+ar(i,j)))
       enddo

     enddo
   endif

   lmt = iord - 3

   call lmppmv( dm, a6, ar, al, p, im, lmt, jan, ja, ilow, ihigh,    &
                jlow, jhigh, jl11, jh11 )

 endif

 jbtn = 0
 jbfn = 0
!dir$ concurrent
 do jj = 1, jan
   j = ja(jj)
   if( ffsl(j) ) then
     jbtn = jbtn + 1
     jbt(jbtn) = j
   else
     jbfn = jbfn + 1
     jbf(jbfn) = j
   endif
 enddo

!dir$ concurrent
 do jj = 1, jbtn
   j = jbt(jj)

   do i=iuw(j), 0
     al(i,j) = al(im+i,j)
     ar(i,j) = ar(im+i,j)
     a6(i,j) = a6(im+i,j)
   enddo

   do i=im+1, iue(j)
     al(i,j) = al(i-im,j)
     ar(i,j) = ar(i-im,j)
     a6(i,j) = a6(i-im,j)
   enddo

   do i=1,im
     iu = c(i,j)
     ru = c(i,j) - iu
     if(c(i,j) .gt. D0_0) then
       itmp = i - iu - 1
       isave(i,j) = itmp + 1
       fx(i,j) = ru*(ar(itmp,j)+D0_5*ru*(al(itmp,j)-ar(itmp,j) +     &
                 a6(itmp,j)*(D1_0-r23*ru)) )
     else
       itmp = i - iu
       isave(i,j) = itmp - 1
       fx(i,j) = ru*(al(itmp,j)-D0_5*ru*(ar(itmp,j)-al(itmp,j) +     &
                 a6(itmp,j)*(D1_0+r23*ru)) )
     endif
   enddo

 enddo

!dir$ concurrent
 do jj = 1, jbfn
   j = jbf(jj)

   al(0,j) = al(im,j)
   ar(0,j) = ar(im,j)
   a6(0,j) = a6(im,j)
   do i=1,im
     if(c(i,j) .gt. D0_0) then
       fx(i,j) = ar(i-1,j) + D0_5*c(i,j)*(al(i-1,j) - ar(i-1,j) +   &
                 a6(i-1,j)*(D1_0-r23*c(i,j)) )
     else
       fx(i,j) = al(i,j) - D0_5*c(i,j)*(ar(i,j) - al(i,j) +         &
                 a6(i,j)*(D1_0+r23*c(i,j)))
     endif
     fx(i,j) = mfx(i,j) * fx(i,j)
   enddo

 enddo

 return
!EOC
 end subroutine fxppmv
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: steepxv
!
! !INTERFACE: 
 subroutine steepxv(im, p, al, dm, jan, ja, jlow, jhigh, jl11, jh11 )
!-----------------------------------------------------------------------

! !USES:
 implicit none

! !INPUT PARAMETERS:
 integer im
 integer jan, ja(jan), jlow, jhigh, jl11, jh11
 real (r8)  p(-im/3:im+im/3,jl11:jh11)
 real (r8) dm(-im/3:im+im/3,jlow:jhigh)

! !INPUT/OUTPUT PARAMETERS:
 real (r8) al(im,jlow:jhigh)

! !DESCRIPTION:
!   
!
! !REVISION HISTORY:
!   99.01.01 Lin      Creation
!   01.03.27 Sawyer   Additional ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
 integer i, jj, j
 real (r8) r3
 parameter ( r3 = D1_0/D3_0 )

 real (r8) dh(0:im,jlow:jhigh)
 real (r8) d2(0:im+1,jlow:jhigh)
 real (r8) eta(0:im,jlow:jhigh)
 real (r8) xxx, bbb, ccc

!dir$ concurrent
 do jj = 1, jan
   j = ja(jj)

   do i=0,im
      dh(i,j) = p(i+1,j) - p(i,j)
   enddo

! Needs dh(0:im,j)
   do i=1,im
      d2(i,j) = dh(i,j) - dh(i-1,j)
   enddo
   d2(0,j) = d2(im,j)
   d2(im+1,j) = d2(1,j)

! needs p(-1:im+2,j), d2(0:im+1,j)
   do i=1,im
      if( d2(i+1,j)*d2(i-1,j).lt.D0_0 .and. p(i+1,j).ne.p(i-1,j) ) then
          xxx    = D1_0 - D0_5 * ( p(i+2,j) - p(i-2,j) ) / ( p(i+1,j) - p(i-1,j) )
          eta(i,j) = max(D0_0, min(xxx, D0_5) )
      else
          eta(i,j) = D0_0
      endif
   enddo

   eta(0,j) = eta(im,j)

! needs eta(0:im,j), dh(0:im-1,j), dm(0:im,j)
   do i=1,im
      bbb = ( D2_0*eta(i,j  ) - eta(i-1,j) ) * dm(i-1,j) 
      ccc = ( D2_0*eta(i-1,j) - eta(i,j  ) ) * dm(i,j  ) 
      al(i,j) = al(i,j) + D0_5*( eta(i-1,j) - eta(i,j)) * dh(i-1,j) + (bbb - ccc) * r3
   enddo

 enddo

 return
!EOC
 end subroutine steepxv
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: huynhv --- Enforce Huynh's 2nd constraint in 1D periodic domain
!
! !INTERFACE: 
 subroutine huynhv(im, ar, al, p, d2, d1, jan, ja, jlow, jhigh, jl11, jh11)
!-----------------------------------------------------------------------

! !USES:

 implicit none

! !INPUT PARAMETERS:
 integer im
 integer jan, ja(jan), jlow, jhigh, jl11, jh11
 real(r8)  p(im,jl11:jh11)

! !OUTPUT PARAMETERS:
 real(r8) ar(im,jlow:jhigh)
 real(r8) al(im,jlow:jhigh)
 real(r8) d2(im,jlow:jhigh)
 real(r8) d1(im,jlow:jhigh)

! !DESCRIPTION:
!
!   Enforce Huynh's 2nd constraint in 1D periodic domain
!
! !REVISION HISTORY:
!   99.01.01 Lin      Creation
!   01.03.27 Sawyer   Additional ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
 integer  i, jj, j
 real(r8) pmp
 real(r8) lac
 real(r8) pmin
 real(r8) pmax

!dir$ concurrent
   do jj = 1, jan
      j = ja(jj)

! Compute d1 and d2
      d1(1,j) = p(1,j) - p(im,j)
      do i=2,im
         d1(i,j) = p(i,j) - p(i-1,j)
      enddo

      do i=1,im-1
         d2(i,j) = d1(i+1,j) - d1(i,j)
      enddo
      d2(im,j) = d1(1,j) - d1(im,j)

! Constraint for AR
!     i = 1
      pmp   = p(1,j) + D2_0 * d1(1,j)
      lac   = p(1,j) + D0_5 * (d1(1,j)+d2(im,j)) + d2(im,j) 
      pmin  = min(p(1,j), pmp, lac)
      pmax  = max(p(1,j), pmp, lac)
      ar(1,j) = min(pmax, max(ar(1,j), pmin))

      do i=2, im
         pmp   = p(i,j) + D2_0*d1(i,j)
         lac   = p(i,j) + D0_5*(d1(i,j)+d2(i-1,j)) + d2(i-1,j)
         pmin  = min(p(i,j), pmp, lac)
         pmax  = max(p(i,j), pmp, lac)
         ar(i,j) = min(pmax, max(ar(i,j), pmin))
      enddo

! Constraint for AL
      do i=1, im-1
         pmp   = p(i,j) - D2_0*d1(i+1,j)
         lac   = p(i,j) + D0_5*(d2(i+1,j)-d1(i+1,j)) + d2(i+1,j)
         pmin  = min(p(i,j), pmp, lac)
         pmax  = max(p(i,j), pmp, lac)
         al(i,j) = min(pmax, max(al(i,j), pmin))
      enddo

! i=im
      i = im
      pmp    = p(im,j) - D2_0*d1(1,j)
      lac    = p(im,j) + D0_5*(d2(1,j)-d1(1,j)) + d2(1,j)
      pmin   = min(p(im,j), pmp, lac)
      pmax   = max(p(im,j), pmp, lac)
      al(im,j) = min(pmax, max(al(im,j), pmin))

! compute A6 (d2)
      do i=1, im
         d2(i,j) = D3_0*(p(i,j)+p(i,j)  - (al(i,j)+ar(i,j)))
      enddo

   enddo

   return
!EOC
 end subroutine huynhv
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: lmppmv
!
! !INTERFACE: 
 subroutine lmppmv(dm, a6, ar, al, p, im, lmt, jan, ja,           &
                  ilow, ihigh, jlow, jhigh, jl11, jh11)
!-----------------------------------------------------------------------

 implicit none

! !INPUT PARAMETERS:
 integer im   ! Total longitudes
 integer jan, ja(jan), ilow, ihigh, jlow, jhigh, jl11, jh11
 integer lmt  ! LMT = 0: full monotonicity
              ! LMT = 1: Improved and simplified full monotonic constraint
              ! LMT = 2: positive-definite constraint
              ! LMT = 3: Quasi-monotone constraint
 real(r8) p(ilow:ihigh,jl11:jh11)
 real(r8) dm(ilow:ihigh,jlow:jhigh)

! !OUTPUT PARAMETERS:
 real(r8) a6(ilow:ihigh,jlow:jhigh)
 real(r8) ar(ilow:ihigh,jlow:jhigh)
 real(r8) al(ilow:ihigh,jlow:jhigh)

! !DESCRIPTION:
!
!
! !REVISION HISTORY:
!   99.01.01 Lin      Creation
!   01.03.27 Sawyer   Additional ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
 real (r8) r12
 parameter ( r12 = D1_0/D12_0 )

 real (r8) da1, da2, fmin, a6da
 real (r8) dr, dl

 integer i, jj, j

! LMT = 0: full monotonicity
! LMT = 1: Improved and simplified full monotonic constraint
! LMT = 2: positive-definite constraint
! LMT = 3: Quasi-monotone constraint

  if( lmt == 0 ) then

! Full constraint

!dir$ concurrent
    do jj = 1, jan
      j = ja(jj)

        do i=1,im
           if(dm(i,j) .eq. D0_0) then
               ar(i,j) = p(i,j)
               al(i,j) = p(i,j)
               a6(i,j) = D0_0
           else
               da1  = ar(i,j) - al(i,j)
               da2  = da1**2
               a6da = a6(i,j)*da1
               if(a6da .lt. -da2) then
                  a6(i,j) = D3_0*(al(i,j)-p(i,j))
                  ar(i,j) = al(i,j) - a6(i,j)
               elseif(a6da .gt. da2) then
                  a6(i,j) = D3_0*(ar(i,j)-p(i,j))
                  al(i,j) = ar(i,j) - a6(i,j)
               endif
           endif
        enddo

    enddo

  elseif( lmt == 1 ) then

! Improved (Lin 2001?) full constraint

!dir$ concurrent
    do jj = 1, jan
      j = ja(jj)

        do i=1,im
            da1 = dm(i,j) + dm(i,j)
            dl = sign(min(abs(da1),abs(al(i,j)-p(i,j))), da1)
            dr = sign(min(abs(da1),abs(ar(i,j)-p(i,j))), da1)
            ar(i,j) = p(i,j) + dr
            al(i,j) = p(i,j) - dl
            a6(i,j) = D3_0*(dl-dr)
        enddo

    enddo

  elseif( lmt == 2 ) then

! Positive definite constraint

!dir$ concurrent
    do jj = 1, jan
      j = ja(jj)

      do i=1,im
        if(abs(ar(i,j)-al(i,j)) .lt. -a6(i,j)) then
          fmin = p(i,j) + D0_25*(ar(i,j)-al(i,j))**2/a6(i,j) + a6(i,j)*r12
          if(fmin.lt.D0_0) then
            if(p(i,j).lt.ar(i,j) .and. p(i,j).lt.al(i,j)) then
                ar(i,j) = p(i,j)
                al(i,j) = p(i,j)
                a6(i,j) = D0_0
            elseif(ar(i,j) .gt. al(i,j)) then
                a6(i,j) = D3_0*(al(i,j)-p(i,j))
                ar(i,j) = al(i,j) - a6(i,j)
            else
                a6(i,j) = D3_0*(ar(i,j)-p(i,j))
                al(i,j) = ar(i,j) - a6(i,j)
            endif
          endif
        endif
      enddo

    enddo

  elseif(lmt .eq. 3) then

! Quasi-monotone constraint

!dir$ concurrent
    do jj = 1, jan
      j = ja(jj)

      do i=1,im
         da1 = D4_0*dm(i,j)
         dl = sign(min(abs(da1),abs(al(i,j)-p(i,j))), da1)
         dr = sign(min(abs(da1),abs(ar(i,j)-p(i,j))), da1)
         ar(i,j) = p(i,j) + dr
         al(i,j) = p(i,j) - dl
         a6(i,j) = D3_0*(dl-dr)
      enddo

    enddo

  endif
  return
!EOC
 end subroutine lmppmv
!-----------------------------------------------------------------------
#endif

end module tp_core
