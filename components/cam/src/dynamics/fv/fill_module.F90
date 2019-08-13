module fill_module
!-----------------------------------------------------------------------
!  $Id$
!BOP
!
! !MODULE: fill_module --- utilities for filling in "bad" data

 use shr_kind_mod, only: r8 => shr_kind_r8
 use cam_logfile,  only: iulog

#ifdef NO_R16
   integer,parameter :: r16= selected_real_kind(12) ! 8 byte real
#else
   integer,parameter :: r16= selected_real_kind(24) ! 16 byte real
#endif

!
! !PUBLIC MEMBER FUNCTIONS:
      public filew, fillxy, fillz, filns, pfix

!
! !DESCRIPTION:
!
!    This module provides the basic utilities to fill in regions
!    with bad "data", for example slightly negative values in fields
!    which must be positive, like mixing ratios.  Generally this 
!    means borrowing positive values from neighboring cells.
!
! !REVISION HISTORY:
!   99.03.01   Lin        Creation
!   01.02.14   Lin        Routines coalesced into this module
!   01.03.26   Sawyer     Added ProTeX documentation
!   05.05.25   Sawyer     Merged CAM and GEOS5 versions
!
!EOP
!-----------------------------------------------------------------------

private
real(r8), parameter ::  D0_0                    =  0.0_r8
real(r8), parameter ::  D0_5                    =  0.5_r8
real(r8), parameter ::  D1_0                    =  1.0_r8
real(r8), parameter ::  D1_5                    =  1.5_r8

contains

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: filew --- Fill from east and west neighbors; essentially
!                      performing local flux adjustment
!
! !INTERFACE: 
 subroutine filew(q, im, jm, jfirst, jlast, acap, ipx, tiny, cosp2)

! !USES:

 implicit none

! !INPUT PARAMETERS:
 integer im                  ! Longitudes
 integer jm                  ! Total latitudes
 integer jfirst              ! Starting latitude
 integer jlast               ! Finishing latitude

 real(r8) tiny               ! A small number to pump up value
 real(r8) acap               ! 1/(polar cap area)
 real(r8) cosp2              ! cosine(lat) at j=2

! !INPUT/OUTPUT PARAMETERS:
 real(r8) q(im,jfirst:jlast) ! Field to adjust

! !OUTPUT PARAMETERS:
 integer ipx                 ! Flag:  0 if Q not change, 1 if changed

! !DESCRIPTION:
!   Check for "bad" data and fill from east and west neighbors
!
! !REVISION HISTORY:
!   01.99.10   Lin        Creation
!   01.07.30   Lin        Improvement
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
 real(r8) d0, d1, d2
 real(r8) qtmp(jfirst:jlast,im)
 real(r8) tinyl   ! local tiny mixing ratio
 real(r8) qmin

 integer i, j, jm1, ip2
 integer j1, j2
 integer imin, jmin

 j1 = max( jfirst,   2 )
 j2 = min( jlast, jm-1 )
 jm1 = jm-1
 ipx = 0

! Copy & swap direction for vectorization.
  do j=j1,j2
     do i=1,im
        qtmp(j,i) = q(i,j)
     enddo
  enddo
 
  do i=2,im-1
     do j=j1,j2
        if(qtmp(j,i) < D0_0) then
           tinyl = max(D0_0,qtmp(j,i-1),qtmp(j,i+1))*tiny
           ipx =  1
! west
           d0 = max(D0_0,qtmp(j,i-1))
           d1 = min(-qtmp(j,i),d0)
           qtmp(j,i-1) = qtmp(j,i-1) - d1
           qtmp(j,i) = qtmp(j,i) + d1
! east
           d0 = max(D0_0,qtmp(j,i+1))
           d2 = min(-qtmp(j,i),d0)
           qtmp(j,i+1) = qtmp(j,i+1) - d2
           qtmp(j,i) = qtmp(j,i) + d2 + tinyl
        endif
    enddo
  enddo
 
     i=1
  do j=j1,j2
     if(qtmp(j,i) < D0_0) then
        ipx =  1
        tinyl = max(D0_0,qtmp(j,im),qtmp(j,i+1))*tiny
! west
        d0 = max(D0_0,qtmp(j,im))
        d1 = min(-qtmp(j,i),d0)
        qtmp(j,im) = qtmp(j,im) - d1
        qtmp(j,i) = qtmp(j,i) + d1
! east
        d0 = max(D0_0,qtmp(j,i+1))
        d2 = min(-qtmp(j,i),d0)
        qtmp(j,i+1) = qtmp(j,i+1) - d2
        qtmp(j,i) = qtmp(j,i) + d2 + tinyl
      endif
  enddo

     i=im
  do j=j1,j2
     if(qtmp(j,i) < D0_0) then
        ipx =  1
        tinyl = max(D0_0,qtmp(j,i-1),qtmp(j,1))*tiny
! west
        d0 = max(D0_0,qtmp(j,i-1))
        d1 = min(-qtmp(j,i),d0)
        qtmp(j,i-1) = qtmp(j,i-1) - d1
        qtmp(j,i) = qtmp(j,i) + d1
! east
        d0 = max(D0_0,qtmp(j,1))
        d2 = min(-qtmp(j,i),d0)
        qtmp(j,1) = qtmp(j,1) - d2
        qtmp(j,i) = qtmp(j,i) + d2 + tinyl
     endif
  enddo

 if(ipx .ne. 0) then

!-----------
! Final pass
!-----------
    do i=1,im-1
       do j=j1,j2
          if (qtmp(j,i) < D0_0 ) then
! Take mass from east (essentially adjusting fx(i+1,j))
              qtmp(j,i+1) = qtmp(j,i+1) + qtmp(j,i)
              qtmp(j,i) = D0_0
          endif
       enddo
    enddo

    do i=im,2,-1
       do j=j1,j2
          if (qtmp(j,i) < D0_0 ) then
! Take mass from west (essentially adjusting fx(i,j))
              qtmp(j,i-1) = qtmp(j,i-1) + qtmp(j,i)
              qtmp(j,i) = D0_0
          endif
       enddo
    enddo

    do j=j1,j2

       qmin = D0_0
       do i=1, im
          if (qtmp(j,i) < qmin) then
             qmin = qtmp(j,i)
             imin = i
             jmin = j
          endif
       enddo

       if ( qmin < D0_0 ) then
          write(iulog,*) ' filew failed, worst i, j, qtmp, q = ', imin, jmin, qtmp(jmin,imin), q(imin,jmin)
       end if

       do i=1,im
          q(i,j) = qtmp(j,i)
       enddo
    enddo

 endif
 
! Check Poles.

 if ( jfirst == 1 ) then
      if(q(1,1) < D0_0) then
         call pfix(q(1,2),q(1,1),im,ipx,acap,cosp2)
      else
!            Check j=2
             ip2 = 0
         do i=1,im
            if(q(i,2).lt.D0_0) then
               ip2 = 1
               go to 322
            endif
         enddo
322      continue
         if(ip2.ne.0) call pfix(q(1,2),q(1,1),im,ipx,acap,cosp2)
      endif
 endif
 
 if ( jlast == jm ) then
      if(q(1,jm) < D0_0) then
         call pfix(q(1,jm1),q(1,jm),im,ipx,acap,cosp2)
      else
!             Check j=jm1
              ip2 = 0
         do i=1,im
            if(q(i,jm1) < D0_0) then
               ip2 = 1
               go to 323
            endif
         enddo
323      continue
         if(ip2.ne.0) call pfix(q(1,jm1),q(1,jm),im,ipx,acap,cosp2)
      endif
 endif

!EOC
 end subroutine filew
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: fillxy --- Fill from east, west, north and south neighbors
!
! !INTERFACE: 
 subroutine fillxy(q, im, jm, jfirst, jlast, acap, cosp, acosp)

! !USES:

 implicit none

 integer im                  ! Longitudes
 integer jm                  ! Total latitudes
 integer jfirst              ! Starting latitude
 integer jlast               ! Finishing latitude

 real(r8) acap               ! ???
 real(r8) cosp(jm)           ! ???
 real(r8) acosp(jm)          ! ???
!
! !INPUT/OUTPUT PARAMETERS:
 real(r8) q(im,jfirst:jlast) ! Field to adjust

! !DESCRIPTION:
!   Check for "bad" data and fill from east and west neighbors
!
! !BUGS:
!   Currently this routine only performs the east-west fill algorithm.
!   This is because the N-S fill is very hard to do in a reproducible
!   fashion when the problem is decomposed by latitudes.
!
! !REVISION HISTORY:
!   99.03.01   Lin        Creation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
  integer ipx, ipy, j1, j2
  real(r8) tiny
  parameter( tiny = 1.e-20_r8 )

    call filew(q,im,jm,jfirst,jlast,acap,ipx,tiny,cosp(2))

! WS 99.08.03 : S.-J. can you clean up the j1, j2 stuff here?
   if(ipx.ne.0) then

      j1 = max( 2,    jfirst )
      j2 = min( jm-1, jlast )
!
! WS 99.08.03 : see comments in "BUGS" above
!!!      call filns(q,im,jm,j1,j2,cosp,acosp,ipy,tiny)

!     if(ipy .ne. 0) then
! do fill zonally
! xfx is problematic
!     call xfix(q,IM,JM,tiny,qt)
!     endif

   endif

!EOC
 end subroutine fillxy
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: fillz --- Fill from neighbors below and above
!
! !INTERFACE: 
 subroutine fillz(im, i1, i2, km, nq, q, dp)

! !USES:

 implicit none

! !INPUT PARAMETERS:
   integer, intent(in) :: im                ! No. of longitudes
   integer, intent(in) :: km                ! No. of levels
   integer, intent(in) :: i1                ! Starting longitude
   integer, intent(in) :: i2                ! Finishing longitude
   integer, intent(in) :: nq                ! Total number of tracers
   real(r8), intent(in) ::  dp(im,km)       ! pressure thickness

! !INPUT/OUTPUT PARAMETERS:
   real(r8), intent(inout) :: q(im,km,nq)   ! tracer mixing ratio

! !DESCRIPTION:
!   Check for "bad" data and fill from east and west neighbors
!
! !BUGS:
!   Currently this routine only performs the east-west fill algorithm.
!   This is because the N-S fill is very hard to do in a reproducible
!   fashion when the problem is decomposed by latitudes.
!
! !REVISION HISTORY:
!   00.04.01   Lin        Creation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
   integer i, k, ic
   real(r8) qup, qly, dup

   do ic=1,nq
! Top layer
      do i=i1,i2
         if( q(i,1,ic) < D0_0) then
             q(i,2,ic) = q(i,2,ic) + q(i,1,ic)*dp(i,1)/dp(i,2)
             q(i,1,ic) = D0_0
          endif
      enddo

! Interior
      do k=2,km-1
         do i=i1,i2
         if( q(i,k,ic) < D0_0 ) then
! Borrow from above
             qup =  q(i,k-1,ic)*dp(i,k-1)
             qly = -q(i,k  ,ic)*dp(i,k  )
             dup =  min( D0_5*qly, qup )        !borrow no more than 50%
             q(i,k-1,ic) = q(i,k-1,ic) - dup/dp(i,k-1) 
! Borrow from below: q(i,k,ic) is still negative at this stage
             q(i,k+1,ic) = q(i,k+1,ic) + (dup-qly)/dp(i,k+1) 
             q(i,k  ,ic) = D0_0
          endif
          enddo
      enddo
 
! Bottom layer
      k = km
      do i=i1,i2
         if( q(i,k,ic) < D0_0) then
! Borrow from above
             qup =  q(i,k-1,ic)*dp(i,k-1)
             qly = -q(i,k  ,ic)*dp(i,k  )
             dup =  min( qly, qup )
             q(i,k-1,ic) = q(i,k-1,ic) - dup/dp(i,k-1) 
             q(i,k,ic) = D0_0
          endif
      enddo
   enddo
!EOC
end subroutine fillz
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: filns --- Fill from north and south neighbors
!
! !INTERFACE: 
 subroutine filns(q,im,jm,j1,j2,cosp,acosp,ipy,tiny)

! !USES:

 implicit none

! !INPUT PARAMETERS:
 integer im                  ! Longitudes
 integer jm                  ! Total latitudes
 integer j1                  ! Starting latitude
 integer j2                  ! Finishing latitude


 real(r8) tiny               ! A small number to pump up value
 real(r8) cosp(*)            ! ???
 real(r8) acosp(*)           ! ???

! !INPUT/OUTPUT PARAMETERS:
 real(r8) q(im,*)            ! Field to adjust

! !OUTPUT PARAMETERS:
 integer  ipy                ! Flag: 0 if no fill-in, 1 if fill-in

! !DESCRIPTION:
!   Check for "bad" data and fill from north and south neighbors
!
! !BUGS:
!   Currently this routine can only be used performs when the
!   problem is *not* distributed in latitude (i.e. j1=1, j2=jm).
!   This is because the N-S fill is very hard to do in a reproducible
!   fashion when the problem is decomposed by latitudes.
!
! !REVISION HISTORY:
!   99.03.01   Lin        Creation
!   05.06.30   Sawyer     Removed SAVE attribute for cap1 (recalculated)
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
 integer  i, j
! This definition of PI as opposed to 4._r16*atan(1._r16) does not 
! appear to generate non-zero differences in GEOS5 checkpoint files
 real(R16),parameter :: pi     = 3.1415926535897932384626433832795028841971_R16
 real(r8) :: dp, cap1, dq, dn, ds, d0, d1, d2
 
    dp = pi/real(jm-1,r16)
    cap1 = im*(D1_0-cos((j1-D1_5)*dp))/dp
 
    ipy = 0
    do j=j1+1,j2-1
      do i=1,im
      if(q(i,j).lt.D0_0) then
         ipy =  1
         dq  = - q(i,j)*cosp(j)
! North
         dn = q(i,j+1)*cosp(j+1)
         d0 = max(D0_0,dn)
         d1 = min(dq,d0)
         q(i,j+1) = (dn - d1)*acosp(j+1)
         dq = dq - d1
! South
         ds = q(i,j-1)*cosp(j-1)
         d0 = max(D0_0,ds)
         d2 = min(dq,d0)
         q(i,j-1) = (ds - d2)*acosp(j-1)
         q(i,j) = (d2 - dq)*acosp(j) + tiny
      endif
      enddo
    enddo
 
      do i=1,im
      if(q(i,j1).lt.D0_0) then
      ipy =  1
      dq  = - q(i,j1)*cosp(j1)
! North
      dn = q(i,j1+1)*cosp(j1+1)
      d0 = max(D0_0,dn)
      d1 = min(dq,d0)
      q(i,j1+1) = (dn - d1)*acosp(j1+1)
      q(i,j1) = (d1 - dq)*acosp(j1) + tiny
      endif
      enddo
 
      j = j2
      do i=1,im
      if(q(i,j).lt.D0_0) then
      ipy =  1
      dq  = - q(i,j)*cosp(j)
! South
      ds = q(i,j-1)*cosp(j-1)
      d0 = max(D0_0,ds)
      d2 = min(dq,d0)
      q(i,j-1) = (ds - d2)*acosp(j-1)
      q(i,j) = (d2 - dq)*acosp(j) + tiny
      endif
      enddo
 
! Check Poles.
      if(q(1,1).lt.D0_0) then
      dq = q(1,1)*cap1/real(im,r8)*acosp(j1)
      do i=1,im
      q(i,1) = tiny
      q(i,j1) = q(i,j1) + dq
      q(i,j1) = max(tiny, q(i,j1) + dq )
      enddo
      endif
 
      if(q(1,jm).lt.D0_0) then
      dq = q(1,jm)*cap1/real(im,r8)*acosp(j2)
      do i=1,im
      q(i,jm) = tiny
      q(i,j2) = max(tiny,  q(i,j2) + dq )
      enddo
      endif
!EOC 
 end subroutine filns
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: pfix --- fix an individual latitude-level
!
! !INTERFACE: 
 subroutine pfix(q, qp, im, ipx, acap, cosp2)

! !USES:
 implicit none

! !INPUT PARAMETERS:
 integer im                  ! Longitudes
 real(r8) acap               ! ???
 real(r8) cosp2              ! ???

! !INPUT/OUTPUT PARAMETERS:
 real(r8) q(im)              ! Latitude-level field to adjust
 real(r8) qp(im)             ! Second latitude-level field to adjust (usually pole)

! !OUTPUT PARAMETERS:
 integer ipx                 ! Flag:  0 if Q not change, 1 if changed


! !DESCRIPTION:
!   Fill one latitude-level from east and west neighbors
!
! !REVISION HISTORY:
!   99.03.01   Lin        Creation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
 integer i
 real(r8) summ, sump, pmean
 
   summ = D0_0
   sump = D0_0
   do i=1,im
     summ = summ + q(i)
     sump = sump + qp(i)
   enddo
 
   sump = sump/im
   pmean = (sump*acap + summ*cosp2) / (acap + cosp2*im)
 
   do i=1,im
      q(i) = pmean
      qp(i) = pmean
   enddo
 
   if( qp(1) < D0_0 ) then
      ipx = 1
   endif

!EOC
 end subroutine pfix
!-----------------------------------------------------------------------

end module fill_module
