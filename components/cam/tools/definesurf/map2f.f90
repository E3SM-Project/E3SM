      subroutine map2f(im, jm, qm, iord, jord, pfilter) 
!
! This is a stand alone 2-Grid-Wave filter for filtering the terrain for
! the finite-volume dynamical core 
! Developed and coded by S.-J. Lin
! Data Assimilation Office, NASA/GSFC
!
      implicit none
! Input
      integer, intent(in):: im          ! E-W diimension (e.g., 144 for 2.5 deg)
      integer, intent(in):: jm          ! N-S dimension (S pole to N pole; 91 for 2 deg)
      integer, intent(in):: iord        ! Mapping accuracy for E-W; recommended value=7
      integer, intent(in):: jord        ! Mapping accuracy for N-S; recommended value=3
      logical, intent(in):: pfilter     ! Polar filter (set to .T. for normal application)

! Input/Output
      real*8, intent(inout):: qm(im,jm) ! array to be filtered

! Local
      integer im2, jm2
      integer ndeg
      real*8, allocatable::  q2(:,:)
      real*8, allocatable::  lon1(:)
      real*8, allocatable::  lon2(:)
      real*8, allocatable::  sin1(:)
      real*8, allocatable::  sin2(:)
      real*8, allocatable::  qt1(:,:), qt2(:,:)

      real*8 dx1, dx2
      real*8 dy1, dy2

      integer i, j
      real*8 pi

      ndeg = 45              ! starting latitude for polar filter
      pi = 4.d0 * datan(1.d0)

      im2 = im / 2
      if (im2*2 /= im) then
          write(*,*) 'Stop in map2f; im=', im
          stop
      endif

      jm2 = (jm-1) / 2 + 1

      allocate ( qt1(im2,jm) )
      allocate ( qt2(im2,jm2) )

      allocate ( q2(im2,jm2) )
      allocate ( lon1(im+1)  )
      allocate ( lon2(im2+1) )
      allocate ( sin1(jm+1)  )
      allocate ( sin2(jm2+1) )

      dx1 = 360./im
      dx2 = 360./im2

      dy1 = pi/(jm-1)
      dy2 = pi/(jm2-1)

      do i=1,im+1
         lon1(i) = dx1 * (-0.5 + (i-1) )
      enddo

      do i=1,im2+1
         lon2(i) = dx2 * (-0.5 + (i-1) )
      enddo

         sin1(1) = -1.
         sin2(1) = -1.

         sin1(jm +1) =  1.
         sin2(jm2+1) =  1.

      do j=2,jm
         sin1(j) = dsin( -0.5*pi + dy1*(-0.5+(j-1)) )
      enddo

      do j=2,jm2
         sin2(j) = dsin( -0.5*pi + dy2*(-0.5+(j-1)) )
      enddo

      call polavg(qm, im, jm, 1, jm)
      if( pfilter ) call plft2d(im, jm, qm, 2, jm-1, ndeg)

!==============================
! From full --> half resolution
!==============================

      call xmap(iord, im, jm, sin1, lon1, qm, im2, lon2, qt1 )
      call ymap(im2, jm, sin1, qt1, jm2, sin2, qt2, 0, jord)

!==============================
! From half --> full resolution
!==============================

      call ymap(im2, jm2, sin2, qt2, jm, sin1, qt1, 0, jord)
      call xmap(iord, im2, jm, sin1, lon2, qt1, im, lon1, qm )

! Apply Monotonicity preserving polar filter
      if( pfilter ) call plft2d(im, jm, qm, 2, jm-1, ndeg)
      call polavg(qm, im, jm, 1, jm)

      deallocate ( q2 )
      deallocate ( lon1 )
      deallocate ( lon2 )
      deallocate ( sin1 )
      deallocate ( sin2 )

      deallocate ( qt1 )
      deallocate ( qt2 )

      return
      end

      subroutine polavg(p, im, jm, jfirst, jlast)

      implicit none

      integer im, jm, jfirst, jlast
      real*8 p(im,jfirst:jlast)
      real*8 sum1
      integer i

      if ( jfirst == 1 ) then 
          sum1 = 0.
        do i=1,im
          sum1 = sum1 + p(i,1)
        enddo
          sum1 = sum1/im

        do i=1,im
          p(i,1) = sum1
        enddo
      endif

      if ( jlast == jm ) then
          sum1 = 0.
        do i=1,im
          sum1 = sum1 + p(i,jm)
        enddo
          sum1 = sum1/im

        do i=1,im
          p(i,jm) = sum1
        enddo
      endif

      return
      end

      subroutine setrig(im, jm, dp, dl, cosp, cose, sinp, sine)

      implicit none

      integer im, jm
      integer j, jm1
      real*8 sine(jm),cosp(jm),sinp(jm),cose(jm)
      real*8 dp, dl
      real*8 pi, ph5

      jm1 = jm - 1
      pi  = 4.d0 * datan(1.d0)
      dl  = (pi+pi)/dble(im)
      dp  = pi/dble(jm1)

      do 10 j=2,jm
         ph5  = -0.5d0*pi + (dble(j-1)-0.5d0)*(pi/dble(jm1))
10    sine(j) = dsin(ph5)

      cosp( 1) =  0.
      cosp(jm) =  0.

      do 80 j=2,jm1
80    cosp(j) = (sine(j+1)-sine(j)) / dp

! Define cosine at edges..

      do 90 j=2,jm
90    cose(j) = 0.5 * (cosp(j-1) + cosp(j))
      cose(1) = cose(2)

      sinp( 1) = -1.
      sinp(jm) =  1.

      do 100 j=2,jm1
100   sinp(j) = 0.5 * (sine(j) + sine(j+1))

      return
      end

      subroutine ymap(im, jm, sin1, q1, jn, sin2, q2, iv, jord)

! Routine to perform area preserving mapping in N-S from an arbitrary
! resolution to another.
!
! sin1 (1) = -1 must be south pole; sin1(jm+1)=1 must be N pole.
!
! sin1(1) < sin1(2) < sin1(3) < ... < sin1(jm) < sin1(jm+1)
! sin2(1) < sin2(2) < sin2(3) < ... < sin2(jn) < sin2(jn+1)
!
! Developer: S.-J. Lin
! First version: piece-wise constant mapping
! Apr 1, 2000
! Last modified:

      implicit none

! Input
      integer im              ! original E-W dimension
      integer jm              ! original N-S dimension
      integer jn              ! Target N-S dimension
      integer jord
      integer iv              ! iv=0 scalar; iv=1: vector
      real*8  sin1(jm+1)      ! original southern edge of the cell
                              ! sin(lat1)
      real*8  sin2(jn+1)      ! Target cell's southern edge
      real*8  q1(im,jm)       ! original data at center of the cell
                              ! sin(lat2)
! Output
      real*8  q2(im,jn)       ! Mapped data at the target resolution

! Local
      integer i, j0, m, mm
      integer j

! PPM related arrays
      real*8 al(im,jm)
      real*8 ar(im,jm)
      real*8 a6(im,jm)
      real*8 dy1(jm)

      real*8 r3, r23
      parameter ( r3 = 1./3., r23 = 2./3. )
      real*8 pl, pr, qsum, esl
      real*8 dy, sum

      do j=1,jm
         dy1(j) = sin1(j+1) - sin1(j)
      enddo

! ***********************
! Area preserving mapping
! ***********************

! Construct subgrid PP distribution
      if ( jord == 1 ) then

      do j=1,jm
         do i=1,im
            a6(i,j) = 0.
            ar(i,j) = q1(i,j)
            al(i,j) = q1(i,j)
         enddo
      enddo

      else

      call ppm_lat(im, jm, q1, al, ar, a6, jord, iv)
      do i=1,im
! SP
         a6(i, 1) = 0.
         ar(i, 1) = q1(i,1)
         al(i, 1) = q1(i,1)
! NP
         a6(i,jm) = 0.
         ar(i,jm) = q1(i,jm)
         al(i,jm) = q1(i,jm)
      enddo
      endif

      do 1000 i=1,im
         j0 = 1
      do 555 j=1,jn
      do 100 m=j0,jm
!
! locate the southern edge: sin2(i)
!
      if(sin2(j) .ge. sin1(m) .and. sin2(j) .le. sin1(m+1)) then
         pl = (sin2(j)-sin1(m)) / dy1(m)
         if(sin2(j+1) .le. sin1(m+1)) then
! entire new cell is within the original cell
            pr = (sin2(j+1)-sin1(m)) / dy1(m)
            q2(i,j) = al(i,m) + 0.5*(a6(i,m)+ar(i,m)-al(i,m)) &
                          *(pr+pl)-a6(i,m)*r3*(pr*(pr+pl)+pl**2)
               j0 = m
               goto 555
          else
! South most fractional area
            qsum = (sin1(m+1)-sin2(j))*(al(i,m)+0.5*(a6(i,m)+ &
                    ar(i,m)-al(i,m))*(1.+pl)-a6(i,m)* &
                     (r3*(1.+pl*(1.+pl))))
              do mm=m+1,jm
! locate the eastern edge: sin2(j+1)
                 if(sin2(j+1) .gt. sin1(mm+1) ) then
! Whole layer
                     qsum = qsum + dy1(mm)*q1(i,mm)
                 else
! North most fractional area
                     dy = sin2(j+1)-sin1(mm)
                    esl = dy / dy1(mm)
                   qsum = qsum + dy*(al(i,mm)+0.5*esl*  &
                         (ar(i,mm)-al(i,mm)+a6(i,mm)*(1.-r23*esl)))
                     j0 = mm
                     goto 123
                 endif
              enddo
              goto 123
           endif
      endif
100   continue
123   q2(i,j) = qsum / ( sin2(j+1) - sin2(j) )
555   continue
1000  continue

! Final processing for poles

      if ( iv == 0 ) then

! South pole
           sum = 0.
         do i=1,im
           sum = sum + q2(i,1)
         enddo

           sum = sum / im
         do i=1,im
           q2(i,1) = sum
         enddo

! North pole:
           sum = 0.
         do i=1,im
           sum = sum + q2(i,jn)
         enddo

           sum = sum / im
         do i=1,im
           q2(i,jn) = sum
         enddo

      endif

      return
      end

      subroutine ppm_lat(im, jm, q, al, ar, a6, jord, iv)
      implicit none

!INPUT
      integer im, jm                      !  Dimensions
      real*8  q(im,jm)
      real*8 al(im,jm)
      real*8 ar(im,jm)
      real*8 a6(im,jm)
      integer jord
      integer iv                             ! iv=0 scalar
                                             ! iv=1 vector
! Local
      real*8 dm(im,jm)
      real*8    r3
      parameter ( r3 = 1./3. )
      integer i, j, im2, iop, jm1
      real*8 tmp, qmax, qmin
      real*8 qop

! Compute dm: linear slope

      do j=2,jm-1
         do i=1,im
            dm(i,j) = 0.25*(q(i,j+1) - q(i,j-1))
            qmax = max(q(i,j-1),q(i,j),q(i,j+1)) - q(i,j)
            qmin = q(i,j) - min(q(i,j-1),q(i,j),q(i,j+1))
            dm(i,j) = sign(min(abs(dm(i,j)),qmin,qmax),dm(i,j))
         enddo
      enddo

      im2 = im/2
      jm1 = jm - 1

!Poles:
      if (iv == 1 ) then
! SP
          do i=1,im
              if( i .le. im2) then
                  qop = -q(i+im2,2)
              else
                  qop = -q(i-im2,2)
              endif
              tmp = 0.25*(q(i,2) - qop)
              qmax = max(q(i,2),q(i,1), qop) - q(i,1)
              qmin = q(i,1) - min(q(i,2),q(i,1), qop)
              dm(i,1) = sign(min(abs(tmp),qmax,qmin),tmp)
           enddo
! NP
           do i=1,im
              if( i .le. im2) then
                  qop = -q(i+im2,jm1)
              else
                  qop = -q(i-im2,jm1)
              endif
              tmp = 0.25*(qop - q(i,jm1))
              qmax = max(qop,q(i,jm), q(i,jm1)) - q(i,jm)
              qmin = q(i,jm) - min(qop,q(i,jm), q(i,jm1))
              dm(i,jm) = sign(min(abs(tmp),qmax,qmin),tmp)
           enddo
      else
!
!*********
! Scalar:
!*********
! SP
          do i=1,im2
            tmp = 0.25*(q(i,2)-q(i+im2,2))
            qmax = max(q(i,2),q(i,1), q(i+im2,2)) - q(i,1)
            qmin = q(i,1) - min(q(i,2),q(i,1), q(i+im2,2))
            dm(i,1) = sign(min(abs(tmp),qmax,qmin),tmp)
          enddo

          do i=im2+1,im
            dm(i, 1) =  - dm(i-im2, 1)
          enddo
! NP
          do i=1,im2
            tmp = 0.25*(q(i+im2,jm1)-q(i,jm1))
            qmax = max(q(i+im2,jm1),q(i,jm), q(i,jm1)) - q(i,jm)
            qmin = q(i,jm) - min(q(i+im2,jm1),q(i,jm), q(i,jm1))
            dm(i,jm) = sign(min(abs(tmp),qmax,qmin),tmp)
          enddo

          do i=im2+1,im
            dm(i,jm) =  - dm(i-im2,jm)
          enddo
      endif

      do j=2,jm
        do i=1,im
          al(i,j) = 0.5*(q(i,j-1)+q(i,j)) + r3*(dm(i,j-1) - dm(i,j))
        enddo
      enddo

      do j=1,jm-1
        do i=1,im
          ar(i,j) = al(i,j+1)
        enddo
      enddo

      do j=2,jm-1
         do i=1,im
            a6(i,j) = 3.*(q(i,j)+q(i,j) - (al(i,j)+ar(i,j)))
         enddo

        call lmppm(dm(1,j), a6(1,j), ar(1,j),  &
                   al(1,j),  q(1,j), im, jord-3)
      enddo

      return
      end

      subroutine xmap(iord, im, jm, sin1, lon1, q1, in, lon2, q2)

! Routine to perform area preserving mapping in E-W from an arbitrary
! resolution to another.
! Periodic domain will be assumed, i.e., the eastern wall bounding cell
! im is lon1(im+1) = lon1(1); Note the equal sign is true geographysically.
!
! lon1(1) < lon1(2) < lon1(3) < ... < lon1(im) < lon1(im+1)
! lon2(1) < lon2(2) < lon2(3) < ... < lon2(in) < lon2(in+1)
!
! Developer: S.-J. Lin
! First version: piece-wise constant mapping
! Apr 1, 2000
! Last modified:

      implicit none

! Input
      integer iord
      integer im              ! original E-W dimension
      integer in              ! Target E-W dimension
      integer jm              ! original N-S dimension
      real*8    lon1(im+1)      ! original western edge of the cell
      real*8    sin1(jm+1)
      real*8    q1(im,jm)       ! original data at center of the cell
      real*8    lon2(in+1)      ! Target cell's western edge

! Output
      real*8    q2(in,jm)       ! Mapped data at the target resolution

! Local
      integer i1, i2
      integer i, i0, m, mm
      integer j
      integer ird

! PPM related arrays
      real*8 qtmp(-im:im+im)
      real*8   al(-im:im+im)
      real*8   ar(-im:im+im)
      real*8   a6(-im:im+im)
      real*8   x1(-im:im+im+1)
      real*8  dx1(-im:im+im)
      real*8  r3, r23
      parameter ( r3 = 1./3., r23 = 2./3. )
      real*8 pl, pr, qsum, esl
      real*8 dx
      logical found

      do i=1,im+1
         x1(i) = lon1(i)
      enddo

      do i=1,im
         dx1(i) = x1(i+1) - x1(i)
      enddo

! check to see if ghosting is necessary

!**************
! Western edge:
!**************
          found = .false.
          i1 = 1
      do while ( .not. found )
         if( lon2(1) .ge. x1(i1) ) then
             found = .true.
         else
                  i1 = i1 - 1
             if (i1 .lt. -im) then
                 write(6,*) 'failed in xmap'
                 stop
             else
                 x1(i1) = x1(i1+1) - dx1(im+i1)
                dx1(i1) = dx1(im+i1)
             endif
         endif
      enddo

!**************
! Eastern edge:
!**************
          found = .false.
          i2 = im+1
      do while ( .not. found )
         if( lon2(in+1) .le. x1(i2) ) then
             found = .true.
         else
                  i2 = i2 + 1
             if (i2 .gt. 2*im) then
                 write(6,*) 'failed in xmap'
                 stop
             else
                dx1(i2-1) = dx1(i2-1-im)
                 x1(i2) = x1(i2-1) + dx1(i2-1)
             endif
         endif
      enddo

      do 1000 j=1,jm

! ***********************
! Area preserving mapping
! ***********************

! Construct subgrid PP distribution
      if ( abs(sin1(j)+sin1(j+1)) > 1.5 ) then
           ird = 3
      elseif ( abs(sin1(j)+sin1(j+1)) < 1.0 ) then
           ird = 8
      else
           ird = iord
      endif

      if ( iord == 1 ) then
           do i=1,im
              qtmp(i) = q1(i,j)
                al(i) = q1(i,j)
                ar(i) = q1(i,j)
                a6(i) = 0.
           enddo
              qtmp(0   ) = q1(im,j)
              qtmp(im+1) = q1(1, j)
      else
      call ppm_cycle(im, q1(1,j), al(1), ar(1), a6(1), qtmp, ird)
      endif

! check to see if ghosting is necessary

! Western edge
          if ( i1 .le. 0 ) then
               do i=i1,0
                  qtmp(i) = qtmp(im+i)
                    al(i) = al(im+i)
                    ar(i) = ar(im+i)
                    a6(i) = a6(im+i)
               enddo
          endif

! Eastern edge:
          if ( i2 .gt. im+1 ) then
             do i=im+1,i2-1
                qtmp(i) = qtmp(i-im)
                  al(i) =   al(i-im)
                  ar(i) =   ar(i-im)
                  a6(i) =   a6(i-im)
             enddo
          endif

         i0 = i1

      do 555 i=1,in
      do 100 m=i0,i2-1
!
! locate the western edge: lon2(i)
!
      if(lon2(i) .ge. x1(m) .and. lon2(i) .le. x1(m+1)) then
         pl = (lon2(i)-x1(m)) / dx1(m)
         if(lon2(i+1) .le. x1(m+1)) then
! entire new grid is within the original grid
            pr = (lon2(i+1)-x1(m)) / dx1(m)
            q2(i,j) = al(m) + 0.5*(a6(m)+ar(m)-al(m))  &
                          *(pr+pl)-a6(m)*r3*(pr*(pr+pl)+pl**2)
               i0 = m
               goto 555
          else
! Left most fractional area
            qsum = (x1(m+1)-lon2(i))*(al(m)+0.5*(a6(m)+ &
                    ar(m)-al(m))*(1.+pl)-a6(m)*         &
                     (r3*(1.+pl*(1.+pl))))
              do mm=m+1,i2-1
! locate the eastern edge: lon2(i+1)
                 if(lon2(i+1) .gt. x1(mm+1) ) then
! Whole layer
                     qsum = qsum + dx1(mm)*qtmp(mm)
                 else
! Right most fractional area
                     dx = lon2(i+1)-x1(mm)
                    esl = dx / dx1(mm)
                   qsum = qsum + dx*(al(mm)+0.5*esl*  &
                         (ar(mm)-al(mm)+a6(mm)*(1.-r23*esl)))
                     i0 = mm
                     goto 123
                 endif
              enddo
              goto 123
           endif
      endif
100   continue
123   q2(i,j) = qsum / ( lon2(i+1) - lon2(i) )
555   continue
1000  continue

      return
      end

      subroutine ppm_cycle(im, q, al, ar, a6, p, iord)
      implicit none

      real*8 r3
      parameter ( r3 = 1./3. )

! Input 
      integer im, iord
      real*8  q(1)
! Output
      real*8 al(1)
      real*8 ar(1)
      real*8 a6(1)
      real*8 p(-im:im+im)

! local
      real*8  dm(0:im)
      integer i, lmt
      real*8 tmp, qmax, qmin

         p(0) = q(im)
      do i=1,im
         p(i) = q(i)
      enddo
         p(im+1) = q(1)

! 2nd order slope
      do i=1,im
         tmp = 0.25*(p(i+1) - p(i-1))
         qmax = max(p(i-1), p(i), p(i+1)) - p(i)
         qmin = p(i) - min(p(i-1), p(i), p(i+1))
         dm(i) = sign(min(abs(tmp),qmax,qmin), tmp)
      enddo
         dm(0) = dm(im)

      do i=1,im
         al(i) = 0.5*(p(i-1)+p(i)) + (dm(i-1) - dm(i))*r3
      enddo

      do i=1,im-1
         ar(i) = al(i+1)
      enddo
         ar(im) = al(1)

         do i=1,im
            a6(i) = 3.*(p(i)+p(i)  - (al(i)+ar(i)))
         enddo

      if(iord <= 6) then
         lmt = iord - 3
         if(lmt <= 2) call lmppm(dm(1),a6(1),ar(1),al(1),p(1),im,lmt)
      else
         call huynh(im, ar(1), al(1), p(1), a6(1), dm(1))
         call lmppm(dm(1),a6(1),ar(1),al(1),p(1),im,2)
      endif

      return
      end

      subroutine lmppm(dm, a6, ar, al, p, im, lmt)
      implicit none
      real*8 r12
      parameter ( r12 = 1./12. )

      integer im, lmt
      integer i
      real*8 a6(im),ar(im),al(im),p(im),dm(im)
      real*8 da1, da2, fmin, a6da

! LMT = 0: full monotonicity
! LMT = 1: semi-monotonic constraint (no undershoot)
! LMT = 2: positive-definite constraint

      if(lmt.eq.0) then

! Full constraint
      do 100 i=1,im
      if(dm(i) .eq. 0.) then
         ar(i) = p(i)
         al(i) = p(i)
         a6(i) = 0.
      else
         da1  = ar(i) - al(i)
         da2  = da1**2
         a6da = a6(i)*da1
         if(a6da .lt. -da2) then
            a6(i) = 3.*(al(i)-p(i))
            ar(i) = al(i) - a6(i)
         elseif(a6da .gt. da2) then
            a6(i) = 3.*(ar(i)-p(i))
            al(i) = ar(i) - a6(i)
         endif
      endif
100   continue

      elseif(lmt == 1) then
! Semi-monotonic constraint
      do 150 i=1,im
      if(abs(ar(i)-al(i)) .ge. -a6(i)) go to 150
      if(p(i).lt.ar(i) .and. p(i).lt.al(i)) then
            ar(i) = p(i)
            al(i) = p(i)
            a6(i) = 0.
      elseif(ar(i) .gt. al(i)) then
            a6(i) = 3.*(al(i)-p(i))
            ar(i) = al(i) - a6(i)
      else
            a6(i) = 3.*(ar(i)-p(i))
            al(i) = ar(i) - a6(i)
      endif
150   continue
      elseif(lmt == 2) then
! Positive definite constraint
      do 250 i=1,im
      if(abs(ar(i)-al(i)) >= -a6(i)) go to 250
      fmin = p(i) + 0.25*(ar(i)-al(i))**2/a6(i) + a6(i)*r12
      if(fmin >= 0.) go to 250
      if(p(i).lt.ar(i) .and. p(i).lt.al(i)) then
            ar(i) = p(i)
            al(i) = p(i)
            a6(i) = 0.
      elseif(ar(i) .gt. al(i)) then
            a6(i) = 3.*(al(i)-p(i))
            ar(i) = al(i) - a6(i)
      else
            a6(i) = 3.*(ar(i)-p(i))
            al(i) = ar(i) - a6(i)
      endif
250   continue
      endif
      return
      end

      subroutine huynh(im, ar, al, p, d2, d1)

! Enforce Huynh's 2nd constraint in 1D periodic domain

      implicit none
      integer im, i
      real*8 ar(im)
      real*8 al(im)
      real*8  p(im)
      real*8 d2(im)
      real*8 d1(im)

! Local scalars:
      real*8 pmp
      real*8 lac
      real*8 pmin
      real*8 pmax

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
         pmp   = p(1) + 2.0 * d1(1)
         lac   = p(1) + 0.5 * (d1(1)+d2(im)) + d2(im) 
         pmin  = min(p(1), pmp, lac)
         pmax  = max(p(1), pmp, lac)
         ar(1) = min(pmax, max(ar(1), pmin))

      do i=2, im
         pmp   = p(i) + 2.0*d1(i)
         lac   = p(i) + 0.5*(d1(i)+d2(i-1)) + d2(i-1)
         pmin  = min(p(i), pmp, lac)
         pmax  = max(p(i), pmp, lac)
         ar(i) = min(pmax, max(ar(i), pmin))
      enddo

! Constraint for AL
      do i=1, im-1
         pmp   = p(i) - 2.0*d1(i+1)
         lac   = p(i) + 0.5*(d2(i+1)-d1(i+1)) + d2(i+1)
         pmin  = min(p(i), pmp, lac)
         pmax  = max(p(i), pmp, lac)
         al(i) = min(pmax, max(al(i), pmin))
      enddo

! i=im
         i = im
         pmp    = p(im) - 2.0*d1(1)
         lac    = p(im) + 0.5*(d2(1)-d1(1)) + d2(1)
         pmin   = min(p(im), pmp, lac)
         pmax   = max(p(im), pmp, lac)
         al(im) = min(pmax, max(al(im), pmin))

! compute A6 (d2)
      do i=1, im
         d2(i) = 3.*(p(i)+p(i)  - (al(i)+ar(i)))
      enddo
      return
      end

 subroutine plft2d(im, jm, p, JS, JN, ndeg)
!
! This is a weak LOCAL polar filter.
! Developer: Shian-Jiann Lin

      implicit none

      integer im
      integer jm
      integer js, jn, ndeg
      real*8 p(im,jm)

      integer i, j, n, ideg, jj, jc
      real*8 cosp(jm),cose(jm)
      real*8 a(0:im/2+1)

      real*8 sine(jm),sinp(jm)
      real*8, allocatable, save :: se(:), sc(:)

      real*8 pi, dp, dl, e0, ycrit, coszc, smax, rn, rn2, esl, tmp

      data IDEG /0/

      if(IDEG .ne. ndeg) then
	IDEG = ndeg
!  (e0 = 2.6)
      e0 = 0.5 * sqrt(27.)
      PI = 4. * ATAN(1.)

      allocate( sc(jm), se(jm))

      call setrig(im, jm, dp, dl, cosp, cose, sinp, sine)

      ycrit = IDEG*PI/180.
      coszc = cos(ycrit)

      smax = (jm-1)/2
      write(6,*) 'Critical latitude in local pft = ',ndeg

         a(0) = 1.
      do n=1,im/2+1
         rn = n
         rn2 = 2*n
         a(n) = sqrt(rn2+1.) * ((rn2+1.)/rn2)**rn
      enddo

      do j=2,jm-1
      sc(j) = coszc / cosp(j)

      IF(sc(j) > 1. .and. sc(j) <= 1.5 ) THEN
         esl = 1./ sc(j)
         sc(j) =  1. +  (1.-esl) / (1.+esl)
      ELSEIF(sc(j) > 1.5 .and. sc(j) <= e0 ) THEN
         esl = 1./ sc(j)
         sc(j) =  1. + 2./ (27.*esl**2 - 2.)
      ELSEIF(sc(j) > e0) THEN
! Search
      do jj=1,im/2
      if(sc(j) <= a(jj)) then
      jc = jj
!     write(*,*) 'jc=', jc
      goto 111
      endif
      enddo
      jc = im/2 + 1
111   continue

      tmp = ((sc(j) - a(jc-1))/(a(jc) - a(jc-1)))**0.25
      sc(j) =  jc + min(1.d0, tmp)
!     sc(j) =  min(smax,sc(j))
      ENDIF
      enddo
! ====================================================
      do j=2,jm
      se(j) = coszc / cose(j)
      IF(se(j) > 1. .and. se(j) <= 1.5 ) THEN
         esl = 1./ se(j)
         se(j) =  1. + (1.-esl) / (1.+esl)
      ELSEIF(se(j) > 1.5 .and. se(j) <= e0 ) THEN
         esl = 1./ se(j)
         se(j) =  1. + 2./ (27.*esl**2 - 2.)
      ELSEIF(se(j) > e0) THEN
! Search
      do jj=1,im/2
      if(se(j) <= a(jj)) then
      jc = jj
      goto 222
      endif
      enddo

      jc = im/2 + 1
222   continue
      tmp = ((se(j) - a(jc-1))/(a(jc) - a(jc-1)))**0.25
      se(j) =  jc + min(1.d0, tmp)
!     se(j) =  min(smax,se(j))
      ENDIF
      enddo

      do i=1,im
        se( 2) = sc(2)
        se(jm) = sc(jm-1)
      enddo

      do j=2,jm-1
!        write(*,*) j,sc(j)
      enddo
      ENDIF

      if( JN == (jm-1) ) then
! Cell-centered variables
         call lpft(im, jm, p, 2, jm-1, Sc)
      else
! Cell-edge variables
         call lpft(im, jm, p, 2, jm, Se)
      endif
      return
 end


 subroutine lpft(im, jm, p, j1, j2, s)
      implicit none

      integer im, jm, j1, j2
      real*8  p(im,jm)
      real*8  s(jm)

! Local
      integer i, j, n, nt

      real*8 ptmp(0:im+1)
      real*8 q(0:im+1)
      real*8 frac, rsc, bt

      do 2500 j=j1,j2
      if(s(j) > 1.02) then

        NT  = INT(S(j))
        frac = S(j) - NT
        NT = NT-1

        rsc = 1. / (1.+frac)
        bt = 0.5 * frac

        do i=1,im
           ptmp(i) = p(i,j)
        enddo

        ptmp(0)    = p(im,j)
        ptmp(im+1) = p(1 ,j)

        if( NT < 1 ) then
          do i=1,im
             p(i,j) = rsc * (ptmp(i) + bt*(ptmp(i-1)+ptmp(i+1)))
          enddo
        else
          do i=1,im
             q(i) = rsc * (ptmp(i) + bt*(ptmp(i-1)+ptmp(i+1)))
          enddo

           do 500 N=1,NT
              q(0) = q(im)
              do i=1,im
                 ptmp(i) = q(i) + q(i-1)
              enddo
                 ptmp(im+1) = ptmp(1)

             if ( n == nt ) then
               do i=1,im
                  p(i,j) = 0.25*(ptmp(i) + ptmp(i+1))
               enddo
             else
               do i=1,im
                  q(i) = 0.25*(ptmp(i) + ptmp(i+1))
               enddo
             endif
500        continue
        endif
      endif
2500  continue

      return
 end
