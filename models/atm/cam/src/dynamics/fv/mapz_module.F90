module mapz_module

 use shr_kind_mod,  only : r8 => shr_kind_r8
 use FVperf_module, only : FVstartclock, FVstopclock
 use abortutils,    only : endrun
 use cam_logfile,   only : iulog

 public map1_cubic_te, map1_ppm, mapn_ppm, mapn_ppm_tracer, ppme

 private

 real(r8), parameter ::  D0_0                    =  0.0_r8
 real(r8), parameter ::  D1EM14                  =  1.0e-14_r8
 real(r8), parameter ::  D0_125                  =  0.125_r8
 real(r8), parameter ::  D0_1875                 =  0.1875_r8
 real(r8), parameter ::  D0_25                   =  0.25_r8
 real(r8), parameter ::  D0_5                    =  0.5_r8
 real(r8), parameter ::  D1_0                    =  1.0_r8
 real(r8), parameter ::  D1_5                    =  1.5_r8
 real(r8), parameter ::  D2_0                    =  2.0_r8
 real(r8), parameter ::  D3_0                    =  3.0_r8
 real(r8), parameter ::  D4_0                    =  4.0_r8
 real(r8), parameter ::  D5_0                    =  5.0_r8
 real(r8), parameter ::  D8_0                    =  8.0_r8
 real(r8), parameter ::  D12_0                   = 12.0_r8

contains

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  map1_cubic_te --- Cubic Interpolation for TE mapping
!
! !INTERFACE:
  subroutine map1_cubic_te ( km,   pe1,   q1,  kn,   pe2,   q2,         &
                             ng_s, ng_n, itot, i1, i2,                  &
                             j, jfirst, jlast, iv, kord)
      implicit none

! !INPUT PARAMETERS:
      integer, intent(in) :: i1                ! Starting longitude
      integer, intent(in) :: i2                ! Finishing longitude
      integer, intent(in) :: itot              ! Total latitudes
      integer, intent(in) :: iv                ! Mode: 0 ==  constituents  1 == ???
      integer, intent(in) :: kord              ! Method order
      integer, intent(in) :: j                 ! Current latitude
      integer, intent(in) :: jfirst            ! Starting latitude
      integer, intent(in) :: jlast             ! Finishing latitude
      integer, intent(in) :: ng_s              ! Ghosted latitudes south
      integer, intent(in) :: ng_n              ! Ghosted latitudes north
      integer, intent(in) :: km                ! Original vertical dimension
      integer, intent(in) :: kn                ! Target vertical dimension

      real(r8), intent(in) ::  pe1(itot,km+1)  ! pressure at layer edges 
                                               ! (from model top to bottom surface)
                                               ! in the original vertical coordinate
      real(r8), intent(in) ::  pe2(itot,kn+1)  ! pressure at layer edges 
                                               ! (from model top to bottom surface)
                                               ! in the new vertical coordinate

      real(r8), intent(in)   ::  q1(itot,jfirst-ng_s:jlast+ng_n,km) ! Field input

! !INPUT/OUTPUT PARAMETERS:
      real(r8), intent(inout)::  q2(itot,jfirst-ng_s:jlast+ng_n,kn) ! Field output

! !DESCRIPTION:
!
!     Perform Cubic Interpolation a given latitude
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
!
! !REVISION HISTORY:
!    05.11.14   Takacs    Initial Code
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      real(r8)       qx(i1:i2,km)
      real(r8)   logpl1(i1:i2,km)
      real(r8)   logpl2(i1:i2,kn)
      real(r8)   dlogp1(i1:i2,km)
      real(r8)    vsum1(i1:i2)
      real(r8)    vsum2(i1:i2)
      real(r8)   am2,am1,ap0,ap1,P,PLP1,PLP0,PLM1,PLM2,DLP0,DLM1,DLM2

      integer i, k, LM2,LM1,LP0,LP1

! Initialization
! --------------
      do k=1,km
          qx(:,k) = q1(:,j,k)
      logpl1(:,k) = log( D0_5*(pe1(:,k)+pe1(:,k+1)) )
      enddo
      do k=1,kn
      logpl2(:,k) = log( D0_5*(pe2(:,k)+pe2(:,k+1)) )
      enddo

      do k=1,km-1
      dlogp1(:,k) = logpl1(:,k+1)-logpl1(:,k)
      enddo

! Compute vertical integral of Input TE
! -------------------------------------
      vsum1(:) = D0_0
      do i=i1,i2
      do k=1,km
      vsum1(i) = vsum1(i) + qx(i,k)*( pe1(i,k+1)-pe1(i,k) )
      enddo
      vsum1(i) = vsum1(i) / ( pe1(i,km+1)-pe1(i,1) )
      enddo

! Interpolate TE onto target Pressures
! ------------------------------------
      do i=i1,i2
      do k=1,kn
         LM1 = 1
         LP0 = 1
         do while( logpl1(i,LP0).lt.logpl2(i,k) .and. LP0.le.km )
         LP0 = LP0+1
         enddo
         LM1 = max(LP0-1,1)
         LP0 = min(LP0, km)

! Extrapolate Linearly in LogP above first model level
! ----------------------------------------------------
         if( LM1.eq.1 .and. LP0.eq.1 ) then
             q2(i,j,k) = qx(i,1) + ( qx(i,2)-qx(i,1) )*( logpl2(i,k)-logpl1(i,1) ) &
                                                      /( logpl1(i,2)-logpl1(i,1) )

! Extrapolate Linearly in LogP below last model level
! ---------------------------------------------------
         else if( LM1.eq.km .and. LP0.eq.km ) then
             q2(i,j,k) = qx(i,km) + ( qx(i,km)-qx(i,km-1) )*( logpl2(i,k )-logpl1(i,km  ) ) &
                                                           /( logpl1(i,km)-logpl1(i,km-1) )

! Interpolate Linearly in LogP between levels 1 => 2 and km-1 => km
! -----------------------------------------------------------------
         else if( LM1.eq.1 .or. LP0.eq.km ) then
             q2(i,j,k) = qx(i,LP0) + ( qx(i,LM1)-qx(i,LP0) )*( logpl2(i,k  )-logpl1(i,LP0) ) &
                                                            /( logpl1(i,LM1)-logpl1(i,LP0) )

! Interpolate Cubicly in LogP between other model levels
! ------------------------------------------------------
         else
              LP1 = LP0+1
              LM2 = LM1-1
             P    = logpl2(i,k)
             PLP1 = logpl1(i,LP1)
             PLP0 = logpl1(i,LP0)
             PLM1 = logpl1(i,LM1)
             PLM2 = logpl1(i,LM2)
             DLP0 = dlogp1(i,LP0)
             DLM1 = dlogp1(i,LM1)
             DLM2 = dlogp1(i,LM2)

              ap1 = (P-PLP0)*(P-PLM1)*(P-PLM2)/( DLP0*(DLP0+DLM1)*(DLP0+DLM1+DLM2) )
              ap0 = (PLP1-P)*(P-PLM1)*(P-PLM2)/( DLP0*      DLM1 *(     DLM1+DLM2) )
              am1 = (PLP1-P)*(PLP0-P)*(P-PLM2)/( DLM1*      DLM2 *(DLP0+DLM1     ) )
              am2 = (PLP1-P)*(PLP0-P)*(PLM1-P)/( DLM2*(DLM1+DLM2)*(DLP0+DLM1+DLM2) )

             q2(i,j,k) = ap1*qx(i,LP1) + ap0*qx(i,LP0) + am1*qx(i,LM1) + am2*qx(i,LM2)

         endif

      enddo
      enddo

! Compute vertical integral of Output TE
! --------------------------------------
      vsum2(:) = D0_0
      do i=i1,i2
      do k=1,kn
      vsum2(i) = vsum2(i) + q2(i,j,k)*( pe2(i,k+1)-pe2(i,k) )
      enddo
      vsum2(i) = vsum2(i) / ( pe2(i,kn+1)-pe2(i,1) )
      enddo

! Adjust Final TE to conserve
! ---------------------------
      do i=i1,i2
      do k=1,kn
         q2(i,j,k) = q2(i,j,k) + vsum1(i)-vsum2(i)
!        q2(i,j,k) = q2(i,j,k) * vsum1(i)/vsum2(i)
      enddo
      enddo

      return
!EOC
 end subroutine map1_cubic_te
!----------------------------------------------------------------------- 

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  map1_ppm --- Piecewise parabolic mapping, variant 1
!
! !INTERFACE:
  subroutine map1_ppm( km,   pe1,   q1,  kn,   pe2,   q2,                &
                       ng_s, ng_n, itot, i1, i2,                         &
                       j, jfirst, jlast, iv, kord)

      implicit none

! !INPUT PARAMETERS:
      integer, intent(in) :: i1                ! Starting longitude
      integer, intent(in) :: i2                ! Finishing longitude
      integer, intent(in) :: itot              ! Total latitudes
      integer, intent(in) :: iv                ! Mode: 0 ==  constituents  1 == ???
      integer, intent(in) :: kord              ! Method order
      integer, intent(in) :: j                 ! Current latitude
      integer, intent(in) :: jfirst            ! Starting latitude
      integer, intent(in) :: jlast             ! Finishing latitude
      integer, intent(in) :: ng_s              ! Ghosted latitudes south
      integer, intent(in) :: ng_n              ! Ghosted latitudes north
      integer, intent(in) :: km                ! Original vertical dimension
      integer, intent(in) :: kn                ! Target vertical dimension

      real(r8), intent(in) ::  pe1(itot,km+1)  ! pressure at layer edges 
                                               ! (from model top to bottom surface)
                                               ! in the original vertical coordinate
      real(r8), intent(in) ::  pe2(itot,kn+1)  ! pressure at layer edges 
                                               ! (from model top to bottom surface)
                                               ! in the new vertical coordinate
      real(r8), intent(in) ::  q1(itot,jfirst-ng_s:jlast+ng_n,km) ! Field input

! !INPUT/OUTPUT PARAMETERS:
      real(r8), intent(inout)::  q2(itot,jfirst-ng_s:jlast+ng_n,kn) ! Field output

! !DESCRIPTION:
!
!     Perform piecewise parabolic method on a given latitude    
! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
!
! !REVISION HISTORY: 
!    00.04.24   Lin       Last modification
!    01.03.26   Sawyer    Added ProTeX documentation
!    02.04.04   Sawyer    incorporated latest FVGCM version
!    02.06.20   Sawyer    made Q2 inout since the args for Q1/Q2 same
!    03.07.22   Parks     Cleaned main loop, removed gotos
!    05.05.25   Sawyer    Merged CAM and GEOS5 versions
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      real(r8)       r3, r23
      parameter (r3 = D1_0/D3_0, r23 = D2_0/D3_0)
      real(r8)   dp1(i1:i2,km)
      real(r8)  q4(4,i1:i2,km)

      integer i, k, kk, kl, k0(i1:i2,0:kn+1), k0found
      real(r8)    pl, pr, qsum, qsumk(i1:i2,kn), delp, esl

      do k=1,km
         do i=i1,i2
             dp1(i,k) = pe1(i,k+1) - pe1(i,k)
            q4(1,i,k) = q1(i,j,k)
         enddo
      enddo

! Mapping

! Compute vertical subgrid distribution
      call ppm2m( q4, dp1, km, i1, i2, iv, kord )

! For each pe2(i,k), determine lowest pe1 interval = smallest k0 (= k0(i,k))
!   such that pe1(i,k0) <= pe2(i,k) <= pe1(i,k0+1)
!   Note that pe2(i,1)==pe1(i,1) and pe2(i,kn+1)==pe1(i,kn+1)
!   Note also that pe1, pe2 are assumed to be monotonically increasing
#if defined( UNICOSMP ) || defined ( NEC_SX )
      do kk = km, 1, -1
         do k = 1, kn+1
!dir$ prefervector
            do i = i1, i2
               if (pe2(i,k) <= pe1(i,kk+1)) then
                  k0(i,k) = kk
               endif
            enddo
         enddo
      enddo
#else
      do i = i1, i2
         k0(i,0) = 1
         do k = 1, kn+1
            k0found = -1
            do kk = k0(i,k-1), km
               if (pe2(i,k) <= pe1(i,kk+1)) then
                  k0(i,k) = kk
                  k0found = kk
                  exit
               endif
            enddo
            if (k0found .lt. 0) then
               write(iulog,*) 'mapz error - k0found i j k (kk,pe1,pe2) = ',   &
                  k0found, i, j, k, (kk,pe1(i,kk),pe2(i,kk),kk=1,km+1)
               call endrun('MAPZ_MODULE')
               return
            endif
         enddo
      enddo
#endif

! Interpolate
      do k = 1, kn

! Prepare contribution between pe1(i,ko(i,k)+1) and pe1(i,k0(i,k+1))
         qsumk(:,k) = D0_0
         do i = i1, i2
            do kl = k0(i,k)+1, k0(i,k+1)-1
               qsumk(i,k) = qsumk(i,k) + dp1(i,kl)*q4(1,i,kl)
            enddo
         enddo

         do i = i1, i2
            kk = k0(i,k)
! Consider contribution between pe1(i,kk) and pe2(i,k)
            pl = (pe2(i,k)-pe1(i,kk)) / dp1(i,kk)
! Check to see if pe2(i,k+1) and pe2(i,k) are in same pe1 interval
            if (k0(i,k+1) == k0(i,k)) then
               pr = (pe2(i,k+1)-pe1(i,kk)) / dp1(i,kk)
               q2(i,j,k) = q4(2,i,kk) + D0_5*(q4(4,i,kk)+q4(3,i,kk)-q4(2,i,kk))  &
                  *(pr+pl) - q4(4,i,kk)*r3*(pr*(pr+pl)+pl**2)
            else
! Consider contribution between pe2(i,k) and pe1(i,kk+1)
               qsum = (pe1(i,kk+1)-pe2(i,k))*(q4(2,i,kk)+D0_5*(q4(4,i,kk)+       &
                  q4(3,i,kk)-q4(2,i,kk))*(D1_0+pl)-q4(4,i,kk)*                    &
                  (r3*(D1_0+pl*(D1_0+pl))))
! Next consider contribution between pe1(i,kk+1) and pe1(i,k0(i,k+1))
               qsum = qsum + qsumk(i,k)
! Now consider contribution between pe1(i,k0(i,k+1)) and pe2(i,k+1)
               kl = k0(i,k+1)
               delp = pe2(i,k+1)-pe1(i,kl)
               esl = delp / dp1(i,kl)
               qsum = qsum + delp*(q4(2,i,kl)+D0_5*esl*                          &
                  (q4(3,i,kl)-q4(2,i,kl)+q4(4,i,kl)*(D1_0-r23*esl)))
               q2(i,j,k) = qsum / ( pe2(i,k+1) - pe2(i,k) )
            endif
         enddo
      enddo

      return
!EOC
 end subroutine map1_ppm
!----------------------------------------------------------------------- 

!-----------------------------------------------------------------------
!BOP
! !ROUTINE:  mapn_ppm --- Piecewise parabolic mapping, variant 1
!
! !INTERFACE:
 subroutine mapn_ppm(km,   pe1,   q1, nq,                  &
                     kn,   pe2,   q2, ng_s, ng_n,          &
                     itot, i1, i2, j,                      &
                     jfirst, jlast, iv, kord)

! !USES:
      implicit none

! !INPUT PARAMETERS:
      integer, intent(in) :: i1                ! Starting longitude
      integer, intent(in) :: i2                ! Finishing longitude
      integer, intent(in) :: itot              ! Total latitudes
      integer, intent(in) :: iv                ! Mode: 0 ==  constituents  1 == ???
      integer, intent(in) :: kord              ! Method order
      integer, intent(in) :: j                 ! Current latitude
      integer, intent(in) :: jfirst            ! Starting latitude
      integer, intent(in) :: jlast             ! Finishing latitude
      integer, intent(in) :: ng_s              ! Ghosted latitudes south
      integer, intent(in) :: ng_n              ! Ghosted latitudes north
      integer, intent(in) :: km                ! Original vertical dimension
      integer, intent(in) :: kn                ! Target vertical dimension
      integer, intent(in) :: nq                ! Number of tracers

      real(r8), intent(in) :: pe1(itot,km+1)   ! pressure at layer edges 
                                               ! (from model top to bottom surface)
                                               ! in the original vertical coordinate
      real(r8), intent(in) :: pe2(itot,kn+1)   ! pressure at layer edges 
                                               ! (from model top to bottom surface)
                                               ! in the new vertical coordinate
      real(r8), intent(in) ::  q1(itot,jfirst-ng_s:jlast+ng_n,km,nq) ! Field input
! !INPUT/OUTPUT PARAMETERS:
      real(r8), intent(inout)::  q2(itot,jfirst-ng_s:jlast+ng_n,kn,nq) ! Field output

! !DESCRIPTION:
!
!     Perform piecewise parabolic method on a given latitude    
! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
!
! !REVISION HISTORY: 
!    02.04.04   Sawyer    incorporated latest FVGCM version, ProTeX
!    02.06.20   Sawyer    made Q2 inout since the args for Q1/Q2 same
!    03.07.22   Parks     Cleaned main loop, removed gotos
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      real(r8)       r3, r23
      parameter (r3 = D1_0/D3_0, r23 = D2_0/D3_0)
      real(r8)   dp1(i1:i2,km) 
     real(r8)  q4(4,i1:i2,km)

      integer i, k, kk, kl, k0(i1:i2,0:kn+1), iq
      real(r8)    pl, pr, qsum, qsumk(i1:i2,kn), delp, esl

      do k=1,km
         do i=i1,i2
             dp1(i,k) = pe1(i,k+1) - pe1(i,k)
         enddo
      enddo

! Mapping

! For each pe2(i,k), determine lowest pe1 interval = smallest k0 (= k0(i,k))
!   such that pe1(i,k0) <= pe2(i,k) <= pe1(i,k0+1)
!   Note that pe2(i,1)==pe1(i,1) and pe2(i,kn+1)==pe1(i,kn+1)
!   Note also that pe1, pe2 are assumed to be monotonically increasing
#if defined( UNICOSMP ) || defined ( NEC_SX )
      do kk = km, 1, -1
         do k = 1, kn+1
!dir$ prefervector
            do i = i1, i2
               if (pe2(i,k) <= pe1(i,kk+1)) then
                  k0(i,k) = kk
               endif
            enddo
         enddo
      enddo
#else
      do i = i1, i2
         k0(i,0) = 1
         do k = 1, kn+1
            do kk = k0(i,k-1), km
               if (pe2(i,k) <= pe1(i,kk+1)) then
                  k0(i,k) = kk
                  exit
               endif
            enddo
         enddo
      enddo
#endif

    do iq=1,nq

      do k=1,km
         do i=i1,i2
            q4(1,i,k) = q1(i,j,k,iq)
         enddo
      enddo

! Compute vertical subgrid distribution
      call ppm2m( q4, dp1, km, i1, i2, iv, kord )
! Interpolate
      do k = 1, kn

! Prepare contribution between pe1(i,ko(i,k)+1) and pe1(i,k0(i,k+1))
         qsumk(:,k) = D0_0
         do i = i1, i2
            do kl = k0(i,k)+1, k0(i,k+1)-1
               qsumk(i,k) = qsumk(i,k) + dp1(i,kl)*q4(1,i,kl)
            enddo
         enddo

         do i = i1, i2
            kk = k0(i,k)
! Consider contribution between pe1(i,kk) and pe2(i,k)
            pl = (pe2(i,k)-pe1(i,kk)) / dp1(i,kk)
! Check to see if pe2(i,k+1) and pe2(i,k) are in same pe1 interval
            if (k0(i,k+1) == k0(i,k)) then
               pr = (pe2(i,k+1)-pe1(i,kk)) / dp1(i,kk)
               q2(i,j,k,iq) = q4(2,i,kk) + D0_5*(q4(4,i,kk)+q4(3,i,kk)-q4(2,i,kk))  &
                  *(pr+pl) - q4(4,i,kk)*r3*(pr*(pr+pl)+pl**2)
            else
! Consider contribution between pe2(i,k) and pe1(i,kk+1)
               qsum = (pe1(i,kk+1)-pe2(i,k))*(q4(2,i,kk)+D0_5*(q4(4,i,kk)+       &
                  q4(3,i,kk)-q4(2,i,kk))*(D1_0+pl)-q4(4,i,kk)*                    &
                  (r3*(D1_0+pl*(D1_0+pl))))
! Next consider contribution between pe1(i,kk+1) and pe1(i,k0(i,k+1))
               qsum = qsum + qsumk(i,k)
! Now consider contribution between pe1(i,k0(i,k+1)) and pe2(i,k+1)
               kl = k0(i,k+1)
               delp = pe2(i,k+1)-pe1(i,kl)
               esl = delp / dp1(i,kl)
               qsum = qsum + delp*(q4(2,i,kl)+D0_5*esl*                          &
                  (q4(3,i,kl)-q4(2,i,kl)+q4(4,i,kl)*(D1_0-r23*esl)))
               q2(i,j,k,iq) = qsum / ( pe2(i,k+1) - pe2(i,k) )
            endif
         enddo
      enddo

    enddo

    return
!EOC
 end subroutine mapn_ppm
!----------------------------------------------------------------------- 


!-----------------------------------------------------------------------
!BOP
! !ROUTINE:  mapn_ppm_tracer --- Piecewise parabolic mapping, multiple tracers
!
! !INTERFACE:
 subroutine mapn_ppm_tracer(km,   pe1,   tracer, nq,                    &
                            kn,   pe2,   i1, i2, j,                     &
                            ifirst, ilast, jfirst, jlast, iv, kord)

! !USES:
      implicit none

! !INPUT PARAMETERS:
      integer, intent(in) :: i1                ! Starting longitude
      integer, intent(in) :: i2                ! Finishing longitude
      integer, intent(in) :: iv                ! Mode: 0 ==  constituents  1 == ???
      integer, intent(in) :: kord              ! Method order
      integer, intent(in) :: j                 ! Current latitude
      integer, intent(in) :: ifirst            ! Starting segment
      integer, intent(in) :: ilast             ! Finishing segment
      integer, intent(in) :: jfirst            ! Starting latitude
      integer, intent(in) :: jlast             ! Finishing latitude
      integer, intent(in) :: km                ! Original vertical dimension
      integer, intent(in) :: kn                ! Target vertical dimension
      integer, intent(in) :: nq                ! Number of tracers

      real(r8), intent(in) :: pe1(ifirst:ilast,km+1) ! pressure at layer edges 
                                               ! (from model top to bottom surface)
                                               ! in the original vertical coordinate
      real(r8), intent(in) :: pe2(ifirst:ilast,kn+1) ! pressure at layer edges 
                                               ! (from model top to bottom surface)
                                               ! in the new vertical coordinate
! !INPUT/OUTPUT PARAMETERS:
      real (r8), intent(inout)::  tracer(ifirst:ilast,jfirst:jlast,km,nq) ! Field output

! !DESCRIPTION:
!
!     Perform piecewise parabolic method on a given latitude    
! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
!
! !REVISION HISTORY: 
!    05.03.20   Sawyer    Created from mapn_ppm
!    05.04.04   Sawyer    Simplified indexing, removed ifirst
!    05.04.12   Sawyer    Added r4/r8 distinction
!    05.10.12   Worley    Made mapn_ppm_tracer vector-friendly
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

      real(r8)       r3, r23
      parameter (r3 = D1_0/D3_0, r23 = D2_0/D3_0)

      real(r8)   dp1(i1:i2,km)
      real(r8)  q4(4,i1:i2,km)

      integer i, k, kk, kl, k0(i1:i2,0:kn+1), iq
      real(r8)    pl, pr, qsum, qsumk(i1:i2,kn), delp, esl

      do k=1,km
         do i=i1,i2
             dp1(i,k) = pe1(i,k+1) - pe1(i,k)
         enddo
      enddo

! Mapping

! For each pe2(i,k), determine lowest pe1 interval = smallest k0 (= k0(i,k))
!   such that pe1(i,k0) <= pe2(i,k) <= pe1(i,k0+1)
!   Note that pe2(i,1)==pe1(i,1) and pe2(i,kn+1)==pe1(i,kn+1)
!   Note also that pe1, pe2 are assumed to be monotonically increasing
#if defined( UNICOSMP ) || defined ( NEC_SX )
      do kk = km, 1, -1
         do k = 1, kn+1
!dir$ prefervector
            do i = i1, i2
               if (pe2(i,k) <= pe1(i,kk+1)) then
                  k0(i,k) = kk
               endif
            enddo
         enddo
      enddo
#else
      do i = i1, i2
         k0(i,0) = 1
         do k = 1, kn+1
            do kk = k0(i,k-1), km
               if (pe2(i,k) <= pe1(i,kk+1)) then
                  k0(i,k) = kk
                  exit
               endif
            enddo
         enddo
      enddo
#endif

    do iq=1,nq
      do k=1,km
        do i=i1,i2
          q4(1,i,k) = tracer(i,j,k,iq)
        enddo
      enddo

! Compute vertical subgrid distribution
      call ppm2m( q4, dp1, km, i1, i2, iv, kord )

! Interpolate
      do k = 1, kn

! Prepare contribution between pe1(i,ko(i,k)+1) and pe1(i,k0(i,k+1))
         qsumk(:,k) = D0_0
         do i = i1, i2
            do kl = k0(i,k)+1, k0(i,k+1)-1
               qsumk(i,k) = qsumk(i,k) + dp1(i,kl)*q4(1,i,kl)
            enddo
         enddo

         do i = i1, i2
            kk = k0(i,k)
! Consider contribution between pe1(i,kk) and pe2(i,k)
            pl = (pe2(i,k)-pe1(i,kk)) / dp1(i,kk)
! Check to see if pe2(i,k+1) and pe2(i,k) are in same pe1 interval
            if (k0(i,k+1) == k0(i,k)) then
               pr = (pe2(i,k+1)-pe1(i,kk)) / dp1(i,kk)
               tracer(i,j,k,iq) = q4(2,i,kk)          &
                          + D0_5*(q4(4,i,kk)+q4(3,i,kk)-q4(2,i,kk))     &
                          *(pr+pl)-q4(4,i,kk)*r3*(pr*(pr+pl)+pl**2)
            else
! Consider contribution between pe2(i,k) and pe1(i,kk+1)
               qsum = (pe1(i,kk+1)-pe2(i,k))*(q4(2,i,kk)+D0_5*(q4(4,i,kk)+       &
                  q4(3,i,kk)-q4(2,i,kk))*(D1_0+pl)-q4(4,i,kk)*                    &
                  (r3*(D1_0+pl*(D1_0+pl))))
! Next consider contribution between pe1(i,kk+1) and pe1(i,k0(i,k+1))
               qsum = qsum + qsumk(i,k)
! Now consider contribution between pe1(i,k0(i,k+1)) and pe2(i,k+1)
               kl = k0(i,k+1)
               delp = pe2(i,k+1)-pe1(i,kl)
               esl = delp / dp1(i,kl)
               qsum = qsum + delp*(q4(2,i,kl)+D0_5*esl*                          &
                  (q4(3,i,kl)-q4(2,i,kl)+q4(4,i,kl)*(D1_0-r23*esl)))
               tracer(i,j,k,iq) = qsum / ( pe2(i,k+1) - pe2(i,k) )
            endif
         enddo
      enddo

    enddo   ! do iq=1,nq

    return
!EOC
 end subroutine mapn_ppm_tracer
!----------------------------------------------------------------------- 


!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  ppm2m --- Piecewise parabolic method for fields
!
! !INTERFACE:
 subroutine ppm2m(a4, delp, km, i1, i2, iv, kord)

! !USES:
 implicit none

! !INPUT PARAMETERS:
 integer, intent(in):: iv      ! iv =-1: winds
                               ! iv = 0: positive definite scalars
                               ! iv = 1: others
 integer, intent(in):: i1      ! Starting longitude
 integer, intent(in):: i2      ! Finishing longitude
 integer, intent(in):: km      ! vertical dimension
 integer, intent(in):: kord    ! Order (or more accurately method no.):
                               ! 
 real (r8), intent(in):: delp(i1:i2,km)     ! layer pressure thickness

! !INPUT/OUTPUT PARAMETERS:
 real (r8), intent(inout):: a4(4,i1:i2,km)  ! Interpolated values

! !DESCRIPTION:
!
!   Perform the piecewise parabolic method 
! 
! !REVISION HISTORY: 
!   ??.??.??    Lin        Creation
!   02.04.04    Sawyer     Newest release from FVGCM
!   02.04.23    Sawyer     Incorporated minor algorithmic change to 
!                          maintain CAM zero diffs (see comments inline)
! 
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! local arrays:
      real(r8)  dc(i1:i2,km)
      real(r8)  h2(i1:i2,km)
      real(r8) delq(i1:i2,km)
      real(r8) df2(i1:i2,km)
      real(r8) d4(i1:i2,km)

! local scalars:
      real(r8) fac
      real(r8) a1, a2, c1, c2, c3, d1, d2
      real(r8) qmax, qmin, cmax, cmin
      real(r8) qm, dq, tmp

      integer i, k, km1, lmt
      real(r8) qmp, pmp
      real(r8) lac
      integer it

      km1 = km - 1
       it = i2 - i1 + 1

      do k=2,km
         do i=i1,i2
            delq(i,k-1) =   a4(1,i,k) - a4(1,i,k-1)
              d4(i,k  ) = delp(i,k-1) + delp(i,k)
         enddo
      enddo

      do k=2,km1
         do i=i1,i2
            c1  = (delp(i,k-1)+D0_5*delp(i,k))/d4(i,k+1)
            c2  = (delp(i,k+1)+D0_5*delp(i,k))/d4(i,k)
            tmp = delp(i,k)*(c1*delq(i,k) + c2*delq(i,k-1)) /      &
                                    (d4(i,k)+delp(i,k+1))
            qmax = max(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1)) - a4(1,i,k)
            qmin = a4(1,i,k) - min(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1))
             dc(i,k) = sign(min(abs(tmp),qmax,qmin), tmp)
            df2(i,k) = tmp
         enddo
      enddo

!****6***0*********0*********0*********0*********0*********0**********72
! 4th order interpolation of the provisional cell edge value
!****6***0*********0*********0*********0*********0*********0**********72

      do k=3,km1
      do i=i1,i2
      c1 = delq(i,k-1)*delp(i,k-1) / d4(i,k)
      a1 = d4(i,k-1) / (d4(i,k) + delp(i,k-1))
      a2 = d4(i,k+1) / (d4(i,k) + delp(i,k))
      a4(2,i,k) = a4(1,i,k-1) + c1 + D2_0/(d4(i,k-1)+d4(i,k+1)) *    &
                ( delp(i,k)*(c1*(a1 - a2)+a2*dc(i,k-1)) -          &
                                delp(i,k-1)*a1*dc(i,k  ) )
      enddo
      enddo

      call steepz(i1, i2, km, a4, df2, dc, delq, delp, d4)

! Area preserving cubic with 2nd deriv. = 0 at the boundaries
! Top
      do i=i1,i2
      d1 = delp(i,1)
      d2 = delp(i,2)
      qm = (d2*a4(1,i,1)+d1*a4(1,i,2)) / (d1+d2)
      dq = D2_0*(a4(1,i,2)-a4(1,i,1)) / (d1+d2)
      c1 = D4_0*(a4(2,i,3)-qm-d2*dq) / ( d2*(D2_0*d2*d2+d1*(d2+D3_0*d1)) )
      c3 = dq - D0_5*c1*(d2*(D5_0*d1+d2)-D3_0*d1**2)
      a4(2,i,2) = qm - D0_25*c1*d1*d2*(d2+D3_0*d1)
      a4(2,i,1) = d1*(D2_0*c1*d1**2-c3) + a4(2,i,2)
      dc(i,1) =  a4(1,i,1) - a4(2,i,1)
! No over- and undershoot condition
      cmax = max(a4(1,i,1), a4(1,i,2))
      cmin = min(a4(1,i,1), a4(1,i,2))
      a4(2,i,2) = max(cmin,a4(2,i,2))
      a4(2,i,2) = min(cmax,a4(2,i,2))
      enddo

      if( iv == 0 ) then
         do i=i1,i2
!
! WS: 02.04.23  Algorithmic difference with FVGCM.  FVGCM does this:
!
!!!            a4(2,i,1) = a4(1,i,1)
!!!            a4(3,i,1) = a4(1,i,1)
!
!     CAM does this:
!
            a4(2,i,1) = max(D0_0,a4(2,i,1))
            a4(2,i,2) = max(D0_0,a4(2,i,2))
         enddo
      elseif ( iv == -1 ) then
! Winds:
        if( km > 32 ) then
          do i=i1,i2
! More dampping: top layer as the sponge
             a4(2,i,1) = a4(1,i,1)
             a4(3,i,1) = a4(1,i,1)
          enddo
        else
          do i=i1,i2
             if( a4(1,i,1)*a4(2,i,1) <=  D0_0 ) then
                 a4(2,i,1) = D0_0
             else
                 a4(2,i,1) = sign(min(abs(a4(1,i,1)),    &
                                      abs(a4(2,i,1))),   &
                                          a4(1,i,1)  )
            endif
          enddo
        endif
      endif

! Bottom
! Area preserving cubic with 2nd deriv. = 0 at the surface
      do i=i1,i2
         d1 = delp(i,km)
         d2 = delp(i,km1)
         qm = (d2*a4(1,i,km)+d1*a4(1,i,km1)) / (d1+d2)
         dq = D2_0*(a4(1,i,km1)-a4(1,i,km)) / (d1+d2)
         c1 = (a4(2,i,km1)-qm-d2*dq) / (d2*(D2_0*d2*d2+d1*(d2+D3_0*d1)))
         c3 = dq - D2_0*c1*(d2*(D5_0*d1+d2)-D3_0*d1**2)
         a4(2,i,km) = qm - c1*d1*d2*(d2+D3_0*d1)
         a4(3,i,km) = d1*(D8_0*c1*d1**2-c3) + a4(2,i,km)
         dc(i,km) = a4(3,i,km) -  a4(1,i,km)
! No over- and under-shoot condition
         cmax = max(a4(1,i,km), a4(1,i,km1))
         cmin = min(a4(1,i,km), a4(1,i,km1))
         a4(2,i,km) = max(cmin,a4(2,i,km))
         a4(2,i,km) = min(cmax,a4(2,i,km))
      enddo

! Enforce constraint at the surface

      if ( iv == 0 ) then
! Positive definite scalars:
           do i=i1,i2
              a4(3,i,km) = max(D0_0, a4(3,i,km))
           enddo
      elseif ( iv == -1 ) then
! Winds:
           do i=i1,i2
              if( a4(1,i,km)*a4(3,i,km) <=  D0_0 ) then
                  a4(3,i,km) = D0_0
              else
                  a4(3,i,km) = sign( min(abs(a4(1,i,km)),   &
                                         abs(a4(3,i,km))),  &
                                             a4(1,i,km)  )
              endif
           enddo
      endif

      do k=1,km1
         do i=i1,i2
            a4(3,i,k) = a4(2,i,k+1)
         enddo
      enddo
 
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
 
! Top 2 and bottom 2 layers always use monotonic mapping
      do k=1,2
         do i=i1,i2
            a4(4,i,k) = D3_0*(D2_0*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call kmppm(dc(i1,k), a4(1,i1,k), it, 0)
      enddo

      if(kord >= 7) then
!****6***0*********0*********0*********0*********0*********0**********72
! Huynh's 2nd constraint
!****6***0*********0*********0*********0*********0*********0**********72
      do k=2, km1
         do i=i1,i2
! Method#1
!           h2(i,k) = delq(i,k) - delq(i,k-1)
! Method#2
!           h2(i,k) = D2_0*(dc(i,k+1)/delp(i,k+1) - dc(i,k-1)/delp(i,k-1))
!    &               / ( delp(i,k)+D0_5*(delp(i,k-1)+delp(i,k+1)) )
!    &               * delp(i,k)**2
! Method#3
            h2(i,k) = dc(i,k+1) - dc(i,k-1)
         enddo
      enddo

      if( kord == 7 ) then
         fac = D1_5           ! original quasi-monotone
      else
         fac = D0_125         ! full monotone
      endif

      do k=3, km-2
        do i=i1,i2
! Right edges
!        qmp   = a4(1,i,k) + D2_0*delq(i,k-1)
!        lac   = a4(1,i,k) + fac*h2(i,k-1) + D0_5*delq(i,k-1)
!
         pmp   = D2_0*dc(i,k)
         qmp   = a4(1,i,k) + pmp
         lac   = a4(1,i,k) + fac*h2(i,k-1) + dc(i,k)
         qmin  = min(a4(1,i,k), qmp, lac)
         qmax  = max(a4(1,i,k), qmp, lac)
         a4(3,i,k) = min(max(a4(3,i,k), qmin), qmax)
! Left  edges
!        qmp   = a4(1,i,k) - D2_0*delq(i,k)
!        lac   = a4(1,i,k) + fac*h2(i,k+1) - D0_5*delq(i,k)
!
         qmp   = a4(1,i,k) - pmp
         lac   = a4(1,i,k) + fac*h2(i,k+1) - dc(i,k)
         qmin  = min(a4(1,i,k), qmp, lac)
         qmax  = max(a4(1,i,k), qmp, lac)
         a4(2,i,k) = min(max(a4(2,i,k), qmin), qmax)
! Recompute A6
         a4(4,i,k) = D3_0*(D2_0*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
        enddo
! Additional constraint to prevent negatives when kord=7
         if (iv == 0 .and. kord == 7) then
             call kmppm(dc(i1,k), a4(1,i1,k), it, 2)
         endif
      enddo

      else
 
         lmt = kord - 3
         lmt = max(0, lmt)
         if (iv == 0) lmt = min(2, lmt)

      do k=3, km-2
      if( kord /= 4) then
         do i=i1,i2
            a4(4,i,k) = D3_0*(D2_0*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
      endif
         call kmppm(dc(i1,k), a4(1,i1,k), it, lmt)
      enddo
      endif

      do k=km1,km
         do i=i1,i2
            a4(4,i,k) = D3_0*(D2_0*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call kmppm(dc(i1,k), a4(1,i1,k), it, 0)
      enddo

      return
!EOC
 end subroutine ppm2m
!-----------------------------------------------------------------------


!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  ppme --- PPM scheme at vertical edges
!
! !INTERFACE:
      subroutine ppme(p,qe,delp,im,km)
! !USES:
      use shr_kind_mod, only : r8 => shr_kind_r8, i4 => shr_kind_i4 
      implicit none

! !INPUT PARAMETERS:
      integer,  intent(in)     ::  im, km
      real(r8), intent(in)     ::  p(im,km), delp(im,km)

! !INPUT/OUTPUT PARAMETERS:
      real(r8), intent(out)    ::  qe(im,km+1)

! !DESCRIPTION:
!
!
! !REVISION HISTORY:
!    05.06.13   Sawyer    Inserted file ppme.F90 here, added ProTeX
!
!EOP
!-----------------------------------------------------------------------
!BOC

      integer(i4)  km1
      integer(i4)  i, k
! local arrays.
      real(r8) dc(im,km),delq(im,km), a6(im,km)
      real(r8) c1, c2, c3, tmp, qmax, qmin
      real(r8) a1, a2, s1, s2, s3, s4, ss3, s32, s34, s42
      real(r8) a3, b2, sc, dm, d1, d2, f1, f2, f3, f4
      real(r8) qm, dq

      km1 = km - 1

      do 500 k=2,km
      do 500 i=1,im
500   a6(i,k) = delp(i,k-1) + delp(i,k)

      do 1000 k=1,km1
      do 1000 i=1,im
      delq(i,k) = p(i,k+1) - p(i,k)
1000  continue

      do 1220 k=2,km1
      do 1220 i=1,im
      c1 = (delp(i,k-1)+D0_5*delp(i,k))/a6(i,k+1)
      c2 = (delp(i,k+1)+D0_5*delp(i,k))/a6(i,k)
      tmp = delp(i,k)*(c1*delq(i,k) + c2*delq(i,k-1)) /    &
                                    (a6(i,k)+delp(i,k+1))
      qmax = max(p(i,k-1),p(i,k),p(i,k+1)) - p(i,k)
      qmin = p(i,k) - min(p(i,k-1),p(i,k),p(i,k+1))
      dc(i,k) = sign(min(abs(tmp),qmax,qmin), tmp)
1220  continue

!****6***0*********0*********0*********0*********0*********0**********72
! 4th order interpolation of the provisional cell edge value
!****6***0*********0*********0*********0*********0*********0**********72

      do 12 k=3,km1
      do 12 i=1,im
      c1 = delq(i,k-1)*delp(i,k-1) / a6(i,k)
      a1 = a6(i,k-1) / (a6(i,k) + delp(i,k-1))
      a2 = a6(i,k+1) / (a6(i,k) + delp(i,k))
      qe(i,k) = p(i,k-1) + c1 + D2_0/(a6(i,k-1)+a6(i,k+1)) *        &
                ( delp(i,k)*(c1*(a1 - a2)+a2*dc(i,k-1)) -         &
                                delp(i,k-1)*a1*dc(i,k  ) )
12    continue

! three-cell parabolic subgrid distribution at model top

      do 10 i=1,im
! three-cell PP-distribution
! Compute a,b, and c of q = aP**2 + bP + c using cell averages and delp
! a3 = a / 3
! b2 = b / 2
      s1 = delp(i,1)
      s2 = delp(i,2) + s1
!
      s3 = delp(i,2) + delp(i,3)
      s4 = s3 + delp(i,4)
      ss3 =  s3 + s1
      s32 = s3*s3
      s42 = s4*s4
      s34 = s3*s4
! model top
      a3 = (delq(i,2) - delq(i,1)*s3/s2) / (s3*ss3)
!
      if(abs(a3) .gt. D1EM14) then
         b2 =  delq(i,1)/s2 - a3*(s1+s2)
         sc = -b2/(D3_0*a3)
         if(sc .lt. D0_0 .or. sc .gt. s1) then
             qe(i,1) = p(i,1) - s1*(a3*s1 + b2)
         else
             qe(i,1) = p(i,1) - delq(i,1)*s1/s2
         endif
      else
! Linear
         qe(i,1) = p(i,1) - delq(i,1)*s1/s2
      endif
      dc(i,1) = p(i,1) - qe(i,1)
! compute coef. for the off-centered area preserving cubic poly.
      dm = delp(i,1) / (s34*ss3*(delp(i,2)+s3)*(s4+delp(i,1)))
      f1 = delp(i,2)*s34 / ( s2*ss3*(s4+delp(i,1)) )
      f2 = (delp(i,2)+s3) * (ss3*(delp(i,2)*s3+s34+delp(i,2)*s4)   &
            + s42*(delp(i,2)+s3+s32/s2))
      f3 = -delp(i,2)*( ss3*(s32*(s3+s4)/(s4-delp(i,2))            &
            + (delp(i,2)*s3+s34+delp(i,2)*s4))                     &
            + s42*(delp(i,2)+s3) )
      f4 = ss3*delp(i,2)*s32*(delp(i,2)+s3) / (s4-delp(i,2))
      qe(i,2) = f1*p(i,1)+(f2*p(i,2)+f3*p(i,3)+f4*p(i,4))*dm
10    continue

! Bottom
! Area preserving cubic with 2nd deriv. = 0 at the surface
      do 15 i=1,im
      d1 = delp(i,km)
      d2 = delp(i,km1)
      qm = (d2*p(i,km)+d1*p(i,km1)) / (d1+d2)
      dq = D2_0*(p(i,km1)-p(i,km)) / (d1+d2)
      c1 = (qe(i,km1)-qm-d2*dq) / (d2*(D2_0*d2*d2+d1*(d2+D3_0*d1)))
      c3 = dq - D2_0*c1*(d2*(D5_0*d1+d2)-D3_0*d1**2)
      qe(i,km  ) = qm - c1*d1*d2*(d2+D3_0*d1)
      qe(i,km+1) = d1*(D8_0*c1*d1**2-c3) + qe(i,km)
15    continue
      return
!EOC
      end subroutine ppme
!----------------------------------------------------------------------- 

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  kmppm --- Perform piecewise parabolic method in vertical
!
! !INTERFACE:
 subroutine kmppm(dm, a4, itot, lmt)

! !USES:
 implicit none

! !INPUT PARAMETERS:
      real(r8), intent(in)::     dm(*)     ! ??????
      integer, intent(in) ::     itot      ! Total Longitudes
      integer, intent(in) ::     lmt       ! 0: Standard PPM constraint
                                           ! 1: Improved full monotonicity constraint (Lin)
                                           ! 2: Positive definite constraint
                                           ! 3: do nothing (return immediately)

! !INPUT/OUTPUT PARAMETERS:
      real(r8), intent(inout) :: a4(4,*)   ! ???????
                                           ! AA <-- a4(1,i)
                                           ! AL <-- a4(2,i)
                                           ! AR <-- a4(3,i)
                                           ! A6 <-- a4(4,i)

! !DESCRIPTION:
!
!    Writes a standard set of data to the history buffer. 
!
! !REVISION HISTORY: 
!    00.04.24   Lin       Last modification
!    01.03.26   Sawyer    Added ProTeX documentation
!    02.04.04   Sawyer    Incorporated newest FVGCM version
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

      real(r8)       r12
      parameter (r12 = D1_0/D12_0)

      real(r8) qmp
      integer i
      real(r8) da1, da2, a6da
      real(r8) fmin

! Developer: S.-J. Lin, NASA-GSFC
! Last modified: Apr 24, 2000

      if ( lmt == 3 ) return

      if(lmt == 0) then
! Standard PPM constraint
      do i=1,itot
      if(dm(i) == D0_0) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = D0_0
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da < -da2) then
            a4(4,i) = D3_0*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da > da2) then
            a4(4,i) = D3_0*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
      enddo

      elseif (lmt == 1) then

! Improved full monotonicity constraint (Lin)
! Note: no need to provide first guess of A6 <-- a4(4,i)
      do i=1, itot
           qmp = D2_0*dm(i)
         a4(2,i) = a4(1,i)-sign(min(abs(qmp),abs(a4(2,i)-a4(1,i))), qmp)
         a4(3,i) = a4(1,i)+sign(min(abs(qmp),abs(a4(3,i)-a4(1,i))), qmp)
         a4(4,i) = D3_0*( D2_0*a4(1,i) - (a4(2,i)+a4(3,i)) )
      enddo

      elseif (lmt == 2) then

! Positive definite constraint
      do i=1,itot
      if( abs(a4(3,i)-a4(2,i)) < -a4(4,i) ) then
      fmin = a4(1,i)+D0_25*(a4(3,i)-a4(2,i))**2/a4(4,i)+a4(4,i)*r12
         if( fmin < D0_0 ) then
         if(a4(1,i)<a4(3,i) .and. a4(1,i)<a4(2,i)) then
            a4(3,i) = a4(1,i)
            a4(2,i) = a4(1,i)
            a4(4,i) = D0_0
         elseif(a4(3,i) > a4(2,i)) then
            a4(4,i) = D3_0*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         else
            a4(4,i) = D3_0*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
         endif
      endif
      enddo

      endif

      return
!EOC
 end subroutine kmppm
!-----------------------------------------------------------------------

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  steepz --- Calculate attributes for PPM
!
! !INTERFACE:
 subroutine steepz(i1, i2, km, a4, df2, dm, dq, dp, d4)

! !USES:
   implicit none

! !INPUT PARAMETERS:
      integer, intent(in) :: km                   ! Total levels
      integer, intent(in) :: i1                   ! Starting longitude
      integer, intent(in) :: i2                   ! Finishing longitude
      real(r8), intent(in) ::  dp(i1:i2,km)       ! grid size
      real(r8), intent(in) ::  dq(i1:i2,km)       ! backward diff of q
      real(r8), intent(in) ::  d4(i1:i2,km)       ! backward sum:  dp(k)+ dp(k-1) 
      real(r8), intent(in) :: df2(i1:i2,km)       ! first guess mismatch
      real(r8), intent(in) ::  dm(i1:i2,km)       ! monotonic mismatch

! !INPUT/OUTPUT PARAMETERS:
      real(r8), intent(inout) ::  a4(4,i1:i2,km)  ! first guess/steepened

!
! !DESCRIPTION:
!   This is complicated stuff related to the Piecewise Parabolic Method
!   and I need to read the Collela/Woodward paper before documenting
!   thoroughly.
!
! !REVISION HISTORY: 
!   ??.??.??    Lin?       Creation
!   01.03.26    Sawyer     Added ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, k
      real(r8) alfa(i1:i2,km)
      real(r8)    f(i1:i2,km)
      real(r8)  rat(i1:i2,km)
      real(r8)  dg2

! Compute ratio of dq/dp
      do k=2,km
         do i=i1,i2
            rat(i,k) = dq(i,k-1) / d4(i,k)
         enddo
      enddo

! Compute F
      do k=2,km-1
         do i=i1,i2
            f(i,k) = (rat(i,k+1) - rat(i,k))                             &
                     / ( dp(i,k-1)+dp(i,k)+dp(i,k+1) )
         enddo
      enddo

      do k=3,km-2
         do i=i1,i2
         if(f(i,k+1)*f(i,k-1)<D0_0 .and. df2(i,k)/=D0_0) then
            dg2 = (f(i,k+1)-f(i,k-1))*((dp(i,k+1)-dp(i,k-1))**2          &
                   + d4(i,k)*d4(i,k+1) )
            alfa(i,k) = max(D0_0, min(D0_5, -D0_1875*dg2/df2(i,k))) 
         else
            alfa(i,k) = D0_0
         endif
         enddo
      enddo

      do k=4,km-2
         do i=i1,i2
            a4(2,i,k) = (D1_0-alfa(i,k-1)-alfa(i,k)) * a4(2,i,k) +         &
                        alfa(i,k-1)*(a4(1,i,k)-dm(i,k))    +             &
                        alfa(i,k)*(a4(1,i,k-1)+dm(i,k-1))
         enddo
      enddo

      return
!EOC
 end subroutine steepz
!----------------------------------------------------------------------- 

end module mapz_module
