module pft_module
!BOP
!
! !MODULE: pft_module --- polar filters
!
! !USES:

 use shr_kind_mod,   only: r8 => shr_kind_r8
 use fv_control_mod, only: fft_flt

#ifdef NO_R16
   integer,parameter :: r16= selected_real_kind(12) ! 8 byte real
#else
   integer,parameter :: r16= selected_real_kind(24) ! 16 byte real
#endif

!
! !PUBLIC MEMBER FUNCTIONS:
      public pft2d, pft_cf, fftfax, pftinit, fftrans
!
! !DESCRIPTION:
!
!      This module provides fast-Fourier transforms
!
!      \begin{tabular}{|l|l|} \hline \hline
!         pftinit   &  \\ \hline
!         pft2d     &  \\ \hline
!         pft\_cf   &  \\ \hline
!         fftfax    &  \\ \hline
!         fftrans   &  \\ \hline
!                                \hline
!      \end{tabular}
!
! !REVISION HISTORY:
!   01.01.30   Lin        Integrated into this module
!   01.03.26   Sawyer     Added ProTeX documentation
!   05.05.25   Sawyer     Merged CAM and GEOS5 versions (CAM vectorization)
!   05.07.26   Worley     Revised module using for Cray X1 version
!   06.09.08   Sawyer     Magic numbers isolated in F90 parameters
!
!EOP
!-----------------------------------------------------------------------
      private
      real(r8), parameter ::  D0_0                    =  0.0_r8
      real(r8), parameter ::  D1EM20                  =  1.0e-20_r8
      real(r8), parameter ::  D0_5                    =  0.5_r8
      real(r8), parameter ::  D1_0                    =  1.0_r8
      real(r8), parameter ::  D1_01                   =  1.01_r8
      real(r8), parameter ::  D2_0                    =  2.0_r8
      real(r8), parameter ::  D4_0                    =  4.0_r8
      real(r8), parameter ::  D8_0                    =  8.0_r8
      real(r8), parameter ::  D180_0                  =180.0_r8

      integer, save :: ifax(13)                      !ECMWF fft
      real(r8), allocatable, save :: trigs(:)        ! reentrant code??

CONTAINS

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: pftinit --- Two-dimensional FFT initialization
!
! !INTERFACE:
 subroutine pftinit(im)

! !USES:
 implicit none

! !INPUT PARAMETERS:
      integer im                   ! Total X dimension

! !DESCRIPTION:
!
!   Perform a two-dimensional FFT initialization
!
! !REVISION HISTORY:
!   05.05.15   Mirin        Put into this module
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer icffta
      real(r8) rcffta

#if defined( LIBSCI_FFT )
          allocate( trigs(2*im+100) )
          icffta = 0
          rcffta = D0_0
          call dzfftm(0, im, icffta, rcffta, rcffta, icffta,     &
               rcffta, icffta, trigs, rcffta, icffta)
#else
          allocate( trigs(3*im/2+1) )
          call fftfax(im, ifax, trigs)
#endif

      return
!EOC
 end subroutine pftinit
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: pft2d --- Two-dimensional fast Fourier transform
!
! !INTERFACE: 
 subroutine pft2d(p, s, damp, im,  jp, q1, q2)

! !USES:
 implicit none

! !INPUT PARAMETERS:
      integer im                   ! Total X dimension
      integer jp                   ! Total Y dimension
      real(r8)   s(jp)             ! 3-point algebraic filter
      real(r8)  damp(im,jp)        ! FFT damping coefficients

! !INPUT/OUTPUT PARAMETERS:
      real(r8) q1( im+2, *)        ! Work array
      real(r8) q2(*)               ! Work array
      real(r8)  p(im,jp)           ! Array to be polar filtered

! !DESCRIPTION:
!
!   Perform a two-dimensional fast Fourier transformation.
!
! !REVISION HISTORY:
!   01.01.30   Lin          Put into this module
!   01.03.26   Sawyer       Added ProTeX documentation
!   02.04.05   Sawyer       Integrated newest FVGCM version
!   05.05.17   Sawyer       Merged CAM and GEOS-5
!   05.07.26   Worley       Removed ifax, trigs from arg list
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      real(r8) rsc, bt 
      integer  i, j, n, nj

!Local Auto arrays:
      real(r8) ptmp(0:im+1)
!!!      real(r8) q1(  im+2, jp)
!!!      real(r8) q2( (im+1)*jp )
      integer  jf(jp)

      nj = 0

      do 200 j=1,jp

      if(s(j) > D1_01 ) then
       if(fft_flt .eq. 0 .and. s(j) <= D4_0) then

         rsc = D1_0/s(j)
         bt  = D0_5*(s(j)-D1_0)

         do i=1,im
            ptmp(i) = p(i,j)
         enddo
           ptmp(   0) = p(im,j)
           ptmp(im+1) = p(1 ,j)

         do i=1,im
            p(i,j) = rsc * ( ptmp(i) + bt*(ptmp(i-1)+ptmp(i+1)) )
         enddo

       else

! Packing for FFT 
         nj  = nj + 1
         jf(nj) = j

         do i=1,im
            q1(i,nj) = p(i,j)
         enddo
            q1(im+1,nj) = D0_0
            q1(im+2,nj) = D0_0

       endif
      endif
200   continue

      if( nj == 0) return

      call fftrans(damp, im, jp, nj, jf, q1, q2)

      do n=1,nj
         do i=1,im
            p(i,jf(n)) = q1(i,n)
         enddo
      enddo

      return
!EOC
 end subroutine pft2d
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: fftrans --- Two-dimensional fast Fourier transform
!
! !INTERFACE:
 subroutine fftrans(damp, im, jp, nj, jf, q1, q2)

! !USES:
 implicit none

! !INPUT PARAMETERS:
      integer im                   ! Total X dimension
      integer jp                   ! Total Y dimension
      integer nj                   ! Number of transforms
      integer jf(jp)               ! J index versus transform number
      real(r8)  damp(im,jp)        ! FFT damping coefficients

! !INPUT/OUTPUT PARAMETERS:
      real(r8) q1( im+2, *)        ! Work array
      real(r8) q2(*)               ! Work array

! !DESCRIPTION:
!
!   Perform a two-dimensional fast Fourier transformation.
!
! !REVISION HISTORY:
!   05.05.15   Mirin        Initial combined version
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, n
      real (r8) ooim

!Local Auto arrays:

#if defined( LIBSCI_FFT )
      real (r8) qwk(2*im+4, jp)
      complex(r8) cqf(im/2+1, jp)
      integer imo2p
#elif defined( SGI_FFT )
      integer*4 im_4, nj_4, imp2_4
#endif

#if defined( LIBSCI_FFT )
      imo2p = im/2 + 1
      ooim = D1_0/real(im,r8)

      call dzfftm(-1, im, nj, D1_0, q1, im+2, cqf, imo2p,      &
                  trigs, qwk, 0)

      do n=1,nj
         do i=3,imo2p
            cqf(i,n) = cqf(i,n) * damp(2*i-2,jf(n))
         enddo
      enddo

      call zdfftm( 1, im, nj, ooim, cqf, imo2p, q1, im+2,    &
                  trigs, qwk, 0)
#elif defined( SGI_FFT )
      im_4 = im
      nj_4 = nj
      imp2_4 = im+2
      call dzfftm1du (-1, im_4, nj_4, q1, 1, imp2_4, trigs)
      do n=1,nj
         do i=5,im+2
            q1(i,n) = q1(i,n) * damp(i-2,jf(n))
         enddo
      enddo
      call dzfftm1du (1, im_4, nj_4, q1, 1, imp2_4, trigs)
      ooim = D1_0/real(im,r8)
      do n=1,nj
        do i=1,im+2
          q1(i,n) = ooim*q1(i,n)
        enddo
      enddo
#else
      call fft991 (q1, q2, trigs, ifax, 1, im+2, im, nj, -1)
      do n=1,nj
         do i=5,im+2
            q1(i,n) = q1(i,n) * damp(i-2,jf(n))
         enddo
      enddo
      call fft991 (q1, q2, trigs, ifax, 1, im+2, im, nj, 1)
#endif

      return
!EOC
 end subroutine fftrans
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: pft_cf --- Calculate algebraic and FFT polar filters
!
! !INTERFACE: 
 subroutine pft_cf(im, jm, js2g0, jn2g0, jn1g1, sc, se, dc, de,          &
                   cosp, cose, ycrit)

! !USES:
 implicit none

! !INPUT PARAMETERS:
      integer im                      ! Total X dimension
      integer jm                      ! Total Y dimension
      integer js2g0                   ! j south limit ghosted 0 (SP: from 2)
      integer jn2g0                   ! j north limit ghosted 0 (NP: from jm-1)
      integer jn1g1                   ! j north limit ghosted 1 (starts jm)
      real (r8)   cosp(jm)            ! cosine array
      real (r8)   cose(jm)            ! cosine array
      real (r8)   ycrit               ! critical value

! !OUTPUT PARAMETERS:
      real (r8)   sc(js2g0:jn2g0)     ! Algebric filter at center
      real (r8)   se(js2g0:jn1g1)     ! Algebric filter at edge
      real (r8)   dc(im,js2g0:jn2g0)  ! FFT filter at center
      real (r8)   de(im,js2g0:jn1g1)  ! FFT filter at edge

! !DESCRIPTION:
!
!   Compute coefficients for the 3-point algebraic and the FFT
!   polar filters.
!
! !REVISION HISTORY:
!
!   99.01.01   Lin          Creation
!   99.08.20   Sawyer/Lin   Changes for SPMD mode
!   01.01.30   Lin          Put into this module
!   01.03.26   Sawyer       Added ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      real (r8), parameter ::  pi = 3.14159265358979323846_R8
      integer i, j
      real (r8)   dl, coszc, cutoff, phi, damp

      coszc = cos(ycrit*pi/D180_0)

! INIT fft polar coefficients:
      dl = pi/real(im,r8)
      cutoff = D1EM20

      do j=js2g0,jn2g0
         do i=1,im
            dc(i,j) = D1_0
         enddo
      enddo

      do j=js2g0,jn1g1
         do i=1,im
            de(i,j) = D1_0
         enddo
      enddo

!     write(iulog,*) '3-point polar filter coefficients:'

!************
! Cell center
!************
      do j=js2g0,jn2g0
            sc(j) = (coszc/cosp(j))**2

         if(sc(j) > D1_0 ) then
          if(fft_flt .eq. 0 .and. sc(j) <= D2_0) then
            sc(j) =  D1_0 +  (sc(j)-D1_0)/(sc(j)+D1_0)
          elseif(fft_flt .eq. 0 .and. sc(j) <= D4_0) then
            sc(j) =  D1_0 +  sc(j)/(D8_0-sc(j))
          else

! FFT filter
            do i=1,im/2
               phi = dl * i
               damp = min((cosp(j)/coszc)/sin(phi),D1_0)**2
               if(damp < cutoff) damp = D0_0
               dc(2*i-1,j) = damp
               dc(2*i  ,j) = damp
            enddo

          endif
         endif
      enddo

!************
! Cell edges
!************
      do j=js2g0,jn1g1
            se(j) = (coszc/cose(j))**2

         if(se(j) > D1_0 ) then
          if(fft_flt .eq. 0 .and. se(j) <= D2_0) then
            se(j) =  D1_0 +  (se(j)-D1_0)/(se(j)+D1_0)
          elseif(fft_flt .eq. 0 .and. se(j) <= D4_0) then
            se(j) =  D1_0 +  se(j)/(D8_0-se(j))
          else
! FFT
            do i=1,im/2
               phi = dl * i
               damp = min((cose(j)/coszc)/sin(phi), D1_0)**2
               if(damp < cutoff) damp = D0_0
               de(2*i-1,j) = damp
               de(2*i  ,j) = damp
            enddo
          endif
         endif
      enddo
      return
!EOC
 end subroutine pft_cf
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: fftfax --- Initialize FFT
!
! !INTERFACE: 
 subroutine fftfax (n, ifaxx, trigss)

! !USES:
      implicit none

! !DESCRIPTION:
!
!   Initialize the fast Fourier transform.  If CPP token SGI_FFT is
!   set, SGI libraries will be used.  Otherwise the Fortran code
!   is inlined.
!
! !REVISION HISTORY:
!
!   99.11.24   Sawyer       Added wrappers for SGI
!   01.03.26   Sawyer       Added ProTeX documentation
!   05.07.26   Worley       Modified version for Cray X1
!
!EOP
!-----------------------------------------------------------------------
!BOC

      integer n

#if defined( SGI_FFT )
      real(r8)    trigss(1)
      integer ifaxx(*)
! local
      integer*4 nn

      nn=n
      call dzfftm1dui (nn,trigss)
#else
      integer ifaxx(13)
      real(r8) trigss(3*n/2+1)
      call set99(trigss,ifaxx,n)
#endif
      return
!EOC
 end subroutine fftfax
!-----------------------------------------------------------------------

end module pft_module
