
      module mo_xsections

      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

      private
      public :: r08_inti, r44_inti
      public :: r01, r04, r06, r08, r10, r11
      public :: r14, r15, r17, r18, r44, xs_mvk

      save

      real(r8), allocatable :: a(:)
      real(r8), allocatable :: b(:)
      real(r8), allocatable :: suma(:)
      real(r8), allocatable :: sumb(:)

      contains

      subroutine r01( nw, wl, wc, tlev, airlev, jlabel, xs )
!-----------------------------------------------------------------------------
!   purpose:
!   provide the product of (cross section) x (quantum yield) for the two
!   o3 photolysis reactions:
!              (a) o3 + hv -> o2 + o(1d)
!              (b) o3 + hv -> o2 + o(3p)
!   cross section:  combined data from wmo 85 ozone assessment (use 273k
!                   value from 175.439-847.5 nm) and data from molina and
!                   molina (use in hartley and huggins bans (240.5-350 nm)
!   quantum yield:  choice between
!                    (1) data from michelsen et al, 1994
!                    (2) jpl 87 recommendation
!                    (3) jpl 90/92 recommendation (no "tail")
!                    (4) data from shetter et al., 1996
!                    (5) jpl 97 recommendation
!-----------------------------------------------------------------------------
!   parameters:
!   nw     - integer, number of specified intervals + 1 in working        (i)
!            wavelength grid
!   wl     - real(r8), vector of lower limits of wavelength intervals in      (i)
!            working wavelength grid
!   wc     - real(r8), vector of center points of wavelength intervals in     (i)
!            working wavelength grid
!   nz     - integer, number of altitude levels in working altitude grid  (i)
!   tlev   - real(r8), temperature (k) at each specified altitude level       (i)
!   airlev - real(r8), air density (molec/cc) at each altitude level          (i)
!   j      - integer, counter for number of weighting functions defined  (io)
!   sq     - real(r8), cross section x quantum yield (cm^2) for each          (o)
!            photolysis reaction defined, at each defined wavelength and
!            at each defined altitude level
!   jlabel - character*40, string identifier for each photolysis reaction (o)
!            defined
!-----------------------------------------------------------------------------
!   edit history:
!   05/98  original, adapted from former jspec1 subroutine
!-----------------------------------------------------------------------------
!  this program is free software;  you can redistribute it and/or modify
!  it under the terms of the gnu general public license as published by the
!  free software foundation;  either version 2 of the license, or (at your
!  option) any later version.
!  the tuv package is distributed in the hope that it will be useful, but
!  without any warranty;  without even the implied warranty of merchantibi-
!  lity or fitness for a particular purpose.  see the gnu general public
!  license for more details.
!  free software foundation, inc., 675 mass ave, cambridge, ma 02139, usa.
!-----------------------------------------------------------------------------

      use mo_params,    only : kw
      use ppgrid,       only : pverp
      use mo_waveo3
      use mo_waveall,   only : r01g1, r01g2

      implicit none

!-----------------------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------------------
      integer, intent(in)          :: nw
      real(r8), intent(in)         :: wl(kw)
      real(r8), intent(in)         :: wc(kw)
      real(r8), intent(in)         :: tlev(pverp)
      real(r8), intent(in)         :: airlev(pverp)
      real(r8), intent(inout)      :: xs(:,:)
      character(len=*), intent(in) :: jlabel

!-----------------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------------
      real(r8), parameter ::     c0  = 12._r8/19._r8
      real(r8), parameter ::     a1  = 0.887_r8
      real(r8), parameter ::     a2  = 2.35_r8
      real(r8), parameter ::     a3  = 57.0_r8
      real(r8), parameter ::     wc1 = 302._r8
      real(r8), parameter ::     wc2 = 311.1_r8
      real(r8), parameter ::     wc3 = 313.9_r8
      real(r8), parameter ::     v2  = 820.0_r8
      real(r8), parameter ::     v3  = 1190.0_r8
      real(r8), parameter ::     w1  = 7.9_r8
      real(r8), parameter ::     w2  = 2.2_r8
      real(r8), parameter ::     w3  = 7.4_r8
      real(r8), parameter ::     xk  = 0.695_r8

      integer  :: i, wn
      integer  :: myld
      real(r8) :: qy1d, qy3p, so3, wrk, tinv
      logical  :: tl

!-----------------------------------------------------------------------------
!       NEW O3 qy 2002 
!       Matsumi et al., 2002
!-----------------------------------------------------------------------------
      integer, parameter :: kmats = 7

      myld  = kmats

level_loop : &
      do i = 1,pverp
         tinv = 1._r8/tlev(i)
         tl   = tlev(i) < 263._r8
         if( tl ) then
            wrk = (tlev(i) - 226._r8)/(263._r8 - 226._r8)
         else
            wrk = (tlev(i) - 263._r8)/(298._r8 - 263._r8)
         end if
wave_loop : &
         do wn = 1,nw
            if( wl(wn) > 240.5_r8  .and. wl(wn+1) < 350._r8 ) then
               if( tl ) then
                  so3 = s226(wn) + (s263(wn) - s226(wn)) * wrk
               else
                  so3 = s263(wn) + (s298(wn) - s263(wn)) * wrk
               end if
            else
               so3 = xso3(wn)
            end if
!-----------------------------------------------------------------------------
! 	... from jpl97
!-----------------------------------------------------------------------------
             if( wc(wn) < 271._r8 ) then
                qy1d = .87_r8
             else if( wc(wn) >= 271._r8 .and. wc(wn) < 290._r8 ) then
                qy1d = .87_r8 + (wc(wn) - 271._r8)*c0
             else if( wc(wn) >= 290._r8 .and. wc(wn) < 305._r8 ) then
                qy1d = .95_r8
             else if( wc(wn) >= 305._r8 .and. wc(wn) <= 325._r8 ) then
                qy1d = r01g1(wn) * exp ( -r01g2(wn)*tinv )
             else
                qy1d = 0._r8
             end if
!-------------------------------------------------------------------------------
!	... from jpl2000
!-------------------------------------------------------------------------------
             if( wc(wn) < 300._r8 ) then
                qy1d = 0.95_r8
             else if( wc(wn) >= 300._r8 .and. wc(wn) < 331._r8 ) then
                qy1d = a1*exp( -((wc(wn) - wc1 )/w1)**4 ) &
                     + a2*(tlev(i)/300._r8)**4*exp( -v2/xk*tinv ) &
                    * exp( -((wc(wn) - wc2)/w2)**2 ) &
                     + a3*exp( -v3/xk*tinv ) * exp( -((wc(wn) - wc3)/w3)**2 ) &
                     + 0.06_r8
             else if( wc(wn) >= 331._r8 .and. wc(wn) <= 345._r8 ) then
                qy1d = 0.06_r8
             else
                qy1d = 0._r8
             end if

             if( myld == kmats ) then
                qy1d = fo3qy( wc(wn), tlev(i) )
             end if

             if( (trim( jlabel ) == 'jo1d') .or. (trim( jlabel ) == 'j2oh') ) then
                xs(wn,i) = qy1d*so3
             else
                qy3p     = 1._r8 - qy1d
                xs(wn,i) = qy3p*so3
             end if
         end do wave_loop
      end do level_loop

      end subroutine r01

      function fo3qy( w, t )
!-----------------------------------------------------------------------------
!   PURPOSE:
! function to calculate the quantum yield O3 + hv -> O(1D) + O2,
! according to:                                                             
! Matsumi, Y., F. J. Comes, G. Hancock, A. Hofzumanhays, A. J. Hynes,
! M. Kawasaki, and A. R. Ravishankara, QUantum yields for production of O(1D)
! in the ultraviolet photolysis of ozone:  Recommendation based on evaluation
! of laboratory data, J. Geophys. Res., 107, 10.1029/2001JD000510, 2002.
!-----------------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------------------
      real(r8), intent(in)  :: w
      real(r8), intent(in)  :: t

!-----------------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------------
      real(r8) :: kt
      real(r8) :: q1
      real(r8) :: q2 
      real(r8) :: a(3)  = (/ 0.8036_r8, 8.9061_r8, 0.1192_r8 /)
      real(r8) :: x(3)  = (/ 304.225_r8, 314.957_r8, 310.737_r8 /)
      real(r8) :: om(3) = (/ 5.576_r8, 6.601_r8, 2.187_r8 /)

!-----------------------------------------------------------------------------
!	... function declarations
!-----------------------------------------------------------------------------
      real(r8) :: fo3qy
      
      fo3qy = 0._r8
      kt = 0.695_r8 * t
      q1 = 1._r8
      q2 = exp( -825.518_r8/kt )
      
      if( w <= 305._r8 ) then
         fo3qy = .90_r8
      else if( w > 305._r8 .and. w <= 328._r8 ) then
         fo3qy = 0.0765_r8 &
                + a(1)*             (q1/(q1+q2))*exp( -((x(1)-w)/om(1))**4 ) &
                + a(2)*(t/300._r8)**2 *(q2/(q1+q2))*exp( -((x(2)-w)/om(2))**2 ) &
                + a(3)*(t/300._r8)**1.5_r8            *exp( -((x(3)-w)/om(3))**2 )
      else if( w > 328._r8 .and. w <= 340._r8 ) then
         fo3qy = 0.08_r8
      else if( w > 340._r8 ) then
         fo3qy = 0._r8
      end if

      end function fo3qy

      subroutine r04( nw, wl, wc, tlev, airlev, jlabel, xs )
!-----------------------------------------------------------------------------
!   purpose:
!   provide product of (cross section) x (quantum yiels) for n2o5 photolysis
!   reactions:
!        (a) n2o5 + hv -> no3 + no + o(3p)
!        (b) n2o5 + hv -> no3 + no2
!   cross section from jpl97: use tabulated values up to 280 nm, use expon.
!                             expression for >285nm, linearly interpolate
!                             between s(280) and s(285,t) in between
!   quantum yield: analysis of data in jpl94 (->dataj1/yld/n2o5.qy)
!-----------------------------------------------------------------------------
!   parameters:
!   nw     - integer, number of specified intervals + 1 in working        (i)
!            wavelength grid
!   wl     - real(r8), vector of lower limits of wavelength intervals in      (i)
!            working wavelength grid
!   wc     - real(r8), vector of center points of wavelength intervals in     (i)
!            working wavelength grid
!   nz     - integer, number of altitude levels in working altitude grid  (i)
!   tlev   - real(r8), temperature (k) at each specified altitude level       (i)
!   airlev - real(r8), air density (molec/cc) at each altitude level          (i)
!   j      - integer, counter for number of weighting functions defined  (io)
!   sq     - real(r8), cross section x quantum yield (cm^2) for each          (o)
!            photolysis reaction defined, at each defined wavelength and
!            at each defined altitude level
!   jlabel - character*40, string identifier for each photolysis reaction (o)
!            defined
!-----------------------------------------------------------------------------
!   edit history:
!   05/98  original, adapted from former jspec1 subroutine
!-----------------------------------------------------------------------------

      use mo_params
      use mo_waveall
      use ppgrid,     only : pverp

      implicit none

!-----------------------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------------------
      integer, intent(in)          :: nw
      real(r8), intent(in)         :: wl(kw)
      real(r8), intent(in)         :: wc(kw)
      real(r8), intent(in)         :: tlev(pverp)
      real(r8), intent(in)         :: airlev(pverp)
      real(r8), intent(inout)      :: xs(:,:)
      character(len=*), intent(in) :: jlabel

!-----------------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------------
      real(r8), parameter :: xs280 = 1.16999993e-19_r8

      integer  :: k, wn
      real(r8) :: qy
      real(r8) :: xsect, xst285
      real(r8) :: t

!-----------------------------------------------------------------------------
!	... n2o5 photodissociation
!-----------------------------------------------------------------------------
! 	... cross section from jpl97, table up to 280 nm
!           quantum yield : see dataj1/yld/n2o5.qy for explanation
!           correct for t-dependence of cross section
!-----------------------------------------------------------------------------
level_loop : &
      do k = 1,pverp
!-----------------------------------------------------------------------------
! 	... temperature dependence only valid for 225 - 300 k.
!-----------------------------------------------------------------------------
         t = 1._r8/max( 225._r8,min( tlev(k),300._r8 ) )
wave_loop : &
         do wn = 1,nw
            qy = max( 0._r8,min( 1._r8, 3.832441_r8 - 0.012809638_r8 * wc(wn) ) )
!-----------------------------------------------------------------------------
! 	... evaluate exponential
!-----------------------------------------------------------------------------
            if( wl(wn) >= 285._r8 .and. wl(wn+1) <= 380._r8 ) then
!!$               xs(wn,k) = qy        * 1.e-20_r8*exp( 2.735_r8 + (4728.5_r8 - 17.127_r8*wc(wn)) * t )
! fvitt - made to match the old jn2o5
               xs(wn,k) = (1._r8 - qy) * 1.e-20_r8*exp( 2.735_r8 + (4728.5_r8 - 17.127_r8*wc(wn)) * t )
!-----------------------------------------------------------------------------
! 	... between 280 and 285, interpolate between temperature evaluated exponential
!           at 285 nm and the tabulated value at 280 nm.
!-----------------------------------------------------------------------------
            else if( wl(wn) >= 280._r8 .and. wl(wn+1) <= 286._r8 ) then
               xst285 = 1.e-20_r8* exp( 2.735_r8 + (4728.5_r8 - 17.127_r8*286._r8)*t )
               xsect  = xs280 + (wc(wn) - 280._r8)*(xst285 - xs280)/(286._r8 - 280._r8)
!!$               xs(wn,k) = qy * xsect 
               xs(wn,k) = (1._r8-qy) * xsect 
!-----------------------------------------------------------------------------
! 	... use tabulated values
!-----------------------------------------------------------------------------
            else if (wl(wn) <= 280._r8 ) then
!!$               xs(wn,k) = qy * r04g(wn)
               xs(wn,k) = (1._r8-qy) * r04g(wn)
!-----------------------------------------------------------------------------
! 	... beyond 380 nm, set to zero
!-----------------------------------------------------------------------------
            else
               xs(wn,k) = 0._r8
            end if
         end do wave_loop
      end do level_loop

      end subroutine r04

      subroutine r44_inti( nw, wc )
!-----------------------------------------------------------------------------
!	... initialize subroutine r44
!-----------------------------------------------------------------------------

      use mo_params,   only : kw
      use abortutils,  only : endrun
      use cam_logfile, only : iulog

      implicit none

!-----------------------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------------------
      integer, intent(in)  :: nw
      real(r8), intent(in) :: wc(kw)

!-----------------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!	... cross sections according to jpl97 recommendation (identical to 94 rec.)
!           see file dataj1/abs/n2o_jpl94.abs for detail
!-----------------------------------------------------------------------------
      real(r8), parameter :: a0 = 68.21023_r8
      real(r8), parameter :: a1 = -4.071805_r8
      real(r8), parameter :: a2 = 4.301146e-02_r8
      real(r8), parameter :: a3 = -1.777846e-04_r8
      real(r8), parameter :: a4 = 2.520672e-07_r8

      real(r8), parameter :: b0 = 123.4014_r8
      real(r8), parameter :: b1 = -2.116255_r8
      real(r8), parameter :: b2 = 1.111572e-02_r8
      real(r8), parameter :: b3 = -1.881058e-05_r8

      integer  :: wn
      integer  :: astat
      real(r8) :: lambda

      allocate( a(nw), b(nw), stat=astat )
      if( astat /= 0 ) then
         write(iulog,*) 'r44_inti: a,b allocate failed; error = ',astat
         call endrun
      end if
      do wn = 1,nw
         lambda = wc(wn)   
         if( lambda >= 173._r8 .and. lambda <= 240._r8 ) then
            a(wn) = (((a4*lambda + a3)*lambda + a2)*lambda + a1)*lambda + a0
            b(wn) = ((b3*lambda + b2)*lambda + b1)*lambda + b0
         end if
      end do

      end subroutine r44_inti

      subroutine r44( nw, wl, wc, tlev, airlev, jlabel, xs )
!-----------------------------------------------------------------------------
!   purpose:
!   provide product (cross section) x (quantum yield) for n2o photolysis:
!               n2o + hv -> n2 + o(1d)
!   cross section: from jpl 97 recommendation
!   quantum yield: assumed to be unity, based on greenblatt and ravishankara
!-----------------------------------------------------------------------------
!   parameters:
!   nw     - integer, number of specified intervals + 1 in working        (i)
!            wavelength grid
!   wl     - real(r8), vector of lower limits of wavelength intervals in      (i)
!            working wavelength grid
!   wc     - real(r8), vector of center points of wavelength intervals in     (i)
!            working wavelength grid
!   nz     - integer, number of altitude levels in working altitude grid  (i)
!   tlev   - real(r8), temperature (k) at each specified altitude level       (i)
!   airlev - real(r8), air density (molec/cc) at each altitude level          (i)
!   j      - integer, counter for number of weighting functions defined  (io)
!   sq     - real(r8), cross section x quantum yield (cm^2) for each          (o)
!            photolysis reaction defined, at each defined wavelength and
!            at each defined altitude level
!   jlabel - character*40, string identifier for each photolysis reaction (o)
!            defined
!-----------------------------------------------------------------------------
!   edit history:
!   05/98  original, adapted from former jspec1 subroutine
!-----------------------------------------------------------------------------

      use mo_params
      use ppgrid,    only : pverp

      implicit none

!-----------------------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------------------
      integer, intent(in)          :: nw
      real(r8), intent(in)         :: wl(kw), wc(kw)
      real(r8), intent(in)         :: tlev(pverp)
      real(r8), intent(in)         :: airlev(pverp)
      real(r8), intent(inout)      :: xs(:,:)
      character(len=*), intent(in) :: jlabel

!-----------------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------------
      integer  :: k
      real(r8) :: t

!-----------------------------------------------------------------------------
!	... n2o photodissociation
!-----------------------------------------------------------------------------
!	... quantum yield of n(4s) and no(2pi) is less than 1% (greenblatt and
!           ravishankara), so quantum yield of o(1d) is assumed to be unity
!-----------------------------------------------------------------------------
      do k = 1,pverp
         t = max( 194._r8,min( tlev(k),320._r8 ) ) - 300._r8
         where( wc(:nw) >= 173._r8 .and. wc(:nw) <= 240._r8 )
            xs(:nw,k) = exp( a(:nw) + t*exp( b(:nw) ) )
         elsewhere
            xs(:nw,k) = 0._r8
         endwhere
      end do

      end subroutine r44

      subroutine r08_inti( nw, wl, wc )
!-----------------------------------------------------------------------------
!	... initialize subroutine r08
!-----------------------------------------------------------------------------

      use mo_params,   only : kw
      use abortutils,  only : endrun
      use cam_logfile, only : iulog

      implicit none

!-----------------------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------------------
      integer, intent(in)   :: nw
      real(r8), intent(in)  :: wl(kw)
      real(r8), intent(in)  :: wc(kw)

!-----------------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------------
      real(r8), parameter :: a0 = 6.4761e+04_r8
      real(r8), parameter :: a1 = -9.2170972e+02_r8
      real(r8), parameter :: a2 = 4.535649_r8
      real(r8), parameter :: a3 = -4.4589016e-03_r8
      real(r8), parameter :: a4 = -4.035101e-05_r8
      real(r8), parameter :: a5 = 1.6878206e-07_r8
      real(r8), parameter :: a6 = -2.652014e-10_r8
      real(r8), parameter :: a7 = 1.5534675e-13_r8

      real(r8), parameter :: b0 = 6.8123e+03_r8
      real(r8), parameter :: b1 = -5.1351e+01_r8
      real(r8), parameter :: b2 = 1.1522e-01_r8
      real(r8), parameter :: b3 = -3.0493e-05_r8
      real(r8), parameter :: b4 = -1.0924e-07_r8

      integer  :: astat
      integer  :: wn
      real(r8) :: lambda

      allocate( suma(nw), sumb(nw), stat=astat )
      if( astat /= 0 ) then
         write(iulog,*) 'r08_inti: suma,sumb allocate failed; error = ',astat
         call endrun
      end if
      do wn = 1,nw
         if( wl(wn) >= 260._r8 .and. wl(wn) < 350._r8 ) then
            lambda = wc(wn)
            suma(wn) = ((((((a7*lambda + a6)*lambda + a5)*lambda + a4)*lambda +a3)*lambda + a2)*lambda + a1)*lambda + a0
            sumb(wn) = (((b4*lambda + b3)*lambda + b2)*lambda + b1)*lambda + b0
         end if
      end do

      end subroutine r08_inti

      subroutine r08( nw, wl, wc, tlev, airlev, jlabel, xs )
!-----------------------------------------------------------------------------
!   purpose:
!   provide product of (cross section) x (quantum yield) for h2o2 photolysis
!          h2o2 + hv -> 2 oh
!   cross section:  from jpl97, tabulated values @ 298k for <260nm, t-depend.
!                   parameterization for 260-350nm
!   quantum yield:  assumed to be unity
!-----------------------------------------------------------------------------
!   parameters:
!   nw     - integer, number of specified intervals + 1 in working        (i)
!            wavelength grid
!   wl     - real(r8), vector of lower limits of wavelength intervals in      (i)
!            working wavelength grid
!   wc     - real(r8), vector of center points of wavelength intervals in     (i)
!            working wavelength grid
!   nz     - integer, number of altitude levels in working altitude grid  (i)
!   tlev   - real(r8), temperature (k) at each specified altitude level       (i)
!   airlev - real(r8), air density (molec/cc) at each altitude level          (i)
!   j      - integer, counter for number of weighting functions defined  (io)
!   sq     - real(r8), cross section x quantum yield (cm^2) for each          (o)
!            photolysis reaction defined, at each defined wavelength and
!            at each defined altitude level
!   jlabel - character*40, string identifier for each photolysis reaction (o)
!            defined
!-----------------------------------------------------------------------------
!   edit history:
!   05/98  original, adapted from former jspec1 subroutine
!-----------------------------------------------------------------------------

      use mo_params
      use mo_waveall
      use ppgrid,       only : pverp

      implicit none

!-----------------------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------------------
      integer, intent(in)          :: nw
      real(r8), intent(in)         :: wl(kw)
      real(r8), intent(in)         :: wc(kw)
      real(r8), intent(in)         :: tlev(pverp)
      real(r8), intent(in)         :: airlev(pverp)
      real(r8), intent(inout)      :: xs(:,:)
      character(len=*), intent(in) :: jlabel

!-----------------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------------
      integer  :: k
      real(r8) :: t
      real(r8) :: chi

!-----------------------------------------------------------------------------
!	... h2o2 photodissociation
!           cross section from lin et al. 1978
!-----------------------------------------------------------------------------
      do k = 1,pverp
         t = 1._r8/min( max( tlev(k),200._r8 ),400._r8 )            
         chi = 1._r8/(1._r8 + exp( -1265._r8*t ))
!-----------------------------------------------------------------------------
! 	... parameterization (jpl94)
!           range 260-350 nm; 200-400 k
!-----------------------------------------------------------------------------
         where( wl(:nw) > 260._r8 .and. wl(:nw) < 350._r8 )
            xs(:nw,k) = (chi * suma(:nw) + (1._r8 - chi)*sumb(:nw))*1.e-21_r8
         elsewhere
            xs(:nw,k) = r08g(:nw)
         endwhere
      end do

      end subroutine r08

      subroutine r06( nw, wl, wc, tlev, airlev, jlabel, xs )
!-----------------------------------------------------------------------------
!   purpose:
!   provide product of (cross section) x (quantum yield) for hno3 photolysis 
!         hno3 + hv -> oh + no2
!   cross section: burkholder et al., 1993
!   quantum yield: assumed to be unity
!-----------------------------------------------------------------------------
!   parameters:
!   nw     - integer, number of specified intervals + 1 in working        (i)
!            wavelength grid
!   wl     - real(r8), vector of lower limits of wavelength intervals in      (i)
!            working wavelength grid
!   wc     - real(r8), vector of center points of wavelength intervals in     (i)
!            working wavelength grid
!   nz     - integer, number of altitude levels in working altitude grid  (i)
!   tlev   - real(r8), temperature (k) at each specified altitude level       (i)
!   airlev - real(r8), air density (molec/cc) at each altitude level          (i)
!   j      - integer, counter for number of weighting functions defined  (io)
!   sq     - real(r8), cross section x quantum yield (cm^2) for each          (o)
!            photolysis reaction defined, at each defined wavelength and
!            at each defined altitude level
!   jlabel - character*40, string identifier for each photolysis reaction (o)
!            defined
!-----------------------------------------------------------------------------

      use mo_params
      use mo_waveall
      use ppgrid,       only : pverp
      
      implicit none

!-----------------------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------------------
      integer, intent(in)          :: nw
      real(r8), intent(in)         :: wl(kw)
      real(r8), intent(in)         :: wc(kw)
      real(r8), intent(in)         :: tlev(pverp)
      real(r8), intent(in)         :: airlev(pverp)
      real(r8), intent(inout)      :: xs(:,:)
      character(len=*), intent(in) :: jlabel

!-----------------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------------
      integer  :: k
      real(r8) :: wrk

!-----------------------------------------------------------------------------
!	... hno3 photodissociation
!-----------------------------------------------------------------------------
! 	... hno3 cross section parameters from burkholder et al. 1993
!           quantum yield = 1
!           correct for temperature dependence
!-----------------------------------------------------------------------------
      do k = 1,pverp
         wrk = 1.e-3_r8*(tlev(k) - 298._r8)
         xs(:nw,k) = r06g1(:nw) * 1.e-20_r8 * exp( r06g2(:nw)*wrk )
      end do

      end subroutine r06

      subroutine r10( nw, wl, wc, tlev, airlev, jlabel, xs )
!-----------------------------------------------------------------------------
!   purpose:
!   provide product of (cross section) x (quantum yield) for ch2o photolysis *
!         (a) ch2o + hv -> h + hco
!         (b) ch2o + hv -> h2 + co
!   cross section: choice between
!                  1) bass et al., 1980 (resolution: 0.025 nm)
!                  2) moortgat and schneider (resolution: 1 nm)
!                  3) cantrell et al. (orig res.) for > 301 nm,
!                     iupac 92, 97 elsewhere
!                  4) cantrell et al. (2.5 nm res.) for > 301 nm,
!                     iupac 92, 97 elsewhere
!                  5) rogers et al., 1990
!                  6) new ncar recommendation, based on averages of
!                     cantrell et al., moortgat and schneider, and rogers
!                     et al.
!   quantum yield: choice between
!                  1) evaluation by madronich 1991 (unpublished)
!                  2) iupac 89, 92, 97
!                  3) madronich, based on 1), updated 1998.
!-----------------------------------------------------------------------------
!   parameters:
!   nw     - integer, number of specified intervals + 1 in working        (i) 
!            wavelength grid
!   wl     - real(r8), vector of lower limits of wavelength intervals in      (i) 
!            working wavelength grid
!   wc     - real(r8), vector of center points of wavelength intervals in     (i) 
!            working wavelength grid
!   nz     - integer, number of altitude levels in working altitude grid  (i) 
!   tlev   - real(r8), temperature (k) at each specified altitude level       (i) 
!   airlev - real(r8), air density (molec/cc) at each altitude level          (i) 
!   j      - integer, counter for number of weighting functions defined  (io) 
!   sq     - real(r8), cross section x quantum yield (cm^2) for each          (o) 
!            photolysis reaction defined, at each defined wavelength and
!            at each defined altitude level
!   jlabel - character*40, string identifier for each photolysis reaction (o)
!            defined
!-----------------------------------------------------------------------------

      use mo_params
      use mo_waveall
      use ppgrid,       only : pverp

      implicit none

!-----------------------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------------------
      integer, intent(in)          :: nw
      real(r8), intent(in)         :: wl(kw)
      real(r8), intent(in)         :: wc(kw)
      real(r8), intent(in)         :: tlev(pverp)
      real(r8), intent(in)         :: airlev(pverp)
      real(r8), intent(inout)      :: xs(:,:)
      character(len=*), intent(in) :: jlabel


!-----------------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------------
      integer, parameter :: mopt1 = 6
      integer, parameter :: mopt2 = 1
      real(r8), parameter :: c0 = 1._r8/70._r8

      integer   :: k, wn
      real(r8)  :: phi1, phi2, phi20, ak300, akt
      real(r8)  :: qy, qy1, qy2, qy3, t, t1
      real(r8)  :: sigma, sig, slope

!-----------------------------------------------------------------------------
!	... ch2o photodissociatation
!-----------------------------------------------------------------------------
! 	... combine
!           y1 = xsect
!           y2 = xsect(223), cantrell et al.
!           y3 = xsect(293), cantrell et al.
!           y4 = qy for radical channel
!           y5 = qy for molecular channel
!           pressure and temperature dependent for w > 330.
!-----------------------------------------------------------------------------
level_loop : &
      do k = 1,pverp
         t = max( 223.15_r8,min( tlev(k),293.15_r8 ) )
         t1 = airlev(k)
wave_loop : &
         do wn = 1,nw
            if( mopt1 == 6 ) then
               sig = r10g2(wn)
            else
               sig = r10g1(wn)
            end if
!-----------------------------------------------------------------------------
! 	... correct cross section for temperature dependence for > 301. nm
!-----------------------------------------------------------------------------
            if( wl(wn) >= 301._r8 ) then 
               if( mopt1 == 3 .or. mopt1 == 6 ) then
                  sig = r10g2(wn) + r10g3(wn) * (t - 273.15_r8)
               else if( mopt1 == 4 ) then
                  slope = (r10g3(wn) - r10g2(wn)) * c0
                  slope = (r10g3(wn) - r10g2(wn)) * c0
                  sig = r10g2(wn) + slope * (t - 223._r8)
               end if
            end if
            sig = max( sig,0._r8 )
!-----------------------------------------------------------------------------
! 	... quantum yields:
!           temperature and pressure dependence beyond 330 nm
!-----------------------------------------------------------------------------
            qy1 = r10g4(wn)
            if( trim(jlabel) == 'jch2o_a' ) then
               xs(wn,k) = sig * qy1
            else
               if( wc(wn) >= 330._r8 .and. r10g5(wn) > 0._r8 ) then
                  phi1 = r10g4(wn)
                  phi2 = r10g5(wn)
                  phi20 = 1._r8 - phi1
                  ak300 = (phi20 - phi2)/(phi20*phi2*2.54e+19_r8)
                  akt = ak300*(1._r8 + 61.69_r8*(1._r8 - tlev(k)/300._r8)*(wc(wn)/329._r8 - 1._r8))
                  qy2 = 1._r8 / ((1._r8/phi20) + t1*akt)
               else
                  qy2 = r10g5(wn)
               end if
               qy2 = min( 1._r8,max( 0._r8,qy2 ) )
               xs(wn,k) = sig * qy2
            end if
         end do wave_loop
      end do level_loop

      end subroutine r10

      subroutine r11( nw, wl, wc, tlev, airlev, jlabel, xs )
!-----------------------------------------------------------------------------
!   purpose:
!   provide product (cross section) x (quantum yield) for ch3cho photolysis:
!       (a)  ch3cho + hv -> ch3 + hco
!       (b)  ch3cho + hv -> ch4 + co
!       (c)  ch3cho + hv -> ch3co + h
!   cross section:  choice between
!                    (1) iupac 97 data, from martinez et al.
!                    (2) calvert and pitts
!                    (3) martinez et al., table 1 scanned from paper
!                    (4) kfa tabulations
!   quantum yields: choice between
!                    (1) iupac 97, pressure correction using horowith and
!                                  calvert, 1982
!                    (2) ncar data file, from moortgat, 1986
!-----------------------------------------------------------------------------
!   parameters:
!   nw     - integer, number of specified intervals + 1 in working        (i) 
!            wavelength grid
!   wl     - real(r8), vector of lower limits of wavelength intervals in      (i) 
!            working wavelength grid
!   wc     - real(r8), vector of center points of wavelength intervals in     (i) 
!            working wavelength grid
!   nz     - integer, number of altitude levels in working altitude grid  (i) 
!   tlev   - real(r8), temperature (k) at each specified altitude level       (i) 
!   airlev - real(r8), air density (molec/cc) at each altitude level          (i) 
!   j      - integer, counter for number of weighting functions defined  (io) 
!   sq     - real(r8), cross section x quantum yield (cm^2) for each          (o) 
!            photolysis reaction defined, at each defined wavelength and      
!            at each defined altitude level
!   jlabel - character*40, string identifier for each photolysis reaction (o) 
!            defined
!-----------------------------------------------------------------------------

      use mo_params
      use mo_waveall
      use ppgrid,       only : pverp

      implicit none

!-----------------------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------------------
      integer, intent(in)          :: nw
      real(r8), intent(in)         :: wl(kw)
      real(r8), intent(in)         :: wc(kw)
      real(r8), intent(in)         :: tlev(pverp)
      real(r8), intent(in)         :: airlev(pverp)
      real(r8), intent(inout)      :: xs(:,:)
      character(len=*), intent(in) :: jlabel

!-----------------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------------
      real(r8), parameter :: c0 = 1._r8/2.465e19_r8
      integer, parameter  :: mabs = 3
      integer, parameter  :: myld = 1

      integer  :: k, wn
      real(r8) :: qy1
      real(r8) :: sig, t

!-----------------------------------------------------------------------------
!	... ch3cho photolysis
!           1:  ch3 + hco
!           2:  ch4 + co
!           3:  ch3co + h
!-----------------------------------------------------------------------------
! 	... options
!           mabs for cross sections
!           myld for quantum yields
!
!           absorption:
!           1:  iupac-97 data, from martinez et al.
!           2:  calvert and pitts
!           3:  martinez et al., table 1 scanned from paper
!           4:  kfa tabulations, 6 choices, see file open statements
!
!           quantum yield
!           1:  dataj1/ch3cho/ch3cho_iup.yld
!               pressure correction using horowitz and calvert 1982, 
!               based on slope/intercept of stern-volmer plots
!
!           2:  ncar data file, from moortgat 1986.
!-----------------------------------------------------------------------------

      do k = 1,pverp
         t = airlev(k)*c0
         do wn = 1,nw
            sig = r11g(wn)
!-----------------------------------------------------------------------------
! 	... pressure correction for channel 1, ch3 + cho
!           based on horowitz and calvert 1982.
!-----------------------------------------------------------------------------
            if( trim( jlabel ) == 'jch3cho_a' ) then
               qy1 = r11g1(wn) * (1._r8 + r11g4(wn))/(1._r8 + r11g4(wn)*t)
               qy1 = min( 1._r8,max( 0._r8,qy1 ) )
               xs(wn,k) = sig * qy1
            else if( trim( jlabel ) == 'jch3cho_b' ) then
               xs(wn,k) = sig * r11g2(wn)
            else if( trim( jlabel ) == 'jch3cho_c' ) then
               xs(wn,k) = sig * r11g3(wn)
            end if
         end do
      end do

      end subroutine r11

      subroutine r14( nw, wl, wc, tlev, airlev, jlabel, xs )
!-----------------------------------------------------------------------------
!   purpose:
!   provide the product (cross section) x (quantum yield) for ch3cocho
!   photolysis:
!            ch3cocho + hv -> products
!
!   cross section: choice between
!                   (1) from meller et al., 1991, as tabulated by iupac 97
!                          5 nm resolution (table 1) for < 402 nm
!                          2 nm resolution (table 2) for > 402 nm
!                   (2) average at 1 nm of staffelbach et al., 1995, and
!                       meller et al., 1991
!                   (3) plum et al., 1983, as tabulated by kfa
!                   (4) meller et al., 1991 (0.033 nm res.), as tab. by kfa
!                   (5) meller et al., 1991 (1.0 nm res.), as tab. by kfa
!                   (6) staffelbach et al., 1995, as tabulated by kfa
!   quantum yield: choice between
!                   (1) plum et al., fixed at 0.107
!                   (2) plum et al., divided by 2, fixed at 0.0535
!                   (3) staffelbach et al., 0.45 for < 300 nm, 0 for > 430 nm
!                       linear interp. in between
!                   (4) koch and moortgat, prv. comm., 1997
!-----------------------------------------------------------------------------
!   parameters:
!   nw     - integer, number of specified intervals + 1 in working        (i) 
!            wavelength grid
!   wl     - real(r8), vector of lower limits of wavelength intervals in      (i) 
!            working wavelength grid
!   wc     - real(r8), vector of center points of wavelength intervals in     (i) 
!            working wavelength grid
!   nz     - integer, number of altitude levels in working altitude grid  (i) 
!   tlev   - real(r8), temperature (k) at each specified altitude level       (i) 
!   airlev - real(r8), air density (molec/cc) at each altitude level          (i) 
!   j      - integer, counter for number of weighting functions defined  (io) 
!   sq     - real(r8), cross section x quantum yield (cm^2) for each          (o) 
!            photolysis reaction defined, at each defined wavelength and      
!            at each defined altitude level
!   jlabel - character*40, string identifier for each photolysis reaction (o) 
!            defined
!-----------------------------------------------------------------------------

      use mo_params
      use mo_waveall
      use ppgrid,       only : pverp

      implicit none

!-----------------------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------------------
      integer, intent(in)          :: nw
      real(r8), intent(in)         :: wl(kw)
      real(r8), intent(in)         :: wc(kw)
      real(r8), intent(in)         :: tlev(pverp)
      real(r8), intent(in)         :: airlev(pverp)
      real(r8), intent(inout)      :: xs(:,:)
      character(len=*), intent(in) :: jlabel

!-----------------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------------
      real(r8), parameter :: sig2 = .45_r8
      real(r8), parameter :: sig3 = 1._r8/130._r8
      real(r8), parameter :: expfac = 760._r8/2.456e19_r8

      integer  :: k, wn
      real(r8) :: qy
      real(r8) :: phi0, kq

!-----------------------------------------------------------------------------
!	... ch3cocho photolysis
!-----------------------------------------------------------------------------
! 	... options
!           mabs for cross sections
!           myld for quantum yields
!
!           absorption:
!           1:  from meller et al. (1991), as tabulated by iupac-97
!               for wc < 402, use coarse data (5 nm, table 1)
!               for wc > 402, use finer data (2 nm, table 2)
!           2: average at 1nm of  staffelbach et al. 1995 and meller et al. 1991
!               cross section from kfa tables:
!           3: ch3cocho.001 - plum et al. 1983
!           4: ch3cocho.002 - meller et al. 1991, 0.033 nm resolution
!           5: ch3cocho.003 - meller et al. 1991, 1.0   nm resolution
!           6: ch3cocho.004 - staffelbach et al. 1995
!
!           quantum yield
!           1:  plum et al., 0.107
!           2:  plum et al., divided by two = 0.0535
!           3:  staffelbach et al., 0.45 at wc .le. 300, 0 for wc .gt. 430, linear 
!               interp in between
!           4:  koch and moortgat, prv. comm. 1997. - pressure-dependent
!         * 5:  Chen, Y., W. Wang, and L. Zhu, Wavelength-dependent photolysis of methylglyoxal
!         *      in the 290-440 nm region, J Phys Chem A, 104, 11126-11131, 2000
!-----------------------------------------------------------------------------
      do k = 1,pverp
         do wn = 1,nw
            phi0 = 1._r8 - (wc(wn) - 380._r8)/60._r8
            phi0 = max( 0._r8,min( phi0,1._r8 ) )
            kq = 1.36e8_r8 * exp( -8793._r8/wc(wn) )
            if( phi0 > 0._r8 ) then
               if( wc(wn) >= 380._r8 .and. wc(wn) <= 440._r8 ) then
                  qy = phi0 / (phi0 + kq * airlev(k) * expfac )
               else
                  qy = phi0
               end if
            else
               qy = 0._r8
            end if
            xs(wn,k) = r14g(wn) * qy
         end do
      end do

      end subroutine r14

      subroutine r15( nw, wl, wc, tlev, airlev, jlabel, xs )
!-----------------------------------------------------------------------------
!   purpose:
!   provide product (cross section) x (quantum yield) for ch3coch3 photolysis
!           ch3coch3 + hv -> products
!
!   cross section:  choice between
!                    (1) calvert and pitts
!                    (2) martinez et al., 1991, alson in iupac 97
!                    (3) noaa, 1998, unpublished as of 01/98
!   quantum yield:  choice between
!                    (1) gardiner et al, 1984
!                    (2) iupac 97
!                    (3) mckeen et al., 1997
!-----------------------------------------------------------------------------
!   parameters:
!   nw     - integer, number of specified intervals + 1 in working        (i) 
!            wavelength grid
!   wl     - real(r8), vector of lower limits of wavelength intervals in      (i) 
!            working wavelength grid
!   wc     - real(r8), vector of center points of wavelength intervals in     (i) 
!            working wavelength grid
!   nz     - integer, number of altitude levels in working altitude grid  (i) 
!   tlev   - real(r8), temperature (k) at each specified altitude level       (i) 
!   airlev - real(r8), air density (molec/cc) at each altitude level          (i) 
!   j      - integer, counter for number of weighting functions defined  (io) 
!   sq     - real(r8), cross section x quantum yield (cm^2) for each          (o) 
!            photolysis reaction defined, at each defined wavelength and      
!            at each defined altitude level
!   jlabel - character*40, string identifier for each photolysis reaction (o) 
!            defined
!-----------------------------------------------------------------------------

      use mo_params
      use mo_waveall
      use ppgrid,       only : pverp

      implicit none

!-----------------------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------------------
      integer, intent(in)          :: nw
      real(r8), intent(in)         :: wl(kw)
      real(r8), intent(in)         :: wc(kw)
      real(r8), intent(in)         :: tlev(pverp)
      real(r8), intent(in)         :: airlev(pverp)
      real(r8), intent(inout)      :: xs(:,:)
      character(len=*), intent(in) :: jlabel

!-----------------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------------
      integer, parameter :: mabs = 2
      integer, parameter :: myld = 3

      integer  :: k, wn
      real(r8) :: qy
      real(r8) :: sig
      real(r8) :: a, b, t, t1
      real(r8) :: m, fco, fac, w

!-----------------------------------------------------------------------------
!	... ch3coch3 photodissociation
!-----------------------------------------------------------------------------
! 	... options
!           mabs for cross sections
!           myld for quantum yields
!
!           absorption:
!           1:  cross section from calvert and  pitts
!           2:  martinez et al. 1991, also in iupac97
!           3:  noaa 1998, unpublished as of jan 98.
!
!           quantum yield
!           1:  gardiner et al. 1984
!           2:  iupac 97
!           3:  mckeen, s. a., t. gierczak, j. b. burkholder, p. o. wennberg, t. f. hanisco,
!               e. r. keim, r.-s. gao, s. c. liu, a. r. ravishankara, and d. w. fahey, 
!               the photochemistry of acetone in the upper troposphere:  a source of 
!               odd-hydrogen radicals, geophys. res. lett., 24, 3177-3180, 1997.
!           4:  Blitz, M. A., D. E. Heard, M. J. Pilling, S. R. Arnold, and M. P. Chipperfield 
!              (2004), Pressure and temperature-dependent quantum yields for the 
!               photodissociation of acetone between 279 and 327.5 nm, Geophys. 
!               Res. Lett., 31, L06111, doi:10.1029/2003GL018793.
!
!-----------------------------------------------------------------------------
      do k = 1,pverp
         m = airlev(k)
         t = tlev(k)
         do wn = 1,nw
            sig = r15g(wn)
            w   = wc(wn)
            call qyacet( w, t, m, fco, fac )
            qy = min( 1._r8,max( 0._r8,fac ) )
            xs(wn,k) = sig*qy
         end do
      end do

      end subroutine r15

      subroutine qyacet( w, t, m, fco, fac )
!-----------------------------------------------------------------------------
! Compute acetone quantum yields according to the parameterization of:
! Blitz, M. A., D. E. Heard, M. J. Pilling, S. R. Arnold, and M. P. Chipperfield 
!       (2004), Pressure and temperature-dependent quantum yields for the 
!       photodissociation of acetone between 279 and 327.5 nm, Geophys. 
!       Res. Lett., 31, L06111, doi:10.1029/2003GL018793.
!-----------------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------------------
      real(r8), intent(in)  :: w            ! w = wavelength (nm)
      real(r8), intent(in)  :: t            ! T = temperature (K)
      real(r8), intent(in)  :: m            ! m = air number density (molec/cm^3)
      real(r8), intent(out) :: fco          ! fco = quantum yield for product CO
      real(r8), intent(out) :: fac          ! fac = quantum yield for product CH3CO (acetyl radical)

!-----------------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------------
      real(r8) :: a0, a1, a2, a3, a4
      real(r8) :: b0, b1, b2, b3, b4
      real(r8) :: c3
      real(r8) :: cA0, cA1, cA2, cA3, cA4
      real(r8) :: tratio
      real(r8) :: wrk

!-----------------------------------------------------------------------------
!** set out-of-range values:
! use low pressure limits for shorter wavelengths
! set to zero beyound 327.5
!-----------------------------------------------------------------------------
      if( w < 279._r8 ) then
         fco = 0.05_r8
         fac = 0.95_r8
      else if( w > 327.5_r8 ) then
         fco = 0._r8
         fac = 0._r8
      else
!-----------------------------------------------------------------------------
!	... CO (carbon monoxide) quantum yields
!-----------------------------------------------------------------------------
      tratio = t/295._r8
      a0  = .350_r8 * tratio**(-1.28_r8)
      b0  = .068_r8 * tratio**(-2.65_r8)
      cA0 = exp( b0*(w - 248._r8) ) * a0 / (1._r8 - a0)
      fco = 1._r8 / (1._r8 + cA0)
!-----------------------------------------------------------------------------
!	... CH3CO (acetyl radical) quantum yields:
!-----------------------------------------------------------------------------
      if( w >= 279._r8 .and. w < 302._r8 ) then
         a1  = 1.600e-19_r8 * tratio**(-2.38_r8)
         b1  = 0.55e-3_r8   * tratio**(-3.19_r8)
         cA1 = a1 * exp( -b1*((1.e7_r8/w) - 33113._r8) )
         fac = (1._r8 - fco) / (1._r8 + cA1 * m)
      else if( w >= 302._r8 .and. w < 327.5_r8 ) then
         a2  = 1.62e-17_r8 * tratio**(-10.03_r8)
         b2  = 1.79e-3_r8  * tratio**(-1.364_r8)
         wrk = 1.e7_r8/w
         cA2 = a2 * exp( -b2*(wrk - 30488._r8) )
         a3  = 26.29_r8   * tratio**(-6.59_r8)
         b3  = 5.72e-7_r8 * tratio**(-2.93_r8)
         c3  = 30006._r8  * tratio**(-0.064_r8)
         ca3 = a3 * exp( -b3*(wrk - c3)**2 )
         a4  = 1.67e-15_r8 * tratio**(-7.25_r8)
         b4  = 2.08e-3_r8  * tratio**(-1.16_r8)
         cA4 = a4 * exp( -b4*(wrk - 30488._r8) )
         fac = (1._r8 - fco) * (1._r8 + cA3 + cA4 * m) &
               /((1._r8 + cA3 + cA2 * M)*(1._r8 + cA4 * m))
         end if
      end if

      end subroutine qyacet

      subroutine r17( nw, wl, wc, tlev, airlev, jlabel, xs )
!-----------------------------------------------------------------------------
!   purpose:
!   provide product (cross section) x (quantum yield) for ch3ono2
!   photolysis:
!           ch3ono2 + hv -> ch3o + no2
!
!   cross section: choice between
!                   (1) calvert and pitts, 1966
!                   (2) talukdar, burkholder, hunter, gilles, roberts,
!                       ravishankara, 1997
!                   (3) iupac 97, table of values for 198k
!                   (4) iupac 97, temperature-dependent equation
!                   (5) taylor et al, 1980
!                   (6) fit from roberts and fajer, 1989
!                   (7) rattigan et al., 1992
!                   (8) libuda and zabel, 1995
!   quantum yield: assumed to be unity
!-----------------------------------------------------------------------------
!   parameters:
!   nw     - integer, number of specified intervals + 1 in working        (i) 
!            wavelength grid
!   wl     - real(r8), vector of lower limits of wavelength intervals in      (i) 
!            working wavelength grid
!   wc     - real(r8), vector of center points of wavelength intervals in     (i) 
!            working wavelength grid
!   nz     - integer, number of altitude levels in working altitude grid  (i) 
!   tlev   - real(r8), temperature (k) at each specified altitude level       (i) 
!   airlev - real(r8), air density (molec/cc) at each altitude level          (i) 
!   j      - integer, counter for number of weighting functions defined  (io) 
!   sq     - real(r8), cross section x quantum yield (cm^2) for each          (o) 
!            photolysis reaction defined, at each defined wavelength and      
!            at each defined altitude level
!   jlabel - character*40, string identifier for each photolysis reaction (o) 
!            defined
!-----------------------------------------------------------------------------

      use mo_params
      use mo_waveall
      use ppgrid,       only : pverp

      implicit none

!-----------------------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------------------
      integer, intent(in)          :: nw
      real(r8), intent(in)         :: wl(kw)
      real(r8), intent(in)         :: wc(kw)
      real(r8), intent(in)         :: tlev(pverp)
      real(r8), intent(in)         :: airlev(pverp)
      real(r8), intent(inout)      :: xs(:,:)
      character(len=*), intent(in) :: jlabel

!-----------------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------------
      integer  :: k
      real(r8) :: t

!-----------------------------------------------------------------------------
!	... ch3ono2 photodissociation
!-----------------------------------------------------------------------------
! 	... mabs: absorption cross section options:
!           1:  calvert and  pitts 1966
!           2:  talukdar, burkholder, hunter, gilles, roberts, ravishankara, 1997.
!           3:  iupac-97, table of values for 298k.
!           4:  iupac-97, temperature-dependent equation
!           5:  taylor et al. 1980
!           6:  fit from roberts and fajer, 1989
!           7:  rattigan et al. 1992
!           8:  libuda and zabel 1995
!-----------------------------------------------------------------------------
      do k = 1,pverp
         t = tlev(k) - 298._r8
         xs(:nw,k) = r17g(:nw) * exp( r17g1(:nw)*t )
      end do

      end subroutine r17

      subroutine r18( nw, wl, wc, tlev, airlev, jlabel, xs )
!-----------------------------------------------------------------------------
!   purpose:
!   provide product (cross section) x (quantum yield) for pan photolysis:
!        pan + hv -> products
!
!   cross section: from talukdar et al., 1995
!   quantum yield: assumed to be unity
!-----------------------------------------------------------------------------
!   parameters:
!   nw     - integer, number of specified intervals + 1 in working        (i) 
!            wavelength grid
!   wl     - real(r8), vector of lower limits of wavelength intervals in      (i) 
!            working wavelength grid
!   wc     - real(r8), vector of center points of wavelength intervals in     (i) 
!            working wavelength grid
!   nz     - integer, number of altitude levels in working altitude grid  (i) 
!   tlev   - real(r8), temperature (k) at each specified altitude level       (i) 
!   airlev - real(r8), air density (molec/cc) at each altitude level          (i) 
!   j      - integer, counter for number of weighting functions defined  (io) 
!   sq     - real(r8), cross section x quantum yield (cm^2) for each          (o) 
!            photolysis reaction defined, at each defined wavelength and      
!            at each defined altitude level
!   jlabel - character*40, string identifier for each photolysis reaction (o) 
!            defined
!-----------------------------------------------------------------------------

      use mo_params
      use mo_waveall
      use ppgrid,       only : pverp

      implicit none

!-----------------------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------------------
      integer, intent(in)          :: nw
      real(r8), intent(in)         :: wl(kw)
      real(r8), intent(in)         :: wc(kw)
      real(r8), intent(in)         :: tlev(pverp)
      real(r8), intent(in)         :: airlev(pverp)
      real(r8), intent(inout)      :: xs(:,:)
      character(len=*), intent(in) :: jlabel

!-----------------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------------
      integer  :: k
      real(r8) :: t

!-----------------------------------------------------------------------------
!	... pan photodissociation
!-----------------------------------------------------------------------------
! 	... cross section from senum et al., 1984, j.phys.chem. 88/7, 1269-1270
!           quantum yield
!           yet unknown, but assumed to be 1.0 (talukdar et al., 1995)
!-----------------------------------------------------------------------------
      do k = 1,pverp
         t = tlev(k) - 298._r8
         xs(:nw,k) = r18g(:nw) * exp( r18g2(:nw)*t )
      end do 

      end subroutine r18

      subroutine xs_mvk( nw, wl, wc, tlev, airlev, xs )
!-----------------------------------------------------------------------------
!   purpose:
!   provide product (cross section) x (quantum yield) for mvk photolysis:
!        mvk + hv -> products
!-----------------------------------------------------------------------------
!   parameters:
!   nw     - integer, number of specified intervals + 1 in working        (i) 
!            wavelength grid
!   wl     - real(r8), vector of lower limits of wavelength intervals in      (i) 
!            working wavelength grid
!   wc     - real(r8), vector of center points of wavelength intervals in     (i) 
!            working wavelength grid
!   nz     - integer, number of altitude levels in working altitude grid  (i) 
!   tlev   - real(r8), temperature (k) at each specified altitude level       (i) 
!   airlev - real(r8), air density (molec/cc) at each altitude level          (i) 
!   sq     - real(r8), cross section x quantum yield (cm^2) for each          (o) 
!            photolysis reaction defined, at each defined wavelength and      
!            at each defined altitude level
!-----------------------------------------------------------------------------

      use mo_params
      use mo_waveall
      use ppgrid,       only : pverp

      implicit none

!-----------------------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------------------
      integer, intent(in)      :: nw
      real(r8), intent(in)     :: wl(kw)
      real(r8), intent(in)     :: wc(kw)
      real(r8), intent(in)     :: tlev(pverp)
      real(r8), intent(in)     :: airlev(pverp)
      real(r8), intent(inout)  :: xs(:,:)

!-----------------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------------
      integer  :: k
      integer  :: w
      real(r8) :: denomi
      real(r8) :: qy(nw)

      do k = 1,pverp
         denomi = 1._r8/(5.5_r8 + 9.2e-19_r8*airlev(k))
         qy(:)  = exp( -.055_r8*(wc(:nw) - 308._r8) )*denomi
         xs(:nw,k) = min( qy(:nw),1._r8 ) * xs(:nw,k)
      end do

      end subroutine xs_mvk

      end module mo_xsections
