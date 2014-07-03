
      module mo_ps2str

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: ps2str

      contains

      subroutine ps2str( nw, zen, rsfc, tauu, omu, &
                         gu, dsdh, nid, radfld )
!-----------------------------------------------------------------------------
!   purpose:
!   solve two-stream equations for multiple layers.  the subroutine is based
!   on equations from:  toon et al., j.geophys.res., v94 (d13), nov 20, 1989.
!   it contains 9 two-stream methods to choose from.  a pseudo-spherical
!   correction has also been added.
!-----------------------------------------------------------------------------
!   parameters:
!   nlevel  - integer, number of specified altitude levels in the working (i)
!             grid
!   zen     - real(r8), solar zenith angle (degrees)                          (i)
!   rsfc    - real(r8), surface albedo at current wavelength                  (i)
!   tauu    - real(r8), unscaled optical depth of each layer                  (i)
!   omu     - real(r8), unscaled single scattering albedo of each layer       (i)
!   gu      - real(r8), unscaled asymmetry parameter of each layer            (i)
!   dsdh    - real(r8), slant path of direct beam through each layer crossed  (i)
!             when travelling from the top of the atmosphere to layer i;
!             dsdh(i,j), i = 0..nz-1, j = 1..nz-1
!   nid     - integer, number of layers crossed by the direct beam when   (i)
!             travelling from the top of the atmosphere to layer i;
!             nid(i), i = 0..nz-1
!   delta   - logical, switch to use delta-scaling                        (i)
!             .true. -> apply delta-scaling
!             .false.-> do not apply delta-scaling
!   fdr     - real(r8), contribution of the direct component to the total     (o)
!             actinic flux at each altitude level
!   fup     - real(r8), contribution of the diffuse upwelling component to    (o)
!             the total actinic flux at each altitude level
!   fdn     - real(r8), contribution of the diffuse downwelling component to  (o)
!             the total actinic flux at each altitude level
!   edr     - real(r8), contribution of the direct component to the total     (o)
!             spectral irradiance at each altitude level
!   eup     - real(r8), contribution of the diffuse upwelling component to    (o)
!             the total spectral irradiance at each altitude level
!   edn     - real(r8), contribution of the diffuse downwelling component to  (o)
!             the total spectral irradiance at each altitude level
!-----------------------------------------------------------------------------

      use mo_params,    only : smallest, largest
      use mo_constants, only : d2r
      use ppgrid,       only : pver, pverp
      use mo_trislv,    only : tridec, trislv

      implicit none

!-----------------------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------------------
      integer, intent(in)     :: nw
      integer, intent(in)     :: nid(0:pver)
      real(r8), intent(in)    :: zen
      real(r8), intent(in)    :: rsfc(nw)
      real(r8), intent(in)    :: tauu(pver,nw)
      real(r8), intent(in)    :: omu(pver,nw)
      real(r8), intent(in)    :: gu(pver,nw)
      real(r8), intent(in)    :: dsdh(0:pver,pver)
      real(r8), intent(out)   :: radfld(pverp,nw)

!-----------------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! 	... mu = cosine of solar zenith angle
!           rsfc = surface albedo
!           tauu =  unscaled optical depth of each layer
!           omu  =  unscaled single scattering albedo
!           gu   =  unscaled asymmetry factor
!           klev = max dimension of number of layers in atmosphere
!           nlayer = number of layers in the atmosphere
!           nlevel = nlayer + 1 = number of levels
!-----------------------------------------------------------------------------
      integer, parameter  :: mrows = 2*pver
      integer, parameter  :: pverm = pver - 1
      real(r8), parameter :: eps = 1.e-3_r8
      real(r8), parameter :: pifs = 1._r8
      real(r8), parameter :: fdn0 = 0._r8

      integer :: row
      integer :: lev
      integer :: i, ip1, wn
      integer :: j, jl, ju

      real(r8) :: precis, wrk
      real(r8) :: tempg
      real(r8) :: mu, suma
      real(r8) :: g, om
      real(r8) :: gam1, gam2, gam3, gam4
      real(r8), dimension(pver)    :: f, gi, omi
      real(r8), dimension(0:pverp) :: tauc, mu2
      real(r8), dimension(pver)    :: lam, taun, bgam
      real(r8), dimension(pver)    :: cdn
      real(r8), dimension(0:pverp,nw) :: tausla
      real(r8), dimension(pver,nw)  :: cup, cuptn, cdntn
      real(r8), dimension(pver,nw)  :: e1, e2, e3, e4
      real(r8), dimension(mrows)    :: a, b, d, e
      real(r8), dimension(nw,mrows) :: sub, main, super, y
!-----------------------------------------------------------------------------
! 	... for calculations of associated legendre polynomials for gama1,2,3,4
!           in delta-function, modified quadrature, hemispheric constant,
!           hybrid modified eddington-delta function metods, p633,table1.
!           w.e.meador and w.r.weaver, gas,1980,v37,p.630
!           w.j.wiscombe and g.w. grams, gas,1976,v33,p2440, 
!           uncomment the following two lines and the appropriate statements
!           further down.
!-----------------------------------------------------------------------------
      real(r8) :: expon, expon0, expon1, divisr, temp, up, dn
      real(r8) :: ssfc

!-----------------------------------------------------------------------------
! 	... initial conditions:  pi*solar flux = 1;  diffuse incidence = 0
!-----------------------------------------------------------------------------
      precis = epsilon( precis )

      mu = cos( zen*d2r )
wave_loop : &
      do wn = 1,nw
!-----------------------------------------------------------------------------
!	... compute coefficients for each layer:
!           gam1 - gam4 = 2-stream coefficients, different for different approximations
!           expon0 = calculation of e when tau is zero
!           expon1 = calculation of e when tau is taun
!           cup and cdn = calculation when tau is zero
!           cuptn and cdntn = calc. when tau is taun
!           divisr = prevents division by zero
!-----------------------------------------------------------------------------
         tauc(0:pverp)      = 0._r8
         tausla(0:pverp,wn) = 0._r8
         mu2(0:pverp)       = sqrt( smallest )

!-----------------------------------------------------------------------------
! 	... delta-scaling. have to be done for delta-eddington approximation, 
!           delta discrete ordinate, practical improved flux method, delta function,
!           and hybrid modified eddington-delta function methods approximations
!-----------------------------------------------------------------------------
         f(1:pver)    = gu(:,wn)*gu(:,wn)
         gi(1:pver)   = (gu(:,wn) - f(1:pver))/(1._r8 - f(1:pver))
         omi(1:pver)  = (1._r8 - f(1:pver))*omu(1:pver,wn)/(1._r8 - omu(1:pver,wn)*f(1:pver))       
         taun(1:pver) = (1._r8 - omu(1:pver,wn)*f(1:pver))*tauu(1:pver,wn)

!-----------------------------------------------------------------------------
! 	... calculate slant optical depth at the top of the atmosphere when zen>90.
!           in this case, higher altitude of the top layer is recommended which can 
!           be easily changed in gridz.f.
!-----------------------------------------------------------------------------
         if( zen > 90._r8 ) then
            if( nid(0) < 0 ) then
               tausla(0,wn) = largest
            else
	       ju = nid(0)
               tausla(0,wn) = 2._r8*dot_product( taun(1:ju),dsdh(0,1:ju) )
            end if
         end if
level_loop : &  
         do i = 1,pver
            g = gi(i)
            om = omi(i)
            tauc(i) = tauc(i-1) + taun(i)
!-----------------------------------------------------------------------------
! 	... stay away from 1 by precision.  for g, also stay away from -1
!-----------------------------------------------------------------------------
            tempg = min( abs(g),1._r8 - precis )
            g     = sign( tempg,g )
            om    = min( om,1._r8 - precis )
!-----------------------------------------------------------------------------
! 	... calculate slant optical depth
!-----------------------------------------------------------------------------
            if( nid(i) < 0 ) then
               tausla(i,wn) = largest
            else
               ju = min( nid(i),i )
               suma = dot_product( taun(1:ju),dsdh(i,1:ju) )
               jl = min( nid(i),i ) + 1
               tausla(i,wn) = suma + 2._r8*dot_product( taun(jl:nid(i)),dsdh(i,jl:nid(i)) )
               if( tausla(i,wn) == tausla(i-1,wn) ) then
                 mu2(i) = sqrt( largest )
               else
                 mu2(i) = (tauc(i) - tauc(i-1))/(tausla(i,wn) - tausla(i-1,wn))
                 mu2(i) = sign( max( abs(mu2(i)),sqrt(smallest) ),mu2(i) )
               end if
            end if
!-----------------------------------------------------------------------------
!	... the following gamma equations are from pg 16,289, table 1
!           eddington approximation(joseph et al., 1976, jas, 33, 2452):
!-----------------------------------------------------------------------------
            gam1 =  .25_r8*(7._r8 - om*(4._r8 + 3._r8*g))
            gam2 = -.25_r8*(1._r8 - om*(4._r8 - 3._r8*g))
            gam3 = .25_r8*(2._r8 - 3._r8*g*mu)
            gam4 = 1._r8 - gam3
!-----------------------------------------------------------------------------
! 	... lambda = pg 16,290 equation 21
!           big gamma = pg 16,290 equation 22
!-----------------------------------------------------------------------------
            lam(i) = sqrt( gam1*gam1 - gam2*gam2 )
            bgam(i) = (gam1 - lam(i))/gam2
	    wrk = lam(i)*taun(i)
	    if( wrk < 500._r8 ) then
               expon = exp( -wrk )
	    else
               expon = 0._r8
	    end if
!-----------------------------------------------------------------------------
! 	... e1 - e4 = pg 16,292 equation 44
!-----------------------------------------------------------------------------
            e1(i,wn) = 1._r8 + bgam(i)*expon
            e2(i,wn) = 1._r8 - bgam(i)*expon
            e3(i,wn) = bgam(i) + expon
            e4(i,wn) = bgam(i) - expon
!-----------------------------------------------------------------------------
! 	... the following sets up for the c equations 23, and 24
!           found on page 16,290
!           prevent division by zero (if lambda=1/mu, shift 1/mu^2 by eps = 1.e-3
!           which is approx equiv to shifting mu by 0.5*eps* (mu)**3
!-----------------------------------------------------------------------------
	    if( tausla(i-1,wn) < 500._r8 ) then
               expon0 = exp( -tausla(i-1,wn) )
	    else
               expon0 = 0._r8
	    end if
	    if( tausla(i,wn) < 500._r8 ) then
               expon1 = exp( -tausla(i,wn) )
	    else
               expon1 = 0._r8
	    end if
            divisr = lam(i)*lam(i) - 1._r8/(mu2(i)*mu2(i))
            temp = max( eps,abs(divisr) )
            divisr = 1._r8/sign( temp,divisr )
            up = om*pifs*((gam1 - 1._r8/mu2(i))*gam3 + gam4*gam2)*divisr
            dn = om*pifs*((gam1 + 1._r8/mu2(i))*gam4 + gam2*gam3)*divisr
!-----------------------------------------------------------------------------
! 	... cup and cdn are when tau is equal to zero
!           cuptn and cdntn are when tau is equal to taun
!-----------------------------------------------------------------------------
            cup(i,wn) = up*expon0
            cdn(i)    = dn*expon0
            cuptn(i,wn) = up*expon1
            cdntn(i,wn) = dn*expon1
         end do level_loop

!-----------------------------------------------------------------------------
!	... set up matrix
!           ssfc = pg 16,292 equation 37  where pi fs is one (unity).
!-----------------------------------------------------------------------------
	if( tausla(pver,wn) < 500._r8 ) then
           ssfc = rsfc(wn)*mu*exp( -tausla(pver,wn) )*pifs
	else
           ssfc = 0._r8
	end if

!-----------------------------------------------------------------------------
! 	... the following are from pg 16,292  equations 39 - 43.
!           set up first row of matrix:
!-----------------------------------------------------------------------------
        a(1) = 0._r8
        b(1) = e1(1,wn)
        d(1) = -e2(1,wn)
        e(1) = fdn0 - cdn(1)

!-----------------------------------------------------------------------------
! 	... set up odd rows 3 thru (mrows - 1):
!-----------------------------------------------------------------------------
        a(3:mrows-1:2) = e2(1:pverm,wn)*e3(1:pverm,wn) - e4(1:pverm,wn)*e1(1:pverm,wn)
        b(3:mrows-1:2) = e1(1:pverm,wn)*e1(2:pver,wn) - e3(1:pverm,wn)*e3(2:pver,wn)
        d(3:mrows-1:2) = e3(1:pverm,wn)*e4(2:pver,wn) - e1(1:pverm,wn)*e2(2:pver,wn)
        e(3:mrows-1:2) = e3(1:pverm,wn)*(cup(2:pver,wn) - cuptn(1:pverm,wn)) + e1(1:pverm,wn)*(cdntn(1:pverm,wn) - cdn(2:pver))

!-----------------------------------------------------------------------------
! 	... set up even rows 2 thru (mrows - 2): 
!-----------------------------------------------------------------------------
        a(2:mrows-2:2) = e2(2:pver,wn)*e1(1:pverm,wn) - e3(1:pverm,wn)*e4(2:pver,wn)
        b(2:mrows-2:2) = e2(1:pverm,wn)*e2(2:pver,wn) - e4(1:pverm,wn)*e4(2:pver,wn)
        d(2:mrows-2:2) = e1(2:pver,wn)*e4(2:pver,wn) - e2(2:pver,wn)*e3(2:pver,wn)
        e(2:mrows-2:2) = (cup(2:pver,wn) - cuptn(1:pverm,wn))*e2(2:pver,wn) - (cdn(2:pver) - cdntn(1:pverm,wn))*e4(2:pver,wn)

!-----------------------------------------------------------------------------
! 	... set up last row of matrix at mrows:
!-----------------------------------------------------------------------------
        a(mrows) = e1(pver,wn) - rsfc(wn)*e3(pver,wn)
        b(mrows) = e2(pver,wn) - rsfc(wn)*e4(pver,wn)
        d(mrows) = 0._r8
        e(mrows) = ssfc - cuptn(pver,wn) + rsfc(wn)*cdntn(pver,wn)

	sub(wn,1:mrows)   = a(1:mrows)
	main(wn,1:mrows)  = b(1:mrows)
	super(wn,1:mrows) = d(1:mrows)
	y(wn,1:mrows)     = e(1:mrows)
      end do wave_loop

!-----------------------------------------------------------------------------
! 	... solve the system
!-----------------------------------------------------------------------------
      call tridec( nw, mrows, sub, main, super )
      call trislv( nw, mrows, sub, main, super, y )

!-----------------------------------------------------------------------------
!	... unfold solution of matrix, compute output fluxes
!-----------------------------------------------------------------------------
      do wn = 1,nw
!-----------------------------------------------------------------------------
! 	... the following equations are from pg 16,291  equations 31 & 32
!-----------------------------------------------------------------------------
	 e(:mrows) = y(wn,:mrows)
	 if( tausla(0,wn) < 500._r8 ) then
            radfld(1,wn) = 2._r8*(fdn0 +  e(1)*e3(1,wn) - e(2)*e4(1,wn) + cup(1,wn)) + exp( -tausla(0,wn) )
	 else
            radfld(1,wn) = 2._r8*(fdn0 +  e(1)*e3(1,wn) - e(2)*e4(1,wn) + cup(1,wn))
	 end if
	 where( tausla(1:pver,wn) < 500._r8 )
	    cdn(1:pver) = exp( -tausla(1:pver,wn) )
	 elsewhere
	    cdn(1:pver) = 0._r8
	 endwhere
         radfld(2:pverp,wn) = 2._r8*(e(1:mrows-1:2)*(e3(1:pver,wn) + e1(1:pver,wn)) &
                            + e(2:mrows:2)*(e4(1:pver,wn) + e2(1:pver,wn))       &
                            + cdntn(1:pver,wn) + cuptn(1:pver,wn)) + cdn(1:pver)
      end do

      end subroutine ps2str

      end module mo_ps2str
