
      module mo_aurora
!-----------------------------------------------------------------------
!
! Auroral oval parameterization. See reference:
! R.G. Roble, E.C. Ridley
! An auroral model for the NCAR thermospheric general circulation model (TGCM)
! Annales Geophysicae,5A, (6), 369-382, 1987. 
!
! The aurora oval is a circle in auroral circle coordinates.  Auroral circle
!  coordinates are offset from magnetic coordinates by offa degrees (radians)
!  towards 0 MLT and by dskofa degrees (radians) towards dusk (18 MLT).
! The aurora assumes a Maxwellian in energy, so that the characteristic
!  energy is half of the mean energy (or mean energy = 2*alfa, where alfa
!  is the characteristic energy).  The Maxwellian is approximated in the
!  aion and bion subroutines.
! The aurora oval is assumed to be a Gaussian in auroral latitude, with
!  peak values on the day (=1) and night (=2) sides that change from one to
!  the other using cosines of the auroral longitude coordinate.
! There is provision for a low energy (~75 eV) aurora at the location of the
!  regular (~1-6 keV) aurora in order to simulate the energy flux found
!  at higher altitudes that is non-Maxwellian, but the flux is usually
!  set to zero (1.e-80).
! There is provision for a proton (MeV) aurora, but the flux is usually
!  set to zero (1.e-20).
! The drizzle is a constant low energy electron flux over the polar cap,
!  which goes to 1/e over twice the half-width of the aurora at the
!  radius of the aurora.
! The cusp is a low energy electron flux centered over the dayside convection
!  entrance at phid at the convection reversal boundary theta0.  The cusp
!  falls off over 5 degrees in latitude and over 20 degrees in longitude
!  to 1/e values of the peak at the center.
! 1.e-20 and 1.e-80 are used to give a near zero answer.
!
! The polar drizzle and cusp electron energies are low, and soft particles
!  have great influence on the over-all thermospheric and ionospheric
!  structure, especially on the electron density profiles at mid-latitudes
!  and in winter since low energy electrons produce ionization at high
!  altitudes where loss rates are very low.  (Comment by Wenbin Wang.)
! The original energies for drizzle and cusp were alfad=0.75, alfac=0.5 keV.
! The original guess at energy fluxes were: ed=0.1+2.0*power/100.,ec=0.1+0.9*power/100.
! The next guess at energy fluxes were: ed=0.01+0.2*power/100., ec=0.01+0.09*power/100.
! The values below reflect higher estimates for the electron energy (lower alt)
!
! Calling sequence (all subs in mo_aurora, mo_aurora.F):
!   1) sub aurora_cons called once per time step from advance.
!   2) sub aurora called from dynamics, inside parallel latitude scan.
!   3) subs aurora_cusp and aurora_heat called from sub aurora.
!   4) sub aurora_ions called from sub aurora. 
!
!-----------------------------------------------------------------------

      use shr_kind_mod,  only: r8 => shr_kind_r8
      use mo_constants,  only: pi, &
                               avo => avogadro, &
                               boltz_cgs, &
                               gask => rgas_cgs
      use cam_logfile,   only: iulog
      use spmd_utils,    only: masterproc

      implicit none

      interface aurora
         module procedure aurora_prod
         module procedure aurora_hrate
      end interface

      save

      integer, parameter  :: isouth = 1
      integer, parameter  :: inorth = 2

      ! g = 8.7 m/s^2? Because this is 400 km up?
      real(r8), parameter :: grav   = 870._r8          ! (cm/s^2)

      integer  :: lev1 = 1
      real(r8) :: twopi
      real(r8) :: rmass_o1
      real(r8) :: rmass_o2
      real(r8) :: rmass_n2
      real(r8) :: rmassinv_o1
      real(r8) :: rmassinv_o2
      real(r8) :: rmassinv_n2
      real(r8) :: dtr

!-----------------------------------------------------------------------
! 	... polar drizzle parameters:
!   alfad: Characteristic Maxwellian energy of drizzle electrons (keV)
!   ed   : Column energy input of drizzle electrons (ergs/cm**2/s)
!   fd   : Electron particle flux of drizzle electrons (particles/cm**2/s)
!-----------------------------------------------------------------------
      real(r8), parameter :: alfad = 2.0_r8, &
                             ed    = 0.5_r8     
      real(r8) :: fd                     ! set in sub aurora_ions

!-----------------------------------------------------------------------
! 	... polar cusp parameters:
!   alfac: Characteristic Maxwellian energy of polar cusp electons (keV)
!   ec   : Column energy input of polar cusp electrons (ergs/cm**2/s)
!   fc   : Electron particle flux of polar cusp electrons (particles/cm**2/s)
!-----------------------------------------------------------------------
      real(r8), parameter :: alfac = 1.0_r8, &
                             ec    = 0.5_r8
      real(r8) :: fc                     ! set in sub aurora_ions

!-----------------------------------------------------------------------
! e1: Peak energy flux in noon sector of the aurora (ergs/cm**2/s)
! e2: Peak energy flux in midnight sector of the aurora (ergs/cm**2/s)
! h1: Gaussian half-width of the noon auroral oval in degrees
! h2: Gaussian half-width of the midnight auroral oval in degrees
!-----------------------------------------------------------------------
      real(r8) :: &
        e1, e2, &                        ! set in sub aurora_cons (function of hem power)
        h1, h2                           ! set in sub aurora_cons (function of hem power)

!-----------------------------------------------------------------------
! 	... solar proton parameters for the polar cap (time-gcm only)
!   alfa_sp: Characteristic Maxwellian energy of solar protons (MeV) (was alfad2)
!   e_sp   : Column energy input of solar protons (ergs/cm**2/s) (was ed2)
!   flx_sp : e_sp/1.602e-6, for input to sub bion                (was fd2)
! Add solar protons to ionization if add_sproton is true (time-gcm only)
!-----------------------------------------------------------------------
      logical :: add_sproton = .false.
      real(r8), parameter :: &
        alfa_sp = 10._r8, &
        e_sp    = 1.e-20_r8
      real(r8) :: flx_sp

!-----------------------------------------------------------------------
! 	... high energy electron parameters in the auroral oval (time-gcm only):
!-----------------------------------------------------------------------
      logical :: add_helectron = .false.
      real(r8), parameter :: alfa30 = 40._r8, &                  ! Characteristic energy of auroral electrons (Mev)
                             e30    = .05_r8                     ! Column energy of auroral electrons (ergs/cm**2/s)

!-----------------------------------------------------------------------
! 	... additional auroral parameters
!-----------------------------------------------------------------------
      real(r8) :: &
        alfa0, &        ! average of noon and midnight characteristic Maxw energies
        ralfa,ralfa2, & ! difference ratios of characteristic energies
        rrote, &        ! clockwise rotation from noon of peak dayside energy flux (e1)
        rroth, &        ! clockwise rotation from noon of dayside h1 Gaussian half-width
        h0, &           ! average of noon and midnight Gaussian half-widths
        rh, &           ! difference ratio of half-widths (rh=(h2-h1)/(h2+h1))
        e0,e20, &       ! e0 = average of noon and midnight electrons
        ree,re2, &      ! difference ratios of peak energy fluxes (ree=(e2-e1)/(e2+e1))
        alfa20          ! average of noon and midnight char energies for high alt aurora
      real(r8) :: &
        theta0(2), &    ! convection reversal boundary in radians
        offa(2), &      ! offset of oval towards 0 MLT relative to magnetic pole (rad)
        dskofa(2), &    ! offset of oval in radians towards 18 MLT (f(By))
        phid(2), &      ! dayside convection entrance in MLT converted to radians (f(By))
        rrad(2)         ! radius of auroral circle in radians
      real(r8) :: ctpoten    ! cross-cap potential (kV)
      real(r8) :: byimf      ! BY component of IMF (nT)


      private
      public :: aurora_inti, aurora_timestep_init, aurora
      public :: aurora_register

      logical :: has_ions = .false.
      integer :: indxAIPRS = -1

      contains

        
      !----------------------------------------------------------------------
      !----------------------------------------------------------------------
      subroutine aurora_register
        use ppgrid,       only : pver,pcols
        use physics_buffer, only : pbuf_add_field, dtype_r8

        ! add ionization rates to phys buffer for waccmx ionosphere module

        call pbuf_add_field('AurIPRateSum' , 'physpkg', dtype_r8, (/pcols,pver/), indxAIPRS)     ! Sum of ion auroral production rates for O2

      endsubroutine aurora_register

      subroutine aurora_inti
!-----------------------------------------------------------------------
! 	... initialize aurora module
!-----------------------------------------------------------------------

      use ppgrid,       only : pver
      use pmgrid,       only : plev, plevp
      use constituents, only : cnst_get_ind, cnst_mw
      use chem_mods,    only : adv_mass
      use ref_pres,     only : pref_mid
      use mo_chem_utls, only : get_spc_ndx
      use cam_history,  only : addfld, phys_decomp

      implicit none

!-----------------------------------------------------------------------
! 	... local variables
!-----------------------------------------------------------------------
      integer             :: k, m
      real(r8), parameter :: e       = 1.e-10_r8
      real(r8), parameter :: convert = 3.1211e8_r8

      real(r8) :: plb
      real(r8) :: alfa_1, alfa_2, alfa21, alfa22
      real(r8) :: e21, e22

      integer :: op_ndx,o2p_ndx,np_ndx,n2p_ndx,e_ndx

      op_ndx   = get_spc_ndx( 'Op' )
      o2p_ndx  = get_spc_ndx( 'O2p' )
      np_ndx   = get_spc_ndx( 'Np' )
      n2p_ndx  = get_spc_ndx( 'N2p' )
      e_ndx    = get_spc_ndx( 'e' )

      has_ions = op_ndx > 0 .and. o2p_ndx > 0 .and. np_ndx > 0 .and. n2p_ndx > 0 .and. e_ndx > 0

      if (.not. has_ions) return

!-----------------------------------------------------------------------
!	... initialize module variables
!-----------------------------------------------------------------------
      twopi  = 2._r8*pi
      dtr    = pi/180._r8

!-----------------------------------------------------------------------
!	... set molecular weights
!-----------------------------------------------------------------------
      call cnst_get_ind( 'O2', m )
      rmass_o2    = cnst_mw(m)
      rmassinv_o2 = 1._r8/rmass_o2
      call cnst_get_ind( 'O', m )
      rmass_o1    = cnst_mw(m)
      rmassinv_o1 = 1._r8/rmass_o1
      call cnst_get_ind( 'N', m )
      rmass_n2    = 2._r8*cnst_mw(m)
      rmassinv_n2 = 1._r8/rmass_n2

      offa(isouth)   = 4.3_r8*dtr
      offa(inorth)   = 3.7_r8*dtr
      phid(isouth)   = 0._r8
      phid(inorth)   = 0._r8
      alfa_1         = 2._r8
      alfa_2         = 2._r8

!-----------------------------------------------------------------------
! Values from 10/05/94 HPI estimates (50% or more higher than old estimates):
!     alfa_1 = amin1(1.5,1.25+0.05*plevel)
!     alfa_2 = 1.2 + 0.095*plevel
!-----------------------------------------------------------------------
      alfa0  = 0.5_r8*(alfa_1 + alfa_2)
      ralfa  = (alfa_2 - alfa_1) / (alfa_1 + alfa_2 + e)
      alfa21 = 0.075_r8
      alfa22 = 0.075_r8
      alfa20 = 0.5_r8 * (alfa21 + alfa22)
      ralfa2 = (alfa22 - alfa21) / (alfa21 + alfa22 + e)
      e21    = 1.e-80_r8
      e22    = 1.e-80_r8
      e20    = 0.5_r8 * (e21 + e22)
      re2    = (e22 - e21) / (e21 + e22)

!-----------------------------------------------------------------------
! Set cusp and drizzle parameters:
! (conversion between particle number density and characteristic
!  energy and column energy input)
!-----------------------------------------------------------------------
      fc = convert * ec / alfac
      fd = convert * ed / alfad

!-----------------------------------------------------------------------
! Solar proton flux:
!-----------------------------------------------------------------------
      flx_sp = e_sp/1.602e-6_r8

!-----------------------------------------------------------------------
! 	... set auroral lower bndy index
!-----------------------------------------------------------------------
      plb = 5.e-4_r8*exp( 7._r8 ) * .1_r8             ! Pa
      do k = 1,pver
	 if( pref_mid(k) >= plb ) then
	    lev1 = k-1
	    exit
	 end if
      end do

      if (masterproc) write(iulog,*) ' '
      if (masterproc) write(iulog,*) 'aurora_inti: aurora will go down to lev,p = ',lev1,pref_mid(lev1)
      if (masterproc) write(iulog,*) ' '

!-----------------------------------------------------------------------
! Report to stdout:
!-----------------------------------------------------------------------
#ifdef AURORA_DIAGS
        write(iulog,"(/,'aurora_cons:')")
        write(iulog,"('  cusp:    alfac=',f8.3,' ec=',f8.3,' fc=',e10.4)") &
          alfac,ec,fc
        write(iulog,"('  drizzle: alfad=',f8.3,' ed=',f8.3,' fd=',e10.4)") &
          alfad,ed,fd
        write(iulog,"('  half-widths = h1,h2=',2f10.3)") h1,h2
        write(iulog,"('  energy flux = e1,e2=',2f10.3)") e1,e2
        write(iulog,"('  add_sproton = ',l1)") add_sproton
        if( add_sproton ) then
           write(iulog,"('  solar protons: alfa_sp=',f8.3, &
          ' e_sp = ',e10.4,' flx_sp = ',e10.4)") alfa_sp,e_sp,flx_sp
        end if
        if( add_helectron ) then
          write(iulog,"('  high-energy electrons: alfa30 = ',f8.3, &
          ' e30 = ',f8.3)") alfa30,e30
        end if
        write(iulog,"(' ')")
#endif

        call addfld( 'QSUM', '/s ', pver, 'I', 'total ion production', phys_decomp )

      end subroutine aurora_inti

      subroutine aurora_timestep_init
!-----------------------------------------------------------------------
! 	... per timestep initialization
!-----------------------------------------------------------------------

      use mag_parms,  only : get_mag_parms
      use spmd_utils, only : masterproc

!-----------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------
      real(r8) :: power, plevel
      real(r8) :: roth, rote, rcp, rhp
      real(r8) :: arad

      if (.not. has_ions) return

!-----------------------------------------------------------------------
!	... get hemispheric power
!-----------------------------------------------------------------------
      call get_mag_parms( by = byimf, hpower = power, ctpoten = ctpoten )
#ifdef AURORA_DIAGS
      if( masterproc ) then
         write(iulog,*) '----------------------------------------'
         write(iulog,*) 'aurora_timestep_init: by,power,ctpoten = ',byimf,power,ctpoten
         write(iulog,*) '----------------------------------------'
      end if
#endif

      if( power >= .01_r8 ) then
         plevel = 2.09_r8*log( power )
      else
         plevel = 0._r8
      end if

! fvitt -- moved the calc of h1, h2, rh, and h0 from aurora_inti
!          This was done to for bit-for-bit restarts.
!          These aurora oval dimension quantities power dependent and should be updated.
!-----------------------------------------------------------------------
! h1 = Gaussian half-width of the noon auroral oval in degrees
! h2 = Gaussian half-width of the midnight auroral oval in degrees
!-----------------------------------------------------------------------
!      h1 = 3._r8
!      h2 = 10._r8
! modified by LQIAN, 2007
! produce realistic oval compared to NOAA empirical auroral oval and TIMED/GUVI
! h1 formula given by Wenbin base on POLARVIS image;
! h2 formula based on Emery et al original auroral parameterization report
      h1 = min(2.35_r8, 0.83_r8 + 0.33_r8*plevel)
      h2 = 2.5_r8+0.025_r8*max(power,55._r8)+0.01_r8*min(0._r8,power-55._r8)
!-----------------------------------------------------------------------
! Values from corrections to Emery et al Parameterization report:
!     h1 = amin1(2.35, 0.83 + 0.33*plevel)
!     h2 = 2.87 + 0.15*plevel
!-----------------------------------------------------------------------

      rh = (h2 - h1) / (h1 + h2)
      h0     = 0.5_r8 * (h1 + h2) * dtr

      theta0(isouth) = (-3.80_r8 + 8.48_r8*(ctpoten**.1875_r8))*dtr
      theta0(inorth) = theta0(isouth)
      dskofa(isouth) = (-1.26_r8 + 0.15_r8 * byimf)*dtr
      dskofa(inorth) = dskofa(isouth)
      roth           = (12.18_r8 - 0.89_r8 * plevel)
      rote           = ( 2.62_r8 - 0.55_r8 * plevel)
      rroth          = roth * dtr
      rrote          = rote * dtr

!-----------------------------------------------------------------------
! e1 = energy flux in the noon sector of the aurora (ergs/cm**2/s)
! e2 = energy flux in the midnight sector of the aurora (ergs/cm**2/s)
!-----------------------------------------------------------------------
!      e1 = (0.5_r8 + 0.15_r8 * power)
!      e2 = (1.5_r8 + 0.25_r8 * power)
! modified by LQIAN, 2008
! produce realistic oval compared to NOAA empirical auroral oval and TIMED/GUVI
! e1 formula given by Wenbin base on POLARVIS image;
! e2 formula based on Emery et al original auroral parameterization report
      e1 = max(0.50_r8, -2.15_r8 + 0.62_r8*plevel)
      e2=1._r8+0.11_r8*power
!-----------------------------------------------------------------------
! Values from corrections to Emery et al Parameterization report:
!-----------------------------------------------------------------------
      e0  = 0.5_r8 * (e1 + e2)
      ree = (e2 - e1) / (e1 + e2)

      rhp          = 14.20_r8 + 0.96_r8*plevel
      rcp          = -0.43_r8 + 9.69_r8 * (ctpoten**.1875_r8)
      arad         = max( rcp,rhp )
      rrad(isouth) = arad*dtr 
      rrad(inorth) = arad*dtr

      end subroutine aurora_timestep_init

      subroutine aurora_prod( tn, o2, o1, mbar, rlats, &
                              qo2p, qop, qn2p, qnp, pmid, &
                              lchnk, calday,  ncol, rlons, pbuf )
!-----------------------------------------------------------------------
! 	... auroral parameterization driver
!-----------------------------------------------------------------------

      use mo_apex,     only : alatm, alonm                      ! magnetic latitude,longitude grid (radians)
      use ppgrid,      only : pcols, pver
      use cam_history, only : outfld
      use physics_buffer,only: physics_buffer_desc

      implicit none

!-----------------------------------------------------------------------
! 	... dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) ::  &
       ncol, &                           ! column count
       lchnk                             ! chunk index
      real(r8), intent(in) :: &
        calday                           ! calendar day of year
      real(r8), intent(in) :: &
        tn(pcols,pver), &                ! neutral gas temperature (K)
        o2(ncol,pver), &                 ! O2 concentration (kg/kg)
        o1(ncol,pver), &                 ! O concentration (kg/kg)
        mbar(ncol,pver)                  ! mean molecular weight (g/mole)
      real(r8), intent(in) :: &
        pmid(pcols,pver)                 ! midpoint pressure (Pa)
      real(r8), intent(in) :: &
        rlats(ncol), &                   ! column latitudes (radians)
        rlons(ncol)
      real(r8), intent(out) :: &
        qo2p(ncol,pver), &               ! o2+ production
        qop(ncol,pver), &                ! o+ production
        qn2p(ncol,pver), &               ! n2+ production
        qnp(ncol,pver)                   ! n+ production

      type(physics_buffer_desc),pointer :: pbuf(:)

!-----------------------------------------------------------------------
! 	... local variables
!-----------------------------------------------------------------------
      integer  :: i, k
      integer  :: hemis(ncol)
      real(r8) :: r2d
      real(r8) :: ofda, cosofa, sinofa, aslona
      real(r8) :: sunlons(ncol)                     ! sun's mag longitudes
      real(r8) :: dlat_aur(ncol)
      real(r8) :: dlon_aur(ncol)
      real(r8) :: colat(ncol)
      real(r8) :: sinlat(ncol)
      real(r8) :: coslat(ncol)
      real(r8) :: coslon(ncol)
      real(r8) :: sinlon(ncol)
      real(r8) :: alon(ncol)
      real(r8) :: cusp(ncol)
      real(r8) :: alfa(ncol)
      real(r8) :: alfa2(ncol)
      real(r8) :: alfa3(ncol)
      real(r8) :: flux(ncol)
      real(r8) :: flux2(ncol)
      real(r8) :: flux3(ncol)
      real(r8) :: drizl(ncol)
      real(r8) :: qteaur(ncol)                         ! for electron temperature
      logical  :: do_aurora(ncol)

      if (.not. has_ions) return

!-----------------------------------------------------------------------
! 	... initialize ion production
!-----------------------------------------------------------------------
      do k = 1,pver
        qo2p(:,k) = 0._r8
        qop(:,k)  = 0._r8
        qn2p(:,k) = 0._r8
        qnp(:,k)  = 0._r8
      end do

      r2d = 180._r8/pi

!-----------------------------------------------------------------------
! 	... output mag lons, lats
!-----------------------------------------------------------------------
      call outfld( 'ALONM', r2d*alonm(:ncol,lchnk), pcols, lchnk )
      call outfld( 'ALATM', r2d*alatm(:ncol,lchnk), pcols, lchnk )

!-----------------------------------------------------------------------
! 	... check latitudes, and return if all below 32.5 deg
!-----------------------------------------------------------------------
      do_aurora(:) = abs( rlats(:) ) > pi/6._r8
      if( all( .not. do_aurora(:) ) ) then
         return
      end if

!-----------------------------------------------------------------------
! 	... set sun location
!-----------------------------------------------------------------------
      call sunloc( calday, sunlons, lchnk, ncol )

      do i = 1,ncol
        if( do_aurora(i) ) then
          dlat_aur(i) = alatm(i,lchnk)
          dlon_aur(i) = alonm(i,lchnk) - sunlons(i)
          if( dlon_aur(i) > pi ) then
             dlon_aur(i) = dlon_aur(i) - twopi
          else if( dlon_aur(i) < -pi ) then
             dlon_aur(i) = dlon_aur(i) + twopi
          end if
          if( dlat_aur(i) > 0._r8 ) then
            hemis(i) = 2
          else
            hemis(i) = 1
          end if
!-----------------------------------------------------------------------
! 	... find auroral circle coordinates
!-----------------------------------------------------------------------
            ofda      = sqrt( offa(hemis(i))**2 + dskofa(hemis(i))**2)
            cosofa    = cos( ofda )
            sinofa    = sin( ofda )
            aslona    = asin( dskofa(hemis(i))/ofda )
            sinlat(i) = sin( abs( dlat_aur(i) ) )
            coslat(i) = cos( dlat_aur(i) )
            sinlon(i) = sin( dlon_aur(i) + aslona )
            coslon(i) = cos( dlon_aur(i) + aslona )
            colat(i)  = acos( cosofa*sinlat(i) - sinofa*coslat(i)*coslon(i))
            alon(i)   = mod( atan2( sinlon(i)*coslat(i),sinlat(i)*sinofa &
                                    + cosofa*coslat(i)*coslon(i) ) - aslona + 3._r8*pi,twopi) - pi
        end if
      end do
#ifdef AURORA_DIAGS
      write(iulog,*) '-----------------------------------------------------'
      write(iulog,*) 'aurora: diagnostics for lchnk = ',lchnk
      write(iulog,*) '        geo lats'
      write(iulog,'(1p,5g15.7)') r2d*rlats(:ncol)
      write(iulog,*) '        geo lons'
      write(iulog,'(1p,5g15.7)') r2d*rlons(:ncol)
      write(iulog,*) '        mag lats'
      write(iulog,'(1p,5g15.7)') r2d*dlat_aur(:ncol)
      write(iulog,*) '        mag lons'
      write(iulog,'(1p,5g15.7)') r2d*alonm(:ncol,lchnk)
      write(iulog,*) '        mag table lons'
      write(iulog,'(1p,5g15.7)') r2d*dlon_aur(:ncol)
      write(iulog,*) '        sunlons'
      write(iulog,'(1p,5g15.7)') r2d*sunlons(:ncol)
      write(iulog,*) '     min,max mag lons = ',r2d*minval(alonm(:ncol,lchnk)),r2d*maxval(alonm(:ncol,lchnk))
      write(iulog,*) '-----------------------------------------------------'
#endif

!-----------------------------------------------------------------------
! 	... make cusp
!-----------------------------------------------------------------------
      call aurora_cusp( cusp, do_aurora, hemis, colat, alon, ncol )

!-----------------------------------------------------------------------
! 	... make alfa, flux, and drizzle
!-----------------------------------------------------------------------
      call aurora_heat( flux, flux2, flux3, alfa, alfa2, &
                        alfa3, qteaur, drizl, do_aurora, hemis, &
                        alon, colat, ncol )

!-----------------------------------------------------------------------
! 	... auroral additions to ionization rates
!-----------------------------------------------------------------------
      call aurora_ions( drizl, cusp, alfa, alfa2, alfa3, &
                        flux, flux2, flux3, tn, o2, &
                        o1, mbar, qo2p, qop, qn2p, &
                        qnp, pmid, do_aurora, ncol, lchnk, pbuf )

      end subroutine aurora_prod

      subroutine aurora_hrate( tn, o2, o1, mbar, rlats, &
                               aur_hrate, cpair, pmid, lchnk, calday, &
                               ncol, rlons )
!-----------------------------------------------------------------------
! 	... auroral parameterization driver
!-----------------------------------------------------------------------

      use mo_apex, only : alatm, alonm                      ! magnetic latitude,longitude grid (radians)
      use ppgrid,  only : pcols, pver

      implicit none

!-----------------------------------------------------------------------
! 	... dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) ::  &
       ncol, &                           ! column count
       lchnk                             ! chunk index
      real(r8), intent(in) :: &
        calday                           ! calendar day of year
      real(r8), intent(in) :: &
        tn(pcols,pver), &                ! neutral gas temperature (K)
        o2(ncol,pver), &                 ! O2 concentration (kg/kg)
        o1(ncol,pver), &                 ! O concentration (kg/kg)
        mbar(ncol,pver)                  ! mean molecular weight (g/mole)
      real(r8), intent(in) :: &
        cpair(ncol,pver)                 ! specific heat capacity (J/K/kg)
      real(r8), intent(in) :: &
        pmid(pcols,pver)                 ! midpoint pressure (Pa)
      real(r8), intent(in) :: &
        rlats(ncol), &                   ! column latitudes (radians)
        rlons(ncol)
      real(r8), intent(out) :: &
        aur_hrate(ncol,pver)             ! auroral heating rate

!-----------------------------------------------------------------------
! 	... local variables
!-----------------------------------------------------------------------
      real(r8), parameter :: aur_therm     = 807._r8
      real(r8), parameter :: jkcal         = 4184._r8
      real(r8), parameter :: aur_heat_eff  = .05_r8
      real(r8), parameter :: aur_hconst    = 1.e3_r8*jkcal*aur_therm*aur_heat_eff

      integer  :: i, k
      integer  :: hemis(ncol)
      real(r8) :: r2d
      real(r8) :: ofda, cosofa, sinofa, aslona
      real(r8) :: sunlons(ncol)                     ! sun's mag longitudes
      real(r8) :: dlat_aur(ncol)
      real(r8) :: dlon_aur(ncol)
      real(r8) :: colat(ncol)
      real(r8) :: sinlat(ncol)
      real(r8) :: coslat(ncol)
      real(r8) :: coslon(ncol)
      real(r8) :: sinlon(ncol)
      real(r8) :: alon(ncol)
      real(r8) :: cusp(ncol)
      real(r8) :: alfa(ncol)
      real(r8) :: alfa2(ncol)
      real(r8) :: alfa3(ncol)
      real(r8) :: flux(ncol)
      real(r8) :: flux2(ncol)
      real(r8) :: flux3(ncol)
      real(r8) :: drizl(ncol)
      real(r8) :: qteaur(ncol)                         ! for electron temperature
      real(r8) :: qsum(ncol,pver)                      ! total ion production (1/s)
      logical  :: do_aurora(ncol)

!-----------------------------------------------------------------------
! 	... initialize ion production
!-----------------------------------------------------------------------
      do k = 1,pver
        aur_hrate(:,k) = 0._r8
      end do

      if (.not. has_ions) return

      r2d = 180._r8/pi

!-----------------------------------------------------------------------
! 	... check latitudes, and return if all below 32.5 deg
!-----------------------------------------------------------------------
      do_aurora(:) = abs( rlats(:) ) > pi/6._r8
      if( all( .not. do_aurora(:) ) ) then
         return
      end if

!-----------------------------------------------------------------------
! 	... set sun location
!-----------------------------------------------------------------------
      call sunloc( calday, sunlons, lchnk, ncol )

      do i = 1,ncol
        if( do_aurora(i) ) then
          dlat_aur(i) = alatm(i,lchnk)
          dlon_aur(i) = alonm(i,lchnk) - sunlons(i)
          if( dlon_aur(i) > pi ) then
             dlon_aur(i) = dlon_aur(i) - twopi
          else if( dlon_aur(i) < -pi ) then
             dlon_aur(i) = dlon_aur(i) + twopi
          end if
          if( dlat_aur(i) > 0._r8 ) then
            hemis(i) = 2
          else
            hemis(i) = 1
          end if
!-----------------------------------------------------------------------
! 	... find auroral circle coordinates
!-----------------------------------------------------------------------
            ofda      = sqrt( offa(hemis(i))**2 + dskofa(hemis(i))**2)
            cosofa    = cos( ofda )
            sinofa    = sin( ofda )
            aslona    = asin( dskofa(hemis(i))/ofda )
            sinlat(i) = sin( abs( dlat_aur(i) ) )
            coslat(i) = cos( dlat_aur(i) )
            sinlon(i) = sin( dlon_aur(i) + aslona )
            coslon(i) = cos( dlon_aur(i) + aslona )
            colat(i)  = acos( cosofa*sinlat(i) - sinofa*coslat(i)*coslon(i))
            alon(i)   = mod( atan2( sinlon(i)*coslat(i),sinlat(i)*sinofa &
                                    + cosofa*coslat(i)*coslon(i) ) - aslona + 3._r8*pi,twopi) - pi
        end if
      end do
#ifdef AURORA_DIAGS
      write(iulog,*) '-----------------------------------------------------'
      write(iulog,*) 'aurora: diagnostics for lchnk = ',lchnk
      write(iulog,*) '        geo lats'
      write(iulog,'(1p,5g15.7)') r2d*rlats(:ncol)
      write(iulog,*) '        geo lons'
      write(iulog,'(1p,5g15.7)') r2d*rlons(:ncol)
      write(iulog,*) '        mag lats'
      write(iulog,'(1p,5g15.7)') r2d*dlat_aur(:ncol)
      write(iulog,*) '        mag lons'
      write(iulog,'(1p,5g15.7)') r2d*alonm(:ncol,lchnk)
      write(iulog,*) '        mag table lons'
      write(iulog,'(1p,5g15.7)') r2d*dlon_aur(:ncol)
      write(iulog,*) '        sunlons'
      write(iulog,'(1p,5g15.7)') r2d*sunlons(:ncol)
      write(iulog,*) '     min,max mag lons = ',r2d*minval(alonm(:ncol,lchnk)),r2d*maxval(alonm(:ncol,lchnk))
      write(iulog,*) '-----------------------------------------------------'
#endif

!-----------------------------------------------------------------------
! 	... make cusp
!-----------------------------------------------------------------------
      call aurora_cusp( cusp, do_aurora, hemis, colat, alon, ncol )

!-----------------------------------------------------------------------
! 	... make alfa, flux, and drizzle
!-----------------------------------------------------------------------
      call aurora_heat( flux, flux2, flux3, alfa, alfa2, &
                        alfa3, qteaur, drizl, do_aurora, hemis, &
                        alon, colat, ncol )

!-----------------------------------------------------------------------
! 	... auroral additions to ionization rates
!-----------------------------------------------------------------------
      call total_ion_prod( drizl, cusp, alfa, alfa2, alfa3, &
                           flux, flux2, flux3, tn, o2, &
                           o1, mbar, qsum, pmid, do_aurora, &
                           ncol, lchnk )
!-----------------------------------------------------------------------
! 	... form auroral heating rate
!-----------------------------------------------------------------------
      do k = 1,pver
         aur_hrate(:,k) = aur_hconst * qsum(:,k) / (cpair(:,k) * mbar(:,k))
      end do

      end subroutine aurora_hrate

      subroutine aurora_cusp( cusp, do_aurora, hemis, colat, alon, ncol )
!-----------------------------------------------------------------------
! 	... calculate horizontal variation of polar cusp heating
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
! 	... dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in)     :: ncol
      integer, intent(in)     :: hemis(ncol)
      real(r8), intent(in)    :: colat(ncol)
      real(r8), intent(in)    :: alon(ncol)
      real(r8), intent(out)   :: cusp(ncol)
      logical, intent(in)     :: do_aurora(ncol)

!-----------------------------------------------------------------------
! 	... local variables
!-----------------------------------------------------------------------
      real(r8), parameter :: s5  =.08726646_r8, &
                             s20 =.34906585_r8

      where( do_aurora(:) )
         cusp(:) = (exp( -((theta0(hemis(:)) - colat(:))/s5)**2 ) &
                      + exp( -((pi - theta0(hemis(:)) - colat(:))/s5)**2) ) &
                        *exp( -(atan2( sin(alon(:) - phid(hemis(:))), cos(alon(:) - phid(hemis(:))) )/s20)**2 )
      elsewhere
         cusp(:) = 0._r8
      endwhere

      end subroutine aurora_cusp 

      subroutine aurora_heat( flux, flux2, flux3, alfa, alfa2, &
                              alfa3, qteaur, drizl, do_aurora, hemis, &
                              alon, colat, ncol )
!-----------------------------------------------------------------------
! 	... calculate alfa, flux, and drizzle
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
! 	... dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in)     :: ncol
      integer, intent(in)     :: hemis(ncol)
      real(r8), intent(in)    :: colat(ncol)
      real(r8), intent(in)    :: alon(ncol)
      real(r8), intent(inout) :: flux(ncol)
      real(r8), intent(inout) :: flux2(ncol)
      real(r8), intent(inout) :: flux3(ncol)
      real(r8), intent(inout) :: drizl(ncol)
      real(r8), intent(inout) :: qteaur(ncol)
      real(r8), intent(inout) :: alfa(ncol)
      real(r8), intent(inout) :: alfa2(ncol)
      real(r8), intent(inout) :: alfa3(ncol)
      logical, intent(in)     :: do_aurora(ncol)

!-----------------------------------------------------------------------
! 	... local variables
!-----------------------------------------------------------------------
      real(r8), dimension(ncol) :: &
        coslamda, &                             ! cos(angle from throat)
        halfwidth, &                            ! oval half-width
        wrk, &                                  ! temp wrk array
        dtheta                                  ! latitudinal variation (Gaussian)

!-----------------------------------------------------------------------
! Low-energy protons:
!
!     alfap0 = 0.5*(alfap1+alfap2)
!     e0p = 0.5*(pe1+pe2)
!
! coslamda  = cos(lamda)
! halfwidth = auroral half width
! dtheta    = colat-theta0(ihem)
! alfa      = electron energy
!-----------------------------------------------------------------------
      where( do_aurora(:) )
         coslamda(:) = cos( atan2( sin( alon(:) - rrote ),cos( alon(:) - rrote ) ) )
!-----------------------------------------------------------------------
! 	... auroral oval half-width (equation (1) in Roble,1987):
!-----------------------------------------------------------------------
         halfwidth(:) = h0*(1._r8 - rh*cos( atan2( sin(alon(:) - rroth),cos( alon(:) - rroth ) ) ) )
         dtheta(:)    = colat(:) - rrad(hemis(:))
      endwhere
!-----------------------------------------------------------------------
! 	... characteristic energy (equation (2) in Roble,1987):
!-----------------------------------------------------------------------
      if( alfa0 > .01_r8 ) then
         where( do_aurora(:) )
            alfa(:) =  alfa0*(1._r8 - ralfa*coslamda(:))
         endwhere
      else
         alfa(:) =  0._r8
      end if

      where( do_aurora(:) )
!-----------------------------------------------------------------------
! 	... flux, drizzle, alfa2, flux2
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 	... energy flux (equation (3) in Roble,1987):
!-----------------------------------------------------------------------
         wrk(:)   = exp( -(dtheta(:)/halfwidth(:))**2 )
         flux(:)  = e0*(1._r8 - ree*coslamda(:))*wrk(:) / (2._r8*alfa(:)*1.602e-9_r8)
         drizl(:) = exp( -((dtheta(:) + abs(dtheta(:)))/(2._r8*h0))**2 )
         alfa2(:) = alfa20*(1._r8 - ralfa2*coslamda(:))
         flux2(:) = e20*(1._r8 - re2*coslamda(:))*wrk(:) / (2._r8*alfa2(:)*1.602e-9_r8)
!-----------------------------------------------------------------------
! 	... alfa3, flux3 for high energy electrons:
!-----------------------------------------------------------------------
         alfa3(:) = alfa30
         flux3(:) = e30*wrk(:) / 1.602e-6_r8
!-----------------------------------------------------------------------
! 	... for electron temperature (used in settei):  
!-----------------------------------------------------------------------
         qteaur(:) = -7.e8_r8*wrk(:)
      endwhere

      end subroutine aurora_heat

      subroutine aurora_ions( drizl, cusp, alfa1, alfa2, alfa3, &
                              flux1, flux2, flux3, tn, o2, &
                              o1, mbar, qo2p, qop, qn2p, &
                              qnp, pmid, do_aurora, ncol, lchnk, pbuf )
!-----------------------------------------------------------------------
! 	... calculate auroral additions to ionization rates
!-----------------------------------------------------------------------

      use ppgrid,      only : pcols, pver
      use cam_history, only : outfld

      use physics_buffer,only: physics_buffer_desc, pbuf_set_field

      implicit none

!-----------------------------------------------------------------------
! 	... dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: ncol
      integer, intent(in) :: lchnk
      real(r8), intent(in), dimension(ncol) :: &
                             drizl, &
                             cusp, &
                             alfa1, &
                             alfa2, &
                             alfa3, &
                             flux1, &
                             flux2, &
                             flux3
      real(r8), dimension(pcols,pver), intent(in) :: &
                             tn, &                     ! midpoint neutral temperature (K)
                             pmid                      ! midpoint pressure (Pa)
      real(r8), dimension(ncol,pver), intent(in) :: &
                             o2, &                     ! midpoint o2 concentration (kg/kg)
                             o1, &                     ! midpoint o  concentration (kg/kg)
                             mbar                      ! mean molecular mass (g/mole)
      real(r8), dimension(ncol,pver), intent(inout) :: &
                             qo2p, &                   ! o2p prod from aurora (molecules/cm^3/s)
                             qop, &                    ! op prod from aurora (molecules/cm^3/s)
                             qn2p, &                   ! n2p prod from aurora (molecules/cm^3/s)
                             qnp                       ! np prod from aurora (molecules/cm^3/s)
      logical, intent(in) :: do_aurora(ncol)

      type(physics_buffer_desc),pointer :: pbuf(:)

!-----------------------------------------------------------------------
! 	... local variables
!-----------------------------------------------------------------------
      real(r8), parameter :: const0        = 1.e-20_r8

      integer  :: i, k
      real(r8), dimension(ncol) :: &
        p0ez, &
        press, &                                   ! pressure at interface levels (dyne/cm^2)
        tempi, &                                   ! temperature at interface levels (K)
        xalfa1, &
        xalfa2, &
        xcusp, &
        xdrizl, &                                  ! input to sub aion
        xalfa_sp, &
        xalfa3, &
        flux1_ion, &
        flux2_ion, &
        cusp_ion, &
        drizl_ion, &                               ! output from sub aion
        alfa1_ion, &
        alfa2_ion, &
        alfa3_ion, &                               ! output from sub aion
        alfasp_bion, &                             ! output from sub bion
        barm_t, &
        qsum, &
        denom, &
        p0ez_mbar, &
        tk_mbar, &
        barm, &
        falfa1, &
        falfa2, &
        fcusp, &
        fdrizl, &
        falfa_sp, &
        xn2, &
        falfa3
      real(r8), dimension(ncol) :: &
        qo2p_aur, &
        qop_aur, &
        qn2p_aur                                   ! auroral ionization for O2+, O+, N2+
      real(r8) :: qia(5)                           ! low energy proton source (not in use, 1/02)
      real(r8) :: wrk(ncol,pver)

      qia(:) = 0._r8
      wrk(:,:) = 0._r8

level_loop : &
      do k = 1,lev1
          where( do_aurora(:) )
             press(:ncol) = 10._r8*pmid(:ncol,k)              ! from Pa to dyne/cm^2
             tempi(:ncol) = tn(:ncol,k)
             barm(:)      = mbar(:,k)
             p0ez(:)      = (press(:)/(grav*4.e-6_r8))**.606_r8
             xalfa1(:)    = p0ez(:)/alfa1(:)
             xalfa2(:)    = p0ez(:)/alfa2(:)
             xcusp (:)    = p0ez(:)/alfac
             xdrizl(:)    = p0ez(:)/alfad

!-----------------------------------------------------------------------
! 	... initialize (whole array operations):
!-----------------------------------------------------------------------
             flux1_ion(:) = const0
             flux2_ion(:) = const0
             alfa1_ion(:) = const0
             alfa2_ion(:) = const0
             cusp_ion(:)  = const0
             drizl_ion(:) = const0
          endwhere
!-----------------------------------------------------------------------
! 	... auroral electrons
!-----------------------------------------------------------------------
          call aion( xalfa1, alfa1_ion, do_aurora, ncol )
          call aion( xalfa2, alfa2_ion, do_aurora, ncol )
          call aion( xcusp , cusp_ion, do_aurora, ncol  )
          call aion( xdrizl, drizl_ion, do_aurora, ncol )
          where( do_aurora(:) )
             falfa1(:) = alfa1(:)*flux1(:)  ! s7
             falfa2(:) = alfa2(:)*flux2(:)  ! s8
             fcusp (:) = cusp(:)*alfac*fc   ! s9
             fdrizl(:) = drizl(:)*alfad*fd  ! s10
             qsum(:)   = falfa1(:)*alfa1_ion(:) &    ! s7*s3
                       + falfa2(:)*alfa2_ion(:) &    ! s8*s4
                       + fcusp(:)*cusp_ion (:) &     ! s9*s5
                       + fdrizl(:)*drizl_ion(:)       ! s10*s6
          endwhere
!-----------------------------------------------------------------------
! 	... include solar protons if add_sproton is set, 
!           and high energy electrons if add_helectron is set
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 	... solar protons
!-----------------------------------------------------------------------
          if( add_sproton ) then
             if( flx_sp > 1.e-19_r8 ) then
                where( do_aurora(:) )
                   p0ez_mbar(:) = press(:)/(boltz_cgs*tempi(:)*barm(:))/avo
                   tk_mbar(:)   = gask*tempi(:)/(grav*barm(:))*p0ez_mbar(:)
                   xalfa_sp(:)  = ((tk_mbar(:)/.00271_r8)**.58140_r8)/alfa_sp
                endwhere
                call bion( xalfa_sp, alfasp_bion, do_aurora, ncol )
                where( do_aurora(:) )
                   falfa_sp(:) = drizl(:)*p0ez_mbar(:)*flx_sp*1.e6_r8/(tk_mbar(:)*35._r8)
                   qsum(:)     = qsum(:) + falfa_sp(:)*alfasp_bion(:)
                endwhere
             end if
          end if
!-----------------------------------------------------------------------
! 	... high energy electrons
!-----------------------------------------------------------------------
          if( add_helectron ) then
             if( e30 > 1.e-19_r8 ) then
                where( do_aurora(:) )
                   xalfa3(:) = p0ez(:)/alfa3(:)        ! alfa3(:)==alfa30
                endwhere
                call aion( xalfa3, alfa3_ion, do_aurora, ncol )
                where( do_aurora(:) )
                   falfa3(:) = alfa3(:)*flux3(:)  ! s13 (high energy electrons)
                   qsum(:)   = qsum(:) + falfa3(:)*alfa3_ion(:)
                endwhere
             end if
          end if

!-----------------------------------------------------------------------
! 	... form production
!-----------------------------------------------------------------------
          where( do_aurora(:) )
             barm_t(:) = grav*barm(:)/(35.e-3_r8*gask*tempi(:))
             qsum(:)   = qsum(:)*barm_t(:)               ! s1 = s1*s11
             wrk(:,k)  = qsum(:)
!-----------------------------------------------------------------------
! 	... denominator of equations (13-16) in Roble,1987.
!-----------------------------------------------------------------------
             xn2(:)   = max( (1._r8 - o2(:,k) - o1(:,k)),1.e-8_r8 )
             denom(:) = 0.92_r8*xn2(:)*rmassinv_n2 &
                      + 1.5_r8*o2(:,k) *rmassinv_o2 + 0.56_r8*o1(:,k) *rmassinv_o1
!-----------------------------------------------------------------------
! 	... production of O2+ (equation (15) in Roble,1987):
!-----------------------------------------------------------------------
             qo2p_aur(:) = qsum(:)*o2(:,k)/(rmass_o2*denom(:)) + qia(2)
!-----------------------------------------------------------------------
! 	... production of O+ (equation (16) in Roble,1987):
!-----------------------------------------------------------------------
             qop_aur(:) = qsum(:)*(.5_r8 *o2(:,k)*rmassinv_o2 &
                                   + .56_r8*o1(:,k)*rmassinv_o1)/denom(:) + qia(3)
!-----------------------------------------------------------------------
! 	... production of N2+ (equation (13) in Roble,1987)
!-----------------------------------------------------------------------
             qn2p_aur(:) = qsum(:)*.7_r8*xn2(:)/(rmass_n2*denom(:)) + qia(1)
             qo2p(:,k)   = qo2p(:,k) + qo2p_aur(:)
             qop(:,k)    = qop(:,k) + qop_aur(:)
             qn2p(:,k)   = qn2p(:,k) + qn2p_aur(:)
             qnp(:,k)    = qnp (:,k) + .22_r8/.7_r8 * qn2p_aur(:)
          endwhere
      end do level_loop

!--------------------------------------------------------------------------------------------
!  Save sum of auroral ion production rates to be accessed in ionosphere module in physics 
!  buffer using pointer for WACCM-X
!--------------------------------------------------------------------------------------------
      if (indxAIPRS>0) then
         call pbuf_set_field(pbuf, indxAIPRS,  wrk(1:ncol,1:pver))
      endif

      call outfld( 'QSUM', wrk, ncol, lchnk )

      end subroutine aurora_ions

      subroutine total_ion_prod( drizl, cusp, alfa1, alfa2, alfa3, &
                                 flux1, flux2, flux3, tn, o2, &
                                 o1, mbar, tpions, pmid, do_aurora, &
                                 ncol, lchnk )
!-----------------------------------------------------------------------
! 	... calculate auroral additions to ionization rates
!-----------------------------------------------------------------------

      use ppgrid,      only : pcols, pver
      use cam_history, only : outfld

      implicit none

!-----------------------------------------------------------------------
! 	... dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: ncol
      integer, intent(in) :: lchnk
      real(r8), intent(in), dimension(ncol) :: &
                             drizl, &
                             cusp, &
                             alfa1, &
                             alfa2, &
                             alfa3, &
                             flux1, &
                             flux2, &
                             flux3
      real(r8), dimension(pcols,pver), intent(in) :: &
                             tn, &                     ! midpoint neutral temperature (K)
                             pmid                      ! midpoint pressure (Pa)
      real(r8), dimension(ncol,pver), intent(in) :: &
                             o2, &                     ! midpoint o2 concentration (kg/kg)
                             o1, &                     ! midpoint o  concentration (kg/kg)
                             mbar                      ! mean molecular mass (g/mole)
      real(r8), dimension(ncol,pver), intent(inout) :: &
                             tpions                    ! total ion production (1/s)
      logical, intent(in) :: do_aurora(ncol)

!-----------------------------------------------------------------------
! 	... local variables
!-----------------------------------------------------------------------
      real(r8), parameter :: const0        = 1.e-20_r8

      integer  :: i, k
      real(r8), dimension(ncol) :: &
        p0ez, &
        press, &                                   ! pressure at interface levels (dyne/cm^2)
        tempi, &                                   ! temperature at interface levels (K)
        xalfa1, &
        xalfa2, &
        xcusp, &
        xdrizl, &                                  ! input to sub aion
        xalfa_sp, &
        xalfa3, &
        flux1_ion, &
        flux2_ion, &
        cusp_ion, &
        drizl_ion, &                               ! output from sub aion
        alfa1_ion, &
        alfa2_ion, &
        alfa3_ion, &                               ! output from sub aion
        alfasp_bion, &                             ! output from sub bion
        barm_t, &
        qsum, &
        denom, &
        p0ez_mbar, &
        tk_mbar, &
        barm, &
        falfa1, &
        falfa2, &
        fcusp, &
        fdrizl, &
        falfa_sp, &
        xn2, &
        falfa3
      real(r8), dimension(ncol) :: &
        qo2p_aur, &
        qop_aur, &
        qn2p_aur                                   ! auroral ionization for O2+, O+, N2+
      real(r8) :: qia(5)                           ! low energy proton source (not in use, 1/02)

      qia(:)      = 0._r8
      tpions(:,:) = 0._r8

level_loop : &
      do k = 1,lev1
          where( do_aurora(:) )
             press(:ncol) = 10._r8*pmid(:ncol,k)              ! from Pa to dyne/cm^2
             tempi(:ncol) = tn(:ncol,k)
             barm(:)      = mbar(:,k)
             p0ez(:)      = (press(:)/(grav*4.e-6_r8))**.606_r8
             xalfa1(:)    = p0ez(:)/alfa1(:)
             xalfa2(:)    = p0ez(:)/alfa2(:)
             xcusp (:)    = p0ez(:)/alfac
             xdrizl(:)    = p0ez(:)/alfad

!-----------------------------------------------------------------------
! 	... initiliaze (whole array operations):
!-----------------------------------------------------------------------
             flux1_ion(:) = const0
             flux2_ion(:) = const0
             alfa1_ion(:) = const0
             alfa2_ion(:) = const0
             cusp_ion(:)  = const0
             drizl_ion(:) = const0
          endwhere
!-----------------------------------------------------------------------
! 	... auroral electrons
!-----------------------------------------------------------------------
          call aion( xalfa1, alfa1_ion, do_aurora, ncol )
          call aion( xalfa2, alfa2_ion, do_aurora, ncol )
          call aion( xcusp , cusp_ion, do_aurora, ncol  )
          call aion( xdrizl, drizl_ion, do_aurora, ncol )
          where( do_aurora(:) )
             falfa1(:) = alfa1(:)*flux1(:)  ! s7
             falfa2(:) = alfa2(:)*flux2(:)  ! s8
             fcusp (:) = cusp(:)*alfac*fc   ! s9
             fdrizl(:) = drizl(:)*alfad*fd  ! s10
             qsum(:)   = falfa1(:)*alfa1_ion(:) &    ! s7*s3
                       + falfa2(:)*alfa2_ion(:) &    ! s8*s4
                       + fcusp(:)*cusp_ion (:) &     ! s9*s5
                       + drizl(:)*drizl_ion(:)       ! s10*s6
          endwhere
!-----------------------------------------------------------------------
! 	... include solar protons if add_sproton is set, 
!           and high energy electrons if add_helectron is set
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 	... solar protons
!-----------------------------------------------------------------------
          if( add_sproton ) then
             if( flx_sp > 1.e-19_r8 ) then
                where( do_aurora(:) )
                   p0ez_mbar(:) = press(:)/(boltz_cgs*tempi(:)*barm(:))/avo
                   tk_mbar(:)   = gask*tempi(:)/(grav*barm(:))*p0ez_mbar(:)
                   xalfa_sp(:)  = ((tk_mbar(:)/.00271_r8)**.58140_r8)/alfa_sp
                endwhere
                call bion( xalfa_sp, alfasp_bion, do_aurora, ncol )
                where( do_aurora(:) )
                   falfa_sp(:) = drizl(:)*p0ez_mbar(:)*flx_sp*1.e6_r8/(tk_mbar(:)*35._r8)
                   qsum(:)     = qsum(:) + falfa_sp(:)*alfasp_bion(:)
                endwhere
             end if
          end if
!-----------------------------------------------------------------------
! 	... high energy electrons
!-----------------------------------------------------------------------
          if( add_helectron ) then
             if( e30 > 1.e-19_r8 ) then
                where( do_aurora(:) )
                   xalfa3(:) = p0ez(:)/alfa3(:)        ! alfa3(:)==alfa30
                endwhere
                call aion( xalfa3, alfa3_ion, do_aurora, ncol )
                where( do_aurora(:) )
                   falfa3(:) = alfa3(:)*flux3(:)  ! s13 (high energy electrons)
                   qsum(:)   = qsum(:) + falfa3(:)*alfa3_ion(:)
                endwhere
             end if
          end if

!-----------------------------------------------------------------------
! 	... form production
!-----------------------------------------------------------------------
          where( do_aurora(:) )
             barm_t(:)   = grav*barm(:)/(35.e-3_r8*gask*tempi(:))
             tpions(:,k) = qsum(:)*barm_t(:)               ! s1 = s1*s11
          endwhere
      end do level_loop

      end subroutine total_ion_prod

      subroutine aion( si, so, do_aurora, ncol )
!-----------------------------------------------------------------------
! Calculates integrated f(x) needed for total auroral ionization.
! See equations (10-12) in Roble,1987.
! Coefficients for equation (12) of Roble,1987 are in variable cc 
! (revised since 1987):
! Uses the identity x**y = exp(y*ln(x)) for performance 
! (fewer (1/2) trancendental functions are required).
!------------------------------------------------------------------------

      implicit none

!------------------------------------------------------------------------
! 	... dummy arguments
!------------------------------------------------------------------------
      integer,  intent(in)  :: ncol
      real(r8), intent(in)  :: si(ncol)
      real(r8), intent(out) :: so(ncol)
      logical,  intent(in)  :: do_aurora(ncol)

!------------------------------------------------------------------------
! 	... local variables
!------------------------------------------------------------------------
      real(r8), parameter :: cc(8) = &
       (/ 3.2333134511131_r8 ,  2.5658873458085_r8 ,  2.2540957232641_r8 , &
          0.72971983372673_r8,  1.1069072431948_r8 ,  1.7134937681128_r8 , &
          1.8835442312993_r8 ,  0.86472135072090_r8 /)

      real(r8) :: xlog(ncol)

      where( do_aurora(:) )
         xlog(:) = log( si(:) )
         so(:)   = cc(1)*exp( cc(2)*xlog(:) - cc(3)*exp( cc(4)*xlog(:) ) ) &
                   + cc(5)*exp( cc(6)*xlog(:) - cc(7)*exp( cc(8)*xlog(:) ) )
      elsewhere
         so(:) = 0._r8
      endwhere

      end subroutine aion

      subroutine bion( si, so, do_aurora, ncol )
!-----------------------------------------------------------------------
! Calculates integrated f(x) needed for total auroral ionization.
! See equations (10-12) in Roble,1987.
! Use the identity x**y = exp(y*ln(x)) for performance 
! (fewer (1/2) trancendental functions are required).
!-----------------------------------------------------------------------

      implicit none

!------------------------------------------------------------------------
! 	... dummy arguments
!------------------------------------------------------------------------
      integer, intent(in)   :: ncol
      real(r8), intent(in)  :: si(ncol)
      real(r8), intent(out) :: so(ncol)
      logical,  intent(in)  :: do_aurora(ncol)

!------------------------------------------------------------------------
! 	... local variables
!------------------------------------------------------------------------
      real(r8), parameter :: cc(8) = &
        (/ 0.12718_r8, 4.9119_r8, 1.8429_r8, 0.99336_r8, 0.52472_r8, &
           1.5565_r8,  .85732_r8, 1.4116_r8 /)

      real(r8) :: xlog(ncol)

      where( do_aurora(:) )
         xlog(:) = log( si(:) )
         so(:)   = cc(1)*exp( cc(2)*xlog(:) - cc(3)*exp( cc(4)*xlog(:) ) ) &
                   + cc(5)*exp( cc(6)*xlog(:) - cc(7)*exp( cc(8)*xlog(:) ) )
      elsewhere
         so(:) = 0._r8
      endwhere

      end subroutine bion

      subroutine sunloc( calday, sunlons, lchnk, ncol )
!-----------------------------------------------------------------------
! 	... calculate sun's longitude in dipole coordinates, defining sunlon
!-----------------------------------------------------------------------
      use dyn_grid,   only : get_dyn_grid_parm, get_horiz_grid_d
      use spmd_utils, only : masterproc
      use  phys_grid, only : get_lat_all_p, get_rlat_all_p, get_rlon_all_p
      use  mo_apex,   only : glonm                                              ! magnetic longitude grid (radians)

      implicit none

!-----------------------------------------------------------------------
! 	... dummy arguments
!-----------------------------------------------------------------------
      integer,  intent(in)   :: ncol
      integer,  intent(in)   :: lchnk
      real(r8), intent(in)   :: calday  ! calendar day of year
      real(r8), intent(out)  :: sunlons(ncol)

!-----------------------------------------------------------------------
! 	... local variables
!-----------------------------------------------------------------------
      integer  :: col
      integer  :: latndx(ncol)
      integer  :: sunlon_ndx
      integer  :: sunlon_ndxp1
      integer :: plon, plat
      real(r8), allocatable :: clon(:)
      real(r8) :: wght1, wght2, dellon, r2d
      real(r8) :: rlats(ncol), rlons(ncol)

      r2d = 180._r8/pi
!-----------------------------------------------------------------------
! 	... sun's geographic coordinates
!-----------------------------------------------------------------------
      plon = get_dyn_grid_parm('plon')
      plat = get_dyn_grid_parm('plat')
      allocate(clon(plon))
      call get_horiz_grid_d(plon,clon_d_out=clon)


      dellon       = .5_r8 - (calday - int(calday))
      sunlon_ndx   = mod( nint( dellon*plon ) - 1,plon ) + 1
      if( sunlon_ndx < 1 ) then
         sunlon_ndx = plon + sunlon_ndx
      end if
      sunlon_ndxp1 = mod( sunlon_ndx,plon ) + 1
      wght2        = min( 1._r8, max( (dellon*twopi - clon(sunlon_ndx))*plon/twopi,0._r8 ) )
      wght1        = 1._r8 - wght2
      deallocate(clon)

!-----------------------------------------------------------------------      
!        ... get chunck latitudes
!-----------------------------------------------------------------------      
      call get_lat_all_p( lchnk, ncol, latndx )
      call get_rlat_all_p( lchnk, ncol, rlats )
      call get_rlon_all_p( lchnk, ncol, rlons )

      do col = 1,ncol
!        sunlons(col) = wght1*glonm(sunlon_ndx,latndx(col)) + wght2*glonm(sunlon_ndxp1,latndx(col))
	 sunlons(col) = wght1*glonm(sunlon_ndx,plat/2) + wght2*glonm(sunlon_ndxp1,plat/2)
      end do

#ifdef AURORA_DIAGS
      write(iulog,*) '-----------------------------------------------------'
      write(iulog,*) 'sunloc: diagnostics for lchnk = ',lchnk
      write(iulog,*) '        dellon,sunlon_ndx,sunlon_ndxp1,wght1,wght2'
      write(iulog,'(1p,g15.7,2i6,2g15.7)') dellon,sunlon_ndx,sunlon_ndxp1,wght1,wght2
      write(iulog,*) '   phys lons'
      write(iulog,'(1p,5g15.7)') r2d*rlons(:ncol)
      write(iulog,*) '  phys lats'
!     write(iulog,'(10i5)') latndx(:ncol)
      write(iulog,'(1p,5g15.7)') r2d*rlats(:ncol)
      write(iulog,*) '    mag lons'
      write(iulog,'(1p,5g15.7)') r2d*glonm(sunlon_ndx,latndx(:ncol))
      write(iulog,*) '    mag lons'
      write(iulog,'(1p,5g15.7)') r2d*glonm(sunlon_ndxp1,latndx(:ncol))
      write(iulog,*) '-----------------------------------------------------'
      write(iulog,"('sunloc: sunlons=',/,(e12.4))") sunlons
#endif

      end subroutine sunloc

      end module mo_aurora
