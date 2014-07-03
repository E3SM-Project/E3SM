
      module efield
!------------------------------------------------------------------------------ 
! description: calculates the electric potential for a given year,
!      day of year,UT, F10.7, B_z(K_p)
!     - low/midlatitudes electric potential is from an empirical model from
!       L.Scherliess ludger@gaim.cass.usu.edu
!     - high latitude electric potential is from Weimer96 model
!     - the transition zone is smoothed
!     - output is the horizontal global electric field in magnetic coordinates direction
!      at every magnetic local time grid point expressed in degrees (0 deg-0MLT; 360deg 24 MLT)
!
! input 
!      integer :: iday,     ! day number of year
!                 iyear     ! year
!      real(r8):: ut,       ! universal time 
!                 F10.7,    ! solar flux       (see ionosphere module)
!                 bz        ! component of IMF (see ionosphere module)
! output
!      real(r8) ::               &
!       ed1(0:nmlon,0:nmlat),    &  ! zonal electric field Ed1  [V/m] 
!       ed2(0:nmlon,0:nmlat)        ! meridional electric field Ed2/sin I_m  [V/m]  
!
! notes:
!
! - !to be done (commented out): input S_a F10.7/ Kp from WACCM and calculate B_z 
!    from these inputs
! - assume regular geomagnetic grid 
! - uses average year 365.24 days/year 30.6001 day/mo s. Weimer
! - get_tilt works only for iyear >= 1900
! - Weimer model 1996, Dan Weimer (not with the updates from B.Emery)
! - fixed parameters: B_z, B_y units nT  CHANGE THIS
!                     F10.7
! - we assume that the reference height is 300km for the emperical potential model
! - as a first approximation the electric field is constant in height
!   WATCH what is the upper boundary condition in WACCM
! - for all the calculation done here we set the reference height to the same 
!   value as in tiegcm (hr=130km)
! - 12/15/03 input value iseasav : replaced by day -> month and day of month
! - 12/15/03 S_aM calculated according to Scherliess draft paper and added
!   S_aM(corrected) = 90*(S_aM+1) to get variation in fig 1 Scherliess draft
!
! Author: A. Maute Dec 2003  am 12/30/03 
!------------------------------------------------------------------------------ 

      use shr_kind_mod,  only: r8 => shr_kind_r8
      use physconst,     only: pi
      use abortutils,    only: endrun
      use cam_logfile,   only: iulog
   
      implicit none

      public :: efield_init,  & ! interface routine
                get_efield      ! interface routine
      public :: ed1,          & ! zonal electric field Ed1  [V/m] 
                ed2,          & ! meridional electric field Ed2 [V/m] 
                potent,       & ! electric potential [V]
	        nmlon, nmlat, & ! dimension of mag. grid 
                dlatm, dlonm, & ! grid spacing of mag. grid 
	        ylonm, ylatm    ! magnetic longitudes/latitudes (deg)
      private

      integer ::  &
        iday,     &      ! day number of year
        iyear,    &      ! year
        iday_m,   &      ! day of month
        imo              ! month
      real(r8) ::  ut    ! universal time  

!------------------------------------------------------------------------------ 
! solar parameters
!------------------------------------------------------------------------------ 
      real(r8) ::   f107d           ! 10.7 cm solar flux
      real(r8) ::   by              ! By component of IMF [nT]
      real(r8) ::   bz              ! Bz component of IMF [nT]
!------------------------------------------------------------------------------ 
! mag. grid dimensions (assumed resolution of 2deg)
!------------------------------------------------------------------------------ 
      integer, parameter ::  &
      nmlon = 180,       &  ! mlon 
      nmlat = 90,        &  ! mlat
      nmlath= nmlat/2,   &  ! mlat/2
      nmlonh= nmlon/2,   &  ! mlon/2
      nmlonp1 = nmlon+1, &  ! mlon+1 
      nmlatp1 = nmlat+1     ! mlat+1

      real(r8) ::        &
        ylatm(0:nmlat),  &   ! magnetic latitudes (deg)
        ylonm(0:nmlon),  &   ! magnetic longitudes (deg)
        dlonm,	         &   ! delon lon grid spacing
        dlatm		     ! delat lat grid spacing

!------------------------------------------------------------------------------ 
! array on magnetic grid:    
!------------------------------------------------------------------------------ 
      real(r8) ::                &
        potent(0:nmlon,0:nmlat), &  ! electric potential   [V]  
        ed1(0:nmlon,0:nmlat),    &  ! zonal electric field Ed1  [V/m] 
        ed2(0:nmlon,0:nmlat)        ! meridional electric field Ed2/sin I_m  [V/m]  
       
      real(r8) :: &
       date,  & ! iyear+iday+ut
       day      ! iday+ut

      logical, parameter :: iutav=.false.   ! .true.  means UT-averaging 
                                            ! .false. means no UT-averaging
      real(r8), parameter ::          &
        v_sw = 400._r8 	                   ! solar wind velocity [km/s]

!------------------------------------------------------------------------------ 
! boundary for Weimer
!------------------------------------------------------------------------------ 
      real(r8), parameter :: bnd_wei = 44._r8 ! colat. [deg]
      integer :: nmlat_wei
      
!------------------------------------------------------------------------------ 
! flag for choosing factors for empirical low latitude model      
!------------------------------------------------------------------------------ 
      integer, parameter ::  iseasav = 0  ! flag for season 

!------------------------------------------------------------------------------ 
! constants:
!------------------------------------------------------------------------------ 
      real(r8), parameter ::        & 
     	r_e  =  6.371e6_r8,         &  ! radius_earth [m] (same as for apex.F90)
        h_r  = 130.0e3_r8,          &  ! reference height [m] (same as for apex.F90)
     	dy2yr= 365.24_r8,           &  ! day per avg. year used in Weimer
     	dy2mo= 30.6001_r8  	       ! day per avg. month used in Weimer

      real(r8) :: &
     	rtd ,     &    ! radians -> deg
     	dtr,      &    ! deg -> radians
     	sqr2,     &      
     	hr2rd,    &    ! pi/12 hrs
     	dy2rd,    &    ! 2*pi/365.24  average year
     	deg2mlt,  &    ! for mlon to deg
     	mlt2deg,  &    ! for mlt to mlon
        sinIm_mag(0:nmlat)    ! sinIm

      integer :: jmin, jmax   ! latitude index for interpolation of 
                              ! northward e-field ed2 at mag. equator

!------------------------------------------------------------------------------ 
!  for spherical harmonics
!------------------------------------------------------------------------------ 
      integer, parameter ::   &
     	nm   = 19,     &
     	mm   = 18,     &					
     	nmp  = nm + 1, &					       
     	mmp  = mm + 1	  

      real(r8) :: r(0:nm,0:mm)      ! R_n^m
      real(r8) :: pmopmmo(0:mm)     ! sqrt(1+1/2m)

!------------------------------------------------------------------------------ 
!  index for factors f_m(mlt),f_l(UT),f_-k(d)
!------------------------------------------------------------------------------ 
      integer, parameter :: ni = 1091  ! for n=12 m=-18:18
      integer :: imax                                         ! max number of index
      integer,dimension(0:ni) :: &
     	kf, &
     	lf, &
     	mf, &
     	nf, &
     	jf
      real(r8) :: ft(1:3,0:2)  ! used for f_-k(season,k)

      real(r8) ::  a_klnm(0:ni)        !  A_klm
      real(r8) ::  a_lf(0:ni)          ! A_klmn^lf for minimum &
      real(r8) ::  a_hf(0:ni)          ! A_klmn^hf for maximum

!------------------------------------------------------------------------------ 
! high_latitude boundary
!------------------------------------------------------------------------------ 
      real(r8), parameter ::  & 
     	ef_max  = 0.015_r8,   &  ! max e-field for high latitude boundary location [V/m]
     	lat_sft = 54._r8	 ! shift of highlat_bnd to 54 deg
      integer :: ilat_sft        ! index of shift for high latitude boundary
      integer, parameter :: nmax_sin = 2 ! max. wave number to be represented
      logical, parameter :: debug =.false.

      contains

      subroutine efield_init(efield_lflux_file, efield_hflux_file, efield_wei96_file)
!-----------------------------------------------------------------------
! Purpose: read in and set up coefficients needed for electric field
!          calculation (independent of time & geog. location)
!
! Method:
!
! Author: A. Maute Dec 2003  am 12/17/03 
!-----------------------------------------------------------------------
      character(len=*), intent(in) :: efield_lflux_file
      character(len=*), intent(in) :: efield_hflux_file
      character(len=*), intent(in) :: efield_wei96_file

      call constants	 ! calculate constants
!-----------------------------------------------------------------------
! low/midlatitude potential from Scherliess model
!-----------------------------------------------------------------------
      call read_acoef (efield_lflux_file, efield_hflux_file)	! read in A_klnm for given S_aM
      call index_quiet  ! set up index for f_m(mlt),f_l(UT),f_-k(d)
      call prep_fk	! set up the constant factors for f_k
      call prep_pnm	! set up the constant factors for P_n^m & dP/d phi
!-----------------------------------------------------------------------
! following part should be independent of time & location if IMF constant
!-----------------------------------------------------------------------
      call ReadCoef (efield_wei96_file)

      end subroutine efield_init

      subroutine get_efield
!-----------------------------------------------------------------------
! Purpose: calculates the global electric potential field on the
!          geomagnetic grid (MLT in deg) and derives the electric field 
!
! Method:
!
! Author: A. Maute Dec 2003  am 12/17/03    
!-----------------------------------------------------------------------

      use time_manager,   only : get_curr_calday, get_curr_date
      use mo_solar_parms, only : get_solar_parms
      use mag_parms,      only : get_mag_parms
      use cam_control_mod, only: magfield_fix_year
      use spmd_utils,      only: masterproc

      integer :: idum1, idum2, tod ! time of day [s] 

!-----------------------------------------------------------------------
! get current calendar day of year & date components 
! valid at end of current timestep
!-----------------------------------------------------------------------
      iday = get_curr_calday()                   ! day of year
      call get_curr_date (iyear,imo,iday_m,tod)  ! year, time of day [sec]
      iyear = magfield_fix_year

      if( iyear < 1900 ) then
        write(iulog,"(/,'>>> get_efield: year < 1900 not possible: year=',i5)") iyear
        call endrun
      end if

      ut = tod/3600._r8                   ! UT of day [sec]

!-----------------------------------------------------------------------
! get solar parms
!-----------------------------------------------------------------------
      call get_solar_parms( f107_s = f107d )
!-----------------------------------------------------------------------
! get mag parms
!-----------------------------------------------------------------------
      call get_mag_parms( by = by, bz = bz )
#ifdef EFIELD_DIAGS
      if( masterproc ) then
         write(iulog,*) 'get_efield: f107d,by,bz = ', f107d,by,bz 
      end if
#endif
!-----------------------------------------------------------------------
! ajust S_a
!-----------------------------------------------------------------------
      call adj_S_a
!-----------------------------------------------------------------------
! calculate global electric potential    
!-----------------------------------------------------------------------
      call GlobalElPotential

!-----------------------------------------------------------------------
! calculate derivative of global electric potential 
!-----------------------------------------------------------------------
      call DerivPotential

      end subroutine get_efield

      subroutine GlobalElPotential
!-----------------------------------------------------------------------
! Purpose: calculates the global electric potential field on the
!          geomagnetic grid (MLT in deg) 
!
! Method: rewritten code from Luedger Scherliess (11/20/99 LS)
!     routine to calculate the global electric potential in magnetic
!     Apex coordinates (Latitude and MLT).
!     High Latitude Model is Weimer 1996.
!     Midlatitude model is Scherliess 1999.
!     Interpolation in a transition region at about 60 degree 
!     magnetic apex lat
!
! Author: A. Maute Dec 2003  am 12/17/03 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! local variables
!-----------------------------------------------------------------------
      integer  :: ilon, ilat, idlat
      integer  :: ihlat_bnd(0:nmlon)                  ! high latitude boundary
      integer  :: itrans_width(0:nmlon)               ! width of transition zone
      real(r8) :: mlt, mlon, mlat, mlat_90, pot
      real(r8) :: pot_midlat(0:nmlon,0:nmlat)         ! potential from L. Scherliess model
      real(r8) :: pot_highlat(0:nmlon,0:nmlat)        ! potential from Weimer model
      real(r8) :: pot_highlats(0:nmlon,0:nmlat)	      ! smoothed potential from Weimer model

!-----------------------------------------------------------------------
! Externals
!-----------------------------------------------------------------------
      real(r8),external :: EpotVal        ! in wei96.f

!-----------------------------------------------------------------------
! convert to date and day	
!-----------------------------------------------------------------------
      day  = iday + ut/24._r8
      date = iyear + day/dy2yr

!-----------------------------------------------------------------------
! low/midlatitude electric potential - empirical model Scherliess 1999  
!-----------------------------------------------------------------------
!$omp parallel do private(ilat, ilon, mlat, pot)
      do ilat = 0,nmlath                        ! Calculate only for one magn. hemisphere
	mlat = ylatm(ilat)                      ! mag. latitude
        do ilon = 0,nmlon	                ! lon. loop
          call efield_mid( mlat, ylonm(ilon), pot )
	  pot_midlat(ilon,ilat+nmlath) = pot	! SH/NH symmetry 
	  pot_midlat(ilon,nmlath-ilat) = pot
        end do
      end do

!-----------------------------------------------------------------------
! hight latitude potential from Weimer model
! at the poles Weimer potential is not longitudinal dependent
!-----------------------------------------------------------------------
      call prep_weimer                             ! calculate IMF angle & magnitude, tilt

!$omp parallel do private(ilat, ilon, mlat_90, pot)
      do ilat = 0,nmlat_wei 	                   ! Calculate only for one magn. hemisphere
        mlat_90 = 90._r8 - ylatm(ilat)             ! mag. latitude
        do ilon = 0,nmlon
    	  pot  = 1000._r8*EpotVal( mlat_90, ylonm(ilon)*deg2mlt )     ! calculate potential (kv -> v)
!-----------------------------------------------------------------------
! NH/SH symmetry
!-----------------------------------------------------------------------
    	  pot_highlat(ilon,ilat)        = pot
    	  pot_highlat(ilon,nmlat-ilat)  = pot
    	  pot_highlats(ilon,ilat)       = pot
    	  pot_highlats(ilon,nmlat-ilat) = pot
        end do
      end do     

!-----------------------------------------------------------------------
! weighted smoothing of high latitude potential
!-----------------------------------------------------------------------
      idlat = 2              ! smooth over -2:2 = 5 grid points
      call pot_latsmo( pot_highlats, idlat )
!-----------------------------------------------------------------------
! calculate the height latitude bounday ihl_bnd
! 1. calculate E field from weimar model
!    boundary is set where the total electric field exceeds
!    0.015 V/m (corresp. approx. to 300 m/s)
! 2. moved halfways to 54 deg 
! output : index 0-pole nmlath-equator
!-----------------------------------------------------------------------
      call highlat_getbnd( ihlat_bnd )
!-----------------------------------------------------------------------
! 3. adjust high latitude boundary sinusoidally
!    calculate width of transition zone
!-----------------------------------------------------------------------
      call bnd_sinus( ihlat_bnd, itrans_width ) 
!-----------------------------------------------------------------------
! 4. ajust high latitude potential to low latitude potential      
!-----------------------------------------------------------------------
      call highlat_adjust( pot_highlats, pot_highlat, pot_midlat, ihlat_bnd )
!-----------------------------------------------------------------------
! interpolation of high and low/midlatitude potential in the
! transition zone and put it into global potent array
!-----------------------------------------------------------------------
      call interp_poten( pot_highlats, pot_highlat, pot_midlat, ihlat_bnd, itrans_width) 
!-----------------------------------------------------------------------
! potential weighted smoothing in latitude
!-----------------------------------------------------------------------
      idlat = 2                 ! smooth over -2:2 = 5 grid points
      call pot_latsmo2( potent, idlat )
!-----------------------------------------------------------------------
! potential smoothing in longitude
!-----------------------------------------------------------------------
      idlat = nmlon/48          ! smooth over -idlat:idlat grid points
      call pot_lonsmo( potent, idlat )

!-----------------------------------------------------------------------
! output
!-----------------------------------------------------------------------
! output ( change later to netcdf file)
!      do ilat=0,nmlat
!       do ilon=0,nmlon
!         write(iulog,'(4(x,f12.5))') ylatm(ilat),ylonm(ilon), &
!           potent(ilon,ilat),potent(ilon,nmlat-ilat)
!         write(iulog,'(4(x,f12.5))') ylatm(ilat),ylonm(ilon), &
!           potent(ilon,ilat),potent(ilon,nmlat-ilat)
!	write(iulog,'(f10.3)') potent(ilon,ilat)
!       end do
!      end do

      end subroutine GlobalElPotential

      subroutine ff( ph, mt, f )                                                    
!-----------------------------------------------------------------------
! Purpose: calculate F for normalized associated Legendre polynomial P_n^m
!          Ref.: Richmond J.Atm.Ter.Phys. 1974
!
! Method:  f_m(phi) = sqrt(2) sin(m phi) m > 0
!                   = 1                  m = 0
!                   = sqrt(2) cos(m phi) m < 0
!
! Author: A. Maute Nov 2003  am 11/18/03
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
! dummy arguments
!-----------------------------------------------------------------------
      integer,intent(in)   :: mt
      real(r8),intent(in)  :: ph	! geo. longitude of 0SLT (ut*15)
      real(r8),intent(out) :: f(-mt:mt)

!-----------------------------------------------------------------------
! local variables
!-----------------------------------------------------------------------
      integer  :: m, i, j, mmo
      real(r8) :: sp, cp    

      sp   = sin( ph/rtd )
      cp   = cos( ph/rtd )
      f(0) = 1.e0_r8
                                                                
      f(-1) = sqr2*cp
      f(1)  = sqr2*sp      								 
      do m = 2,mt
        mmo   = m - 1  
        f(m)  = f(-mmo)*sp + cp*f(mmo)
        f(-m) = f(-mmo)*cp - sp*f(mmo)
      end do      

      end subroutine ff                                                                      

      subroutine pnm( ct, p )
!----------------------------------------------------------------------------                                                                   
! Purpose: normalized associated Legendre polynomial P_n^m
!          Ref.: Richmond J.Atm.Ter.Phys. 1974
! Method:
!   P_m^m    = sqrt(1+1/2m)*si*P_m-1^m-1                  m>0
!   P_n^m    = [cos*P_n-1^m - R_n-1^m*P_n-2^m ]/R_n^m     n>m>=0
!   dP/d phi = n*cos*P_n^m/sin-(2*n+1)*R_n^m*P_n-1^m/sin  n>=m>=0
!   R_n^m    = sqrt[ (n^2-m^2)/(4n^2-1) ]
!
! Author: A. Maute Nov 2003  am 11/18/03
!----------------------------------------------------------------------------                                                                   

      implicit none

!-----------------------------------------------------------------------
! dummy arguments
!-----------------------------------------------------------------------
      real(r8), intent(inout) :: ct ! cos(colat)                 
      real(r8), intent(inout) :: p(0:nm,0:mm)

!-----------------------------------------------------------------------
! local variables
!-----------------------------------------------------------------------
      integer  :: mp, m, n, np
      real(r8) :: pm2, st

!      ct = min( ct,.99_r8 )		! cos(colat)
      st = sqrt( 1._r8 - ct*ct ) 	! sin(colat)

      p(0,0) = 1._r8  
      do mp = 1,mmp  ! m+1=1,mm+1
        m = mp - 1
	if( m >= 1 ) then
           p(m,m) = pmopmmo(m)*p(m-1,m-1)*st 			
	end if
	pm2 = 0._r8                                                                  
	do n = mp,nm                    ! n=m+1,N
	   np     = n + 1
	   p(n,m) = (ct*p(n-1,m) - r(n-1,m)*pm2)/r(n,m)
	   pm2    = p(n-1,m)
        end do
      end do

      end subroutine pnm                                                                         

      subroutine prep_pnm
!----------------------------------------------------------------------------                                                                   
! Purpose: constant factors for normalized associated Legendre polynomial P_n^m
!          Ref.: Richmond J.Atm.Ter.Phys. 1974
!
! Method:
!   PmoPmmo(m) = sqrt(1+1/2m)
!   R_n^m      = sqrt[ (n^2-m^2)/(4n^2-1) ]
!
! Author: A. Maute Nov 2003  am 11/18/03
!----------------------------------------------------------------------------                                                                   

      implicit none                

!-----------------------------------------------------------------------
! local variables
!-----------------------------------------------------------------------
      integer  :: mp, m, n
      real(r8) :: xms, xns, den

      do mp = 1, mmp            ! m+1 = 1,mm+1                                     
	m = mp - 1                                               
	xms = m*m                                                
	if( mp /= 1 ) then
           pmopmmo(m) = sqrt( 1._r8 + .5_r8/M )
	end if
	do n = m,nm      ! n = m,N                                     
	  xns    = n*n                                       
	  den    = max(4._r8*xns - 1._r8,1._r8)
	  r(n,m) = sqrt( (xns  - xms)/den )
	end do                 
      end do 

      end subroutine prep_pnm                                                                         

      subroutine index_quiet
!----------------------------------------------------------------------------                                                                   
! Purpose: set up index for factors f_m(mlt),f_l(UT),f_-k(d) to
!    describe the electric potential Phi for the empirical model   
!
! Method:
!    Phi = sum_k sum_l sum_m sum_n [ A_klmn * P_n^m *f_m(mlt)*f_l(UT)*f_-k(d)]
!    - since the electric potential is symmetric about the equator
!      n+m odd terms are set zero resp. not used
!    - in the summation for calculation Phi the index have the following
!      range n=1,12 and m=-n,n, k=0,2 l=-2,2
!
! Author: A. Maute Nov 2003  am 11/18/03
!----------------------------------------------------------------------------                                                                   

      implicit none

!----------------------------------------------------------------------------                                                                   
!	... local variables
!----------------------------------------------------------------------------                                                                   
      integer :: i, j, k, l, n, m

      i = 0 	! initialize
      j = 1 
      do k = 2,0,-1
        do l = -2,2
          if( k == 2 .and. abs(l) == 2 ) then
             cycle
          end if
          do n = 1,12
            do m = -18,18 
              if( abs(m) <= n ) then		    !  |m| < n
                if( (((n-m)/2)*2) == (n-m) ) then   ! only n+m even
             	  if( n-abs(m) <= 9 ) then	    ! n-|m| <= 9 why?
             	    kf(i) = 2-k
             	    lf(i) = l
             	    nf(i) = n
             	    mf(i) = m
             	    jf(i) = j
             	    i	  = i + 1	 ! counter
                  end if
                end if
              end if
            end do ! m
          end do ! n
        end do ! l
      end do ! k

      imax = i - 1  
      if(imax /= ni ) then    ! check if imax == ni 
        write(iulog,'(a19,i5,a18,i5)') 'index_quiet: imax= ',imax, &
          ' not equal to ni =',ni 
        call endrun('index_quiet ERROR')
      end if							
      if(debug) write(iulog,*) 'imax=',imax

      end subroutine index_quiet                                                           

      subroutine read_acoef (efield_lflux_file, efield_hflux_file)
!----------------------------------------------------------------------------                                                                   
! Purpose:  
!    1. read in coefficients A_klmn^lf for solar cycle minimum and
!      A_klmn^hf for maximum 
!    2. adjust S_a (f107d) such that if S_a<80 or S_a > 220 it has reasonable numbers
!      S_aM = [atan{(S_a-65)^2/90^2}-a90]/[a180-a90]
!      a90  = atan [(90-65)/90]^2
!      a180 = atan [(180-65)/90]^2
!    3. inter/extrapolation of the coefficient to the actual flux which is
!      given by the user
!      A_klmn = S_aM [A_klmn^hf-A_klmn^lf]/90. + 2*A_klmn^lf-A_klmn^hf
!
! Method:
!
! Author: A. Maute Nov 2003  am 11/19/03
!----------------------------------------------------------------------------                                                                   

      use ioFileMod,     only : getfil
      use units,         only : getunit, freeunit

      character(len=*), intent(in) :: efield_lflux_file
      character(len=*), intent(in) :: efield_hflux_file

      integer  :: i,ios,unit,istat
      character (len=256):: locfn

!----------------------------------------------------------------------------                                                                   
!  get coefficients file for solar minimum: 
!----------------------------------------------------------------------------                                                                   
      unit     = getunit()
      call getfil( efield_lflux_file, locfn, 0 )

!----------------------------------------------------------------------------                                                                   
! open datafile with coefficients A_klnm
!----------------------------------------------------------------------------                                                                   
      if (debug) write(iulog,*) 'read_acoef: open file ',trim(locfn),' unit ',unit
      open(unit=unit,file=trim(locfn), &
           status = 'old',iostat = ios)
      if(ios.gt.0) then
	write(iulog,*) 'read_acoef: error in opening coeff_lf file',' unit ',unit
        call endrun
      end if

!----------------------------------------------------------------------------                                                                   
! read datafile with coefficients A_klnm
!----------------------------------------------------------------------------                                                                   
      if (debug) write(iulog,*) 'read_acoef: read file ',trim(locfn),' unit ',unit
      read(unit,*,iostat = ios) a_lf
      if(ios.gt.0) then
	write(iulog,*) 'read_acoef: error in reading coeff_lf file',' unit ',unit
        call endrun
      end if

!----------------------------------------------------------------------------                                                                   
! close & free unit      
!----------------------------------------------------------------------------                                                                   
      close(unit)
      call freeunit(unit)
      if (debug) write(iulog,*) 'read_acoef: free unit ',unit

!----------------------------------------------------------------------------                                                                   
!  get coefficients file for solar maximum: 
!----------------------------------------------------------------------------                                                                   
      unit     = getunit()
      call getfil( efield_hflux_file, locfn, 0 )

!----------------------------------------------------------------------------                                                                   
! open datafile with coefficients A_klnm
!----------------------------------------------------------------------------                                                                   
      if (debug) write(iulog,*) 'read_acoef: open file ',trim(locfn),' unit ',unit
      open(unit=unit,file=trim(locfn), &
           status = 'old',iostat = ios)
      if(ios.gt.0) then
	write(iulog,*) 'read_acoef: error in opening coeff_hf file',' unit ',unit
        call endrun
      end if

!----------------------------------------------------------------------------                                                                   
! read datafile with coefficients A_klnm
!----------------------------------------------------------------------------                                                                   
      if (debug) write(iulog,*) 'read_acoef: read file ',trim(locfn)
      read(unit,*,iostat = ios) a_hf
      if(ios.gt.0) then
	write(iulog,*) 'read_acoef: error in reading coeff_hf file',' unit ',unit
        call endrun
      end if

!----------------------------------------------------------------------------                                                                   
! close & free unit      
!----------------------------------------------------------------------------                                                                   
      close(unit)
      call freeunit(unit)
      if (debug) write(iulog,*) 'read_acoef: free unit ',unit

      end subroutine read_acoef

      subroutine adj_S_a
!----------------------------------------------------------------------------                                                                   
! adjust S_a -> S_aM   eqn.8-11 Scherliess draft
!----------------------------------------------------------------------------                                                                   

      implicit none

!----------------------------------------------------------------------------                                                                   
! local variables
!----------------------------------------------------------------------------                                                                   
      integer  :: i
      real(r8) :: x2, y2, a90, a180, S_aM

      x2 = 90._r8*90._r8
      y2 = (90._r8 - 65._r8)
      y2 = y2*y2
      a90  = atan2(y2,x2)
      y2 = (180._r8 - 65._r8)
      y2 = y2*y2
      a180 = atan2(y2,x2)
!     y2 = (S_a-65._r8)
      y2 = (f107d - 65._r8)
      y2 = y2*y2
      S_aM = (atan2(y2,x2) - a90)/(a180 - a90) 
      S_aM = 90._r8*(1._r8 + S_aM)
      if(debug) write(iulog,*) 'f107d=',f107d,' S_aM =',S_aM
      if(debug) write(iulog,*) 'By=',by

!----------------------------------------------------------------------------                                                                   
! inter/extrapolate to S_a (f107d)
!----------------------------------------------------------------------------                                                                   
      do i = 0,ni                       ! eqn.8 Scherliess draft
        a_klnm(i) = S_aM*(a_hf(i) - a_lf(i))/90._r8 + 2._r8*a_lf(i) - a_hf(i)
! for testing like in original code
!        a_klnm(i)=S_a*(a_hf(i)-a_lf(i))/90.+2.*a_lf(i)-a_hf(i)
!        a_klnm(i)=f107d*(a_hf(i)-a_lf(i))/90.+2.*a_lf(i)-a_hf(i)
      end do

      end subroutine adj_S_a

      subroutine constants
!----------------------------------------------------------------------------                                                                   
! Purpose: set up constant values (e.g. magnetic grid, convertion
!      constants etc)
!
! Method:
!
! Author: A. Maute Nov 2003  am 11/19/03
!----------------------------------------------------------------------------                                                                   

!----------------------------------------------------------------------------                                                                   
! local variables
!----------------------------------------------------------------------------                                                                   
      integer  :: i,j
      real(r8) :: fac,lat

      rtd     = 180._r8/pi	        ! radians -> deg
      dtr     = pi/180._r8	        ! deg -> radians
      sqr2    = sqrt(2.e0_r8)
      hr2rd   = pi/12._r8	        ! pi/12 hrs
      dy2rd   = 2._r8*pi/dy2yr          ! 2*pi/365.24  average year
      deg2mlt = 24._r8/360._r8          ! convert degrees to MLT hours
      mlt2deg = 360._r8/24._r8          ! for mlt to mlon       

!----------------------------------------------------------------------------                                                                   
! Set grid deltas:
!----------------------------------------------------------------------------                                                                   
      dlatm = 180._r8/nmlat
      dlonm = 360._r8/nmlon

!----------------------------------------------------------------------------                                                                   
! Set magnetic latitude array 
!----------------------------------------------------------------------------                                                                   
      do j = 0,nmlat
        ylatm(j) = j*dlatm
        lat = (ylatm(j) - 90._r8)*dtr
	fac = cos(lat)                   ! sinIm = 2*sin(lam_m)/sqrt[4-3*cos^2(lam_m)]
	fac = 4._r8 - 3._r8*fac*fac
	fac = 2._r8/sqrt( fac )
	sinIm_mag(j) = fac*sin( lat )
      end do 

!----------------------------------------------------------------------------                                                                   
! Set magnetic longitude array
!----------------------------------------------------------------------------                                                                   
      do i = 0,nmlon
        ylonm(i) = i*dlonm
      end do ! i=1,nmlonp1

!----------------------------------------------------------------------------                                                                   
! find boundary index for weimer
!----------------------------------------------------------------------------                                                                   
      do j = 0,nmlat
        nmlat_wei = j
        if( bnd_wei <= ylatm(j) ) then
           exit
        end if
      end do 

!----------------------------------------------------------------------------                                                                   
! find latitudinal shift
!----------------------------------------------------------------------------                                                                   
      do j = 0,nmlat
        ilat_sft = j
        if( lat_sft <= ylatm(j) ) then
           exit
        end if
      end do 

!----------------------------------------------------------------------------                                                                   
! find index for linear interpolation of ed2 at mag.equator 
! use 12 deg - same as in TIEGCM      
!----------------------------------------------------------------------------                                                                   
      do j = 0,nmlat
        lat = ylatm(j) - 90._r8
        if( lat <= -12._r8 ) then
	  jmin = j
        else if( lat > 12._r8 ) then
	  jmax = j
	  exit
       end if
      end do

      end subroutine constants

      subroutine prep_fk
!----------------------------------------------------------------------------                                                                   
! Purpose: set up constants factors for f_-k(day) used for empirical model
!     to calculate the electric potential
!
! Method:
!
! Author: A. Maute Nov 2003  am 11/19/03
!----------------------------------------------------------------------------                                                                   

      ft(1,0) = .75_r8*sqrt( 6.e0_r8 )/pi			
      ft(1,1) = 2.e0_r8*ft(1,0)					      
      ft(1,2) = 1.e0_r8						      
      ft(2,0) = ft(1,0) 					      
      ft(2,1) = -ft(1,1)					      
      ft(2,2) = 1.e0_r8						      
      ft(3,0) = ft(2,1) 					      
      ft(3,1) = 0._r8						      
      ft(3,2) = 1.e0_r8							   

      end subroutine prep_fk

      subroutine set_fkflfs( fk, fl, fs )
!----------------------------------------------------------------------------                                                                   
! Purpose:  set f_-k(day) depending on seasonal flag used for empirical model
!     to calculate the electric potential
!
! Method:
!
! Author: A. Maute Nov 2003  am 11/20/03
!----------------------------------------------------------------------------                                                                   

!----------------------------------------------------------------------------                                                                   
!	... dummy arguments
!----------------------------------------------------------------------------                                                                   
      real(r8), intent(out) ::  &
     	fk(0:2),  &	                ! f_-k(day) 
     	fl(-2:2), &	                ! f_l(ut)  
     	fs(2)		                ! f_s(f10.7) 
!----------------------------------------------------------------------------                                                                   
! local variables
!----------------------------------------------------------------------------                                                                   
      integer  :: lp
      real(r8) :: ang
      real(r8) :: lon_ut

!----------------------------------------------------------------------------                                                                   
! f_-k(day) 
! use factors for iseasav == 0 - Scherliess had iseasav as an input parameter
!----------------------------------------------------------------------------                                                                   
      lp = iseasav
      if( iseasav == 0 ) then
        ang   = (day + 9._r8)*dy2rd
        fk(0) = sqr2*cos( 2._r8*ang )
        fk(1) = sqr2*cos( ang )
        fk(2) = 1._r8
      else if( iseasav >= 1 .and. iseasav <= 3 ) then
        fk(0) = ft(lp,0)
        fk(1) = ft(lp,1)
        fk(2) = ft(lp,2)
      else if( iseasav == 4 ) then
        fk(0) =0._r8
        fk(1) =0._r8
        fk(2) =1._r8
      end if

!----------------------------------------------------------------------------                                                                   
! f_l(ut) 
!----------------------------------------------------------------------------                                                                   
      lon_ut = 15._r8*ut        ! 15.*mlt - xmlon + 69. 
      call ff( lon_ut, 2, fl )                                                 
      if( iutav ) then  	! UT-averaging
     
	ang   = fl(0)
        fl(:) = 0._r8
        fl(0) = ang
	
      end if

!----------------------------------------------------------------------------                                                                   
! f_s(f10.7)  only fs(1) used  	
!----------------------------------------------------------------------------                                                                   
      fs(1) = 1._r8
!     fs(2) = S_a			  
      fs(2) = f107d			  

      end subroutine set_fkflfs

      subroutine efield_mid( mlat, mlon, pot )
!----------------------------------------------------------------------------                                                                   
! Purpose: calculate the electric potential for low and 
!      midlatitudes from an empirical model (Scherliess 1999)
!
! Method:
!
! Author: A. Maute Nov 2003  am 11/20/03
!----------------------------------------------------------------------------                                                                   

!----------------------------------------------------------------------------                                                                   
!	... dummy arguments
!----------------------------------------------------------------------------                                                                   
      real(r8), intent(in)  :: mlat, mlon
      real(r8), intent(out) :: pot               ! electric potential (V)

!----------------------------------------------------------------------------                                                                   
! local variables
!----------------------------------------------------------------------------                                                                   
      integer  :: i, mp, np, nn
      real(r8) :: mod_mlat, ct, x
      real(r8) :: fk(0:2)      	    ! f_-k(day) 
      real(r8) :: fl(-2:2)          ! f_l(ut)  
      real(r8) :: fs(2)	            ! f_s(f10.7) 
      real(r8) :: f(-18:18)
      real(r8) :: p(0:nm,0:mm)      ! P_n^m	 spherical harmonics

      pot = 0._r8 ! initialize                                        

      mod_mlat = mlat
      if( abs(mlat) <= 0.5_r8 ) then
         mod_mlat = 0.5_r8                     ! avoid geomag.equator
      end if

!----------------------------------------------------------------------------                                                                   
! set f_-k, f_l, f_s depending on seasonal flag
!----------------------------------------------------------------------------                                                                   
      call set_fkflfs( fk, fl, fs ) 
      
!----------------------------------------------------------------------------                                                                   
! spherical harmonics 
!----------------------------------------------------------------------------                                                                   
      ct = cos( (90._r8 - mod_mlat)*dtr )  ! magnetic colatitude 
      call pnm( ct, p )	                   ! calculate P_n^m
      call ff( mlon, 18, f )               ! calculate f_m (phi) why 18 if N=12                              

      do i = 0,imax  
        mp  = mf(i)                                                      
        np  = nf(i)
        nn  = abs(mp)                      !   P_n^m = P_n^-m  
        x   = a_klnm(i)* fl(lf(i)) * fk(kf(i)) * fs(jf(i))
	pot = pot + x*f(mp)*p(np,nn) 
      end do 
      
      end subroutine efield_mid                                              

      subroutine prep_weimer
!----------------------------------------------------------------------------                                                                   
! Purpose:  for Weimer model calculate IMF angle, IMF magnitude
!  tilt of earth
!
! Method: using functions and subroutines from Weimer Model 1996
!     output:  angle, &  ! IMF angle
!     	       bt,    &  ! IMF magnitude
!     	       tilt      ! tilt of earth
!
! Author: A. Maute Nov 2003  am 11/20/03
!----------------------------------------------------------------------------                                                                   

!----------------------------------------------------------------------------                                                                   
!  local variables
!----------------------------------------------------------------------------                                                                   
      real(r8) ::  &
        angle,  &  ! IMF angle
	bt,     &  ! IMF magnitude
	tilt       ! tilt of earth

!----------------------------------------------------------------------------                                                                   
! function declarations
!----------------------------------------------------------------------------                                                                   
      real(r8), external :: get_tilt	 ! in wei96.f

      if( by == 0._r8 .and. bz == 0._r8) then
         angle = 0._r8
      else
         angle = atan2( by,bz )
      end if
      
      angle = angle*rtd
      call adjust( angle )
      bt = sqrt( by*by + bz*bz )
!----------------------------------------------------------------------------                                                                   
! use month and day of month - calculated with average no.of days per month
! as in Weimer
!----------------------------------------------------------------------------                                                                   
      if(debug) write(iulog,*) 'prep_weimer: day->day of month',iday,imo,iday_m,ut
      tilt = get_tilt( iyear, imo, iday_m, ut )

       if(debug) then
        write(iulog,"(/,'efield prep_weimer:')")
        write(iulog,*)  '  Bz   =',bz
        write(iulog,*)  '  By   =',by
        write(iulog,*)  '  Bt   =',bt
        write(iulog,*)  '  angle=',angle
        write(iulog,*)  '  VSW  =',v_sw
        write(iulog,*)  '  tilt =',tilt
       end if

      call SetModel( angle, bt, tilt, v_sw )

      end subroutine prep_weimer

      subroutine pot_latsmo( pot, idlat )  ! pots == pot_highlats
!----------------------------------------------------------------------------                                                                   
! Purpose: smoothing in latitude of  potential
!
! Method: weighted smoothing in latitude 
! assume regular grid spacing
!
! Author: A. Maute Nov 2003  am 11/20/03
!----------------------------------------------------------------------------                                                                   

!----------------------------------------------------------------------------                                                                   
!	... dummy arguments
!----------------------------------------------------------------------------                                                                   
      integer, intent(in)     :: idlat
      real(r8), intent(inout) :: pot(0:nmlon,0:nmlat)

!----------------------------------------------------------------------------                                                                   
! local variables
!----------------------------------------------------------------------------                                                                   
      integer  :: ilon, ilat, id
      real(r8) :: wgt, del
      real(r8) :: w(-idlat:idlat)
!     real(r8) :: pot_smo(0:nmlat) ! temp array for smooth. potential
      real(r8) :: pot_smo(0:nmlon,0:nmlat_wei) ! temp array for smooth. potential

!----------------------------------------------------------------------------                                                                   
! weighting factors (regular grid spacing) 
!----------------------------------------------------------------------------                                                                   
      wgt = 0._r8 
      do id = -idlat,idlat
	del   = abs(id)*dlatm	! delta lat_m
	w(id) = 1._r8/(del + 1._r8)
        wgt   = wgt + w(id)
      end do
      wgt = 1._r8/wgt

!     do ilon = 0,nmlon
!        do ilat = idlat,nmlat_wei-idlat
!       do ilat = idlat,nmlat-idlat
!         pot_smo(ilat) = 0._r8
!         do id = -idlat,idlat	!  org. was degree now grid points
!           pot_smo(ilat) = pot_smo(ilat) + w(id)*pot(ilon,ilat+id)
!           write(iulog,"('pot_latsmo: ilon=',i3,' ilat=',i3,' id=',i3,' pot(ilon,ilat+id)=',e12.4)") ilon,ilat,id,pot(ilon,ilat+id)
!         end do
!         pot_smo(ilat)       = pot_smo(ilat)*wgt
!         pot_smo(nmlat-ilat) = pot_smo(ilat)
!       end do
!      pot(ilon,idlat:nmlat-idlat) =  &        ! copy back into pot
!         pot_smo(idlat:nmlat-idlat)
!        pot(ilon,idlat:nmlat_wei-idlat)       = pot_smo(idlat:nmlat_wei-idlat)
!       pot(ilon,nmlat-nmlat_wei+idlat:nmlat) = pot_smo(nmlat-nmlat_wei+idlat:nmlat)
!        pot(ilon,nmlat-nmlat_wei+idlat:nmlat-idlat) = pot_smo(nmlat-nmlat_wei+idlat:nmlat-idlat)
!     end do

!$omp parallel do private(ilat)
      do ilat = idlat,nmlat_wei-idlat
         pot_smo(:,ilat) = matmul( pot(:,ilat-idlat:ilat+idlat),w )*wgt
      end do

      do ilat = idlat,nmlat_wei-idlat
         pot(:,ilat)       = pot_smo(:,ilat)
         pot(:,nmlat-ilat) = pot_smo(:,ilat)
      end do

      end subroutine pot_latsmo

      subroutine pot_latsmo2( pot, idlat ) 
!----------------------------------------------------------------------------                                                                   
! Purpose:  smoothing in latitude of  potential
!
! Method: weighted smoothing in latitude 
!         assume regular grid spacing
!
! Author: A. Maute Nov 2003  am 11/20/03
!----------------------------------------------------------------------------                                                                   

!----------------------------------------------------------------------------                                                                   
!	... dummy arguments
!----------------------------------------------------------------------------                                                                   
      integer, intent(in)     :: idlat
      real(r8), intent(inout) :: pot(0:nmlon,0:nmlat)

!----------------------------------------------------------------------------                                                                   
! local variables
!----------------------------------------------------------------------------                                                                   
      integer  :: ilon, ilat, id
      real(r8) :: wgt, del
      real(r8) :: w(-idlat:idlat)
!     real(r8) :: pot_smo(0:nmlat) ! temp array for smooth. potential
      real(r8) :: pot_smo(0:nmlon,0:nmlath) ! temp array for smooth. potential

!----------------------------------------------------------------------------                                                                   
! weighting factors (regular grid spacing)  
!----------------------------------------------------------------------------                                                                   
      wgt = 0._r8
      do id = -idlat,idlat
	del   = abs(id)*dlatm	! delta lat_m
	w(id) = 1._r8/(del + 1._r8)
        wgt   = wgt + w(id)
      end do
      wgt = 1._r8/wgt

!     do ilon = 0,nmlon
!       do ilat = idlat,nmlath-idlat  ! ilat = 5:175
!         pot_smo(ilat) = 0._r8
!         do id = -idlat,idlat	!  org. was degree now grid points
!           pot_smo(ilat) = pot_smo(ilat) + w(id)*pot(ilon,ilat+id)
!         end do
!         pot_smo(ilat) = pot_smo(ilat)*wgt
!       end do
!       pot(ilon,idlat:nmlath-idlat) = pot_smo(idlat:nmlath-idlat) ! copy back into pot
!     end do

!$omp parallel do private(ilat)
      do ilat = idlat,nmlath-idlat
         pot_smo(:,ilat) = matmul( pot(:,ilat-idlat:ilat+idlat),w )*wgt
      end do

      do ilat = idlat,nmlath-idlat
         pot(:,ilat) = pot_smo(:,ilat)
      end do

      end subroutine pot_latsmo2

      subroutine pot_lonsmo( pot, idlon ) 
!----------------------------------------------------------------------------                                                                   
! Purpose: smoothing in longitude of potential
!
! Method:  weighted smoothing in longitude
!          assume regular grid spacing
!
! Author: A. Maute Nov 2003  am 11/20/03
!----------------------------------------------------------------------------                                                                   

!----------------------------------------------------------------------------                                                                   
!	... dummy arguments
!----------------------------------------------------------------------------                                                                   
      integer, intent(in)     :: idlon
      real(r8), intent(inout) :: pot(0:nmlon,0:nmlat)

!----------------------------------------------------------------------------                                                                   
! local variables
!----------------------------------------------------------------------------                                                                   
      integer  :: ilon, ilat, id, iabs
      real(r8) :: wgt, del
      real(r8) :: w(-idlon:idlon)
      real(r8) :: pot_smo(0:nmlath) ! temp array for smooth. potential
      real(r8) :: tmp(-idlon:nmlon+idlon) ! temp array for smooth. potential

!----------------------------------------------------------------------------                                                                   
! weighting factors (regular grid spacing) 
!----------------------------------------------------------------------------                                                                   
      wgt = 0._r8
      do id = -idlon,idlon
	del   = abs(id)*dlonm	! delta lon_m
	w(id) = 1._r8/(del + 1._r8)
        wgt   = wgt + w(id)
      end do
     wgt = 1._r8/wgt

!----------------------------------------------------------------------------                                                                   
! averaging     
!----------------------------------------------------------------------------                                                                   
!     do ilon = 0,nmlon
!       do ilat = 0,nmlath
!         pot_smo(ilat) = 0._r8
!         do id = -idlon,idlon	                  !  org. was degree now grid points
!           iabs = ilon + id
!           if( iabs > nmlon ) then
!              iabs = iabs - nmlon ! test if wrap around
!           end if
!           if( iabs < 0 ) then
!              iabs = iabs + nmlon ! test if wrap around
!           end if
!           pot_smo(ilat) = pot_smo(ilat) + w(id)*pot(iabs,ilat)
!         end do
!         pot_smo(ilat)  = pot_smo(ilat)*wgt
!         pot(ilon,ilat) = pot_smo(ilat)       ! copy back into pot 
!         pot(ilon,nmlat-ilat) = pot_smo(ilat) ! copy back into pot    
!       end do   
!     end do

!$omp parallel do private(ilat,ilon,tmp)
      do ilat = 0,nmlath
          tmp(0:nmlon)             = pot(0:nmlon,ilat)
          tmp(-idlon:-1)           = pot(nmlon-idlon:nmlon-1,ilat)
          tmp(nmlon+1:nmlon+idlon) = pot(1:idlon,ilat)
          do ilon = 0,nmlon
	     pot(ilon,ilat) = dot_product( tmp(ilon-idlon:ilon+idlon),w )*wgt
	     pot(ilon,nmlat-ilat) = pot(ilon,ilat)
          end do   
      end do
      
      end subroutine pot_lonsmo

      subroutine highlat_getbnd( ihlat_bnd ) 
!----------------------------------------------------------------------------                                                                   
! Purpose: calculate the height latitude bounday index ihl_bnd
!
! Method:
! 1. calculate E field from weimar model
!    boundary is set where the total electric field exceeds
!    0.015 V/m (corresp. approx. to 300 m/s)
! 2. moved halfways to 54 deg not necessarily equatorwards as in the
!    original comment from L. Scherliess- or?
!
! Author: A. Maute Nov 2003  am 11/20/03
!----------------------------------------------------------------------------                                                                   

!----------------------------------------------------------------------------                                                                   
!	... dummy arguments
!----------------------------------------------------------------------------                                                                   
      integer, intent(out) :: ihlat_bnd(0:nmlon)

!----------------------------------------------------------------------------                                                                   
! local variables
!----------------------------------------------------------------------------                                                                   
      integer  :: ilon, ilat, ilat_sft_rvs
      real(r8) :: mlat, mlt, es, ez, e_tot

      ilat_sft_rvs = nmlath - ilat_sft          ! pole =0, equ=90
!$omp parallel do private(ilat,ilon,mlt,mlat,es,ez,e_tot)
      do ilon = 0,nmlon                         ! long.
	ihlat_bnd(ilon) = 0
        mlt  = ylonm(ilon)*deg2mlt              ! mag.local time ?
        do ilat = nmlat_wei+1,0,-1              ! lat. loop moving torwards pole
	  mlat = 90._r8 - ylatm(ilat)           ! mag. latitude pole = 90 equator = 0
          call gecmp( mlat, mlt, es, ez )	! get electric field
          e_tot = sqrt( es**2 + ez**2 )
          if( abs(e_tot) >= ef_max ) then                        ! e-filed > limit -> boundary
            ihlat_bnd(ilon) = ilat - (ilat - ilat_sft_rvs)/2     ! shift boundary to lat_sft (54deg)
            exit
          end if
        end do
      end do     

!     write(iulog,"('highlat_getbnd: ihlat_bnd=',/,(12i6))") ihlat_bnd

      end subroutine highlat_getbnd

      subroutine bnd_sinus( ihlat_bnd, itrans_width )  
!----------------------------------------------------------------------------                                                                   
! Purpose: 
!   1. adjust high latitude boundary (ihlat_bnd) sinusoidally
!   2. width of transition zone from midlatitude potential to high latitude
!      potential (itrans_width)
!
! Method:
! 1.adjust boundary sinusoidally
!   max. wave number to be represented nmax_sin
!   RHS(mi) = Sum_phi Sum_(mi=-nmax_sin)^_(mi=nmax_sin) f_mi(phi)*hlat_bnd(phi) 
!   U(mi,mk)   = Sum_phi Sum_(mi=-nmax_sin)^_(mi=nmax_sin) f_mi(phi) *
!                Sum_(mk=-nmax_sin)^_(mk=nmax_sin) f_mk(phi)
!   single values decomposition of U
!   solving U*LSG = RHS 
!   calculating hlat_bnd:
!   hlat_bnd = Sum_(mi=-nmax_sin)^_(mi=nmax_sin) f_mi(phi)*LSG(mi)
!
! 2. width of transition zone from midlatitude potential to high latitude
!    potential
!    trans_width(phi)=8.-2.*cos(phi) 
!
! Author: A. Maute Nov 2003  am 11/20/03
!----------------------------------------------------------------------------                                                                   

      use sv_decomp, only : svdcmp, svbksb
     
!----------------------------------------------------------------------------                                                                   
!	... dummy arguments
!----------------------------------------------------------------------------                                                                   
      integer, intent(inout) :: ihlat_bnd(0:nmlon)    ! loaction of boundary
      integer, intent(out)   :: itrans_width(0:nmlon) ! width of transition zone

!----------------------------------------------------------------------------                                                                   
! local variables
!----------------------------------------------------------------------------                                                                   
      integer, parameter :: nmax_a = 2*nmax_sin+1 ! absolute array length
      integer, parameter :: ishf   = nmax_sin+1
      integer  :: ilon, i, i1, j, bnd
      real(r8) :: sum, mlon
      real(r8) :: rhs(nmax_a)
      real(r8) :: lsg(nmax_a)
      real(r8) :: u(nmax_a,nmax_a)
      real(r8) :: v(nmax_a,nmax_a)
      real(r8) :: w(nmax_a,nmax_a)
      real(r8) :: f(-nmax_sin:nmax_sin,0:nmlon)

!----------------------------------------------------------------------------                                                                   
!    Sinusoidal Boundary calculation
!----------------------------------------------------------------------------                                                                   
      rhs(:) = 0._r8
      lsg(:) = 0._r8
      u(:,:) = 0._r8
      v(:,:) = 0._r8
      w(:,:) = 0._r8

      do ilon = 0,nmlon                  ! long.
        bnd  = nmlath - ihlat_bnd(ilon)  ! switch from pole=0 to pole =90
        call ff( ylonm(ilon), nmax_sin, f(-nmax_sin,ilon) )
        do i = -nmax_sin,nmax_sin
	  i1 = i + ishf
          rhs(i1) = rhs(i1) + f(i,ilon) * bnd
!	  write(iulog,*) 'rhs ',ilon,i1,bnd,f(i,ilon),rhs(i+ishf)
          do j = -nmax_sin,nmax_sin 
            u(i1,j+ishf) = u(i1,j+ishf) + f(i,ilon)*f(j,ilon)
!	    write(iulog,*) 'u ',ilon,i1,j+ishf,u(i+ishf,j+ishf)
          end do
        end do
      end do

!     if (debug) write(iulog,*) ' Single Value Decomposition'
      call svdcmp( u, nmax_a, nmax_a, nmax_a, nmax_a, w, v )

!     if (debug) write(iulog,*) ' Solving'
      call svbksb( u, w, v, nmax_a, nmax_a, nmax_a, nmax_a, rhs, lsg )
!      
      do ilon = 0,nmlon  ! long.
!       sum = 0._r8
	sum = dot_product( lsg(-nmax_sin+ishf:nmax_sin+ishf),f(-nmax_sin:nmax_sin,ilon) )
!       do i = -nmax_sin,nmax_sin
!         sum = sum + lsg(i+ishf)*f(i,ilon)  
!       end do
        ihlat_bnd(ilon)    = nmlath - int( sum + .5_r8 )                                ! closest point
        itrans_width(ilon) = int( 8._r8 - 2._r8*cos( ylonm(ilon)*dtr ) + .5_r8 )/dlatm  ! 6 to 10 deg.
      end do      
!     write(iulog,"('bnd_sinus: ihlat_bnd=',/,(12i6))") ihlat_bnd
!     write(iulog,"('bnd_sinus: itrans_width=',/,(12i6))") itrans_width

      end subroutine bnd_sinus

      subroutine highlat_adjust( pot_highlats, pot_highlat, pot_midlat, ihlat_bnd )
!----------------------------------------------------------------------------                                                                   
! Purpose: Adjust mid/low latitude electric potential and high latitude
!          potential such that there are continous across the mid to high 
!          latitude boundary
!
! Method:
! 1. integrate Phi_low/mid(phi,bnd) along the boundary mid to high latitude
! 2. integrate Phi_high(phi,bnd) along the boundary mid to high latitude
! 3. adjust Phi_high by delta =
!    Int_phi Phi_high(phi,bnd) d phi/360. - Int_phi Phi_low/mid(phi,bnd) d phi/360.
!
! Author: A. Maute Nov 2003  am 11/21/03
!----------------------------------------------------------------------------                                                                   

!----------------------------------------------------------------------------                                                                   
!	... dummy arguments
!----------------------------------------------------------------------------                                                                   
      integer, intent(in)     :: ihlat_bnd(0:nmlon)	                  ! boundary mid to high latitude
      real(r8), intent(in)    :: pot_midlat(0:nmlon,0:nmlat)              ! low/mid latitude potentail
      real(r8), intent(inout) :: pot_highlat(0:nmlon,0:nmlat)             ! high_lat potential
      real(r8), intent(inout) :: pot_highlats(0:nmlon,0:nmlat)            ! high_lat potential! smoothed high_lat potential

!----------------------------------------------------------------------------                                                                   
! local:     
!----------------------------------------------------------------------------                                                                   
      integer  :: bnd, ilon, ilat, ilatS, ibnd60, ibnd_hl
      real(r8) :: pot60, pot_hl, del

!----------------------------------------------------------------------------                                                                   
! 1. integrate Phi_low/mid(phi,bnd) along the boundary mid to high latitude
! 2. integrate Phi_high(phi,bnd) along the boundary mid to high latitude
!----------------------------------------------------------------------------                                                                   
      pot60  = 0._r8
      pot_hl = 0._r8
      do ilon = 1,nmlon  ! long.           ! bnd -> eq to pole 0:90
    	ibnd60  = nmlat - ihlat_bnd(ilon)   ! 0:180 pole to pole
    	ibnd_hl = ihlat_bnd(ilon)         ! colatitude
        pot60   = pot60 + pot_midlat(ilon,ibnd60)
        pot_hl  = pot_hl + pot_highlats(ilon,ibnd_hl)
      end do
      pot60  = pot60/(nmlon)
      pot_hl = pot_hl/(nmlon)
      
      if (debug) write(iulog,*) 'Mid-Latitude Boundary Potential =',pot60
      if (debug) write(iulog,*) 'High-Latitude Boundary Potential=',pot_hl

!----------------------------------------------------------------------------                                                                   
! 3. adjust Phi_high by delta =
!    Int_phi Phi_high(phi,bnd) d phi/360. - Int_phi Phi_low/mid(phi,bnd) d phi/360.
!----------------------------------------------------------------------------                                                                   
      del = pot_hl - pot60

!$omp parallel do private(ilat,ilon,ilats)
      do ilat = 0,nmlat_wei      ! colatitude
        ilats = nmlat - ilat
        do ilon = 0,nmlon
	  pot_highlat(ilon,ilat)   = pot_highlat(ilon,ilat)   - del
	  pot_highlat(ilon,ilats)  = pot_highlat(ilon,ilats)  - del
	  pot_highlats(ilon,ilat)  = pot_highlats(ilon,ilat)  - del
	  pot_highlats(ilon,ilats) = pot_highlats(ilon,ilats) - del
        end do
      end do

      end subroutine highlat_adjust

      subroutine interp_poten( pot_highlats, pot_highlat, pot_midlat, &
                               ihlat_bnd, itrans_width ) 
!----------------------------------------------------------------------------                                                                   
! Purpose: construct a smooth global electric potential field 
!
! Method: construct one global potential field
! 1. low/mid latitude: |lam| < bnd-trans_width
!   Phi(phi,lam) = Phi_low(phi,lam)
! 2. high latitude: |lam| > bnd+trans_width
!   Phi(phi,lam) = Phi_hl(phi,lam)
! 3. transition zone: bnd-trans_width <= lam <= bnd+trans_width 
! a. interpolate between high and low/midlatitude potential
!   Phi*(phi,lam) = 1/15*[ 5/(2*trans_width) * {Phi_low(phi,bnd-trans_width)*
!   [-lam+bnd+trans_width] + Phi_hl(phi,bnd+trans_width)*
!   [lam-bnd+trans_width]} + 10/(2*trans_width) {Phi_low(phi,lam)*
!   [-lam+bnd+trans_width] + Phi_hl(phi,lam)*
!   [lam-bnd+trans_width]}]
! b.  Interpolate between just calculated Potential and the high latitude
!    potential in a 3 degree zone poleward of the boundary:
!    bnd+trans_width < lam <= bnd+trans_width+ 3 deg 
!   Phi(phi,lam) = 1/3 { [3-(lam-bnd-trans_width)]* Phi*(phi,lam) +
!   [lam-bnd-trans_width)]* Phi_hl*(phi,lam) }
!
! Author: A. Maute Nov 2003  am 11/21/03      
!----------------------------------------------------------------------------                                                                   

!----------------------------------------------------------------------------                                                                   
!	... dummy arguments
!----------------------------------------------------------------------------                                                                   
      integer, intent(in)  :: ihlat_bnd(0:nmlon)
      integer, intent(in)  :: itrans_width(0:nmlon)
      real(r8), intent(in) :: pot_highlats(0:nmlon,0:nmlat)
      real(r8), intent(in) :: pot_highlat(0:nmlon,0:nmlat)
      real(r8), intent(in) :: pot_midlat(0:nmlon,0:nmlat)

!----------------------------------------------------------------------------                                                                   
! local variables
!----------------------------------------------------------------------------                                                                   
      real(r8), parameter :: fac = 1._r8/3._r8
      integer  :: ilon, ilat
      integer  :: ibnd, tw, hb1, hb2, lat_ind
      integer  :: j1, j2
      real(r8) :: a, b, lat, b1, b2
      real(r8) :: wrk1, wrk2

!$omp parallel do private(ilat,ilon,ibnd,tw)
      do ilon = 0,nmlon
        ibnd = ihlat_bnd(ilon)     ! high latitude boundary index
	tw   = itrans_width(ilon)  ! width of transition zone (index)
!----------------------------------------------------------------------------                                                                   
! 1. low/mid latitude: |lam| < bnd-trans_width
!   Phi(phi,lam) = Phi_low(phi,lam)
!----------------------------------------------------------------------------                                                                   
        do ilat = 0,nmlath-(ibnd+tw+1)
          potent(ilon,nmlath+ilat) = pot_midlat(ilon,nmlath+ilat)
          potent(ilon,nmlath-ilat) = pot_midlat(ilon,nmlath+ilat)
        end do
!----------------------------------------------------------------------------                                                                   
! 2. high latitude: |lam| > bnd+trans_width
!   Phi(phi,lam) = Phi_hl(phi,lam)
!----------------------------------------------------------------------------                                                                   
        do ilat = 0,ibnd-tw-1
          potent(ilon,ilat)       = pot_highlats(ilon,nmlat-ilat)
          potent(ilon,nmlat-ilat) = pot_highlats(ilon,nmlat-ilat)
        end do
      end do
!----------------------------------------------------------------------------                                                                   
! 3. transition zone: bnd-trans_width <= lam <= bnd+trans_width 
!----------------------------------------------------------------------------                                                                   
! a. interpolate between high and low/midlatitude potential
! update only southern hemisphere (northern hemisphere is copied
! after smoothing)
!----------------------------------------------------------------------------                                                                   
!$omp parallel do private(ilat,ilon,ibnd,tw,a,b,b1,b2,hb1,hb2,lat_ind,j1,j2,wrk1,wrk2)
      do ilon = 0,nmlon
        ibnd = ihlat_bnd(ilon)          ! high latitude boundary index
	tw   = itrans_width(ilon)       ! width of transition zone (index)
        a    = 1._r8/(2._r8*tw)
	b1   = (nmlath - ibnd + tw)*a
	b2   = (nmlath - ibnd - tw)*a
	hb1  = nmlath - (ibnd + tw)
	j1   = nmlath - hb1
	hb2  = nmlath - (ibnd - tw)
	j2   = nmlath - hb2
	wrk1 = pot_midlat(ilon,j1)
	wrk2 = pot_highlats(ilon,j2)
!        write(iulog,*) 'pot_all ',ilon,hb1,hb2,nmlath -ibnd,tw
	do ilat = ibnd-tw,ibnd+tw
	  lat_ind = nmlath - ilat
          potent(ilon,ilat) = &
     	     fac*((wrk1 + 2._r8*pot_midlat(ilon,ilat))*(b1 - a*lat_ind) &
     	          + (wrk2 + 2._r8*pot_highlats(ilon,ilat))*(a*lat_ind - b2))
          potent(ilon,nmlat-ilat) = potent(ilon,ilat)
        end do
!----------------------------------------------------------------------------                                                                   
! b.  Interpolate between just calculated Potential and the high latitude
!    potential in a 3 degree zone poleward of the boundary
!----------------------------------------------------------------------------                                                                   
	do ilat = hb2+1,nmlath
	  a = max( 3._r8/dlatm - (ilat - hb2 - 1),0._r8 )
	  b = 3._r8/dlatm - a
          potent(ilon,nmlath-ilat) = (a*potent(ilon,nmlath-ilat)   &
                                      + b*pot_highlat(ilon,nmlath-ilat))/3._r8*dlatm
          potent(ilon,nmlath+ilat) = potent(ilon,nmlath-ilat)
        end do
      end do      

      end subroutine interp_poten

      subroutine DerivPotential
!-----------------------------------------------------------------------
! Purpose: calulates the electric field [V/m] from the electric potential
!
! Method:  Richmond [1995] eqn 5.9-5.10
! ed1(:,:) = Ed1 = - 1/[R cos lam_m] d PHI/d phi_m
! ed2(:,:) = Ed2 = 1/R d PHI/d lam_m /sin I_m
! R = R_e + h_r we assume a reference height of 130 km which is also
! used in the TIEGCM code
!
! Author: A. Maute Dec 2003  am 12/16/03
!-----------------------------------------------------------------------

      integer  :: i, j, ip1f, ip2f, ip3f
      real(r8) :: coslm, r, fac, wrk
      real(r8) :: wrk1d(0:nmlon)

      r = r_e + h_r  ! earth radius + reference height [m]
!-----------------------------------------------------------------------
! ed2= Ed2 is the equatorward/downward component of the electric field, at all 
! geomagnetic grid points (central differencing)
!-----------------------------------------------------------------------
      fac = .5_r8/(dlatm*dtr*r)
!$omp parallel do private(j, i, wrk )
      do j = 1,nmlath-1		! southern hemisphere
        wrk = fac/sinIm_mag(j)
        do i = 0,nmlon
          ed2(i,j) = (potent(i,j+1) - potent(i,j-1))*wrk
        end do
      end do

!$omp parallel do private(j, i, wrk )
      do j = nmlath+1,nmlat-1	! northern hemisphere
        wrk = fac/sinIm_mag(j)
        do i = 0,nmlon
          ed2(i,j) = (potent(i,j+1) - potent(i,j-1))*wrk
        end do
      end do

!-----------------------------------------------------------------------
! Interpolate of ed2 between between -12 <= lam_m <= 12 degrees:
!-----------------------------------------------------------------------
      wrk1d(:) = ed2(:,jmax) - ed2(:,jmin)
      do j = jmin+1,jmax-1
        fac = (ylatm(j) - ylatm(jmin))/(ylatm(jmax) - ylatm(jmin))
        do i = 0,nmlon
	    ed2(i,j) = ed2(i,jmin) + fac*wrk1d(i)
	end do
      end do

!-----------------------------------------------------------------------
! ed1= Ed1 is the zonal component of the electric field, at all 
! geomagnetic grid points (central differencing)
!-----------------------------------------------------------------------
      fac = .5_r8/(dlonm*dtr*r)
!$omp parallel do private(j, i, wrk, coslm )
      do j = 1,nmlat-1
        coslm = ylatm(j) - 90._r8
        coslm = cos( coslm*dtr )
        wrk = fac/coslm
        do i = 1,nmlon-1
          ed1(i,j) = -(potent(i+1,j) - potent(i-1,j))*wrk
        end do
	i = 0
	ed1(i,j)     = -(potent(i+1,j) - potent(nmlon-1,j))*wrk
	ed1(nmlon,j) = ed1(i,j)
      end do

!-----------------------------------------------------------------------
! Poles:
!-----------------------------------------------------------------------
      do i = 0,nmlon
        ip1f = i + nmlon/4
        if( ip1f > nmlon ) then
           ip1f = ip1f - nmlon
        end if
        ip2f = i + nmlon/2
        if( ip2f > nmlon ) then
           ip2f = ip2f - nmlon
        end if
        ip3f = i + 3*nmlon/4
        if( ip3f > nmlon ) then
           ip3f = ip3f - nmlon
        end if
        ed1(i,0)     = .25_r8*(ed1(i,1) - ed1(ip2f,1) + ed2(ip1f,1) - ed2(ip3f,1))
        ed1(i,nmlat) = .25_r8*(ed1(i,nmlat-1) - ed1(ip2f,nmlat-1) &
                               + ed2(ip1f,nmlat-1) - ed2(ip3f,nmlat-1))
        ed2(i,0)     = .25_r8*(ed2(i,1) - ed2(ip2f,1) - ed1(ip1f,1) + ed1(ip3f,1))
        ed2(i,nmlat) = .25_r8*(ed2(i,nmlat-1) - ed2(ip2f,nmlat-1) &
                               - ed1(ip1f,nmlat-1) + ed1(ip3f,nmlat-1))
      end do

      end subroutine DerivPotential

   end module efield
