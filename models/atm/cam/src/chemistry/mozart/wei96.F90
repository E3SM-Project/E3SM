!      
! Purpose: 
! Subroutines to calculate the electric potentials from the Weimer '96 model of
! the polar cap ionospheric electric potentials.
!
! Method:
!
! To use, first call subroutine ReadCoef once.
! Next, call SetModel with the specified input parameters.
! The function EpotVal(gLAT,gMLT) can then be used repeatively to get the
! electric potential at the desired location in geomagnetic coordinates.
! Subroutines to calculate the electric potentials from the Weimer '96 model of
! the polar cap ionospheric electric potentials.
!
!
! Author: A. Maute Dec 2003  
! This code is protected by copyright and is
! distributed for research or educational use only.
! Commerical use without written permission from Dan Weimer/MRC is prohibited.
!
!*********************** Copyright 1996, Dan Weimer/MRC ***********************
!================================================================================================

	FUNCTION EpotVal(gLAT,gMLT)
!
!-----------------------------------------------------------------------
! Return the value of the electric potential in kV at
! corrected geomagnetic coordinates gLAT (degrees) and gMLT (hours).
!
! Must first call ReadCoef and SetModel to set up the model coeficients for
! the desired values of Bt, IMF clock angle, Dipole tilt angle, and SW Vel.
!-----------------------------------------------------------------------
!
        use shr_kind_mod, only: r8 => shr_kind_r8
        implicit none 
!
!-----------------------------Return Value------------------------------
!
        real(r8) EpotVal
!
!-------------------------------Commons---------------------------------
!
	INTEGER ML,MM
	REAL(r8) Coef(0:1,0:8,0:3),pi
	COMMON/SetCoef/Coef,pi,ML,MM
!
!------------------------------Arguments--------------------------------
!
	REAL(r8) gLAT,gMLT
!
!---------------------------Local variables-----------------------------
!
        integer limit,l,m

	Real(r8) Theta,Phi,Z,ct,Phim
        real(r8) r
	REAL(r8) Plm(0:20,0:20)
!
!-----------------------------------------------------------------------
!
	r=90._r8-gLAT
	IF(r .LT. 45._r8)THEN
	  Theta=r*pi/45._r8
          Phi=gMLT*pi/12._r8
	  Z=Coef(0,0,0)
	  ct=COS(Theta)
	  CALL Legendre(ct,ML,MM,Plm)
	  DO l=1,ML
	    Z=Z + Coef(0,l,0)*Plm(l,0)
	    IF(l.LT.MM)THEN
	      limit=l
	    ELSE
	      limit=MM
	    ENDIF
	    DO m=1,limit
	      phim=phi*m
              Z=Z + Coef(0,l,m)*Plm(l,m)*COS(phim) + &
     		   Coef(1,l,m)*Plm(l,m)*SIN(phim)
	    ENDDO
	  ENDDO
	ELSE
	  Z=0._r8
	ENDIF
	EpotVal=Z
	RETURN
	END FUNCTION EpotVal

!================================================================================================

	SUBROUTINE ReadCoef (wei96_file)
!
!-----------------------------------------------------------------------
!
! Read in the data file with the model coefficients
!
!*********************** Copyright 1996, Dan Weimer/MRC ***********************
!
! NCAR addition (Jan 97):  initialize constants used in GECMP
!-----------------------------------------------------------------------
!
      use shr_kind_mod,  only: r8 => shr_kind_r8
      use ioFileMod,     only : getfil
      use units,         only : getunit, freeunit
      use abortutils,    only : endrun
      use cam_logfile,   only : iulog
      implicit none 
!
!-------------------------------Commons---------------------------------
!
      real(r8) alamn, alamx, alamr, stpd, stp2, cstp, sstp
      COMMON /CECMP/ ALAMN,ALAMX,ALAMR,STPD,STP2,CSTP,SSTP
!            ALAMN = Absolute min latitude (deg) of model
!            ALAMX = Absolute max latitude (deg) for normal gradient calc.
!            STPD  = Angular dist (deg) of step @ 300km above earth (r=6371km)
!            STP2  = Denominator in gradient calc

!
!------------------------------Arguments--------------------------------
!
      character(len=*), intent(in) :: wei96_file
!
!-----------------------------Parameters------------------------------
!
      real(r8) d2r, r2d
      PARAMETER ( D2R =  0.0174532925199432957692369076847_r8 , &
                  R2D = 57.2957795130823208767981548147_r8)
!
!---------------------------Local variables-----------------------------
!
      INTEGER udat,unit,ios
      integer ll, mm, k, m, klimit, kk, nn, ii, i, n, ilimit, mlimit, l

      REAL(r8) C(0:3)
      real(r8) stpr, step

      CHARACTER*15 skip

      INTEGER MaxL,MaxM,MaxN
      REAL(r8) Cn( 0:3 , 0:1 , 0:4 , 0:1 , 0:8 , 0:3 )
      COMMON /AllCoefs/Cn,MaxL,MaxM,MaxN

      character(len=256) :: locfn
      logical, parameter :: debug = .false.
!
!-----------------------------------------------------------------------
!
      STEP = 10._r8
      STPR = STEP/6671._r8
      STPD = STPR*R2D
      STP2 = 2._r8*STEP
      CSTP = COS (STPR)
      SSTP = SQRT (1._r8 - CSTP*CSTP)
      ALAMN = 45._r8
      ALAMX = 90._r8 - STPD
      ALAMR = ALAMN*D2R
!          End NCAR addition
! 
!  get coeff_file  
      unit= getunit()
      if (debug) write(iulog,*) 'Weimer: getting file ',trim(wei96_file),' unit ',unit
      call getfil( wei96_file, locfn, 0 )
!      
      if (debug) write(iulog,*) 'Weimer: opening file ',trim(locfn),' unit ',unit	
      OPEN(unit=unit,file=trim(locfn), &
           status = 'old',iostat = ios)
      if(ios.gt.0) then
	write(iulog,*) 'Weimer: error in opening wei96.cofcnts',' unit ',unit
        call endrun
      endif
      
  900 FORMAT(A15)
 1000 FORMAT(3I8)
 2000 FORMAT(3I2)
 3000 FORMAT(2I2,4E15.6)

!     READ(udat,900) skip
      if (debug) write(iulog,*) 'Weimer: reading file ',trim(locfn),' unit ',unit	
      READ(unit,1000,iostat = ios) MaxL,MaxM,MaxN
      if(ios.gt.0) then
	write(iulog,*) 'ReadCoef: error in reading wei96.cofcnts file',' unit ',unit	
        call endrun
      endif
      DO l=0,MaxL
        IF(l.LT.MaxM)THEN
          mlimit=l
        ELSE
          mlimit=MaxM
        ENDIF
        DO m=0,mlimit
          IF(m.LT.1)THEN
            klimit=0
          ELSE
            klimit=1
          ENDIF
          DO k=0,klimit
            READ(unit,2000,iostat = ios) ll,mm,kk
            if(ios.gt.0) then
      	      write(iulog,*) 'ReadCoef: error in reading wei96.cofcnts file',' unit ',unit	
              call endrun
            endif
            IF(ll.NE.l .OR. mm.NE.m .OR. kk.NE.k)THEN
              WRITE(IULOG,*)'Data File Format Error'
              CALL ENDRUN
            ENDIF
            DO n=0,MaxN
              IF(n.LT.1)THEN
        	ilimit=0
              ELSE
        	ilimit=1
              ENDIF
              DO i=0,ilimit
        	READ(unit,3000,iostat = ios) nn,ii,C
                if(ios.gt.0) then
      	          write(iulog,*) 'ReadCoef: error in reading', &
                       ' wei96.cofcnts file',' unit ',unit	
                  call endrun
                endif
        	IF(nn.NE.n .OR. ii.NE.i)THEN
        	  WRITE(IULOG,*)'Data File Format Error'
        	  CALL ENDRUN
        	ENDIF
        	Cn(0,i,n,k,l,m)=C(0)
        	Cn(1,i,n,k,l,m)=C(1)
        	Cn(2,i,n,k,l,m)=C(2)
        	Cn(3,i,n,k,l,m)=C(3)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!      
      close(unit)
      call freeunit(unit)
!    
      RETURN
      END SUBROUTINE ReadCoef

!================================================================================================

	FUNCTION FSVal(omega,MaxN,FSC)
!
!-----------------------------------------------------------------------
! Evaluate a  Sine/Cosine Fourier series for N terms up to MaxN
! at angle omega, given the coefficients in FSC
!
!*********************** Copyright 1996, Dan Weimer/MRC ***********************
!-----------------------------------------------------------------------
!
        use shr_kind_mod, only: r8 => shr_kind_r8
        implicit none 
!
!-----------------------------Return Value------------------------------
!
        real(r8) FSVal
!
!------------------------------Arguments--------------------------------
!
	INTEGER MaxN
	REAL(r8) omega,FSC(0:1,0:*)
!
!---------------------------Local variables-----------------------------
!
	INTEGER n
	REAL(r8) Y,theta
!
!-----------------------------------------------------------------------
!
	Y=0._r8
	DO n=0,MaxN
	  theta=omega*n
	  Y=Y + FSC(0,n)*COS(theta) + FSC(1,n)*SIN(theta)
	ENDDO
	FSVal=Y
	RETURN
	END FUNCTION FSVal

!================================================================================================

	SUBROUTINE SetModel(angle,Bt,Tilt,SWVel)
!
!-----------------------------------------------------------------------
! Calculate the complete set of spherical harmonic coefficients,
! given an arbitrary IMF angle (degrees from northward toward +Y),
! magnitude Bt (nT), dipole tilt angle (degrees),
! and solar wind velocity (km/sec).
! Returns the Coef in the common block SetCoef.
!
!*********************** Copyright 1996, Dan Weimer/MRC ***********************
!-----------------------------------------------------------------------
!
        use shr_kind_mod, only: r8 => shr_kind_r8
        implicit none 
!
!-------------------------------Commons---------------------------------
!
	INTEGER MaxL,MaxM,MaxN
	REAL(r8) Cn( 0:3 , 0:1 , 0:4 , 0:1 , 0:8 , 0:3 )
	COMMON /AllCoefs/Cn,MaxL,MaxM,MaxN

	INTEGER ML,MM
	REAL(r8) Coef(0:1,0:8,0:3),pi
	COMMON/SetCoef/Coef,pi,ML,MM
!
!------------------------------Arguments--------------------------------
!
	REAL(r8) angle,Bt,Tilt,SWVel
!
!---------------------------Local variables-----------------------------
!
        integer n, k, ilimit, i, klimit, l, m, mlimit
	REAL(r8) FSC(0:1,0:4), fsval, omega, sintilt
!
!-----------------------------------------------------------------------
!
	pi=2._r8*ASIN(1._r8)
	ML=MaxL
	MM=MaxM
	SinTilt=SIN(Tilt*pi/180._r8)
!	SinTilt=SIND(Tilt)

	omega=angle*pi/180._r8

        fsc(1,0) = 0._r8
	DO l=0,MaxL
	  IF(l.LT.MaxM)THEN
	    mlimit=l
	  ELSE
	    mlimit=MaxM
	  ENDIF
	  DO m=0,mlimit
	    IF(m.LT.1)THEN
	      klimit=0
	    ELSE
	      klimit=1
	    ENDIF
	    DO k=0,klimit
! Retrieve the regression coefficients and evaluate the function
! as a function of Bt,Tilt,and SWVel to get each Fourier coefficient.
	      DO n=0,MaxN
	        IF(n.LT.1)THEN
	          ilimit=0
	        ELSE
	          ilimit=1
	        ENDIF
		DO i=0,ilimit
		  FSC(i,n)=Cn(0,i,n,k,l,m) + Bt*Cn(1,i,n,k,l,m) + &
     		   SinTilt*Cn(2,i,n,k,l,m) + SWVel*Cn(3,i,n,k,l,m)
		ENDDO
	      ENDDO
! Next evaluate the Fourier series as a function of angle.
      	      Coef(k,l,m)=FSVal(omega,MaxN,FSC)
	    ENDDO
	  ENDDO
	ENDDO
	RETURN
	END SUBROUTINE SetModel

!================================================================================================

	SUBROUTINE LEGENDRE(x,lmax,mmax,Plm)
!
!-----------------------------------------------------------------------
! compute Associate Legendre Function P_l^m(x)
! for all l up to lmax and all m up to mmax.
! returns results in array Plm
! if X is out of range ( abs(x)>1 ) then value is returned as if x=1.
!
!*********************** Copyright 1996, Dan Weimer/MRC ***********************
!-----------------------------------------------------------------------
!
        use shr_kind_mod, only: r8 => shr_kind_r8
        use cam_logfile,  only : iulog

        implicit none 
!
!------------------------------Arguments--------------------------------
!
        integer lmax, mmax
	real(r8) x, Plm(0:20,0:20)
!
!---------------------------Local variables-----------------------------
!
        integer m, lm2, l
        real(r8) xx, fact
!
!-----------------------------------------------------------------------
!
	  DO l=0,20
	    DO m=0,20
		Plm(l,m)=0._r8
	    ENDDO
	  ENDDO
	xx=MIN(x,1._r8)
	xx=MAX(xx,-1._r8)
	IF(lmax .LT. 0 .OR. mmax .LT. 0 .OR. mmax .GT. lmax )THEN
	  write(iulog,*)'Bad arguments to Legendre'
	  RETURN
	ENDIF
! First calculate all Pl0 for l=0 to l
	Plm(0,0)=1._r8
	IF(lmax.GT.0)Plm(1,0)=xx
	IF (lmax .GT. 1 )THEN
	  DO L=2,lmax
	    Plm(L,0)=( (2._r8*L-1)*xx*Plm(L-1,0) - (L-1)*Plm(L-2,0) )/L
	  ENDDO
	ENDIF
	IF (mmax .EQ. 0 )RETURN
	fact=SQRT( (1._r8-xx)*(1._r8+xx) )
	DO M=1,mmax
	  DO L=m,lmax
	    lm2=MAX(L-2,0)
	    Plm(L,M)=Plm(lm2,M) - ( 2*L-1)*fact*Plm(L-1,M-1)
	  ENDDO
	ENDDO
	RETURN
	END SUBROUTINE LEGENDRE

!================================================================================================

!*********************** Copyright 1996, Dan Weimer/MRC ***********************

!CC NCAR MODIFIED (3/96) CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  The following routines (translib.for) were added to return the dipole tilt. C
!  GET_TILT was initially a procedure (TRANS), here it has been changed into   C
!  a function which returns the dipole tilt.                                   C
! Barbara Emery (emery@ncar.ucar.edu) and William Golesorkhi, HAO/NCAR (3/96)  C
!CC NCAR MODIFIED (3/96) CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! COORDINATE TRANSFORMATION UTILITIES
!**********************************************************************        
	FUNCTION GET_TILT(YEAR,MONTH,DAY,HOUR)
!
!-----------------------------------------------------------------------
!CC NCAR MODIFIED (3/96) CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  The following line initially was:                                           C
!       SUBROUTINE TRANS(YEAR,MONTH,DAY,HOUR,IDBUG)                            C
!  It has been changed to return the dipole tilt from this function call.      C
!CC NCAR MODIFIED (3/96) CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!         
!      THIS SUBROUTINE DERIVES THE ROTATION MATRICES AM(I,J,K) FOR 11
!      TRANSFORMATIONS, IDENTIFIED BY K.
!          K=1 TRANSFORMS GSE to GEO
!          K=2     "      GEO to MAG
!          K=3     "      GSE to MAG
!          K=4     "      GSE to GSM
!          K=5     "      GEO to GSM
!          K=6     "      GSM to MAG
!          K=7     "      GSE to GEI
!          K=8     "      GEI to GEO
!          K=9     "      GSM to SM 
!	   K=10    "      GEO to SM 
!	   K=11    "      MAG to SM 
!
!      IF IDBUG IS NOT 0, THEN OUTPUTS DIAGNOSTIC INFORMATION TO
!      FILE UNIT=IDBUG
!
!      The formal names of the coordinate systems are:
!	GSE - Geocentric Solar Ecliptic
!	GEO - Geographic
!	MAG - Geomagnetic
!	GSM - Geocentric Solar Magnetospheric
!	SM  - Solar Magnetic
!	
!      THE ARRAY CX(I) ENCODES VARIOUS ANGLES, STORED IN DEGREES
!      ST(I) AND CT(I) ARE SINES & COSINES.       
!
!      Program author:  D. R. Weimer
!
!      Some of this code has been copied from subroutines which had been
!      obtained from D. Stern, NASA/GSFC.  Other formulas are from "Space 
!      Physics Coordinate Transformations: A User Guide" by M. Hapgood (1991).
!
!      The formulas for the calculation of Greenwich mean sidereal time (GMST)
!      and the sun's location are from "Almanac for Computers 1990",
!      U.S. Naval Observatory.
!
!-----------------------------------------------------------------------
!
        use shr_kind_mod, only: r8 => shr_kind_r8
        implicit none 
!
!-----------------------------Return Value------------------------------
!
        real(r8) get_tilt
!
!-------------------------------Commons---------------------------------
!
        real(r8) cx, st, ct, am
	COMMON/TRANSDAT/CX(9),ST(6),CT(6),AM(3,3,11)

        real(r8) :: th0=11.19_r8
        real(r8) :: ph0=-70.76_r8
!
!------------------------------Arguments--------------------------------
!
	INTEGER YEAR, MONTH, DAY
	REAL(r8) HOUR
!
!-----------------------------Parameters------------------------------
!
	INTEGER GSEGEO,GEOGSE,GEOMAG,MAGGEO
	INTEGER GSEMAG,MAGGSE,GSEGSM,GSMGSE
	INTEGER GEOGSM,GSMGEO,GSMMAG,MAGGSM
	INTEGER GSEGEI,GEIGSE,GEIGEO,GEOGEI
	INTEGER GSMSM,SMGSM,GEOSM,SMGEO,MAGSM,SMMAG

	PARAMETER (GSEGEO= 1,GEOGSE=-1,GEOMAG= 2,MAGGEO=-2)
	PARAMETER (GSEMAG= 3,MAGGSE=-3,GSEGSM= 4,GSMGSE=-4)
	PARAMETER (GEOGSM= 5,GSMGEO=-5,GSMMAG= 6,MAGGSM=-6)
	PARAMETER (GSEGEI= 7,GEIGSE=-7,GEIGEO= 8,GEOGEI=-8)
	PARAMETER (GSMSM = 9,SMGSM =-9,GEOSM =10,SMGEO=-10)
	PARAMETER (MAGSM =11,SMMAG =-11)
!
!---------------------------Local variables-----------------------------
!
        integer IDBUG
        integer j, k, jd, iyr, i, mjd

        REAL(r8) UT, T0, GMSTD, GMSTH, ECLIP, MA, LAMD, SUNLON, pi
        real(r8) b32, b33, b3
!
!-------------------------External Functions----------------------------
!
        integer julday
        external julday
!
!-----------------------------------------------------------------------
!
	pi=2._r8*ASIN(1._r8)
!CC NCAR MODIFICATION (3/96) CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! IDBUG=0 to prevent printing data to the screen or writing data to a file.    C
        IDBUG = 0
!CC NCAR MODIFICATION (3/96) CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	IF(YEAR.LT.1900)THEN
	  IYR=1900+YEAR
	ELSE
	  IYR=YEAR
	ENDIF
	UT=HOUR
	JD=JULDAY(MONTH,DAY,IYR)
	MJD=JD-2400001
	T0=(real(MJD,r8)-51544.5_r8)/36525.0_r8
	GMSTD=100.4606184_r8 + 36000.770_r8*T0 + 3.87933E-4_r8*T0*T0 + &
             15.0410686_r8*UT
	CALL ADJUST(GMSTD)
	GMSTH=GMSTD*24._r8/360._r8
	ECLIP=23.439_r8 - 0.013_r8*T0
        MA=357.528_r8 + 35999.050_r8*T0 + 0.041066678_r8*UT
        CALL ADJUST(MA)
        LAMD=280.460_r8 + 36000.772_r8*T0 + 0.041068642_r8*UT
        CALL ADJUST(LAMD)
        SUNLON=LAMD + (1.915_r8-0.0048_r8*T0)*SIN(MA*pi/180._r8) + 0.020_r8* &
           SIN(2._r8*MA*pi/180._r8)
        CALL ADJUST(SUNLON)
	IF(IDBUG.NE.0)THEN
	  WRITE(IDBUG,*) YEAR,MONTH,DAY,HOUR
	  WRITE(IDBUG,*) 'MJD=',MJD
	  WRITE(IDBUG,*) 'T0=',T0
	  WRITE(IDBUG,*) 'GMSTH=',GMSTH
	  WRITE(IDBUG,*) 'ECLIPTIC OBLIQUITY=',ECLIP
	  WRITE(IDBUG,*) 'MEAN ANOMALY=',MA
	  WRITE(IDBUG,*) 'MEAN LONGITUDE=',LAMD
	  WRITE(IDBUG,*) 'TRUE LONGITUDE=',SUNLON
	ENDIF

	CX(1)= GMSTD
	CX(2) = ECLIP
	CX(3) = SUNLON
	CX(4) = TH0
	CX(5) = PH0
! Derived later:
!       CX(6) = Dipole tilt angle  
!       CX(7) = Angle between sun and magnetic pole
!       CX(8) = Subsolar point latitude
!       CX(9) = Subsolar point longitude

	DO I=1,5
	  ST(I) = SIN(CX(I)*pi/180._r8)
	  CT(I) = COS(CX(I)*pi/180._r8)
	ENDDO
!         
      AM(1,1,GSEGEI) = CT(3)
      AM(1,2,GSEGEI) = -ST(3)
      AM(1,3,GSEGEI) = 0._r8         
      AM(2,1,GSEGEI) = ST(3)*CT(2)
      AM(2,2,GSEGEI) = CT(3)*CT(2)
      AM(2,3,GSEGEI) = -ST(2)
      AM(3,1,GSEGEI) = ST(3)*ST(2)
      AM(3,2,GSEGEI) = CT(3)*ST(2)
      AM(3,3,GSEGEI) = CT(2)      
!         
      AM(1,1,GEIGEO) = CT(1)      
      AM(1,2,GEIGEO) = ST(1)      
      AM(1,3,GEIGEO) = 0._r8         
      AM(2,1,GEIGEO) = -ST(1)     
      AM(2,2,GEIGEO) = CT(1)      
      AM(2,3,GEIGEO) = 0._r8         
      AM(3,1,GEIGEO) = 0._r8         
      AM(3,2,GEIGEO) = 0._r8         
      AM(3,3,GEIGEO) = 1._r8         
!         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,GSEGEO) = AM(I,1,GEIGEO)*AM(1,J,GSEGEI) + &
         AM(I,2,GEIGEO)*AM(2,J,GSEGEI) +  AM(I,3,GEIGEO)*AM(3,J,GSEGEI)
      ENDDO
      ENDDO
!         
      AM(1,1,GEOMAG) = CT(4)*CT(5) 
      AM(1,2,GEOMAG) = CT(4)*ST(5) 
      AM(1,3,GEOMAG) =-ST(4)       
      AM(2,1,GEOMAG) =-ST(5)       
      AM(2,2,GEOMAG) = CT(5)       
      AM(2,3,GEOMAG) = 0._r8
      AM(3,1,GEOMAG) = ST(4)*CT(5) 
      AM(3,2,GEOMAG) = ST(4)*ST(5) 
      AM(3,3,GEOMAG) = CT(4)       
!         
      DO I=1,3   
      DO J=1,3   
       AM(I,J,GSEMAG) = AM(I,1,GEOMAG)*AM(1,J,GSEGEO) + &
        AM(I,2,GEOMAG)*AM(2,J,GSEGEO) +  AM(I,3,GEOMAG)*AM(3,J,GSEGEO)
      ENDDO
      ENDDO
!         
      B32 = AM(3,2,GSEMAG)         
      B33 = AM(3,3,GSEMAG)         
      B3  = SQRT(B32*B32+B33*B33)       
      IF (B33.LE.0._r8) B3 = -B3  
!         
      AM(2,2,GSEGSM) = B33/B3      
      AM(3,3,GSEGSM) = AM(2,2,GSEGSM)   
      AM(3,2,GSEGSM) = B32/B3      
      AM(2,3,GSEGSM) =-AM(3,2,GSEGSM)   
      AM(1,1,GSEGSM) = 1._r8
      AM(1,2,GSEGSM) = 0._r8
      AM(1,3,GSEGSM) = 0._r8
      AM(2,1,GSEGSM) = 0._r8
      AM(3,1,GSEGSM) = 0._r8
!         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,GEOGSM) = AM(I,1,GSEGSM)*AM(J,1,GSEGEO) + &
         AM(I,2,GSEGSM)*AM(J,2,GSEGEO) + AM(I,3,GSEGSM)*AM(J,3,GSEGEO)
      ENDDO
      ENDDO
!         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,GSMMAG) = AM(I,1,GEOMAG)*AM(J,1,GEOGSM) + &
        AM(I,2,GEOMAG)*AM(J,2,GEOGSM) + AM(I,3,GEOMAG)*AM(J,3,GEOGSM)
      ENDDO
      ENDDO 
!
	ST(6) = AM(3,1,GSEMAG)        
	CT(6) = SQRT(1._r8-ST(6)*ST(6))      
	CX(6) = ASIN(ST(6)*pi/180._r8)  

        AM(1,1,GSMSM) = CT(6)
        AM(1,2,GSMSM) = 0._r8
        AM(1,3,GSMSM) = -ST(6)
        AM(2,1,GSMSM) = 0._r8
        AM(2,2,GSMSM) = 1._r8
        AM(2,3,GSMSM) = 0._r8
        AM(3,1,GSMSM) = ST(6)
        AM(3,2,GSMSM) = 0._r8
        AM(3,3,GSMSM) = CT(6)  
!         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,GEOSM) = AM(I,1,GSMSM)*AM(1,J,GEOGSM) + &
         AM(I,2,GSMSM)*AM(2,J,GEOGSM) +  AM(I,3,GSMSM)*AM(3,J,GEOGSM)
      ENDDO
      ENDDO
!         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,MAGSM) = AM(I,1,GSMSM)*AM(J,1,GSMMAG) + &
        AM(I,2,GSMSM)*AM(J,2,GSMMAG) + AM(I,3,GSMSM)*AM(J,3,GSMMAG)
      ENDDO
      ENDDO
      
!
      CX(7)=ATAN2( AM(2,1,11) , AM(1,1,11) )
      
      CX(7)=CX(7)*180._r8/pi
      CX(8)=ASIN( AM(3,1,1)*pi/180._r8 )
      CX(9)=ATAN2( AM(2,1,1) , AM(1,1,1) )
      CX(9)=CX(9)*180._r8/pi

      IF(IDBUG.NE.0)THEN
	  WRITE(IDBUG,*) 'Dipole tilt angle=',CX(6)
	  WRITE(IDBUG,*) 'Angle between sun and magnetic pole=',CX(7)
	  WRITE(IDBUG,*) 'Subsolar point latitude=',CX(8)
	  WRITE(IDBUG,*) 'Subsolar point longitude=',CX(9)

        DO K=1,11
         WRITE(IDBUG,1001) K
         DO I=1,3
           WRITE(IDBUG,1002) (AM(I,J,K),J=1,3)
         ENDDO
        ENDDO
 1001   FORMAT(' ROTATION MATRIX ',I2)
 1002   FORMAT(3F9.5)
      ENDIF

!CC NCAR MODIFICATION (3/96) CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  The next line was added to return the dipole tilt from this function call.  C
      GET_TILT = CX(6)
!CC NCAR MODIFICATION (3/96) CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      RETURN
      END FUNCTION GET_TILT

!================================================================================================

      SUBROUTINE ROTATE (X,Y,Z,I)       
!
!-----------------------------------------------------------------------
!     THIS SUBROUTINE APPLIES TO THE VECTOR (X,Y,Z) THE ITH ROTATION  
!     MATRIX AM(N,M,I) GENERATED BY SUBROUTINE TRANS
!     IF I IS NEGATIVE, THEN THE INVERSE ROTATION IS APPLIED
!-----------------------------------------------------------------------
!
      use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none 
!
!------------------------------Arguments--------------------------------
!
      integer i
      REAL(r8) X,Y,Z
!
!---------------------------Local variables-----------------------------
!
      REAL(r8) A(3)
!
!-----------------------------------------------------------------------
!
      A(1)=X
      A(2)=Y
      A(3)=Z
      CALL ROTATEV(A,A,I)
      X=A(1)
      Y=A(2)
      Z=A(3)
    
      RETURN        
      END SUBROUTINE ROTATE

!================================================================================================

      SUBROUTINE ROTATEV (A,B,I)       
!         
!-----------------------------------------------------------------------
!     THIS SUBROUTINE APPLIES TO THE VECTOR A(3) THE ITH ROTATION  
!     MATRIX AM(N,M,I) GENERATED BY SUBROUTINE TRANS
!     AND OUTPUTS THE CONVERTED VECTOR B(3), WITH NO CHANGE TO A.
!     IF I IS NEGATIVE, THEN THE INVERSE ROTATION IS APPLIED
!-----------------------------------------------------------------------
!
      use shr_kind_mod, only: r8 => shr_kind_r8
      use cam_logfile,  only : iulog
      use abortutils,   only : endrun

      implicit none 
!
!-------------------------------Commons---------------------------------
!
      real(r8) cx, st, ct, am
      COMMON/TRANSDAT/CX(9),ST(6),CT(6),AM(3,3,11)
!
!------------------------------Arguments--------------------------------
!
      integer i
      REAL(r8) A(3),B(3)
!
!---------------------------Local variables-----------------------------
!
      integer id, j
      real(r8) xa, ya, za
!
!-----------------------------------------------------------------------
!
      IF(I.EQ.0 .OR. IABS(I).GT.11)THEN
      WRITE(IULOG,*)'ROTATEV CALLED WITH UNDEFINED TRANSFORMATION'
      CALL ENDRUN
      ENDIF

      XA = A(1)
      YA = A(2)
      ZA = A(3)
      IF(I.GT.0)THEN
	ID=I
        DO J=1,3
          B(J) = XA*AM(J,1,ID) + YA*AM(J,2,ID) + ZA*AM(J,3,ID)
        ENDDO
      ELSE
	ID=-I
        DO J=1,3
          B(J) = XA*AM(1,J,ID) + YA*AM(2,J,ID) + ZA*AM(3,J,ID)
        ENDDO
      ENDIF
      RETURN        
      END SUBROUTINE ROTATEV

!================================================================================================

	SUBROUTINE FROMCART(R,LAT,LONG,POS)
!
!-----------------------------------------------------------------------
! CONVERT CARTESIAN COORDINATES POS(3)
! TO SPHERICAL COORDINATES R, LATITUDE, AND LONGITUDE (DEGREES)
!-----------------------------------------------------------------------
!
        use shr_kind_mod, only: r8 => shr_kind_r8
        implicit none 
!
!------------------------------Arguments--------------------------------
!
	REAL(r8) R, LAT, LONG, POS(3)
!
!---------------------------Local variables-----------------------------
!
        real(r8) pi
!
!-----------------------------------------------------------------------
!
	pi=2._r8*ASIN(1._r8)
	R=SQRT(POS(1)*POS(1) + POS(2)*POS(2) + POS(3)*POS(3))
	IF(R.EQ.0._r8)THEN
	  LAT=0._r8
	  LONG=0._r8
	ELSE
	  LAT=ASIN(POS(3)*pi/180._r8/R)
	  LONG=ATAN2(POS(2),POS(1))
	  LONG=LONG*180._r8/pi
	ENDIF
	RETURN
	END SUBROUTINE FROMCART

!================================================================================================

	SUBROUTINE TOCART(R,LAT,LONG,POS)
!
!-----------------------------------------------------------------------
! CONVERT SPHERICAL COORDINATES R, LATITUDE, AND LONGITUDE (DEGREES)
! TO CARTESIAN COORDINATES POS(3)
!-----------------------------------------------------------------------
!
        use shr_kind_mod, only: r8 => shr_kind_r8
        implicit none 
!
!------------------------------Arguments--------------------------------
!
	REAL(r8) R, LAT, LONG, POS(3)
!
!---------------------------Local variables-----------------------------
!
        real(r8) pi, stc, ctc, sf, cf
!
!-----------------------------------------------------------------------
!
	pi=2._r8*ASIN(1._r8)
        STC = SIN(LAT*pi/180._r8)    
        CTC = COS(LAT*pi/180._r8)    
        SF = SIN(LONG*pi/180._r8)     
        CF = COS(LONG*pi/180._r8)     
        POS(1) = R*CTC*CF        
        POS(2) = R*CTC*SF        
        POS(3) = R*STC
	RETURN
	END SUBROUTINE TOCART

!================================================================================================

	SUBROUTINE ADJUST(ANGLE)
!
!-----------------------------------------------------------------------
!	ADJUST AN ANGLE IN DEGREES TO BE IN RANGE OF 0 TO 360.
!-----------------------------------------------------------------------
!
        use shr_kind_mod, only: r8 => shr_kind_r8
        implicit none 
!
!------------------------------Arguments--------------------------------
!
        real(r8) angle
!
!-----------------------------------------------------------------------
!
 10	CONTINUE
	IF(ANGLE.LT.0._r8)THEN
	  ANGLE=ANGLE+360._r8
	  GOTO 10
	ENDIF
 20	CONTINUE
 	IF(ANGLE.GE.360._r8)THEN
	  ANGLE=ANGLE-360._r8
	  GOTO 20
	ENDIF
	RETURN
	END SUBROUTINE ADJUST

!================================================================================================

      INTEGER FUNCTION JULDAY(MM,ID,IYYY)
!
!-----------------------------------------------------------------------
!
      use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none 
!
!------------------------------Arguments--------------------------------
!
      integer mm, id, iyyy
!
!-----------------------------Parameters------------------------------
!
      integer igreg
      PARAMETER (IGREG=15+31*(10+12*1582))
!
!---------------------------Local variables-----------------------------
!
      integer ja, jm, jy
!
!-----------------------------------------------------------------------
!
      IF (IYYY.EQ.0) PAUSE 'There is no Year Zero.'
      IF (IYYY.LT.0) IYYY=IYYY+1
      IF (MM.GT.2) THEN
        JY=IYYY
        JM=MM+1
      ELSE
        JY=IYYY-1
        JM=MM+13
      ENDIF
      JULDAY=INT(365.25_r8*JY)+INT(30.6001_r8*JM)+ID+1720995
      IF (ID+31*(MM+12*IYYY).GE.IGREG) THEN
        JA=INT(0.01_r8*JY)
        JULDAY=JULDAY+2-JA+INT(0.25_r8*JA)
      ENDIF
      RETURN
      END FUNCTION JULDAY

!================================================================================================

	FUNCTION MLT(MagLong)
!
!-----------------------------------------------------------------------
! given magnetic longitude in degrees, return Magnetic Local Time
! assuming that TRANS has been called with the date & time to calculate
! the rotation matrices.
!
! btf 11/06/03:
! Call sub adjust instead of referencing it as a function
!-----------------------------------------------------------------------
!
        use shr_kind_mod, only: r8 => shr_kind_r8
        implicit none 
!
!-----------------------------Return Value------------------------------
!
        real(r8) mlt
!
!-------------------------------Commons---------------------------------
!
        real(r8) cx, st, ct, am
	COMMON/TRANSDAT/CX(9),ST(6),CT(6),AM(3,3,11)

!
!------------------------------Arguments--------------------------------
!
	REAL(r8) MagLong
!
!---------------------------Local variables-----------------------------
!
	REAL(r8) angle, rotangle
!
!-----------------------------------------------------------------------
!
	RotAngle=CX(7)
!       MLT=ADJUST(Maglong+RotAngle+180.)/15.
        angle = Maglong+RotAngle+180._r8
        call adjust(angle)
        mlt = angle/15._r8
	RETURN
	END FUNCTION MLT

!================================================================================================

	FUNCTION MagLong(MLT)
!
!-----------------------------------------------------------------------
! return magnetic longitude in degrees, given Magnetic Local Time
! assuming that TRANS has been called with the date & time to calculate
! the rotation matrices.
!
! btf 11/06/03:
! Call sub adjust instead of referencing it as a function
!-----------------------------------------------------------------------
!
        use shr_kind_mod, only: r8 => shr_kind_r8
        implicit none 
!
!-----------------------------Return Value------------------------------
!
        real(r8) MagLong
!
!-------------------------------Commons---------------------------------
!
        real(r8) cx, st, ct, am
        COMMON/TRANSDAT/CX(9),ST(6),CT(6),AM(3,3,11)
!
!------------------------------Arguments--------------------------------
!
	REAL(r8) MLT
!
!---------------------------Local variables-----------------------------
!
	REAL(r8) angle, rotangle
!
!-----------------------------------------------------------------------
!
	RotAngle=CX(7)
	angle=MLT*15._r8-RotAngle-180._r8
!       MagLong=ADJUST(angle)
        call adjust(angle)
        MagLong = angle
	RETURN
	END FUNCTION MagLong

!================================================================================================

	SUBROUTINE SunLoc(SunLat,SunLong)
!
!-----------------------------------------------------------------------
! Return latitude and longitude of sub-solar point.
! Assumes that TRANS has previously been called with the
! date & time to calculate the rotation matrices.
!-----------------------------------------------------------------------
!
        use shr_kind_mod, only: r8 => shr_kind_r8
        implicit none 
!
!-------------------------------Commons---------------------------------
!
        real(r8) cx, st, ct, am
        COMMON/TRANSDAT/CX(9),ST(6),CT(6),AM(3,3,11)
!
!------------------------------Arguments--------------------------------
!
	Real(r8) SunLat,SunLong
!
!-----------------------------------------------------------------------
!
	SunLong=CX(9)
	SunLat=CX(8)
	RETURN
	END SUBROUTINE SunLoc

!================================================================================================

      SUBROUTINE GECMP (AMLA,RMLT,ET,EP)
!
!-----------------------------------------------------------------------
!          Get Electric field components for the Weimer electrostatic
!          potential model.  Before use, first load coefficients (CALL
!          READCOEF) and initialize model conditions (CALL SETMODEL).
!
!          INPUTS:
!            AMLA = Absolute value of magnetic latitude (deg)
!            RMLT = Magnetic local time (hours).
!          RETURNS:
!            ET = Etheta (magnetic equatorward*) E field component (V/m)
!            EP = Ephi   (magnetic eastward)     E field component (V/m)
!
!          * ET direction is along the magnetic meridian away from the
!            current hemisphere; i.e., when ET > 0, the direction is
!              southward when RMLA > 0
!              northward when RMLA < 0
!
!          NCAR addition (Jan 97).  R.Barnes
!-----------------------------------------------------------------------
!
      use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none 
!
!-------------------------------Commons---------------------------------
!
!          CECMP contains constants initialized in READCOEF
      real(r8) alamn, alamx, alamr, stpd, stp2, cstp, sstp
      COMMON /CECMP/ ALAMN,ALAMX,ALAMR,STPD,STP2,CSTP,SSTP
!
!------------------------------Arguments--------------------------------
!
      real(r8) amla, rmlt, et, ep
!
!-----------------------------Parameters------------------------------
!
      real(r8) d2r, r2d
      PARAMETER ( D2R =  0.0174532925199432957692369076847_r8 , &
                 R2D = 57.2957795130823208767981548147_r8)
!
!---------------------------Local variables-----------------------------
!
      real(r8) p1, p2
      real(r8) xmlt, xmlt1, kpol, dphi, amla1
!
!-------------------------External Functions----------------------------
!
      real(r8) epotval
      external epotval
!
!-----------------------------------------------------------------------
!
      ET = -99999._r8
      EP = -99999._r8
      IF (AMLA .LT. 0._r8) GO TO 100

!          Calculate -(latitude gradient) by stepping 10 km along the
!          meridian in each direction (flipping coordinates when going
!          over pole to keep lat <= 90).
      KPOL  = 0
      XMLT  = RMLT
   10 XMLT1 = XMLT
      AMLA1 = AMLA + STPD
      IF (AMLA1 .GT. 90._r8) THEN
	AMLA1 = 180._r8 - AMLA1
	XMLT1 = XMLT1 + 12._r8
      ENDIF
      P1 = EPOTVAL (AMLA1    ,XMLT1)
      P2 = EPOTVAL (AMLA-STPD,XMLT )
      IF (KPOL .EQ. 1) GO TO 20
      ET = (P1 - P2) / STP2

!          Calculate -(lon gradient).  For most latitudes, step along a
!          great circle.  However, limit minimum latitude to the model
!          minimum (distorting the path onto a latitude line).  Also,
!          avoid a divide by zero at the pole avoid by using Art's trick
!          where Ephi(90,lon) = Etheta(90,lon+90)
      IF (AMLA .LT. ALAMX) THEN
	AMLA1 = MAX (ASIN(SIN(AMLA*D2R)*CSTP) , ALAMR)
	DPHI  = ASIN (SSTP/SIN(AMLA1))*R2D
	AMLA1 = AMLA1*R2D
	P1 = EPOTVAL (AMLA1,XMLT+DPHI)
	P2 = EPOTVAL (AMLA1,XMLT-DPHI)
      ELSE
	AMLA = 90._r8
	XMLT = XMLT + 6._r8
	KPOL = 1
	GO TO 10
      ENDIF
   20 EP = (P2 - P1) / STP2
      IF (KPOL .EQ. 1) EP = -EP

!          Below model minimum lat, the potential is value at min lat
      IF (AMLA .LT. ALAMN) THEN
	ET = 0._r8
	EP = EP * COS(ALAMR)/COS(AMLA*D2R)
      ENDIF

  100 RETURN
      END SUBROUTINE GECMP

!================================================================================================
