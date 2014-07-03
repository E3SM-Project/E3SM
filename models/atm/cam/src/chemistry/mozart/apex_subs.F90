!
! Original file ~bozo/pgms/apex/apxntrp.f copied 2/25/00.
!
!================================================================================================

      SUBROUTINE APXMKA (MSGUN, EPOCH, GPLAT,GPLON,GPALT,NLAT,NLON, &
                         NALT,WK,LWK, IST)
!
!-----------------------------------------------------------------------
!          This contains ENTRYs that determine quantities related to Apex
!          coordinates.  They are designed to speed up coordinate determination
!          by substituting interpolation from previously calculated tables
!          for direct calculation.  They also permit calculating quantities
!	   involving gradients of Apex coordinates, such as base vectors in
!	   the Modified-Apex and Quasi-Dipole systems.  Options allow table 
!	   creation, storage, and read-back in addition to interpolation.
!
!          Tables (arrays) must be prepared prior to coordinate determination.
!          Tables may be created for one date and held only in memory by calling
!
!            APXMKA - make magnetic arrays.
!
!          Or, they may be created and written, then read back later by
!
!            APXWRA - make and write arrays
!            APXRDA - read stored arrays.
!
!          Tables may be created for multiple times (e.g., the DGRF dates)
!          and written to a file using APXWRA.  However, tables are held in
!          memory for only one time due to the potential to exceed memory
!          capacity.  Thus, when reading back multiples times, APXRDA
!          interpolates/extrapolates to a single time.
!
!          Alternatively, If APXWRA is called to make and write tables for one
!          time only, then it is not necessary to call APXRDA because tables
!          have already been initialized for that time.
!
!          Creation of global lookup tables can be time consuming, thus, it
!          is expected that this is done rarely; it might be considered an
!          installation step.  Once the tables are established, programs may
!          reference them by calling APXRDA.  In this case, the grid coordinates
!          (originally specified during the table creation) may be ascertained
!          by calling
!
!            APXGGC - get grid coordinates.
!
!          After the tables have been loaded in memory for a given date, it
!          is possible to interpolate in space to determinine various magnetic
!          parameters using the following calls:
!
!            APXALL  - get Apex coordinates (radius, lat, lon).
!            APXMALL - get Modified Apex coordinates, Quasi-Dipole
!                      coordinates and associated base vectors.
!            APXQ2G  - convert quasi-dipole to geodetic coordinates.
!
!          Details for each ENTRY introduced above are provided below.
!          Following the ENTRY summaries are comments describing EXTERNALS,
!          INSTALLATION SPECIFICS for different computers, and a HISTORY.
!
!          REFERENCE:
!          Richmond, A. D., Ionospheric Electrodynamics Using Magnetic Apex
!          Coordinates, J. Geomag. Geoelectr., 47, 191-212, 1995.
!
!          ALGORITHM:
!          When arrays are created, APXMKA calls subroutine MAKEXYZV, which
!          in turn calls subroutine APEX at each grid point to get the apex
!          radius A, the apex longitude PHIA, and the magnetic potential
!          VMP.  The cosine (CLP) and sine (SLP) of the quasi-dipole
!          latitude are computed from A.  From these can be computed
!          preliminary quantites defined as
!
!            x = cos(quasi-dipole latitude)*cos(apex longitude)
!            y = cos(quasi-dipole latitude)*sin(apex longitude)
!            z = sin(quasi-dipole latitude)
!            v = (VMP/VP)*((RE+ALT)/RE)**2
!
!          where VP is the magnitude of the magnetic potential of the
!          geomagnetic dipole at a radius of RE; ALT is altitude; and RE is
!          the mean Earth radius.  Note that all of these quantities vary
!          smoothly across the apex poles, unlike apex latitude or
!          quasi-dipole latitude, so that they can be linearly interpolated
!          near these poles.  Corresponding values of x,y,z,v for a dipole
!          field on a spherical Earth are computed analytically and
!          subtracted from the above quantities, and the results are put
!          into the 3D arrays X,Y,Z,V.  When APXALL or APXMALL is called,
!          trilinear interpolations (in latitude, longitude, and inverse
!          altitude) are carried out between the grid points.  Gradients are
!          calculated for APXMALL in order to determine the base vectors.
!          Analytic formulas appropriate for a dipole field on a spherical
!          Earth are used to determine the previously removed dipole
!          components of x,y,z,v, and their gradients and these are added
!          back to the interpolated values obtained for the non-dipole
!          components.  Finally, the apex-based coordinates and their base
!          vectors are calculated from x,y,z,v and their gradients.
!
!------------------------------------------------------------------------------
!
!       ENTRY APXMKA and APXWRA:
!          These create gridded arrays of quantities that can later be
!          linearly interpolated to any desired geographic point within
!          their range by APXALL, APXMALL, and APXQ2G.  The spatial extent
!          and resolution of the arrays in latitude, longitude, and altitude
!          can be tailored for the application at hand.  In one extreme,
!          very high interpolation accuracy can be achieved by creating a
!          2x2x2 array with very small spatial increments surrounding the
!          point of interest. (However, since vector quantities are obtained
!          by taking differences of the quantities at adjacent points,
!          computer roundoff errors will limit the benefits of making the
!          increments too small.) In another extreme, a global array from
!          the ground to an altitude of several Earth radii can be created.
!          (In this case, the spatial resolution is limited by the need to
!          maintain reasonable array sizes.) For most purposes, we have
!          found it adequate to use a global grid with dimensions 91,151,6
!          (latitude,longitude,altitude) and altitude from ground to 1274 km
!          as created by GGRID with NVERT = 30 (see file test.f).
!
!          WARNING: execution time to create these look-up tables can be
!          substantial.  It is dependent on the grid dimensions and the
!          number of Epochs to be computed.  For eight epochs and dimensions
!          (91,151,6), it took 76 minutes on a Sun SPARCstation IPX and the
!          resulting file is 15.8 Mbytes.
!
!          Use APXMKA to create the interpolation tables for a single time
!          without writing them.  Otherwise, use APXWRA create the tables
!          and save them in a file.
!
!              CALL APXMKA (MSGUN, EPOCH, GPLAT,GPLON,GPALT,NLAT,NLON,NALT,
!             +            WK,LWK, IST)
!
!              CALL APXWRA (MSGUN, FILNAM,IUN, EPOCH,NEPOCH,
!             +            GPLAT,GPLON,GPALT,NLAT,NLON,NALT, WK,LWK, IST)
!
!          INPUTS:
!            MSGUN  = Fortran unit number to write diagnostic messages.
!            EPOCH  = Time formatted as UT year and fraction; e.g., 1990.0
!                     is 0 UT 1 Jan 1990.  For APXMKA EPOCH is single valued;
!                     for APXWRA EPOCH contains NEPOCH values.
!            GPLAT,GPLON,GPALT = Grid point latitudes, longitudes, and
!                     altitudes.  NLAT values in GPLAT must be in the range
!                     -90. to +90. degrees North.  NLON values in GPLON
!                     must be in the range -180. to 180. degrees East.  NALT
!                     values in GPALT are km MSL.  Each array must be in
!                     ascending numerical order but they may have arbitrary
!                     spacing (as desired for model grid resolution).
!                     Subroutine GGRID can be used to set up these arrays.
!            NLAT,NLON,NALT = Number of latitudes, longitudes and altitudes
!                     are the single dimensions of the arrays GPLAT,GPLON,
!                     and GPALT.  Subroutine GGRID can be used to select
!                     appropriate values.  In general, increasing the number
!                     of grid points, increases the model resolution.  While
!                     array values at the grid points are exact, values at
!                     locations between grid points are estimated by linear
!                     interpolation.  One can increase the grid resolution to
!                     the point of degraded accuracy very close to the poles
!                     of quantities involving east-west gradients, viz.,
!                     B(1), G, H, W, Bhat(1), D1(1), D2(1), D3, E1, E2, E3,
!                     F1(2), and F2(2).  This is evident with dimensions
!                     301x501x15 for a global grid, 0-1000 km.
!            WK    =  Work array is used internally, i.e., there is no need
!                     to access the contents.  WK should not be altered
!                     between initialization (by APXMKA, APXWRA or APXRDA)
!                     and use (by APXGGC, APXALL, APXMALL, or APXQ2G).
!            LWK   =  Dimension of WK >=  NLAT*NLON*NALT*5 + NLAT+NLON+NALT
!
!          Additional INPUTS for APXWRA only:
!            FILNAM = file name where arrays are stored.
!            IUN    = Fortran unit number to be associated with FILNAM.
!            NEPOCH = Number of times in EPOCH.
!
!          RETURNS:
!            IST = Return status: = 0  okay
!                                 > 0  failed
!
!          Declarations for APXMKA formal arguments:
!            DIMENSION GPLAT(NLAT), GPLON(NLON), GPALT(NALT), WK(LWK)
!
!          Additional declarations for APXWRA formal arguments:
!            DIMENSION EPOCH(NEPOCH)
!            CHARACTER*(*) FILNAM
!
!------------------------------------------------------------------------------
!
!       ENTRY APXRDA:
!          Read back tables (previously created by calling APXWRA) in
!          preparation for magnetic coordinate determination (using APXALL,
!          APXMALL, APXQ2G).  If multiple dates were stored, interpolate
!          in time.
!
!              CALL APXRDA (MSGUN, FILNAM,IUN, DATE, WK,LWK, IST)
!
!          INPUTS:
!            MSGUN  = Fortran unit number to write diagnostic messages.
!            FILNAM = file name where arrays are stored.
!            IUN    = Fortran unit number to be associated with FILNAM.
!            DATE   = Time formatted as UT year and fraction; e.g., 1990.0
!                     is 0 UT 1 Jan 1990.
!            WK,LWK = Same as APXMKA
!
!          RETURNS:
!            IST = Return status: = 0  okay
!                                 = -1 non-fatal date problem
!                                 > 0  failed
!
!          Declarations for APXRDA formal arguments:
!            CHARACTER*(*) FILNAM
!            DIMENSION WK(LWK)
!
!------------------------------------------------------------------------------
!
!       ENTRY APXGGC:
!          Obtain grid coordinates from tables read back using APXRDA.
!
!              CALL APXGGC (MSGUN, WK,LWK, GPLAT,GPLON,GPALT,NLAT,NLON,NALT,
!             +            IST)
!
!          INPUTS:
!            MSGUN  = Fortran unit number to write diagnostic messages.
!            WK,LWK = Same as APXMKA
!
!          RETURNS:
!            GPLAT,GPLON,GPALT = Grid point latitudes, longitudes, and
!                                altitudes.  Units are degrees and kilometers.
!            NLAT,NLON,NALT    = Number of latitudes, longitudes and altitudes.
!            IST               = Return status, where IST=0 implies okay.
!
!          Declarations for APXGGC formal arguments:
!            DIMENSION GPLAT(NLAT), GPLON(NLON), GPALT(NALT), WK(LWK)
!
!
!------------------------------------------------------------------------------
!
!
!       ENTRY APXALL:
!          Determine Apex coordinates by interpolation from precalculated
!          arrays.  First call APXMKA, APXWRA, or APXRDA to load look-up
!          tables in memory.
!
!              CALL APXALL (GLAT,GLON,ALT, WK, A,ALAT,ALON, IST)
!
!          INPUTS:
!            GLAT = Geographic (geodetic) latitude, degrees, must be within
!                   the grid domain (GPLAT(1) <= GLAT <= GPLAT(NLAT)).
!            GLON = Geographic (geodetic) longitude, degrees, must be within
!                   one revolution of the grid domain:
!                     GPLON(1) <= GLON-360.,GLON, or GLON+360. <= GPLON(NLON))
!            ALT  = Altitude, km
!            WK   = same as entry APXMKA
!          RETURNS:
!            A    = Apex radius, normalized by Req
!            ALAT = Apex latitude, degrees
!            ALON = Apex longitude, degrees
!            IST  = Return status:  okay (0); or failure (1).
!
!          Dimensions of non-scalar arguments to APXALL:
!            WK(LWK)
!
!
!------------------------------------------------------------------------------
!
!       ENTRY APXMALL:
!          Compute Modified Apex coordinates, quasi-dipole coordinates,
!          base vectors and other parameters by interpolation from
!          precalculated arrays.  First call APXMKA, APXWRA, or APXRDA
!          to load look-up tables in memory.
!
!              CALL APXMALL (GLAT,GLON,ALT,HR, WK,
!             +             B,BHAT,BMAG,SI,
!             +             ALON,
!             +             XLATM,VMP,W,D,BE3,SIM,D1,D2,D3,E1,E2,E3,
!             +             XLATQD,F,F1,F2 , IST)
!
!          Reference:  Richmond, A. D., Ionospheric Electrodynamics Using
!          Magnetic Apex Coordinates, J. Geomag. Geoelectr., 47, 191-212, 1995.
!
!          INPUTS:
!            GLAT = Geographic (geodetic) latitude, degrees, must be within
!                   the grid domain (GPLAT(1) <= GLAT <= GPLAT(NLAT)).
!            GLON = Geographic (geodetic) longitude, degrees, must be within
!                   one revolution of the grid domain:
!                     GPLON(1) <= GLON-360.,GLON, or GLON+360. <= GPLON(NLON))
!            ALT  = Altitude, km
!            HR   = Reference altitude, km (used only for Modified Apex coords)
!            WK   = same as entry APXMKA
!          RETURNS:
!            B      = magnetic field components (east, north, up), in nT
!            BHAT   = components (east, north, up) of unit vector along
!                     geomagnetic field direction
!            BMAG   = magnitude of magnetic field, in nT
!            SI     = sin(I)
!            ALON   = apex longitude = modified apex longitude = quasi-dipole
!                     longitude, degrees
!            XLATM  = Modified Apex latitude, degrees
!            VMP    = magnetic potential, in T.m
!            W      = W of reference above, in km**2 /nT  (i.e., 10**15 m**2 /T)
!            D      = D of reference above
!            BE3    = B_e3 of reference above (= Bmag/D), in nT
!            SIM    = sin(I_m) described in reference above
!            D1,D2,D3,E1,E2,E3 = components (east, north, up) of base vectors
!                     described in reference above
!            XLATQD = quasi-dipole latitude, degrees
!            F      = F described in reference above for quasi-dipole
!                     coordinates
!            F1,F2  = components (east, north) of base vectors described in
!                     reference above for quasi-dipole coordinates
!            IST    = Return status:  okay (0); or failure (1).
!
!          Dimensions of non-scalar arguments to APXMALL:
!            GPLAT(NLAT),GPLON(NLON),GPALT(NALT),WK(LWK),
!            B(3),BHAT(3),D1(3),D2(3),D3(3), E1(3),E2(3),E3(3), F1(2),F2(2)
!
!
!------------------------------------------------------------------------------
!
!       ENTRY APXQ2G:
!          Convert from quasi-dipole to geodetic coordinates, APXQ2G
!          (input magnetic, output geodetic) is the functional inverse
!          of APXALL or APXMALL (input geodetic, output magnetic).  First
!          call APXMKA, APXWRA, or APXRDA to load look-up tables in memory.
!
!              CALL APXQ2G (QDLAT,QDLON,ALT, WK, GDLAT,GDLON, IST)
!
!          INPUTS:
!            QDLAT = quasi-dipole latitude in degrees
!            QDLON = quasi-dipole longitude in degrees
!            ALT   = altitude in km
!            WK    = same as entry APXMKA
!          RETURNS:
!            GDLAT = geodetic latitude in degrees
!            GDLON = geodetic longitude in degrees
!            IST   = Return status: =  0 okay
!                                   = -1 non-fatal, results are not as
!                                        close as PRECISE
!                                   >  0 failure
!
!          Dimensions of non-scalar arguments to APXQ2G:
!            WK(LWK)
!
!
!------------------------------------------------------------------------------
!
!          EXTERNALS:
!            apex.f   - APEX, etc. source code
!            magfld.f - COFRM, DYPOL, FELDG are the DGRF/IGRF.
!            magloctm.f, subsol.f and cossza.f - are not externals, but they
!                       are related, computing magnetic local time, the sub-
!                       solar point, and the cosine of the solar zenith angle.
!            test.f  -  also not an external, rather it is an example driver
!                       program demonstrating various call sequences.  It
!                       includes subroutine GGRID.
!
!          INSTALLATION SPECIFICS:
!            The only (known) machine dependency is parameter IRLF described
!            below.
!            Compiler default is usually to initialize memory to 0.  If this
!            is not true, an error trap involving KGMA may not work.
!
!          HISTORY:
!            Originally coded 940824 by A. D. Richmond, NCAR.  Components of
!            routines in the package came from previous work by Vincent
!            Wickwar (COFRM, FELDG in magfld.f).
!
!            Modified Sep 95 (Roy Barnes):  Changes were made with the
!            objective to allow the user to completely control definition of
!            the interpolation grid.  While doing this the order of the
!            ENTRYs in the first subroutine and arguments lists were
!            changed:  APXMKA, APXWRA, APXRDA (formerly GETARR) and the other
!            ENTRYs now include a work array (WK) which holds arrays X,Y,Z,V,
!            GPLAT,GPLON and GPALT.  Subroutine SETLIM was removed and the
!            grid creation algorithm based on NVERT originally integral to
!            GETARR and SETLIM has been extracted, but still available as
!            an example (see GGRID in file test.f).  Subroutine TSTDIM has a
!            different role, so it is now CKGP (check grid points).
!            MAKEXYZV was also changed to accomodate explicit grid point
!            arrays instead of computed values from an index, minimum and
!            delta.  Only one format is written now, so that is is possible
!            to concatinate files for different epochs.  This required
!            changing delta time logic from fixed 5 yr epochs to computed
!            differences which may vary.  Also, the ties to DGRF/IGRF dates
!            have been removed.
!
!            Modified Sep 96 (Roy Barnes):  Corrected bug in APXQ2G longitude
!            iteration.  The latitude iteration is now constrained to not
!            step beyond the (partial global) interpolation grid.  NITER was
!            increased from 9 to 14 and code to determine cos(lambda') (CLP)
!            was revised by Art Richmond to reduce truncation errors near
!            the geographic poles.
!
!            Modified Sep 97:  Corrected comments, definition of COLAT.
!
!            Modified Dec 98:  Change GLON input to try +/-360 values
!            before rejecting when out of defined grid (GDLON) range;
!            affects INTRP
!
!            Modified Feb-Mar 99:  (1) Corrected a typo (bad variable name)
!            in diagnostic print in INTRP: GLO -> GLON.  This error was
!            probably introduced Dec 98, so no-one had access to the
!            bad copy and, therefore, no announcement is needed.  (2) Also
!            modified APXMALL and APXQ2G:  When very close to a geographic
!            pole, gradients are recalculated after stepping away; in this
!            situation, the latitude input to INTRP was changed from GLAT
!            to GLATX.  This is affects gradients when within 0.1 degrees
!            of a pole (defined as the larger of GLATLIM, 89.9 deg, or the
!            second largest grid point).  (3) Changed MAKEXYZV to make X,Y,Z,
!            V constant for all longitudes (at each alt) when poleward of
!            POLA; POLA was added to /APXCON/.  (4)  Replaced definition of
!            DYLON in APXQ2G.  (5) Reduced NITER to 10 from 14.  Changes 3-5
!            fix a problem where APXQ2G calculations were failing to satisify
!            PRECISE within NITER iterations at the pole.  (6) Replace
!            XYZ2APX with revised trigonometry to avoid a truncation problem
!            at the magnetic equator.  Most changes were devised by Art
!            Richmond and installed by Roy B.
!
!            Questions about this version should be directed to Roy Barnes
!            NCAR (email: bozo@ucar.edu phone: 303/497-1230).
!
!------------------------------------------------------------------------------
!
      use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none
!
!-------------------------------Commons---------------------------------
!
!          Common APXDIPL is assigned in MAKEXYZV but computed in DYPOL
!            COLAT = geocentric colatitude (degrees) of north geomagnetic pole
!            ELON  = geocentric east longitude (degrees) of north geomagnetic
!                    pole
!            VP    = magnetic potential at 1 RE, south geomagnetic pole
!            CTP   = cos(colat*dtor)
!            STP   = sin(colat*dtor)

      real(r8) colat, elon, vp, ctp, stp
      COMMON /APXDIPL/  COLAT,ELON,VP,CTP,STP

!          Common APXCON is assigned here
!            RTOD, DTOR = radians to degrees (180/pi) and inverse
!            RE, REQ    = 6371.2, 6378.160 (Mean and equatorial Earth radius)
!            MSGU       = MSGUN to be passed to subroutines
!            POLA       = Pole angle (deg); when the geographic latitude is
!                         poleward of POLA, X,Y,Z,V are forced to be constant.
!                         for all longitudes at each altitude
!
! bf: separate ints from reals in commons:
      integer msgu
      COMMON /APXCON_int/ MSGU

      real(r8) rtod, dtor, re, req, pola
      COMMON /APXCON/ RTOD,DTOR,RE,REQ,POLA
!
!------------------------------Arguments--------------------------------
!
      integer msgun, nlat, nlon, nalt, lwk, ist
      real(r8) GPLAT(NLAT),GPLON(NLON),GPALT(NALT), EPOCH(*), WK(LWK)
!
!-----------------------------Parameters------------------------------
!
!          XMISS   = value used to fill returned variables when requested
!                    point is outside array boundaries
!          GLATLIM = Extreme polar latitude allowed before changing east-west
!                    gradient calculation to avoid possible underflow at
!                    poles.  GLATLIM is superseded when the second to last
!                    grid point value is closer to a pole.
!          PRECISE = (angular distance)**2 (radians**2) of precision of
!                    transform from quasi-dipole to geodetic coordinates.
!                    7.6E-11 corresponds to an accuracy of 0.0005 degree.
!          IRLF    = Record length factor required to convert the computed
!                    length of direct access records from words to those
!                    units appropriate to this computer.  IRLF is 1 when
!                    RECL units are words, otherwise it is a function of
!                    the computer word length; e.g.,
!                      IRLF = 1 on a DEC  (32-bit words, RECL units are words)
!                      IRLF = 4 on a PC   (32-bit words, RECL units are bytes)
!                      IRLF = 4 on a Sun  (32-bit words, RECL units are bytes)
!                      IRLF = 8 on a Cray (64-bit words, RECL units are bytes)
!          DATDMX  = maximum time difference (years) allowed between the
!                    requested date and the single epoch in arrays.
!          DATIMX  = maximum time difference (years) allowed between the
!                    requested date and the closest epoch in the stored
!                    arrays (apropos multiple stored dates).
      integer irlf
      real(r8) xmiss, glatlim, precise, datdmx, datimx
      PARAMETER (XMISS=-32767._r8 , GLATLIM=89.9_r8 , PRECISE=7.6E-11_r8, &
                DATDMX=1._r8 , DATIMX=2.5_r8 , IRLF=8 )
!
!---------------------------Local variables-----------------------------
!
      integer jst
      integer nloi, nlai, nali, i1, le, lc, n, it, iter, niter, iun, nlar
      integer nepoch, nla, nal, nlo, nalr, nlor, lbx, ldr, ngp
      integer lfn, nepok, lcn, i, j, i2, il, lbt, ngm1, k, leg
      integer lal, llo, lbz, lby, lla, lbv

      real(r8) sth, cth, fxdum, glatx, dfvdh, dfxdh, dfvdln, dfzdh, dfydh
      real(r8) dmvdth, dmzdth, dmydth, dmydh, dmxdh, fzdum, fydum
      real(r8) dfzdth, dfydth, dfxdln, dfvdth, fz, fy, dfxdth, fv, dmvdh
      real(r8) dmzdh, r3, clm, dfzdln, dfydln, dmxdth, fvdum, fx, ti, to
      real(r8) year1, datd, frac, glon, omf, cola2, year2, vp2, elon2, x0
      real(r8) slp, sad, slm2, clm2, coslm, sal, cad, slm, dylon, cal2
      real(r8) sad2, clp, clp2, cal, xnorm, distlon, z0, y0, ylon, ylat, hgrd2n
      real(r8) hgrd2e, angdist, hgrd2, ydif, xdif, dist2, zdif, t, date, sim, be3
      real(r8) f, xlatqd, glalmx, glalmn, d, glat, gdlon, a, glo
      real(r8) qdlon, qdlat, gdlat, alt, xlatm, si, w, vmp, alon, alat
      real(r8) bmag, hr, year, tb, tl, r3_2

      real(r8) GRADX(3), GRADY(3), GRADZ(3), GRADV(3), &
                GRCLM(3), CLMGRP(3), RGRLP(3), &
                B(3),BHAT(3), D1(3),D2(3),D3(3), E1(3),E2(3),E3(3), &
                F1(2),F2(2)

      CHARACTER*(*) FILNAM

      CHARACTER CALNM*7, CMP*10, EDG*5

      integer kgma
      DATA KGMA /0/

      SAVE KGMA, GLALMN,GLALMX, NLA,NLO,NAL, LBX,LBY,LBZ,LBV,LLA,LLO,LAL
!
!-----------------------------------------------------------------------
!
      CALNM = 'APXMKA'
      LCN   = 6
      KGMA  = 1
      MSGU  = MSGUN
      NEPOK = 1
      IF (NLAT .LT. 2 .OR. NLON .LT. 2 .OR. NALT .LT. 2) GO TO 9100
      NLA = NLAT
      NLO = NLON
      NAL = NALT
      GO TO 40

!      ENTRY APXGGC (MSGUN,WK,LWK, GPLAT,GPLON,GPALT,NLAR,NLOR,NALR,IST)
!          Sep 95 R. Barnes
      CALNM = 'APXGGC'
      LCN   = 6
      MSGU  = MSGUN
      IF (KGMA .LT. 1) GO TO 9300
      J = LLA
      DO 10 I=1,NLA
      GPLAT(I) = WK(J)
   10 J = J + 1
      DO 20 I=1,NLO
      GPLON(I) = WK(J)
   20 J = J + 1
      DO 30 I=1,NAL
      GPALT(I) = WK(J)
   30 J = J + 1
      NLAR = NLA
      NLOR = NLO
      NALR = NAL
      IST = 0
      RETURN

      ENTRY APXWRA (MSGUN, FILNAM,IUN, EPOCH,NEPOCH, &
                   GPLAT,GPLON,GPALT,NLAT,NLON,NALT, WK,LWK, IST)
!          Sep 95 R. Barnes
      CALNM = 'APXWRA'
      LCN   = 6
      KGMA  = 2
      MSGU  = MSGUN
      NEPOK = NEPOCH
      IF (NLAT .LT. 2 .OR. NLON .LT. 2 .OR. NALT .LT. 2) GO TO 9100
      NLA = NLAT
      NLO = NLON
      NAL = NALT
      GO TO 40

      ENTRY APXRDA (MSGUN, FILNAM,IUN, DATE, WK,LWK, IST)
!          Sep 95 R. Barnes
      CALNM = 'APXRDA'
      LCN   = 6
      KGMA  = 3

!          Open the read-back file with a temporary record length, get
!          the grid dimensions from the first values, then close it, (so
!          it can be reopened later with the proper LDR):
      LDR = 7*IRLF
      OPEN (IUN,FILE=FILNAM,ACCESS='direct',RECL=LDR,STATUS='old', &
                                                             IOSTAT=IST)
      MSGU = MSGUN
      IF (IST .NE. 0) GO TO 9110
      READ  (IUN,REC=1,IOSTAT=IST) YEAR,COLAT,ELON,VP,NLA,NLO,NAL
      IF (IST .NE. 0) GO TO 9120
      CLOSE (IUN)

   40 RE   = 6371.2_r8
      REQ  = 6378.160_r8
      RTOD = 45._r8/ATAN(1._r8)
      DTOR = 1._r8/RTOD
      POLA = 90._r8 - SQRT (PRECISE) * RTOD
      LFN = 0
      IF (KGMA .EQ. 1) GO TO 51
      DO 50 I=1,LEN(FILNAM)
      IF (FILNAM(I:I) .EQ. ' ') GO TO 51
   50 LFN = LFN + 1
   51 CONTINUE

!          Save grid dimensions, establish direct access rec length, and
!          determine indices into the work array.  WK contains arrays
!          X,Y,Z,V,temp,GPLAT,GPLON,GPALT where X thru tmp are dimensioned
!          (NLAT,NLON,NALT); tmp is scratch space used during read back.
      NGP = NLA*NLO*NAL
      NGM1= NGP - 1
      LDR = NGP * IRLF
      LBX = 1
      LBY = LBX + NGP
      LBZ = LBY + NGP
      LBV = LBZ + NGP
      LBT = LBV + NGP
      LLA = LBT + NGP
      LLO = LLA + NLA
      LAL = LLO + NLO
      LEG = LAL + NAL-1
      IF (LWK .LT. LEG) GO TO 9130

      IF (KGMA .EQ. 3) GO TO 200

!          Make and optionally write interpolation arrays for NEPOK times
      IF (KGMA .EQ. 2) THEN
	OPEN (IUN,FILE=FILNAM,ACCESS='direct',RECL=LDR,STATUS='new', &
             IOSTAT=IST)
	IF (IST .NE. 0) GO TO 9115
      ENDIF

      CALL CKGP (CALNM(:LCN),MSGUN,NLAT,NLON,NALT,GPLAT,GPLON,GPALT,IST)
      IF (IST .NE. 0) RETURN
      I = LLA - 1
      DO 60 J=1,NLAT
      I = I + 1
   60 WK(I) = GPLAT(J)
      DO 70 J=1,NLON
      I = I + 1
   70 WK(I) = GPLON(J)
      DO 80 J=1,NALT
      I = I + 1
   80 WK(I) = GPALT(J)

      IF (NEPOK .LT. 1) GO TO 9140
      J = 1
      DO 100 I=1,NEPOK
      CALL MAKEXYZV (EPOCH(I),NLAT,NLON,NALT,GPLAT,GPLON,GPALT, &
                     WK(LBX),WK(LBY),WK(LBZ),WK(LBV))
      IF (KGMA .EQ. 1) GO TO 100
      WRITE (IUN,REC=J) EPOCH(I),COLAT,ELON,VP,NLAT,NLON,NALT
      WRITE (IUN,REC=J+1) (WK(K),K=LLA,LEG)
      WRITE (IUN,REC=J+2) (WK(K),K=LBX,LBX+NGM1)
      WRITE (IUN,REC=J+3) (WK(K),K=LBY,LBY+NGM1)
      WRITE (IUN,REC=J+4) (WK(K),K=LBZ,LBZ+NGM1)
      WRITE (IUN,REC=J+5) (WK(K),K=LBV,LBV+NGM1)
  100 J = J + 6
      IF (KGMA .EQ. 2) CLOSE (IUN)
      IST = 0
      GO TO 300

!          Read back interpolation arrays.  When arrays for multiple times
!          are available, interpolate using the pair bounding the desired
!          time (DATE).  Make an initial pass only to identify closest
!          available times and issue any diagnostics, then the second pass
!          to read the stored arrays (GPLAT,GPLON,GPALT,X,Y,Z,V) and do
!          the linear interpolation/extrapolation.
  200 OPEN (IUN,FILE=FILNAM,ACCESS='direct',RECL=LDR,STATUS='old', &
            IOSTAT=IST)
      IF (IST .NE. 0) GO TO 9110

      READ (IUN,REC=1,IOSTAT=IST) TB
      IF (IST .NE. 0) GO TO 9120
      I2 = 1
      TL = TB
      IL = 1
      I = 1
  210 I = I + 6
      READ (IUN,REC=I,IOSTAT=JST) T
!         JST .NE. 0 is assumed to mean read beyond last record
      IF (JST .NE. 0) GO TO 220
      TO = TL
      TL = T
      IL = I
      IF (DATE .GT. TL) GO TO 210

  220 I1 = IL - 6

      IST = 0
      IF (TL .EQ. TB) THEN
	DATD = ABS (DATE-TB)
	IF (DATD .GT. DATDMX) THEN
	  WRITE (MSGU,9150) CALNM(1:LCN),DATE,DATD,TB
	  IF (TB .EQ. 0._r8) WRITE (MSGU,9155) FILNAM(1:LFN)
	  IST = -1
	ENDIF
	I1 = 1
	I2 = 0
      ELSE IF (DATE .LT. TB) THEN
	WRITE (MSGU,9160) CALNM(1:LCN),DATE,TB,FILNAM(1:LFN)
	IST = -1
      ELSE IF (DATE .GT. TL) THEN
	WRITE (MSGU,9170) CALNM(1:LCN),DATE,TL,FILNAM(1:LFN)
	IST = -1
      ELSE
	DATD = min (DATE-TO,TL-DATE)
	IF (DATD .GT. DATIMX) THEN
	  WRITE (MSGU,9180) CALNM(1:LCN),DATE,TB,TL,FILNAM(1:LFN),DATD
	  IST = -1
	ENDIF
      ENDIF

      READ (IUN,REC=I1) YEAR1,COLAT,ELON,VP,NLAI,NLOI,NALI
      TI = YEAR1
      IF (NLAI.NE.NLA .OR. NLOI.NE.NLO .OR. NALI.NE.NAL) GO TO 9190
      READ (IUN,REC=I1+1) (WK(I),I=LLA,LEG)
      READ (IUN,REC=I1+2) (WK(I),I=LBX,LBX+NGM1)
      READ (IUN,REC=I1+3) (WK(I),I=LBY,LBY+NGM1)
      READ (IUN,REC=I1+4) (WK(I),I=LBZ,LBZ+NGM1)
      READ (IUN,REC=I1+5) (WK(I),I=LBV,LBV+NGM1)
      IF (I2 .EQ. 1) THEN
	READ (IUN,REC=I1+6) YEAR2,COLA2,ELON2,VP2,NLAI,NLOI,NALI
	TI = YEAR2
	IF (NLAI.NE.NLA .OR. NLOI.NE.NLO .OR. NALI.NE.NAL) GO TO 9190
	LE = LBT + NLA+NLO+NAL - 1
	READ (IUN,REC=I1+7) (WK(I),I=LBT,LE)
	J = LLA
	DO 230 I=LBT,LE
	IF (WK(J) .NE. WK(I)) GO TO 9200
  230   J = J + 1
	FRAC = (DATE-YEAR1) / (YEAR2-YEAR1)
	OMF  = 1._r8 - FRAC
	LE = LBT + NGM1
	READ (IUN,REC=I1+8) (WK(I),I=LBT,LE)
	J = LBX
	DO 240 I=LBT,LE
	WK(J) =  OMF*WK(J) + FRAC*WK(I)
  240   J = J + 1
	READ (IUN,REC=I1+9) (WK(I),I=LBT,LE)
	DO 250 I=LBT,LE
	WK(J) =  OMF*WK(J) + FRAC*WK(I)
  250   J = J + 1
	READ (IUN,REC=I1+10) (WK(I),I=LBT,LE)
	DO 260 I=LBT,LE
	WK(J) =  OMF*WK(J) + FRAC*WK(I)
  260   J = J + 1
	READ (IUN,REC=I1+11) (WK(I),I=LBT,LE)
	DO 270 I=LBT,LE
	WK(J) =  OMF*WK(J) + FRAC*WK(I)
  270   J = J + 1
	COLAT = OMF*COLAT + FRAC*COLA2
	ELON  = OMF*ELON  + FRAC*ELON2
	VP    = OMF*VP    + FRAC*VP2
      ENDIF

      CTP  = COS (COLAT*DTOR)
      STP  = SIN (COLAT*DTOR)
      YEAR = DATE
      CLOSE (IUN)

!          Establish for this grid polar latitude limits beyond which east-west
!          gradients are computed differently to avoid potential underflow
  300 GLALMX = max ( GLATLIM,WK(LLA+NLA-2))
      GLALMN = min (-GLATLIM,WK(LLA+1))

      RETURN

!*******************************************************************************

      ENTRY APXMALL (GLAT,GLO ,ALT,HR, WK, LWK,                      &  !Inputs
                    B,BHAT,BMAG,SI,                                  &  !Mag Fld
                    ALON,                                            &  !Apx Lon
                    XLATM,VMP,W,D,BE3,SIM,D1,D2,D3,E1,E2,E3,         &  !Mod Apx
                    XLATQD,F,F1,F2 , IST)                               !Qsi-Dpl
!          940822 A. D. Richmond, Sep 95 R. Barnes

!          Test to see if WK has been initialized
      CALNM = 'APXMALL'
      LCN   = 7
      IF (KGMA .LT. 1) GO TO 9300

!          Alias geodetic longitude input in case INTRP adjusts it by +/-360
      GLON = GLO

      CALL INTRP (GLAT,GLON,ALT, WK(LBX),WK(LBY),WK(LBZ),WK(LBV), &
                 NLA,NLO,NAL,    WK(LLA),WK(LLO),WK(LAL), &
                 FX,FY,FZ,FV, &
                 DFXDTH,DFYDTH,DFZDTH,DFVDTH,DFXDLN,DFYDLN,DFZDLN, &
                 DFVDLN,DFXDH,DFYDH,DFZDH,DFVDH, CALNM(1:LCN),IST)

      IF (IST .NE. 0) THEN
	CALL SETMISS (XMISS, XLATM,ALON,VMP,B,BMAG,BE3,SIM,SI,F,D,W, &
                     BHAT,D1,D2,D3,E1,E2,E3,F1,F2)
	RETURN
      ENDIF

      CALL ADPL (GLAT,GLON,CTH,STH,FX,FY,FZ,FV, &
                DFXDTH,DFYDTH,DFZDTH,DFVDTH,DFXDLN,DFYDLN,DFZDLN,DFVDLN)
      CALL GRADXYZV (ALT,CTH,STH, &
                DFXDTH,DFYDTH,DFZDTH,DFVDTH,DFXDLN,DFYDLN,DFZDLN,DFVDLN, &
                DFXDH,DFYDH,DFZDH,DFVDH,GRADX,GRADY,GRADZ,GRADV)

      IF (GLAT .GT. GLALMX .OR. GLAT .LT. GLALMN) THEN
!          If the point is very close to either the North or South
!          geographic pole, recompute the east-west gradients after
!          stepping a small distance from the pole.
	GLATX = GLALMX
	IF (GLAT .LT. 0._r8) GLATX = GLALMN
!990225 CALL INTRP (GLAT ,GLON,ALT, WK(LBX),WK(LBY),WK(LBZ),WK(LBV),
	CALL INTRP (GLATX,GLON,ALT, WK(LBX),WK(LBY),WK(LBZ),WK(LBV), &
                   NLA,NLO,NAL,    WK(LLA),WK(LLO),WK(LAL), &
                   FXDUM,FYDUM,FZDUM,FVDUM, &
                   DMXDTH,DMYDTH,DMZDTH,DMVDTH,DFXDLN,DFYDLN,DFZDLN, &
                   DFVDLN,DMXDH,DMYDH,DMZDH,DMVDH, CALNM(1:LCN),IST)
	CALL ADPL (GLATX,GLON,CTH,STH,FXDUM,FYDUM,FZDUM,FVDUM, DMXDTH, &
                  DMYDTH,DMZDTH,DMVDTH,DFXDLN,DFYDLN,DFZDLN,DFVDLN)
	CALL GRAPXYZV (ALT,CTH,STH, DFXDLN, &
                      DFYDLN,DFZDLN,DFVDLN,GRADX,GRADY,GRADZ,GRADV)
      ENDIF

      CALL GRADLPV (HR,ALT,FX,FY,FZ,FV,GRADX,GRADY,GRADZ,GRADV, &
                   XLATM,ALON,VMP,GRCLM,CLMGRP,XLATQD,RGRLP,B,CLM,R3_2)
      CALL BASVEC (HR,XLATM,GRCLM,CLMGRP,RGRLP,B,CLM,R3_2, &
                  BMAG,SIM,SI,F,D,W,BHAT,D1,D2,D3,E1,E2,E3,F1,F2)
      BE3 = BMAG/D

      IST = 0
      RETURN


!*******************************************************************************

!      ENTRY APXALL (GLAT,GLO ,ALT, WK, A,ALAT,ALON, IST)
!          940802 A. D. Richmond, Sep 95 R. Barnes

!          Test to see if WK has been initialized
      CALNM = 'APXALL'
      LCN   = 6
      IF (KGMA .LT. 1) GO TO 9300

!          Alias geodetic longitude input in case INTRPSC adjusts it by +/-360
      GLON = GLO

      CALL INTRPSC (GLAT,GLON,ALT, WK(LBX),WK(LBY),WK(LBZ), &
                    NLA,NLO,NAL,   WK(LLA),WK(LLO),WK(LAL), &
                    FX,FY,FZ, CALNM(1:LCN), IST)
      IF (IST .NE. 0) GO TO 600

      CALL ADPLSC (GLAT,GLON,FX,FY,FZ)

      CALL XYZ2APX (ALT,FX,FY,FZ,A,ALAT,ALON,IST)
      IF (IST .EQ. 0) GO TO 601

  600 A    = XMISS
      ALAT = XMISS
      ALON = XMISS
  601 CONTINUE

      RETURN

!*******************************************************************************

!      ENTRY APXQ2G (QDLAT,QDLON,ALT, WK, GDLAT,GDLON, IST)
!          940819 A. D. Richmond, Sep 95 R. Barnes, Sep 96 mod A. D. Richmond
!          Input guessed geodetic coordinates (YLAT,YLON) to INTRP and
!          compare the returned magnetic coordinates to those desired.
!          If the guess is not sufficiently close (PRECISE), make another
!          guess by moving in the direction of the gradient of a quantity
!          (DIST2) that approximates the squared angular distance between
!          the returned and desired magnetic coordinates.

!          Test to see if WK has been initialized
      CALNM = 'APXQ2G'
      LCN   = 6
      IF (KGMA .LT. 1) GO TO 9300

!          Determine quasi-cartesian coordinates on a unit sphere of the
!          desired magnetic lat,lon in quasi-dipole coordinates.
      X0 = COS (QDLAT*DTOR) * COS (QDLON*DTOR)
      Y0 = COS (QDLAT*DTOR) * SIN (QDLON*DTOR)
      Z0 = SIN (QDLAT*DTOR)

!          Initial guess:  use centered dipole, convert to geocentric coords
      CALL GM2GC (QDLAT,QDLON,YLAT,YLON)

!          Iterate until (angular distance)**2 (units: radians) is within
!          PRECISE of location (QDLAT,QDLON) on a unit sphere.
!
! 4/00: 10 iters not enough for tgcm14:
!     NITER = 10
      NITER = 20
      DO 1000 ITER=1,NITER
      CALL INTRP (YLAT,YLON,ALT, WK(LBX),WK(LBY),WK(LBZ),WK(LBV), &
                 NLA,NLO,NAL,    WK(LLA),WK(LLO),WK(LAL), &
                 FX,FY,FZ,FV, &
                 DFXDTH,DFYDTH,DFZDTH,DFVDTH,DFXDLN,DFYDLN,DFZDLN, &
                 DFVDLN,DFXDH,DFYDH,DFZDH,DFVDH, CALNM(1:LCN),IST)
      IF (IST .NE. 0) GO TO 9400
      CALL ADPL (YLAT,YLON,CTH,STH,FX,FY,FZ,FV, &
                DFXDTH,DFYDTH,DFZDTH,DFVDTH,DFXDLN,DFYDLN,DFZDLN,DFVDLN)

      DISTLON = COS(YLAT*DTOR)
      IF (YLAT .GT. GLALMX .OR. YLAT .LT. GLALMN) THEN
	GLATX = GLALMX
	IF (YLAT.LT.0._r8) GLATX = GLALMN
	DISTLON = COS (GLATX*DTOR)
!990225 CALL INTRP (YLAT ,YLON,ALT, WK(LBX),WK(LBY),WK(LBZ),WK(LBV),
	CALL INTRP (GLATX,YLON,ALT, WK(LBX),WK(LBY),WK(LBZ),WK(LBV), &
                   NLA,NLO,NAL,     WK(LLA),WK(LLO),WK(LAL), &
                   FXDUM,FYDUM,FZDUM,FVDUM, &
                   DMXDTH,DMYDTH,DMZDTH,DMVDTH,DFXDLN,DFYDLN,DFZDLN, &
                   DFVDLN,DMXDH,DMYDH,DMZDH,DMVDH, CALNM(1:LCN),IST)
	CALL ADPL (GLATX,YLON,CTH,STH,FXDUM,FYDUM,FZDUM,FVDUM, &
                DMXDTH,DMYDTH,DMZDTH,DMVDTH,DFXDLN,DFYDLN,DFZDLN,DFVDLN)
      ENDIF

!          At this point, FX,FY,FZ are approximate quasi-cartesian
!          coordinates on a unit sphere for the quasi-dipole coordinates
!          corresponding to the geodetic coordinates YLAT, YLON.
!          Normalize the vector length of (FX,FY,FZ) to unity using XNORM
!          so that the resultant vector can be directly compared with the
!          target vector (X0,Y0,Z0).
      XNORM = SQRT(FX*FX + FY*FY + FZ*FZ)
      XDIF = FX/XNORM - X0
      YDIF = FY/XNORM - Y0
      ZDIF = FZ/XNORM - Z0
!          DIST2 = square of distance between normalized (FX,FY,FZ) and
!          X0,Y0,Z0.
      DIST2 = XDIF*XDIF + YDIF*YDIF + ZDIF*ZDIF

      IF (DIST2 .LE. PRECISE) then
!       if (iter > 10)
!    |    write(iulog,"('Note apxparm: took > 10 iterations: iter=',i3)")
!    |      iter
        GO TO 1001
      endif
!          HGRD2* = one-half of east or north gradient of DIST2 on unit sphere.
      HGRD2E =  (XDIF*DFXDLN + YDIF*DFYDLN + ZDIF*DFZDLN)/DISTLON
      HGRD2N = -(XDIF*DFXDTH + YDIF*DFYDTH + ZDIF*DFZDTH)
      HGRD2  = SQRT(HGRD2E*HGRD2E + HGRD2N*HGRD2N)
!          ANGDIST = magnitude of angular distance to be moved for new guess
!          of YLAT, YLON.
      ANGDIST = DIST2/HGRD2

!          Following spherical trigonometry moves YLAT, YLON to new location,
!          in direction of grad(DIST2), by amount ANGDIST.
      CAL = -HGRD2N/HGRD2
      SAL = -HGRD2E/HGRD2
      COSLM = COS(YLAT*DTOR)
      SLM = SIN(YLAT*DTOR)
      CAD = COS(ANGDIST)
      SAD = SIN(ANGDIST)
      SLP = SLM*CAD + COSLM*SAD*CAL

!          Old code (below) introduced truncation errors near geographic poles
!     SLP = AMIN1 (SLP,SGLAN) ; sglan = sin(wk(lla+nla-1))*dtor)
!     SLP = AMAX1 (SLP,SGLAS) ; sglas = sin(wk(lla))      *dtor)
!     CLP = SQRT(1. - SLP*SLP)
!     YLAT = ASIN(SLP)*RTOD

      CLM2 = COSLM*COSLM
      SLM2 = SLM*SLM
      SAD2 = SAD*SAD
      CAL2 = CAL*CAL
      CLP2 = CLM2 + SLM2*SAD2 - 2._r8*SLM*CAD*COSLM*SAD*CAL -CLM2*SAD2*CAL2
      CLP = SQRT (max(0._r8,CLP2))
      YLAT = ATAN2(SLP,CLP)*RTOD

!          Restrict latitude iterations to stay within the interpolation grid
!          limits, but let INTRP find any longitude exceedence.  This is only
!          an issue when the interpolation grid does not cover the entire
!          magnetic pole region.
      YLAT = min(YLAT,WK(LLA+NLA-1))
      YLAT = max(YLAT,WK(LLA))

      DYLON = ATAN2 (SAD*SAL,CAD*COSLM-SAD*SLM*CAL)*RTOD
!          Old formula (replaced Mar 99) had truncation problem near poles
!     DYLON = ATAN2 (CLP*SAD*SAL,CAD-SLM*SLP)*RTOD

      YLON  = YLON + DYLON
      IF (YLON .GT.  WK(LLO+NLO-1)) YLON = YLON - 360._r8
      IF (YLON .LT.  WK(LLO)      ) YLON = YLON + 360._r8
 1000 CONTINUE

      WRITE (MSGU,9145) NITER, SQRT(DIST2  )*RTOD , SQRT(PRECISE)*RTOD
       
      EDG = ' '
      IF (YLAT .EQ. WK(LLA+NLA-1)) EDG = 'north'
      IF (YLAT .EQ. WK(LLA))       EDG = 'south'
      IF (EDG .NE. ' ') WRITE (MSGU,9146) EDG

      IST = -1
      GO TO 1010

 1001 IST = 0

 1010 GDLAT = YLAT
      GDLON = YLON

      RETURN

!*******************************************************************************

!          Error Trap diagnostics
 9100 WRITE (MSGU,'(A,'':  NLAT,NLON or NALT < 2 '',3I8)') &
                   CALNM(1:LCN),  NLAT,NLON,NALT
      IST = 1
      RETURN
 9110 WRITE (MSGU,'(A,'': Trouble opening old file "'',A,''"'')') &
       CALNM(1:LCN), FILNAM(1:LFN)
      RETURN
 9115 WRITE (MSGU,'(A,'': Trouble opening new file "'',A,''"'')') &
       CALNM(1:LCN), FILNAM(1:LFN)
      RETURN
 9120 WRITE (MSGU,'(A,'': Trouble reading first record of '',A)') &
       CALNM(1:LCN), FILNAM(1:LFN)
      RETURN
 9130 WRITE (MSGU,&
      '(A,'': LWK is too small; LWK must be at least'',I5,''but LWK ='',I5)') &
      CALNM(1:LCN), LEG, LWK
      IST = 1
      RETURN
 9140 WRITE (MSGU,'(A,'':  NEPOCH must be positive; NEPOCH ='',I8)') &
             CALNM(1:LCN), NEPOK
      IST = 1
      RETURN

 9145 FORMAT('APXQ2G: Warning',I3, &
      ' iterations only reduced the angular difference to',/,8X,F8.5,&
      ' degrees (',F8.5,' degrees is the test criterion)')
 9146 FORMAT( &
      "        Coordinates are on the ",A," edge of the interpolation grid and",/,&
      "        latitude is constrained to stay within grid limits when iterating.") 
 9150 FORMAT (A,': DATE (',F7.2,') differs by',F5.2,' years from the stored EPOCH (',F7.2,')')
 9155 FORMAT ('        A stored date = 0. implies "',A,'" is incorrectly formatted')
 9160 FORMAT (A,': DATE (',F7.2,') is extrapolated before first EPOCH (',F7.2,') in "',A,'"')
 9170 FORMAT (A,': DATE (',F7.2,') is extrapolated after last EPOCH (',F7.2,') in "',A,'"')
 9180 FORMAT (A,': DATE (',F7.2,') minimum difference from the nearest stored ',/,&
      'EPOCHs (',F7.2,', ',F7.2,') in "',A,'"',/,'is',F6.2,' years')

 9190 WRITE (MSGU,9195) CALNM(1:LCN), FILNAM(1:LFN), TI , &
       NLAI,NLOI,NALI,  NLA,NLO,NAL, TB
 9195 FORMAT(A,': Dimensions (lat,lon,alt) read from "',A,&
      '" for EPOCH ',F7.2,/,'        (',I5,',',I5,',',I5,&
      ') do not match (',I5,',',I5,',',I5,') from the',/,&
      '        first EPOCH read (',F7.2,')')
      IST = 1
      RETURN
 9200 CMP = 'latitudes'
      LC  = 9
      I1 = LLA - 1
      IT = LBT - 1
      N  = NLA
      IF (I .LT. LLO) GO TO 9201
      CMP = 'longitudes'
      LC   = 10
      I1 = I1 + NLA
      IT = IT + NLA
      N  = NLO
      IF (I .LT. LAL) GO TO 9201
      CMP = 'altitudes'
      LC  = 9
      I1 = I1 + NLO
      IT = IT + NLO
      N  = NAL
 9201 WRITE (MSGU,9205)  CALNM(1:LCN), &
        CMP(1:LC), FILNAM(1:LFN), TI,TB, (WK(I1+I),WK(IT+I),I=1,N)
 9205 FORMAT(A,': Grid ',A,' read from "',A,'" for EPOCH',F8.2,&
      ' do not match the',/,'        first EPOCH (',F7.2,')',/,&
      '        First    Current',/,(4X,2F10.3))
      IST = 1
      RETURN

 9300 WRITE(MSGU,'(A,'': Must first load lookup tables by calling APXMKA, APXWRA or APXRDA'')') &
      CALNM(1:LCN)
      IST = 1
      RETURN
 9400 WRITE (MSGU,'(''APXQ2G: INTRP failed (maybe coordinates are not within interpolation grid)'')')
      IST = 1
      RETURN

      END subroutine apxmka

!================================================================================================

      SUBROUTINE INTRP (GLAT,GLON,ALT, X,Y,Z,V, NLAT,NLON,NALT, &
                       GPLAT,GPLON,GPALT, FX,FY,FZ,FV, &
                       DFXDTH,DFYDTH,DFZDTH,DFVDTH,DFXDLN,DFYDLN,DFZDLN, &
                       DFVDLN,DFXDH,DFYDH,DFZDH,DFVDH, CALNM,IST)
!
!-----------------------------------------------------------------------
!          Interpolation of x,y,z,v and their derivatives
!          940806 A. D. Richmond
!          INPUTS:
!            GLAT    = latitude, degrees
!            GLON    = longitude, degrees
!            ALT     = altitude, km
!            X,Y,Z,V = gridded arrays
!            NLAT,NLON,NALT = three dimensions of x,y,z,v and respective
!                      dimensions of GP___ arrays
!            GPLAT,GPLON,GPALT = grid point geographic locations
!            CALNM   = Name of calling routine (for error diagnostics)
!          OUTPUT:
!            FX = interpolated value of x
!            FY = interpolated value of y
!            FZ = interpolated value of z
!            FV = interpolated value of v
!            DFXDTH,DFYDTH,DFZDTH,DFVDTH = derivatives of x,y,z,v with
!                  respect to colatitude, in radians-1
!            DFXDLN,DFYDLN,DFZDLN,DFVDLN = derivatives of x,y,z,v with
!                  respect to longitude, in radians-1
!            DFXDH,DFYDH,DFZDH,DFVDH = derivatives of x,y,z,v with
!                  respect to altitude, in km-1
!            IST = Status =  0 = okay.
!-----------------------------------------------------------------------
!
      use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none
!
!-------------------------------Commons---------------------------------
!
! bf: separate ints from reals in commons:
      integer msgu
      COMMON /APXCON_int/ MSGU

      real(r8) rtod, dtor, re, req, pola
      COMMON /APXCON/ RTOD,DTOR,RE,REQ,POLA
!
!------------------------------Arguments--------------------------------
!
      integer nlat, nlon, nalt, ist

      real(r8) glat, glon, alt, fx, fy, fz, fv
      real(r8) dfxdth, dfydth, dfzdth, dfvdth, dfxdln, dfydln, dfzdln, dfvdln
      real(r8) dfxdh, dfydh, dfzdh, dfvdh
      real(r8) X(NLAT,NLON,NALT), Y(NLAT,NLON,NALT), Z(NLAT,NLON,NALT), &
               V(NLAT,NLON,NALT) , GPLAT(NLAT),GPLON(NLON),GPALT(NALT)

      CHARACTER*(*) CALNM
!
!---------------------------Local variables-----------------------------
!
      integer ientry, i, j, k

      real(r8) dfzdd, dfvdn, dfvde, dfzde, dfyde, dfydd, dfzdn
      real(r8) dmdfdn, dmdfde, dmdfdd, dmf, dfvdd, fac, omfac
      real(r8) dfydn, dlon, yj, xi, dlat, dfxdn, dfxde, dfxdd
      real(r8) zk, hti, diht

      integer io, jo, ko
      DATA IO, JO, KO / 1, 1, 1 /

      SAVE IO, JO, KO
!
!-----------------------------------------------------------------------
!
      IENTRY = 0
      GO TO 5

!*******************************************************************************
      ENTRY INTRPSC (GLAT,GLON,ALT, X,Y,Z, NLAT,NLON,NALT, &
                     GPLAT,GPLON,GPALT, FX,FY,FZ , CALNM,IST)
!          Interpolation of x,y,z
!          940803 A. D. Richmond.
!          Inputs and outputs:  for definitions, see above.

      IENTRY = 1
    5 IST = 0

      IF (GLAT .LT. GPLAT(1) .OR. GLAT .GT. GPLAT(NLAT)) GO TO 9100
      IF (ALT  .LT. GPALT(1) .OR. ALT  .GT. GPALT(NALT)) GO TO 9300
!          Accept input longitude range +/- one revolution (Dec 98)
      IF (GLON .LT. GPLON(1)   ) GLON = GLON + 360._r8
      IF (GLON .GT. GPLON(NLON)) GLON = GLON - 360._r8
      IF (GLON .LT. GPLON(1) .OR. GLON .GT. GPLON(NLON)) GO TO 9200

      I = IO
      IF (GLAT .LE. GPLAT(I)) GO TO 15
   12 I = I + 1
      IF (GPLAT(I) .LT. GLAT) GO TO 12
      I = I - 1
      GO TO 16
   14 I = I - 1
   15 IF (GPLAT(I) .GT. GLAT) GO TO 14
   16 IO = I
      DLAT = GPLAT(I+1) - GPLAT(I)
      XI   = (GLAT - GPLAT(I)) / DLAT

      J = JO
      IF (GLON .LE. GPLON(J)) GO TO 25
   22 J = J + 1
      IF (GPLON(J) .LT. GLON) GO TO 22
      J = J - 1
      GO TO 26
   24 J = J - 1
   25 IF (GPLON(J) .GT. GLON) GO TO 24
   26 JO = J
      DLON = GPLON(J+1) - GPLON(J)
      YJ   = (GLON - GPLON(J)) / DLON

      K = KO
      IF (ALT .LE. GPALT(K)) GO TO 35
   32 K = K + 1
      IF (GPALT(K) .LT. ALT) GO TO 32
      K = K - 1
      GO TO 36
   34 K = K - 1
   35 IF (GPALT(K) .GT. ALT) GO TO 34
   36 KO = K
      HTI  = RE/(RE+ALT)
      DIHT = RE/(RE+GPALT(K+1)) - RE/(RE+GPALT(K))
      ZK  = (HTI - RE/(RE+GPALT(K))) / DIHT

!          For intrpsc:
      IF (IENTRY.EQ.1) THEN
	CALL TRILINS (X(I,J,K),NLAT,NLON,XI,YJ,ZK,FX)
	CALL TRILINS (Y(I,J,K),NLAT,NLON,XI,YJ,ZK,FY)
	CALL TRILINS (Z(I,J,K),NLAT,NLON,XI,YJ,ZK,FZ)
        RETURN
      ENDIF
	
!          For intrp:
      CALL TRILIN (X(I,J,K),NLAT,NLON,XI,YJ,ZK,FX,DFXDN,DFXDE,DFXDD)
      DFXDTH = -DFXDN*RTOD/DLAT
      DFXDLN =  DFXDE*RTOD/DLON
      DFXDH  = -HTI*HTI*DFXDD/(RE*DIHT)
      CALL TRILIN (Y(I,J,K),NLAT,NLON,XI,YJ,ZK,FY,DFYDN,DFYDE,DFYDD)
      DFYDTH = -DFYDN*RTOD/DLAT
      DFYDLN =  DFYDE*RTOD/DLON
      DFYDH  = -HTI*HTI*DFYDD/(RE*DIHT)
      CALL TRILIN (Z(I,J,K),NLAT,NLON,XI,YJ,ZK,FZ,DFZDN,DFZDE,DFZDD)
      DFZDTH = -DFZDN*RTOD/DLAT
      DFZDLN =  DFZDE*RTOD/DLON
      DFZDH  = -HTI*HTI*DFZDD/(RE*DIHT)
      CALL TRILIN (V(I,J,K),NLAT,NLON,XI,YJ,ZK,FV,DFVDN,DFVDE,DFVDD)
      DFVDTH = -DFVDN*RTOD/DLAT
      DFVDLN =  DFVDE*RTOD/DLON
      DFVDH  = -HTI*HTI*DFVDD/(RE*DIHT)

      IF (NLAT .LT. 3) RETURN

!          Improve calculation of longitudinal derivatives near poles
      IF (GLAT .GE. DLAT-90._r8) GO TO 40
      FAC = .5_r8*XI
      OMFAC = 1._r8 - FAC
      XI = XI - 1._r8
      I = I + 1
      CALL TRILIN (X(I,J,K),NLAT,NLON,XI,YJ,ZK,DMF,DMDFDN,DMDFDE,DMDFDD)
      DFXDLN = DFXDLN*OMFAC + FAC*DMDFDE*RTOD/DLON
      CALL TRILIN (Y(I,J,K),NLAT,NLON,XI,YJ,ZK,DMF,DMDFDN,DMDFDE,DMDFDD)
      DFYDLN = DFYDLN*OMFAC + FAC*DMDFDE*RTOD/DLON
      CALL TRILIN (V(I,J,K),NLAT,NLON,XI,YJ,ZK,DMF,DMDFDN,DMDFDE,DMDFDD)
      DFVDLN = DFVDLN*OMFAC + FAC*DMDFDE*RTOD/DLON

   40 IF (GLAT .LE. 90._r8-DLAT) GO TO 50
      FAC = .5_r8*(1._r8- XI)
      OMFAC = 1._r8 - FAC
      XI = XI + 1._r8
      I = I - 1
      CALL TRILIN (X(I,J,K),NLAT,NLON,XI,YJ,ZK,DMF,DMDFDN,DMDFDE,DMDFDD)
      DFXDLN = DFXDLN*OMFAC + FAC*DMDFDE*RTOD/DLON
      CALL TRILIN (Y(I,J,K),NLAT,NLON,XI,YJ,ZK,DMF,DMDFDN,DMDFDE,DMDFDD)
      DFYDLN = DFYDLN*OMFAC + FAC*DMDFDE*RTOD/DLON
      CALL TRILIN (V(I,J,K),NLAT,NLON,XI,YJ,ZK,DMF,DMDFDN,DMDFDE,DMDFDD)
      DFVDLN = DFVDLN*OMFAC + FAC*DMDFDE*RTOD/DLON
   50 RETURN

!          Error trap diagnostics
 9100 WRITE (MSGU,'(A,'':  Latitude out of range; GPLAT(1),GLAT,GPLAT(NLAT)='',3F10.3)') &
             CALNM,GPLAT(1),GLAT,GPLAT(NLAT)
      IST = 1
      RETURN
 9200 WRITE (MSGU,'(A,'':  Longitude out of range; GPLON(1),GLON,GPLON(NLON)='',3F10.3)') &
             CALNM,GPLON(1),GLON,GPLON(NLON)
      IST = 1
      RETURN
 9300 WRITE (MSGU,'(A,'':  Altitude out of range; GPALT(1),ALT,GPALT(NALT)='',3F10.3)') &
             CALNM,GPALT(1),ALT,GPALT(NALT)
      IST = 1
      RETURN

      END subroutine intrp

!================================================================================================

      SUBROUTINE TRILIN (U,NLAT,NLON,XI,YJ,ZK,FU,DFUDX,DFUDY,DFUDZ)
!
!-----------------------------------------------------------------------
!  Trilinear interpolation of u and its derivatives
! 940803 A. D. Richmond
! Inputs:
!   u(1,1,1) = address of lower corner of interpolation box 
!   nlat = first dimension of u from calling routine
!   nlon = second dimension of u from calling routine
!   xi = fractional distance across box in x direction 
!   yj = fractional distance across box in y direction 
!   zk = fractional distance across box in z direction 
! Outputs:
!   fu = interpolated value of u
!   dfudx = interpolated derivative of u with respect to i (x direction)
!   dfudy = interpolated derivative of u with respect to j (y direction)
!   dfudz = interpolated derivative of u with respect to k (z direction)
!-----------------------------------------------------------------------
!
      use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none
!
!------------------------------Arguments--------------------------------
!
      integer nlat, nlon

      real(r8) U(NLAT,NLON,2)
      real(r8) xi, yj, zk, fu, dfudx, dfudy, dfudz
!
!---------------------------Local variables-----------------------------
!
      integer ientry

      real(r8) omyj, omzk, omxi
!
!-----------------------------------------------------------------------
!
      IENTRY = 0
      GOTO 5
!*******************************************************************************
      ENTRY TRILINS (U,NLAT,NLON,XI,YJ,ZK,FU)
!  Trilinear interpolation of u only
! 940803 A. D. Richmond
! Inputs and outputs:  see above for definitions
      IENTRY = 1

    5 CONTINUE
      OMXI = 1._r8 - XI
      OMYJ = 1._r8 - YJ
      OMZK = 1._r8 - ZK

      FU = U(1,1,1)*OMXI*OMYJ*OMZK &
         + U(2,1,1)*XI*OMYJ*OMZK &
         + U(1,2,1)*OMXI*YJ*OMZK &
         + U(1,1,2)*OMXI*OMYJ*ZK &
         + U(2,2,1)*XI*YJ*OMZK &
         + U(2,1,2)*XI*OMYJ*ZK &
         + U(1,2,2)*OMXI*YJ*ZK &
         + U(2,2,2)*XI*YJ*ZK

      IF (IENTRY.NE.0) RETURN

      DFUDX = (U(2,1,1)-U(1,1,1))*OMYJ*OMZK &
            + (U(2,2,1)-U(1,2,1))*YJ*OMZK &
            + (U(2,1,2)-U(1,1,2))*OMYJ*ZK &
            + (U(2,2,2)-U(1,2,2))*YJ*ZK
      DFUDY = (U(1,2,1)-U(1,1,1))*OMXI*OMZK &
            + (U(2,2,1)-U(2,1,1))*XI*OMZK &
            + (U(1,2,2)-U(1,1,2))*OMXI*ZK &
            + (U(2,2,2)-U(2,1,2))*XI*ZK
      DFUDZ = (U(1,1,2)-U(1,1,1))*OMXI*OMYJ &
            + (U(2,1,2)-U(2,1,1))*XI*OMYJ &
            + (U(1,2,2)-U(1,2,1))*OMXI*YJ &
            + (U(2,2,2)-U(2,2,1))*XI*YJ
      RETURN
      END subroutine trilin

!================================================================================================

      SUBROUTINE ADPL(GLAT,GLON,CTH,STH,FX,FY,FZ,FV &
        ,DFXDTH,DFYDTH,DFZDTH,DFVDTH,DFXDLN,DFYDLN,DFZDLN,DFVDLN)
!
!-----------------------------------------------------------------------
!  v is used for vr2n
!  Add-back of pseudodipole component to x,y,z,v and their derivatives.
!  940715 A. D. Richmond
! Inputs:
!   glat = latitude, degrees
!   glon = longitude, degrees
!   fx = interpolated value of x
!   fy = interpolated value of y
!   fz = interpolated value of z
!   fv = interpolated value of v
!   dfxdth,dfydth,dfzdth,dfvdth = derivatives of x,y,z,v with respect to 
!	colatitude, in radians-1
!   dfxdln,dfydln,dfzdln,dfvdln = derivatives of x,y,z,v with respect to 
!	longitude, in radians-1
! Output:
!   cth,sth = cos(colatitude), sin(colatitude)
!   fx = interpolated value of x
!   fy = interpolated value of y
!   fz = interpolated value of z
!   fv = interpolated value of v
!   dfxdth,dfydth,dfzdth,dfvdth = derivatives of x,y,z,v with respect to 
!	colatitude, in radians-1
!   dfxdln,dfydln,dfzdln,dfvdln = derivatives of x,y,z,v with respect to 
!	longitude, in radians-1
!-----------------------------------------------------------------------
!
      use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none
!
!-------------------------------Commons---------------------------------
!
! bf: separate ints from reals in commons:
      integer msgu
      COMMON /APXCON_int/ MSGU

      real(r8) rtod, dtor, re, req, pola
      COMMON /APXCON/ RTOD,DTOR,RE,REQ,POLA

      real(r8) colat, elon, vp, ctp, stp
      COMMON /APXDIPL/ COLAT,ELON,VP,CTP,STP
!
!------------------------------Arguments--------------------------------
!
      real(r8) glat, glon, cth, sth, fx, fy, fz, fv
      real(r8) dfxdth, dfydth, dfzdth, dfvdth, dfxdln
      real(r8) dfydln, dfzdln, dfvdln
!
!---------------------------Local variables-----------------------------
!
      real(r8) ctm, sph, cph
!
!-----------------------------------------------------------------------
!
      CPH = COS((GLON-ELON)*DTOR)
      SPH = SIN((GLON-ELON)*DTOR)
      CTH = SIN(GLAT*DTOR)
      STH = COS(GLAT*DTOR)
      CTM = CTP*CTH + STP*STH*CPH
      FX = FX + STH*CTP*CPH - CTH*STP
      FY = FY + STH*SPH
      FZ = FZ + CTM
      FV = FV - CTM
      DFXDTH = DFXDTH + CTP*CTH*CPH + STP*STH
      DFYDTH = DFYDTH + CTH*SPH
      DFZDTH = DFZDTH - CTP*STH + STP*CTH*CPH
      DFVDTH = DFVDTH + CTP*STH - STP*CTH*CPH
      DFXDLN = DFXDLN - CTP*STH*SPH
      DFYDLN = DFYDLN + STH*CPH
      DFZDLN = DFZDLN - STP*STH*SPH
      DFVDLN = DFVDLN + STP*STH*SPH
      RETURN
      END subroutine adpl

!================================================================================================

      SUBROUTINE ADPLSC (GLAT,GLON,FX,FY,FZ)
!
!-----------------------------------------------------------------------
!  Add-back of pseudodipole component to x,y,z 
!  940801 A. D. Richmond
! Inputs:
!   glat = latitude, degrees
!   glon = longitude, degrees
!   fx = interpolated value of x
!   fy = interpolated value of y
!   fz = interpolated value of z
! Output:
!   fx = interpolated value of x
!   fy = interpolated value of y
!   fz = interpolated value of z
!-----------------------------------------------------------------------
!
      use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none
!
!-------------------------------Commons---------------------------------
!
! bf: separate ints from reals in commons:
      integer msgu
      COMMON /APXCON_int/ MSGU

      real(r8) rtod, dtor, re, req, pola
      COMMON /APXCON/ RTOD,DTOR,RE,REQ,POLA

      real(r8) colat, elon, vp, ctp, stp
      COMMON /APXDIPL/ COLAT,ELON,VP,CTP,STP
!
!------------------------------Arguments--------------------------------
!
      real(r8) glat, glon, fx, fy, fz
!
!---------------------------Local variables-----------------------------
!
      real(r8) sth, ctm, cth, cph, sph 
!
!-----------------------------------------------------------------------
!
      CPH = COS((GLON-ELON)*DTOR)
      SPH = SIN((GLON-ELON)*DTOR)
      CTH = SIN(GLAT*DTOR)
      STH = COS(GLAT*DTOR)
      CTM = CTP*CTH + STP*STH*CPH
      FX = FX + STH*CTP*CPH - CTH*STP
      FY = FY + STH*SPH
      FZ = FZ + CTM
      RETURN
      END subroutine adplsc

!================================================================================================

      SUBROUTINE GRADXYZV (ALT,CTH,STH, &
                DFXDTH,DFYDTH,DFZDTH,DFVDTH,DFXDLN,DFYDLN,DFZDLN,DFVDLN, &
                DFXDH,DFYDH,DFZDH,DFVDH,GRADX,GRADY,GRADZ,GRADV)
!
!-----------------------------------------------------------------------
!          Calculates east,north,up components of gradients of x,y,z,v in
!          geodetic coordinates.  All gradients are in inverse km.  Assumes
!          flatness of 1/298.25 and equatorial radius (REQ) of 6378.16 km.
!          940803 A. D. Richmond
!-----------------------------------------------------------------------
!
      use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none
!
!-------------------------------Commons---------------------------------
!
      real(r8) rho, ddisdth
      COMMON /APXGEOD/ RHO,DDISDTH
!
!------------------------------Arguments--------------------------------
!
      real(r8) alt, cth, sth, dfxdth, dfydth, dfzdth, dfvdth
      real(r8) dfxdln, dfydln, dfzdln, dfvdln, dfxdh, dfydh, dfzdh, dfvdh
      real(r8) GRADX(3),GRADY(3),GRADZ(3),GRADV(3)
!
!---------------------------Local variables-----------------------------
!
      real(r8) dddthod, drhodth, dzetdth, d2, d
      integer ientry
!
!-----------------------------------------------------------------------
!
      IENTRY = 0
      GOTO 5
!*******************************************************************************
      ENTRY GRAPXYZV (ALT,CTH,STH, &
                    DFXDLN,DFYDLN,DFZDLN,DFVDLN,GRADX,GRADY,GRADZ,GRADV)
!          Calculates east component of gradient near pole.
!          940803 A. D. Richmond
!          Inputs and outputs:  see above for definitions
      IENTRY = 1

    5 CONTINUE
      D2 = 40680925.E0_r8 - 272340.E0_r8*CTH*CTH
!          40680925. = req**2 (rounded off)
!          272340.   = req**2 * E2, where E2 = (2. - 1./298.25)/298.25
!                      is square of eccentricity of ellipsoid.
      D = SQRT(D2)
      RHO = STH*(ALT + 40680925.E0_r8/D)
      DDDTHOD = 272340.E0_r8*CTH*STH/D2
      DRHODTH = ALT*CTH + (40680925.E0_r8/D)*(CTH-STH*DDDTHOD)
      DZETDTH =-ALT*STH - (40408585.E0_r8/D)*(STH+CTH*DDDTHOD)
      DDISDTH = SQRT(DRHODTH*DRHODTH + DZETDTH*DZETDTH)
      GRADX(1) = DFXDLN/RHO
      GRADY(1) = DFYDLN/RHO
      GRADZ(1) = DFZDLN/RHO
      GRADV(1) = DFVDLN/RHO

      IF (IENTRY .NE. 0) RETURN

      GRADX(2) = -DFXDTH/DDISDTH
      GRADY(2) = -DFYDTH/DDISDTH
      GRADZ(2) = -DFZDTH/DDISDTH
      GRADV(2) = -DFVDTH/DDISDTH
      GRADX(3) = DFXDH
      GRADY(3) = DFYDH
      GRADZ(3) = DFZDH
      GRADV(3) = DFVDH

      RETURN
      end subroutine gradxyzv

!================================================================================================

      SUBROUTINE GRADLPV (HR,ALT,FX,FY,FZ,FV,GRADX,GRADY,GRADZ,GRADV, &
                    XLATM,XLONM,VMP,GRCLM,CLMGRP,QDLAT,RGRLP,B,CLM,R3_2)
!
!-----------------------------------------------------------------------
!          Uses gradients of x,y,z,v to compute geomagnetic field and
!          gradients of apex latitude, longitude.
!          940819 A. D. Richmond
!          INPUT:
!            HR     = reference altitude
!            ALT    = altitude
!            FX,FY,FZ,FV = interpolated values of x,y,z,v, plus pseudodipole
!                     component
!            GRADX,GRADY,GRADZ,GRADV = interpolated gradients of x,y,z,v,
!                     including pseudodipole components (east,north,up)
!          OUTPUT:
!            XLATM  = modified apex latitude (lambda_m), degrees
!            XLONM  = apex longitude (phi_a), degrees
!            VMP    = magnetic potential, in T.m.
!            GRCLM  = grad(cos(lambda_m)), in km-1
!            CLMGRP = cos(lambda_m)*grad(phi_a), in km-1
!            QDLAT  = quasi-dipole latitude, degrees
!            RGRLP  = (re + alt)*grad(lambda')
!            B      = magnetic field, in nT
!            CLM    = cos(lambda_m)
!            R3_2   = ((re + alt)/(re + hr))**(3/2)
!-----------------------------------------------------------------------
!
      use shr_kind_mod, only: r8 => shr_kind_r8
      use abortutils,   only: endrun
      use cam_logfile,  only: iulog
      implicit none
!
!-------------------------------Commons---------------------------------
!
!
! bf: separate ints from reals in commons:
      integer msgu
      COMMON /APXCON_int/ MSGU

      real(r8) rtod, dtor, re, req, pola
      COMMON /APXCON/ RTOD,DTOR,RE,REQ,POLA

      real(r8) colat, elon, vp, ctp, stp
      COMMON /APXDIPL/ COLAT,ELON,VP,CTP,STP
!
!------------------------------Arguments--------------------------------
!
      real(r8) hr, alt, fx, fy, fz, fv, xlatm, xlonm, vmp, qdlat
      real(r8) clm, r3_2
      real(r8) GRADX(3),GRADY(3),GRADZ(3),GRADV(3), &
               GRCLM(3),CLMGRP(3),RGRLP(3),B(3)
!
!---------------------------Local variables-----------------------------
!
      integer i, ierr

      real(r8) xlp, slp, xnorm, rn2, x2py2, clp, clp2, am1, slp2
      real(r8) grclp, rr, alon, a, alat, r, spm, bo, cpm, rn, sqrror
!
!-----------------------------------------------------------------------
!
      RR = RE + HR
      R  = RE + ALT
      RN = R/RE
      SQRROR = SQRT(RR/R)
      R3_2 = 1._r8/SQRROR/SQRROR/SQRROR
      XLONM = ATAN2(FY,FX)
      CPM = COS(XLONM)
      SPM = SIN(XLONM)
      XLONM = RTOD*XLONM
      BO = VP*1.E6_r8
!             1.E6 converts T to nT and km-1 to m-1.
      RN2 = RN*RN
      VMP = VP*FV/RN2
      B(1) = -BO*GRADV(1)/RN2
      B(2) = -BO*GRADV(2)/RN2
      B(3) = -BO*(GRADV(3)-2._r8*FV/R)/RN2

      X2PY2 = FX*FX + FY*FY
      XNORM = SQRT(X2PY2 + FZ*FZ)
      XLP = ATAN2(FZ,SQRT(X2PY2))
      SLP = SIN(XLP)
      CLP = COS(XLP)
      QDLAT = XLP*RTOD
      CLM = SQRROR*CLP
      IF (CLM.LE.1._r8) GOTO 5
      write(iulog,*) 'Stopped in gradlpv because point lies below field line that peaks at reference height.'
      call endrun
    5 XLATM = RTOD*ACOS(CLM)
!  If southern magnetic hemisphere, reverse sign of xlatm
      IF (SLP.LT.0._r8) XLATM = - XLATM
      DO 10 I=1,3
	GRCLP = CPM*GRADX(I) + SPM*GRADY(I)
	RGRLP(I) = R*(CLP*GRADZ(I) - SLP*GRCLP)
	GRCLM(I) = SQRROR*GRCLP
   10   CLMGRP(I) = SQRROR*(CPM*GRADY(I)-SPM*GRADX(I))
      GRCLM(3) = GRCLM(3) - SQRROR*CLP/(2._r8*R)
      RETURN
!*******************************************************************************
      ENTRY XYZ2APX (ALT,FX,FY,FZ,A,ALAT,ALON,IERR)
!          Computes apex latitude, longitude.
!          990309 A. D. Richmond
!          INPUT:
!            ALT      = altitude
!            FX,FY,FZ = interpolated values of x,y,z, plus pseudodipole
!                       component
!          OUTPUT:
!            A    = apex radius, normalized by req
!            ALAT = apex latitude, degrees
!            ALON = apex longitude, degrees
!
!          Mod (Mar 99):  Lines 19-30 are changed from the original in order
!          to avoid a truncation error near the magnetic equator.  What was
!          done is to make use of the identity
!
!                  SIN(ALAT/RTOD)**2 + COS(ALAT/RTOD)**2 = 1,
!
!          so that since
!
!                  COS(ALAT/RTOD)**2 = 1/A (Eq. 1),
!
!          then
!
!                  SIN(ALAT/RTOD)**2 = (A-1)/A (Eq. 2)
!
!          Below AM1 = A-1.  If AM1 is less than 1, use Eq. 1;
!          otherwise use Eq. 2.  Mathematically, both equations should
!          give identical results, but numerically it is better to use
!          that function ASIN or ACOS that has the smaller argument.
!          The jump from one equation to the other occurs at ALAT = 45.

      IERR  = 0
      ALON  = ATAN2(FY,FX)*RTOD
      SLP2  = FZ*FZ
      X2PY2 = FX*FX + FY*FY
      XNORM = SLP2 + X2PY2
      SLP2  = SLP2/XNORM
      CLP2  = X2PY2/XNORM
      AM1   = (RE*SLP2 + ALT)/(REQ*CLP2)
      A = 1._r8 + AM1

      IF (AM1.LT.0._r8) THEN
        IERR = 1
        write(iulog,*) 'Missing alat returned because point lies below field line that peaks at Earth surface.'
        RETURN
      ELSEIF (AM1.LT.1._r8) THEN
        ALAT = RTOD*ASIN(SQRT(AM1/A))
      ELSE
        ALAT = RTOD*ACOS(1._r8/SQRT(A))
      ENDIF
!  If southern magnetic hemisphere, reverse sign of alat
      IF (FZ.LT.0._r8) ALAT = - ALAT
      RETURN
      end subroutine gradlpv

!================================================================================================

      SUBROUTINE BASVEC (HR,XLATM,GRCLM,CLMGRP,RGRLP,B,CLM,R3_2, &
                        BMAG,SIM,SI,F,D,W,BHAT,D1,D2,D3,E1,E2,E3,F1,F2)
!
!-----------------------------------------------------------------------
!          Computes base vectors and other parameters for apex coordinates.
!          Vector components:  east, north, up
!          940801 A. D. Richmond
!          Reference:
!            Richmond, A. D., Ionospheric Electrodynamics Using Magnetic Apex
!            Coordinates, J. Geomag. Geoelectr., 47, 191-212, 1995.
!          INPUTS:
!            HR     = reference altitude
!            XLATM  = modified apex latitude, degrees
!            GRCLM  = grad(cos(lambda_m)), in km-1
!            CLMGRP = cos(lambda_m)*grad(phi_a), in km-1
!            RGRLP  = (re + altitude)*grad(lambda')
!            B      = magnetic field, in nT
!            CLM    = cos(lambda_m)
!            R3_2   = ((re + altitude)/(re + hr))**(3/2)
!          RETURNS:
!            BMAG    = magnitude of magnetic field, in nT
!            SIM     = sin(I_m) of article
!            SI      = sin(I)
!            F       = F of article
!            D       = D of article
!            W       = W of article
!            BHAT    = unit vector along geomagnetic field direction
!            D1...F2 = base vectors of article
!-----------------------------------------------------------------------
!
      use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none
!
!-------------------------------Commons---------------------------------
!
! bf: separate ints from reals in commons:
      integer msgu
      COMMON /APXCON_int/ MSGU

      real(r8) rtod, dtor, re, req, pola
      COMMON /APXCON/ RTOD,DTOR,RE,REQ,POLA
!
!------------------------------Arguments--------------------------------
!
      real(r8) hr, xlatm
      real(r8) GRCLM(3),CLMGRP(3),RGRLP(3),B(3)
      real(r8) clm, r3_2, bmag, sim, si, f, d, w
      real(r8) BHAT(3),D1(3),D2(3),D3(3),E1(3),E2(3),E3(3),F1(2),F2(2)
!
!---------------------------Local variables-----------------------------
!
      integer i

      real(r8) d2db, d1db, rr, simoslm
!
!-----------------------------------------------------------------------
!
      RR = RE + HR
      SIMOSLM = 2._r8/SQRT(4._r8 - 3._r8*CLM*CLM)
      SIM = SIMOSLM*SIN(XLATM*DTOR)
      BMAG = SQRT(B(1)*B(1) + B(2)*B(2) + B(3)*B(3))
      D1DB = 0._r8
      D2DB = 0._r8
      DO 10 I=1,3
        BHAT(I) = B(I)/BMAG
        D1(I) = RR*CLMGRP(I)
        D1DB = D1DB + D1(I)*BHAT(I)
        D2(I) = RR*SIMOSLM*GRCLM(I)
   10   D2DB = D2DB + D2(I)*BHAT(I)
! Ensure that d1,d2 are exactly perpendicular to B:
      DO 15 I=1,3
        D1(I) = D1(I) - D1DB*BHAT(I)
   15   D2(I) = D2(I) - D2DB*BHAT(I)
      E3(1) = D1(2)*D2(3) - D1(3)*D2(2)
      E3(2) = D1(3)*D2(1) - D1(1)*D2(3)
      E3(3) = D1(1)*D2(2) - D1(2)*D2(1)
      D = BHAT(1)*E3(1) + BHAT(2)*E3(2) + BHAT(3)*E3(3)
      DO 20 I=1,3
        D3(I) = BHAT(I)/D
! Following step may be unnecessary, but it ensures that e3 lies along bhat.
        E3(I) = BHAT(I)*D
   20   CONTINUE
      E1(1) = D2(2)*D3(3) - D2(3)*D3(2)
      E1(2) = D2(3)*D3(1) - D2(1)*D3(3)
      E1(3) = D2(1)*D3(2) - D2(2)*D3(1)
      E2(1) = D3(2)*D1(3) - D3(3)*D1(2)
      E2(2) = D3(3)*D1(1) - D3(1)*D1(3)
      E2(3) = D3(1)*D1(2) - D3(2)*D1(1)
      W = RR*RR*CLM*ABS(SIM)/(BMAG*D)
      SI = -BHAT(3)
      F1(1) =  RGRLP(2) 
      F1(2) = -RGRLP(1)
      F2(1) = -D1(2)*R3_2
      F2(2) =  D1(1)*R3_2
      F = F1(1)*F2(2) - F1(2)*F2(1)
      RETURN
      end subroutine basvec

!================================================================================================

      SUBROUTINE CKGP (CALNM,MSGUN,NLAT,NLON,NALT,GPLAT,GPLON,GPALT,IST)
!
!-----------------------------------------------------------------------
!          Check grid point values tests extremes and order of the grid
!          point arrays, producing diagnostics to MSGUN and IST=1 when
!          rules have been broken.
!-----------------------------------------------------------------------
!
      use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none
!
!------------------------------Arguments--------------------------------
!
      integer msgun, nlat, nlon, nalt, ist

      real(r8) GPLAT(NLAT), GPLON(NLON), GPALT(NALT)

      CHARACTER*(*) CALNM
!
!---------------------------Local variables-----------------------------
!
      integer i

      real(r8) olat, oalt, olon
!
!-----------------------------------------------------------------------
!
      IST = 1
      OLAT = -90._r8
      DO 10 I=1,NLAT
      IF (ABS (GPLAT(I)) .GT.  90._r8) GO TO 9100
      IF (     GPLAT(I)  .LT. OLAT) GO TO 9200
   10 OLAT = GPLAT(I)
!
! 2/00: reset longitude limits to -270 -> +270, as per communication
!       with Richmond, 2/25/00.
!
!     OLON = -180.
      OLON = -270._r8
      DO 20 I=1,NLON
!     IF (ABS (GPLON(I)) .GT. 180.) GO TO 9300
      IF (ABS (GPLON(I)) .GT. 270._r8) GO TO 9300
      IF (     GPLON(I)  .LT. OLON) GO TO 9400
   20 OLON = GPLON(I)

      OALT = 0._r8
      DO 30 I=1,NALT
      IF (GPALT(I) .LT.   0._r8) GO TO 9500
      IF (GPALT(I) .LT. OALT) GO TO 9600
   30 OALT = GPALT(I)

      IST = 0

  100 RETURN

 9100 WRITE (MSGUN,'(A,'':  |GPLAT(I)| > 90; I,GPLAT(I)'',I5,F10.3)') &
                  CALNM,                     I,GPLAT(I)
      GO TO 100
 9200 WRITE (MSGUN,'(A,'':  GPLAT(I) < GPLAT(I-1); I,GPLAT(I),GPLAT(I-1) ='',I5,2F10.3)') &
                     CALNM,  I,GPLAT(I),OLAT
      GO TO 100
 9300 WRITE (MSGUN,'(A,'':  |GPLON(I)| > 180; I,GPLON(I)'',I5,F10.3)') &
                                       CALNM, I,GPLON(I)
      GO TO 100
 9400 WRITE (MSGUN,'(A,'':  GPLON(I) < GPLON(I-1); I,GPLON(I),GPLON(I-1) ='',I5,2F10.3)') &
                   CALNM, I,GPLON(I),OLON
      GO TO 100
 9500 WRITE (MSGUN,'(A,'':  GPALT(I) <  0; I,GPALT(I)'',I5,F10.3)') &
                                      CALNM, I,GPALT(I)
      GO TO 100
 9600 WRITE (MSGUN,'(A,'':  GPALT(I) < GPALT(I-1); I,GPALT(I),GPALT(I-1) ='',I5,2F10.3)') &
                     CALNM, I,GPALT(I),OALT
      GO TO 100
      END subroutine ckgp

!================================================================================================

      SUBROUTINE MAKEXYZV (EPOCH,NLAT,NLON,NALT,GPLAT,GPLON,GPALT, &
                           X,Y,Z,V)
!
!-----------------------------------------------------------------------
!          Sets up grid arrays for later interpolation
!          940822 A. D. Richmond, NCAR
!          INPUT:
!            EPOCH = year and fraction (e.g., 1994.50 for 1994 July 2)
!            NLAT,NLON,NALT = triple dimensions of X,Y,Z,V and respective
!                    single dimensions of GP___ arrays
!            GPLAT,GPLON,GPALT = grid point latitudes, longitudes and altitudes
!          OUTPUT:
!            X = array containing cos(lambda')cos(phi_a) less pseudodipole
!                component
!            Y = array containing cos(lambda')sin(phi_a) less pseudodipole
!                component
!            Z = array containing sin(lambda') less pseudodipole component
!            V = array containing ((magnetic potential)/vp)*((re+height)/re)**2,
!                less pseudodipole component
!
!          Modification (99 Mar):  Make X,Y,Z,V constant near the poles
!          for all GPLON(j) at each height.  Add POLA to APXCON
!            POLA       = Pole angle (deg); when the geographic latitude is
!                         poleward of POLA, X,Y,Z,V are forced to be constant.
!                         for all longitudes at each altitude.  POLA is defined
!                         in APXMKA (POLA = 90. - SQRT (PRECISE) * RTOD),
!                         which currently makes POLA = 89.995
!-----------------------------------------------------------------------
!
      use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none
!
!-------------------------------Commons---------------------------------
!
! bf: separate ints from reals in commons:
      integer msgu
      COMMON /APXCON_int/ MSGU

      real(r8) rtod, dtor, re, req, pola
      COMMON /APXCON/ RTOD,DTOR,RE,REQ,POLA

      real(r8) colat, elon, vp, ctp, stp
      COMMON /APXDIPL/ COLAT,ELON,VP,CTP,STP
!
!------------------------------Arguments--------------------------------
!
      integer nlat, nlon, nalt

      real(r8) epoch
      real(r8) X(NLAT,NLON,NALT), Y(NLAT,NLON,NALT), Z(NLAT,NLON,NALT), &
           V(NLAT,NLON,NALT), GPLAT(NLAT), GPLON(NLON), GPALT(NALT)
!
!---------------------------Local variables-----------------------------
!
      integer kpol, i, j, k

      real(r8) xmag, ymag, zdown, alat, phia, bmag, vmp, slp, clp, phiar
      real(r8)  vnor, rp, reqam1, a, ct, st, reqore, rqorm1, ctm, stmcpm
      real(r8)  stmspm, cp, sp
!
!-----------------------------------------------------------------------
!
      CALL COFRM (EPOCH)
      CALL DYPOL (COLAT,ELON,VP)
      CTP = COS (COLAT*DTOR)
      STP = SIN (COLAT*DTOR)
      REQORE = REQ/RE
      RQORM1 = REQORE-1._r8

      DO 100 I=1,NLAT
      CT = SIN (GPLAT(I)*DTOR)
      ST = COS (GPLAT(I)*DTOR)
      KPOL = 0
      IF (ABS (GPLAT(I)) .GT. POLA) KPOL = 1

      DO 100 J=1,NLON
      IF (KPOL .EQ. 0) GO TO 20
      IF (J    .EQ. 1) GO TO 20
!          KPOL = 1 (poleward of POLA) and first lon's XYZV are defined
      DO 10 K=1,NALT
      V(I,J,K) = V(I,1,K)
      X(I,J,K) = X(I,1,K)
      Y(I,J,K) = Y(I,1,K)
   10 Z(I,J,K) = Z(I,1,K)
      GO TO 100

   20 CP  = COS ((GPLON(J)-ELON)*DTOR)
      SP  = SIN ((GPLON(J)-ELON)*DTOR)
!           ctm   is pseudodipole component of z
!          -ctm   is pseudodipole component of v
!          stmcpm is pseudodipole component of x
!          stmspm is pseudodipole component of y
      CTM    = CTP*CT + STP*ST*CP
      STMCPM = ST*CTP*CP - CT*STP
      STMSPM = ST*SP

      DO 30 K=1,NALT
      CALL APEX (EPOCH,GPLAT(I),GPLON(J),GPALT(K), &
                 A,ALAT,PHIA,BMAG,XMAG,YMAG,ZDOWN,VMP)
      VNOR = VMP/VP
      RP = 1._r8 + GPALT(K)/RE
      V(I,J,K) = VNOR*RP*RP + CTM
      REQAM1 = REQ*(A-1._r8)
      SLP = SQRT(max(REQAM1-GPALT(K),0._r8)/(REQAM1+RE))
!          Reverse sign of slp in southern magnetic hemisphere
      IF (ZDOWN.LT.0._r8) SLP = -SLP
      CLP = SQRT (RP/(REQORE*A-RQORM1))
      PHIAR = PHIA*DTOR
      X(I,J,K) = CLP*COS (PHIAR) - STMCPM
      Y(I,J,K) = CLP*SIN (PHIAR) - STMSPM
      Z(I,J,K) = SLP - CTM
   30 CONTINUE

  100 CONTINUE

      RETURN
      end subroutine makexyzv

!================================================================================================

      SUBROUTINE SETMISS (XMISS &
       ,XLATM,XLONM,VMP,B,BMAG,BE3,SIM,SI,F,D,W &
       ,BHAT,D1,D2,D3,E1,E2,E3,F1,F2)
!
!-----------------------------------------------------------------------
!
      use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none
!
!------------------------------Arguments--------------------------------
!
      real(r8) xmiss, xlatm, xlonm, vmp, bmag, be3, sim, si, f, d, w
      real(r8) BHAT(3),D1(3),D2(3),D3(3),E1(3),E2(3),E3(3),F1(2),F2(2)
      real(r8) B(3)
!
!---------------------------Local variables-----------------------------
!
      integer i
!
!-----------------------------------------------------------------------
!
      XLATM = XMISS
      XLONM = XMISS
      VMP = XMISS
      BMAG = XMISS
      BE3 = XMISS
      SIM = XMISS
      SI = XMISS
      F = XMISS
      D = XMISS
      W = XMISS
      DO 5 I=1,3
	B(I) = XMISS
	BHAT(I) = XMISS
	D1(I) = XMISS
	D2(I) = XMISS
	D3(I) = XMISS
	E1(I) = XMISS
	E2(I) = XMISS
    5   E3(I) = XMISS
      DO 6 I=1,2
	F1(I) = XMISS
    6   F2(I) = XMISS
      RETURN
      END subroutine setmiss

!================================================================================================

      SUBROUTINE GM2GC (GMLAT,GMLON,GCLAT,GCLON)
!
!-----------------------------------------------------------------------
!  Converts geomagnetic to geocentric coordinates.
!  940819 A. D. Richmond
!
!  Inputs:
!	gmlat = geomagnetic latitude in degrees
!	gmlon = geomagnetic longitude in degrees
!  Outputs:
!	gclat = geocentric latitude in degrees
!	gclon = geocentric longitude in degrees
!
!  Common/consts/
!	rtod, dtor = 180/pi, pi/180
!	re, req    = 6371.2, 6378.160
!-----------------------------------------------------------------------
!
      use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none
!
!-------------------------------Commons---------------------------------
!
! bf: separate ints from reals in commons:
      integer msgu
      COMMON /APXCON_int/ MSGU

      real(r8) rtod, dtor, re, req, pola
      COMMON /APXCON/ RTOD,DTOR,RE,REQ,POLA

!  Common/dipol/
!       colat = geocentric colatitude of north geomagnetic pole, in degrees
!	elon  = geocentric east longitude of north geomagnetic pole, in degrees
!	vp    = magnetic potential at 1 RE, south geomagnetic pole
!	ctp   = cos(colat*dtor)
!	stp   = sin(colat*dtor)

      real(r8) colat, elon, vp, ctp, stp
      COMMON /APXDIPL/ COLAT,ELON,VP,CTP,STP
!
!------------------------------Arguments--------------------------------
!
      real(r8) gmlat, gmlon, gclat, gclon
!
!---------------------------Local variables-----------------------------
!
      real(r8) ctc, ctm, stm
!
!-----------------------------------------------------------------------
!
      STM = COS(GMLAT*DTOR)
      CTM = SIN(GMLAT*DTOR)
      CTC = CTP*CTM - STP*STM*COS(GMLON*DTOR)
      CTC = min(CTC,1._r8)
      CTC = max(CTC,-1._r8)
      GCLAT = ASIN(CTC)*RTOD
      GCLON = ATAN2(STP*STM*SIN(GMLON*DTOR),CTM-CTP*CTC)
      GCLON = GCLON*RTOD + ELON
      IF (GCLON.LT.-180._r8) GCLON = GCLON + 360._r8
      RETURN
      END subroutine gm2gc

!================================================================================================

      SUBROUTINE APEX (DATE,DLAT,DLON,ALT, &
                      A,ALAT,ALON,BMAG,XMAG,YMAG,ZMAG,V)
!
!-----------------------------------------------------------------------
!
! Original file ~bozo/pgms/apex/apex.f copied 2/25/00.
!
! Calculate apex radius, latitude, longitude; and magnetic field and potential
! 940822 A. D. Richmond
!
!          INPUTS:
!            DATE = Year and fraction (1990.0 = 1990 January 1, 0 UT)
!            DLAT = Geodetic latitude in degrees
!            DLON = Geodetic longitude in degrees
!            ALT = Altitude in km
!
!          RETURNS:
!            A    = (Apex height + REQ)/REQ, where REQ = equatorial Earth
!                    radius.  A is analogous to the L value in invariant
!                    coordinates.
!            ALAT = Apex latitude in degrees (negative in S. magnetic hemisphere)
!            ALON = Apex longitude (Geomagnetic longitude of apex)
!            BMAG = Magnitude of geomagnetic field in nT
!            XMAG,YMAG,ZMAG = North, east, and downward geomagnetic field components in nT
!            V    = Magnetic potential, in T.m
!
!     COMMON /DIPOLE/ COLAT,ELON,VP,CTP,STP
!          COLAT = Geocentric colatitude of geomagnetic dipole north pole (deg)
!          ELON  = East longitude of geomagnetic dipole north pole (deg)
!          VP    = Magnitude, in T.m, of dipole component of magnetic potential at
!                  geomagnetic pole and geocentric radius of 6371.2 km
!          CTP,STP = cosine, sine of COLAT
!
!          MODIFICATIONS:
!          May 1999:  Revise DS calculation in LINAPX to avoid divide by zero.
!-----------------------------------------------------------------------
!
      use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none
!
!-------------------------------Commons---------------------------------
!
      real(r8) colat, elon, vp, ctp, stp
      COMMON /DIPOLE/ COLAT,ELON,VP,CTP,STP
!
!------------------------------Arguments--------------------------------
!
      real(r8) date, dlat, dlon, alt, a, alat, alon, bmag, xmag, ymag, zmag, v
!
!-----------------------------Parameters------------------------------
!
      real(r8) rtod, dtor, re, req
      PARAMETER (RTOD=5.72957795130823E1_r8, DTOR=1.745329251994330E-2_r8, &
                 RE=6371.2_r8, REQ=6378.160_r8)
!
!---------------------------Local variables-----------------------------
!
      real(r8) bx, z, bz, by, y, polon, clatp, x, vpol
!
!-----------------------------------------------------------------------
!
      CALL COFRM (DATE)
      CALL DYPOL (CLATP,POLON,VPOL)
      COLAT = CLATP
      CTP   = COS(CLATP*DTOR)
      STP   = SQRT(1._r8 - CTP*CTP)
      ELON  = POLON
      VP    = VPOL
      CALL LINAPX (DLAT,DLON,ALT,A,ALAT,ALON,XMAG,YMAG,ZMAG,BMAG)
      XMAG = XMAG*1.E5_r8
      YMAG = YMAG*1.E5_r8
      ZMAG = ZMAG*1.E5_r8
      BMAG = BMAG*1.E5_r8
      CALL GD2CART (DLAT,DLON,ALT,X,Y,Z)
      CALL FELDG (3,X/RE,Y/RE,Z/RE,BX,BY,BZ,V)
      RETURN
      END subroutine apex

!================================================================================================

       SUBROUTINE LINAPX(GDLAT,GLON,ALT &
       ,A,ALAT,ALON,XMAG,YMAG,ZMAG,F)
!
!-----------------------------------------------------------------------
!***BEGIN PROLOGUE  LINAPX                                                      
!***DATE WRITTEN   731029   (YYMMDD)                                            
!***AUTHOR  CLARK, W., N.O.A.A. ERL LAB.                                        
!***REVISION DATE  880201   (YYMMDD) Harsh Anand Passi, NCAR
!***LATEST REVISION 940803 A. D. Richmond
!***PURPOSE  Transforms the geographic coordinates to apex coordinates.         
!***DESCRIPTION                                                                 
!     The method used is as follow:                                             
!       1. Calculates step size as a function of the geomagnetic                
!          dipole coordinates of the starting point.                            
!       2. Determine direction of trace                                         
!       3. Convert the geodetic coordinates of the starting point               
!          to the cartesian coordinates for tracing.                            
!       Loop:                                                                   
!       i)   Increment step count, if count > 200,                              
!            assume it is dipole field, call DIPAPX to                          
!            determine Apex coordinates else continue.                          
!       ii)  Get field components in X, Y and Z direction                       
!       iii) Call ITRACE to trace field line.                                   
!       iv)  Test if Apex passed call FNDAPX to determine Apex coordinates      
!            else loop:                                                         
!                                                                               
!   INPUT                                                                       
!     GDLAT  Latitude of starting point (Deg)                                   
!     GLON   Longitude (East=+) of starting point (Deg)                         
!     ALT    Ht. of starting point (Km)                                         
!                                                                               
!  OUTPUT
!          A  (Apex height + REQ)/REQ, where REQ = equatorial Earth radius.
!             A is analogous to the L value in invariant coordinates.
!       ALAT  Apex Lat. (deg)                                                   
!       ALON  Apex Lon. (deg)                                                   
!       XMAG  North component of magnetic field at starting point
!       YMAG  East component of magnetic field at starting point
!       ZMAG  Down component of magnetic field at starting point
!          F  Magnetic field magnitude at starting point
!
!     COMMON Blocks Used                                                        
!     /APXIN/ YAPX(3,3)                                                         
!       YAPX    Matrix of cartesian coordinates (loaded columnwise)             
!               of the 3 points about APEX. Set in subprogram ITRACE.           
!                                                                               
!   /DIPOLE/COLAT,ELON,VP,CTP,STP                           
!     COLAT   Geographic colatitude of the north pole of the                    
!             earth-centered dipole (Deg).                                      
!     ELON    Geographic longitude of the north pole of the                     
!             earth-centered dipole (Deg).                                      
!     VP      Magnetic potential magnitude at geomagnetic pole (T.m)
!     CTP     cos(COLAT*DTOR)
!     STP     sin(COLAT*DTOR)
!                                                                               
!     /ITRA/ NSTP, Y(3), YOLD(3), SGN, DS                                       
!     NSTP      Step count. Incremented in subprogram LINAPX.                   
!     Y         Array containing current tracing point cartesian coordinates.   
!     YOLD      Array containing previous tracing point cartesian coordinates.
!     SGN       Determines direction of trace. Set in subprogram LINAPX         
!     DS        Step size (Km) Computed in subprogram LINAPX.                   
!                                                                               
!     /FLDCOMD/ BX, BY, BZ, BB                                                  
!     BX        X comp. of field vector at the current tracing point (Gauss)    
!     BY        Y comp. of field vector at the current tracing point (Gauss)    
!     BZ        Z comp. of field vector at the current tracing point (Gauss)    
!     BB        Magnitude of field vector at the current tracing point (Gauss)  
!                                                                               
!***REFERENCES  Stassinopoulos E. G. , Mead Gilbert D., X-841-72-17             
!                 (1971) GSFC, Greenbelt, Maryland                              
!                                                                               
!***ROUTINES CALLED  GD2CART,CONVRT,FELDG,
!                    ITRACE,DIPAPX,FNDAPX                                       
!                                                                               
!***END PROLOGUE  LINAPX                                                        
!-----------------------------------------------------------------------
!
      use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none
!
!-------------------------------Commons---------------------------------
!
      real(r8) bx, by, bz, bb
      COMMON /FLDCOMD/ BX, BY, BZ, BB

      real(r8) yapx
      COMMON /APXIN/   YAPX(3,3)

      real(r8) colat, elon, vp, ctp, stp
      COMMON /DIPOLE/  COLAT,ELON,VP,CTP,STP
!
! bf: separate ints and reals into 2 different commons:
!     COMMON /ITRA/    NSTP, Y(3), YP(3),  SGN, DS

      integer nstp
      COMMON /ITRA_int/    NSTP

      real(r8) y, yp, sqn, ds
      COMMON /ITRA/  Y(3), YP(3),  SGN, DS
!
!------------------------------Arguments--------------------------------
!
      real(r8) gdlat, glon, alt, a, alat, alon, xmag, ymag, zmag, f
!
!-----------------------------Parameters------------------------------
!
      integer maxs
      real(r8) rtod, dtor, re, req
      PARAMETER (MAXS = 200, &
                 RTOD=5.72957795130823E1_r8, DTOR=1.745329251994330E-2_r8, &
                 RE=6371.2_r8, REQ=6378.160_r8)
!
!---------------------------Local variables-----------------------------
!
      integer ientry, i, j

      real(r8) rho, bnrth, xlon, ht, beast, iapx, babs, bdown, xlat, r
      real(r8)  gclat, sgn, singml, cgml2
!
!-----------------------------------------------------------------------
!
!          Determine step size as a function of geomagnetic dipole
!          coordinates of the starting point
      CALL CONVRT (2,GDLAT,ALT,GCLAT,R)
!     SINGML = .98*SIN(GCLAT*DTOR) + .199*COS(GCLAT*DTOR)*
      SINGML = CTP*SIN(GCLAT*DTOR) +  STP*COS(GCLAT*DTOR)* &
                                                   COS((GLON-ELON)*DTOR)

!          May 99: avoid possible divide by zero (when SINGML = 1.)
      CGML2 = max (0.25_r8,1._r8-SINGML*SINGML)
      DS = .06_r8*R/CGML2 - 370._r8
!old:      Limit DS to its value at 60 deg gm latitude.
!old: DS = .06*R/(1.-SINGML*SINGML) - 370.
!old: IF (DS .GT. 1186.) DS = 1186.

!          Initialize YAPX array
      DO 4 J=1,3
      DO 4 I=1,3
    4 YAPX(I,J) = 0._r8

!          Convert from geodetic to earth centered cartesian coordinates
      CALL GD2CART (GDLAT,GLON,ALT,Y(1),Y(2),Y(3))
      NSTP = 0

!          Get magnetic field components to determine the direction for
!          tracing the field line
      CALL FELDG (1,GDLAT,GLON,ALT,XMAG,YMAG,ZMAG,F)
      SGN = SIGN (1._r8,-ZMAG)

!          Use cartesian coordinates to get magnetic field components
!          (from which gradients steer the tracing)
   10 IENTRY=2

      CALL FELDG (IENTRY,Y(1)/RE,Y(2)/RE,Y(3)/RE,BX,BY,BZ,BB)

      NSTP = NSTP + 1

!          Quit if too many steps
      IF (NSTP .GE. MAXS) THEN                                                 
        RHO = SQRT(Y(1)*Y(1) + Y(2)*Y(2))
        CALL CONVRT(3,XLAT,HT,RHO,Y(3))
        XLON = RTOD*ATAN2(Y(2),Y(1))
	CALL FELDG (1,XLAT,XLON,HT,BNRTH,BEAST,BDOWN,BABS)
	CALL DIPAPX (XLAT,XLON,HT,BNRTH,BEAST,BDOWN,A,ALON)
        ALAT = -SGN*RTOD*ACOS(SQRT(1._r8/A))
        RETURN                                                                
      END IF                                                                   

!          Find next point using adams algorithm after 7 points
       CALL ITRACE (IAPX)
       IF (IAPX .EQ. 1) GO TO 10

!          Maximum radius just passed.  Find apex coords
      CALL FNDAPX (ALT,ZMAG,A,ALAT,ALON)
      RETURN
      end subroutine linapx

!================================================================================================

      SUBROUTINE ITRACE (IAPX)
!
!-----------------------------------------------------------------------
!***BEGIN PROLOGUE  ITRACE                                                      
!***DATE WRITTEN   731029   (YYMMDD)                                            
!***AUTHOR  CLARK, W. N.O.A.A. ERL LAB.                                         
!***REVISION DATE  880201, H. Passi 
!***PURPOSE  Field line integration routine.                                    
!***DESCRIPTION                                                                 
!     It uses 4-point ADAMS formula after initialization.                       
!     First 7 iterations advance point by 3 steps.                              
!                                                                               
!   INPUT                                                                       
!            Passed in through the common blocks ITRA, FLDCOMD.                 
!            See the description below.                                         
!   OUTPUT                                                                      
!    IAPX    Flag set to 2 when APEX passed, otherwise set to 1.                
!                                                                               
!            Passed out through the common block APXIN.                         
!            See the description below.                                         
!                                                                               
!     COMMON Blocks Used                                                        
!     /APXIN/ YAPX(3,3)                                                         
!     YAPX      Matrix of cartesian coordinates (loaded columnwise)             
!               of the 3 points about APEX. Set in subprogram ITRACE.           
!                                                                               
!     /FLDCOMD/ BX, BY, BZ, BB                                                  
!     BX        X comp. of field vector at the current tracing point (Gauss)    
!     BY        Y comp. of field vector at the current tracing point (Gauss)    
!     BZ        Z comp. of field vector at the current tracing point (Gauss)    
!     BB        Magnitude of field vector at the current tracing point (Gauss)  
!                                                                               
!     /ITRA/ NSTP, Y(3), YOLD(3), SGN, DS                                       
!     NSTP      Step count for line integration.                                
!               Incremented in subprogram LINAPX.                               
!     Y         Array containing current tracing point cartesian coordinates.   
!     YOLD      Array containing previous tracing point cartesian coordinates.  
!     SGN       Determines direction of trace. Set in subprogram LINAPX         
!     DS        Integration step size (arc length Km)                           
!               Computed in subprogram LINAPX.                                  
!                                                                               
!***REFERENCES  reference 1                                                     
!                                                                               
!***ROUTINES CALLED  None                                                       
!***COMMON BLOCKS    APXIN,FLDCOMD,ITRA                                         
!***END PROLOGUE  ITRACE                                                        
! **                                                                            
!-----------------------------------------------------------------------
!
      use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none
!
!-------------------------------Commons---------------------------------
!
! bf: separate ints and reals into 2 different commons:
!     COMMON /ITRA/    NSTP, Y(3), YOLD(3), SGN, DS

      integer nstp
      COMMON /ITRA_int/    NSTP

      real(r8) y, yold, sgn, ds
      COMMON /ITRA/ Y(3), YOLD(3), SGN, DS

      real(r8) bx, by, bz, bb
      COMMON /FLDCOMD/ BX, BY, BZ, BB

      real(r8) yapx
      COMMON /APXIN/   YAPX(3,3)
!
!------------------------------Arguments--------------------------------
!
      real(r8) iapx
!
!---------------------------Local variables-----------------------------
!
      integer i, j

      real(r8)  YP(3, 4)
      real(r8) rc, d24, rp, term, d2, d12, d6

      SAVE
!
!-------------------------Statement Functions----------------------------
!
      real(r8) rdus, d, e, f
      RDUS(D,E,F) = SQRT( D**2 + E**2 + F**2 )
!
!-----------------------------------------------------------------------
!
      IAPX = 1
!          Field line is defined by the following differential equations
!          in cartesian coordinates:
      YP(1, 4) = SGN*BX/BB
      YP(2, 4) = SGN*BY/BB
      YP(3, 4) = SGN*BZ/BB
      IF (NSTP .GT. 7) GO TO 90

!          First seven steps use this block
      DO 80 I = 1, 3
      GO TO (10, 20, 30, 40, 50, 60, 70) NSTP
   10 D2  = DS/2._r8
      D6  = DS/6._r8
      D12 = DS/12._r8
      D24 = DS/24._r8
      YP(I,1)   = YP(I,4)
      YOLD(I)   = Y(I)
      YAPX(I,1) = Y(I)
      Y(I) = YOLD(I) + DS*YP(I, 1)
      GO TO 80

   20 YP(I, 2) = YP(I,4)
      Y(I) = YOLD(I) + D2*(YP(I,2)+YP(I,1))
      GO TO 80

   30 Y(I) = YOLD(I) + D6*(2._r8*YP(I,4)+YP(I,2)+3._r8*YP(I,1))
      GO TO 80

   40 YP(I,2)  = YP(I,4)
      YAPX(I,2)= Y(I)
      YOLD(I)  = Y(I)
      Y(I)     = YOLD(I) + D2*(3._r8*YP(I,2)-YP(I,1))
      GO TO 80

   50 Y(I) = YOLD(I) + D12*(5._r8*YP(I,4)+8._r8*YP(I,2)-YP(I,1))
      GO TO 80

   60 YP(I,3)  = YP(I,4)
      YOLD(I)  = Y(I)
      YAPX(I,3)= Y(I)
      Y(I)     = YOLD(I) + D12*(23._r8*YP(I,3)-16._r8*YP(I,2)+5._r8*YP(I,1))
      GO TO 80

   70 YAPX(I, 1) = YAPX(I, 2)
      YAPX(I, 2) = YAPX(I, 3)
      Y(I) = YOLD(I) + D24*(9._r8*YP(I,4)+19._r8*YP(I,3)-5._r8*YP(I,2)+YP(I,1))
      YAPX(I, 3) = Y(I)
   80 CONTINUE

!          Signal if apex passed
      IF ( NSTP .EQ. 6 .OR. NSTP .EQ. 7) THEN
	RC = RDUS( YAPX(1,3), YAPX(2,3), YAPX(3,3))
        RP = RDUS( YAPX(1,2), YAPX(2,2), YAPX(3,2))                             
	IF (RC .LT. RP) IAPX=2
      ENDIF
      RETURN

!          Stepping block for NSTP .gt. 7
   90 DO 110 I = 1, 3                                                           
      YAPX(I, 1) = YAPX(I, 2)
      YAPX(I, 2) = Y(I)
      YOLD(I) = Y(I)
      TERM = 55._r8*YP(I, 4) - 59._r8*YP(I, 3) + 37._r8*YP(I, 2) - 9._r8*YP(I, 1)
      Y(I) = YOLD(I) + D24*TERM
      YAPX(I, 3) = Y(I)

      DO 100 J = 1, 3
      YP(I, J) = YP(I, J+1)
  100 CONTINUE
  110 CONTINUE                                                                  

      RC = RDUS (   Y(1),    Y(2),    Y(3))
      RP = RDUS (YOLD(1), YOLD(2), YOLD(3))
      IF (RC .LT. RP) IAPX=2

      RETURN
      END subroutine itrace

!================================================================================================

       SUBROUTINE FNDAPX (ALT,ZMAG,A,ALAT,ALON)
!
!-----------------------------------------------------------------------
!***BEGIN PROLOGUE  FNDAPX                                                      
!***DATE WRITTEN   731023   (YYMMDD)                                            
!***AUTHOR  CLARK, W., NOAA BOULDER                                             
!***REVISION DATE  940803, A. D. Richmond, NCAR 
!***PURPOSE  Finds apex coords once tracing has signalled that the apex         
!            has been passed.  
!***DESCRIPTION                                                                 
!                                                                               
!     It uses second degree interpolation routine, FINT, to find                
!     apex latitude and apex longtitude.                                        
!   INPUT                                                                       
!     ALT    Altitude of starting point
!     ZMAG   Downward component of geomagnetic field at starting point
!     NMAX   Order of IGRF coefficients being used
!     G      Array of coefficients from COFRM
!   OUTPUT                                                                      
!          A  (Apex height + REQ)/REQ, where REQ = equatorial Earth radius.
!             A is analogous to the L value in invariant coordinates.
!       ALAT  Apex Lat. (deg)                                                   
!       ALON  Apex Lon. (deg)                                                   
!                                                                               
!***LONG DESCRIPTION                                                            
!                                                                               
!     COMMON Blocks Used                                                        
!     /APXIN/ YAPX(3,3)                                                         
!     YAPX      Matrix of cartesian coordinates (loaded columnwise)             
!               of the 3 points about APEX. Set in subprogram ITRACE.           
!                                                                               
!   /DIPOLE/COLAT,ELON,VP,CTP,STP                           
!     COLAT   Geocentric colatitude of the north pole of the                    
!             earth-centered dipole (Deg).                                      
!     ELON    Geographic longitude of the north pole of the                     
!             earth-centered dipole (Deg).                                      
!     CTP     cos(COLAT*DTOR)
!     STP     sin(COLAT*DTOR)
!                                                                               
!***ROUTINES CALLED  FINT                                                       
!***COMMON BLOCKS    APXIN,DIPOLE                                         
!***END PROLOGUE  FNDAPX                                                        
!-----------------------------------------------------------------------
!
      use shr_kind_mod, only: r8 => shr_kind_r8
      use abortutils,   only: endrun
      use cam_logfile,  only: iulog
      implicit none
!
!-------------------------------Commons---------------------------------
!
      real(r8) yapx
      COMMON /APXIN/  YAPX(3,3)

      real(r8) colat, elon, vp, ctp, stp
      COMMON /DIPOLE/ COLAT,ELON,VP,CTP,STP
!
!------------------------------Arguments--------------------------------
!
      real(r8) alt, zmag, a, alat, alon
!
!-----------------------------Parameters------------------------------
!
      real(r8) rtod, dtor, re, req
      PARAMETER (RTOD=5.72957795130823E1_r8,DTOR=1.745329251994330E-2_r8)
      PARAMETER (RE=6371.2_r8,REQ=6378.160_r8)
!
!---------------------------Local variables-----------------------------
!
      integer i, ientry

      real(r8) Z(3), HT(3), Y(3)
      real(r8) xlon, ang, cang, f, xinter, rasq, ste, stfcpa, stfspa
      real(r8) sang, r, cte, gdln, x, ydum, gdlt, rho
!
!-----------------------------------------------------------------------
!
!       ****
!       ****     GET GEODETIC FIELD COMPONENTS
!       ****
      IENTRY = 1
      DO 2 I = 1,3
	RHO = SQRT(YAPX(1,I)**2+YAPX(2,I)**2)
	CALL CONVRT (3,GDLT,HT(I),RHO,YAPX(3,I))
	GDLN = RTOD*ATAN2 (YAPX(2,I),YAPX(1,I))
	CALL FELDG (IENTRY,GDLT,GDLN,HT(I),X,YDUM,Z(I),F)
    2 CONTINUE
!     ****
!     ****     FIND CARTESIAN COORDINATES AT DIP EQUATOR BY INTERPOLATION
!     ****
      DO 3 I = 1,3
	CALL FINT(Z(1),Z(2),Z(3),YAPX(I,1),YAPX(I,2),YAPX(I,3),0._r8,Y(I))
    3 CONTINUE
!     ****
!     ****     FIND APEX HEIGHT BY INTERPOLATION
!     ****
      CALL FINT(Z(1),Z(2),Z(3),HT(1),HT(2),HT(3),0._r8,XINTER)
!          Ensure that XINTER is not less than original starting altitude (ALT)
      XINTER = max(ALT,XINTER)
      A = (REQ+XINTER)/(REQ)
!     ****
!     ****     FIND APEX COORDINATES , GIVING ALAT SIGN OF DIP AT
!     ****       STARTING POINT.  ALON IS THE VALUE OF THE GEOMAGNETIC
!     ****       LONGITUDE AT THE APEX.
!     ****
      IF (A.LT.1._r8) THEN
	write(iulog,20)
  20    FORMAT (' BOMBED! THIS MAKES A LESS THAN ONE')
	call endrun
      ENDIF
      RASQ = RTOD*ACOS(SQRT(1._r8/A))
      ALAT = SIGN(RASQ,ZMAG)

! Algorithm for ALON:
!   Use spherical coordinates.
!   Let GP be geographic pole.
!   Let GM be geomagnetic pole (colatitude COLAT, east longitude ELON).
!   Let XLON be longitude of apex.
!   Let TE be colatitude of apex.
!   Let ANG be longitude angle from GM to apex.
!   Let TP be colatitude of GM.
!   Let TF be arc length between GM and apex.
!   Let PA = ALON be geomagnetic longitude, i.e., Pi minus angle measured 
!     counterclockwise from arc GM-apex to arc GM-GP.
!   Then, using notation C=cos, S=sin, spherical-trigonometry formulas 
!     for the functions of the angles are as shown below.  Note: STFCPA,
!     STFSPA are sin(TF) times cos(PA), sin(PA), respectively.

      XLON = ATAN2(Y(2),Y(1))
      ANG  = XLON-ELON*DTOR
      CANG = COS(ANG)
      SANG = SIN(ANG)
      R    = SQRT(Y(1)**2+Y(2)**2+Y(3)**2)
      CTE  = Y(3)/R
      STE  = SQRT(1._r8-CTE*CTE)
      STFCPA = STE*CTP*CANG - CTE*STP
      STFSPA = SANG*STE
      ALON = ATAN2(STFSPA,STFCPA)*RTOD
      RETURN
      END subroutine fndapx 

!================================================================================================

      SUBROUTINE DIPAPX(GDLAT,GDLON,ALT,BNORTH,BEAST,BDOWN,A,ALON)
!
!-----------------------------------------------------------------------
! Compute a, alon from local magnetic field using dipole and spherical approx.
! 940501 A. D. Richmond
! Input:
!   GDLAT  = geodetic latitude, degrees
!   GDLON  = geodetic longitude, degrees
!   ALT    = altitude, km
!   BNORTH = geodetic northward magnetic field component (any units)
!   BEAST  = eastward magnetic field component
!   BDOWN  = geodetic downward magnetic field component
! Output:
!   A      = apex radius, 1 + h_A/R_eq
!   ALON   = apex longitude, degrees
!
! Algorithm:
!   Use spherical coordinates.
!   Let GP be geographic pole.
!   Let GM be geomagnetic pole (colatitude COLAT, east longitude ELON).
!   Let G be point at GDLAT,GDLON.
!   Let E be point on sphere below apex of dipolar field line passing through G.
!   Let TD be dipole colatitude of point G, found by applying dipole formula
!     for dip angle to actual dip angle.
!   Let B be Pi plus local declination angle.  B is in the direction 
!     from G to E.
!   Let TG be colatitude of G.
!   Let ANG be longitude angle from GM to G.
!   Let TE be colatitude of E.
!   Let TP be colatitude of GM.
!   Let A be longitude angle from G to E.
!   Let APANG = A + ANG
!   Let PA be geomagnetic longitude, i.e., Pi minus angle measured 
!     counterclockwise from arc GM-E to arc GM-GP.
!   Let TF be arc length between GM and E.
!   Then, using notation C=cos, S=sin, COT=cot, spherical-trigonometry formulas 
!     for the functions of the angles are as shown below.  Note: STFCPA,
!     STFSPA are sin(TF) times cos(PA), sin(PA), respectively.
!-----------------------------------------------------------------------
!
      use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none
!
!-------------------------------Commons---------------------------------
!
      real(r8) colat, elon, vp, ctp, stp
      COMMON/DIPOLE/COLAT,ELON,VP,CTP,STP
!
!------------------------------Arguments--------------------------------
!
      real(r8) gdlat, gdlon, alt, bnorth, beast, bdown, a, alon
!
!-----------------------------Parameters------------------------------
!
      real(r8) rtod, dtor, re, req
      PARAMETER (RTOD=5.72957795130823E1_r8,DTOR=1.745329251994330E-2_r8)
      PARAMETER (RE=6371.2_r8,REQ=6378.160_r8)
!
!---------------------------Local variables-----------------------------
!
      real(r8) bhor, ca, capang, sa, cte, ste, r, ha, stfspa, sapang
      real(r8) stfcpa, sb, cb, ctd, cottd, std, sang, cang, ang, ctg
      real(r8) stg
!
!-----------------------------------------------------------------------
!
      BHOR = SQRT(BNORTH*BNORTH + BEAST*BEAST)
      IF (BHOR.EQ.0._r8) THEN
	ALON = 0._r8
	A = 1.E34_r8
	RETURN
      ENDIF
      COTTD = BDOWN*.5_r8/BHOR
      STD = 1._r8/SQRT(1._r8+COTTD*COTTD)
      CTD = COTTD*STD
      SB = -BEAST/BHOR
      CB = -BNORTH/BHOR
      CTG = SIN(GDLAT*DTOR) 
      STG = COS(GDLAT*DTOR)
      ANG = (GDLON-ELON)*DTOR 
      SANG = SIN(ANG)
      CANG = COS(ANG)
      CTE = CTG*STD + STG*CTD*CB
      STE = SQRT(1._r8 - CTE*CTE)
      SA = SB*CTD/STE
      CA = (STD*STG - CTD*CTG*CB)/STE
      CAPANG = CA*CANG - SA*SANG
      SAPANG = CA*SANG + SA*CANG
      STFCPA = STE*CTP*CAPANG - CTE*STP
      STFSPA = SAPANG*STE
      ALON = ATAN2(STFSPA,STFCPA)*RTOD
      R = ALT + RE
      HA = ALT + R*COTTD*COTTD
      A = 1._r8 + HA/REQ
      RETURN
      end subroutine dipapx

!================================================================================================

       SUBROUTINE FINT (A1, A2, A3, A4, A5, A6, A7, RESULT)
!
!-----------------------------------------------------------------------
!***PURPOSE  Second degree interpolation routine                                
!***REFER TO  FNDAPX                                                            
!-----------------------------------------------------------------------
!
      use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none
!
!------------------------------Arguments--------------------------------
!
      real(r8) a1, a2, a3, a4, a5, a6, a7, result

       RESULT = ((A2-A3)*(A7-A2)*(A7-A3)*A4-(A1-A3)*(A7-A1)*(A7-A3)*A5+          &
        (A1-A2)*(A7-A1)*(A7-A2)*A6)/((A1-A2)*(A1-A3)*(A2-A3))
       RETURN                                                                   
       END subroutine fint 

!================================================================================================

      SUBROUTINE COFRM (DATE)
!
!-----------------------------------------------------------------------
!
! Original file ~bozo/pgms/apex/magfld.f copied 2/25/00.
!
!          Assign DGRF/IGRF spherical harmonic coefficients, to degree and
!          order NMAX, for DATE, yyyy.fraction, into array G.  Coefficients
!          are interpolated from the DGRF dates through the current IGRF year.
!          Coefficients for a later DATE are extrapolated using the IGRF
!          initial value and the secular change coefficients.  A warning
!          message is issued if DATE is later than the last recommended
!          (5 yrs later than the IGRF).  An DATE input earlier than the
!          first DGRF (EPOCH(1)), results in a diagnostic and a STOP.
!
!          Output in COMMON /MAGCOF/ NMAX,GB(144),GV(144),ICHG
!             NMAX = Maximum order of spherical harmonic coefficients used
!             GB   = Coefficients for magnetic field calculation
!             GV   = Coefficients for magnetic potential calculation
!             ICHG = Flag indicating when GB,GV have been changed in COFRM
!
!          HISTORY (blame):
!          COFRM and FELDG originated 15 Apr 83 by Vincent B. Wickwar
!          (formerly at SRI. Int., currently at Utah State).  Although set
!          up to accomodate second order time derivitives, the IGRF
!          (GTT, HTT) have been zero.  The spherical harmonic coefficients
!          degree and order is defined by NMAX (currently 10).
!
!          Jun 86:  Updated coefficients adding DGRF 1980 & IGRF 1985, which
!          were obtained from Eos Vol. 7, No. 24.  Common block MAG was
!          replaced by MAGCOF, thus removing variables not used in subroutine
!          FELDG.  (Roy Barnes)
!
!          Apr 1992 (Barnes):  Added DGRF 1985 and IGRF 1990 as described
!          in EOS Vol 73 Number 16 Apr 21 1992.  Other changes were made so
!          future updates should:
!            (1) Increment NDGY;
!            (2) Append to EPOCH the next IGRF year;
!            (3) Append the next DGRF coefficients to G1DIM and H1DIM; and
!            (4) Replace the IGRF initial values (G0, GT) and rates of
!                change indices (H0, HT).
!
!          Apr 94 (Art Richmond): Computation of GV added, for finding
!          magnetic potential.
!
!          Aug 95 (Barnes):  Added DGRF for 1990 and IGRF for 1995, which were
!          obtained by anonymous ftp geomag.gsfc.nasa.gov (cd pub, mget table*)
!          as per instructions from Bob Langel (langel@geomag.gsfc.nasa.gov),
!          but, problems are to be reported to baldwin@geomag.gsfc.nasa.gov
 
!          Oct 95 (Barnes):  Correct error in IGRF-95 G 7 6 and H 8 7 (see
!          email in folder).  Also found bug whereby coefficients were not being
!          updated in FELDG when IENTY did not change. ICHG was added to flag
!          date changes.  Also, a vestigial switch (IS) was removed from COFRM:
!          It was always 0 and involved 3 branch if statements in the main
!          polynomial construction loop (now numbered 200).
! 
!          Feb 99 (Barnes):  Explicitly initialize GV(1) in COFRM to avoid
!          possibility of compiler or loader options initializing memory
!          to something else (e.g., indefinite).  Also simplify the algebra
!          in COFRM; this does not effect results.
!
!          Mar 99 (Barnes):  Removed three branch if's from FELDG and changed
!          statement labels to ascending order
!
!          Jun 99 (Barnes):  Corrected RTOD definition in GD2CART.
!
!          May 00 (Barnes):  Replace IGRF 1995 with DGRF 1995, add IGRF
!          2000, and extend the earlier DGRF's backward to 1900.  A complete
!          set of coefficients came from a NGDC web page
!
!
!-----------------------------------------------------------------------
!
      use shr_kind_mod,  only: r8 => shr_kind_r8
      use abortutils,    only: endrun
      use cam_logfile,   only: iulog
      implicit none

#ifdef NO_R16
   integer,parameter :: r16= selected_real_kind(12) ! 8 byte real
#else
   integer,parameter :: r16= selected_real_kind(24) ! 16 byte real
#endif

!
!-------------------------------Commons---------------------------------
!
!     COMMON /MAGCOF/ NMAX,G(144),GV(144),ICHG

      integer nmax, ichg
      COMMON /MAGCOF_int/ NMAX,ICHG

      real(r8) g, gv
      COMMON /MAGCOF_real/ G(144),GV(144)
!     DATA NMAX,ICHG /10,-99999/
!
!------------------------------Arguments--------------------------------
!
      real(r8) date
!
!-----------------------------Parameters------------------------------
!
      real(r8) rtod, dtor
      PARAMETER (RTOD=5.72957795130823E1_r8,DTOR=1.745329251994330E-2_r8)

!      PARAMETER (NDGY=6 , NYT = NDGY+1 , NGH = 144*NDGY)

      integer ndgy, nyt, ngh
      PARAMETER (NDGY=20 , NYT = NDGY+1 , NGH = 144*NDGY)
!          NDGY = Number of DGRF years of sets of coefficients
!          NYT  = Add one for the IGRF set (and point to it).
!          NGH  = Dimension of the equivalenced arrays
!
!---------------------------Local variables-----------------------------
!
      integer i, mm, nn, m, i1, n, iy
 
      real(r8) a1, a2, a3, a4, a5, a6, a7, result
      real(r8) rnn, datel, t, time
      real(r16) F,F0

      REAL(r8) GYR(12,12,NYT) , HYR(12,12,NYT), EPOCH(NYT) , &
                G1DIM(NGH)     , H1DIM(NGH) , &
                G0(12,12) , GT(12,12) , GTT(12,12) , &
                H0(12,12) , HT(12,12) ,  HTT(12,12)
      EQUIVALENCE (GYR(1,1,1),G1DIM(1))  , (HYR(1,1,1),H1DIM(1)) , &
                  (GYR(1,1,NYT),G0(1,1)) , (HYR(1,1,NYT),H0(1,1))

      SAVE DATEL,GYR,HYR,EPOCH,G1DIM,H1DIM,G0,H0,GT,HT,GTT,HTT
!     SAVE DATEL,GYR,HYR,EPOCH,            G0,H0,GT,HT,GTT,HTT

      DATA DATEL /-999._r8/
!      DATA EPOCH /1965. , 1970. , 1975. , 1980. , 1985. , 1990. , 1995./
      DATA EPOCH /1900, 1905, 1910, 1915, 1920, 1925, 1930, 1935, 1940, &
      1945, 1950, 1955, 1960, 1965, 1970, 1975, 1980, 1985, 1990, 1995, &
      2000/
!          D_/Dtime2 coefficients are 0
      DATA GTT/144*0._r8/,HTT/144*0._r8/
!          DGRF g(n,m) for 1965:
!          The "column" corresponds to "n" and
!          the "line" corresponds to "m" as indicated in column 6;
!          e.g., for 1965 g(0,3) = 1297. or g(6,6) = -111.
      DATA (G1DIM(I),I=1,144) /0, &
        -31543,  -677,  1022,  876, -184,   63,  70,  11,   8, -3,  2*0, &
         -2298,  2905, -1469,  628,  328,   61, -55,   8,  10, -4,  3*0, &
                  924,  1256,  660,  264,  -11,   0,  -4,   1,  2,  4*0, &
                         572, -361,    5, -217,  34,  -9, -11, -5,  5*0, &
                               134,  -86,  -58, -41,   1,  12, -2,  6*0, &
                                     -16,   59, -21,   2,   1,  6,  7*0, &
                                           -90,  18,  -9,  -2,  4,  8*0, &
                                                  6,   5,   2,  0,  9*0, &
                                                       8,  -1,  2, 10*0, &
                                                           -1,  2, 11*0, &
                                                                0, 13*0/
!          DGRF g(n,m) for 1905:
      DATA (G1DIM(I),I=145,288) /0, &
        -31464,  -728,  1037,  880, -192,   62,  70,  11,   8, -3,  2*0, &
         -2298,  2928, -1494,  643,  328,   60, -54,   8,  10, -4,  3*0, &
                 1041,  1239,  653,  259,  -11,   0,  -4,   1,  2,  4*0, &
                         635, -380,   -1, -221,  33,  -9, -11, -5,  5*0, &
                               146,  -93,  -57, -41,   1,  12, -2,  6*0, &
                                     -26,   57, -20,   2,   1,  6,  7*0, &
                                           -92,  18,  -8,  -2,  4,  8*0, &
                                                  6,   5,   2,  0,  9*0, &
                                                       8,   0,  2, 10*0, &
                                                           -1,  2, 11*0, &
                                                                0, 13*0/
!          DGRF g(n,m) for 1910:
      DATA (G1DIM(I),I=289,432) /0, &
        -31354,  -769,  1058,  884, -201,   62,  71,  11,   8, -3,  2*0, &
         -2297,  2948, -1524,  660,  327,   58, -54,   8,  10, -4,  3*0, &
                 1176,  1223,  644,  253,  -11,   1,  -4,   1,  2,  4*0, &
                         705, -400,   -9, -224,  32,  -9, -11, -5,  5*0, &
                               160, -102,  -54, -40,   1,  12, -2,  6*0, &
                                     -38,   54, -19,   2,   1,  6,  7*0, &
                                           -95,  18,  -8,  -2,  4,  8*0, &
                                                  6,   5,   2,  0,  9*0, &
                                                       8,   0,  2, 10*0, &
                                                           -1,  2, 11*0, &
                                                                0, 13*0/
!          DGRF g(n,m) for 1915:
      DATA (G1DIM(I),I=433,576) /0, &
        -31212,  -802,  1084,  887, -211,   61,  72,  11,   8, -3,  2*0, &
         -2306,  2956, -1559,  678,  327,   57, -54,   8,  10, -4,  3*0, &
                 1309,  1212,  631,  245,  -10,   2,  -4,   1,  2,  4*0, &
                         778, -416,  -16, -228,  31,  -9, -11, -5,  5*0, &
                               178, -111,  -51, -38,   2,  12, -2,  6*0, &
                                     -51,   49, -18,   3,   1,  6,  7*0, &
                                           -98,  19,  -8,  -2,  4,  8*0, &
                                                  6,   6,   2,  0,  9*0, &
                                                       8,   0,  1, 10*0, &
                                                           -1,  2, 11*0, &
                                                                0, 13*0/
!          DGRF g(n,m) for 1920:
      DATA (G1DIM(I),I=577,720) /0, &
        -31060,  -839,  1111,  889, -221,   61,  73,  11,   8, -3,  2*0, &
         -2317,  2959, -1600,  695,  326,   55, -54,   7,  10, -4,  3*0, &
                 1407,  1205,  616,  236,  -10,   2,  -3,   1,  2,  4*0, &
                         839, -424,  -23, -233,  29,  -9, -11, -5,  5*0, &
                               199, -119,  -46, -37,   2,  12, -2,  6*0, &
                                     -62,   44, -16,   4,   1,  6,  7*0, &
                                          -101,  19,  -7,  -2,  4,  8*0, &
                                                  6,   6,   2,  0,  9*0, &
                                                       8,   0,  1, 10*0, &
                                                           -1,  3, 11*0, &
                                                                0, 13*0/
!          DGRF g(n,m) for 1925:
      DATA (G1DIM(I),I=721,864) /0, &
        -30926,  -893,  1140,  891, -230,   61,  73,  11,   8, -3,  2*0, &
         -2318,  2969, -1645,  711,  326,   54, -54,   7,  10, -4,  3*0, &
                 1471,  1202,  601,  226,   -9,   3,  -3,   1,  2,  4*0, &
                         881, -426,  -28, -238,  27,  -9, -11, -5,  5*0, &
                               217, -125,  -40, -35,   2,  12, -2,  6*0, &
                                     -69,   39, -14,   4,   1,  6,  7*0, &
                                          -103,  19,  -7,  -2,  4,  8*0, &
                                                  6,   7,   2,  0,  9*0, &
                                                       8,   0,  1, 10*0, &
                                                           -1,  3, 11*0, &
                                                                0, 13*0/
!          DGRF g(n,m) for 1930:
      DATA (G1DIM(I),I=865,1008) /0, &
        -30805,  -951,  1172,  896, -237,   60,  74,  11,   8, -3,  2*0, &
         -2316,  2980, -1692,  727,  327,   53, -54,   7,  10, -4,  3*0, &
                 1517,  1205,  584,  218,   -9,   4,  -3,   1,  2,  4*0, &
                         907, -422,  -32, -242,  25,  -9, -12, -5,  5*0, &
                               234, -131,  -32, -34,   2,  12, -2,  6*0, &
                                     -74,   32, -12,   5,   1,  6,  7*0, &
                                          -104,  18,  -6,  -2,  4,  8*0, &
                                                  6,   8,   3,  0,  9*0, &
                                                       8,   0,  1, 10*0, &
                                                           -2,  3, 11*0, &
                                                                0, 13*0/
!          DGRF g(n,m) for 1935:
      DATA (G1DIM(I),I=1009,1152) /0, &
        -30715, -1018,  1206,  903, -241,   59,  74,  11,   8, -3,  2*0, &
         -2306,  2984, -1740,  744,  329,   53, -53,   7,  10, -4,  3*0, &
                 1550,  1215,  565,  211,   -8,   4,  -3,   1,  2,  4*0, &
                         918, -415,  -33, -246,  23,  -9, -12, -5,  5*0, &
                               249, -136,  -25, -33,   1,  11, -2,  6*0, &
                                     -76,   25, -11,   6,   1,  6,  7*0, &
                                          -106,  18,  -6,  -2,  4,  8*0, &
                                                  6,   8,   3,  0,  9*0, &
                                                       7,   0,  2, 10*0, &
                                                           -2,  3, 11*0, &
                                                                0, 13*0/
!          DGRF g(n,m) for 1940:
      DATA (G1DIM(I),I=1153,1296) /0, &
        -30654, -1106,  1240,  914, -241,   57,  74,  11,   8, -3,  2*0, &
         -2292,  2981, -1790,  762,  334,   54, -53,   7,  10, -4,  3*0, &
                 1566,  1232,  550,  208,   -7,   4,  -3,   1,  2,  4*0, &
                         916, -405,  -33, -249,  20, -10, -12, -5,  5*0, &
                               265, -141,  -18, -31,   1,  11, -2,  6*0, &
                                     -76,   18,  -9,   6,   1,  6,  7*0, &
                                          -107,  17,  -5,  -2,  4,  8*0, &
                                                  5,   9,   3,  0,  9*0, &
                                                       7,   1,  2, 10*0, &
                                                           -2,  3, 11*0, &
                                                                0, 13*0/
!          DGRF g(n,m) for 1945:
      DATA (G1DIM(I),I=1297,1440) /0, &
        -30594, -1244,  1282,  944, -253,   59,  70,  13,   5, -3,  2*0, &
         -2285,  2990, -1834,  776,  346,   57, -40,   7, -21, 11,  3*0, &
                 1578,  1255,  544,  194,    6,   0,  -8,   1,  1,  4*0, &
                         913, -421,  -20, -246,   0,  -5, -11,  2,  5*0, &
                               304, -142,  -25, -29,   9,   3, -5,  6*0, &
                                     -82,   21, -10,   7,  16, -1,  7*0, &
                                          -104,  15, -10,  -3,  8,  8*0, &
                                                 29,   7,  -4, -1,  9*0, &
                                                       2,  -3, -3, 10*0, &
                                                           -4,  5, 11*0, &
                                                               -2, 13*0/
!          DGRF g(n,m) for 1950:
      DATA (G1DIM(I),I=1441,1584) /0, &
        -30554, -1341,  1297,  954, -240,   54,  65,  22,   3, -8,  2*0, &
         -2250,  2998, -1889,  792,  349,   57, -55,  15,  -7,  4,  3*0, &
                 1576,  1274,  528,  211,    4,   2,  -4,  -1, -1,  4*0, &
                         896, -408,  -20, -247,   1,  -1, -25, 13,  5*0, &
                               303, -147,  -16, -40,  11,  10, -4,  6*0, &
                                     -76,   12,  -7,  15,   5,  4,  7*0, &
                                          -105,   5, -13,  -5, 12,  8*0, &
                                                 19,   5,  -2,  3,  9*0, &
                                                      -1,   3,  2, 10*0, &
                                                            8, 10, 11*0, &
                                                                3, 13*0/
!          DGRF g(n,m) for 1955:
      DATA (G1DIM(I),I=1585,1728) /0, &
        -30500, -1440,  1302,  958, -229,   47,  65,  11,   4, -3,  2*0, &
         -2215,  3003, -1944,  796,  360,   57, -56,   9,   9, -5,  3*0, &
                 1581,  1288,  510,  230,    3,   2,  -6,  -4, -1,  4*0, &
                         882, -397,  -23, -247,  10, -14,  -5,  2,  5*0, &
                               290, -152,   -8, -32,   6,   2, -3,  6*0, &
                                     -69,    7, -11,  10,   4,  7,  7*0, &
                                          -107,   9,  -7,   1,  4,  8*0, &
                                                 18,   6,   2, -2,  9*0, &
                                                       9,   2,  6, 10*0, &
                                                            5, -2, 11*0, &
                                                                0, 13*0/
!          DGRF g(n,m) for 1960:
      DATA (G1DIM(I),I=1729,1872) /0, &
        -30421, -1555,  1302,  957, -222,   46,  67,  15,   4,  1,  2*0, &
         -2169,  3002, -1992,  800,  362,   58, -56,   6,   6, -3,  3*0, &
                 1590,  1289,  504,  242,    1,   5,  -4,   0,  4,  4*0, &
                         878, -394,  -26, -237,  15, -11,  -9,  0,  5*0, &
                               269, -156,   -1, -32,   2,   1, -1,  6*0, &
                                     -63,   -2,  -7,  10,   4,  4,  7*0, &
                                          -113,  17,  -5,  -1,  6,  8*0, &
                                                  8,  10,  -2,  1,  9*0, &
                                                       8,   3, -1, 10*0, &
                                                           -1,  2, 11*0, &
                                                                0, 13*0/
!          DGRF g(n,m) for 1965:
      DATA (G1DIM(I),I=1873,2016) /0, &
        -30334, -1662,  1297,  957, -219,   45,  75,  13,   8, -2,  2*0, &
         -2119,  2997, -2038,  804,  358,   61, -57,   5,  10, -3,  3*0, &
                 1594,  1292,  479,  254,    8,   4,  -4,   2,  2,  4*0, &
                         856, -390,  -31, -228,  13, -14, -13, -5,  5*0, &
                               252, -157,    4, -26,   0,  10, -2,  6*0, &
                                     -62,    1,  -6,   8,  -1,  4,  7*0, &
                                          -111,  13,  -1,  -1,  4,  8*0, &
                                                  1,  11,   5,  0,  9*0, &
                                                       4,   1,  2, 10*0, &
                                                           -2,  2, 11*0, &
                                                                0, 13*0/
!          DGRF g(n,m) for 1970:
      DATA (G1DIM(I),I=2017,2160) /0, &
        -30220, -1781,  1287,  952, -216,   43,  72,  14,   8, -3,  2*0, &
         -2068,  3000, -2091,  800,  359,   64, -57,   6,  10, -3,  3*0, &
                 1611,  1278,  461,  262,   15,   1,  -2,   2,  2,  4*0, &
                         838, -395,  -42, -212,  14, -13, -12, -5,  5*0, &
                               234, -160,    2, -22,  -3,  10, -1,  6*0, &
                                     -56,    3,  -2,   5,  -1,  6,  7*0, &
                                          -112,  13,   0,   0,  4,  8*0, &
                                                 -2,  11,   3,  1,  9*0, &
                                                       3,   1,  0, 10*0, &
                                                           -1,  3, 11*0, &
                                                               -1, 13*0/
!          DGRF g(n,m) for 1975:
      DATA (G1DIM(I),I=2161,2304) /0, &
        -30100, -1902,  1276,  946, -218,   45,  71,  14,   7, -3,  2*0, &
         -2013,  3010, -2144,  791,  356,   66, -56,   6,  10, -3,  3*0, &
                 1632,  1260,  438,  264,   28,   1,  -1,   2,  2,  4*0, &
                         830, -405,  -59, -198,  16, -12, -12, -5,  5*0, &
                               216, -159,    1, -14,  -8,  10, -2,  6*0, &
                                     -49,    6,   0,   4,  -1,  5,  7*0, &
                                          -111,  12,   0,  -1,  4,  8*0, &
                                                 -5,  10,   4,  1,  9*0, &
                                                       1,   1,  0, 10*0, &
                                                           -2,  3, 11*0, &
                                                               -1, 13*0/
!          DGRF g(n,m) for 1980:
      DATA (G1DIM(I),I=2305,2448) /0, &
        -29992, -1997,  1281,  938, -218,   48,  72,  18,   5, -4,  2*0, &
         -1956,  3027, -2180,  782,  357,   66, -59,   6,  10, -4,  3*0, &
                 1663,  1251,  398,  261,   42,   2,   0,   1,  2,  4*0, &
                         833, -419,  -74, -192,  21, -11, -12, -5,  5*0, &
                               199, -162,    4, -12,  -7,   9, -2,  6*0, &
                                     -48,   14,   1,   4,  -3,  5,  7*0, &
                                          -108,  11,   3,  -1,  3,  8*0, &
                                                 -2,   6,   7,  1,  9*0, &
                                                      -1,   2,  2, 10*0, &
                                                           -5,  3, 11*0, &
                                                                0, 13*0/
!          DGRF g(n,m) for 1985:
      DATA (G1DIM(I),I=2449,2592) /0, &
        -29873, -2072,  1296,  936, -214,   53,  74,  21,   5, -4,  2*0, &
         -1905,  3044, -2208,  780,  355,   65, -62,   6,  10, -4,  3*0, &
                 1687,  1247,  361,  253,   51,   3,   0,   1,  3,  4*0, &
                         829, -424,  -93, -185,  24, -11, -12, -5,  5*0, &
                               170, -164,    4,  -6,  -9,   9, -2,  6*0, &
                                     -46,   16,   4,   4,  -3,  5,  7*0, &
                                          -102,  10,   4,  -1,  3,  8*0, &
                                                  0,   4,   7,  1,  9*0, &
                                                      -4,   1,  2, 10*0, &
                                                           -5,  3, 11*0, &
                                                                0, 13*0/
!          DGRF g(n,m) for 1990:
      DATA (G1DIM(I),I=2593,2736) /0, &
        -29775, -2131,  1314,  939, -214,   61,  77,  23,   4, -3,  2*0, &
         -1848,  3059, -2239,  780,  353,   65, -64,   5,   9, -4,  3*0, &
                 1686,  1248,  325,  245,   59,   2,  -1,   1,  2,  4*0, &
                         802, -423, -109, -178,  26, -10, -12, -5,  5*0, &
                               141, -165,    3,  -1, -12,   9, -2,  6*0, &
                                     -36,   18,   5,   3,  -4,  4,  7*0, &
                                           -96,   9,   4,  -2,  3,  8*0, &
                                                  0,   2,   7,  1,  9*0, &
                                                      -6,   1,  3, 10*0, &
                                                           -6,  3, 11*0, &
                                                                0, 13*0/
!          DGRF g(n,m) for 1995:
      DATA (G1DIM(I),I=2737,2880) /0, &
        -29682, -2197,  1329,  941, -210,   66,  78,  24,   4, -3,  2*0, &
         -1789,  3074, -2268,  782,  352,   64, -67,   4,   9, -4,  3*0, &
                 1685,  1249,  291,  237,   65,   1,  -1,   1,  2,  4*0, &
                         769, -421, -122, -172,  29,  -9, -12, -5,  5*0, &
                               116, -167,    2,   4, -14,   9, -2,  6*0, &
                                     -26,   17,   8,   4,  -4,  4,  7*0, &
                                           -94,  10,   5,  -2,  3,  8*0, &
                                                 -2,   0,   7,  1,  9*0, &
                                                      -7,   0,  3, 10*0, &
                                                           -6,  3, 11*0, &
                                                                0, 13*0/
!          DGRF h(n,m) for 1900:
      DATA (H1DIM(I),I=1,144) /13*0, &
          5922, -1061,  -330,  195, -210,   -9, -45,   8, -20,  2,  3*0, &
                 1121,     3,  -69,   53,   83, -13, -14,  14,  1,  4*0, &
                         523, -210,  -33,    2, -10,   7,   5,  2,  5*0, &
                               -75, -124,  -35,  -1, -13,  -3,  6,  6*0, &
                                       3,   36,  28,   5,  -2, -4,  7*0, &
                                           -69, -12,  16,   8,  0,  8*0, &
                                                -22,  -5,  10, -2,  9*0, &
                                                     -18,  -2,  4, 10*0, &
                                                            2,  0, 11*0, &
                                                               -6, 13*0/
!          DGRF h(n,m) for 1905:
      DATA (H1DIM(I),I=145,288) /13*0, &
          5909, -1086,  -357,  203, -193,   -7, -46,   8, -20,  2,  3*0, &
                 1065,    34,  -77,   56,   86, -14, -15,  14,  1,  4*0, &
                         480, -201,  -32,    4, -11,   7,   5,  2,  5*0, &
                               -65, -125,  -32,   0, -13,  -3,  6,  6*0, &
                                      11,   32,  28,   5,  -2, -4,  7*0, &
                                           -67, -12,  16,   8,  0,  8*0, &
                                                -22,  -5,  10, -2,  9*0, &
                                                     -18,  -2,  4, 10*0, &
                                                            2,  0, 11*0, &
                                                               -6, 13*0/
!          DGRF h(n,m) for 1910:
      DATA (H1DIM(I),I=289,432) /13*0, &
          5898, -1128,  -389,  211, -172,   -5, -47,   8, -20,  2,  3*0, &
                 1000,    62,  -90,   57,   89, -14, -15,  14,  1,  4*0, &
                         425, -189,  -33,    5, -12,   6,   5,  2,  5*0, &
                               -55, -126,  -29,   1, -13,  -3,  6,  6*0, &
                                      21,   28,  28,   5,  -2, -4,  7*0, &
                                           -65, -13,  16,   8,  0,  8*0, &
                                                -22,  -5,  10, -2,  9*0, &
                                                     -18,  -2,  4, 10*0, &
                                                            2,  0, 11*0, &
                                                               -6, 13*0/
!          DGRF h(n,m) for 1915:
      DATA (H1DIM(I),I=433,576) /13*0, &
          5875, -1191,  -421,  218, -148,   -2, -48,   8, -20,  2,  3*0, &
                  917,    84, -109,   58,   93, -14, -15,  14,  1,  4*0, &
                         360, -173,  -34,    8, -12,   6,   5,  2,  5*0, &
                               -51, -126,  -26,   2, -13,  -3,  6,  6*0, &
                                      32,   23,  28,   5,  -2, -4,  7*0, &
                                           -62, -15,  16,   8,  0,  8*0, &
                                                -22,  -5,  10, -2,  9*0, &
                                                     -18,  -2,  4, 10*0, &
                                                            2,  0, 11*0, &
                                                               -6, 13*0/
!          DGRF h(n,m) for 1920:
      DATA (H1DIM(I),I=577,720) /13*0, &
          5845, -1259,  -445,  220, -122,    0, -49,   8, -20,  2,  3*0, &
                  823,   103, -134,   58,   96, -14, -15,  14,  1,  4*0, &
                         293, -153,  -38,   11, -13,   6,   5,  2,  5*0, &
                               -57, -125,  -22,   4, -14,  -3,  6,  6*0, &
                                      43,   18,  28,   5,  -2, -4,  7*0, &
                                           -57, -16,  17,   9,  0,  8*0, &
                                                -22,  -5,  10, -2,  9*0, &
                                                     -19,  -2,  4, 10*0, &
                                                            2,  0, 11*0, &
                                                               -6, 13*0/
!          DGRF h(n,m) for 1925:
      DATA (H1DIM(I),I=721,864) /13*0, &
          5817, -1334,  -462,  216,  -96,    3, -50,   8, -20,  2,  3*0, &
                  728,   119, -163,   58,   99, -14, -15,  14,  1,  4*0, &
                         229, -130,  -44,   14, -14,   6,   5,  2,  5*0, &
                               -70, -122,  -18,   5, -14,  -3,  6,  6*0, &
                                      51,   13,  29,   5,  -2, -4,  7*0, &
                                           -52, -17,  17,   9,  0,  8*0, &
                                                -21,  -5,  10, -2,  9*0, &
                                                     -19,  -2,  4, 10*0, &
                                                            2,  0, 11*0, &
                                                               -6, 13*0/
!          DGRF h(n,m) for 1930:
      DATA (H1DIM(I),I=865,1008) /13*0, &
          5808, -1424,  -480,  205,  -72,    4, -51,   8, -20,  2,  3*0, &
                  644,   133, -195,   60,  102, -15, -15,  14,  1,  4*0, &
                         166, -109,  -53,   19, -14,   5,   5,  2,  5*0, &
                               -90, -118,  -16,   6, -14,  -3,  6,  6*0, &
                                      58,    8,  29,   5,  -2, -4,  7*0, &
                                           -46, -18,  18,   9,  0,  8*0, &
                                                -20,  -5,  10, -2,  9*0, &
                                                     -19,  -2,  4, 10*0, &
                                                            2,  0, 11*0, &
                                                               -6, 13*0/
!          DGRF h(n,m) for 1935:
      DATA (H1DIM(I),I=1009,1152) /13*0, &
          5812, -1520,  -494,  188,  -51,    4, -52,   8, -20,  2,  3*0, &
                  586,   146, -226,   64,  104, -17, -15,  15,  1,  4*0, &
                         101,  -90,  -64,   25, -14,   5,   5,  2,  5*0, &
                              -114, -115,  -15,   7, -15,  -3,  6,  6*0, &
                                      64,    4,  29,   5,  -3, -4,  7*0, &
                                           -40, -19,  18,   9,  0,  8*0, &
                                                -19,  -5,  11, -1,  9*0, &
                                                     -19,  -2,  4, 10*0, &
                                                            2,  0, 11*0, &
                                                               -6, 13*0/
!          DGRF h(n,m) for 1940:
      DATA (H1DIM(I),I=1153,1296) /13*0, &
          5821, -1614,  -499,  169,  -33,    4, -52,   8, -21,  2,  3*0, &
                  528,   163, -252,   71,  105, -18, -14,  15,  1,  4*0, &
                          43,  -72,  -75,   33, -14,   5,   5,  2,  5*0, &
                              -141, -113,  -15,   7, -15,  -3,  6,  6*0, &
                                      69,    0,  29,   5,  -3, -4,  7*0, &
                                           -33, -20,  19,   9,  0,  8*0, &
                                                -19,  -5,  11, -1,  9*0, &
                                                     -19,  -2,  4, 10*0, &
                                                            2,  0, 11*0, &
                                                               -6, 13*0/
!          DGRF h(n,m) for 1945:
      DATA (H1DIM(I),I=1297,1440) /13*0, &
          5810, -1702,  -499,  144,  -12,    6, -45,  12, -27,  5,  3*0, &
                  477,   186, -276,   95,  100, -18, -21,  17,  1,  4*0, &
                         -11,  -55,  -67,   16,   2, -12,  29,-20,  5*0, &
                              -178, -119,   -9,   6,  -7,  -9, -1,  6*0, &
                                      82,  -16,  28,   2,   4, -6,  7*0, &
                                           -39, -17,  18,   9,  6,  8*0, &
                                                -22,   3,   6, -4,  9*0, &
                                                     -11,   1, -2, 10*0, &
                                                            8,  0, 11*0, &
                                                               -2, 13*0/
!          DGRF h(n,m) for 1950:
      DATA (H1DIM(I),I=1441,1584) /13*0, &
          5815, -1810,  -476,  136,    3,   -1, -35,   5, -24, 13,  3*0, &
                  381,   206, -278,  103,   99, -17, -22,  19, -2,  4*0, &
                         -46,  -37,  -87,   33,   0,   0,  12,-10,  5*0, &
                              -210, -122,  -12,  10, -21,   2,  2,  6*0, &
                                      80,  -12,  36,  -8,   2, -3,  7*0, &
                                           -30, -18,  17,   8,  6,  8*0, &
                                                -16,  -4,   8, -3,  9*0, &
                                                     -17, -11,  6, 10*0, &
                                                           -7, 11, 11*0, &
                                                                8, 13*0/
!          DGRF h(n,m) for 1955:
      DATA (H1DIM(I),I=1585,1728) /13*0, &
          5820, -1898,  -462,  133,   15,   -9, -50,  10, -11, -4,  3*0, &
                  291,   216, -274,  110,   96, -24, -15,  12,  0,  4*0, &
                         -83,  -23,  -98,   48,  -4,   5,   7, -8,  5*0, &
                              -230, -121,  -16,   8, -23,   6, -2,  6*0, &
                                      78,  -12,  28,   3,  -2, -4,  7*0, &
                                           -24, -20,  23,  10,  1,  8*0, &
                                                -18,  -4,   7, -3,  9*0, &
                                                     -13,  -6,  7, 10*0, &
                                                            5, -1, 11*0, &
                                                               -3, 13*0/
!          DGRF h(n,m) for 1960:
      DATA (H1DIM(I),I=1729,1872) /13*0, &
          5791, -1967,  -414,  135,   16,  -10, -55,  11, -18,  4,  3*0, &
                  206,   224, -278,  125,   99, -28, -14,  12,  1,  4*0, &
                        -130,    3, -117,   60,  -6,   7,   2,  0,  5*0, &
                              -255, -114,  -20,   7, -18,   0,  2,  6*0, &
                                      81,  -11,  23,   4,  -3, -5,  7*0, &
                                           -17, -18,  23,   9,  1,  8*0, &
                                                -17,   1,   8, -1,  9*0, &
                                                     -20,   0,  6, 10*0, &
                                                            5,  0, 11*0, &
                                                               -7, 13*0/
!          DGRF h(n,m) for 1965:
      DATA (H1DIM(I),I=1873,2016) /13*0, &
          5776, -2016,  -404,  148,   19,  -11, -61,   7, -22,  2,  3*0, &
                  114,   240, -269,  128,  100, -27, -12,  15,  1,  4*0, &
                        -165,   13, -126,   68,  -2,   9,   7,  2,  5*0, &
                              -269,  -97,  -32,   6, -16,  -4,  6,  6*0, &
                                      81,   -8,  26,   4,  -5, -4,  7*0, &
                                            -7, -23,  24,  10,  0,  8*0, &
                                                -12,  -3,  10, -2,  9*0, &
                                                     -17,  -4,  3, 10*0, &
                                                            1,  0, 11*0, &
                                                               -6, 13*0/
!          DGRF h(n,m) for 1970:
      DATA (H1DIM(I),I=2017,2160) /13*0, &
          5737, -2047,  -366,  167,   26,  -12, -70,   7, -21,  1,  3*0, &
                   25,   251, -266,  139,  100, -27, -15,  16,  1,  4*0, &
                        -196,   26, -139,   72,  -4,   6,   6,  3,  5*0, &
                              -279,  -91,  -37,   8, -17,  -4,  4,  6*0, &
                                      83,   -6,  23,   6,  -5, -4,  7*0, &
                                             1, -23,  21,  10,  0,  8*0, &
                                                -11,  -6,  11, -1,  9*0, &
                                                     -16,  -2,  3, 10*0, &
                                                            1,  1, 11*0, &
                                                               -4, 13*0/
!          DGRF h(n,m) for 1975:
      DATA (H1DIM(I),I=2161,2304) /13*0, &
          5675, -2067,  -333,  191,   31,  -13, -77,   6, -21,  1,  3*0, &
                  -68,   262, -265,  148,   99, -26, -16,  16,  1,  4*0, &
                        -223,   39, -152,   75,  -5,   4,   7,  3,  5*0, &
                              -288,  -83,  -41,  10, -19,  -4,  4,  6*0, &
                                      88,   -4,  22,   6,  -5, -4,  7*0, &
                                            11, -23,  18,  10, -1,  8*0, &
                                                -12, -10,  11, -1,  9*0, &
                                                     -17,  -3,  3, 10*0, &
                                                            1,  1, 11*0, &
                                                               -5, 13*0/
!          DGRF h(n,m) for 1980:
      DATA (H1DIM(I),I=2305,2448) /13*0, &
          5604, -2129,  -336,  212,   46,  -15, -82,   7, -21,  1,  3*0, &
                 -200,   271, -257,  150,   93, -27, -18,  16,  0,  4*0, &
                        -252,   53, -151,   71,  -5,   4,   9,  3,  5*0, &
                              -297,  -78,  -43,  16, -22,  -5,  6,  6*0, &
                                      92,   -2,  18,   9,  -6, -4,  7*0, &
                                            17, -23,  16,   9,  0,  8*0, &
                                                -10, -13,  10, -1,  9*0, &
                                                     -15,  -6,  4, 10*0, &
                                                            2,  0, 11*0, &
                                                               -6, 13*0/
!          DGRF h(n,m) for 1985:
      DATA (H1DIM(I),I=2449,2592) /13*0, &
          5500, -2197,  -310,  232,   47,  -16, -83,   8, -21,  1,  3*0, &
                 -306,   284, -249,  150,   88, -27, -19,  15,  0,  4*0, &
                        -297,   69, -154,   69,  -2,   5,   9,  3,  5*0, &
                              -297,  -75,  -48,  20, -23,  -6,  6,  6*0, &
                                      95,   -1,  17,  11,  -6, -4,  7*0, &
                                            21, -23,  14,   9,  0,  8*0, &
                                                 -7, -15,   9, -1,  9*0, &
                                                     -11,  -7,  4, 10*0, &
                                                            2,  0, 11*0, &
                                                               -6, 13*0/
!          DGRF h(n,m) for 1990:
      DATA (H1DIM(I),I=2593,2736) /13*0, &
          5406, -2279,  -284,  247,   46,  -16, -80,  10, -20,  2,  3*0, &
                 -373,   293, -240,  154,   82, -26, -19,  15,  1,  4*0, &
                        -352,   84, -153,   69,   0,   6,  11,  3,  5*0, &
                              -299,  -69,  -52,  21, -22,  -7,  6,  6*0, &
                                      97,    1,  17,  12,  -7, -4,  7*0, &
                                            24, -23,  12,   9,  0,  8*0, &
                                                 -4, -16,   8, -2,  9*0, &
                                                     -10,  -7,  3, 10*0, &
                                                            2, -1, 11*0, &
                                                               -6, 13*0/
!          DGRF h(n,m) for 1995:
      DATA (H1DIM(I),I=2737,2880) /13*0, &
          5318, -2356,  -263,  262,   44,  -16, -77,  12, -19,  2,  3*0, &
                 -425,   302, -232,  157,   77, -25, -20,  15,  1,  4*0, &
                        -406,   98, -152,   67,   3,   7,  11,  3,  5*0, &
                              -301,  -64,  -57,  22, -21,  -7,  6,  6*0, &
                                      99,    4,  16,  12,  -7, -4,  7*0, &
                                            28, -23,  10,   9,  0,  8*0, &
                                                 -3, -17,   7, -2,  9*0, &
                                                     -10,  -8,  3, 10*0, &
                                                            1, -1, 11*0, &
                                                               -6, 13*0/
!          Initial coefficients g0 (IGRF for 2000):
      DATA G0 /0, &
        -29615, -2267,  1341,  935, -217,   72,  79,  25,   5, -2,  2*0, &
         -1728,  3072, -2290,  787,  351,   68, -74,   6,   9, -6,  3*0, &
                 1672,  1253,  251,  222,   74,   0,  -9,   3,  2,  4*0, &
                         715, -405, -131, -161,  33,  -8,  -8, -3,  5*0, &
                               110, -169,   -5,   9, -17,   6,  0,  6*0, &
                                     -12,   17,   7,   9,  -9,  4,  7*0, &
                                           -91,   8,   7,  -2,  1,  8*0, &
                                                 -2,  -8,   9,  2,  9*0, &
                                                      -7,  -4,  4, 10*0, &
                                                           -8,  0, 11*0, &
                                                               -1, 13*0/
!          D_/Dtime coefficients gt (IGRF for 2000-2005):
      DATA GT /0, &
          14.6_r8, -12.4_r8,   0.7_r8, -1.3_r8,  0.0_r8,  1.0_r8,-0.4_r8,-0.3_r8, 0.0_r8,0.0_r8,  2*0, &
          10.7_r8,   1.1_r8,  -5.4_r8,  1.6_r8, -0.7_r8, -0.4_r8,-0.4_r8, 0.2_r8, 0.0_r8,0.0_r8,  3*0, &
                 -1.1_r8,   0.9_r8, -7.3_r8, -2.1_r8,  0.9_r8,-0.3_r8,-0.3_r8, 0.0_r8,0.0_r8,  4*0, &
                        -7.7_r8,  2.9_r8, -2.8_r8,  2.0_r8, 1.1_r8, 0.4_r8, 0.0_r8,0.0_r8,  5*0, &
                              -3.2_r8, -0.8_r8, -0.6_r8, 1.1_r8,-1.0_r8, 0.0_r8,0.0_r8,  6*0, &
                                     2.5_r8, -0.3_r8,-0.2_r8, 0.3_r8, 0.0_r8,0.0_r8,  7*0, &
                                           1.2_r8, 0.6_r8,-0.5_r8, 0.0_r8,0.0_r8,  8*0, &
                                               -0.9_r8,-0.7_r8, 0.0_r8,0.0_r8,  9*0, &
                                                    -0.4_r8, 0.0_r8,0.0_r8, 10*0, &
                                                          0.0_r8,0.0_r8, 11*0, &
                                                              0.0_r8, 13*0/
!          Initial coefficients h0 (IGRF for 2000):
      DATA H0 /13*0, &
          5186, -2478,  -227,  272,   44,  -17, -65,  12, -20,  1,  3*0, &
                 -458,   296, -232,  172,   64, -24, -22,  13,  0,  4*0, &
                        -492,  119, -134,   65,   6,   8,  12,  4,  5*0, &
                              -304,  -40,  -61,  24, -21,  -6,  5,  6*0, &
                                     107,    1,  15,  15,  -8, -6,  7*0, &
                                            44, -25,   9,   9, -1,  8*0, &
                                                 -6, -16,   4, -3,  9*0, &
                                                      -3,  -8,  0, 10*0, &
                                                            5, -2, 11*0, &
                                                               -8, 13*0/
!          D_/Dtime coefficients ht (IGRF for 2000-2005):
      DATA HT /13*0, &
         -22.5_r8, -20.6_r8,   6.0_r8,  2.1_r8, -0.1_r8, -0.2_r8, 1.1_r8, 0.1_r8, 0.0_r8,0.0_r8,  3*0, &
                 -9.6_r8,  -0.1_r8,  1.3_r8,  0.6_r8, -1.4_r8, 0.0_r8, 0.0_r8, 0.0_r8,0.0_r8,  4*0, &
                       -14.2_r8,  5.0_r8,  1.7_r8,  0.0_r8, 0.3_r8, 0.0_r8, 0.0_r8,0.0_r8,  5*0, &
                               0.3_r8,  1.9_r8, -0.8_r8,-0.1_r8, 0.3_r8, 0.0_r8,0.0_r8,  6*0, &
                                     0.1_r8,  0.0_r8,-0.6_r8, 0.6_r8, 0.0_r8,0.0_r8,  7*0, &
                                           0.9_r8,-0.7_r8,-0.4_r8, 0.0_r8,0.0_r8,  8*0, &
                                                0.2_r8, 0.3_r8, 0.0_r8,0.0_r8,  9*0, &
                                                     0.7_r8, 0.0_r8,0.0_r8, 10*0, &
                                                          0.0_r8,0.0_r8, 11*0, &
                                                              0.0_r8, 13*0/
!
!-----------------------------------------------------------------------
!
!          Do not need to load new coefficients if date has not changed
      ICHG = 0
      IF (DATE .EQ. DATEL) GO TO 300
      DATEL = DATE
      ICHG = 1
 
!          Trap out of range date:
      IF (DATE .LT. EPOCH(1)) GO TO 9100
      IF (DATE .GT. EPOCH(NYT)+5._r8) write(iulog,9200) DATE
 
      DO 100 I=1,NYT
      IF (DATE .LT. EPOCH(I)) GO TO 110
      IY = I
  100 CONTINUE
  110 CONTINUE
 
      TIME = DATE
      T = TIME-EPOCH(IY)
      G(1)  = 0.0_r8
      GV(1) = 0.0_r8
      I = 2
      F0 = -1.0E-5_r16
      DO 200 N=1,NMAX
      F0 = F0 * REAL(N,r8)/2._r8
      F  = F0 / SQRT(2.0_r8)
      NN = N+1
      MM = 1
      IF (IY .EQ. NYT) THEN
!          Extrapolate coefficients
        G(I) = ((GTT(NN,MM)*T + GT(NN,MM))*T + G0(NN,MM)) * F0
      ELSE
!          Interpolate coefficients
        G(I) = (GYR(NN,MM,IY) + &
               T/5.0_r8 * (GYR(NN,MM,IY+1)-GYR(NN,MM,IY))) * F0
      ENDIF
      GV(I) = G(I) / REAL(NN,r8)
      I = I+1
      DO 200 M=1,N
      F = F / SQRT( REAL(N-M+1,r8) / REAL(N+M,r8) )
      NN = N+1
      MM = M+1
      I1 = I+1
      IF (IY .EQ. NYT) THEN
!          Extrapolate coefficients
        G(I)  = ((GTT(NN,MM)*T + GT(NN,MM))*T + G0(NN,MM)) * F
        G(I1) = ((HTT(NN,MM)*T + HT(NN,MM))*T + H0(NN,MM)) * F
      ELSE
!          Interpolate coefficients
        G(I)  = (GYR(NN,MM,IY) + &
                T/5.0_r8 * (GYR(NN,MM,IY+1)-GYR(NN,MM,IY))) * F
        G(I1) = (HYR(NN,MM,IY) + &
                T/5.0_r8 * (HYR(NN,MM,IY+1)-HYR(NN,MM,IY))) * F
      ENDIF
      RNN = REAL(NN,r8)
      GV(I)  = G(I)  / RNN
      GV(I1) = G(I1) / RNN
  200 I = I+2
 
  300 CONTINUE

      RETURN
 
!          Error trap diagnostics:
 9100 write(iulog,'('' '',/,'' COFRM:  DATE'',F9.3,'' preceeds DGRF coefficients'','' presently coded.'')') DATE
!     STOP 'mor cod'
      call endrun
 9200 FORMAT(' ',/,' COFRM:  DATE',F9.3,' is after the maximum',' recommended for extrapolation.')
      END SUBROUTINE COFRM
 
!================================================================================================

      SUBROUTINE DYPOL (COLAT,ELON,VP)
!
!-----------------------------------------------------------------------
!          Computes parameters for dipole component of geomagnetic field.
!          COFRM must be called before calling DYPOL!
!          940504 A. D. Richmond
!
!          INPUT from COFRM through COMMON /MAGCOF/ NMAX,GB(144),GV(144),ICHG
!            NMAX = Maximum order of spherical harmonic coefficients used
!            GB   = Coefficients for magnetic field calculation
!            GV   = Coefficients for magnetic potential calculation
!            ICHG = Flag indicating when GB,GV have been changed
!
!          RETURNS:
!            COLAT = Geocentric colatitude of geomagnetic dipole north pole
!                    (deg)
!            ELON  = East longitude of geomagnetic dipole north pole (deg)
!            VP    = Magnitude, in T.m, of dipole component of magnetic
!                    potential at geomagnetic pole and geocentric radius
!                    of 6371.2 km
!-----------------------------------------------------------------------
!
      use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none
!
!-------------------------------Commons---------------------------------
!
!     COMMON /MAGCOF/ NMAX,G(144),GV(144),ICHG

      integer nmax, ichg
      COMMON /MAGCOF_int/ NMAX,ICHG

      real(r8) g, gv
      COMMON /MAGCOF_real/ G(144),GV(144)
!
!------------------------------Arguments--------------------------------
!
      real(r8) colat, elon, vp
!
!-----------------------------Parameters------------------------------
!
      real(r8) rtod, dtor, re, req
      PARAMETER (RTOD=5.72957795130823E1_r8 , DTOR=1.745329251994330E-2_r8, &
                 RE=6371.2_r8, REQ=6378.160_r8)
!
!---------------------------Local variables-----------------------------
!
      real(r8) gpl, stp, ctp
!
!-----------------------------------------------------------------------
!
!          Compute geographic colatitude and longitude of the north pole of
!          earth centered dipole
      GPL = SQRT( G(2  )**2+ G(3  )**2+ G(4  )**2)
      CTP = G(2  )/GPL
      STP = SQRT(1._r8 - CTP*CTP)
      COLAT = (ACOS(CTP))*RTOD
      ELON = ATAN2( G(4  ), G(3  ))*RTOD
 
!          Compute magnitude of magnetic potential at pole, radius Re.
      VP = .2_r8*GPL*RE
!          .2 = 2*(10**-4 T/gauss)*(1000 m/km) (2 comes through F0 in COFRM).
 
      RETURN
      END SUBROUTINE DYPOL
 
!================================================================================================

      SUBROUTINE FELDG (IENTY,GLAT,GLON,ALT, BNRTH,BEAST,BDOWN,BABS)
!
!-----------------------------------------------------------------------
!          Compute the DGRF/IGRF field components at the point GLAT,GLON,ALT.
!          COFRM must be called to establish coefficients for correct date
!          prior to calling FELDG.
!
!          IENTY is an input flag controlling the meaning and direction of the
!                remaining formal arguments:
!          IENTY = 1
!            INPUTS:
!              GLAT = Latitude of point (deg)
!              GLON = Longitude (east=+) of point (deg)
!              ALT  = Ht of point (km)
!            RETURNS:
!              BNRTH  north component of field vector (Gauss)
!              BEAST  east component of field vector  (Gauss)
!              BDOWN  downward component of field vector (Gauss)
!              BABS   magnitude of field vector (Gauss)
!
!          IENTY = 2
!            INPUTS:
!              GLAT = X coordinate (in units of earth radii 6371.2 km )
!              GLON = Y coordinate (in units of earth radii 6371.2 km )
!              ALT  = Z coordinate (in units of earth radii 6371.2 km )
!            RETURNS:
!              BNRTH = X component of field vector (Gauss)
!              BEAST = Y component of field vector (Gauss)
!              BDOWN = Z component of field vector (Gauss)
!              BABS  = Magnitude of field vector (Gauss)
!          IENTY = 3
!            INPUTS:
!              GLAT = X coordinate (in units of earth radii 6371.2 km )
!              GLON = Y coordinate (in units of earth radii 6371.2 km )
!              ALT  = Z coordinate (in units of earth radii 6371.2 km )
!            RETURNS:
!              BNRTH = Dummy variable
!              BEAST = Dummy variable
!              BDOWN = Dummy variable
!              BABS  = Magnetic potential (T.m)
!
!          INPUT from COFRM through COMMON /MAGCOF/ NMAX,GB(144),GV(144),ICHG
!            NMAX = Maximum order of spherical harmonic coefficients used
!            GB   = Coefficients for magnetic field calculation
!            GV   = Coefficients for magnetic potential calculation
!            ICHG = Flag indicating when GB,GV have been changed
!
!          HISTORY:
!          COFRM and FELDG originated 15 Apr 83 by Vincent B. Wickwar
!          (formerly at SRI. Int., currently at Utah State).
!
!          May 94 (A.D. Richmond): Added magnetic potential calculation
!
!          Oct 95 (Barnes): Added ICHG
!-----------------------------------------------------------------------
!
      use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none
!
!-------------------------------Commons---------------------------------
!
!     COMMON /MAGCOF/ NMAX,GB(144),GV(144),ICHG

      integer nmax, ichg
      COMMON /MAGCOF_int/ NMAX,ICHG

      real(r8) gb, gv
      COMMON /MAGCOF_real/ GB(144),GV(144)
!
!------------------------------Arguments--------------------------------
!
      integer ienty

      real(r8) glat, glon, alt, bnrth, beast, bdown, babs
!
!-----------------------------Parameters------------------------------
!
      real(r8) rtod, dtor, re, req
      PARAMETER (RTOD=57.2957795130823_r8, DTOR=0.01745329251994330_r8, &
                 RE=6371.2_r8, REQ=6378.160_r8)
!
!---------------------------Local variables-----------------------------
!
      integer m, ih, k, il, ihm, ilm, is, imax, last, mk, i, ihmax

      real(r8) y, x, z, f, byyy, bxxx, brho, bzzz, t, s, cp, rlon, xxx
      real(r8) sp, rlat, st, ct, zzz, yyy, rq
      real(r8) XI(3),H(144),G(144)
      SAVE G

      integer ientyp
      SAVE IENTYP
      DATA IENTYP/-10000/
!
!-----------------------------------------------------------------------
!
      IF (IENTY .EQ. 1) THEN
        IS=1
        RLAT = GLAT*DTOR
        CT   = SIN(RLAT)
        ST   = COS(RLAT)
        RLON = GLON*DTOR
        CP   = COS(RLON)
        SP   = SIN(RLON)
        CALL GD2CART (GLAT,GLON,ALT,XXX,YYY,ZZZ)
        XXX = XXX/RE
        YYY = YYY/RE
        ZZZ = ZZZ/RE
      ELSE
        IS   = 2
        XXX  = GLAT
        YYY  = GLON
        ZZZ  = ALT
      ENDIF
      RQ    = 1._r8/(XXX**2+YYY**2+ZZZ**2)
      XI(1) = XXX*RQ
      XI(2) = YYY*RQ
      XI(3) = ZZZ*RQ
      IHMAX = NMAX*NMAX+1
      LAST  = IHMAX+NMAX+NMAX
      IMAX  = NMAX+NMAX-1
 
      IF (IENTY .NE. IENTYP .OR. ICHG .EQ. 1) THEN
        IENTYP = IENTY
	ICHG = 0
        IF (IENTY .NE. 3) THEN
	  DO 10 I=1,LAST
   10     G(I) = GB(I)
        ELSE
	  DO 20 I=1,LAST
   20     G(I) = GV(I)
        ENDIF
      ENDIF
 
      DO 30 I=IHMAX,LAST
   30 H(I) = G(I)

      MK = 3
      IF (IMAX .EQ. 1) MK=1

      DO 100 K=1,MK,2
      I  = IMAX
      IH = IHMAX

   60 IL = IH-I
      F = 2._r8/real(I-K+2,r8)
      X = XI(1)*F
      Y = XI(2)*F
      Z = XI(3)*(F+F)

      I = I-2
      IF (I .LT. 1) GO TO 90
      IF (I .EQ. 1) GO TO 80

      DO 70 M=3,I,2
      IHM = IH+M
      ILM = IL+M
      H(ILM+1) = G(ILM+1)+ Z*H(IHM+1) + X*(H(IHM+3)-H(IHM-1)) &
                                              -Y*(H(IHM+2)+H(IHM-2))
   70 H(ILM)   = G(ILM)  + Z*H(IHM)   + X*(H(IHM+2)-H(IHM-2)) &
                                              +Y*(H(IHM+3)+H(IHM-1))

   80 H(IL+2) = G(IL+2) + Z*H(IH+2) + X*H(IH+4) - Y*(H(IH+3)+H(IH))
      H(IL+1) = G(IL+1) + Z*H(IH+1) + Y*H(IH+4) + X*(H(IH+3)-H(IH))

   90 H(IL)   = G(IL)   + Z*H(IH)   + 2._r8*(X*H(IH+1)+Y*H(IH+2))
      IH = IL
      IF (I .GE. K) GO TO 60
  100 CONTINUE
 
      S = .5_r8*H(1)+2._r8*(H(2)*XI(3)+H(3)*XI(1)+H(4)*XI(2))
      T = (RQ+RQ)*SQRT(RQ)
      BXXX = T*(H(3)-S*XXX)
      BYYY = T*(H(4)-S*YYY)
      BZZZ = T*(H(2)-S*ZZZ)
      BABS = SQRT(BXXX**2+BYYY**2+BZZZ**2)
      IF (IS .EQ. 1) THEN            ! (convert back to geodetic)
        BEAST = BYYY*CP-BXXX*SP
        BRHO  = BYYY*SP+BXXX*CP
        BNRTH =  BZZZ*ST-BRHO*CT
        BDOWN = -BZZZ*CT-BRHO*ST
      ELSEIF (IS .EQ. 2) THEN        ! (leave in earth centered cartesian)
        BNRTH = BXXX
        BEAST = BYYY
        BDOWN = BZZZ
      ENDIF
 
!          Magnetic potential computation makes use of the fact that the
!          calculation of V is identical to that for r*Br, if coefficients
!          in the latter calculation have been divided by (n+1) (coefficients
!          GV).  Factor .1 converts km to m and gauss to tesla.
      IF (IENTY.EQ.3) BABS = (BXXX*XXX + BYYY*YYY + BZZZ*ZZZ)*RE*.1_r8
 
      RETURN
      END SUBROUTINE FELDG
 
!================================================================================================

      SUBROUTINE GD2CART (GDLAT,GLON,ALT,X,Y,Z)
!
!-----------------------------------------------------------------------
!          Convert geodetic to cartesian coordinates by calling CONVRT
!          940503 A. D. Richmond
!-----------------------------------------------------------------------
!
      use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none
!
!------------------------------Arguments--------------------------------
!
      real(r8) gdlat, glon, alt, x, y, z
!
!-----------------------------Parameters------------------------------
!
      real(r8) rtod, dtor
      PARAMETER (RTOD=57.2957795130823_r8, DTOR=0.01745329251994330_r8)
!
!---------------------------Local variables-----------------------------
!
      real(r8) ang, rho
!
!-----------------------------------------------------------------------
!
      CALL CONVRT (1,GDLAT,ALT,RHO,Z)
      ANG = GLON*DTOR
      X = RHO*COS(ANG)
      Y = RHO*SIN(ANG)
      RETURN
      END SUBROUTINE GD2CART
 
!================================================================================================

      SUBROUTINE CONVRT (I,GDLAT,ALT,X1,X2)
!
!-----------------------------------------------------------------------
!          Convert space point from geodetic to geocentric or vice versa.
!
!          I is an input flag controlling the meaning and direction of the
!            remaining formal arguments:
!
!          I = 1  (convert from geodetic to cylindrical)
!            INPUTS:
!              GDLAT = Geodetic latitude (deg)
!              ALT   = Altitude above reference ellipsoid (km)
!            RETURNS:
!              X1    = Distance from Earth's rotation axis (km)
!              X2    = Distance above (north of) Earth's equatorial plane (km)
!
!          I = 2  (convert from geodetic to geocentric spherical)
!            INPUTS:
!              GDLAT = Geodetic latitude (deg)
!              ALT   = Altitude above reference ellipsoid (km)
!            RETURNS:
!              X1    = Geocentric latitude (deg)
!              X2    = Geocentric distance (km)
!
!          I = 3  (convert from cylindrical to geodetic)
!            INPUTS:
!              X1    = Distance from Earth's rotation axis (km)
!              X2    = Distance from Earth's equatorial plane (km)
!            RETURNS:
!              GDLAT = Geodetic latitude (deg)
!              ALT   = Altitude above reference ellipsoid (km)
!
!          I = 4  (convert from geocentric spherical to geodetic)
!            INPUTS:
!              X1    = Geocentric latitude (deg)
!              X2    = Geocentric distance (km)
!            RETURNS:
!              GDLAT = Geodetic latitude (deg)
!              ALT   = Altitude above reference ellipsoid (km)
!
!
!          HISTORY:
!          940503 (A. D. Richmond):  Based on a routine originally written
!          by V. B. Wickwar.
!
!          REFERENCE:  ASTRON. J. VOL. 66, p. 15-16, 1961.
!-----------------------------------------------------------------------
!
      use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none
!
!------------------------------Arguments--------------------------------
!
      integer i

      real(r8) gdlat, alt, x1, x2
!
!-----------------------------Parameters------------------------------
!
      real(r8) rtod, dtor, re, req, fltnvrs, e2, e4, e6, e8, ome2req
      real(r8) a21, a22, a23, a41, a42, a43, a44, a61, a62, a63, a81, a82
      real(r8) a83, a84
      PARAMETER (RTOD=57.2957795130823_r8, DTOR=0.01745329251994330_r8 , &
        RE=6371.2_r8 , REQ=6378.160_r8 , FLTNVRS=298.25_r8 , &
        E2=(2._r8-1._r8/FLTNVRS)/FLTNVRS , E4=E2*E2 , E6=E4*E2 , E8=E4*E4 , &
        OME2REQ = (1._r8-E2)*REQ , &
           A21 =     (512._r8*E2 + 128._r8*E4 + 60._r8*E6 + 35._r8*E8)/1024._r8 , &
           A22 =     (                        E6 +     E8)/  32._r8 , &
           A23 = -3._r8*(                     4._r8*E6 +  3._r8*E8)/ 256._r8 , &
           A41 =    -(           64._r8*E4 + 48._r8*E6 + 35._r8*E8)/1024._r8 , &
           A42 =     (            4._r8*E4 +  2._r8*E6 +     E8)/  16._r8 , &
           A43 =                                   15._r8*E8 / 256._r8 , &
           A44 =                                      -E8 /  16._r8 , &
           A61 =  3._r8*(                     4._r8*E6 +  5._r8*E8)/1024._r8 , &
           A62 = -3._r8*(                        E6 +     E8)/  32._r8 , &
           A63 = 35._r8*(                     4._r8*E6 +  3._r8*E8)/ 768._r8 , &
           A81 =                                   -5._r8*E8 /2048._r8 , &
           A82 =                                   64._r8*E8 /2048._r8 , &
           A83 =                                 -252._r8*E8 /2048._r8 , &
           A84 =                                  320._r8*E8 /2048._r8 )
!          E2 = Square of eccentricity
!
!---------------------------Local variables-----------------------------
!
      real(r8) a2, a4, a6, a8, ccl, c2cl, c4cl, dltcl, gclat, gclatm, &
           coslat, sinlat, s2cl, s6cl, sgl, s4cl, s8cl, ri, &
           rkm, scl, d, z, rho
!
!-----------------------------------------------------------------------
!
      IF (I .GE. 3) GO TO 300
 
!          Geodetic to geocentric
 
!          Compute RHO,Z
      SINLAT = SIN(GDLAT*DTOR)
      COSLAT = SQRT(1._r8-SINLAT*SINLAT)
      D      = SQRT(1._r8-E2*SINLAT*SINLAT)
      Z      = (ALT+OME2REQ/D)*SINLAT
      RHO    = (ALT+REQ/D)*COSLAT
      X1 = RHO
      X2 = Z
      IF (I .EQ. 1) RETURN
 
!          Compute GCLAT,RKM
      RKM   = SQRT(Z*Z + RHO*RHO)
      GCLAT = RTOD*ATAN2(Z,RHO)
      X1 = GCLAT
      X2 = RKM
      RETURN
 
!          Geocentric to geodetic
  300 IF (I .EQ. 3) THEN
         RHO = X1
         Z = X2
         RKM = SQRT(Z*Z+RHO*RHO)
         SCL = Z/RKM
         GCLAT = ASIN(SCL)*RTOD
      ELSEIF (I .EQ. 4) THEN
         GCLAT = X1
         RKM = X2
         SCL = SIN(GCLAT*DTOR)
      ELSE
         RETURN
      ENDIF
 
      RI = REQ/RKM
      A2 = RI*(A21+RI*(A22+RI* A23))
      A4 = RI*(A41+RI*(A42+RI*(A43+RI*A44)))
      A6 = RI*(A61+RI*(A62+RI* A63))
      A8 = RI*(A81+RI*(A82+RI*(A83+RI*A84)))
      CCL = SQRT(1._r8-SCL*SCL)
      S2CL = 2._r8*SCL*CCL
      C2CL = 2._r8*CCL*CCL-1._r8
      S4CL = 2._r8*S2CL*C2CL
      C4CL = 2._r8*C2CL*C2CL-1._r8
      S8CL = 2._r8*S4CL*C4CL
      S6CL = S2CL*C4CL+C2CL*S4CL
      DLTCL = S2CL*A2+S4CL*A4+S6CL*A6+S8CL*A8
      GDLAT = DLTCL*RTOD+GCLAT
      SGL = SIN(GDLAT*DTOR)
      ALT = RKM*COS(DLTCL)-REQ*SQRT(1._r8-E2*SGL*SGL)
      RETURN
      END SUBROUTINE CONVRT

!================================================================================================

      block data blkmagcof
!
!-----------------------------------------------------------------------
!
      use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none
!
!-------------------------------Commons---------------------------------
!
!     COMMON /MAGCOF/ NMAX,G(144),GV(144),ICHG

      integer :: nmax,ichg
      COMMON /MAGCOF_int/ NMAX,ICHG

      real(r8) :: g,gv
      COMMON /MAGCOF_real/ G(144),GV(144)
!
!-----------------------------------------------------------------------
!
      DATA NMAX,ICHG /10,-99999/
      end block data blkmagcof

!================================================================================================

!-----------------------------------------------------------------------
! the following part is used to calculate magnetic local time
! the orginal code can be found at HAO machines ~bozo/apex/src 
! (Roy Barnes) 
!  A.Maute copied 2004/2/6 
!
!-----------------------------------------------------------------------
      SUBROUTINE MAGLOCTM (ALON,SBSLLAT,SBSLLON,CLATP,POLON,MLT)
!
!-----------------------------------------------------------------------
!  Computes magnetic local time from magnetic longitude, subsolar coordinates,
!   and geomagnetic pole coordinates.
!  950302 A. D. Richmond, NCAR
!  Algorithm:  MLT is calculated from the difference of the apex longitude,
!   alon, and the geomagnetic dipole longitude of the subsolar point.
!
!   Inputs:
!    ALON    = apex magnetic longitude of the point (deg)
!    SBSLLAT = geographic latitude of subsolar point (degrees)
!    SBSLLON = geographic longitude of subsolar point (degrees)
!    CLATP   = Geocentric colatitude of geomagnetic dipole north pole (deg)
!    POLON   = East longitude of geomagnetic dipole north pole (deg)
!
!   Output:
!    mlt (real) = magnetic local time for the apex longitude alon (hours)
!
!-----------------------------------------------------------------------
!
        use shr_kind_mod, only: r8 => shr_kind_r8
        implicit none
!
!------------------------------Arguments--------------------------------
!
	REAL(r8) alon, sbsllat, sbsllon, clatp, polon, MLT
!
!---------------------------Local variables-----------------------------
!
        real(r8) smlon
!
!-----------------------------------------------------------------------
!
	CALL SOLGMLON (SBSLLAT,SBSLLON,CLATP,POLON,SMLON)
	MLT = (ALON - SMLON)/15.0_r8 + 12.0_r8
	IF (MLT .GE. 24.0_r8) MLT = MLT - 24.0_r8
	IF (MLT .LT.   0._r8) MLT = MLT + 24.0_r8
	RETURN
        END SUBROUTINE MAGLOCTM

!================================================================================================

      SUBROUTINE SOLGMLON (XLAT,XLON,COLAT,ELON,MLON)
!
!-----------------------------------------------------------------------
! Computes geomagnetic longitude of the point with geocentric spherical
!  latitude and longitude of XLAT and XLON, respectively.
! 940719 A. D. Richmond, NCAR
! Inputs:
!   XLAT  = geocentric spherical latitude (deg)
!   XLON  = geocentric spherical longitude (deg)
!   COLAT = Geocentric colatitude of geomagnetic dipole north pole (deg)
!   ELON  = East longitude of geomagnetic dipole north pole (deg)
! Output:
!   MLON  = Geomagnetic dipole longitude of the point (deg, -180. to 180.)
!-----------------------------------------------------------------------
!
      use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none
!
!------------------------------Arguments--------------------------------
!
      REAL(r8) xlat, xlon, colat, elon, MLON
!
!-----------------------------Parameters------------------------------
!
      real(r8) rtod, dtor
      PARAMETER (RTOD=5.72957795130823E1_r8,DTOR=1.745329251994330E-2_r8)
!
!---------------------------Local variables-----------------------------
!
      real(r8) ste, cte, stfspa, stfcpa, sang, stp, ctp, cang, ang
!
!-----------------------------------------------------------------------
!
! Algorithm:
!   Use spherical coordinates.
!   Let GP be geographic pole.
!   Let GM be geomagnetic pole (colatitude COLAT, east longitude ELON).
!   Let XLON be longitude of point P.
!   Let TE be colatitude of point P.
!   Let ANG be longitude angle from GM to P.
!   Let TP be colatitude of GM.
!   Let TF be arc length between GM and P.
!   Let PA = MLON be geomagnetic longitude, i.e., Pi minus angle measured
!     counterclockwise from arc GM-P to arc GM-GP.
!   Then, using notation C=cos, S=sin, spherical-trigonometry formulas
!     for the functions of the angles are as shown below.  Note: STFCPA,
!     STFSPA are sin(TF) times cos(PA), sin(PA), respectively.

      CTP = COS(COLAT*DTOR)
      STP = SQRT(1._r8 - CTP*CTP)
      ANG = (XLON-ELON)*DTOR
      CANG = COS(ANG)
      SANG = SIN(ANG)
      CTE = SIN(XLAT*DTOR)
      STE = SQRT(1._r8-CTE*CTE)
      STFCPA = STE*CTP*CANG - CTE*STP
      STFSPA = SANG*STE
      MLON = ATAN2(STFSPA,STFCPA)*RTOD
      RETURN
      END SUBROUTINE SOLGMLON

!================================================================================================

      SUBROUTINE SUBSOL (IYR,IDAY,IHR,IMN,SEC,SBSLLAT,SBSLLON)
!
!-----------------------------------------------------------------------
!          Find subsolar geographic latitude and longitude given the
!          date and time (Universal Time).
!
!          This is based on formulas in Astronomical Almanac for the
!          year 1996, p.  C24. (U.S.  Government Printing Office,
!          1994).  According to the Almanac, results are good to at
!          least 0.01 degree latitude and 0.025 degree longitude
!          between years 1950 and 2050.  Accuracy for other years has
!          not been tested although the algorithm has been designed to
!          accept input dates from 1601 to 2100.  Every day is assumed
!          to have exactly 86400 seconds; thus leap seconds that
!          sometimes occur on June 30 and December 31 are ignored:
!          their effect is below the accuracy threshold of the algorithm.
!
!          961026 A. D. Richmond, NCAR
!
!          INPUTS:
!            IYR  = Year (e.g., 1994). IYR must be in the range: 1601 to 2100.
!            IDAY = Day number of year (e.g., IDAY = 32 for Feb 1)
!            IHR  = Hour of day    (e.g., 13 for 13:49)
!            IMN  = Minute of hour (e.g., 49 for 13:49)
!            SEC  = Second and fraction after the hour/minute.
!          Note:  While IYR is bounds tested, there is no constraint
!                 placed on values: IDAY,IHR,IMN,SEC; e.g., IHR=25 is
!                 properly interpreted.
!          RETURNS:
!            SBSLLAT = geographic latitude of subsolar point (degrees)
!            SBSLLON = geographic longitude of subsolar point (degrees,
!                      between -180 and +180)
!-----------------------------------------------------------------------
!
      use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none
!
!------------------------------Arguments--------------------------------
!
      INTEGER IYR,IDAY,IHR,IMN

      REAL(r8) SEC,SBSLLAT,SBSLLON
!
!-----------------------------Parameters------------------------------
!
      integer msgun
      PARAMETER (MSGUN=6)

      real(r8) d2r, r2d
      PARAMETER (D2R=0.0174532925199432957692369076847_r8 , &
                 R2D=57.2957795130823208767981548147_r8)
!
!---------------------------Local variables-----------------------------
!
      INTEGER YR,NLEAP,NCENT,NROT

      REAL(r8) L0,G0,DF,LF,GF,L,G,LAMBDA,EPSILON,N &
          ,GRAD,LAMRAD,SINLAM,EPSRAD,DELTA,UT,ETDEG
      real(r8) aptime, alpha
!
!-----------------------------------------------------------------------
!
! Number of years from 2000 to IYR (negative if IYR < 2000):
      YR = IYR - 2000
!
! NLEAP (final) = number of leap days from (2000 January 1) to (IYR January 1)
!                 (negative if IYR is before 1997)
      NLEAP = (IYR-1601)/4
      NLEAP = NLEAP - 99
      IF (IYR.LE.1900) THEN
	IF (IYR.LE.1600) THEN
	 WRITE(MSGUN,*) 'SUBSOLR INVALID BEFORE 1601: INPUT YEAR = ',IYR       
	 STOP
	ENDIF
	NCENT = (IYR-1601)/100
	NCENT = 3 - NCENT 
	NLEAP = NLEAP + NCENT
      ENDIF
      IF (IYR.GE.2101) THEN
	WRITE(MSGUN,*) 'SUBSOLR INVALID AFTER 2100:  INPUT YEAR = ',IYR
	STOP
      ENDIF
!
! L0 = Mean longitude of Sun at 12 UT on January 1 of IYR:
!     L0 = 280.461 + .9856474*(365*(YR-NLEAP) + 366*NLEAP) 
!	   - (ARBITRARY INTEGER)*360.
!        = 280.461 + .9856474*(365*(YR-4*NLEAP) + (366+365*3)*NLEAP) 
!	   - (ARBITRARY INTEGER)*360.
!        = (280.461 - 360.) + (.9856474*365 - 360.)*(YR-4*NLEAP) 
!	   + (.9856474*(366+365*3) - 4*360.)*NLEAP,
!  where ARBITRARY INTEGER = YR+1.  This gives:
      L0 = -79.549_r8 + (-.238699_r8*(YR-4*NLEAP) + 3.08514E-2_r8*NLEAP)
!
! G0 = Mean anomaly at 12 UT on January 1 of IYR:
!     G0 = 357.528 + .9856003*(365*(YR-NLEAP) + 366*NLEAP) 
!	   - (ARBITRARY INTEGER)*360.
!        = 357.528 + .9856003*(365*(YR-4*NLEAP) + (366+365*3)*NLEAP) 
!	   - (ARBITRARY INTEGER)*360.
!        = (357.528 - 360.) + (.9856003*365 - 360.)*(YR-4*NLEAP) 
!	   + (.9856003*(366+365*3) - 4*360.)*NLEAP,
!  where ARBITRARY INTEGER = YR+1.  This gives:
      G0 = -2.472_r8 + (-.2558905_r8*(YR-4*NLEAP) - 3.79617E-2_r8*NLEAP)
!
! Universal time in seconds:
      UT = real(IHR*3600 + IMN*60,r8) + SEC
!
! Days (including fraction) since 12 UT on January 1 of IYR:
      DF = (UT/86400._r8 - 1.5_r8) + IDAY
!
! Addition to Mean longitude of Sun since January 1 of IYR:
      LF = .9856474_r8*DF
!
! Addition to Mean anomaly since January 1 of IYR:
      GF = .9856003_r8*DF
!
! Mean longitude of Sun:
      L = L0 + LF
!
! Mean anomaly:
      G = G0 + GF
      GRAD = G*D2R
!
! Ecliptic longitude:
      LAMBDA = L + 1.915_r8*SIN(GRAD) + .020_r8*SIN(2._r8*GRAD)
      LAMRAD = LAMBDA*D2R
      SINLAM = SIN(LAMRAD)
!
! Days (including fraction) since 12 UT on January 1 of 2000:
      N = DF + real(365*YR + NLEAP,r8)
!
! Obliquity of ecliptic:
      EPSILON = 23.439_r8 - 4.E-7_r8*N
      EPSRAD = EPSILON*D2R
!
! Right ascension:
      ALPHA = ATAN2(COS(EPSRAD)*SINLAM,COS(LAMRAD))*R2D
!
! Declination:
      DELTA = ASIN(SIN(EPSRAD)*SINLAM)*R2D
!
! Subsolar latitude:
      SBSLLAT = DELTA
!
! Equation of time (degrees):
      ETDEG = L - ALPHA
      NROT = NINT(ETDEG/360._r8)
      ETDEG = ETDEG - real(360*NROT,r8)
!
! Apparent time (degrees):
      APTIME = UT/240._r8 + ETDEG
!          Earth rotates one degree every 240 s.
!
! Subsolar longitude:
      SBSLLON = 180._r8 - APTIME
      NROT = NINT(SBSLLON/360._r8)
      SBSLLON = SBSLLON - real(360*NROT,r8)
!
      RETURN
      END SUBROUTINE SUBSOL

!================================================================================================

