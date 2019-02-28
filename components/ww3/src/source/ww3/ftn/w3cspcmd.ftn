#include "w3macros.h"
!/ ------------------------------------------------------------------- /
      MODULE W3CSPCMD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         01-Nov-2012 |
!/                  +-----------------------------------+
!/
!/    19-Sep-2005 : Origination.                        ( version 3.08 )
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/    01-Nov-2012 : Minor code clean-up (tabs & coments)( version 4.08 )
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS. 
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Convert spectra to new discrete spectral grid.
!
!  2. Variables and types :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      NCASES    Int.  Private  Number of cases for which interpol.
!                               data is stored.
!      IDATA     CASE  Private  Interpolation data.
!     ----------------------------------------------------------------
!
!     Elements of the data structure CASE are given below. The middle
!     block pf parameters has pointer aliasses with the same name in
!     the subroutine.
!  
!      Name      Type  Description
!     ----------------------------------------------------------------
!      ICASE     Int.  Number of case.
!      NFR1, NTH1, NFR2, NTH2, XF1, FR1, TH1, XF2, FR2, TH2
!                      Same as in parameter list of routine.
!
!      DTH1      Real  Directional increment.
!      DTH2      Real  Directional increment.
!      IDTH      I.A.  Index information for redistribution of
!                      energy in direction space.
!      RDTH      R.A.  Factors corresponding to IDTH.
!      FRQ1      R.A.  Frequencies.
!      FRQ2      R.A.  Frequencies.
!      XDF1      Real  Factor for increments.
!      XDF2      Real  Factor for increments.
!      NFR2T     Int.  Frequency to start the tail.
!      IDFR      I.A.  Idem for frequencies.
!      RDFR      R.A.  Factors corresponding to IDFR.
!
!      NEXT      CASE  Pointer to next data set stored.
!     ----------------------------------------------------------------
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      W3CSPC    Subr. Public   Perform conversion for vector of
!                               spectra.
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!     See subroutine W3CSPC.
!
!  5. Remarks :
!
!     - Conversion data are sored in an endless linked chain, which
!       is tested at the beginning of the routine.
!
!  6. Switches :
!
!     See subroutine.
!
!  7. Source code :
!
!/ ------------------------------------------------------------------- /
      PUBLIC
!/
      TYPE CASE
      INTEGER               :: ICASE, NFR1, NTH1, NFR2, NTH2, NFR2T
      REAL                  :: XF1, FR1, TH1, XF2, FR2, TH2,       &
                               DTH1, DTH2, XDF1, XDF2
      INTEGER, POINTER      :: IDTH(:,:), IDFR(:,:)
      REAL, POINTER         :: RDTH(:,:), FRQ1(:), FRQ2(:), RDFR(:,:)
      TYPE(CASE), POINTER   :: NEXT
      END TYPE CASE
!/
      INTEGER, PRIVATE        :: NCASES = 0
      TYPE(CASE), PRIVATE, POINTER :: IDATA
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3CSPC ( SP1, NFR1, NTH1, XF1, FR1, TH1,             &
                          SP2, NFR2, NTH2, XF2, FR2, TH2,             &
                          NSP, NDST, NDSE, FTL )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         01-Nov-2012 !
!/                  +-----------------------------------+
!/
!/    19-Sep-2005 : Origination.                        ( version 3.08 )
!/    01-Nov-2012 : code clean up (tab spaces, comments)( version 4.08 ) 
!/
!  1. Purpose :
!
!     Convert a set of spectra to a new spectral grid.
!
!  2. Method :
!
!     Conservative distribution of input energies over new grid.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       SP1     R.A.   I   Input spectra.
!       NFR1    Int.   I   Input number of frequencies.
!       NTH1    Int.   I   Input number of directions.
!       XFR     Real   I   Input frequency increment factor.
!       FR1     Real   I   First input frequency.
!       TH1     Real   I   First input direction.
!       SP2     R.A.   O   Output spectra.
!       NFR2, NTH2, XF2, FR2, TH2
!                      !   Specral description for output spectra.
!       NSP     Int.   I   Number of spectra.
!       NDST    int.   I   Unit number for test output.
!       NDSE    int.   I   Unit number for error output.
!       FTAIL   Real   I   Factor for tail description = XF2**N
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Sur.  W3SERVMD Subroutine tracing.
!      EXTCDE    Sur.    Id     program abort.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3IOBC    Subr. W3IOBCMD Updating boundary conditions.
!                Subr           Multi scale model bound. data input.
!     ---------------------------------------------------------------- 
!
!  6. Error messages :
!
!     - Check on input parameters.
!
!  7. Remarks :
!
!     - The inner loop of the actual redistribution is over the 
!       individual spectra, optimizing this routine for large numbers
!       of conversions in a single call.
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S    Enable subroutine tracing.
!
!     !/T    Enable test output.
!     !/T1   Test output for searching in stored data.
!     !/T2   Test output for redistribution data.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE CONSTANTS
!
      USE W3SERVMD, ONLY: EXTCDE
!/S      USE W3SERVMD, ONLY: STRACE
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: NSP, NFR1, NTH1, NFR2, NTH2, NDST, NDSE
      REAL, INTENT(IN)        :: SP1(NTH1,NFR1,NSP), XF1, FR1, TH1,   &
                                 XF2, FR2, TH2, FTL
      REAL, INTENT(OUT)       :: SP2(NTH2,NFR2,NSP)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: I, NRMAX, J, I1, L1, J1, I2, L2, J2, &
                                 ISP
!/S      INTEGER, SAVE           :: IENT = 0
      REAL                    :: LOW, HGH, RLOW, RHGH, BLOW, BHGH,    &
                                 FRAC, AUX1, AUX2, R1, R2, FACT
      LOGICAL                 :: FOUND
      TYPE(CASE), POINTER     :: CURRENT
!/
!/ ------------------------------------------------------------------- /
!/ Pointers for aliases
!/ 
      INTEGER, POINTER        :: IDTH(:,:), IDFR(:,:), NFR2T
      REAL, POINTER           :: DTH1, DTH2, RDTH(:,:), FRQ1(:),      &
                                 FRQ2(:), XDF1, XDF2, RDFR(:,:)
!/
!/S      CALL STRACE (IENT, 'W3CSPC')
!
! -------------------------------------------------------------------- /
! 0.  Initializations
! 0.a Check input
!
      IF ( NFR1.LT.3 .OR. NTH1.LT.4 .OR. XF1.LE.1. .OR. FR1.LE.0. .OR.&
           NFR2.LT.3 .OR. NTH2.LT.4 .OR. XF2.LE.1. .OR. FR2.LE.0. ) THEN
        WRITE (NDSE,900) NFR1, NTH1, XF1, FR1, NFR2, NTH2, XF2, FR2
        CALL EXTCDE ( 1 )
        END IF
!
      IF ( NSP .LT. 0 ) THEN
        WRITE (NDSE,901)   
        CALL EXTCDE ( 2 )
        END IF
!
      IF ( NSP .EQ. 0 ) THEN
        WRITE (NDSE,902) 
        RETURN
        END IF
!
! 0.b Test output
!
!/T      WRITE (NDST,9000) NSP, NFR1, NTH1, XF1, FR1, TH1*RADE,       &
!/T                             NFR2, NTH2, XF2, FR2, TH2*RADE, FTL
!
! -------------------------------------------------------------------- /
! 1.  Search stored interpolation data for match
!
      FOUND  = .FALSE.
!
      DO I=1, NCASES
!
        IF ( I .EQ. 1 ) THEN
          CURRENT => IDATA
        ELSE
          CURRENT => CURRENT%NEXT
          END IF
!
!/T1        WRITE (NDST,9010) I, CURRENT%NFR1, CURRENT%NTH1,          &
!/T1                     CURRENT%XF1, CURRENT%FR1, CURRENT%TH1*RADE,  &
!/T1                             CURRENT%NFR2, CURRENT%NTH2,          &
!/T1                     CURRENT%XF2, CURRENT%FR2, CURRENT%TH2*RADE
!
        FOUND = CURRENT%NFR1.EQ.NFR1 .AND. CURRENT%NFR2.EQ.NFR2 .AND. &
                CURRENT%NTH1.EQ.NTH1 .AND. CURRENT%NTH2.EQ.NTH2 .AND. &
                CURRENT%XF1 .EQ.XF1  .AND. CURRENT%XF2 .EQ.XF2  .AND. &
                CURRENT%FR1 .EQ.FR1  .AND. CURRENT%FR2 .EQ.FR2  .AND. &
                CURRENT%TH1 .EQ.TH1  .AND. CURRENT%TH2 .EQ.TH2 
        IF ( FOUND ) EXIT
!
        END DO
!
! -------------------------------------------------------------------- /
! 2.  Link or compute interpolation data
! 2.a Link
!
      IF ( FOUND ) THEN
!
!/T          WRITE (NDST,9020) I
!
        DTH1   => CURRENT%DTH1
        DTH2   => CURRENT%DTH2
        IDTH   => CURRENT%IDTH
        RDTH   => CURRENT%RDTH
!
        FRQ1   => CURRENT%FRQ1
        FRQ2   => CURRENT%FRQ2
        XDF1   => CURRENT%XDF1
        XDF2   => CURRENT%XDF2
        NFR2T  => CURRENT%NFR2T
        IDFR   => CURRENT%IDFR
        RDFR   => CURRENT%RDFR
!
! 2.b Compute
!
      ELSE
!
        NCASES = NCASES + 1
!/T          WRITE (NDST,9021) NCASES
!
! 2.b.1 Point and allocate as necessary
!
        IF ( NCASES .EQ. 1 ) THEN
          ALLOCATE ( IDATA )
          CURRENT => IDATA
        ELSE
          ALLOCATE ( CURRENT%NEXT )
          CURRENT => CURRENT%NEXT
          END IF 
!
! 2.b.2 Store test data
!
        CURRENT%ICASE = NCASES
        CURRENT%NFR1  = NFR1
        CURRENT%NTH1  = NTH1
        CURRENT%XF1   = XF1
        CURRENT%FR1   = FR1
        CURRENT%TH1   = TH1
        CURRENT%NFR2  = NFR2
        CURRENT%NTH2  = NTH2
        CURRENT%XF2   = XF2
        CURRENT%FR2   = FR2
        CURRENT%TH2   = TH2
!
! 2.b.3 Directional redistribution data
!
        DTH1   => CURRENT%DTH1
        DTH1   = TPI / REAL(NTH1)
        DTH2   => CURRENT%DTH2
        DTH2   = TPI / REAL(NTH2)
!
        IF ( DTH1 .LE. DTH2 ) THEN
          NRMAX  = 2
        ELSE
          NRMAX  = 2 + INT(DTH1/DTH2)
          END IF
!
        ALLOCATE (CURRENT%IDTH(0:NRMAX,NTH1),CURRENT%RDTH(NRMAX,NTH1))
        IDTH   => CURRENT%IDTH
        RDTH   => CURRENT%RDTH
        IDTH   = 0
        RDTH   = 0.
!
        DO I=1, NTH1
          LOW    = TH1 + REAL(I-1)*DTH1 - 0.5*DTH1
          HGH    = LOW + DTH1
          RLOW   = 1. + (LOW-TH2)/DTH2
          RHGH   = 1. + (HGH-TH2)/DTH2
          DO J=NINT(RLOW), NINT(RLOW)+NRMAX-1
            BLOW   = TH2 + REAL(J-1)*DTH2 - 0.5*DTH2
            BHGH   = BLOW + DTH2
            FRAC   = (MIN(BHGH,HGH)-MAX(BLOW,LOW)) / (HGH-LOW)
            IF ( FRAC .GT. 1.E-5 ) THEN
               IDTH(0,I) = IDTH(0,I) + 1
               IDTH(IDTH(0,I),I) = 1 + MOD(J-1+NTH2,NTH2)
               RDTH(IDTH(0,I),I) = FRAC
               END IF
            END DO
          END DO
!
! 2.b.4 Frequency redistribution data
!
        ALLOCATE ( CURRENT%FRQ1(NFR1), CURRENT%FRQ2(NFR2) )
        FRQ1   => CURRENT%FRQ1
        FRQ2   => CURRENT%FRQ2
!
        FRQ1(1) = FR1
        DO I=2, NFR1
          FRQ1(I) = XF1 * FRQ1(I-1)
          END DO
!
        FRQ2(1) = FR2
        DO I=2, NFR2
          FRQ2(I) = XF2 * FRQ2(I-1)
          END DO
!
        XDF1   => CURRENT%XDF1
        XDF1   = 0.5 * ( XF1 - 1./XF1 )
        XDF2   => CURRENT%XDF2
        XDF2   = 0.5 * ( XF2 - 1./XF2 )
!
        IF ( XDF1 .LE. XDF2 ) THEN
          NRMAX  = 2
        ELSE
          NRMAX  = 1
          AUX1   = XDF1
          AUX2   = XDF2
          DO
            NRMAX  = NRMAX + 1
            AUX1   = AUX1 - AUX2
            AUX2   = AUX2 / XF2
            IF ( AUX1 .LT. 0. ) EXIT
            END DO
          END IF
!
        ALLOCATE (CURRENT%IDFR(0:NRMAX,NFR1),CURRENT%RDFR(NRMAX,NFR1))
        IDFR   => CURRENT%IDFR
        RDFR   => CURRENT%RDFR
        IDFR   = 0
        RDFR   = 0.
!
        DO I=1, NFR1
          IF ( I .EQ. 1 ) THEN
            HGH    = 0.5 * ( FRQ1(I) + FRQ1(I+1) )
            LOW    = HGH - XDF1*FRQ1(I)
          ELSE
            LOW    = 0.5 * ( FRQ1(I) + FRQ1(I-1) )
            HGH    = LOW + XDF1*FRQ1(I)
            END IF
          DO J=1, NFR2
            IF ( J .EQ. 1 ) THEN
              BHGH   = 0.5 * ( FRQ2(J) + FRQ2(J+1) )
              BLOW   = BHGH - XDF2*FRQ2(J)
            ELSE
              BLOW   = 0.5 * ( FRQ2(J) + FRQ2(J-1) )
              BHGH   = BLOW + XDF2*FRQ2(J)
              END IF
            IF ( BHGH .LE. LOW ) CYCLE
            IF ( BLOW .GE. HGH ) EXIT
            FRAC   = (MIN(BHGH,HGH)-MAX(BLOW,LOW)) / (HGH-LOW)
            IF ( FRAC .LT. 1.E-5 ) CYCLE
            IDFR(0,I) = IDFR(0,I) + 1
            IDFR(IDFR(0,I),I) = J
            RDFR(IDFR(0,I),I) = FRAC
            END DO
          END DO
!
        NFR2T  => CURRENT%NFR2T
        NFR2T  = NFR2 + 1
        DO J=NFR2, 1, -1
          IF ( J .EQ. 1 ) THEN
            BHGH   = 0.5 * ( FRQ2(J) + FRQ2(J+1) )
          ELSE
            BLOW   = 0.5 * ( FRQ2(J) + FRQ2(J-1) )
            BHGH   = BLOW + XDF2*FRQ2(J)
            END IF
          IF ( BHGH .GT. HGH ) THEN
            NFR2T  = J
          ELSE
            EXIT
            END IF
          END DO
!
        END IF
!
! 2.c Test output
!
!/T2      WRITE (NDST,9022)
!/T2      DO I=1, NTH1
!/T2        WRITE (NDST,9024) I, IDTH(0,I),                           &
!/T2                          (IDTH(J,I),RDTH(J,I),J=1,IDTH(0,I))
!/T2        END DO
!/T2      WRITE (NDST,9023) NFR2T
!/T2      DO I=1, NFR1
!/T2        WRITE (NDST,9024) I, IDFR(0,I),                           &
!/T2                          (IDFR(J,I),RDFR(J,I),J=1,IDFR(0,I))
!/T2        END DO
!
! -------------------------------------------------------------------- /
! 3.  Convert
! 3.a Discrete energies
!
!/T      WRITE (NDST,9030)
!
      SP2    = 0.
!
      DO I2=1, NFR1
      DO L2=1, IDFR(0,I2)
        J2   = IDFR(L2,I2)
          R2   = RDFR(L2,I2)
          DO I1=1,NTH1
            DO L1=1, IDTH( 0,I1)
              J1   = IDTH(L1,I1)
              R1   = RDTH(L1,I1)
              FRAC   = R2 * FRQ1(I2) * XDF1 * R1 * DTH1
              SP2(J1,J2,:) = SP2(J1,J2,:) + FRAC * SP1(I1,I2,:)
              END DO
            END DO
          END DO
        END DO
!
! 3.b Energy densities
!
!/T      WRITE (NDST,9031)
!
      DO J2=1, NFR2
        DO J1=1, NTH2
          FACT   = 1. / ( FRQ2(J2) * XDF2 * DTH2 )
          SP2(J1,J2,:) = FACT * SP2(J1,J2,:)
          END DO
        END DO
!
! 3.c Add the tail
!
!/T      WRITE (NDST,9032)
!
      DO J2=NFR2T, NFR2
        SP2(:,J2,:) = FTL * SP2(:,J2-1,:)
        END DO
!
      RETURN
!
! Formats
!
  900 FORMAT (/' *** ERROR W3CSPC: ILLEGAL INPUT PARAMETERS ***'/     &
               '                   INPUT  : ',2I8,2F10.4/             &
               '                   OUTPUT : ',2I8,2F10.4)
  901 FORMAT (/' *** ERROR W3CSPC: NEGATIVE NUMBER OF SPECTRA ***'/)
  902 FORMAT (/' *** WARNING W3CSPC: NO SPECTRA ***'/)
!
!/T 9000 FORMAT ( ' TEST W3CSPC : NR. OF SPECTRA : ',I8/             &
!/T               '               INPUT SPECTRA  : ',2I4,2F8.4,F6.1/ &
!/T               '               OUTPUT SPECTRA : ',2I4,2F8.4,F6.1/ &
!/T               '               TAIL FACTOR    : ',F8.5)
!
!/T1 9010 FORMAT ( ' TEST W3CSPC : TEST INFO CASE : ',I8/             &
!/T1               '               INPUT SPECTRA  : ',2I4,2F8.4,F6.1/ &
!/T1               '               OUTPUT SPECTRA : ',2I4,2F8.4,F6.1)
!
!/T 9020 FORMAT ( ' TEST W3CSPC : USING STORED DATA FOR CASE',I4)
!/T 9021 FORMAT ( ' TEST W3CSPC : COMPUTING DATA FOR CASE',I4)
!/T2 9022 FORMAT ( ' TEST W3CSPC : DIRECTIONAL DISTRIBUTION DATA')
!/T2 9023 FORMAT ( ' TEST W3CSPC : FREQUENCY DISTRIBUTION DATA, ',    &
!/T2                              'TAIL AT',I4)
!/T2 9024 FORMAT ( '           ',I4,I4,' :',10(I4,F5.2) )
!
!/T 9030 FORMAT ( ' TEST W3CSPC : STARTING CONVERSION')
!/T 9031 FORMAT ( ' TEST W3CSPC : ENERGIES TO DENSITIES')
!/T 9032 FORMAT ( ' TEST W3CSPC : ADD TAIL')
!/
!/ End of W3CSPC ----------------------------------------------------- /
!/
      END SUBROUTINE W3CSPC
!/
!/ End of module W3CSPCMD -------------------------------------------- /
!/
      END MODULE W3CSPCMD
