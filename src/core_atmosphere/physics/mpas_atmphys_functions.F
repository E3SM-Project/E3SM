! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_functions

 implicit none
 private
 public:: gammln,gammp,wgamma,rslf,rsif


 contains


!=================================================================================================================
!NOTE: functions rslf and rsif are taken from module_mp_thompson temporarily for computing
!      the diagnostic relative humidity. These two functions will be removed from this module
!      when the Thompson cloud microphysics scheme will be restored to MPAS-Dev.
!      Laura D. Fowler (birch.mmm.ucar.edu) / 2013-07-11.

!+---+-----------------------------------------------------------------+
      SUBROUTINE GCF(GAMMCF,A,X,GLN)
!     --- RETURNS THE INCOMPLETE GAMMA FUNCTION Q(A,X) EVALUATED BY ITS
!     --- CONTINUED FRACTION REPRESENTATION AS GAMMCF.  ALSO RETURNS
!     --- LN(GAMMA(A)) AS GLN.  THE CONTINUED FRACTION IS EVALUATED BY
!     --- A MODIFIED LENTZ METHOD.
!     --- USES GAMMLN
      IMPLICIT NONE
      INTEGER, PARAMETER:: ITMAX=100
      REAL, PARAMETER:: gEPS=3.E-7
      REAL, PARAMETER:: FPMIN=1.E-30
      REAL, INTENT(IN):: A, X
      REAL:: GAMMCF,GLN
      INTEGER:: I
      REAL:: AN,B,C,D,DEL,H
      GLN=GAMMLN(A)
      B=X+1.-A
      C=1./FPMIN
      D=1./B
      H=D
      DO 11 I=1,ITMAX
        AN=-I*(I-A)
        B=B+2.
        D=AN*D+B
        IF(ABS(D).LT.FPMIN)D=FPMIN
        C=B+AN/C
        IF(ABS(C).LT.FPMIN)C=FPMIN
        D=1./D
        DEL=D*C
        H=H*DEL
        IF(ABS(DEL-1.).LT.gEPS)GOTO 1
 11   CONTINUE
      PRINT *, 'A TOO LARGE, ITMAX TOO SMALL IN GCF'
 1    GAMMCF=EXP(-X+A*LOG(X)-GLN)*H
      END SUBROUTINE GCF
!  (C) Copr. 1986-92 Numerical Recipes Software 2.02
!+---+-----------------------------------------------------------------+
      SUBROUTINE GSER(GAMSER,A,X,GLN)
!     --- RETURNS THE INCOMPLETE GAMMA FUNCTION P(A,X) EVALUATED BY ITS
!     --- ITS SERIES REPRESENTATION AS GAMSER.  ALSO RETURNS LN(GAMMA(A)) 
!     --- AS GLN.
!     --- USES GAMMLN
      IMPLICIT NONE
      INTEGER, PARAMETER:: ITMAX=100
      REAL, PARAMETER:: gEPS=3.E-7
      REAL, INTENT(IN):: A, X
      REAL:: GAMSER,GLN
      INTEGER:: N
      REAL:: AP,DEL,SUM
      GLN=GAMMLN(A)
      IF(X.LE.0.)THEN
        IF(X.LT.0.) PRINT *, 'X < 0 IN GSER'
        GAMSER=0.
        RETURN
      ENDIF
      AP=A
      SUM=1./A
      DEL=SUM
      DO 11 N=1,ITMAX
        AP=AP+1.
        DEL=DEL*X/AP
        SUM=SUM+DEL
        IF(ABS(DEL).LT.ABS(SUM)*gEPS)GOTO 1
 11   CONTINUE
      PRINT *,'A TOO LARGE, ITMAX TOO SMALL IN GSER'
 1    GAMSER=SUM*EXP(-X+A*LOG(X)-GLN)
      END SUBROUTINE GSER
!  (C) Copr. 1986-92 Numerical Recipes Software 2.02
!+---+-----------------------------------------------------------------+
      REAL FUNCTION GAMMLN(XX)
!     --- RETURNS THE VALUE LN(GAMMA(XX)) FOR XX > 0.
      IMPLICIT NONE
      REAL, INTENT(IN):: XX
      DOUBLE PRECISION, PARAMETER:: STP = 2.5066282746310005D0
      DOUBLE PRECISION, DIMENSION(6), PARAMETER:: &
               COF = (/76.18009172947146D0, -86.50532032941677D0, &
                       24.01409824083091D0, -1.231739572450155D0, &
                      .1208650973866179D-2, -.5395239384953D-5/)
      DOUBLE PRECISION:: SER,TMP,X,Y
      INTEGER:: J

      X=XX
      Y=X
      TMP=X+5.5D0
      TMP=(X+0.5D0)*LOG(TMP)-TMP
      SER=1.000000000190015D0
      DO 11 J=1,6
        Y=Y+1.D0
        SER=SER+COF(J)/Y
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER/X)
      END FUNCTION GAMMLN
!  (C) Copr. 1986-92 Numerical Recipes Software 2.02
!+---+-----------------------------------------------------------------+
      REAL FUNCTION GAMMP(A,X)
!     --- COMPUTES THE INCOMPLETE GAMMA FUNCTION P(A,X)
!     --- SEE ABRAMOWITZ AND STEGUN 6.5.1
!     --- USES GCF,GSER
      IMPLICIT NONE
      REAL, INTENT(IN):: A,X
      REAL:: GAMMCF,GAMSER,GLN
      GAMMP = 0.
      IF((X.LT.0.) .OR. (A.LE.0.)) THEN
        PRINT *, 'BAD ARGUMENTS IN GAMMP'
        RETURN
      ELSEIF(X.LT.A+1.)THEN
        CALL GSER(GAMSER,A,X,GLN)
        GAMMP=GAMSER
      ELSE
        CALL GCF(GAMMCF,A,X,GLN)
        GAMMP=1.-GAMMCF
      ENDIF
      END FUNCTION GAMMP
!  (C) Copr. 1986-92 Numerical Recipes Software 2.02
!+---+-----------------------------------------------------------------+
      REAL FUNCTION WGAMMA(y)

      IMPLICIT NONE
      REAL, INTENT(IN):: y

      WGAMMA = EXP(GAMMLN(y))

      END FUNCTION WGAMMA
!+---+-----------------------------------------------------------------+
! THIS FUNCTION CALCULATES THE LIQUID SATURATION VAPOR MIXING RATIO AS
! A FUNCTION OF TEMPERATURE AND PRESSURE
!
      REAL FUNCTION RSLF(P,T)

      IMPLICIT NONE
      REAL, INTENT(IN):: P, T
      REAL:: ESL,X
      REAL, PARAMETER:: C0= .611583699E03
      REAL, PARAMETER:: C1= .444606896E02
      REAL, PARAMETER:: C2= .143177157E01
      REAL, PARAMETER:: C3= .264224321E-1
      REAL, PARAMETER:: C4= .299291081E-3
      REAL, PARAMETER:: C5= .203154182E-5
      REAL, PARAMETER:: C6= .702620698E-8
      REAL, PARAMETER:: C7= .379534310E-11
      REAL, PARAMETER:: C8=-.321582393E-13

      X=MAX(-80.,T-273.16)

!     ESL=612.2*EXP(17.67*X/(T-29.65))
      ESL=C0+X*(C1+X*(C2+X*(C3+X*(C4+X*(C5+X*(C6+X*(C7+X*C8)))))))
      RSLF=.622*ESL/(P-ESL)

!    ALTERNATIVE
!  ; Source: Murphy and Koop, Review of the vapour pressure of ice and
!             supercooled water for atmospheric applications, Q. J. R.
!             Meteorol. Soc (2005), 131, pp. 1539-1565.
!    ESL = EXP(54.842763 - 6763.22 / T - 4.210 * ALOG(T) + 0.000367 * T
!        + TANH(0.0415 * (T - 218.8)) * (53.878 - 1331.22
!        / T - 9.44523 * ALOG(T) + 0.014025 * T))

      END FUNCTION RSLF
!+---+-----------------------------------------------------------------+
! THIS FUNCTION CALCULATES THE ICE SATURATION VAPOR MIXING RATIO AS A
! FUNCTION OF TEMPERATURE AND PRESSURE
!
      REAL FUNCTION RSIF(P,T)

      IMPLICIT NONE
      REAL, INTENT(IN):: P, T
      REAL:: ESI,X
      REAL, PARAMETER:: C0= .609868993E03
      REAL, PARAMETER:: C1= .499320233E02
      REAL, PARAMETER:: C2= .184672631E01
      REAL, PARAMETER:: C3= .402737184E-1
      REAL, PARAMETER:: C4= .565392987E-3
      REAL, PARAMETER:: C5= .521693933E-5
      REAL, PARAMETER:: C6= .307839583E-7
      REAL, PARAMETER:: C7= .105785160E-9
      REAL, PARAMETER:: C8= .161444444E-12

      X=MAX(-80.,T-273.16)
      ESI=C0+X*(C1+X*(C2+X*(C3+X*(C4+X*(C5+X*(C6+X*(C7+X*C8)))))))
      RSIF=.622*ESI/(P-ESI)

!    ALTERNATIVE
!  ; Source: Murphy and Koop, Review of the vapour pressure of ice and
!             supercooled water for atmospheric applications, Q. J. R.
!             Meteorol. Soc (2005), 131, pp. 1539-1565.
!     ESI = EXP(9.550426 - 5723.265/T + 3.53068*ALOG(T) - 0.00728332*T)

      END FUNCTION RSIF

!=================================================================================================================
 end module mpas_atmphys_functions
!=================================================================================================================
