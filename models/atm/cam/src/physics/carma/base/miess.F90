! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This subroutine computes mie scattering by a stratified sphere,
!! i.e. a particle consisting of a spherical core surrounded by a
!! Spherical shell.  The basic code used was that described in the
!! report: " Subroutines for computing the parameters of the
!! electromagnetic radiation scattered by a sphere " J.V. Dave,
!! IBM Scientific Center, Palo Alto , California.
!! Report No. 320 - 3236 .. May 1968 .
!!
!! The modifications for stratified spheres are described in
!!   Toon and Ackerman, Appl. Optics, in press, 1981
!!
!! The definitions for the output parameters can be found in "Light
!! scattering by small particles, H.C.Van de Hulst, John Wiley '
!! Sons, Inc., New York, 1957".
!!
!! Also the subroutine computes the capital A function by making use of
!! downward recurrence relationship.
!!
!! @author Brian Toon
!! @version 1981?
SUBROUTINE miess(carma,RO,RFR,RFI,THETD,JX,QEXT,QSCAT,QBS,CTBRQS,R,RE2,TMAG2,WVNO,rc)

	! types
  use carma_precision_mod
  use carma_constants_mod, only : IT, DEG2RAD
  use carma_enums_mod, only     : RC_ERROR
  use carma_types_mod, only     : carma_type
  use carma_mod

	implicit none

  type(carma_type), intent(in)         :: carma   !! the carma object
  real(kind=f), intent(in)             :: RO         !! OUTER (SHELL) RADIUS
  real(kind=f), intent(in)             :: RFR        !! REAL PART OF THE SHELL INDEX OF REFRACTION
  real(kind=f), intent(in)             :: RFI        !! IMAGINARY PART OF THE SHELL INDEX OF REFRACTION
  real(kind=f), intent(in)             :: R          !! CORE RADIUS
  real(kind=f), intent(in)             :: RE2        !! REAL PART OF THE CORE INDEX OF REFRACTION
  real(kind=f), intent(in)             :: TMAG2      !! IMAGINARY PART OF THE CORE INDEX OF REFRACTION

  !! ANGLE IN DEGREES BETWEEN THE DIRECTIONS OF THE INCIDENT
  !! AND THE SCATTERED RADIATION.  THETD(J) IS< OR= 90.0
  !! IF THETD(J) SHOULD HAPPEN TO BE GREATER THAN 90.0, ENTER WITH
  !! SUPPLEMENTARY VALUE, SEE COMMENTS ON ELTRMX.
  real(kind=f), intent(inout)          :: THETD(IT)
  
  !! TOTAL NUMBER OF THETD FOR WHICH THE COMPUTATIONS ARE
  !! REQUIRED.  JX SHOULD NOT EXCEED IT UNLESS THE DIMENSIONS
  !! STATEMENTS ARE APPROPRIATEDLY MODIFIED.
  integer, intent(in)                  :: JX
  
  real(kind=f), intent(out)            :: QEXT       !! EFFICIENCY FACTOR FOR EXTINCTION,VAN DE HULST,P.14 ' 127.
  real(kind=f), intent(out)            :: QSCAT      !! EFFICIENCY FACTOR FOR SCATTERING,V.D. HULST,P.14 ' 127.
  real(kind=f), intent(out)            :: QBS        !! BACK SCATTER CROSS SECTION.
  real(kind=f), intent(out)            :: CTBRQS     !! AVERAGE(COSINE THETA) * QSCAT,VAN DE HULST,P.128.
  real(kind=f), intent(in)             :: WVNO       !! 2*PI / WAVELENGTH.
  integer, intent(inout)               :: rc         !! return code, negative indicates failure

  ! Local declarations
  real(kind=f), parameter              :: EPSILON_MIE = 1.e-14_f
  
  integer                              :: I
  integer                              :: J
  integer                              :: K
  integer                              :: M
  integer                              :: N
  integer                              :: NN
  integer                              :: NMX1
  integer                              :: NMX2
  integer                              :: IFLAG
  integer                              :: IACAP
  
  complex(kind=f) :: FNAP,   FNBP,   &
                     FNA,    FNB,    RF,       RRF, &
                     RRFX,   WM1,    FN1,      FN2, &
                     TC1,    TC2,    WFN(2),   Z(4), &
                     K1,     K2,     K3,       &
                     RCR,    U(8),   DH1, &
                     DH2,    DH4,    P24H24,   P24H21, &
                     PSTORE, HSTORE, DUMMY,    DUMSQ

  complex(kind=f), allocatable :: ACAP(:), W(:,:)

 
  ! TA(1): REAL PART OF WFN(1).  TA(2): IMAGINARY PART OF WFN(1).
  ! TA(3): REAL PART OF WFN(2).  TA(4): IMAGINARY PART OF WFN(2).
  ! TB(1): REAL PART OF FNA.     TB(2): IMAGINARY PART OF FNA.
  ! TC(1): REAL PART OF FNB.     TC(2): IMAGINARY PART OF FNB.
  ! TD(1): REAL PART OF FNAP.    TD(2): IMAGINARY PART OF FNAP.
  ! TE(1): REAL PART OF FNBP.    TE(2): IMAGINARY PART OF FNBP.
  ! FNAP, FNBP  ARE THE PRECEDING VALUES OF FNA, FNB RESPECTIVELY.
  real(kind=f) :: T(5), TA(4), TB(2), TC(2), TD(2), TE(2), &
                  PI(3,IT), TAU(3,IT), CSTHT(IT), SI2THT(IT), &  
                  X, X1, X4, Y1, Y4, RX, SINX1, SINX4, COSX1, COSX4, &
                  EY1, E2Y1, EY4, EY1MY4, EY1PY4, AA, BB, CC, DD, DENOM, &
                  REALP, AMAGP, QBSR, QBSI, RMM, PIG, RXP4

  equivalence (FNA,TB(1)),(FNB,TC(1)),(FNAP,TD(1)),(FNBP,TE(1))
  
  !! ELTRMX(I,J,K): ELEMENTS OF THE TRANSFORMATION MATRIX F,V.D.HULST,P.34,45 ' 125.
  !!   I=1: ELEMENT M SUB 2..I=2: ELEMENT M SUB 1..
  !!   I = 3: ELEMENT S SUB 21.. I = 4: ELEMENT D SUB 21..
  !! ELTRMX(I,J,1) REPRESENTS THE ITH ELEMENT OF THE MATRIX FOR
  !!   THE ANGLE THETD(J).. ELTRMX(I,J,2) REPRESENTS THE ITH ELEMENT
  !!   OF THE MATRIX FOR THE ANGLE 180.0 - THETD(J) ..
  real(kind=f)                         :: ELTRMX(4,IT,2)


  ! IF THE CORE IS SMALL SCATTERING IS COMPUTED FOR THE SHELL ONLY
  IFLAG = 1
  if ( R/RO .LT. 1.e-6_f ) IFLAG = 2
  
  if ( JX .gt. IT ) then
    if (do_print) write(LUNOPRT, '(a,i3,a)') "miess:: The value of the argument JX=", JX, " is greater than IT."
    rc = RC_ERROR
    return
  endif
      
  RF =  CMPLX( RFR,  -RFI )
  RCR =  CMPLX( RE2, -TMAG2 )
  X  =  RO * WVNO
  K1 =  RCR * WVNO
  K2 =  RF * WVNO
  K3 =  CMPLX( WVNO, 0.0_f )
  Z(1) =  K2 * RO
  Z(2) =  K3 * RO
  Z(3) =  K1 * R
  Z(4) =  K2 * R
  X1   =  REAL( Z(1) )
  X4   =  REAL( Z(4) )
  Y1   =  aimag( Z(1) )
  Y4   =  aimag( Z(4) )
  RRF  =  1.0_f / RF
  RX   =  1.0_f / X
  RRFX =  RRF * RX
  T(1) =  ( X**2 ) * ( RFR**2 + RFI**2 )
  T(1) =  SQRT( T(1) )
  NMX1 =  1.10_f * T(1)
    
  ! The dimension of ACAP.
  !
  ! In the original program the dimension of ACAP was 7000.
  ! For conserving space this should be not much higher than
  ! The value, NMX1=1.1*(NREAL**2 + NIMAG**2)**.5 * X + 1  
  IACAP = max(7000, int(1.5_f * NMX1))
  allocate(ACAP(IACAP))
  allocate(W(3,IACAP))
      
  NMX2 = T(1)
 
  if ( NMX1 .le. 150 ) then
    NMX1 = 150
    NMX2 = 135
  endif

  ACAP( NMX1+1 )  =  ( 0.0_f, 0.0_f )
 
  if ( IFLAG .ne. 2 ) then
    do N = 1,3
      W( N,NMX1+1 )  =  ( 0.0_f, 0.0_f )
    enddo
  endif
 
  do N = 1,NMX1
    NN = NMX1 - N + 1
    ACAP(NN) = (NN+1) * RRFX - 1.0_f / ( (NN+1) * RRFX + ACAP(NN+1) )
    if ( IFLAG .ne. 2 ) then
      do M = 1,3
        W( M,NN ) = (NN+1) / Z(M+1)  - &
                    1.0_f / (  (NN+1) / Z(M+1)  +  W( M,NN+1 )  )
      enddo
    endif
  enddo

  do J = 1,JX
    if ( THETD(J) .lt. 0.0 )  THETD(J) =  ABS( THETD(J) )

    if ( THETD(J) .le. 0.0 )  then
      CSTHT(J)  = 1.0_f
      SI2THT(J) = 0.0_f
    else if ( THETD(J) .lt. 90.0_f ) then
      T(1)      =  THETD(J) * DEG2RAD
      CSTHT(J)  =  COS( T(1) )
      SI2THT(J) =  1.0_f - CSTHT(J)**2
    else if ( THETD(J) .le. 90.0_f ) then
      CSTHT(J)  =  0.0_f
      SI2THT(J) =  1.0_f
    else 
      if (do_print) write(LUNOPRT, '(a,i3)') "miess:: The value of the scattering angle is greater than 90.0 Degrees. It is .", THETD(J)
      rc = RC_ERROR
      return
    end if
  enddo
    
  do J = 1,JX
    PI(1,J)  =  0.0_f
    PI(2,J)  =  1.0_f
    TAU(1,J) =  0.0_f
    TAU(2,J) =  CSTHT(J)
  enddo

  ! INITIALIZATION OF HOMOGENEOUS SPHERE
  T(1)   =  COS(X)
  T(2)   =  SIN(X)
  WM1    =  CMPLX( T(1),-T(2) )
  WFN(1) =  CMPLX( T(2), T(1) )
  TA(1)  =  T(2)
  TA(2)  =  T(1)
  WFN(2) =  RX * WFN(1) - WM1
  TA(3)  =  REAL(WFN(2))
  TA(4)  =  aimag(WFN(2))

  if ( IFLAG .ne. 2 ) then
    N = 1

    ! INITIALIZATION PROCEDURE FOR STRATIFIED SPHERE BEGINS HERE
    SINX1   =  SIN( X1 )
    SINX4   =  SIN( X4 )
    COSX1   =  COS( X1 )
    COSX4   =  COS( X4 )
    EY1     =  EXP( Y1 )
    E2Y1    =  EY1 * EY1
    EY4     =  EXP( Y4 )
    EY1MY4  =  EXP( Y1 - Y4 )
    EY1PY4  =  EY1 * EY4
    EY1MY4  =  EXP( Y1 - Y4 )
    AA  =  SINX4 * ( EY1PY4 + EY1MY4 )
    BB  =  COSX4 * ( EY1PY4 - EY1MY4 )
    CC  =  SINX1 * ( E2Y1 + 1.0 )
    DD  =  COSX1 * ( E2Y1 - 1.0 )
    DENOM   =  1.0_f  +  E2Y1 * ( 4.0_f * SINX1 * SINX1 - 2.0_f + E2Y1 )
    REALP   =  ( AA * CC  +  BB * DD ) / DENOM
    AMAGP   =  ( BB * CC  -  AA * DD ) / DENOM
    DUMMY   =  CMPLX( REALP, AMAGP )
    AA  =  SINX4 * SINX4 - 0.5_f
    BB  =  COSX4 * SINX4
    P24H24  =  0.5_f + CMPLX( AA,BB ) * EY4 * EY4
    AA  =  SINX1 * SINX4  -  COSX1 * COSX4
    BB  =  SINX1 * COSX4  +  COSX1 * SINX4
    CC  =  SINX1 * SINX4  +  COSX1 * COSX4
    DD  = -SINX1 * COSX4  +  COSX1 * SINX4
    P24H21  =  0.5_f * CMPLX( AA,BB ) * EY1 * EY4  + 0.5_f * CMPLX( CC,DD ) * EY1MY4
    DH4  =  Z(4) / ( 1.0_f + ( 0.0_f, 1.0_f ) * Z(4) )  -  1.0_f / Z(4)
    DH1  =  Z(1) / ( 1.0_f + ( 0.0_f, 1.0_f ) * Z(1) )  -  1.0_f / Z(1)
    DH2  =  Z(2) / ( 1.0_f + ( 0.0_f, 1.0_f ) * Z(2) )  -  1.0_f / Z(2)
    PSTORE  =  ( DH4 + N / Z(4) )  *  ( W(3,N) + N / Z(4) )
    P24H24  =  P24H24 / PSTORE
    HSTORE  =  ( DH1 + N / Z(1) )  *  ( W(3,N) + N / Z(4) )
    P24H21  =  P24H21 / HSTORE
    PSTORE  =  ( ACAP(N) + N / Z(1) )  /  ( W(3,N) + N / Z(4) )
    DUMMY   =  DUMMY * PSTORE
    DUMSQ   =  DUMMY * DUMMY

    ! NOTE:  THE DEFINITIONS OF U(I) IN THIS PROGRAM ARE NOT THE SAME AS
    !        THE USUBI DEFINED IN THE ARTICLE BY TOON AND ACKERMAN.  THE
    !        CORRESPONDING TERMS ARE:
    !          USUB1 = U(1)                       USUB2 = U(5)
    !          USUB3 = U(7)                       USUB4 = DUMSQ
    !          USUB5 = U(2)                       USUB6 = U(3)
    !          USUB7 = U(6)                       USUB8 = U(4)
    !          RATIO OF SPHERICAL BESSEL FTN TO SPHERICAL HENKAL FTN = U(8)

    U(1) =  K3 * ACAP(N)  -  K2 * W(1,N)
    U(2) =  K3 * ACAP(N)  -  K2 * DH2
    U(3) =  K2 * ACAP(N)  -  K3 * W(1,N)
    U(4) =  K2 * ACAP(N)  -  K3 * DH2
    U(5) =  K1 *  W(3,N)  -  K2 * W(2,N)
    U(6) =  K2 *  W(3,N)  -  K1 * W(2,N)
    U(7) =  ( 0.0_f, -1.0_f )  *  ( DUMMY * P24H21 - P24H24 )
    U(8) =  TA(3) / WFN(2)

    FNA  =  U(8) * ( U(1)*U(5)*U(7)  +  K1*U(1)  -  DUMSQ*K3*U(5) ) / &
                   ( U(2)*U(5)*U(7)  +  K1*U(2)  -  DUMSQ*K3*U(5) )
    FNB  =  U(8) * ( U(3)*U(6)*U(7)  +  K2*U(3)  -  DUMSQ*K2*U(6) ) / &
                   ( U(4)*U(6)*U(7)  +  K2*U(4)  -  DUMSQ*K2*U(6) )
  else
    TC1  =  ACAP(1) * RRF  +  RX
    TC2  =  ACAP(1) * RF   +  RX
    FNA  =  ( TC1 * TA(3)  -  TA(1) ) / ( TC1 * WFN(2)  -  WFN(1) )
    FNB  =  ( TC2 * TA(3)  -  TA(1) ) / ( TC2 * WFN(2)  -  WFN(1) )
  endif

  FNAP = FNA
  FNBP = FNB
  T(1) = 1.50_f

  ! FROM HERE TO THE STATMENT NUMBER 90, ELTRMX(I,J,K) HAS
  ! FOLLOWING MEANING:
  ! ELTRMX(1,J,K): REAL PART OF THE FIRST COMPLEX AMPLITUDE.
  ! ELTRMX(2,J,K): IMAGINARY PART OF THE FIRST COMPLEX AMPLITUDE.
  ! ELTRMX(3,J,K): REAL PART OF THE SECOND COMPLEX AMPLITUDE.
  ! ELTRMX(4,J,K): IMAGINARY PART OF THE SECOND COMPLEX AMPLITUDE.
  ! K = 1 : FOR THETD(J) AND K = 2 : FOR 180.0 - THETD(J)
  ! DEFINITION OF THE COMPLEX AMPLITUDE: VAN DE HULST,P.125.
  TB(1) = T(1) * TB(1)
  TB(2) = T(1) * TB(2)
  TC(1) = T(1) * TC(1)
  TC(2) = T(1) * TC(2)
  
  do J = 1,JX
    ELTRMX(1,J,1) = TB(1) * PI(2,J) + TC(1) * TAU(2,J)
    ELTRMX(2,J,1) = TB(2) * PI(2,J) + TC(2) * TAU(2,J)
    ELTRMX(3,J,1) = TC(1) * PI(2,J) + TB(1) * TAU(2,J)
    ELTRMX(4,J,1) = TC(2) * PI(2,J) + TB(2) * TAU(2,J)
    ELTRMX(1,J,2) = TB(1) * PI(2,J) - TC(1) * TAU(2,J)
    ELTRMX(2,J,2) = TB(2) * PI(2,J) - TC(2) * TAU(2,J)
    ELTRMX(3,J,2) = TC(1) * PI(2,J) - TB(1) * TAU(2,J)
    ELTRMX(4,J,2) = TC(2) * PI(2,J) - TB(2) * TAU(2,J)
  enddo

  QEXT   = 2.0_f * ( TB(1) + TC(1))
  QSCAT  = ( TB(1)**2 + TB(2)**2 + TC(1)**2 + TC(2)**2 ) / 0.75_f
  CTBRQS = 0.0_f
  QBSR   = -2.0_f*(TC(1) - TB(1))
  QBSI   = -2.0_f*(TC(2) - TB(2))
  RMM    = -1.0_f
  N = 2
  
  ! Iterate until the answer converges.
  T(4) = EPSILON_MIE
  
  do while ( T(4) .ge. EPSILON_MIE )
  
    T(1) = 2*N - 1
    T(2) =   N - 1
    T(3) = 2*N + 1
    
    do J = 1,JX
      PI(3,J)  = ( T(1) * PI(2,J) * CSTHT(J) - N * PI(1,J) ) / T(2)
      TAU(3,J) = CSTHT(J) * ( PI(3,J) - PI(1,J) )  - &
                T(1) * SI2THT(J) * PI(2,J)  +  TAU(1,J)
    end do
  
    ! HERE SET UP HOMOGENEOUS SPHERE
    WM1    =  WFN(1)
    WFN(1) =  WFN(2)
    TA(1)  =  REAL(WFN(1))
    TA(2)  =  aimag(WFN(1))
    TA(4)  =  aimag(WFN(2))
    WFN(2) =  T(1) * RX * WFN(1)  -  WM1
    TA(3)  =  REAL(WFN(2))
  
    if ( IFLAG .ne. 2 ) then
  
      ! HERE SET UP STRATIFIED SPHERE
      DH2  =  - N / Z(2)  +  1.0_f / ( N / Z(2) - DH2 )
      DH4  =  - N / Z(4)  +  1.0_f / ( N / Z(4) - DH4 )
      DH1  =  - N / Z(1)  +  1.0_f / ( N / Z(1) - DH1 )
      PSTORE  =  ( DH4 + N / Z(4) )  *  ( W(3,N) + N / Z(4) )
      P24H24  =  P24H24 / PSTORE
      HSTORE  =  ( DH1 + N / Z(1) )  *  ( W(3,N) + N / Z(4) )
      P24H21  =  P24H21 / HSTORE
      PSTORE  =  ( ACAP(N) + N / Z(1) )  /  ( W(3,N) + N / Z(4) )
      DUMMY   =  DUMMY * PSTORE
      DUMSQ   =  DUMMY * DUMMY
  
      U(1) =  K3 * ACAP(N)  -  K2 * W(1,N)
      U(2) =  K3 * ACAP(N)  -  K2 * DH2
      U(3) =  K2 * ACAP(N)  -  K3 * W(1,N)
      U(4) =  K2 * ACAP(N)  -  K3 * DH2
      U(5) =  K1 *  W(3,N)  -  K2 * W(2,N)
      U(6) =  K2 *  W(3,N)  -  K1 * W(2,N)
      U(7) =  ( 0.0_f, -1.0_f )  *  ( DUMMY * P24H21 - P24H24 )
      U(8) =  TA(3) / WFN(2)
  
      FNA  =  U(8) * ( U(1)*U(5)*U(7)  +  K1*U(1)  -  DUMSQ*K3*U(5) ) / &
                     ( U(2)*U(5)*U(7)  +  K1*U(2)  -  DUMSQ*K3*U(5) )
      FNB  =  U(8) * ( U(3)*U(6)*U(7)  +  K2*U(3)  -  DUMSQ*K2*U(6) ) / &
                     ( U(4)*U(6)*U(7)  +  K2*U(4)  -  DUMSQ*K2*U(6) )
    endif
    
    TC1  =  ACAP(N) * RRF  +  N * RX
    TC2  =  ACAP(N) * RF   +  N * RX
    FN1  =  ( TC1 * TA(3)  -  TA(1) ) /  ( TC1 * WFN(2) - WFN(1) )
    FN2  =  ( TC2 * TA(3)  -  TA(1) ) /  ( TC2 * WFN(2) - WFN(1) )
    M    =  WVNO * R
    
    if ( N .ge. M ) then
      if ( IFLAG .ne. 2 ) then
        if ( abs(  ( FN1-FNA ) / FN1  )  .LT.  EPSILON_MIE .AND. &
             abs(  ( FN2-FNB ) / FN2  )  .LT. EPSILON_MIE  ) IFLAG = 2
             
        if ( IFLAG .ne. 1 ) then
          FNA  =  FN1
          FNB  =  FN2
        endif
      else
        FNA  =  FN1
        FNB  =  FN2
     endif
    endif
    
    T(5)  =  N
    T(4)  =  T(1) / ( T(5) * T(2) )
    T(2)  =  (  T(2) * ( T(5) + 1.0_f )  ) / T(5)
  
    CTBRQS  =  CTBRQS  +  T(2) * ( TD(1) * TB(1)  +  TD(2) * TB(2) &
                      +           TE(1) * TC(1)  +  TE(2) * TC(2) ) &
                      +  T(4) * ( TD(1) * TE(1)  +  TD(2) * TE(2) )
    QEXT    =   QEXT  +  T(3) * ( TB(1) + TC(1) )
    
    !     $        T(3), TB(1), TC(1), QEXT
    T(4)    =  TB(1)**2 + TB(2)**2 + TC(1)**2 + TC(2)**2
    QSCAT   =  QSCAT  +  T(3) * T(4)
    RMM     =  -RMM
    QBSR    =  QBSR + T(3)*RMM*(TC(1) - TB(1))
    QBSI    =  QBSI + T(3)*RMM*(TC(2) - TB(2))
  
    T(2)    =  N * (N+1)
    T(1)    =  T(3) / T(2)
    K = (N/2)*2
    
    do J = 1,JX
      ELTRMX(1,J,1) = ELTRMX(1,J,1)+T(1)*(TB(1)*PI(3,J)+TC(1)*TAU(3,J))
      ELTRMX(2,J,1) = ELTRMX(2,J,1)+T(1)*(TB(2)*PI(3,J)+TC(2)*TAU(3,J))
      ELTRMX(3,J,1) = ELTRMX(3,J,1)+T(1)*(TC(1)*PI(3,J)+TB(1)*TAU(3,J))
      ELTRMX(4,J,1) = ELTRMX(4,J,1)+T(1)*(TC(2)*PI(3,J)+TB(2)*TAU(3,J))
      
      IF ( K .EQ. N )  THEN
       ELTRMX(1,J,2) = ELTRMX(1,J,2)+T(1)*(-TB(1)*PI(3,J)+TC(1)*TAU(3,J))
       ELTRMX(2,J,2) = ELTRMX(2,J,2)+T(1)*(-TB(2)*PI(3,J)+TC(2)*TAU(3,J))
       ELTRMX(3,J,2) = ELTRMX(3,J,2)+T(1)*(-TC(1)*PI(3,J)+TB(1)*TAU(3,J))
       ELTRMX(4,J,2) = ELTRMX(4,J,2)+T(1)*(-TC(2)*PI(3,J)+TB(2)*TAU(3,J))
      ELSE
       ELTRMX(1,J,2) = ELTRMX(1,J,2)+T(1)*(TB(1)*PI(3,J)-TC(1)*TAU(3,J))
       ELTRMX(2,J,2) = ELTRMX(2,J,2)+T(1)*(TB(2)*PI(3,J)-TC(2)*TAU(3,J))
       ELTRMX(3,J,2) = ELTRMX(3,J,2)+T(1)*(TC(1)*PI(3,J)-TB(1)*TAU(3,J))
       ELTRMX(4,J,2) = ELTRMX(4,J,2)+T(1)*(TC(2)*PI(3,J)-TB(2)*TAU(3,J))
      END IF
    enddo
  
    if ( T(4) .ge. EPSILON_MIE ) then
      N = N + 1
      
      do J = 1,JX
        PI(1,J)   =   PI(2,J)
        PI(2,J)   =   PI(3,J)
        TAU(1,J)  =  TAU(2,J)
        TAU(2,J)  =  TAU(3,J)
      enddo
      
      FNAP  =  FNA
      FNBP  =  FNB
      
      if ( N .gt. NMX2 ) then
        if (do_print) write(LUNOPRT, '(a)') "miess:: The upper limit for acap is not enough."
        rc = RC_ERROR
        return
      endif
    endif
  enddo
  
  ! Calculate the results.
  do J = 1,JX
    do K = 1,2
      do I= 1,4
        T(I)  =  ELTRMX(I,J,K)
      enddo
      
      ELTRMX(2,J,K)  =  T(1)**2  +  T(2)**2
      ELTRMX(1,J,K)  =  T(3)**2  +  T(4)**2
      ELTRMX(3,J,K)  =  T(1) * T(3)  +  T(2) * T(4)
      ELTRMX(4,J,K)  =  T(2) * T(3)  -  T(4) * T(1)
    enddo
  enddo
  
  T(1)    = 2.0_f * RX**2
  QEXT    = QEXT * T(1)
  QSCAT   = QSCAT * T(1)
  CTBRQS  = 2.0_f * CTBRQS * T(1)

  ! QBS IS THE BACK SCATTER CROSS SECTION
  PIG   = ACOS(-1.0_f)
  RXP4  = RX*RX/(4.0_f*PIG)
  QBS   = RXP4*(QBSR**2 + QBSI**2)

  deallocate(ACAP)
  deallocate(W)

  return
end
