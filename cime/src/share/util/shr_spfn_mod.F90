
! Define flags for compilers supporting Fortran 2008 intrinsics
! HAVE_GAMMA_INTRINSICS: gamma and log_gamma
! HAVE_ERF_INTRINSICS: erf, erfc, and erfc_scaled
! erfc_scaled(x) = (exp(x**2)*erfc(x))

! Use this flag for compilers that don't have real intrinsics, but link in
! a library for you.
! HAVE_ERF_EXTERNALS: erf and erfc

! These compilers have the intrinsics.
! Intel also has them (and Cray), but as of mid-2015, our implementation is
! actually faster, in part because they do not properly vectorize, so we
! pretend that the compiler version doesn't exist.
#if defined CPRIBM || defined __GFORTRAN__
#define HAVE_GAMMA_INTRINSICS
#define HAVE_ERF_INTRINSICS
#endif

! PGI has external erf/derf and erfc/derfc, and will link them for you, but
! it does not consider them "intrinsics" right now.
#if defined CPRPGI
#define HAVE_ERF_EXTERNALS
#endif

! As of 5.3.1, NAG does not have any of these.

module shr_spfn_mod
! Module for common mathematical functions

! This #ifdef is to allow the module to be compiled with no dependencies,
! even on shr_kind_mod.
#ifndef NO_CSM_SHARE
use shr_kind_mod, only: &
     r4 => shr_kind_r4, &
     r8 => shr_kind_r8
use shr_const_mod, only: &
     pi => shr_const_pi
#endif

implicit none
private
save

#ifdef NO_CSM_SHARE
integer, parameter :: r4 = selected_real_kind(6)  ! 4 byte real
integer, parameter :: r8 = selected_real_kind(12) ! 8 byte real
real(r8), parameter :: pi = 3.1415926535897932384626434E0_r8
#endif

! Error functions
public :: shr_spfn_erf
public :: shr_spfn_erfc
public :: shr_spfn_erfc_scaled

interface shr_spfn_erf
   module procedure shr_spfn_erf_r4
   module procedure shr_spfn_erf_r8
end interface

interface shr_spfn_erfc
   module procedure shr_spfn_erfc_r4
   module procedure shr_spfn_erfc_r8
end interface

interface shr_spfn_erfc_scaled
   module procedure shr_spfn_erfc_scaled_r4
   module procedure shr_spfn_erfc_scaled_r8
end interface

! Gamma functions
! Note that we lack an implementation of log_gamma, but we do have an
! implementation of the upper incomplete gamma function, which is not in
! Fortran 2008.

! Note also that this gamma function is only for double precision. We
! haven't needed an r4 version yet.

public :: shr_spfn_gamma
public :: shr_spfn_igamma

interface shr_spfn_gamma
   module procedure shr_spfn_gamma_r8
end interface

! Mathematical constants
! sqrt(pi)
real(r8), parameter :: sqrtpi = 1.77245385090551602729_r8

! Define machine-specific constants needed in this module.
! These were used by the original gamma and calerf functions to guarantee
! safety against overflow, and precision, on many different machines.

! By defining the constants in this way, we assume that 1/xmin is
! representable (i.e. does not overflow the real type). This assumption was
! not in the original code, but is valid for IEEE single and double
! precision.

! Double precision
!---------------------------------------------------------------------
! Machine epsilon
real(r8), parameter :: epsr8 = epsilon(1._r8)
! "Huge" value is returned when actual value would be infinite.
real(r8), parameter :: xinfr8 = huge(1._r8)
! Smallest normal value.
real(r8), parameter :: xminr8 = tiny(1._r8)
! Largest number that, when added to 1., yields 1.
real(r8), parameter :: xsmallr8 = epsr8/2._r8
! Largest argument for which erfcx > 0.
real(r8), parameter :: xmaxr8 = 1._r8/(sqrtpi*xminr8)

! Single precision
!---------------------------------------------------------------------
! Machine epsilon
real(r4), parameter :: epsr4 = epsilon(1._r4)
! "Huge" value is returned when actual value would be infinite.
real(r4), parameter :: xinfr4 = huge(1._r4)
! Smallest normal value.
real(r4), parameter :: xminr4 = tiny(1._r4)
! Largest number that, when added to 1., yields 1.
real(r4), parameter :: xsmallr4 = epsr4/2._r4
! Largest argument for which erfcx > 0.
real(r4), parameter :: xmaxr4 = 1._r4/(real(sqrtpi,r4)*xminr4)


! For gamma/igamma
! Approximate value of largest acceptable argument to gamma,
! for IEEE double-precision.
real(r8), parameter :: xbig_gamma = 171.624_r8

contains

! Wrapper functions for erf
function shr_spfn_erf_r4(x) result(res)
  real(r4), intent(in) :: x
  real(r4) :: res

#ifdef HAVE_ERF_EXTERNALS
  ! If erf is provided as an external, provide
  ! explicit interface here.
  interface
     function erf(x)
       import :: r4
       real(r4) :: x, erf
     end function erf
  end interface
#endif

#ifdef HAVE_ERF_INTRINSICS

  ! Call intrinsic erf.
  intrinsic erf
  res = erf(x)
#else

#ifdef HAVE_ERF_EXTERNALS
  ! Call compiler-provided external erf.
  res = erf(x)
#else
  ! No compiler-provided erf, so call local version.
  call calerf_r4(x, res, 0)
#endif

#endif

end function shr_spfn_erf_r4

function shr_spfn_erf_r8(x) result(res)
  real(r8), intent(in) :: x
  real(r8) :: res

#ifdef HAVE_ERF_EXTERNALS
  ! If erf is provided as an external, provide
  ! explicit interface here.
  interface
     function derf(x)
       import :: r8
       real(r8) :: x, derf
     end function derf
  end interface
#endif

#ifdef HAVE_ERF_INTRINSICS
  ! Call intrinsic erf.
  intrinsic erf
  res = erf(x)
#else

#ifdef HAVE_ERF_EXTERNALS
  ! Call compiler-provided external erf.
  res = derf(x)
#else
  ! No compiler-provided erf, so call local version.
  call calerf_r8(x, res, 0)
#endif

#endif

end function shr_spfn_erf_r8

! Wrapper functions for erfc
function shr_spfn_erfc_r4(x) result(res)
  real(r4), intent(in) :: x
  real(r4) :: res

#ifdef HAVE_ERF_EXTERNALS
  ! If erfc is provided as an external, provide
  ! explicit interface here.
  interface
     function erfc(x)
       import :: r4
       real(r4) :: x, erfc
     end function erfc
  end interface
#endif

#ifdef HAVE_ERF_INTRINSICS
  ! Call intrinsic erfc.
  intrinsic erfc
  res = erfc(x)
#else

#ifdef HAVE_ERF_EXTERNALS
  ! Call compiler-provided external erfc.
  res = erfc(x)
#else
  ! No compiler-provided erfc, so call local version.
  call calerf_r4(x, res, 1)
#endif

#endif

end function shr_spfn_erfc_r4

function shr_spfn_erfc_r8(x) result(res)
  real(r8), intent(in) :: x
  real(r8) :: res

#ifdef HAVE_ERF_EXTERNALS
  ! If erfc is provided as an external, provide
  ! explicit interface here.
  interface
     function derfc(x)
       import :: r8
       real(r8) :: x, derfc
     end function derfc
  end interface
#endif

#ifdef HAVE_ERF_INTRINSICS
  ! Call intrinsic erfc.
  intrinsic erfc
  res = erfc(x)
#else

#ifdef HAVE_ERF_EXTERNALS
  ! Call compiler-provided external erfc.
  res = derfc(x)
#else
  ! No compiler-provided erfc, so call local version.
  call calerf_r8(x, res, 1)
#endif

#endif

end function shr_spfn_erfc_r8

! Wrapper functions for erfc_scaled
function shr_spfn_erfc_scaled_r4(x) result(res)
  real(r4), intent(in) :: x
  real(r4) :: res

#if defined HAVE_ERF_INTRINSICS
  ! Call intrinsic erfc_scaled.
  intrinsic erfc_scaled
  res = erfc_scaled(x)
#else
  ! No intrinsic.
  call calerf_r4(x, res, 2)
#endif

end function shr_spfn_erfc_scaled_r4

function shr_spfn_erfc_scaled_r8(x) result(res)
  real(r8), intent(in) :: x
  real(r8) :: res

#if defined HAVE_ERF_INTRINSICS
  ! Call intrinsic erfc_scaled.
  intrinsic erfc_scaled
  res = erfc_scaled(x)
#else
  ! No intrinsic.
  call calerf_r8(x, res, 2)
#endif

end function shr_spfn_erfc_scaled_r8

elemental function shr_spfn_gamma_r8(x) result(res)
  real(r8), intent(in) :: x
  real(r8) :: res

#if defined HAVE_GAMMA_INTRINSICS
  ! Call intrinsic gamma.
  intrinsic gamma
  res = gamma(x)
#else
  ! No intrinsic
  res = shr_spfn_gamma_nonintrinsic_r8(x)
#endif

end function shr_spfn_gamma_r8

!------------------------------------------------------------------
!
! 6 December 2006 -- B. Eaton
! The following comments are from the original version of CALERF.
! The only changes in implementing this module are that the function
! names previously used for the single precision versions have been
! adopted for the new generic interfaces.  To support these interfaces
! there is now both a single precision version (calerf_r4) and a
! double precision version (calerf_r8) of CALERF below.  These versions
! are hardcoded to use IEEE arithmetic.
!
!------------------------------------------------------------------
!
! This packet evaluates  erf(x),  erfc(x),  and  exp(x*x)*erfc(x)
!   for a real argument  x.  It contains three FUNCTION type
!   subprograms: ERF, ERFC, and ERFCX (or ERF_R8, ERFC_R8, and ERFCX_R8),
!   and one SUBROUTINE type subprogram, CALERF.  The calling
!   statements for the primary entries are:
!
!                   Y=ERF(X)     (or   Y=ERF_R8(X)),
!
!                   Y=ERFC(X)    (or   Y=ERFC_R8(X)),
!   and
!                   Y=ERFCX(X)   (or   Y=ERFCX_R8(X)).
!
!   The routine  CALERF  is intended for internal packet use only,
!   all computations within the packet being concentrated in this
!   routine.  The function subprograms invoke  CALERF  with the
!   statement
!
!          CALL CALERF(ARG,RESULT,JINT)
!
!   where the parameter usage is as follows
!
!      Function                     Parameters for CALERF
!       call              ARG                  Result          JINT
!
!     ERF(ARG)      ANY REAL ARGUMENT         ERF(ARG)          0
!     ERFC(ARG)     ABS(ARG) .LT. XBIG        ERFC(ARG)         1
!     ERFCX(ARG)    XNEG .LT. ARG .LT. XMAX   ERFCX(ARG)        2
!
!   The main computation evaluates near-minimax approximations
!   from "Rational Chebyshev approximations for the error function"
!   by W. J. Cody, Math. Comp., 1969, PP. 631-638.  This
!   transportable program uses rational functions that theoretically
!   approximate  erf(x)  and  erfc(x)  to at least 18 significant
!   decimal digits.  The accuracy achieved depends on the arithmetic
!   system, the compiler, the intrinsic functions, and proper
!   selection of the machine-dependent constants.
!
!*******************************************************************
!*******************************************************************
!
! Explanation of machine-dependent constants
!
!   XMIN   = the smallest positive floating-point number.
!   XINF   = the largest positive finite floating-point number.
!   XNEG   = the largest negative argument acceptable to ERFCX;
!            the negative of the solution to the equation
!            2*exp(x*x) = XINF.
!   XSMALL = argument below which erf(x) may be represented by
!            2*x/sqrt(pi)  and above which  x*x  will not underflow.
!            A conservative value is the largest machine number X
!            such that   1.0 + X = 1.0   to machine precision.
!   XBIG   = largest argument acceptable to ERFC;  solution to
!            the equation:  W(x) * (1-0.5/x**2) = XMIN,  where
!            W(x) = exp(-x*x)/[x*sqrt(pi)].
!   XHUGE  = argument above which  1.0 - 1/(2*x*x) = 1.0  to
!            machine precision.  A conservative value is
!            1/[2*sqrt(XSMALL)]
!   XMAX   = largest acceptable argument to ERFCX; the minimum
!            of XINF and 1/[sqrt(pi)*XMIN].
!
!   Approximate values for some important machines are:
!
!                          XMIN       XINF        XNEG     XSMALL
!
!  CDC 7600      (S.P.)  3.13E-294   1.26E+322   -27.220  7.11E-15
!  CRAY-1        (S.P.)  4.58E-2467  5.45E+2465  -75.345  7.11E-15
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)  1.18E-38    3.40E+38     -9.382  5.96E-8
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)  2.23D-308   1.79D+308   -26.628  1.11D-16
!  IBM 195       (D.P.)  5.40D-79    7.23E+75    -13.190  1.39D-17
!  UNIVAC 1108   (D.P.)  2.78D-309   8.98D+307   -26.615  1.73D-18
!  VAX D-Format  (D.P.)  2.94D-39    1.70D+38     -9.345  1.39D-17
!  VAX G-Format  (D.P.)  5.56D-309   8.98D+307   -26.615  1.11D-16
!
!
!                          XBIG       XHUGE       XMAX
!
!  CDC 7600      (S.P.)  25.922      8.39E+6     1.80X+293
!  CRAY-1        (S.P.)  75.326      8.39E+6     5.45E+2465
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)   9.194      2.90E+3     4.79E+37
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)  26.543      6.71D+7     2.53D+307
!  IBM 195       (D.P.)  13.306      1.90D+8     7.23E+75
!  UNIVAC 1108   (D.P.)  26.582      5.37D+8     8.98D+307
!  VAX D-Format  (D.P.)   9.269      1.90D+8     1.70D+38
!  VAX G-Format  (D.P.)  26.569      6.71D+7     8.98D+307
!
!*******************************************************************
!*******************************************************************
!
! Error returns
!
!  The program returns  ERFC = 0      for  ARG .GE. XBIG;
!
!                       ERFCX = XINF  for  ARG .LT. XNEG;
!      and
!                       ERFCX = 0     for  ARG .GE. XMAX.
!
!
! Intrinsic functions required are:
!
!     ABS, AINT, EXP
!
!
!  Author: W. J. Cody
!          Mathematics and Computer Science Division
!          Argonne National Laboratory
!          Argonne, IL 60439
!
!  Latest modification: March 19, 1990
!
!------------------------------------------------------------------

SUBROUTINE CALERF_r8(ARG, RESULT, JINT)

   !------------------------------------------------------------------
   !  This version uses 8-byte reals
   !------------------------------------------------------------------
   integer, parameter :: rk = r8

   ! arguments
   real(rk), intent(in)  :: arg
   integer,  intent(in)  :: jint
   real(rk), intent(out) :: result

   ! local variables
   INTEGER :: I

   real(rk) :: X, Y, YSQ, XNUM, XDEN, DEL

   !------------------------------------------------------------------
   !  Mathematical constants
   !------------------------------------------------------------------
   real(rk), parameter :: ZERO   = 0.0E0_rk
   real(rk), parameter :: FOUR   = 4.0E0_rk
   real(rk), parameter :: ONE    = 1.0E0_rk
   real(rk), parameter :: HALF   = 0.5E0_rk
   real(rk), parameter :: TWO    = 2.0E0_rk
   ! 1/sqrt(pi)
   real(rk), parameter :: SQRPI  = 5.6418958354775628695E-1_rk
   real(rk), parameter :: THRESH = 0.46875E0_rk
   real(rk), parameter :: SIXTEN = 16.0E0_rk

   !------------------------------------------------------------------
   !  Machine-dependent constants: IEEE double precision values
   !------------------------------------------------------------------
   real(rk), parameter :: XNEG   = -26.628E0_r8
   real(rk), parameter :: XBIG   =  26.543E0_r8
   real(rk), parameter :: XHUGE  =   6.71E7_r8

   !------------------------------------------------------------------
   !  Coefficients for approximation to  erf  in first interval
   !------------------------------------------------------------------
   real(rk), parameter :: A(5) = (/ 3.16112374387056560E00_rk, 1.13864154151050156E02_rk, &
                                    3.77485237685302021E02_rk, 3.20937758913846947E03_rk, &
                                    1.85777706184603153E-1_rk /)
   real(rk), parameter :: B(4) = (/ 2.36012909523441209E01_rk, 2.44024637934444173E02_rk, &
                                    1.28261652607737228E03_rk, 2.84423683343917062E03_rk /)

   !------------------------------------------------------------------
   !  Coefficients for approximation to  erfc  in second interval
   !------------------------------------------------------------------
   real(rk), parameter :: C(9) = (/ 5.64188496988670089E-1_rk, 8.88314979438837594E00_rk, &
                                    6.61191906371416295E01_rk, 2.98635138197400131E02_rk, &
                                    8.81952221241769090E02_rk, 1.71204761263407058E03_rk, &
                                    2.05107837782607147E03_rk, 1.23033935479799725E03_rk, &
                                    2.15311535474403846E-8_rk /)
   real(rk), parameter :: D(8) = (/ 1.57449261107098347E01_rk, 1.17693950891312499E02_rk, &
                                    5.37181101862009858E02_rk, 1.62138957456669019E03_rk, &
                                    3.29079923573345963E03_rk, 4.36261909014324716E03_rk, &
                                    3.43936767414372164E03_rk, 1.23033935480374942E03_rk /)

   !------------------------------------------------------------------
   !  Coefficients for approximation to  erfc  in third interval
   !------------------------------------------------------------------
   real(rk), parameter :: P(6) = (/ 3.05326634961232344E-1_rk, 3.60344899949804439E-1_rk, &
                                    1.25781726111229246E-1_rk, 1.60837851487422766E-2_rk, &
                                    6.58749161529837803E-4_rk, 1.63153871373020978E-2_rk /)
   real(rk), parameter :: Q(5) = (/ 2.56852019228982242E00_rk, 1.87295284992346047E00_rk, &
                                    5.27905102951428412E-1_rk, 6.05183413124413191E-2_rk, &
                                    2.33520497626869185E-3_rk /)

   !------------------------------------------------------------------
   X = ARG
   Y = ABS(X)
   IF (Y .LE. THRESH) THEN
      !------------------------------------------------------------------
      !  Evaluate  erf  for  |X| <= 0.46875
      !------------------------------------------------------------------
      YSQ = ZERO
      IF (Y .GT. XSMALLR8) YSQ = Y * Y
      XNUM = A(5)*YSQ
      XDEN = YSQ
      DO I = 1, 3
         XNUM = (XNUM + A(I)) * YSQ
         XDEN = (XDEN + B(I)) * YSQ
      end do
      RESULT = X * (XNUM + A(4)) / (XDEN + B(4))
      IF (JINT .NE. 0) RESULT = ONE - RESULT
      IF (JINT .EQ. 2) RESULT = EXP(YSQ) * RESULT
      GO TO 80
   ELSE IF (Y .LE. FOUR) THEN
      !------------------------------------------------------------------
      !  Evaluate  erfc  for 0.46875 <= |X| <= 4.0
      !------------------------------------------------------------------
      XNUM = C(9)*Y
      XDEN = Y
      DO I = 1, 7
         XNUM = (XNUM + C(I)) * Y
         XDEN = (XDEN + D(I)) * Y
      end do
      RESULT = (XNUM + C(8)) / (XDEN + D(8))
      IF (JINT .NE. 2) THEN
         YSQ = AINT(Y*SIXTEN)/SIXTEN
         DEL = (Y-YSQ)*(Y+YSQ)
         RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT
      END IF
   ELSE
      !------------------------------------------------------------------
      !  Evaluate  erfc  for |X| > 4.0
      !------------------------------------------------------------------
      RESULT = ZERO
      IF (Y .GE. XBIG) THEN
         IF ((JINT .NE. 2) .OR. (Y .GE. XMAXR8)) GO TO 30
         IF (Y .GE. XHUGE) THEN
            RESULT = SQRPI / Y
            GO TO 30
         END IF
      END IF
      YSQ = ONE / (Y * Y)
      XNUM = P(6)*YSQ
      XDEN = YSQ
      DO I = 1, 4
         XNUM = (XNUM + P(I)) * YSQ
         XDEN = (XDEN + Q(I)) * YSQ
      end do
      RESULT = YSQ *(XNUM + P(5)) / (XDEN + Q(5))
      RESULT = (SQRPI -  RESULT) / Y
      IF (JINT .NE. 2) THEN
         YSQ = AINT(Y*SIXTEN)/SIXTEN
         DEL = (Y-YSQ)*(Y+YSQ)
         RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT
      END IF
   END IF
30 continue
   !------------------------------------------------------------------
   !  Fix up for negative argument, erf, etc.
   !------------------------------------------------------------------
   IF (JINT .EQ. 0) THEN
      RESULT = (HALF - RESULT) + HALF
      IF (X .LT. ZERO) RESULT = -RESULT
   ELSE IF (JINT .EQ. 1) THEN
      IF (X .LT. ZERO) RESULT = TWO - RESULT
   ELSE
      IF (X .LT. ZERO) THEN
         IF (X .LT. XNEG) THEN
            RESULT = XINFR8
         ELSE
            YSQ = AINT(X*SIXTEN)/SIXTEN
            DEL = (X-YSQ)*(X+YSQ)
            Y = EXP(YSQ*YSQ) * EXP(DEL)
            RESULT = (Y+Y) - RESULT
         END IF
      END IF
   END IF
80 continue
end SUBROUTINE CALERF_r8

!------------------------------------------------------------------------------------------

SUBROUTINE CALERF_r4(ARG, RESULT, JINT)

   !------------------------------------------------------------------
   !  This version uses 4-byte reals
   !------------------------------------------------------------------
   integer, parameter :: rk = r4

   ! arguments
   real(rk), intent(in)  :: arg
   integer,  intent(in)  :: jint
   real(rk), intent(out) :: result

   ! local variables
   INTEGER :: I

   real(rk) :: X, Y, YSQ, XNUM, XDEN, DEL

   !------------------------------------------------------------------
   !  Mathematical constants
   !------------------------------------------------------------------
   real(rk), parameter :: ZERO   = 0.0E0_rk
   real(rk), parameter :: FOUR   = 4.0E0_rk
   real(rk), parameter :: ONE    = 1.0E0_rk
   real(rk), parameter :: HALF   = 0.5E0_rk
   real(rk), parameter :: TWO    = 2.0E0_rk
   ! 1/sqrt(pi)
   real(rk), parameter :: SQRPI  = 5.6418958354775628695E-1_rk
   real(rk), parameter :: THRESH = 0.46875E0_rk
   real(rk), parameter :: SIXTEN = 16.0E0_rk

   !------------------------------------------------------------------
   !  Machine-dependent constants: IEEE single precision values
   !------------------------------------------------------------------
   real(rk), parameter :: XNEG   = -9.382E0_r4
   real(rk), parameter :: XBIG   =  9.194E0_r4
   real(rk), parameter :: XHUGE  =  2.90E3_r4

   !------------------------------------------------------------------
   !  Coefficients for approximation to  erf  in first interval
   !------------------------------------------------------------------
   real(rk), parameter :: A(5) = (/ 3.16112374387056560E00_rk, 1.13864154151050156E02_rk, &
                                    3.77485237685302021E02_rk, 3.20937758913846947E03_rk, &
                                    1.85777706184603153E-1_rk /)
   real(rk), parameter :: B(4) = (/ 2.36012909523441209E01_rk, 2.44024637934444173E02_rk, &
                                    1.28261652607737228E03_rk, 2.84423683343917062E03_rk /)

   !------------------------------------------------------------------
   !  Coefficients for approximation to  erfc  in second interval
   !------------------------------------------------------------------
   real(rk), parameter :: C(9) = (/ 5.64188496988670089E-1_rk, 8.88314979438837594E00_rk, &
                                    6.61191906371416295E01_rk, 2.98635138197400131E02_rk, &
                                    8.81952221241769090E02_rk, 1.71204761263407058E03_rk, &
                                    2.05107837782607147E03_rk, 1.23033935479799725E03_rk, &
                                    2.15311535474403846E-8_rk /)
   real(rk), parameter :: D(8) = (/ 1.57449261107098347E01_rk, 1.17693950891312499E02_rk, &
                                    5.37181101862009858E02_rk, 1.62138957456669019E03_rk, &
                                    3.29079923573345963E03_rk, 4.36261909014324716E03_rk, &
                                    3.43936767414372164E03_rk, 1.23033935480374942E03_rk /)

   !------------------------------------------------------------------
   !  Coefficients for approximation to  erfc  in third interval
   !------------------------------------------------------------------
   real(rk), parameter :: P(6) = (/ 3.05326634961232344E-1_rk, 3.60344899949804439E-1_rk, &
                                    1.25781726111229246E-1_rk, 1.60837851487422766E-2_rk, &
                                    6.58749161529837803E-4_rk, 1.63153871373020978E-2_rk /)
   real(rk), parameter :: Q(5) = (/ 2.56852019228982242E00_rk, 1.87295284992346047E00_rk, &
                                    5.27905102951428412E-1_rk, 6.05183413124413191E-2_rk, &
                                    2.33520497626869185E-3_rk /)

   !------------------------------------------------------------------
   X = ARG
   Y = ABS(X)
   IF (Y .LE. THRESH) THEN
      !------------------------------------------------------------------
      !  Evaluate  erf  for  |X| <= 0.46875
      !------------------------------------------------------------------
      YSQ = ZERO
      IF (Y .GT. XSMALLR4) YSQ = Y * Y
      XNUM = A(5)*YSQ
      XDEN = YSQ
      DO I = 1, 3
         XNUM = (XNUM + A(I)) * YSQ
         XDEN = (XDEN + B(I)) * YSQ
      end do
      RESULT = X * (XNUM + A(4)) / (XDEN + B(4))
      IF (JINT .NE. 0) RESULT = ONE - RESULT
      IF (JINT .EQ. 2) RESULT = EXP(YSQ) * RESULT
      GO TO 80
   ELSE IF (Y .LE. FOUR) THEN
      !------------------------------------------------------------------
      !  Evaluate  erfc  for 0.46875 <= |X| <= 4.0
      !------------------------------------------------------------------
      XNUM = C(9)*Y
      XDEN = Y
      DO I = 1, 7
         XNUM = (XNUM + C(I)) * Y
         XDEN = (XDEN + D(I)) * Y
      end do
      RESULT = (XNUM + C(8)) / (XDEN + D(8))
      IF (JINT .NE. 2) THEN
         YSQ = AINT(Y*SIXTEN)/SIXTEN
         DEL = (Y-YSQ)*(Y+YSQ)
         RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT
      END IF
   ELSE
      !------------------------------------------------------------------
      !  Evaluate  erfc  for |X| > 4.0
      !------------------------------------------------------------------
      RESULT = ZERO
      IF (Y .GE. XBIG) THEN
         IF ((JINT .NE. 2) .OR. (Y .GE. XMAXR4)) GO TO 30
         IF (Y .GE. XHUGE) THEN
            RESULT = SQRPI / Y
            GO TO 30
         END IF
      END IF
      YSQ = ONE / (Y * Y)
      XNUM = P(6)*YSQ
      XDEN = YSQ
      DO I = 1, 4
         XNUM = (XNUM + P(I)) * YSQ
         XDEN = (XDEN + Q(I)) * YSQ
      end do
      RESULT = YSQ *(XNUM + P(5)) / (XDEN + Q(5))
      RESULT = (SQRPI -  RESULT) / Y
      IF (JINT .NE. 2) THEN
         YSQ = AINT(Y*SIXTEN)/SIXTEN
         DEL = (Y-YSQ)*(Y+YSQ)
         RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT
      END IF
   END IF
30 continue
   !------------------------------------------------------------------
   !  Fix up for negative argument, erf, etc.
   !------------------------------------------------------------------
   IF (JINT .EQ. 0) THEN
      RESULT = (HALF - RESULT) + HALF
      IF (X .LT. ZERO) RESULT = -RESULT
   ELSE IF (JINT .EQ. 1) THEN
      IF (X .LT. ZERO) RESULT = TWO - RESULT
   ELSE
      IF (X .LT. ZERO) THEN
         IF (X .LT. XNEG) THEN
            RESULT = XINFR4
         ELSE
            YSQ = AINT(X*SIXTEN)/SIXTEN
            DEL = (X-YSQ)*(X+YSQ)
            Y = EXP(YSQ*YSQ) * EXP(DEL)
            RESULT = (Y+Y) - RESULT
         END IF
      END IF
   END IF
80 continue
end SUBROUTINE CALERF_r4

!------------------------------------------------------------------------------------------

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

pure function shr_spfn_gamma_nonintrinsic_r8(X) result(gamma)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! 7 Feb 2013 -- S. Santos
! The following comments are from the original version. Changes have
! been made to update syntax and allow inclusion into this module.
!
!----------------------------------------------------------------------
!
! THIS ROUTINE CALCULATES THE GAMMA FUNCTION FOR A REAL ARGUMENT X.
!   COMPUTATION IS BASED ON AN ALGORITHM OUTLINED IN REFERENCE 1.
!   THE PROGRAM USES RATIONAL FUNCTIONS THAT APPROXIMATE THE GAMMA
!   FUNCTION TO AT LEAST 20 SIGNIFICANT DECIMAL DIGITS.  COEFFICIENTS
!   FOR THE APPROXIMATION OVER THE INTERVAL (1,2) ARE UNPUBLISHED.
!   THOSE FOR THE APPROXIMATION FOR X .GE. 12 ARE FROM REFERENCE 2.
!   THE ACCURACY ACHIEVED DEPENDS ON THE ARITHMETIC SYSTEM, THE
!   COMPILER, THE INTRINSIC FUNCTIONS, AND PROPER SELECTION OF THE
!   MACHINE-DEPENDENT CONSTANTS.
!
!
!*******************************************************************
!*******************************************************************
!
! EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
!
! BETA   - RADIX FOR THE FLOATING-POINT REPRESENTATION
! MAXEXP - THE SMALLEST POSITIVE POWER OF BETA THAT OVERFLOWS
! XBIG   - THE LARGEST ARGUMENT FOR WHICH GAMMA(X) IS REPRESENTABLE
!          IN THE MACHINE, I.E., THE SOLUTION TO THE EQUATION
!                  GAMMA(XBIG) = BETA**MAXEXP
! XINF   - THE LARGEST MACHINE REPRESENTABLE FLOATING-POINT NUMBER;
!          APPROXIMATELY BETA**MAXEXP
! EPS    - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1.0+EPS .GT. 1.0
! XMININ - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1/XMININ IS MACHINE REPRESENTABLE
!
!     APPROXIMATE VALUES FOR SOME IMPORTANT MACHINES ARE:
!
!                            BETA       MAXEXP        XBIG
!
! CRAY-1         (S.P.)        2         8191        966.961
! CYBER 180/855
!   UNDER NOS    (S.P.)        2         1070        177.803
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)        2          128        35.040
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)        2         1024        171.624
! IBM 3033       (D.P.)       16           63        57.574
! VAX D-FORMAT   (D.P.)        2          127        34.844
! VAX G-FORMAT   (D.P.)        2         1023        171.489
!
!                            XINF         EPS        XMININ
!
! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
! CYBER 180/855
!   UNDER NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
! VAX D-FORMAT   (D.P.)   1.70D+38     1.39D-17    5.88D-39
! VAX G-FORMAT   (D.P.)   8.98D+307    1.11D-16    1.12D-308
!
!*******************************************************************
!*******************************************************************
!
! ERROR RETURNS
!
!  THE PROGRAM RETURNS THE VALUE XINF FOR SINGULARITIES OR
!     WHEN OVERFLOW WOULD OCCUR.  THE COMPUTATION IS BELIEVED
!     TO BE FREE OF UNDERFLOW AND OVERFLOW.
!
!
!  INTRINSIC FUNCTIONS REQUIRED ARE:
!
!     INT, DBLE, EXP, LOG, REAL, SIN
!
!
! REFERENCES:  AN OVERVIEW OF SOFTWARE DEVELOPMENT FOR SPECIAL
!              FUNCTIONS   W. J. CODY, LECTURE NOTES IN MATHEMATICS,
!              506, NUMERICAL ANALYSIS DUNDEE, 1975, G. A. WATSON
!              (ED.), SPRINGER VERLAG, BERLIN, 1976.
!
!              COMPUTER APPROXIMATIONS, HART, ET. AL., WILEY AND
!              SONS, NEW YORK, 1968.
!
!  LATEST MODIFICATION: OCTOBER 12, 1989
!
!  AUTHORS: W. J. CODY AND L. STOLTZ
!           APPLIED MATHEMATICS DIVISION
!           ARGONNE NATIONAL LABORATORY
!           ARGONNE, IL 60439
!
!----------------------------------------------------------------------

  real(r8), intent(in) :: x
  real(r8) :: gamma
  real(r8) :: fact, res, sum, xden, xnum, y, y1, ysq, z

  integer :: i, n
  logical :: negative_odd

  ! log(2*pi)/2
  real(r8), parameter :: logsqrt2pi = 0.9189385332046727417803297E0_r8

!----------------------------------------------------------------------
!  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
!     APPROXIMATION OVER (1,2).
!----------------------------------------------------------------------
  real(r8), parameter :: P(8) = &
       (/-1.71618513886549492533811E+0_r8, 2.47656508055759199108314E+1_r8, &
         -3.79804256470945635097577E+2_r8, 6.29331155312818442661052E+2_r8, &
          8.66966202790413211295064E+2_r8,-3.14512729688483675254357E+4_r8, &
         -3.61444134186911729807069E+4_r8, 6.64561438202405440627855E+4_r8 /)
  real(r8), parameter :: Q(8) = &
       (/-3.08402300119738975254353E+1_r8, 3.15350626979604161529144E+2_r8, &
         -1.01515636749021914166146E+3_r8,-3.10777167157231109440444E+3_r8, &
          2.25381184209801510330112E+4_r8, 4.75584627752788110767815E+3_r8, &
         -1.34659959864969306392456E+5_r8,-1.15132259675553483497211E+5_r8 /)
!----------------------------------------------------------------------
!  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
!----------------------------------------------------------------------
  real(r8), parameter :: C(7) = &
       (/-1.910444077728E-03_r8,          8.4171387781295E-04_r8, &
         -5.952379913043012E-04_r8,       7.93650793500350248E-04_r8, &
         -2.777777777777681622553E-03_r8, 8.333333333333333331554247E-02_r8, &
          5.7083835261E-03_r8 /)

  negative_odd = .false.
  fact = 1._r8
  n = 0
  y = x
  if (y <= 0._r8) then
!----------------------------------------------------------------------
!  ARGUMENT IS NEGATIVE
!----------------------------------------------------------------------
     y = -x
     y1 = aint(y)
     res = y - y1
     if (res /= 0._r8) then
        negative_odd = (y1 /= aint(y1*0.5_r8)*2._r8)
        fact = -pi/sin(pi*res)
        y = y + 1._r8
     else
        gamma = xinfr8
        return
     end if
  end if
!----------------------------------------------------------------------
!  ARGUMENT IS POSITIVE
!----------------------------------------------------------------------
  if (y < epsr8) then
!----------------------------------------------------------------------
!  ARGUMENT .LT. EPS
!----------------------------------------------------------------------
     if (y >= xminr8) then
        res = 1._r8/y
     else
        gamma = xinfr8
        return
     end if
  elseif (y < 12._r8) then
     y1 = y
     if (y < 1._r8) then
!----------------------------------------------------------------------
!  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
        z = y
        y = y + 1._r8
     else
!----------------------------------------------------------------------
!  1.0 .LT. ARGUMENT .LT. 12.0, REDUCE ARGUMENT IF NECESSARY
!----------------------------------------------------------------------
        n = int(y) - 1
        y = y - real(n, r8)
        z = y - 1._r8
     end if
!----------------------------------------------------------------------
!  EVALUATE APPROXIMATION FOR 1.0 .LT. ARGUMENT .LT. 2.0
!----------------------------------------------------------------------
     xnum = 0._r8
     xden = 1._r8
     do i=1,8
        xnum = (xnum+P(i))*z
        xden = xden*z + Q(i)
     end do
     res = xnum/xden + 1._r8
     if (y1 < y) then
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
        res = res/y1
     elseif (y1 > y) then
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  2.0 .LT. ARGUMENT .LT. 12.0
!----------------------------------------------------------------------
        do i = 1,n
           res = res*y
           y = y + 1._r8
        end do
     end if
  else
!----------------------------------------------------------------------
!  EVALUATE FOR ARGUMENT .GE. 12.0,
!----------------------------------------------------------------------
     if (y <= xbig_gamma) then
        ysq = y*y
        sum = C(7)
        do i=1,6
           sum = sum/ysq + C(i)
        end do
        sum = sum/y - y + logsqrt2pi
        sum = sum + (y-0.5_r8)*log(y)
        res = exp(sum)
     else
        gamma = xinfr8
        return
     end if
  end if
!----------------------------------------------------------------------
!  FINAL ADJUSTMENTS AND RETURN
!----------------------------------------------------------------------
  if (negative_odd)  res = -res
  if (fact /= 1._r8) res = fact/res
  gamma = res
! ---------- LAST LINE OF GAMMA ----------
end function shr_spfn_gamma_nonintrinsic_r8

!! Incomplete Gamma function
!!
!! @author  Tianyi Fan
!! @version August-2010
real(r8) elemental function shr_spfn_igamma(a, x)
  ! Upper incomplete gamma function.
  ! Modified for inclusion in this module and made
  ! pure elemental, September 2012

  real(r8), intent(in) ::      a
  real(r8), intent(in) ::      x

  ! local variable
  real(r8) :: xam, gin, s, r, t0
  integer  :: k


  if (x == 0.0_r8) then
     shr_spfn_igamma = shr_spfn_gamma(a)
     return
  end if

  xam = -x + a * log(x)

  if ((xam > 700.0_r8) .or. (a > xbig_gamma)) then
     ! Out of bounds
     ! Return "huge" value.
     shr_spfn_igamma = xinfr8
     return

  else if (x <= (1.0_r8 + a)) then
     s = 1.0_r8 / a
     r = s

     do  k = 1,60
        r = r * x / (a+k)
        s = s + r

        if (abs(r/s) < 1.0e-15_r8) exit
     end do

     gin = exp(xam) * s
     shr_spfn_igamma = shr_spfn_gamma(a) - gin

  else
     t0 = 0.0_r8

     do k = 60,1,-1
        t0 = (k - a) / (1.0_r8 + k / (x + t0))
     end do

     shr_spfn_igamma = exp(xam) / (x + t0)
  endif

end function shr_spfn_igamma


end module shr_spfn_mod
