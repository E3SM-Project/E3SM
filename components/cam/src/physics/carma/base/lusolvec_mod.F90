!! 
!! this module: lusolvec_mod   Numerical solution of a set of linear
!!                             Equations / a matrix equation A * x = b
!!                             using LU decomposition, matrix A and
!!                             vectors b and x being double complex,
!!                             and inversion of A.
!! ******************************************************************
!! Usage:
!! ====== 
!! given a complex matrix A, a right hand side vector b and a 
!! matrix equation      A * x = b      to solve for vector x.
!! 
!!                                                                     
!! First, call LUDCMPC(A,N,NP,INDX,D). The original Matrix A is lost   
!! and substituted by its LU decomposition.                            
!!                                                                     
!! Second, call LUBKSBC(A,N,NP,INDX,B). The original right-hand-side   
!! vector ib in B is lost and replaced/returned as the solution        
!! vector x  ( x(i) = B(i) ).                                          
!! Use same kind of call to solve for successive right-hand-sides.     
!!                                                                     
!! For Inversion of matrix A, call LUBKSBC() subsequently for each     
!! column vector:                                                      
!!   1) Initialize matrix AINV(i,j) to be equal to the                 
!!      identity matrix (AINV(i,j)=1 for i=j; =0 otherwise)            
!!   2) DO jj=1,n                                                      
!!         CALL LUBKSBC(A,N,NP,INDX,AINV(1,jj))                        
!!      END DO                                                         
!! (see textbook for further details).                                 
!! ****************************************************************** 

module lusolvec_mod

  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod

  implicit none
  private

  ! public subroutines
  public :: LUDCMPC
  public :: LUBKSBC

  contains

  !!
  !! SUBROUTINE LUDCMPC(A,N,NP,INDX,D)                                   
  !!
  !! Given a general complex matrix A, this routine replaces it by its   
  !! LU decomposition of a rowwise permutation of itself.                
  !! This routine is used in combination with LUBKSBC(), a complex       
  !! extension of the routine LUBKSB() (DOUBLE COMPLEX).                 
  !! For further details, refer to textbook (see below).                 
  !!
  !! Source: Own adaption/extension to complex matrix of the             
  !!         Subroutine LUDCMP() taken from                              
  !!         Press et al, "Numerical Recipes in Fortran"                 
  !!         The adaption follows the statements given in section 2.3    
  !!         of the textbook "N.R. in C", following Eq.(2.3.16):         
  !!         - definition of variables, vector and matrix elements       
  !!           as complex variables (use of complex arithmetic does      
  !!           not necessitate any adaption in fortran).                 
  !!         - complex modulus instead of absolute values in the         
  !!           construction of the vector vv and in the search for the   
  !!           largest pivot elements.                                   
  !! ******************************************************************  
  !! Version: 28.08.2000                                                 
  !! ******************************************************************
  SUBROUTINE LUDCMPC(A,N,NP,INDX,D)

    INTEGER :: NP
    COMPLEX(kind=f) :: A(NP,NP)
    INTEGER :: N
    INTEGER :: INDX(N)
    REAL(kind=f) :: D

    INTEGER, PARAMETER :: NMAX=100
    REAL(kind=f), PARAMETER :: TINY=1.0e-20_f
    REAL(kind=f) :: VV(NMAX)
    REAL(kind=f) :: DUM,AAMAX
    COMPLEX(kind=f) :: SUM,DUMC,ZEROC,TINYC
    INTEGER I,J,K,IMAX

    D=1._f
    TINYC=cmplx(TINY,0.0_f,kind=f)
    ZEROC=cmplx(0.0_f,0.0_f,kind=f)
    DO I=1,N
      AAMAX=0._f
      DO J=1,N
        IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
      END DO
!      IF (AAMAX.EQ.0._f) PAUSE 'Singular matrix.'
      IF (AAMAX.EQ.0._f) STOP 'Singular matrix.'
      VV(I)=1./AAMAX
    END DO
    DO J=1,N
      IF (J.GT.1) THEN
        DO I=1,J-1
          SUM=A(I,J)
          IF (I.GT.1)THEN
            DO K=1,I-1
              SUM=SUM-A(I,K)*A(K,J)
            END DO
            A(I,J)=SUM
          ENDIF
        END DO
      ENDIF
      AAMAX=0._f
      DO I=J,N
        SUM=A(I,J)
        IF (J.GT.1)THEN
          DO K=1,J-1
            SUM=SUM-A(I,K)*A(K,J)
          END DO
          A(I,J)=SUM
        ENDIF
        DUM=VV(I)*ABS(SUM)
        IF (DUM.GE.AAMAX) THEN
          IMAX=I
          AAMAX=DUM
        ENDIF
      END DO
      IF (J.NE.IMAX)THEN
        DO K=1,N
          DUMC=A(IMAX,K)
          A(IMAX,K)=A(J,K)
          A(J,K)=DUMC
        END DO
        D=-D
        VV(IMAX)=VV(J)
      ENDIF
      INDX(J)=IMAX
      IF(J.NE.N)THEN
        IF(A(J,J).EQ.ZEROC)A(J,J)=TINYC
        DUMC=1./A(J,J)
        DO I=J+1,N
          A(I,J)=A(I,J)*DUMC
        END DO
      ENDIF
    END DO
    IF (A(N,N).EQ.ZEROC) A(N,N)=TINYC
    RETURN
  END SUBROUTINE LUDCMPC

  !!
  !! SUBROUTINE LUBKSBC(A,N,NP,INDX,B)
  !!
  !! Solution of the set of linear equations A' * x = b where 
  !! A is input not as the original matrix, but as a LU decomposition 
  !! of some original matrix A' as determined by the subroutine
  !! LUDCMPC() (matrix and vectors being of type DOUBLE COMPLEX). 
  !! INDX() is input as the permutation vactor returned by LUDCMPC().
  !! B() is input as the right hand side vector b of the Eqn. to solve
  !! and returns with the solution vector x. 
  !! A, N and INDX are not modified by this routine and can be left in
  !! place for successive calls with different right-hand-sides b.
  !! For further details, refer to textbook (see below).
  !!
  !! Source: Own adaption/extension to complex matrix of the
  !!         Subroutine LUBKSB() taken from
  !!         Press et al, "Numerical Recipes in Fortran"
  !!         The adaption follows the statements given in section 2.3
  !!         of the textbook "N.R. in C", following Eq.(2.3.16).
  !! ******************************************************************
  !! Version: 28.08.2000
  !! ******************************************************************
  SUBROUTINE LUBKSBC(A,N,NP,INDX,B)

    INTEGER :: NP
    COMPLEX(kind=f) :: A(NP,NP)
    INTEGER :: N
    INTEGER :: INDX(N)
    COMPLEX(kind=f) :: B(N)

    INTEGER, PARAMETER :: NMAX=100
    REAL(kind=f), PARAMETER :: TINY=1.0e-20_f

    COMPLEX(kind=f) :: SUM,ZEROC
    INTEGER :: II,LL,I,J

    II=0
    ZEROC=cmplx(0.0_f,0.0_f,kind=f)
    DO I=1,N
      LL=INDX(I)
      SUM=B(LL)
      B(LL)=B(I)
      IF (II.NE.0)THEN
        DO J=II,I-1
          SUM=SUM-A(I,J)*B(J)
        END DO
      ELSE IF (SUM.NE.ZEROC) THEN
        II=I
      ENDIF
      B(I)=SUM
    END DO
    DO I=N,1,-1
      SUM=B(I)
      IF(I.LT.N)THEN
        DO J=I+1,N
          SUM=SUM-A(I,J)*B(J)
        END DO
      ENDIF
      B(I)=SUM/A(I,I)
    END DO
    RETURN
  END SUBROUTINE LUBKSBC


end module lusolvec_mod
