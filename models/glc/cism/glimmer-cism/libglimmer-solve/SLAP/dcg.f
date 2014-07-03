*DECK DCG
      SUBROUTINE DCG(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE, 
     $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, DZ,
     $     RWORK, IWORK )
C***BEGIN PROLOGUE  DCG
C***DATE WRITTEN   890404   (YYMMDD)
C***REVISION DATE  890404   (YYMMDD)
C***CATEGORY NO.  D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=DOUBLE PRECISION(DCG-D),
C             Symmetric Linear system, Sparse, Iterative Precondition
C***AUTHOR  Greenbaum, Anne, Courant Institute
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  Preconditioned Conjugate Gradient iterative Ax=b solver.
C            Routine to  solve a  symmetric positive definite linear
C            system    Ax = b    using the Preconditioned  Conjugate
C            Gradient method.
C***DESCRIPTION
C *Usage:
C     INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX
C     INTEGER  ITER, IERR, IUNIT, IWORK(USER DEFINABLE)
C     DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, R(N), Z(N)
C     DOUBLE PRECISION P(N), DZ(N), RWORK(USER DEFINABLE)
C     EXTERNAL MATVEC, MSOLVE
C
C     CALL DCG(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSLOVE, 
C    $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, DZ, 
C    $     RWORK, IWORK )
C         
C *Arguments:
C N      :IN       Integer.
C         Order of the Matrix.
C B      :IN       Double Precision B(N).
C         Right-hand side vector.
C X      :INOUT    Double Precision X(N).
C         On input X is your initial guess for solution vector.
C         On output X is the final approximate solution.
C NELT   :IN       Integer.
C         Number of Non-Zeros stored in A.
C IA     :IN       Integer IA(NELT).
C JA     :IN       Integer JA(NELT).
C A      :IN       Integer A(NELT).
C         These arrays contain the matrix data structure for A.
C         It could take any form.  See ``Description'', below
C         for more late breaking details...
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all nonzero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the upper
C         or lower triangle of the matrix is stored.
C MATVEC :EXT      External.
C         Name of a routine which performs the matrix vector multiply
C         Y = A*X given A and X.  The name of the MATVEC routine must 
C         be declared external in the calling program.  The calling 
C         sequence to MATVEC is:
C                  
C             CALL MATVEC( N, X, Y, NELT, IA, JA, A, ISYM )
C                     
C         Where N is the number of unknowns, Y is the product A*X
C         upon return X is an input vector, NELT is the number of 
C         non-zeros in the SLAP IA, JA, A storage for the matrix A.  
C         ISYM is a flag which, if non-zero, denotest that A is 
C         symmetric and only the lower or upper triangle is stored.
C MSOLVE :EXT      External.
C         Name of a routine which solves a linear system MZ = R for
C         Z given R with the preconditioning matrix M (M is supplied via
C         RWORK and IWORK arrays).  The name of the MSOLVE routine must 
C         be declared external in the calling program.  The calling 
C         sequence to MSLOVE is:
C                                                                
C             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
C                         
C         Where N is the number of unknowns, R is the right-hand side 
C         vector, and Z is the solution upon return.  RWORK is a double
C         precision 
C         array that can be used to pass necessary preconditioning 
C         information and/or workspace to MSOLVE.  IWORK is an integer 
C         work array for the same purpose as RWORK.
C ITOL   :IN       Integer.
C         Flag to indicate type of convergence criterion.
C         If ITOL=1, iteration stops when the 2-norm of the residual 
C         divided by the 2-norm of the right-hand side is less than TOL.
C         If ITOL=2, iteration stops when the 2-norm of M-inv times the 
C         residual divided by the 2-norm of M-inv times the right hand 
C         side is less than TOL, where M-inv is the inverse of the 
C         diagonal of A.
C         ITOL=11 is often useful for checking and comparing different 
C         routines.  For this case, the user must supply the "exact" 
C         solution or a very accurate approximation (one with an error 
C         much less than TOL) through a common block,
C                     COMMON /SOLBLK/ SOLN(1)
C         if ITOL=11, iteration stops when the 2-norm of the difference 
C         between the iterative approximation and the user-supplied
C         solution divided by the 2-norm of the user-supplied solution 
C         is less than TOL.  Note that this requires the user to set up
C         the "COMMON /SOLBLK/ SOLN(LENGTH)" in the calling routine. 
C         The routine with this declaration should be loaded before the
C         stop test so that the correct length is used by the loader.  
C         This procedure is not standard Fortran and may not work 
C         correctly on your system (although it has worked on every
C         system the authors have tried).  If ITOL is not 11 then this
C         common block is indeed standard Fortran.
C TOL    :IN       Double Precision.
C         Convergence criterion, as described above.
C ITMAX  :IN       Integer.
C         Maximum number of iterations.
C ITER   :OUT      Integer.
C         Number of iterations required to reach convergence, or 
C         ITMAX+1 if convergence criterion could not be achieved in 
C         ITMAX iterations.
C ERR    :OUT      Double Precision.
C         Error estimate of error in final approximate solution, as 
C         defined by ITOL.
C IERR   :OUT      Integer.
C         Return error flag.
C           IERR = 0 => All went well.
C           IERR = 1 => Insufficient storage allocated 
C                       for WORK or IWORK.
C           IERR = 2 => Method failed to converge in 
C                       ITMAX steps.
C           IERR = 3 => Error in user input.  Check input
C                       value of N, ITOL.
C           IERR = 4 => User error tolerance set too tight.
C                       Reset to 500.0*R1MACH(3).  Iteration proceeded.
C           IERR = 5 => Preconditioning matrix, M,  is not 
C                       Positive Definite.  $(r,z) < 0.0$.
C           IERR = 6 => Matrix A is not Positive Definite.  
C                       $(p,Ap) < 0.0$.
C IUNIT  :IN       Integer.
C         Unit number on which to write the error at each iteration, 
C         if this is desired for monitoring convergence.  If unit 
C         number is 0, no writing will occur.
C R      :WORK     Double Precision R(N).
C Z      :WORK     Double Precision Z(N).
C P      :WORK     Double Precision P(N).
C DZ     :WORK     Double Precision DZ(N).
C RWORK  :WORK     Double Precision RWORK(USER DEFINABLE).
C         Double Precision array that can be used by  MSOLVE.
C IWORK  :WORK     Integer IWORK(USER DEFINABLE).
C         Integer array that can be used by  MSOLVE.
C
C *Description
C       This routine does  not care  what matrix data   structure is
C       used for  A and M.  It simply   calls  the MATVEC and MSOLVE
C       routines, with  the arguments as  described above.  The user
C       could write any type of structure and the appropriate MATVEC
C       and MSOLVE routines.  It is assumed  that A is stored in the
C       IA, JA, A  arrays in some fashion and  that M (or INV(M)) is
C       stored  in  IWORK  and  RWORK   in  some fashion.   The SLAP
C       routines DSDCG and DSICCG are examples of this procedure.
C       
C       Two  examples  of  matrix  data structures  are the: 1) SLAP
C       Triad  format and 2) SLAP Column format.
C       
C       =================== S L A P Triad format ===================
C               
C       In  this   format only the  non-zeros are  stored.  They may
C       appear  in *ANY* order.   The user  supplies three arrays of
C       length NELT, where  NELT  is the number  of non-zeros in the
C       matrix:  (IA(NELT), JA(NELT),  A(NELT)).  For each  non-zero
C       the  user puts   the row  and  column index   of that matrix
C       element in the IA and JA arrays.  The  value of the non-zero
C       matrix  element is  placed in  the corresponding location of
C       the A  array.  This is  an extremely easy data  structure to
C       generate.  On  the other hand it  is  not too  efficient  on
C       vector  computers   for the  iterative  solution  of  linear
C       systems.  Hence, SLAP  changes this input  data structure to
C       the SLAP   Column  format for the  iteration (but   does not
C       change it back).
C       
C       Here is an example of the  SLAP Triad   storage format for a
C       5x5 Matrix.  Recall that the entries may appear in any order.
C
C           5x5 Matrix      SLAP Triad format for 5x5 matrix on left.
C                              1  2  3  4  5  6  7  8  9 10 11
C       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21
C       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2
C       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C
C       =================== S L A P Column format ==================
C       This routine requires  that the  matrix  A be  stored in the
C       SLAP Column format.  In this format the non-zeros are stored
C       counting down columns (except for  the diagonal entry, which
C       must appear first in each  "column")  and are stored  in the
C       double precision array A.   In other words,  for each column
C       in the matrix put the diagonal entry in  A.  Then put in the
C       other non-zero  elements going down  the column (except  the
C       diagonal) in order.   The  IA array holds the  row index for
C       each non-zero.  The JA array holds the offsets  into the IA,
C       A arrays  for  the  beginning  of each   column.   That  is,
C       IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the
C       ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1),
C       A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column.
C       Note that we always have  JA(N+1) = NELT+1,  where N is  the
C       number of columns in  the matrix and NELT  is the number  of
C       non-zeros in the matrix.
C       
C       Here is an example of the  SLAP Column  storage format for a
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a 
C       column):
C       
C       5x5 Matrix      SLAP Column format for 5x5 matrix on left.
C       1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C       
C *Precision:           Double Precision
C *See Also:
C       DSDCG, DSICCG
C***REFERENCES  1. Louis Hageman \& David Young, ``Applied Iterative 
C                 Methods'', Academic Press, New York (1981) ISBN 
C                 0-12-313340-8.
C
C               2. Concus, Golub \& O'Leary, ``A Generalized Conjugate 
C                 Gradient Method for the Numerical Solution of 
C                 Elliptic Partial Differential Equations,'' in Sparse 
C                 Matrix Computations (Bunch \& Rose, Eds.), Academic 
C                 Press, New York (1979).
C***ROUTINES CALLED  MATVEC, MSOLVE, ISDCG, DCOPY, DDOT, DAXPY, D1MACH
C***END PROLOGUE  DCG
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX, ITER
      INTEGER IUNIT, IERR, IWORK(*)
      DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, R(N), Z(N), P(N)
      DOUBLE PRECISION DZ(N), RWORK(*)
      EXTERNAL MATVEC, MSOLVE
C
C         Check some of the input data.
C***FIRST EXECUTABLE STATEMENT  DCG
      ITER = 0
      IERR = 0
      IF( N.LT.1 ) THEN
         IERR = 3
         RETURN
      ENDIF
      TOLMIN = 500.0*D1MACH(3)
      IF( TOL.LT.TOLMIN ) THEN
         TOL = TOLMIN
         IERR = 4
      ENDIF
C         
C         Calculate initial residual and pseudo-residual, and check
C         stopping criterion.
      CALL MATVEC(N, X, R, NELT, IA, JA, A, ISYM)
      DO 10 I = 1, N
         R(I) = B(I) - R(I)
 10   CONTINUE
      CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
C         
      IF( ISDCG(N, B, X, NELT, IA, JA, A, ISYM, MSOLVE, ITOL, TOL,
     $     ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, DZ,
     $     RWORK, IWORK, AK, BK, BNRM, SOLNRM) .NE. 0 ) GO TO 200
      IF( IERR.NE.0 ) RETURN
C         
C         ***** Iteration loop *****
C         
      DO 100 K=1,ITMAX
         ITER = K
C         
C         Calculate coefficient bk and direction vector p.
         BKNUM = DDOT(N, Z, 1, R, 1)
         IF( BKNUM.LE.0.0D0 ) THEN
            IERR = 5
            RETURN
         ENDIF
         IF(ITER .EQ. 1) THEN
            CALL DCOPY(N, Z, 1, P, 1)
         ELSE
            BK = BKNUM/BKDEN
            DO 20 I = 1, N
               P(I) = Z(I) + BK*P(I)
 20         CONTINUE
         ENDIF
         BKDEN = BKNUM
C         
C         Calculate coefficient ak, new iterate x, new residual r,
C         and new pseudo-residual z.
         CALL MATVEC(N, P, Z, NELT, IA, JA, A, ISYM)
         AKDEN = DDOT(N, P, 1, Z, 1)
         IF( AKDEN.LE.0.0D0 ) THEN
            IERR = 6
            RETURN
         ENDIF
         AK = BKNUM/AKDEN
         CALL DAXPY(N, AK, P, 1, X, 1)
         CALL DAXPY(N, -AK, Z, 1, R, 1)
         CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
C         
C         check stopping criterion.
         IF( ISDCG(N, B, X, NELT, IA, JA, A, ISYM, MSOLVE, ITOL, TOL,
     $        ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, DZ, RWORK,
     $        IWORK, AK, BK, BNRM, SOLNRM) .NE. 0 ) GO TO 200
C         
 100  CONTINUE
C         
C         *****   end of loop  *****
C         
C         stopping criterion not satisfied.
      ITER = ITMAX + 1
      IERR = 2
C         
 200  RETURN
      END
*DECK DSDCG
      SUBROUTINE DSDCG(N, B, X, NELT, IA, JA, A, ISYM, ITOL, TOL,
     $     ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW )
C***BEGIN PROLOGUE  DSDCG
C***DATE WRITTEN   890404   (YYMMDD)
C***REVISION DATE  890404   (YYMMDD)
C***CATEGORY NO.  D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=DOUBLE PRECISION(DSDCG-D),
C             Symmetric Linear system, Sparse, Iterative Precondition
C***AUTHOR  Greenbaum, Anne, Courant Institute
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  Diagonally Scaled Conjugate Gradient Sparse Ax=b Solver.
C            Routine to solve a  symmetric positive definite  linear
C            system Ax  =  b  using   the  Preconditioned  Conjugate
C            Gradient   method.    The preconditioner  is   diagonal
C            scaling.
C***DESCRIPTION
C *Usage:
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX
C     INTEGER ITER, IERR, IUNIT, LENW, IWORK(10), LENIW
C     DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, RWORK(5*N)
C
C     CALL DSDCG(N, B, X, NELT, IA, JA, A, ISYM, ITOL, TOL,
C    $     ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW )
C
C *Arguments:
C N      :IN       Integer.
C         Order of the Matrix.
C B      :IN       Double Precision B(N).
C         Right-hand side vector.
C X      :INOUT    Double Precision X(N).
C         On input X is your initial guess for solution vector.
C         On output X is the final approximate solution.
C NELT   :IN       Integer.
C         Number of Non-Zeros stored in A.
C IA     :INOUT    Integer IA(NELT).
C JA     :INOUT    Integer JA(NELT).
C A      :INOUT    Integer A(NELT).
C         These arrays should hold the matrix A in either the SLAP
C         Triad format or the SLAP Column format.  See ``Description'', 
C         below.  If the SLAP Triad format is chosen it is changed 
C         internally to the SLAP Column format.
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all nonzero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the upper
C         or lower triangle of the matrix is stored.
C ITOL   :IN       Integer.
C         Flag to indicate type of convergence criterion.
C         If ITOL=1, iteration stops when the 2-norm of the residual 
C         divided by the 2-norm of the right-hand side is less than TOL.
C         If ITOL=2, iteration stops when the 2-norm of M-inv times the 
C         residual divided by the 2-norm of M-inv times the right hand 
C         side is less than TOL, where M-inv is the inverse of the 
C         diagonal of A.
C         ITOL=11 is often useful for checking and comparing different 
C         routines.  For this case, the user must supply the "exact" 
C         solution or a very accurate approximation (one with an error 
C         much less than TOL) through a common block,
C                     COMMON /SOLBLK/ SOLN(1)
C         if ITOL=11, iteration stops when the 2-norm of the difference 
C         between the iterative approximation and the user-supplied
C         solution divided by the 2-norm of the user-supplied solution 
C         is less than TOL.  Note that this requires the user to set up
C         the "COMMON /SOLBLK/ SOLN(LENGTH)" in the calling routine. 
C         The routine with this declaration should be loaded before the
C         stop test so that the correct length is used by the loader.  
C         This procedure is not standard Fortran and may not work 
C         correctly on your system (although it has worked on every
C         system the authors have tried).  If ITOL is not 11 then this
C         common block is indeed standard Fortran.
C TOL    :IN       Double Precision.
C         Convergence criterion, as described above.
C ITMAX  :IN       Integer.
C         Maximum number of iterations.
C ITER   :OUT      Integer.
C         Number of iterations required to reach convergence, or 
C         ITMAX+1 if convergence criterion could not be achieved in 
C         ITMAX iterations.
C ERR    :OUT      Double Precision.
C         Error estimate of error in final approximate solution, as 
C         defined by ITOL.
C IERR   :OUT      Integer.
C         Return error flag.
C           IERR = 0 => All went well.
C           IERR = 1 => Insufficient storage allocated 
C                       for WORK or IWORK.
C           IERR = 2 => Method failed to converge in 
C                       ITMAX steps.
C           IERR = 3 => Error in user input.  Check input
C                       value of N, ITOL.
C           IERR = 4 => User error tolerance set too tight.
C                       Reset to 500.0*D1MACH(3).  Iteration proceeded.
C           IERR = 5 => Preconditioning matrix, M,  is not 
C                       Positive Definite.  $(r,z) < 0.0$.
C           IERR = 6 => Matrix A is not Positive Definite.  
C                       $(p,Ap) < 0.0$.
C IUNIT  :IN       Integer.
C         Unit number on which to write the error at each iteration, 
C         if this is desired for monitoring convergence.  If unit 
C         number is 0, no writing will occur.
C RWORK  :WORK     Double Precision RWORK(LENW).
C         Double Precision array used for workspace.
C LENW   :IN       Integer.
C         Length of the double precision workspace, RWORK.  LENW >= 5*N.
C IWORK  :WORK     Integer IWORK(LENIW).
C         Used to hold pointers into the double precision workspace, 
C         RWORK.  Upon return the following locations of IWORK hold 
C         information which may be of use to the user:
C         IWORK(9)  Amount of Integer workspace actually used.
C         IWORK(10) Amount of Double Precision workspace actually used.
C LENIW  :IN       Integer.
C         Length of the integer workspace, IWORK.  LENIW >= 10.
C
C *Description:
C       This  routine   performs preconditioned conjugate   gradient
C       method on  the  symmetric positive definite   linear  system
C       Ax=b.   The preconditioner is  M = DIAG(A), the  diagonal of
C       the matrix A.  This is the  simplest of preconditioners  and
C       vectorizes very well.   This routine is  simply a driver for
C       the DCG routine.  It  calls the DSDS  routine to  set up the
C       preconditioning  and  then  calls  DCG  with the appropriate 
C       MATVEC and MSOLVE routines.
C
C       The Sparse Linear Algebra Package (SLAP) utilizes two matrix
C       data structures: 1) the  SLAP Triad  format or  2)  the SLAP
C       Column format.  The user can hand this routine either of the
C       of these data structures and SLAP  will figure out  which on
C       is being used and act accordingly.
C       
C       =================== S L A P Triad format ===================
C
C       This routine requires that the  matrix A be   stored in  the
C       SLAP  Triad format.  In  this format only the non-zeros  are
C       stored.  They may appear in  *ANY* order.  The user supplies
C       three arrays of  length NELT, where  NELT is  the number  of
C       non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For
C       each non-zero the user puts the row and column index of that
C       matrix element  in the IA and  JA arrays.  The  value of the
C       non-zero   matrix  element is  placed  in  the corresponding
C       location of the A array.   This is  an  extremely  easy data
C       structure to generate.  On  the  other hand it   is  not too
C       efficient on vector computers for  the iterative solution of
C       linear systems.  Hence,   SLAP changes   this  input    data
C       structure to the SLAP Column format  for  the iteration (but
C       does not change it back).
C       
C       Here is an example of the  SLAP Triad   storage format for a
C       5x5 Matrix.  Recall that the entries may appear in any order.
C
C           5x5 Matrix       SLAP Triad format for 5x5 matrix on left.
C                              1  2  3  4  5  6  7  8  9 10 11
C       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21
C       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2
C       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C       
C       =================== S L A P Column format ==================
C       This routine requires  that the  matrix  A be  stored in the
C       SLAP Column format.  In this format the non-zeros are stored
C       counting down columns (except for  the diagonal entry, which
C       must appear first in each  "column")  and are stored  in the
C       double precision array A.   In other words,  for each column
C       in the matrix put the diagonal entry in  A.  Then put in the
C       other non-zero  elements going down  the column (except  the
C       diagonal) in order.   The  IA array holds the  row index for
C       each non-zero.  The JA array holds the offsets  into the IA,
C       A arrays  for  the  beginning  of each   column.   That  is,
C       IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the
C       ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1),
C       A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column.
C       Note that we always have  JA(N+1) = NELT+1,  where N is  the
C       number of columns in  the matrix and NELT  is the number  of
C       non-zeros in the matrix.
C       
C       Here is an example of the  SLAP Column  storage format for a
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a 
C       column):
C       
C       5x5 Matrix      SLAP Column format for 5x5 matrix on left.
C       1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C
C *Precision:           Double Precision
C *Side Effects:
C       The SLAP Triad format (IA, JA, A) is modified internally  to 
C       be the   SLAP   Column format.  See above.
C       
C *See Also:
C       DCG, DSICCG
C***REFERENCES  1. Louis Hageman \& David Young, ``Applied Iterative 
C                 Methods'', Academic Press, New York (1981) ISBN 
C                 0-12-313340-8.
C               2. Concus, Golub \& O'Leary, ``A Generalized Conjugate 
C                 Gradient Method for the Numerical Solution of 
C                 Elliptic Partial Differential Equations,'' in Sparse 
C                 Matrix Computations (Bunch \& Rose, Eds.), Academic 
C                 Press, New York (1979).
C***ROUTINES CALLED  DS2Y, DCHKW, DSDS, DCG
C***END PROLOGUE  DSDCG
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, ITOL
      INTEGER ITMAX, ITER, IERR, IUNIT, LENW, IWORK(LENIW), LENIW
      DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, RWORK(LENW)
      EXTERNAL DSMV, DSDI
      PARAMETER (LOCRB=1, LOCIB=11)
C
C         Modify the SLAP matrix data structure to YSMP-Column.
C***FIRST EXECUTABLE STATEMENT  DSDCG
      IERR = 0
      IF( N.LT.1 .OR. NELT.LT.1 ) THEN
         IERR = 3
         RETURN
      ENDIF
      CALL DS2Y( N, NELT, IA, JA, A, ISYM )
C         
C         Set up the work arrays.
C         Compute the inverse of the diagonal of the matrix.  This 
C         will be used as the preconditioner.
      LOCIW = LOCIB
C
      LOCD = LOCRB
      LOCR = LOCD + N
      LOCZ = LOCR + N
      LOCP = LOCZ + N
      LOCDZ = LOCP + N
      LOCW  = LOCDZ + N
C
C         Check the workspace allocations.
      CALL DCHKW( 'DSDCG', LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR )
      IF( IERR.NE.0 ) RETURN
C
      IWORK(4) = LOCD
      IWORK(9) = LOCIW
      IWORK(10) = LOCW
C
      CALL DSDS(N, NELT, IA, JA, A, ISYM, RWORK(LOCD))
C
C         Do the Preconditioned Conjugate Gradient.
      CALL DCG(N, B, X, NELT, IA, JA, A, ISYM, DSMV, DSDI, 
     $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK(LOCR), 
     $     RWORK(LOCZ), RWORK(LOCP), RWORK(LOCDZ), RWORK, IWORK)
      RETURN
C------------- LAST LINE OF DSDCG FOLLOWS -----------------------------
      END
*DECK DSICCG
      SUBROUTINE DSICCG(N, B, X, NELT, IA, JA, A, ISYM, ITOL, TOL,
     $     ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW )
C***BEGIN PROLOGUE  DSICCG
C***DATE WRITTEN   890404   (YYMMDD)
C***REVISION DATE  890404   (YYMMDD)
C***CATEGORY NO.  D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=DOUBLE PRECISION(DSICCG-D),
C             Symmetric Linear system, Sparse, 
C             Iterative Precondition, Incomplete Cholesky
C***AUTHOR  Greenbaum, Anne, Courant Institute
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  Incomplete Cholesky Conjugate Gradient Sparse Ax=b Solver.
C            Routine to  solve a symmetric  positive definite linear
C            system   Ax    =  b  using  the    incomplete  Cholesky
C            Preconditioned Conjugate Gradient method.
C
C***DESCRIPTION
C *Usage:
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX
C     INTEGER ITER, IERR, IUNIT, LENW, IWORK(NEL+2*n+1), LENIW
C     DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, RWORK(NEL+5*N)
C
C     CALL DSICCG(N, B, X, NELT, IA, JA, A, ISYM, ITOL, TOL,
C    $     ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW )
C
C *Arguments:
C N      :IN       Integer.
C         Order of the Matrix.
C B      :IN       Double Precision B(N).
C         Right-hand side vector.
C X      :INOUT    Double Precision X(N).
C         On input X is your initial guess for solution vector.
C         On output X is the final approximate solution.
C NELT   :IN       Integer.
C         Number of Non-Zeros stored in A.
C IA     :INOUT    Integer IA(NELT).
C JA     :INOUT    Integer JA(NELT).
C A      :INOUT    Integer A(NELT).
C         These arrays should hold the matrix A in either the SLAP
C         Triad format or the SLAP Column format.  See ``Description'', 
C         below.  If the SLAP Triad format is chosen it is changed 
C         internally to the SLAP Column format.
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all nonzero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the upper
C         or lower triangle of the matrix is stored.
C ITOL   :IN       Integer.
C         Flag to indicate type of convergence criterion.
C         If ITOL=1, iteration stops when the 2-norm of the residual 
C         divided by the 2-norm of the right-hand side is less than TOL.
C         If ITOL=2, iteration stops when the 2-norm of M-inv times the 
C         residual divided by the 2-norm of M-inv times the right hand 
C         side is less than TOL, where M-inv is the inverse of the 
C         diagonal of A.
C         ITOL=11 is often useful for checking and comparing different 
C         routines.  For this case, the user must supply the "exact" 
C         solution or a very accurate approximation (one with an error 
C         much less than TOL) through a common block,
C                     COMMON /SOLBLK/ SOLN(1)
C         if ITOL=11, iteration stops when the 2-norm of the difference 
C         between the iterative approximation and the user-supplied
C         solution divided by the 2-norm of the user-supplied solution 
C         is less than TOL.  Note that this requires the user to set up
C         the "COMMON /SOLBLK/ SOLN(LENGTH)" in the calling routine. 
C         The routine with this declaration should be loaded before the
C         stop test so that the correct length is used by the loader.  
C         This procedure is not standard Fortran and may not work 
C         correctly on your system (although it has worked on every
C         system the authors have tried).  If ITOL is not 11 then this
C         common block is indeed standard Fortran.
C TOL    :IN       Double Precision.
C         Convergence criterion, as described above.
C ITMAX  :IN       Integer.
C         Maximum number of iterations.
C ITER   :OUT      Integer.
C         Number of iterations required to reach convergence, or 
C         ITMAX+1 if convergence criterion could not be achieved in 
C         ITMAX iterations.
C ERR    :OUT      Double Precision.
C         Error estimate of error in final approximate solution, as 
C         defined by ITOL.
C IERR   :OUT      Integer.
C         Return error flag.
C           IERR = 0 => All went well.
C           IERR = 1 => Insufficient storage allocated 
C                       for WORK or IWORK.
C           IERR = 2 => Method failed to converge in 
C                       ITMAX steps.
C           IERR = 3 => Error in user input.  Check input
C                       value of N, ITOL.
C           IERR = 4 => User error tolerance set too tight.
C                       Reset to 500.0*D1MACH(3).  Iteration proceeded.
C           IERR = 5 => Preconditioning matrix, M,  is not 
C                       Positive Definite.  $(r,z) < 0.0$.
C           IERR = 6 => Matrix A is not Positive Definite.  
C                       $(p,Ap) < 0.0$.
C           IERR = 7 => Incomplete factorization broke down
C                       and was fudged.  Resulting preconditioning may
C                       be less than the best.
C IUNIT  :IN       Integer.
C         Unit number on which to write the error at each iteration, 
C         if this is desired for monitoring convergence.  If unit 
C         number is 0, no writing will occur.
C RWORK  :WORK     Double Precision RWORK(LENW).
C         Double Precision array used for workspace.  NEL is the 
C         number of non-
C         zeros in the lower triangle of the matrix (including the
C         diagonal)
C LENW   :IN       Integer.
C         Length of the double precision workspace, RWORK.  
C         LENW >= NEL+5*N.
C IWORK  :WORK     Integer IWORK(LENIW).
C         Integer array used for workspace.  NEL is the number of non-
C         zeros in the lower triangle of the matrix (including the
C         diagonal).
C         Upon return the following locations of IWORK hold information
C         which may be of use to the user:
C         IWORK(9)  Amount of Integer workspace actually used.
C         IWORK(10) Amount of Double Precision workspace actually used.
C LENIW  :IN       Integer.
C         Length of the integer workspace, IWORK.  LENIW >= NEL+N+11.
C
C *Description:
C       This routine  performs  preconditioned  conjugate   gradient
C       method on the   symmetric positive  definite  linear  system
C       Ax=b.   The preconditioner  is  the incomplete Cholesky (IC)
C       factorization of the matrix A.  See  DSICS for details about
C       the incomplete   factorization algorithm.  One   should note
C       here however, that the  IC factorization is a  slow  process
C       and  that  one should   save  factorizations  for  reuse, if
C       possible.  The   MSOLVE operation (handled  in  DSLLTI) does
C       vectorize on machines  with  hardware  gather/scatter and is
C       quite fast.
C
C       The Sparse Linear Algebra Package (SLAP) utilizes two matrix
C       data structures: 1) the  SLAP Triad  format or  2)  the SLAP
C       Column format.  The user can hand this routine either of the
C       of these data structures and SLAP  will figure out  which on
C       is being used and act accordingly.
C       
C       =================== S L A P Triad format ===================
C
C       This routine requires that the  matrix A be   stored in  the
C       SLAP  Triad format.  In  this format only the non-zeros  are
C       stored.  They may appear in  *ANY* order.  The user supplies
C       three arrays of  length NELT, where  NELT is  the number  of
C       non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For
C       each non-zero the user puts the row and column index of that
C       matrix element  in the IA and  JA arrays.  The  value of the
C       non-zero   matrix  element is  placed  in  the corresponding
C       location of the A array.   This is  an  extremely  easy data
C       structure to generate.  On  the  other hand it   is  not too
C       efficient on vector computers for  the iterative solution of
C       linear systems.  Hence,   SLAP changes   this  input    data
C       structure to the SLAP Column format  for  the iteration (but
C       does not change it back).
C       
C       Here is an example of the  SLAP Triad   storage format for a
C       5x5 Matrix.  Recall that the entries may appear in any order.
C
C           5x5 Matrix       SLAP Triad format for 5x5 matrix on left.
C                              1  2  3  4  5  6  7  8  9 10 11
C       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21
C       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2
C       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C       
C       =================== S L A P Column format ==================
C       This routine requires  that the  matrix  A be  stored in the
C       SLAP Column format.  In this format the non-zeros are stored
C       counting down columns (except for  the diagonal entry, which
C       must appear first in each  "column")  and are stored  in the
C       double precision array A.   In other words,  for each column
C       in the matrix put the diagonal entry in  A.  Then put in the
C       other non-zero  elements going down  the column (except  the
C       diagonal) in order.   The  IA array holds the  row index for
C       each non-zero.  The JA array holds the offsets  into the IA,
C       A arrays  for  the  beginning  of each   column.   That  is,
C       IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the
C       ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1),
C       A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column.
C       Note that we always have  JA(N+1) = NELT+1,  where N is  the
C       number of columns in  the matrix and NELT  is the number  of
C       non-zeros in the matrix.
C       
C       Here is an example of the  SLAP Column  storage format for a
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a 
C       column):
C       
C       5x5 Matrix      SLAP Column format for 5x5 matrix on left.
C       1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C       
C *Precision:           Double Precision
C *Side Effects:
C       The SLAP Triad format (IA, JA, A) is modified internally to be
C       the SLAP Column format.  See above.
C       
C *See Also:
C       DCG, DSLLTI
C***REFERENCES  1. Louis Hageman \& David Young, ``Applied Iterative 
C                 Methods'', Academic Press, New York (1981) ISBN 
C                 0-12-313340-8.
C               2. Concus, Golub \& O'Leary, ``A Generalized Conjugate 
C                 Gradient Method for the Numerical Solution of 
C                 Elliptic Partial Differential Equations,'' in Sparse 
C                 Matrix Computations (Bunch \& Rose, Eds.), Academic 
C                 Press, New York (1979).
C***ROUTINES CALLED  DS2Y, DCHKW, DSICS, XERRWV, DCG
C***END PROLOGUE  DSICCG
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, ITOL
      INTEGER ITMAX, ITER, IUNIT, LENW, IWORK(LENIW), LENIW
      DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, RWORK(LENW)
      EXTERNAL DSMV, DSLLTI
      PARAMETER (LOCRB=1, LOCIB=11)
C         
C         Change the SLAP input matrix IA, JA, A to SLAP-Column format.
C***FIRST EXECUTABLE STATEMENT  DSICCG
      IERR = 0
      IF( N.LT.1 .OR. NELT.LT.1 ) THEN
         IERR = 3
         RETURN
      ENDIF
      CALL DS2Y( N, NELT, IA, JA, A, ISYM )
C         
C         Count number of elements in lower triangle of the matrix.
C         Then set up the work arrays.
      IF( ISYM.EQ.0 ) THEN
         NEL = (NELT + N)/2
      ELSE
         NEL = NELT
      ENDIF
C
      LOCJEL = LOCIB
      LOCIEL = LOCJEL + NEL
      LOCIW = LOCIEL + N + 1
C
      LOCEL = LOCRB
      LOCDIN = LOCEL + NEL
      LOCR = LOCDIN + N
      LOCZ = LOCR + N
      LOCP = LOCZ + N
      LOCDZ = LOCP + N
      LOCW = LOCDZ + N
C
C         Check the workspace allocations.
      CALL DCHKW( 'DSICCG', LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR )
      IF( IERR.NE.0 ) RETURN
C
      IWORK(1) = NEL
      IWORK(2) = LOCJEL
      IWORK(3) = LOCIEL
      IWORK(4) = LOCEL
      IWORK(5) = LOCDIN
      IWORK(9) = LOCIW
      IWORK(10) = LOCW
C
C         Compute the Incomplete Cholesky decomposition.
C
      CALL DSICS(N, NELT, IA, JA, A, ISYM, NEL, IWORK(LOCIEL),
     $     IWORK(LOCJEL), RWORK(LOCEL), RWORK(LOCDIN),
     $     RWORK(LOCR), IERR )
      IF( IERR.NE.0 ) THEN
         CALL XERRWV('DSICCG: Warning...IC factorization broke down '//
     $        'on step i1.  Diagonal was set to unity and '//
     $        'factorization proceeded.', 113, 1, 1, 1, IERR, 0,
     $        0, 0.0, 0.0 )
         IERR = 7
      ENDIF
C
C         Do the Preconditioned Conjugate Gradient.
      CALL DCG(N, B, X, NELT, IA, JA, A, ISYM, DSMV, DSLLTI, 
     $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK(LOCR),
     $     RWORK(LOCZ), RWORK(LOCP), RWORK(LOCDZ), RWORK(1),
     $     IWORK(1))
      RETURN
C------------- LAST LINE OF DSICCG FOLLOWS ----------------------------
      END
*DECK ISDCG
      FUNCTION ISDCG(N, B, X, NELT, IA, JA, A, ISYM, MSOLVE, ITOL,
     $     TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, DZ, 
     $     RWORK, IWORK, AK, BK, BNRM, SOLNRM)
C***BEGIN PROLOGUE  ISDCG
C***REFER TO  DCG, DSDCG, DSICCG
C***DATE WRITTEN   890404   (YYMMDD)
C***REVISION DATE  890404   (YYMMDD)
C***CATEGORY NO.  D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=DOUBLE PRECISION(ISDCG-D),
C             Linear system, Sparse, Stop Test
C***AUTHOR  Greenbaum, Anne, Courant Institute
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  Preconditioned Conjugate Gradient Stop Test.
C            This routine calculates the stop test for the Conjugate
C            Gradient iteration scheme.  It returns a nonzero if the
C            error estimate (the type of which is determined by ITOL)
C            is less than the user specified tolerance TOL.
C***DESCRIPTION
C *Usage:
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX, ITER
C     INTEGER IERR, IUNIT, IWORK(USER DEFINED)
C     DOUBLE PRECISION B(N), X(N), A(N), TOL, ERR, R(N), Z(N)
C     DOUBLE PRECISION P(N), DZ(N), RWORK(USER DEFINED), AK, BK
C     DOUBLE PRECISION BNRM, SOLNRM
C     EXTERNAL MSOLVE
C
C     IF( ISDCG(N, B, X, NELT, IA, JA, A, ISYM, MSOLVE, ITOL, TOL,
C    $     ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, DZ, RWORK, IWORK,
C    $     AK, BK, BNRM, SOLNRM) .NE. 0 ) THEN ITERATION DONE
C
C *Arguments:
C N      :IN       Integer.
C         Order of the Matrix.
C B      :IN       Double Precision B(N).
C         Right-hand side vector.
C X      :IN       Double Precision X(N).
C         The current approximate solution vector.
C NELT   :IN       Integer.
C         Number of Non-Zeros stored in A.
C IA     :IN       Integer IA(NELT).
C JA     :IN       Integer JA(NELT).
C A      :IN       Double Precision A(NELT).
C         These arrays should hold the matrix A in either the SLAP
C         Triad format or the SLAP Column format.  See ``Description''
C         in the DCG, DSDCG or DSICCG routines.
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all nonzero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the upper
C         or lower triangle of the matrix is stored.
C MSOLVE :EXT      External.
C         Name of a routine which solves a linear system MZ = R for
C         Z given R with the preconditioning matrix M (M is supplied via
C         RWORK and IWORK arrays).  The name of the MSOLVE routine must 
C         be declared external in the calling program.  The calling 
C         sequence to MSLOVE is:
C             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
C         Where N is the number of unknowns, R is the right-hand side 
C         vector, and Z is the solution upon return.  RWORK is a double 
C         precision 
C         array that can be used to pass necessary preconditioning 
C         information and/or workspace to MSOLVE.  IWORK is an integer 
C         work array for the same purpose as RWORK.
C ITOL   :IN       Integer.
C         Flag to indicate type of convergence criterion.
C         If ITOL=1, iteration stops when the 2-norm of the residual 
C         divided by the 2-norm of the right-hand side is less than TOL.
C         If ITOL=2, iteration stops when the 2-norm of M-inv times the 
C         residual divided by the 2-norm of M-inv times the right hand 
C         side is less than tol, where M-inv is the inverse of the 
C         diagonal of A.
C         ITOL=11 is often useful for checking and comparing different 
C         routines.  For this case, the user must supply the ``exact'' 
C         solution or a very accurate approximation (one with an error 
C         much less than tol) through a common block,
C         COMMON /SOLBLK/ SOLN( )
C         if ITOL=11, iteration stops when the 2-norm of the difference 
C         between the iterative approximation and the user-supplied
C         solution divided by the 2-norm of the user-supplied solution 
C         is less than tol.
C TOL    :IN       Double Precision.
C         Convergence criterion, as described above.
C ITMAX  :IN       Integer.
C         Maximum number of iterations.
C ITER   :IN       Integer.
C         The iteration for which to check for convergence.
C ERR    :OUT      Double Precision.
C         Error estimate of error in the X(N) approximate solution, as 
C         defined by ITOL.
C IERR   :OUT      Integer.
C         Error flag.  IERR is set to 3 if ITOL is not on of the 
C         acceptable values, see above. 
C IUNIT  :IN       Integer.
C         Unit number on which to write the error at each iteration, 
C         if this is desired for monitoring convergence.  If unit 
C         number is 0, no writing will occur.
C R      :IN       Double Precision R(N).
C         The residual R = B-AX.
C Z      :WORK     Double Precision Z(N).
C         Workspace used to hold the pseudo-residual M Z = R.
C P      :IN       Double Precision P(N).
C         The conjugate direction vector.
C DZ     :WORK     Double Precision DZ(N).
C         Workspace used to hold temporary vector(s).
C RWORK  :WORK     Double Precision RWORK(USER DEFINABLE).
C         Double Precision array that can be used by MSOLVE.
C IWORK  :WORK     Integer IWORK(USER DEFINABLE).
C         Integer array that can be used by MSOLVE.
C BNRM   :INOUT    Double Precision.
C         Norm of the right hand side.  Type of norm depends on ITOL.
C         Calculated only on the first call.
C SOLNRM :INOUT    Double Precision.
C         2-Norm of the true solution, SOLN.  Only computed and used
C         if ITOL = 11.
C
C *Function Return Values:
C       0 : Error estimate (determined by ITOL) is *NOT* less than the 
C           specified tolerance, TOL.  The iteration must continue.
C       1 : Error estimate (determined by ITOL) is less than the 
C           specified tolerance, TOL.  The iteration can be considered
C           complete.
C
C *Precision:           Double Precision
C *See Also:
C       DCG, DSDCG, DSICCG
C
C *Cautions:
C     This routine will attempt to write to the fortran logical output 
C     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that
C     this  logical  unit  must  be  attached  to  a  file or terminal
C     before calling this routine with a non-zero  value  for   IUNIT.
C     This routine does not check for the validity of a non-zero IUNIT
C     unit number.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  MSOLVE, DNRM2
C***COMMON BLOCKS    SOLBLK
C***END PROLOGUE  ISDCG
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX
      INTEGER ITER, IERR, IUNIT, IWORK(*)
      DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, R(N)
      DOUBLE PRECISION Z(N), P(N), DZ(N), RWORK(*)
      EXTERNAL MSOLVE
      COMMON /SOLBLK/ SOLN(1)
C         
C***FIRST EXECUTABLE STATEMENT  ISDCG
      ISDCG = 0
C         
      IF( ITOL.EQ.1 ) THEN
C         err = ||Residual||/||RightHandSide|| (2-Norms).
         IF(ITER .EQ. 0) BNRM = DNRM2(N, B, 1)
         ERR = DNRM2(N, R, 1)/BNRM
      ELSE IF( ITOL.EQ.2 ) THEN
C                  -1              -1
C         err = ||M  Residual||/||M  RightHandSide|| (2-Norms).
         IF(ITER .EQ. 0) THEN
            CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
            BNRM = DNRM2(N, DZ, 1)
         ENDIF
         ERR = DNRM2(N, Z, 1)/BNRM
      ELSE IF( ITOL.EQ.11 ) THEN
C         err = ||x-TrueSolution||/||TrueSolution|| (2-Norms).
         IF(ITER .EQ. 0) SOLNRM = DNRM2(N, SOLN, 1)
         DO 10 I = 1, N
            DZ(I) = X(I) - SOLN(I)
 10      CONTINUE
         ERR = DNRM2(N, DZ, 1)/SOLNRM
      ELSE
C
C         If we get here ITOL is not one of the acceptable values.
         ERR = 1.0E10
         IERR = 3
      ENDIF
C         
      IF(IUNIT .NE. 0) THEN
         IF( ITER.EQ.0 ) THEN
            WRITE(IUNIT,1000) N, ITOL
         ENDIF
         WRITE(IUNIT,1010) ITER, ERR, AK, BK
      ENDIF     
      IF(ERR .LE. TOL) ISDCG = 1
      RETURN
 1000 FORMAT(' Preconditioned Conjugate Gradient for ',
     $     'N, ITOL = ',I5, I5,
     $     /' ITER','   Error Estimate','            Alpha',
     $     '             Beta')
 1010 FORMAT(1X,I4,1X,E16.7,1X,E16.7,1X,E16.7)
C------------- LAST LINE OF ISDCG FOLLOWS ------------------------------
      END
