*DECK DSDS
      SUBROUTINE DSDS(N, NELT, IA, JA, A, ISYM, DINV)
C***BEGIN PROLOGUE  DSDS
C***DATE WRITTEN   890404   (YYMMDD)
C***REVISION DATE  890404   (YYMMDD)
C***CATEGORY NO.  D2A4, D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=DOUBLE PRECISION(DSDS-D),
C             SLAP Sparse, Diagonal
C***AUTHOR  Greenbaum, Anne, Courant Institute
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  Diagonal Scaling Preconditioner SLAP Set Up.
C            Routine to compute the inverse of the diagonal of a matrix
C            stored in the SLAP Column format.
C***DESCRIPTION
C *Usage:
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
C     DOUBLE PRECISION A(NELT), DINV(N)
C
C     CALL DSDS( N, NELT, IA, JA, A, ISYM, DINV )
C
C *Arguments:
C N      :IN       Integer.
C         Order of the Matrix.
C NELT   :IN       Integer.
C         Number of elements in arrays IA, JA, and A.
C IA     :INOUT    Integer IA(NELT).
C JA     :INOUT    Integer JA(NELT).
C A      :INOUT    Double Precision A(NELT).
C         These arrays should hold the matrix A in the SLAP Column
C         format.  See "Description", below. 
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all nonzero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the upper
C         or lower triangle of the matrix is stored.
C DINV   :OUT      Double Precision DINV(N).
C         Upon return this array holds 1./DIAG(A).
C
C *Description
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
C       With the SLAP  format  all  of  the   "inner  loops" of this
C       routine should vectorize  on  machines with hardware support
C       for vector   gather/scatter  operations.  Your compiler  may
C       require a compiler directive to  convince it that  there are
C       no  implicit  vector  dependencies.  Compiler directives for
C       the Alliant    FX/Fortran and CRI   CFT/CFT77 compilers  are
C       supplied with the standard SLAP distribution.
C
C *Precision:           Double Precision
C
C *Cautions:
C       This routine assumes that the diagonal of A is all  non-zero
C       and that the operation DINV = 1.0/DIAG(A) will not underflow
C       or overflow.    This  is done so that the  loop  vectorizes.
C       Matricies with zero or near zero or very  large entries will
C       have numerical difficulties  and  must  be fixed before this 
C       routine is called.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DSDS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
      DOUBLE PRECISION A(NELT), DINV(N)
C
C         Assume the Diagonal elements are the first in each column.
C         This loop should *VECTORIZE*.  If it does not you may have
C         to add a compiler directive.  We do not check for a zero
C         (or near zero) diagonal element since this would interfere 
C         with vectorization.  If this makes you nervous put a check
C         in!  It will run much slower.
C***FIRST EXECUTABLE STATEMENT  DSDS
 1    CONTINUE
      DO 10 ICOL = 1, N
         DINV(ICOL) = 1.0D0/A(JA(ICOL))
 10   CONTINUE
C         
      RETURN
C------------- LAST LINE OF DSDS FOLLOWS ----------------------------
      END
*DECK DSDSCL
      SUBROUTINE DSDSCL( N, NELT, IA, JA, A, ISYM, X, B, DINV, JOB,
     $     ITOL )
C***BEGIN PROLOGUE  DSDSCL
C***DATE WRITTEN   890404   (YYMMDD)
C***REVISION DATE  890404   (YYMMDD)
C***CATEGORY NO.  D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=DOUBLE PRECISION(DSDSCL-D),
C             SLAP Sparse, Diagonal
C***AUTHOR  Greenbaum, Anne, Courant Institute
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  Diagonal Scaling of system Ax = b.
C            This routine scales (and unscales) the system Ax = b 
C            by symmetric diagonal scaling.  The new system is:
C             -1/2  -1/2  1/2      -1/2
C            D    AD    (D   x) = D    b
C            when scaling is selected with the JOB parameter.  When
C            unscaling is selected this process is reversed.
C            The true solution is also scaled or unscaled if ITOL is set
C            appropriately, see below.
C***DESCRIPTION
C *Usage:
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, JOB, ITOL
C     DOUBLE PRECISION A(NELT), DINV(N)
C
C     CALL DSDSCL( N, NELT, IA, JA, A, ISYM, X, B, DINV, JOB, ITOL )
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C NELT   :IN       Integer.
C         Number of elements in arrays IA, JA, and A.
C IA     :IN       Integer IA(NELT).
C JA     :IN       Integer JA(NELT).
C A      :IN       Double Precision A(NELT).
C         These arrays should hold the matrix A in the SLAP Column
C         format.  See "Description", below. 
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all nonzero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the upper
C         or lower triangle of the matrix is stored.
C X      :INOUT    Double Precision X(N).
C         Initial guess that will be later used in the iterative 
C         solution.
C         of the scaled system.
C B      :INOUT    Double Precision B(N).
C         Right hand side vector.
C DINV   :OUT      Double Precision DINV(N).
C         Upon return this array holds 1./DIAG(A).
C JOB    :IN       Integer.
C         Flag indicating weather to scale or not.  JOB nonzero means
C         do scaling.  JOB = 0 means do unscaling.
C ITOL   :IN       Integer.
C         Flag indicating what type of error estimation to do in the 
C         iterative method.  When ITOL = 11 the exact solution from
C         common block solblk will be used.  When the system is scaled
C         then the true solution must also be scaled.  If ITOL is not
C         11 then this vector is not referenced.
C
C *Common Blocks:
C SOLN    :INOUT   Double Precision SOLN(N).  COMMON BLOCK /SOLBLK/
C         The true solution, SOLN, is scaled (or unscaled) if ITOL is
C         set to 11, see above.
C
C *Description
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
C       With the SLAP  format  all  of  the   "inner  loops" of this
C       routine should vectorize  on  machines with hardware support
C       for vector   gather/scatter  operations.  Your compiler  may
C       require a compiler directive to  convince it that  there are
C       no  implicit  vector  dependencies.  Compiler directives for
C       the Alliant    FX/Fortran and CRI   CFT/CFT77 compilers  are
C       supplied with the standard SLAP distribution.
C
C *Precision:           Double Precision
C
C *Cautions:
C       This routine assumes that the diagonal of A is all  non-zero
C       and that the operation DINV = 1.0/DIAG(A)  will  not  under-
C       flow or overflow. This is done so that the loop  vectorizes.
C       Matricies with zero or near zero or very  large entries will
C       have numerical difficulties  and  must  be fixed before this 
C       routine is called.
C
C *See Also:
C       DSDCG
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DSDSCL
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, JOB, ITOL
      DOUBLE PRECISION A(NELT), X(N), B(N), DINV(N)
      COMMON /SOLBLK/ SOLN(1)
C
C         SCALING...
C
      IF( JOB.NE.0 ) THEN
         DO 10 ICOL = 1, N
            DINV(ICOL) = 1.0D0/SQRT( A(JA(ICOL)) )
 10      CONTINUE
      ELSE
C
C         UNSCALING...
C
         DO 15 ICOL = 1, N
            DINV(ICOL) = 1.0D0/DINV(ICOL)
 15      CONTINUE
      ENDIF
C
      DO 30 ICOL = 1, N
         JBGN = JA(ICOL)
         JEND = JA(ICOL+1)-1
         DI = DINV(ICOL)
         DO 20 J = JBGN, JEND
            A(J) = DINV(IA(J))*A(J)*DI
 20      CONTINUE
 30   CONTINUE
C
      DO 40 ICOL = 1, N
         B(ICOL) = B(ICOL)*DINV(ICOL)
         X(ICOL) = X(ICOL)/DINV(ICOL)
 40   CONTINUE
C
C         Check to see if we need to scale the "true solution" as well.
C
      IF( ITOL.EQ.11 ) THEN
         DO 50 ICOL = 1, N
            SOLN(ICOL) = SOLN(ICOL)/DINV(ICOL)
 50      CONTINUE
      ENDIF
C
      RETURN
      END
*DECK DSD2S
      SUBROUTINE DSD2S(N, NELT, IA, JA, A, ISYM, DINV)
C***BEGIN PROLOGUE  DSD2S
C***DATE WRITTEN   890404   (YYMMDD)
C***REVISION DATE  890404   (YYMMDD)
C***CATEGORY NO.  D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=DOUBLE PRECISION(DSD2S-D),
C             SLAP Sparse, Diagonal
C***AUTHOR  Greenbaum, Anne, Courant Institute
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  Diagonal Scaling Preconditioner SLAP Normal Eqns Set Up.
C            Routine to compute the inverse of the diagonal of the 
C            matrix A*A'.  Where A is stored in SLAP-Column format.
C***DESCRIPTION
C *Usage:
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
C     DOUBLE PRECISION A(NELT), DINV(N)
C
C     CALL DSD2S( N, NELT, IA, JA, A, ISYM, DINV )
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C NELT   :IN       Integer.
C         Number of elements in arrays IA, JA, and A.
C IA     :IN       Integer IA(NELT).
C JA     :IN       Integer JA(NELT).
C A      :IN       Double Precision A(NELT).
C         These arrays should hold the matrix A in the SLAP Column
C         format.  See "Description", below. 
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all nonzero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the upper
C         or lower triangle of the matrix is stored.
C DINV   :OUT      Double Precision DINV(N).
C         Upon return this array holds 1./DIAG(A*A').
C
C *Description
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
C       With the SLAP  format  all  of  the   "inner  loops" of this
C       routine should vectorize  on  machines with hardware support
C       for vector   gather/scatter  operations.  Your compiler  may
C       require a compiler directive to  convince it that  there are
C       no  implicit  vector  dependencies.  Compiler directives for
C       the Alliant    FX/Fortran and CRI   CFT/CFT77 compilers  are
C       supplied with the standard SLAP distribution.
C
C *Precision:           Double Precision
C
C *Cautions:
C       This routine assumes that the diagonal of A is all  non-zero
C       and that the operation DINV = 1.0/DIAG(A*A') will not under-
C       flow or overflow. This is done so that the loop  vectorizes.
C       Matricies with zero or near zero or very  large entries will
C       have numerical difficulties  and  must  be fixed before this 
C       routine is called.
C
C *See Also:
C       DSDCGN
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DSD2S
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
      DOUBLE PRECISION A(NELT), DINV(N)
C
C***FIRST EXECUTABLE STATEMENT  DSD2S
      DO 10 I = 1, N
         DINV(I) = 0.
 10   CONTINUE
C
C         Loop over each column.
CVD$R NOCONCUR
      DO 40 I = 1, N
         KBGN = JA(I)
         KEND = JA(I+1) - 1
C
C         Add in the contributions for each row that has a non-zero
C         in this column.
CLLL. OPTION ASSERT (NOHAZARD)
CDIR$ IVDEP
CVD$ NODEPCHK
         DO 20 K = KBGN, KEND
            DINV(IA(K)) = DINV(IA(K)) + A(K)**2
 20      CONTINUE
         IF( ISYM.EQ.1 ) THEN
C
C         Lower triangle stored by columns => upper triangle stored by
C         rows with Diagonal being the first entry.  Loop across the
C         rest of the row.
            KBGN = KBGN + 1
            IF( KBGN.LE.KEND ) THEN
               DO 30 K = KBGN, KEND
                  DINV(I) = DINV(I) + A(K)**2
 30            CONTINUE
            ENDIF
         ENDIF
 40   CONTINUE
      DO 50 I=1,N
         DINV(I) = 1./DINV(I)
 50   CONTINUE
C         
      RETURN
C------------- LAST LINE OF DSD2S FOLLOWS ----------------------------
      END
*DECK DS2LT
      SUBROUTINE DS2LT( N, NELT, IA, JA, A, ISYM, NEL, IEL, JEL, EL )
C***BEGIN PROLOGUE  DS2LT
C***DATE WRITTEN   890404   (YYMMDD)
C***REVISION DATE  890404   (YYMMDD)
C***CATEGORY NO.  D2A4, D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=DOUBLE PRECISION(DS2LT-D),
C             Linear system, SLAP Sparse, Lower Triangle
C***AUTHOR  Greenbaum, Anne, Courant Institute
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  Lower Triangle Preconditioner SLAP Set Up.
C            Routine to store the lower triangle of a matrix stored
C            in the Slap Column format.
C***DESCRIPTION
C *Usage:
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
C     INTEGER NEL, IEL(N+1), JEL(NEL), NROW(N)
C     DOUBLE PRECISION A(NELT), EL(NEL)
C
C     CALL DS2LT( N, NELT, IA, JA, A, ISYM, NEL, IEL, JEL, EL )
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C NELT   :IN       Integer.
C         Number of non-zeros stored in A.
C IA     :IN       Integer IA(NELT).
C JA     :IN       Integer JA(NELT).
C A      :IN       Double Precision A(NELT).
C         These arrays should hold the matrix A in the SLAP Column
C         format.  See "Description", below. 
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all nonzero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the lower
C         triangle of the matrix is stored.
C NEL    :OUT      Integer.
C         Number of non-zeros in the lower triangle of A.   Also 
C         coresponds to the length of the JEL, EL arrays.
C IEL    :OUT      Integer IEL(N+1).
C JEL    :OUT      Integer JEL(NEL).
C EL     :OUT      Double Precision     EL(NEL).
C         IEL, JEL, EL contain the lower triangle of the A matrix
C         stored in SLAP Column format.  See "Description", below
C         for more details bout the SLAP Column format.
C
C *Description
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
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DS2LT
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM 
      INTEGER NEL, IEL(NEL), JEL(NEL)
      DOUBLE PRECISION A(NELT), EL(NELT)
C***FIRST EXECUTABLE STATEMENT  DS2LT
      IF( ISYM.EQ.0 ) THEN
C
C         The matrix is stored non-symmetricly.  Pick out the lower
C         triangle.
C
         NEL = 0
         DO 20 ICOL = 1, N
            JEL(ICOL) = NEL+1
            JBGN = JA(ICOL)
            JEND = JA(ICOL+1)-1
CVD$ NOVECTOR
            DO 10 J = JBGN, JEND
               IF( IA(J).GE.ICOL ) THEN
                  NEL = NEL + 1
                  IEL(NEL) = IA(J)
                  EL(NEL)  = A(J)
               ENDIF
 10         CONTINUE
 20      CONTINUE
         JEL(N+1) = NEL+1
      ELSE
C
C         The matrix is symmetric and only the lower triangle is 
C         stored.  Copy it to IEL, JEL, EL.
C
         NEL = NELT
         DO 30 I = 1, NELT
            IEL(I) = IA(I)
            EL(I) = A(I)
 30      CONTINUE
         DO 40 I = 1, N+1
            JEL(I) = JA(I)
 40      CONTINUE
      ENDIF
      RETURN
C------------- LAST LINE OF DS2LT FOLLOWS ----------------------------
      END
*DECK DSICS
      SUBROUTINE DSICS(N, NELT, IA, JA, A, ISYM, NEL, IEL, JEL,
     $     EL, D, R, IWARN )
C***BEGIN PROLOGUE  DSICS
C***DATE WRITTEN   890404   (YYMMDD)
C***REVISION DATE  890404   (YYMMDD)
C***CATEGORY NO.  D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=DOUBLE PRECISION(DSICS-D),
C             Linear system, SLAP Sparse, Iterative Precondition
C             Incomplete Cholesky Factorization.
C***AUTHOR  Greenbaum, Anne, Courant Institute
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  Incompl Cholesky Decomposition Preconditioner SLAP Set Up.
C            Routine to generate the Incomplete Cholesky decomposition,
C            L*D*L-trans, of  a symmetric positive definite  matrix, A,
C            which  is stored  in  SLAP Column format.  The  unit lower
C            triangular matrix L is  stored by rows, and the inverse of
C            the diagonal matrix D is stored.
C***DESCRIPTION
C *Usage:
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
C     INTEGER NEL, IEL(NEL), JEL(N+1), IWARN
C     DOUBLE PRECISION A(NELT), EL(NEL), D(N), R(N)
C
C     CALL DSICS( N, NELT, IA, JA, A, ISYM, NEL, IEL, JEL, EL, D, R,
C    $    IWARN )
C
C *Arguments:
C N      :IN       Integer.
C         Order of the Matrix.
C NELT   :IN       Integer.
C         Number of elements in arrays IA, JA, and A.
C IA     :INOUT    Integer IA(NELT).
C JA     :INOUT    Integer JA(NELT).
C A      :INOUT    Double Precision A(NELT).
C         These arrays should hold the matrix A in the SLAP Column
C         format.  See "Description", below. 
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all nonzero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the lower 
C         triangle of the matrix is stored.
C NEL    :OUT      Integer.
C         Number of non-zeros in the lower triangle of A.   Also 
C         coresponds to the length of the JEL, EL arrays.
C IEL    :OUT      Integer IEL(N+1).
C JEL    :OUT      Integer JEL(NEL).
C EL     :OUT      Double Precision EL(NEL).
C         IEL, JEL, EL contain the unit lower triangular factor  of the
C         incomplete decomposition   of the A  matrix  stored  in  SLAP
C         Row format.   The Diagonal of   ones   *IS*   stored.     See
C         "Description", below for more details about the SLAP Row fmt.
C D      :OUT      Double Precision D(N)
C         Upon return this array holds D(I) = 1./DIAG(A).
C R      :WORK     Double Precision R(N).
C         Temporary double precision workspace needed for the 
C         factorization.
C IWARN  :OUT      Integer.
C         This is a warning variable and is zero if the IC factoriza-
C         tion goes well.  It is set to the row index corresponding to
C         the last zero pivot found.  See "Description", below.
C
C *Description
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
C       ==================== S L A P Row format ====================
C       This routine requires  that the matrix A  be  stored  in the
C       SLAP  Row format.   In this format  the non-zeros are stored
C       counting across  rows (except for the diagonal  entry, which
C       must appear first in each "row") and  are stored in the 
C       double precision
C       array A.  In other words, for each row in the matrix put the
C       diagonal entry in  A.   Then   put  in the   other  non-zero
C       elements   going  across the  row (except   the diagonal) in
C       order.   The  JA array  holds   the column   index for  each
C       non-zero.   The IA  array holds the  offsets into  the JA, A
C       arrays  for   the   beginning  of   each  row.   That    is,
C       JA(IA(IROW)),  A(IA(IROW)) points  to  the beginning  of the
C       IROW-th row in JA and A.   JA(IA(IROW+1)-1), A(IA(IROW+1)-1)
C       points to the  end of the  IROW-th row.  Note that we always
C       have IA(N+1) =  NELT+1, where  N  is  the number of rows  in
C       the matrix  and NELT  is the  number   of  non-zeros in  the
C       matrix.
C       
C       Here is an example of the SLAP Row storage format for a  5x5
C       Matrix (in the A and JA arrays '|' denotes the end of a row):
C
C           5x5 Matrix         SLAP Row format for 5x5 matrix on left.
C                              1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 12 15 | 22 21 | 33 35 | 44 | 55 51 53
C       |21 22  0  0  0|  JA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  IA:  1  4  6    8  9   12
C       | 0  0  0 44  0|  
C       |51  0 53  0 55|  
C
C       With the SLAP  format some  of  the   "inner  loops" of this
C       routine should vectorize  on  machines with hardware support
C       for vector   gather/scatter  operations.  Your compiler  may
C       require a compiler directive to  convince it that  there are
C       no  implicit  vector  dependencies.  Compiler directives for
C       the Alliant    FX/Fortran and CRI   CFT/CFT77 compilers  are
C       supplied with the standard SLAP distribution.
C
C       The IC  factorization is not  alway exist for SPD matricies.
C       In the event that a zero pivot is found it is set  to be 1.0
C       and the factorization procedes.   The integer variable IWARN
C       is set to the last row where the Diagonal was fudged.  This 
C       eventuality hardly ever occurs in practice
C       
C *Precision:           Double Precision
C
C *See Also:
C       SCG, DSICCG
C***REFERENCES  1. Gene Golub & Charles Van Loan, "Matrix Computations",
C                 John Hopkins University Press; 3 (1983) IBSN 
C                 0-8018-3010-9.
C***ROUTINES CALLED  XERRWV
C***END PROLOGUE  DSICS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
      INTEGER NEL, IEL(NEL), JEL(NEL)
      DOUBLE PRECISION A(NELT), EL(NEL), D(N), R(N)
C         
C         Set the lower triangle in IEL, JEL, EL
C***FIRST EXECUTABLE STATEMENT  DSICS
      IWARN = 0
C         
C         All matrix elements stored in IA, JA, A.  Pick out the lower 
C         triangle (making sure that the Diagonal of EL is one) and
C         store by rows.
C         
      NEL = 1
      IEL(1) = 1
      JEL(1) = 1
      EL(1) = 1.0D0
      D(1) = A(1)
CVD$R NOCONCUR
      DO 30 IROW = 2, N
C         Put in the Diagonal.
         NEL = NEL + 1
         IEL(IROW) = NEL
         JEL(NEL) = IROW
         EL(NEL) = 1.0D0
         D(IROW) = A(JA(IROW))
C
C         Look in all the lower triangle columns for a matching row.
C         Since the matrix is symmetric, we can look across the 
C         irow-th row by looking down the irow-th column (if it is 
C         stored ISYM=0)...
         IF( ISYM.EQ.0 ) THEN
            ICBGN = JA(IROW)
            ICEND = JA(IROW+1)-1
         ELSE
            ICBGN = 1
            ICEND = IROW-1
         ENDIF
         DO 20 IC = ICBGN, ICEND
            IF( ISYM.EQ.0 ) THEN
               ICOL = IA(IC)
               IF( ICOL.GE.IROW ) GOTO 20
            ELSE
               ICOL = IC
            ENDIF
            JBGN = JA(ICOL)+1
            JEND = JA(ICOL+1)-1
            IF( JBGN.LE.JEND .AND. IA(JEND).GE.IROW ) THEN
CVD$ NOVECTOR
               DO 10 J = JBGN, JEND
                  IF( IA(J).EQ.IROW ) THEN
                     NEL = NEL + 1
                     JEL(NEL) = ICOL
                     EL(NEL)  = A(J)
                     GOTO 20
                  ENDIF
 10            CONTINUE
            ENDIF
 20      CONTINUE
 30   CONTINUE
      IEL(N+1) = NEL+1
C
C         Sort ROWS of lower triangle into descending order (count out 
C         along rows out from Diagonal).
C
      DO 60 IROW = 2, N
         IBGN = IEL(IROW)+1
         IEND = IEL(IROW+1)-1
         IF( IBGN.LT.IEND ) THEN
            DO 50 I = IBGN, IEND-1
CVD$ NOVECTOR
               DO 40 J = I+1, IEND
                  IF( JEL(I).GT.JEL(J) ) THEN
                     JELTMP = JEL(J)
                     JEL(J) = JEL(I)
                     JEL(I) = JELTMP
                     ELTMP = EL(J)
                     EL(J) = EL(I)
                     EL(I) = ELTMP
                  ENDIF
 40            CONTINUE
 50         CONTINUE
         ENDIF
 60   CONTINUE
C         
C         Perform the Incomplete Cholesky decomposition by looping 
C         over the rows.
C         Scale the first column.  Use the structure of A to pick out
C         the rows with something in column 1.
C
      IRBGN = JA(1)+1
      IREND = JA(2)-1
      DO 65 IRR = IRBGN, IREND
         IR = IA(IRR)
C         Find the index into EL for EL(1,IR).  
C         Hint: it's the second entry.
         I = IEL(IR)+1
         EL(I) = EL(I)/D(1)
 65   CONTINUE
C
      DO 110 IROW = 2, N
C
C         Update the IROW-th diagonal.
C
         DO 66 I = 1, IROW-1
            R(I) = 0.0D0
 66      CONTINUE
         IBGN = IEL(IROW)+1
         IEND = IEL(IROW+1)-1
         IF( IBGN.LE.IEND ) THEN
CLLL. OPTION ASSERT (NOHAZARD)
CDIR$ IVDEP
CVD$ NODEPCHK
            DO 70 I = IBGN, IEND
               R(JEL(I)) = EL(I)*D(JEL(I))
               D(IROW) = D(IROW) - EL(I)*R(JEL(I))
 70         CONTINUE
C
C         Check to see if we gota problem with the diagonal.
C
            IF( D(IROW).LE.0.0D0 ) THEN
               IF( IWARN.EQ.0 ) IWARN = IROW
               D(IROW) = 1.0D0
            ENDIF
         ENDIF
C
C         Update each EL(IROW+1:N,IROW), if there are any.
C         Use the structure of A to determine the Non-zero elements
C         of the IROW-th column of EL.
C
         IRBGN = JA(IROW)
         IREND = JA(IROW+1)-1
         DO 100 IRR = IRBGN, IREND
            IR = IA(IRR)
            IF( IR.LE.IROW ) GOTO 100
C         Find the index into EL for EL(IR,IROW)
            IBGN = IEL(IR)+1
            IEND = IEL(IR+1)-1
            IF( JEL(IBGN).GT.IROW ) GOTO 100
            DO 90 I = IBGN, IEND
               IF( JEL(I).EQ.IROW ) THEN
                  ICEND = IEND
 91               IF( JEL(ICEND).GE.IROW ) THEN
                     ICEND = ICEND - 1
                     GOTO 91
                  ENDIF
C         Sum up the EL(IR,1:IROW-1)*R(1:IROW-1) contributions.
CLLL. OPTION ASSERT (NOHAZARD)
CDIR$ IVDEP
CVD$ NODEPCHK
                  DO 80 IC = IBGN, ICEND
                     EL(I) = EL(I) - EL(IC)*R(JEL(IC))
 80               CONTINUE
                  EL(I) = EL(I)/D(IROW)
                  GOTO 100
               ENDIF
 90         CONTINUE
C
C         If we get here, we have real problems...

C WHL: Usually this means the matrix A is not symmetric.
C      The preconditioner is fragile in the sense that it can fail with
C       very small departures from symmetry (due to roundoff errors).
 
            CALL XERRWV('DSICS -- A and EL data structure mismatch'//
     $           ' in row (i1)',53,1,2,1,IROW,0,0,0.0,0.0)
 100     CONTINUE
 110  CONTINUE
C        
C         Replace diagonals by their inverses.
C
CVD$ CONCUR
      DO 120 I =1, N
         D(I) = 1.0D0/D(I)
 120  CONTINUE
      RETURN
C------------- LAST LINE OF DSICS FOLLOWS ----------------------------
      END
*DECK DSILUS
      SUBROUTINE DSILUS(N, NELT, IA, JA, A, ISYM, NL, IL, JL,
     $     L, DINV, NU, IU, JU, U, NROW, NCOL)
C***BEGIN PROLOGUE  DSILUS
C***DATE WRITTEN   890404   (YYMMDD)
C***REVISION DATE  890404   (YYMMDD)
C***CATEGORY NO.  D2A4, D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=DOUBLE PRECISION(DSILUS-D),
C             Non-Symmetric Linear system, Sparse, 
C             Iterative Precondition, Incomplete LU Factorization
C***AUTHOR  Greenbaum, Anne, Courant Institute
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  Incomplete LU Decomposition Preconditioner SLAP Set Up.
C            Routine to generate the incomplete LDU decomposition of a 
C            matrix.  The  unit lower triangular factor L is stored by 
C            rows and the  unit upper triangular factor U is stored by 
C            columns.  The inverse of the diagonal matrix D is stored.
C            No fill in is allowed.
C***DESCRIPTION
C *Usage:
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
C     INTEGER NL, IL(N+1), JL(NL), NU, IU(N+1), JU(NU)
C     INTEGER NROW(N), NCOL(N)
C     DOUBLE PRECISION A(NELT), L(NL), U(NU), DINV(N)
C
C     CALL DSILUS( N, NELT, IA, JA, A, ISYM, NL, IL, JL, L, 
C    $    DINV, NU, IU, JU, U, NROW, NCOL )
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C NELT   :IN       Integer.
C         Number of elements in arrays IA, JA, and A.
C IA     :IN       Integer IA(NELT).
C JA     :IN       Integer JA(NELT).
C A      :IN       Double Precision A(NELT).
C         These arrays should hold the matrix A in the SLAP Column
C         format.  See "Description", below. 
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all nonzero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the lower 
C         triangle of the matrix is stored.
C NL     :OUT      Integer.
C         Number of non-zeros in the EL array.
C IL     :OUT      Integer IL(N+1).
C JL     :OUT      Integer JL(NL).
C L      :OUT      Double Precision L(NL).
C         IL, JL, L  contain the unit ower  triangular factor of  the
C         incomplete decomposition  of some  matrix stored  in   SLAP
C         Row format.     The   Diagonal  of ones  *IS*  stored.  See
C         "DESCRIPTION", below for more details about the SLAP format.
C NU     :OUT      Integer.
C         Number of non-zeros in the U array.     
C IU     :OUT      Integer IU(N+1).
C JU     :OUT      Integer JU(NU).
C U      :OUT      Double Precision     U(NU).
C         IU, JU, U contain   the unit upper triangular factor of the
C         incomplete  decomposition    of some matrix  stored in SLAP
C         Column  format.   The Diagonal of ones   *IS*  stored.  See 
C         "Description", below  for  more  details  about  the   SLAP 
C         format.
C NROW   :WORK     Integer NROW(N).
C         NROW(I) is the number of non-zero elements in the I-th row
C         of L.
C NCOL   :WORK     Integer NCOL(N).
C         NCOL(I) is the number of non-zero elements in the I-th 
C         column of U.
C
C *Description
C       IL, JL, L should contain the unit  lower triangular factor of
C       the incomplete decomposition of the A matrix  stored in SLAP
C       Row format.  IU, JU, U should contain  the unit upper factor
C       of the  incomplete decomposition of  the A matrix  stored in
C       SLAP Column format This ILU factorization can be computed by
C       the DSILUS routine.  The diagonals (which is all one's) are
C       stored.
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
C       ==================== S L A P Row format ====================
C       This routine requires  that the matrix A  be  stored  in the
C       SLAP  Row format.   In this format  the non-zeros are stored
C       counting across  rows (except for the diagonal  entry, which
C       must appear first in each "row") and  are stored in the 
C       double precision
C       array A.  In other words, for each row in the matrix put the
C       diagonal entry in  A.   Then   put  in the   other  non-zero
C       elements   going  across the  row (except   the diagonal) in
C       order.   The  JA array  holds   the column   index for  each
C       non-zero.   The IA  array holds the  offsets into  the JA, A
C       arrays  for   the   beginning  of   each  row.   That    is,
C       JA(IA(IROW)),  A(IA(IROW)) points  to  the beginning  of the
C       IROW-th row in JA and A.   JA(IA(IROW+1)-1), A(IA(IROW+1)-1)
C       points to the  end of the  IROW-th row.  Note that we always
C       have IA(N+1) =  NELT+1, where  N  is  the number of rows  in
C       the matrix  and NELT  is the  number   of  non-zeros in  the
C       matrix.
C       
C       Here is an example of the SLAP Row storage format for a  5x5
C       Matrix (in the A and JA arrays '|' denotes the end of a row):
C
C           5x5 Matrix         SLAP Row format for 5x5 matrix on left.
C                              1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 12 15 | 22 21 | 33 35 | 44 | 55 51 53
C       |21 22  0  0  0|  JA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  IA:  1  4  6    8  9   12
C       | 0  0  0 44  0|  
C       |51  0 53  0 55|  
C
C *Precision:           Double Precision
C *See Also:
C       SILUR
C***REFERENCES  1. Gene Golub & Charles Van Loan, "Matrix Computations",
C                 John Hopkins University Press; 3 (1983) IBSN 
C                 0-8018-3010-9.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DSILUS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, NL, IL(NL), JL(NL)
      INTEGER NU, IU(NU), JU(NU), NROW(N), NCOL(N)
      DOUBLE PRECISION A(NELT), L(NL), DINV(N), U(NU)
C         
C         Count number of elements in each row of the lower triangle.
C***FIRST EXECUTABLE STATEMENT  DSILUS
      DO 10 I=1,N
         NROW(I) = 0
         NCOL(I) = 0
 10   CONTINUE
CVD$R NOCONCUR
CVD$R NOVECTOR
      DO 30 ICOL = 1, N
         JBGN = JA(ICOL)+1
         JEND = JA(ICOL+1)-1
         IF( JBGN.LE.JEND ) THEN
            DO 20 J = JBGN, JEND
               IF( IA(J).LT.ICOL ) THEN
                  NCOL(ICOL) = NCOL(ICOL) + 1
               ELSE
                  NROW(IA(J)) = NROW(IA(J)) + 1
                  IF( ISYM.NE.0 ) NCOL(IA(J)) = NCOL(IA(J)) + 1
               ENDIF
 20         CONTINUE
         ENDIF
 30   CONTINUE
      JU(1) = 1
      IL(1) = 1
      DO 40 ICOL = 1, N
         IL(ICOL+1) = IL(ICOL) + NROW(ICOL)
         JU(ICOL+1) = JU(ICOL) + NCOL(ICOL)
         NROW(ICOL) = IL(ICOL)
         NCOL(ICOL) = JU(ICOL)
 40   CONTINUE
C         
C         Copy the matrix A into the L and U structures.
      DO 60 ICOL = 1, N
         DINV(ICOL) = A(JA(ICOL))
         JBGN = JA(ICOL)+1
         JEND = JA(ICOL+1)-1
         IF( JBGN.LE.JEND ) THEN
            DO 50 J = JBGN, JEND
               IROW = IA(J)
               IF( IROW.LT.ICOL ) THEN
C         Part of the upper triangle.
                  IU(NCOL(ICOL)) = IROW
                  U(NCOL(ICOL)) = A(J)
                  NCOL(ICOL) = NCOL(ICOL) + 1
               ELSE
C         Part of the lower triangle (stored by row).
                  JL(NROW(IROW)) = ICOL
                  L(NROW(IROW)) = A(J)
                  NROW(IROW) = NROW(IROW) + 1
                  IF( ISYM.NE.0 ) THEN
C         Symmetric...Copy lower triangle into upper triangle as well.
                     IU(NCOL(IROW)) = ICOL
                     U(NCOL(IROW)) = A(J)
                     NCOL(IROW) = NCOL(IROW) + 1
                  ENDIF
               ENDIF
 50         CONTINUE
         ENDIF
 60   CONTINUE
C
C         Sort the rows of L and the columns of U.
      DO 110 K = 2, N
         JBGN = JU(K)
         JEND = JU(K+1)-1
         IF( JBGN.LT.JEND ) THEN
            DO 80 J = JBGN, JEND-1
               DO 70 I = J+1, JEND
                  IF( IU(J).GT.IU(I) ) THEN
                     ITEMP = IU(J)
                     IU(J) = IU(I)
                     IU(I) = ITEMP
                     TEMP = U(J)
                     U(J) = U(I)
                     U(I) = TEMP
                  ENDIF
 70            CONTINUE
 80         CONTINUE
         ENDIF
         IBGN = IL(K)
         IEND = IL(K+1)-1
         IF( IBGN.LT.IEND ) THEN
            DO 100 I = IBGN, IEND-1
               DO 90 J = I+1, IEND
                  IF( JL(I).GT.JL(J) ) THEN
                     JTEMP = JU(I)
                     JU(I) = JU(J)
                     JU(J) = JTEMP
                     TEMP = L(I)
                     L(I) = L(J)
                     L(J) = TEMP
                  ENDIF
 90            CONTINUE
 100        CONTINUE
         ENDIF
 110  CONTINUE
C
C         Perform the incomplete LDU decomposition.
      DO 300 I=2,N
C         
C           I-th row of L
         INDX1 = IL(I)
         INDX2 = IL(I+1) - 1
         IF(INDX1 .GT. INDX2) GO TO 200
         DO 190 INDX=INDX1,INDX2
            IF(INDX .EQ. INDX1) GO TO 180
            INDXR1 = INDX1
            INDXR2 = INDX - 1
            INDXC1 = JU(JL(INDX))
            INDXC2 = JU(JL(INDX)+1) - 1
            IF(INDXC1 .GT. INDXC2) GO TO 180
 160        KR = JL(INDXR1)
 170        KC = IU(INDXC1)
            IF(KR .GT. KC) THEN
               INDXC1 = INDXC1 + 1
               IF(INDXC1 .LE. INDXC2) GO TO 170
            ELSEIF(KR .LT. KC) THEN
               INDXR1 = INDXR1 + 1
               IF(INDXR1 .LE. INDXR2) GO TO 160
            ELSEIF(KR .EQ. KC) THEN
               L(INDX) = L(INDX) - L(INDXR1)*DINV(KC)*U(INDXC1)
               INDXR1 = INDXR1 + 1
               INDXC1 = INDXC1 + 1
               IF(INDXR1 .LE. INDXR2 .AND. INDXC1 .LE. INDXC2) GO TO 160
            ENDIF
 180        L(INDX) = L(INDX)/DINV(JL(INDX))
 190     CONTINUE
C         
C         ith column of u
 200     INDX1 = JU(I)
         INDX2 = JU(I+1) - 1
         IF(INDX1 .GT. INDX2) GO TO 260
         DO 250 INDX=INDX1,INDX2
            IF(INDX .EQ. INDX1) GO TO 240
            INDXC1 = INDX1
            INDXC2 = INDX - 1
            INDXR1 = IL(IU(INDX))
            INDXR2 = IL(IU(INDX)+1) - 1
            IF(INDXR1 .GT. INDXR2) GO TO 240
 210        KR = JL(INDXR1)
 220        KC = IU(INDXC1)
            IF(KR .GT. KC) THEN
               INDXC1 = INDXC1 + 1
               IF(INDXC1 .LE. INDXC2) GO TO 220
            ELSEIF(KR .LT. KC) THEN
               INDXR1 = INDXR1 + 1
               IF(INDXR1 .LE. INDXR2) GO TO 210
            ELSEIF(KR .EQ. KC) THEN
               U(INDX) = U(INDX) - L(INDXR1)*DINV(KC)*U(INDXC1)
               INDXR1 = INDXR1 + 1
               INDXC1 = INDXC1 + 1
               IF(INDXR1 .LE. INDXR2 .AND. INDXC1 .LE. INDXC2) GO TO 210
            ENDIF
 240        U(INDX) = U(INDX)/DINV(IU(INDX))
 250     CONTINUE
C         
C         ith diagonal element
 260     INDXR1 = IL(I)
         INDXR2 = IL(I+1) - 1
         IF(INDXR1 .GT. INDXR2) GO TO 300
         INDXC1 = JU(I)
         INDXC2 = JU(I+1) - 1
         IF(INDXC1 .GT. INDXC2) GO TO 300
 270     KR = JL(INDXR1)
 280     KC = IU(INDXC1)
         IF(KR .GT. KC) THEN
            INDXC1 = INDXC1 + 1
            IF(INDXC1 .LE. INDXC2) GO TO 280
         ELSEIF(KR .LT. KC) THEN
            INDXR1 = INDXR1 + 1
            IF(INDXR1 .LE. INDXR2) GO TO 270
         ELSEIF(KR .EQ. KC) THEN
            DINV(I) = DINV(I) - L(INDXR1)*DINV(KC)*U(INDXC1)
            INDXR1 = INDXR1 + 1
            INDXC1 = INDXC1 + 1
            IF(INDXR1 .LE. INDXR2 .AND. INDXC1 .LE. INDXC2) GO TO 270
         ENDIF
C         
 300  CONTINUE
C         
C         replace diagonal lts by their inverses.
CVD$ VECTOR
      DO 430 I=1,N
         DINV(I) = 1./DINV(I)
 430  CONTINUE
C         
      RETURN
C------------- LAST LINE OF DSILUS FOLLOWS ----------------------------
      END
