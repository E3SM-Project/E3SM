module e4d_mat_inv_module
  
  public :: MIGS, ELGS

contains

! Updated 10/24/2001.
!
!cccccccccccccccccccccccc     Program 4.4     cccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                      c
! Please Note:                                                         c
!                                                                      c
! (1) This computer program is part of the book, "An Introduction to   c
!     Computational Physics," written by Tao Pang and published and    c
!     copyrighted by Cambridge University Press in 1997.               c
!                                                                      c
! (2) No warranties, express or implied, are made for this program.    c
!                                                                      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      SUBROUTINE MIGS(A,N,X,INDX)
!
! Subroutine to invert matrix A(N,N) with the inverse stored
! in X(N,N) in the output.

!
      DIMENSION A(N,N),X(N,N),INDX(N),B(N,N)
!

      DO I = 1, N
        DO J = 1, N
          B(I,J) = 0.0
        ENDDO
      ENDDO
      DO I = 1, N
        B(I,I) = 1.0
      ENDDO
!
      CALL ELGS(A,N,INDX)
!
      DO I = 1, N-1
        DO J = I+1, N
          DO K = 1, N
            B(INDX(J),K) = B(INDX(J),K) &
                          -A(INDX(J),I)*B(INDX(I),K)
          ENDDO
        ENDDO
      ENDDO
!
      DO I = 1, N
        X(N,I) = B(INDX(N),I)/A(INDX(N),N)
        DO J = N-1, 1, -1
          X(J,I) = B(INDX(J),I)
          DO K = J+1, N
            X(J,I) = X(J,I)-A(INDX(J),K)*X(K,I)
          ENDDO
          X(J,I) =  X(J,I)/A(INDX(J),J)
        ENDDO
      ENDDO
!
      RETURN
      END SUBROUTINE
!
      SUBROUTINE ELGS(A,N,INDX)
!
! Subroutine to perform the partial-pivoting Gaussian elimination.
! A(N,N) is the original matrix in the input and transformed
! matrix plus the pivoting element ratios below the diagonal in
! the output.  INDX(N) records the pivoting order.

!
      DIMENSION A(N,N),INDX(N),C(N)
!
! Initialize the index
!
      DO I = 1, N
        INDX(I) = I
      ENDDO
!
! Find the rescaling factors, one from each row
!
        DO I = 1, N
          C1= 0.0
          DO J = 1, N
            C1 = AMAX1(C1,ABS(A(I,J)))
          ENDDO
          C(I) = C1
        ENDDO
!
! Search the pivoting (largest) element from each column
!
      DO J = 1, N-1
        PI1 = 0.0
        DO I = J, N
          PI = ABS(A(INDX(I),J))/C(INDX(I))
          IF (PI.GT.PI1) THEN
            PI1 = PI
            K   = I
          ELSE
          ENDIF
        ENDDO
!
! Interchange the rows via INDX(N) to record pivoting order
!
        ITMP    = INDX(J)
        INDX(J) = INDX(K)
        INDX(K) = ITMP
        DO I = J+1, N
          PJ  = A(INDX(I),J)/A(INDX(J),J)
!
! Record pivoting ratios below the diagonal
!
          A(INDX(I),J) = PJ
!
! Modify other elements accordingly
!
          DO K = J+1, N
            A(INDX(I),K) = A(INDX(I),K)-PJ*A(INDX(J),K)
          ENDDO
        ENDDO
      ENDDO
!
      RETURN
      END SUBROUTINE

end module e4d_mat_inv_module
