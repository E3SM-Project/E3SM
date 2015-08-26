Module m_mrgrnk
Integer, Parameter :: kdp = selected_real_kind(15)
public :: mrgrnk
private :: kdp
private :: R_mrgrnk, I_mrgrnk, D_mrgrnk
interface mrgrnk
  module procedure D_mrgrnk, R_mrgrnk, I_mrgrnk
end interface mrgrnk
contains

Subroutine D_mrgrnk (XDONT, IRNGT)
! __________________________________________________________
!   MRGRNK = Merge-sort ranking of an array
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
! __________________________________________________________
! __________________________________________________________
      Real (kind=kdp), Dimension (:), Intent (In) :: XDONT
      Integer*4, Dimension (:), Intent (Out) :: IRNGT
! __________________________________________________________
      Real (kind=kdp) :: XVALA, XVALB
!
      Integer, Dimension (SIZE(IRNGT)) :: JWRKT
      Integer :: LMTNA, LMTNC, IRNG1, IRNG2
      Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
!
      NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
      Select Case (NVAL)
      Case (:0)
         Return
      Case (1)
         IRNGT (1) = 1
         Return
      Case Default
         Continue
      End Select
!
!  Fill-in the index array, creating ordered couples
!
      Do IIND = 2, NVAL, 2
         If (XDONT(IIND-1) <= XDONT(IIND)) Then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         Else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         End If
      End Do
      If (Modulo(NVAL, 2) /= 0) Then
         IRNGT (NVAL) = NVAL
      End If
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      LMTNA = 2
      LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      Do
         If (NVAL <= 2) Exit
!
!   Loop on merges of A and B into C
!
         Do IWRKD = 0, NVAL - 1, 4
            If ((IWRKD+4) > NVAL) Then
               If ((IWRKD+2) >= NVAL) Exit
!
!   1 2 3
!
               If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
!
!   1 3 2
!
               If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
!
!   3 1 2
!
               Else
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               End If
               Exit
            End If
!
!   1 2 3 4
!
            If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
!
!   1 3 x x
!
            If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   1 3 2 4
                  IRNGT (IWRKD+3) = IRNG2
               Else
!   1 3 4 2
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               End If
!
!   3 x x x
!
            Else
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                  IRNGT (IWRKD+2) = IRNG1
                  If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   3 1 2 4
                     IRNGT (IWRKD+3) = IRNG2
                  Else
!   3 1 4 2
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  End If
               Else
!   3 4 1 2
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               End If
            End If
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 4
         Exit
      End Do
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is doubled.
!
      Do
         If (LMTNA >= NVAL) Exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
!
!   Loop on merges of A and B into C
!
         Do
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JINDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            If (IWRKF >= NVAL) Then
               If (JINDA >= NVAL) Exit
               IWRKF = NVAL
            End If
            IINDA = 1
            IINDB = JINDA + 1
!
!   Shortcut for the case when the max of A is smaller
!   than the min of B. This line may be activated when the
!   initial set is already close to sorted.
!
!          IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
!
!  One steps in the C subset, that we build in the final rank array
!
!  Make a copy of the rank array for the merge iteration
!
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
!
            XVALA = XDONT (JWRKT(IINDA))
            XVALB = XDONT (IRNGT(IINDB))
!
            Do
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               If (XVALA > XVALB) Then
                  IRNGT (IWRK) = IRNGT (IINDB)
                  IINDB = IINDB + 1
                  If (IINDB > IWRKF) Then
!  Only A still with unprocessed values
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                     Exit
                  End If
                  XVALB = XDONT (IRNGT(IINDB))
               Else
                  IRNGT (IWRK) = JWRKT (IINDA)
                  IINDA = IINDA + 1
                  If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                  XVALA = XDONT (JWRKT(IINDA))
               End If
!
            End Do
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      End Do
!
      Return
!
End Subroutine D_mrgrnk

Subroutine R_mrgrnk (XDONT, IRNGT)
! __________________________________________________________
!   MRGRNK = Merge-sort ranking of an array
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
! __________________________________________________________
! _________________________________________________________
      Real, Dimension (:), Intent (In) :: XDONT
      Integer*4, Dimension (:), Intent (Out) :: IRNGT
! __________________________________________________________
      Real :: XVALA, XVALB
!
      Integer, Dimension (SIZE(IRNGT)) :: JWRKT
      Integer :: LMTNA, LMTNC, IRNG1, IRNG2
      Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
!
      NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
      Select Case (NVAL)
      Case (:0)
         Return
      Case (1)
         IRNGT (1) = 1
         Return
      Case Default
         Continue
      End Select
!
!  Fill-in the index array, creating ordered couples
!
      Do IIND = 2, NVAL, 2
         If (XDONT(IIND-1) <= XDONT(IIND)) Then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         Else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         End If
      End Do
      If (Modulo(NVAL, 2) /= 0) Then
         IRNGT (NVAL) = NVAL
      End If
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      LMTNA = 2
      LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      Do
         If (NVAL <= 2) Exit
!
!   Loop on merges of A and B into C
!
         Do IWRKD = 0, NVAL - 1, 4
            If ((IWRKD+4) > NVAL) Then
               If ((IWRKD+2) >= NVAL) Exit
!
!   1 2 3
!
               If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
!
!   1 3 2
!
               If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
!
!   3 1 2
!
               Else
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               End If
               Exit
            End If
!
!   1 2 3 4
!
            If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
!
!   1 3 x x
!
            If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   1 3 2 4
                  IRNGT (IWRKD+3) = IRNG2
               Else
!   1 3 4 2
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               End If
!
!   3 x x x
!
            Else
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                  IRNGT (IWRKD+2) = IRNG1
                  If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   3 1 2 4
                     IRNGT (IWRKD+3) = IRNG2
                  Else
!   3 1 4 2
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  End If
               Else
!   3 4 1 2
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               End If
            End If
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 4
         Exit
      End Do
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is doubled.
!
      Do
         If (LMTNA >= NVAL) Exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
!
!   Loop on merges of A and B into C
!
         Do
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JINDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            If (IWRKF >= NVAL) Then
               If (JINDA >= NVAL) Exit
               IWRKF = NVAL
            End If
            IINDA = 1
            IINDB = JINDA + 1
!
!   Shortcut for the case when the max of A is smaller
!   than the min of B. This line may be activated when the
!   initial set is already close to sorted.
!
!          IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
!
!  One steps in the C subset, that we build in the final rank array
!
!  Make a copy of the rank array for the merge iteration
!
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
!
            XVALA = XDONT (JWRKT(IINDA))
            XVALB = XDONT (IRNGT(IINDB))
!
            Do
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               If (XVALA > XVALB) Then
                  IRNGT (IWRK) = IRNGT (IINDB)
                  IINDB = IINDB + 1
                  If (IINDB > IWRKF) Then
!  Only A still with unprocessed values
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                     Exit
                  End If
                  XVALB = XDONT (IRNGT(IINDB))
               Else
                  IRNGT (IWRK) = JWRKT (IINDA)
                  IINDA = IINDA + 1
                  If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                  XVALA = XDONT (JWRKT(IINDA))
               End If
!
            End Do
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      End Do
!
      Return
!
End Subroutine R_mrgrnk
Subroutine I_mrgrnk (XDONT, IRNGT)
! __________________________________________________________
!   MRGRNK = Merge-sort ranking of an array
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
! __________________________________________________________
! __________________________________________________________
      Integer, Dimension (:), Intent (In)  :: XDONT
      Integer*4, Dimension (:), Intent (Out) :: IRNGT
! __________________________________________________________
      Integer :: XVALA, XVALB
!
      Integer, Dimension (SIZE(IRNGT)) :: JWRKT
      Integer :: LMTNA, LMTNC, IRNG1, IRNG2
      Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
!
      NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
      Select Case (NVAL)
      Case (:0)
         Return
      Case (1)
         IRNGT (1) = 1
         Return
      Case Default
         Continue
      End Select
!
!  Fill-in the index array, creating ordered couples
!
      Do IIND = 2, NVAL, 2
         If (XDONT(IIND-1) <= XDONT(IIND)) Then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         Else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         End If
      End Do
      If (Modulo(NVAL, 2) /= 0) Then
         IRNGT (NVAL) = NVAL
      End If
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      LMTNA = 2
      LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      Do
         If (NVAL <= 2) Exit
!
!   Loop on merges of A and B into C
!
         Do IWRKD = 0, NVAL - 1, 4
            If ((IWRKD+4) > NVAL) Then
               If ((IWRKD+2) >= NVAL) Exit
!
!   1 2 3
!
               If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
!
!   1 3 2
!
               If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
!
!   3 1 2
!
               Else
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               End If
               Exit
            End If
!
!   1 2 3 4
!
            If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
!
!   1 3 x x
!
            If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   1 3 2 4
                  IRNGT (IWRKD+3) = IRNG2
               Else
!   1 3 4 2
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               End If
!
!   3 x x x
!
            Else
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                  IRNGT (IWRKD+2) = IRNG1
                  If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   3 1 2 4
                     IRNGT (IWRKD+3) = IRNG2
                  Else
!   3 1 4 2
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  End If
               Else
!   3 4 1 2
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               End If
            End If
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 4
         Exit
      End Do
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is doubled.
!
      Do
         If (LMTNA >= NVAL) Exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
!
!   Loop on merges of A and B into C
!
         Do
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JINDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            If (IWRKF >= NVAL) Then
               If (JINDA >= NVAL) Exit
               IWRKF = NVAL
            End If
            IINDA = 1
            IINDB = JINDA + 1
!
!   Shortcut for the case when the max of A is smaller
!   than the min of B. This line may be activated when the
!   initial set is already close to sorted.
!
!          IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
!
!  One steps in the C subset, that we build in the final rank array
!
!  Make a copy of the rank array for the merge iteration
!
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
!
            XVALA = XDONT (JWRKT(IINDA))
            XVALB = XDONT (IRNGT(IINDB))
!
            Do
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               If (XVALA > XVALB) Then
                  IRNGT (IWRK) = IRNGT (IINDB)
                  IINDB = IINDB + 1
                  If (IINDB > IWRKF) Then
!  Only A still with unprocessed values
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                     Exit
                  End If
                  XVALB = XDONT (IRNGT(IINDB))
               Else
                  IRNGT (IWRK) = JWRKT (IINDA)
                  IINDA = IINDA + 1
                  If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                  XVALA = XDONT (JWRKT(IINDA))
               End If
!
            End Do
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      End Do
!
      Return
!
End Subroutine I_mrgrnk
end module m_mrgrnk
