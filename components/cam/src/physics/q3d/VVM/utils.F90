MODULE utils
!  a collection of utility routines (functions & subroutines)

   USE shr_kind_mod, only: dbl_kind => shr_kind_r8
   USE constld,      only: d0_0,d1_0,d2_0,d3_0,d6_0,d0_5
   
   IMPLICIT NONE
   PRIVATE

!    D0_0  = 0.0_dbl_kind
!    D1_0  = 1.0_dbl_kind
!    D2_0  = 2.0_dbl_kind
!    D3_0  = 3.0_dbl_kind
!    D6_0  = 6.0_dbl_kind
!    D0_5  = 0.5_dbl_kind

!------------------------------------------------------------------------------
! public member functions
!------------------------------------------------------------------------------
   PUBLIC ::    &
     es,        &  ! saturation vapor pressure (Pa)
     indexr,    &  ! finds the index of the element of a given table
     fintrp,    &  ! 1-D interpolation
     ran2,      &  ! random number generator
     cubicint,  &  ! 1-D cubic spline interpolation
     find_face, &  ! find the cube face number 
     latval,    &  ! find the latitude [rad]
     lonval,    &  ! find the longitude [rad]
     alphav,    &  ! find the coordinate alpha [rad]
     betav,     &  ! find the coordinate beta [rad]
     jaco,      &  ! find Jacobian of transformation (sqrt(G))
     gcont11,   &  ! find the contravariant metric tensor: g11
     gcont12,   &  ! find the contravariant metric tensor: g12
     gcont22,   &  ! find the contravariant metric tensor: g22
     amtx11,    &  ! find the transformation matrix A11
     amtx12,    &  ! find the transformation matrix A12
     amtx21,    &  ! find the transformation matrix A21
     amtx22,    &  ! find the transformation matrix A22 
     imatrix,   &  ! calculation of inverse matrix 
     cangle        ! calculation of central angle      
        
   CONTAINS

!========================================================================   
      FUNCTION es ( t ) RESULT(es_result)
!========================================================================      
         REAL (KIND=dbl_kind), INTENT(IN) :: t ! temperature (K)
         REAL(KIND=dbl_kind) :: es_result
         
         REAL (KIND=dbl_kind) :: esatw

         es_result = 100._dbl_kind * esatw ( t )

      END FUNCTION es

!========================================================================
      INTEGER FUNCTION INDEXR( V, NE, TE, LF )
!========================================================================      
!     FINDS THE INDEX OF THE ELEMENT OF TABLE TE THAT IS
!     JUST LESS THAN OR EQUAL TO V.
!
!     V  = VARIABLE
!     NE = NUMBER OF ELEMENTS IN TABLE TE
!     TE = TABLE OF ELEMENTS ARRANGED IN ASCENDING ORDER
!     LF = LOGICAL FLAG

      REAL(KIND=dbl_kind), INTENT(IN) :: V
      INTEGER, INTENT(IN) :: NE
      REAL(KIND=dbl_kind), INTENT(IN) :: TE(NE)
      
      LOGICAL, INTENT(INOUT) :: LF
      
      INTEGER :: I,J,ND

!     ORDER TEST
      IF ( TE(1) .GT. TE(2) ) GO TO 7
!
!     EXTREME TESTS
      IF ( V .GE. TE(1) ) GO TO 1
      LF = .TRUE.
      I = 1
      GO TO 14

   1  IF ( V .LT. TE(NE) ) GO TO 2
      LF = .TRUE.
      I = NE
      GO TO 14

!     INITIALIZATIONS
   2  LF = .FALSE.
      ND = 1

   3  ND = ND + ND
      IF ( ND .GE. NE ) GO TO 4
      GO TO 3

   4  ND = ND / 2
      I = ND

!     BISECTION LOOP
   5  ND = ND / 2
      IF ( ND .LE. 0 ) GO TO 6
      J = MIN( NE, I )
      IF ( V .GT. TE(J) ) I = I + ND
      IF ( V .LT. TE(J) ) I = I - ND
      GO TO 5

   6  IF ( I .GE. NE ) GO TO 14
      IF ( V .GT. TE(I+1) ) I = I + 1
      IF ( V .LT. TE(I) ) I = I - 1

      GO TO 14

!     EXTREME TESTS
   7  IF ( V .GE. TE(NE) ) GO TO 8
      LF = .TRUE.
      I = NE
      GO TO 14

   8  IF ( V .LT. TE(1) ) GO TO 9
      LF = .TRUE.
      I = 1
      GO TO 14

!     INITIALIZATIONS
   9  LF = .FALSE.
      ND = 1

  10  ND = ND + ND
      IF ( ND .GE. NE ) GO TO 11
      GO TO 10

  11  ND = ND / 2
      I = ND

!     BISECTION LOOP
  12  ND = ND / 2
      IF ( ND .LE. 0 ) GO TO 13
      J = MIN( NE, I )
      IF ( V .GT. TE(J) ) I = I - ND
      IF ( V .LT. TE(J) ) I = I + ND
      GO TO 12

  13  IF ( I .GE. NE ) GO TO 14
      IF ( I .LE. 1 ) GO TO 14
      IF ( V .GT. TE(I-1) ) I = I - 1
      IF ( V .LT. TE(I) ) I = I + 1
      I = MAX( 1, I - 1 )

  14  INDEXR = MIN0( NE - 1, I )

      END FUNCTION indexr

!========================================================================      
      FUNCTION FINTRP ( MODE, X, X1, F1, X2, F2 ) RESULT(fintrp_result)
!========================================================================      
!         *FINTRP* INTERPOLATES BETWEEN X1 AND X2 USING A FUNCTIONAL
!         FORM THAT DEPENDS UPON THE SPECIFIED MODE. MODES - -
!     1 = LINEAR - - A * X + B
!     2 = EXPONENTIAL - - A * EXP( B * X )
!     3 = POWER LAW - - A * X**B

         INTEGER, INTENT(IN) :: mode
         REAL (KIND=dbl_kind), INTENT(IN) :: X,X1,F1,X2,F2
         
         REAL (KIND=dbl_kind) :: fintrp_result

         SELECT CASE (mode)
         
!         LINEAR INTERPOLATION
            CASE(1)
               fintrp_result = F1 + ( F1 - F2 ) * ( X - X1 ) / ( X1 - X2 )

!         EXPONENTIAL INTERPOLATION
            CASE(2)
               fintrp_result = F1 * ( F1 / F2 )**( ( X - X1 ) / ( X1 - X2 ) )

!         POWER LAW INTERPOLATION
            CASE(3)
              IF ( X1 .GT. D0_0 ) then
                 fintrp_result = F1 * ( X / X1 )**( LOG( F1 / F2 ) / LOG( X1 / X2 ) )
              ELSE
                 fintrp_result = F1 * ( F1 / F2 )**( ( X - X1 ) / ( X1 - X2 ) )
              ENDIF
            CASE DEFAULT
              WRITE(6,*) mode, 'is an improper value of mode, stopping'
              STOP
         END SELECT

      END FUNCTION fintrp

!========================================================================
      FUNCTION RAN2(IDUM) RESULT(ran2_result)
!========================================================================      
!     Returns a random number between -1.0 and 1.0. 

         INTEGER, INTENT(INOUT) :: idum
         REAL (KIND=dbl_kind) :: ran2_result
         
         ! Local variables
         INTEGER, parameter :: M=714025,IA=1366,IB=3454,IC=150889
         INTEGER :: m1, m2, m3, mfac
         REAL (KIND=dbl_kind) :: tmp, rfac
      
         m1 = mod(IA*IDUM,m)
         m2 = mod(IB*IDUM,m)
         m3 = mod(IC*IDUM,m)
         mfac = (m1+m2+m3)/3
         if(mfac .eq. 0) mfac=1
         tmp = mod(m1*m2+m3, mfac)
         rfac = mfac
         ran2_result = tmp/rfac
         IDUM = IDUM+1

      END FUNCTION ran2

!========================================================================
      SUBROUTINE CUBICINT (NSIZE,NSIZE_NEW,XVAL,YVAL,XNEW,YNEW)
!========================================================================      
      INTEGER,INTENT(IN) :: NSIZE,NSIZE_NEW
      
      REAL (KIND=dbl_kind),INTENT(IN)  :: XVAL(NSIZE),YVAL(NSIZE)
      REAL (KIND=dbl_kind),INTENT(IN)  :: XNEW(NSIZE_NEW)
      REAL (KIND=dbl_kind),INTENT(OUT) :: YNEW(NSIZE_NEW)
      
      ! Local variables
      REAL (KIND=dbl_kind) :: YP1,YPN    ! First derivative of the interpolating function 
                                         ! at points 1 and n, respectively.
      
      REAL (KIND=dbl_kind),DIMENSION(NSIZE) :: Y2
      REAL (KIND=dbl_kind) :: X_INT,Y_INT
      INTEGER :: I
      
      YP1 = (YVAL(2)-YVAL(1))/(XVAL(2)-XVAL(1))
      YPN = (YVAL(NSIZE)-YVAL(NSIZE-1))/(XVAL(NSIZE)-XVAL(NSIZE-1)) 
         
      CALL SPLINE(XVAL,YVAL,NSIZE,YP1,YPN,Y2)

      DO I = 1, NSIZE_NEW
       X_INT =  XNEW(I)
       CALL SPLINT(XVAL,YVAL,Y2,NSIZE,X_INT,Y_INT)
       YNEW(I) =  Y_INT
      ENDDO

      CONTAINS

      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
!     This routine returns an array Y2 (= second derivatives)      
!     yp1 or ypn=0: second derivative is zero (natural cubic spline) 

      INTEGER,INTENT(IN) :: N
      REAL (KIND=dbl_kind),INTENT(IN) :: YP1,YPN
      REAL (KIND=dbl_kind),DIMENSION(N),INTENT(IN) :: X,Y
      REAL (KIND=dbl_kind),DIMENSION(N),INTENT(OUT) :: Y2
      
      REAL (KIND=dbl_kind),DIMENSION(N):: U(N)
      REAL (KIND=dbl_kind) :: SIG,P,QN,UN
      INTEGER :: I,K
      
      IF (YP1 .GT. .99E30) THEN
        Y2(1) = D0_0
        U(1)  = D0_0
      ELSE
        Y2(1) = -D0_5
        U(1)  = (D3_0/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      
      DO I = 2, N-1
        SIG = (X(I)-X(I-1))/(X(I+1)-X(I-1))
        P = SIG*Y2(I-1)+2.
        Y2(I) = (SIG-1.)/P
        U(I) = (D6_0*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1)) &
              /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
      ENDDO
      
      IF (YPN .GT. .99E30) THEN
        QN = D0_0
        UN = D0_0
      ELSE
        QN = D0_5
        UN = (D3_0/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      
      Y2(N) = (UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO K = N-1, 1, -1
        Y2(K) = Y2(K)*Y2(K+1)+U(K)
      ENDDO
      
      END SUBROUTINE spline

      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
!     Returns a cubic-spline interpolated value Y at X-point 

      INTEGER,INTENT(IN) :: N
      REAL (KIND=dbl_kind),DIMENSION(N),INTENT(IN) :: XA,YA,Y2A
      REAL (KIND=dbl_kind),INTENT(IN)  :: X
      REAL (KIND=dbl_kind),INTENT(OUT) :: Y
      
      REAL (KIND=dbl_kind) :: H,A,B
      INTEGER :: I,KLO,KHI,K
      
      KLO = 1
      KHI = N
1     IF (KHI-KLO .GT. 1) THEN
        K = (KHI+KLO)/2
        IF(XA(K) .GT. X)THEN
          KHI = K
        ELSE
          KLO = K
        ENDIF
      GOTO 1
      ENDIF
      
      H = XA(KHI) - XA(KLO)
      IF (H .EQ. D0_0) PAUSE 'Bad XA input.'
      A = (XA(KHI)-X)/H
      B = (X-XA(KLO))/H
      Y = A*YA(KLO)+B*YA(KHI) &
        +((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/D6_0
      
      END SUBROUTINE splint 
             
      END SUBROUTINE cubicint 

!======================================================================== 
      FUNCTION FIND_FACE ( A, B ) RESULT(FIND_FACE_result)
!========================================================================
!     Find the cube-face on which the point (A,B) is located.
 
      REAL (KIND=dbl_kind), INTENT(IN) ::  A   ! longitude [rad]
      REAL (KIND=dbl_kind), INTENT(IN) ::  B   ! latitude  [rad]
      INTEGER :: FIND_FACE_result
      
      REAL (KIND=dbl_kind) :: X,Y,Z,AX,AY,AZ,PM
      
      X = DCOS(A)*DCOS(B)
      Y = DSIN(A)*DCOS(B)
      Z = DSIN(B)

      AX = ABS(X)
      AY = ABS(Y)
      AZ = ABS(Z)
      
      PM = MAX(AX,AY,AZ)
      
      IF (PM .EQ. AX) THEN
         if (x .gt. d0_0) then
           FIND_FACE_result = 1      ! F1 
         else
           FIND_FACE_result = 3      ! F3
         endif
      ELSE IF (PM .EQ. AY) THEN   
         if (y .gt. d0_0) then
           FIND_FACE_result = 2      ! F2
         else
           FIND_FACE_result = 4      ! F4
         endif      
      ELSE
         if (z .gt. d0_0) then
           FIND_FACE_result = 5      ! F5 (top face) 
         else
           FIND_FACE_result = 6      ! F6 (bottom face)  
         endif            
      ENDIF  
       
      END FUNCTION find_face

!========================================================================          
      FUNCTION LATVAL ( N, A, B ) RESULT(LATVAL_result)
!========================================================================
!     INPUT : N (Face number), A (alpha), B (beta)
!     OUTPUT: latitude in [rad]
!
!     FACE1: Latitude = -45 ~  45 (alpha=0)
!     FACE2: Latitude = -45 ~  45 (alpha=0)
!     FACE3: Latitude = -45 ~  45 (alpha=0)
!     FACE4: Latitude = -45 ~  45 (alpha=0)
!     FACE5: Latitude =  45 ~  90 (alpha=0)
!     FACE6: Latitude = -45 ~ -90 (alpha=0)
!========================================================================  
      INTEGER, INTENT(IN) :: N
      REAL (KIND=dbl_kind), INTENT(IN) ::  A,B
      REAL (KIND=dbl_kind) :: LATVAL_result
      
      REAL (KIND=dbl_kind) :: COS_A,TAN_A,TAN_B

      COS_A  = DCOS(A)
      TAN_A  = DTAN(A)
      TAN_B  = DTAN(B)
         
      SELECT CASE (N)
        CASE(1:4)
          LATVAL_result = DATAN(TAN_B*COS_A)    
        CASE(5)
          LATVAL_result =  DATAN(D1_0/SQRT(TAN_A**2 + TAN_B**2))
        CASE(6)
          LATVAL_result = -DATAN(D1_0/SQRT(TAN_A**2 + TAN_B**2))
      END SELECT   
       
      END FUNCTION latval

!========================================================================       
      FUNCTION LONVAL ( N, A, B, PIVAL ) RESULT(LONVAL_result)
!========================================================================
!     INPUT : N (Face number), A (alpha), B (beta)
!     OUTPUT: longitude in [rad]
!
!     FACE1: Longitude = 225 ~ 315 deg 
!     FACE2: Longitude = -45 ~ 45  deg 
!     FACE3: Longitude =  45 ~ 135 deg 
!     FACE4: Longitude = 135 ~ 225 deg 
!     FACE5: Longitude = 0 ~ +/-180 deg 
!            At NP (alpha=beta=0), not determined: DATAN2(0,0)=0   
!     FACE6: Longitude = 0 ~ +/-180 deg 
!            At SP (alpha=beta=0), not determined: DATAN2(0,0)=0  
!======================================================================== 
      INTEGER, INTENT(IN) :: N
      REAL (KIND=dbl_kind), INTENT(IN) ::  A,B,PIVAL
      REAL (KIND=dbl_kind) :: LONVAL_result
      
      REAL (KIND=dbl_kind) :: TAN_A,TAN_B

      TAN_A = DTAN(A)
      TAN_B = DTAN(B)
         
      SELECT CASE (N)
        CASE(1:4)
          LONVAL_result = A + D0_5*(N-1)*PIVAL  
        CASE(5)
          LONVAL_result = DATAN2(TAN_A,-TAN_B) 
        CASE(6)
          LONVAL_result = DATAN2(TAN_A, TAN_B) 
      END SELECT   
       
      END FUNCTION lonval    

!======================================================================== 
      FUNCTION ALPHAV ( N, LONV, LATV, PIVAL ) RESULT(ALPHAV_result)
!========================================================================      
!     INPUT : N (face number), LONV (longitude [rad]), LATV (latitude [rad]), PI
!     OUTPUT: alpha value in [rad]
!========================================================================
      INTEGER, INTENT(IN) :: N
      REAL (KIND=dbl_kind), INTENT(IN) ::  LONV,LATV,PIVAL
      REAL (KIND=dbl_kind) :: ALPHAV_result
      
      REAL (KIND=dbl_kind) :: TAN_LAT,SIN_LON

      TAN_LAT  = DTAN(LATV)
      SIN_LON  = DSIN(LONV)
         
      SELECT CASE (N)
        CASE(1:4)
          ALPHAV_result = LONV - D0_5*(N-1)*PIVAL 
        CASE(5)
          ALPHAV_result = DATAN(SIN_LON/TAN_LAT)
        CASE(6)
          ALPHAV_result = DATAN(-SIN_LON/TAN_LAT)
      END SELECT   
       
      END FUNCTION alphav

!========================================================================             
      FUNCTION BETAV ( N, LONV, LATV, PIVAL ) RESULT(BETAV_result)
!========================================================================       
!     INPUT : N (face number), LONV (longitude [rad]), LATV (latitude [rad]), PI
!     OUTPUT: beta value in [rad]
!========================================================================  
      INTEGER, INTENT(IN) :: N
      REAL (KIND=dbl_kind), INTENT(IN) ::  LONV,LATV,PIVAL
      REAL (KIND=dbl_kind) :: BETAV_result
      
      REAL (KIND=dbl_kind) :: TAN_LAT,COS_LON,TEMP,SECV

      TAN_LAT  = DTAN(LATV)
      COS_LON  = DCOS(LONV)
      
      TEMP = LONV - D0_5*(N-1)*PIVAL
      SECV = D1_0/DCOS(TEMP)
         
      SELECT CASE (N)
        CASE(1:4)
          BETAV_result = DATAN(TAN_LAT*SECV)    
        CASE(5)
          BETAV_result = DATAN(-COS_LON/TAN_LAT)
        CASE(6)
          BETAV_result = DATAN(-COS_LON/TAN_LAT)
      END SELECT   
       
      END FUNCTION betav

!========================================================================               
      FUNCTION JACO ( A, B ) RESULT(JACO_result)
!========================================================================
!     CALCULATION OF Jacobian of Transformation: SQRT(G)
!     INPUT: A (alpha value), B (beta value)
!========================================================================       
 
      REAL (KIND=dbl_kind), INTENT(IN) ::  A,B
      REAL (KIND=dbl_kind) :: JACO_result
      
      REAL (KIND=dbl_kind) :: COS_A2,COS_B2,TAN_A2,TAN_B2
      REAL (KIND=dbl_kind) :: GAMMA,GAMMA3

      COS_A2 = DCOS(A)**2
      COS_B2 = DCOS(B)**2
      
      TAN_A2 = DTAN(A)**2
      TAN_B2 = DTAN(B)**2
      
      GAMMA  = SQRT(D1_0+ TAN_A2 + TAN_B2)
      GAMMA3 = GAMMA**3
         
      JACO_result = D1_0/(GAMMA3*COS_A2*COS_B2)   
       
      END FUNCTION jaco

! CALCULATION OF CONTRAVARIANT METRIC TENSOR: g11, g12=g21, g22
!======================================================================== 
      FUNCTION GCONT11 ( A, B ) RESULT(GCONT11_result)  
!========================================================================       
!     INPUT : A (alpha value), B (beta value)
!     OUTPUT: G11
 
      REAL (KIND=dbl_kind), INTENT(IN) ::  A,B
      REAL (KIND=dbl_kind) :: GCONT11_result
      
      REAL (KIND=dbl_kind) :: COS_A2,COS_B2,TAN_A2,TAN_B2,GAMMA2

      COS_A2 = DCOS(A)*DCOS(A)
      COS_B2 = DCOS(B)*DCOS(B)
      TAN_A2 = DTAN(A)**2
      TAN_B2 = DTAN(B)**2
      GAMMA2 = D1_0 + TAN_A2 + TAN_B2
         
      GCONT11_result = GAMMA2*COS_A2*COS_B2*(D1_0+TAN_B2)
       
      END FUNCTION gcont11      

!======================================================================== 
      FUNCTION GCONT12 ( A, B ) RESULT(GCONT12_result)
!========================================================================       
!     INPUT : A (alpha value), B (beta value)
!     OUTPUT: G12 (=G21)
 
      REAL (KIND=dbl_kind), INTENT(IN) ::  A,B
      REAL (KIND=dbl_kind) :: GCONT12_result
      
      REAL (KIND=dbl_kind) :: COS_A2,COS_B2,TAN_A,TAN_B,GAMMA2

      COS_A2 = DCOS(A)*DCOS(A)
      COS_B2 = DCOS(B)*DCOS(B)
      TAN_A  = DTAN(A)
      TAN_B  = DTAN(B)
      GAMMA2 = D1_0 + TAN_A**2 + TAN_B**2
         
      GCONT12_result = GAMMA2*COS_A2*COS_B2*TAN_A*TAN_B
       
      END FUNCTION gcont12      

!======================================================================== 
      FUNCTION GCONT22 ( A, B ) RESULT(GCONT22_result)
!========================================================================       
!     INPUT : A (alpha value), B (beta value)
!     OUTPUT: G22
 
      REAL (KIND=dbl_kind), INTENT(IN) ::  A,B
      REAL (KIND=dbl_kind) :: GCONT22_result
      
      REAL (KIND=dbl_kind) :: COS_A2,COS_B2,TAN_A2,TAN_B2,GAMMA2

      COS_A2 = DCOS(A)*DCOS(A)
      COS_B2 = DCOS(B)*DCOS(B)
      TAN_A2 = DTAN(A)**2
      TAN_B2 = DTAN(B)**2
      GAMMA2 = D1_0+ TAN_A2 + TAN_B2
         
      GCONT22_result = GAMMA2*COS_A2*COS_B2*(D1_0+TAN_A2)
       
      END FUNCTION gcont22      

! CALCULATION OF MATRIX A (EACH COMPONENTS)
!======================================================================== 
      FUNCTION AMTX11 ( N, A, B ) RESULT(AMTX11_result)
!========================================================================      
!     INPUT : N (Face number), A (alpha value), B (beta value)
!     OUTPUT: a11
 
      INTEGER, INTENT(IN) :: N
      REAL (KIND=dbl_kind), INTENT(IN) ::  A,B
      REAL (KIND=dbl_kind) :: AMTX11_result
      
      REAL (KIND=dbl_kind) :: COS_A,TAN_A,TAN_B
      REAL (KIND=dbl_kind) :: COS_A2,TAN_A2,TAN_B2,GAMMA,RLON
      ! RLON (Longitude)

      COS_A = DCOS(A)  
      TAN_A = DTAN(A)
      TAN_B = DTAN(B)

      COS_A2 = COS_A*COS_A
      TAN_A2 = TAN_A*TAN_A
      TAN_B2 = TAN_B*TAN_B 
            
      GAMMA = DSQRT(D1_0 + TAN_A2 + TAN_B2)
         
      SELECT CASE (N)
        CASE(1:4)
          AMTX11_result = D1_0/(GAMMA*COS_A)     
        CASE(5)
          RLON = DATAN2(TAN_A,-TAN_B)
          AMTX11_result = DCOS(RLON)/(GAMMA*COS_A2)
        CASE(6)
          RLON = DATAN2(TAN_A,TAN_B) 
          AMTX11_result = DCOS(RLON)/(GAMMA*COS_A2) 
      END SELECT   
       
      END FUNCTION amtx11

!========================================================================      
      FUNCTION AMTX12 ( N, A, B ) RESULT(AMTX12_result)
!========================================================================      
!     INPUT : N (Face number), A (alpha value), B (beta value)
!     OUTPUT: a12
 
      INTEGER, INTENT(IN) :: N
      REAL (KIND=dbl_kind), INTENT(IN) ::  A,B
      REAL (KIND=dbl_kind) :: AMTX12_result
      
      REAL (KIND=dbl_kind) :: COS_B,TAN_A,TAN_B
      REAL (KIND=dbl_kind) :: COS_B2,TAN_A2,TAN_B2,GAMMA,RLON

      COS_B = DCOS(B)
      TAN_A = DTAN(A)
      TAN_B = DTAN(B)

      COS_B2 = COS_B*COS_B
      TAN_A2 = TAN_A*TAN_A
      TAN_B2 = TAN_B*TAN_B 
          
      GAMMA = DSQRT(D1_0 + TAN_A2 + TAN_B2)
               
      SELECT CASE (N)
        CASE(1:4)
          AMTX12_result = D0_0 
        CASE(5)
          RLON = DATAN2(TAN_A,-TAN_B)
          AMTX12_result =  DSIN(RLON)/(GAMMA*COS_B2)
        CASE(6)
          RLON = DATAN2(TAN_A,TAN_B) 
          AMTX12_result = -DSIN(RLON)/(GAMMA*COS_B2)
      END SELECT   
       
      END FUNCTION amtx12

!========================================================================
      FUNCTION AMTX21 ( N, A, B ) RESULT(AMTX21_result)
!========================================================================      
!     INPUT : N (Face number), A (alpha value), B (beta value)
!     OUTPUT: a21
 
      INTEGER, INTENT(IN) :: N
      REAL (KIND=dbl_kind), INTENT(IN) ::  A,B
      REAL (KIND=dbl_kind) :: AMTX21_result
      
      REAL (KIND=dbl_kind) :: SIN_B,COS_A,COS_B,TAN_A,TAN_B
      REAL (KIND=dbl_kind) :: COS_A2,TAN_A2,TAN_B2,GAMMA2,RLON
      
      SIN_B  = DSIN(B)
      COS_A  = DCOS(A)
      COS_B  = DCOS(B)
      TAN_A  = DTAN(A)
      TAN_B  = DTAN(B)
      
      COS_A2 = COS_A*COS_A
      TAN_A2 = TAN_A*TAN_A
      TAN_B2 = TAN_B*TAN_B 
             
      GAMMA2 = D1_0 + TAN_A2 + TAN_B2
                     
      SELECT CASE (N)
        CASE(1:4)
          AMTX21_result = -(TAN_A*SIN_B)/(GAMMA2*COS_A*COS_B)  
        CASE(5)
          RLON = DATAN2(TAN_A,-TAN_B)
          AMTX21_result = -DSIN(RLON)/(GAMMA2*COS_A2)
        CASE(6)
          RLON = DATAN2(TAN_A,TAN_B) 
          AMTX21_result =  DSIN(RLON)/(GAMMA2*COS_A2)
      END SELECT   
       
      END FUNCTION amtx21

!========================================================================      
      FUNCTION AMTX22 ( N, A, B ) RESULT(AMTX22_result)
!========================================================================      
!     INPUT : N (Face number), A (alpha value), B (beta value)
!     OUTPUT: a22
 
      INTEGER, INTENT(IN) :: N
      REAL (KIND=dbl_kind), INTENT(IN) ::  A,B
      REAL (KIND=dbl_kind) :: AMTX22_result
      
      REAL (KIND=dbl_kind) :: COS_A,COS_B,TAN_A,TAN_B
      REAL (KIND=dbl_kind) :: COS_B2,TAN_A2,TAN_B2,GAMMA2,RLON

      COS_A  = DCOS(A)
      COS_B  = DCOS(B)
      TAN_A  = DTAN(A)
      TAN_B  = DTAN(B)
      
      COS_B2 = COS_B*COS_B
      TAN_A2 = TAN_A*TAN_A
      TAN_B2 = TAN_B*TAN_B 
      
      GAMMA2 = D1_0 + TAN_A2 + TAN_B2
               
      SELECT CASE (N)
        CASE(1:4)
          AMTX22_result = D1_0/(GAMMA2*COS_A*COS_B2)    
        CASE(5)
          RLON = DATAN2(TAN_A,-TAN_B)
          AMTX22_result = DCOS(RLON)/(GAMMA2*COS_B2)
        CASE(6)
          RLON = DATAN2(TAN_A,TAN_B) 
          AMTX22_result = DCOS(RLON)/(GAMMA2*COS_B2)
      END SELECT   
       
      END FUNCTION amtx22                

!======================================================================== 
      FUNCTION IMATRIX ( MAT ) RESULT(IMATRIX_result) 
!========================================================================
!     CALCULATION OF INVERSE MATRIX for (2x2) matrix
!     INPUT : (2x2) matrix A
!     OUTPUT: (2x2) inverse matrix of A
!======================================================================== 
      
      REAL (KIND=dbl_kind), DIMENSION(2,2), INTENT(IN) :: MAT
      REAL (KIND=dbl_kind), DIMENSION(2,2) :: IMATRIX_result
      
      REAL (KIND=dbl_kind), DIMENSION(2,4) :: A
      REAL (KIND=dbl_kind) :: PivElt, TarElt,TempRow(4)
      INTEGER :: N,PivRow,TarRow,I,J,K 
      
      A = D0_0
      DO J = 1,2
       DO I = 1,2
        A(I,J) = MAT(I,J)
       ENDDO
      ENDDO
      DO I = 1,2
       A(I,2+I) = D1_0
      ENDDO
      
      DO PivRow = 1,2 
       PivElt = A(PivRow,PivRow)
       IF (PivElt == D0_0) THEN
        K = PivRow + 1
        DO WHILE (PivElt == D0_0 .AND. K <= 2)
          PivElt = A(K,PivRow)
          K = K+1
        ENDDO
        IF (PivElt == D0_0) THEN
          PRINT*, "Couldn't find a non-zero pivot: solution rubbish"
          RETURN 
        ELSE
          TempRow = A(PivRow,1:4)
          K = K-1
          A(PivRow,1:4) = A(K,1:4)
          A(K,1:4) = TempRow
        ENDIF
       ENDIF 
        A(PivRow,1:4) = A(PivRow,1:4)/PivElt
        
        DO TarRow = 1,2
          IF (TarRow /= PivRow) THEN
            TarElt = A(TarRow,PivRow)
            A(TarRow,1:4) = A(TarRow,1:4) - A(PivRow,1:4)*TarElt
          ENDIF
        ENDDO
       ENDDO
       
       IMATRIX_result = A(1:2,3:4)          
        
      END FUNCTION imatrix  

!========================================================================       
      FUNCTION CANGLE ( PI,LAT1,LON1,LAT2,LON2 ) RESULT(CANGLE_result)
!========================================================================      
!     Calculate the central angle between two points on the sphere using the Vincenty 
!     formula that is accurate for all distances (the usual spherical law of cosines 
!     can give round errors for the small distances).  
!
!     INPUT : PI, latitudes [rad], longitudes [rad]
!     OUTPUT: central angle between two points     
!========================================================================      
      REAL (KIND=dbl_kind), INTENT(IN) :: PI,LAT1,LON1,LAT2,LON2
      REAL (KIND=dbl_kind) :: CANGLE_result
      
      REAL (KIND=dbl_kind) :: LON1new,LON2new,DLAMDA
      REAL (KIND=dbl_kind) :: TERM1,TERM2,TERM3,TERM4
      
      LON1new = LON1
      LON2new = LON2
      
      IF (LON1.lt.D0_0) LON1new = LON1new + D2_0*PI
      IF (LON2.lt.D0_0) LON2new = LON2new + D2_0*PI
       
      DLAMDA   = DABS(LON2new-LON1new)
      
      TERM1 = DCOS(LAT2)*DSIN(DLAMDA)
      TERM2 = DCOS(LAT1)*DSIN(LAT2) - DSIN(LAT1)*DCOS(LAT2)*DCOS(DLAMDA)
      TERM3 = DSIN(LAT1)*DSIN(LAT2) + DCOS(LAT1)*DCOS(LAT2)*DCOS(DLAMDA)
      TERM4 = SQRT(TERM1**2 + TERM2**2)
      
      CANGLE_result = DATAN2(TERM4,TERM3) 
      
      END FUNCTION cangle 

END MODULE utils
