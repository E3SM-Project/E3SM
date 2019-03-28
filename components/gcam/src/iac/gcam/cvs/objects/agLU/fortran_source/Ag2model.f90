SUBROUTINE Ag2Model( i, T, price2, qdem2, qsup2 )

! This version of the model subroutine is used when AgLU is a module
! This is called for each region at each time period

USE Ag2Global8

IMPLICIT NONE

! arrays set up just to pass the information out to Ag2Link
REAL(8) price2(10,NLP),qdem2(8,NLP,NMP),qsup2(8,NLP,NMP)

INTEGER i,T

price(:,i) = price2(:,i)


! Initialize Forests for the next time period for each region (already done for 1990)
IF (T.GT.2) THEN	
	Forest(3,i,T) = Forest(2,i,T-1)
	Forest(2,i,T) = Forest(1,i,T-1)
	Forest(1,i,T) = ForestN(i,T-1)
END IF

CALL Ag2Demand(i,T)
CALL Ag2Supply(i,T)

! Adjust for fixed trade in beefprds in base year and then phase out trade
!IF (T.EQ.2) qsup(3,i,T) = qsup(3,i,T) - NetExport(i,8)
!IF (T.EQ.3) qsup(3,i,T) = qsup(3,i,T) - NetExport(i,8)/2.0d0

qsup2(:,i,t) = qsup(:,i,t)
qdem2(:,i,t) = qdem(:,i,t)

!if( qsup( 8, i, t ) > 0 ) then
!write( *, * ) "QSUP2: ", qsup2( 8, i, t ),  " QSUP1: ", qsup( 8, i, t )
!endif

! Save prices for output
sprice(:,i,T) = price(:,i)

! IF (T .EQ. 7 .and. MODL1 .EQ. 21174) &
	! BREAK = 1.0d0

END SUBROUTINE Ag2Model
