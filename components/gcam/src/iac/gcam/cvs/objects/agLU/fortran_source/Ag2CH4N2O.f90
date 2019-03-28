	SUBROUTINE AG2CH4N2O(T,AG2CH4ACT,AG2N2OACT)

!	Calculation of Ag methane and nitrous oxide emissions

!   this routine is passed T (the period) and returns "activity levels"
!   for the ag emissions.  All cost curves, user control coefficients, etc.
!   are processed there using emiss_commoncode.f90.  This common code block could
!   be used in a standalone ag model as well if desired, although input cases would
!   need to be defined in the ag model for the cost curves and additional coefficients.

!   rewritten 7/02  mj
!   the "activity level" that is returned is really the same as the old way of getting
!   ag emissions -- the erb just processes additional controls (gdp based, abatement cost curves,
!   scaling exponent on activity, etc.)

!   IMPORTANT NOTE:
!	these arrays store some sort of "base" emissions, NOT final emissions which are
!   processed by the erb and contain additional calibration and adjustments   mj

!	AG2CH4(1,L,M)  :  Ruminant Enteric fermentation
!	AG2CH4(2,L,M)  :  Animal waste (managed)
!	AG2CH4(3,L,M)  :  Rice cultivation/other ag
!	AG2CH4(4,L,M)  :  Total CH4 by region

!	AG2N2O(1,L,M)  :  Soils/Fertilizer
!	AG2N2O(2,L,M)  :  Managed manure
!	AG2N2O(3,L,M)  :  empty 
!	AG2N2O(4,L,M)  :  Total N2O by region

USE Ag2Global8

IMPLICIT NONE

INTEGER mkt,i,T

! local variables: productivity of crop land multiplier, productivity of
!  pasture multiplier, mitigation cost curves,variables to pass the emissions
!  back to the MiniCAM
REAL(8) ferteff, meateff, AG2CH4ACT(1:3,1:NLP), AG2N2OACT(1:3,1:NLP)


AG2CH4ACT(:,:) = 0.0d0  ! initialize returned parameters
AG2N2OACT(:,:) = 0.0d0


IF (T .EQ. 1) RETURN	! no data for 1975 -- pass back 0's


DO i=1,NLP   ! note that i will index regions

	! Initialize the efficiency multipliers	   
	meateff = 1
	ferteff = 1
	! estimate efficiency improvements in meat prod and fertilizer
	! do we need any of these uncommented again?
	! note that they must be set to equal 1 in the base period, or base
	! period emissions will be wrong  mj

	!    IF (KJIMP(2,i,T) .NE. 0) THEN
	!      meateff = meateff * ((1+(KJIMP(2,i,T)*CYRATIO(T))) ** 15)

	!    IF (KJIMP(1,i,T) .NE. 0) THEN
	!      ferteff = ferteff * ((1+(KJIMP(1,i,T)*FYRATIO(T))) ** 15)


    ! *** Begin CH4 -- three sources ***

	! CH4 From Ruminant Enteric Fermentation
	mkt = 3 ! beef/mutton products market
	AG2CH4(1,i,T) = BaseCH4N2O(i,1) * (qsup(mkt,i,T)/qsup(mkt,i,2)) * meateff
	

	! CH4 From Animal Waste
	mkt = 3
    totanimprd(i,T) = qsup(mkt,i,T)+Porkprod(i,T)+poultryprod(i,T) ! total animal products
	AG2CH4(2,i,T) = BaseCH4N2O(i,2) * (totanimprd(i,T)/totanimprd(i,2)) * meateff
 
	! CH4 From Rice Fields
	mkt = 4 ! food grains market
	AG2CH4(3,i,T) = BaseCH4N2O(i,3) * (qsup(mkt,i,T)/qsup(mkt,i,2))


	! ****** Begin N2O -- three sources ******

	! N2O From Soils and Fertilizer
	totcropprd(i,T) = qsup(4,i,T)+qsup(5,i,T)+qsup(6,i,T)+qsup(7,i,T)  ! total crops
	AG2N2O(1,i,T) = BaseCH4N2O(i,4) * (totcropprd(i,T)/totcropprd(i,2)) * ferteff 

	! N2O From Managed Animal Waste
    AG2N2O(2,i,T) = BaseCH4N2O(i,5) * (totanimprd(i,T)/totanimprd(i,2)) * meateff

	! N2O From Unmanaged Animal Waste (not used, should we combine with soils? - mj)
    AG2N2O(3,i,T) = 0.0d0


	! REGIONAL TOTALS -- set to 0 since adding these "activity" numbers doesn't make sense  mj
	AG2CH4(4,i,T) = 0.0d0
	AG2N2O(4,i,T) = 0.0d0
	
	! Pass back the activities back to the MiniCAM
	AG2CH4ACT(1:3,i) = AG2CH4(1:3,i,T)
	AG2N2OACT(1:3,i) = AG2N2O(1:3,i,T)

  END DO

END SUBROUTINE Ag2CH4N2O

