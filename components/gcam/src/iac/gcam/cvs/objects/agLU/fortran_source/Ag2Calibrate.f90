SUBROUTINE Ag2Calibrate

! Procedure incorporating the calculations formerly done in the
! input spreadsheet in sheets parameters1, parameters2, Land2,
! and in the 11 regional sheets.  This subroutine returns data
! into GJ, RJ0, PASTOUT and CALP.
! As well, the parameters of the animal prds production fuction
! are calculated
! There are 3 levels of organization.  Level 0 divides up unmanaged
! and managed land.  Level 1 divides up the different land uses.
! Level 2 divides up the individual crops in the cropland category
! written 5/01 ktg

! Modified 7/01 to change how biomass is calibrated - ktg
! biomass variable cost (GJ) and Biomass avg observed yield are now read-in

USE Ag2Global8

IMPLICIT NONE

! local parameters: average observed yield of unmanaged land, 
!  observed profit rate per hectare, Animal products starting
!  intrinsic yield, start of variable cost 3D array, number of sectors in var cost
REAL, PARAMETER :: Unmanyield=1.0d0, obsprof=100d0, Animyield=1.0d0
INTEGER, PARAMETER :: MSTART = 2, numsectors = 9

! local variables: rate of return factor,
!  biomass area harvested, Total Crop Area (biomass+crops), Total Mananged Land
!  Total Unmanaged Land, foodgr ratio of avg obs yield to intrinsic yield
REAL(8) Factor, CropArea(NLP), ManLand(NLP), UnmanLand(NLP), yieldratio(NLP)

! local integers: sectors in variable cost, 3rd dimension in GJ 3D array
INTEGER	sector, M

! local 2D arrays: Adjusted stdev, unobserved profit rate per hectare calc factor, share of
!  land in nest in the 3 levels, profit per output for the 3 levels, profit per ha for the
!  3 levels, other cost per unit of output for level 1, price per unit output,
!  average yield for levels 1 and 2, forest yield for 2035,forest size requirement
REAL(8) AdjSD(NLP,3), proffactor(NLP,3), shareland0(NLP,2), shareland1(NLP,3), &
	shareland2(NLP,5), profperout0(NLP,1), profperout1(NLP,2), profperout2(NLP,5), &
	profperha0(NLP,2), profperha1(NLP,3), profperha2(NLP,5), othcostperout1(NLP,3), &
	priceperout1(NLP,4), avgyield1(NLP,2), avgyield2(NLP,5), &
	fyield(NLP), forestreq(NLP)

! i is the region we are working with
INTEGER i

! variables for the calibration of beef/mutton substitution curve
INTEGER k,ncrops

REAL(8) ShareFeed(6),p_feed,p_crop(4),TempSum,p_AnimPrds,p_Pasture


! First calibrate the forest products (down at the bottom of calibration.f90)
CALL ForestCal(fyield,forestreq)

Factor = InterestR / ((1.0d0 + InterestR)**45.0d0 - 1.0d0)

  DO i = 1,NLP

  ! Land Area breakup calculations

	! Calc area harvested numbers for Level 2
	CropArea(i) = AreaHarv(i,1)+AreaHarv(i,2)+AreaHarv(i,3)+AreaHarv(i,4)+AreaHarv(i,5)
	
	! calc land area numbers for Level 1
	ManLand(i) = forestreq(i) + LAND(i,2) + CropArea(i)

	! calc land area numbers for Level 0
	UnmanLand(i) = LAND(i,6) - ManLand(i)
	
  ! Calculate share of land in nest

	! calc for Level 0
	shareland0(i,1) = UnmanLand(i) / LAND(i,6)
	shareland0(i,2) = ManLand(i) / LAND(i,6)

	! calc for Level 1 (man forest,pasture & crop/biomass)
	shareland1(i,1) = Forestreq(i) / ManLand(i)
	shareland1(i,2) = LAND(i,2) / ManLand(i)
	shareland1(i,3) = CropArea(i) / ManLand(i)

	! calc for Level 2
	shareland2(i,1) = AreaHarv(i,5) / CropArea(i)
	shareland2(i,2) = AreaHarv(i,1) / CropArea(i)
	shareland2(i,3) = AreaHarv(i,2) / CropArea(i)
	shareland2(i,4) = AreaHarv(i,3) / CropArea(i)
	shareland2(i,5) = AreaHarv(i,4) / CropArea(i)


  ! Level 0 Calculations (Managed vs. Unmanaged land)

	! Calc profit factor
	AdjSD(i,1) = SIG(i) * SQRT(1.0d0-CORRS(i,1))
	proffactor(i,1) = obsprof ** (1.0d0 / AdjSD(i,1))

	! Calc profit per hectare for unmanaged and managed land
	profperha0(i,1) = (proffactor(i,1) * shareland0(i,1)) ** AdjSD(i,1)
	profperha0(i,2) = (proffactor(i,1) * shareland0(i,2)) ** AdjSD(i,1)

	! Calc profit per unit of output for unmanaged land
	profperout0(i,1) = obsprof / Unmanyield

	! Calc intrinsic yield for unmanaged land
	RJ0(1,i) = profperha0(i,1) / profperout0(i,1)

	! Calc unmanaged land's calibration price
	CALP(i,9) = obsprof * RJ0(1,i)

! BREAK = 1.0d0

  ! Level 1 Calculations (land uses)
	
	! Calc profit factor
	AdjSD(i,2) = SIG(i) * SQRT(1.0d0-CORRS(i,2))
	proffactor(i,2) = profperha0(i,2) ** (1.0d0 / AdjSD(i,2))

	! Calc profit per hectare for man forest,pasture & crops/biomass
	profperha1(i,1) = (proffactor(i,2) * shareland1(i,1)) ** AdjSD(i,2)
	profperha1(i,2) = (proffactor(i,2) * shareland1(i,2)) ** AdjSD(i,2)
	profperha1(i,3) = (proffactor(i,2) * shareland1(i,3)) ** AdjSD(i,2)

	! Set average yield for managed forest (in cm/ha)
	avgyield1(i,1) = fyield(i) 

	! Calc profit per unit of output for managed forest (pasture needs lvl 2 data)
	profperout1(i,1) = obsprof / (Factor * avgyield1(i,1))

	! Calc intrinsic yield for managed forest
	RJ0(2,i) = profperha1(i,1) / (Factor * profperout1(i,1))

	! Calc other cost per unit output for managed forest (variable cost)
	priceperout1(i,1) = CALP(i,2)
	othcostperout1(i,1) = priceperout1(i,1) - profperout1(i,1)
	GJ(2,i,2) = othcostperout1(i,1)


  ! Level 2 Calculations (to provide the data to finish level 1)
	
	! Calc profit factor
	AdjSD(i,3) = SIG(i) * SQRT(1.0d0-CORRS(i,3))
	proffactor(i,3) = profperha1(i,3) ** (1.0d0 / AdjSD(i,3))

	! Calc profit per hectare for biomass and the 4 food crops
!	profperha2(i,1) = (proffactor(i,3) * shareland2(i,1)) ** AdjSD(i,3) - old biomass calculation
	profperha2(i,2) = (proffactor(i,3) * shareland2(i,2)) ** AdjSD(i,3)
	profperha2(i,3) = (proffactor(i,3) * shareland2(i,3)) ** AdjSD(i,3)
	profperha2(i,4) = (proffactor(i,3) * shareland2(i,4)) ** AdjSD(i,3)
	profperha2(i,5) = (proffactor(i,3) * shareland2(i,5)) ** AdjSD(i,3)

	! Calc average yield for the 4 food crops (formerly in the Land2 wksheet)
	avgyield2(i,2) = CropSupply(i,1) / AreaHarv (i,1)
	avgyield2(i,3) = CropSupply(i,2) / AreaHarv (i,2)
	avgyield2(i,4) = CropSupply(i,3) / AreaHarv (i,3)
	avgyield2(i,5) = CropSupply(i,4) / AreaHarv (i,4)

	! Calc profit per unit of output for the 4 food crops (in $/ha)
	profperout2(i,2) = obsprof / avgyield2(i,2)
	profperout2(i,3) = obsprof / avgyield2(i,3)
	profperout2(i,4) = obsprof / avgyield2(i,4)
	profperout2(i,5) = obsprof / avgyield2(i,5)	
	
	! Calc intrinsic yield for the 4 food crops (in Gcal/ha)
	RJ0(6,i) = profperha2(i,2) / profperout2(i,2)
	RJ0(7,i) = profperha2(i,3) / profperout2(i,3)
	RJ0(8,i) = profperha2(i,4) / profperout2(i,4)
	RJ0(9,i) = profperha2(i,5) / profperout2(i,5)
		
	! Calc food grains ratio of avg observed yield to intrinsic yield
	yieldratio(i) = avgyield2(i,2) / RJ0(6,i)

	! Calc biomass avg observed yield in Gcal per ha
	avgyield2(i,1) = Bioyield(i)*GJperTON(i)/GJperGcal

	! Calc biomass intrinsic yield in Gcal per ha
	RJ0(5,i) = avgyield2(i,1) / yieldratio(i)

	! Set the intrinsic yield for biomass = that of coarse grains
!	RJ0(5,i) = RJ0(7,i)

	! Biomass profit per unit of output (Gcal)
!	profperout2(i,1) = profperha2(i,1) / RJ0(5,i)

	! Biomass average (observed) yield - NOT USED
!	avgyield2(i,1) = obsprof / profperout2(i,1)

	! Calc other cost per unit of output for 5 crops (same as variable cost)
!	GJ(5,i,2) = CALP(i,8) - profperout2(i,1) ! could use BIOPRICE instead of CALP - biomass calc no longer used
	GJ(6,i,2) = CALP(i,4) - profperout2(i,2)
	GJ(7,i,2) = CALP(i,5) - profperout2(i,3)
	GJ(8,i,2) = CALP(i,6) - profperout2(i,4)
	GJ(9,i,2) = CALP(i,7) - profperout2(i,5)

  
  ! Level 1 Final Calculations

	! Set Animal Products Intrinsic Yield = 1
	RJ0(4,i) = Animyield

	! Set pasture intrinsic yield = food grain intrinsic yield
	RJ0(3,i) = RJ0(6,i)

	! Calc profit per unit of output for pasture
	profperout1(i,2) = profperha1(i,2) / RJ0(3,i)

	! Calc pasture average (observed) yield (in Gcal/ha)
	avgyield1(i,2) = obsprof / profperout1(i,2)

	! Set other cost per unit of output of pasture and animal products
	othcostperout1(i,2) = 0.0d0
	othcostperout1(i,3) = 0.0d0
	
	! Calc the price per unit output of pasture (calibration price)
	CALP(i,10) = profperout1(i,2) + othcostperout1(i,2)

	! Calc the pasture output ratio by dividing supply of pasture by supply of beefprds 
	PASTOUT(i,2) = (avgyield1(i,2) * LAND(i,2)) / FeedOutCalc(i,1)
	
! BREAK = 1.0d0

	! Calc the calibration price of animal products
	CALP(i,3)=AnimFeed(i,5) * FEEDOUT(i,2) + CALP(i,10) * PASTOUT(i,2) + othcostperout1(i,3)

  ! Extend the variable cost data throughout the 3D Array
	DO sector = 1,numsectors
		DO M=MSTART, NZ(i)
		  GJ(sector,i,M) = GJ(sector,i,MSTART) 
		END DO  
	END DO

  END DO

! BREAK = 1.0d0

! Procedure to calculate parameters of Cobb-Douglas production
! function used to create animal products from feed and pasture.

DO i=1,NLP

  ncrops = NC(i)

  TempSum = InputAnimPrds(i,1) + InputAnimPrds(i,2) + &
            InputAnimPrds(i,3) + InputAnimPrds(i,4)

  ShareFeed(1) = InputAnimPrds(i,1)/TempSum
  ShareFeed(2) = InputAnimPrds(i,2)/TempSum
  ShareFeed(3) = InputAnimPrds(i,3)/TempSum
  ShareFeed(4) = InputAnimPrds(i,4)/TempSum

  p_crop(1) = CALP(i,4)
  p_crop(2) = CALP(i,5)
  p_crop(3) = CALP(i,6)
  p_crop(4) = CALP(i,7)

  p_feed = 0.0d0
  do k=1,ncrops
    p_feed = p_feed + ShareFeed(k) * p_crop(k)
  end do

  p_AnimPrds = CALP(i,3)
  p_Pasture  = CALP(i,10)
  
! Cobb-Douglas version

  alpha(i) = FEEDOUT(i,2) * p_feed / p_AnimPrds

  beta(i)  = PASTOUT(i,2) * p_Pasture / p_AnimPrds

  constant(i) = 1.0d0 / (FEEDOUT(i,2)**alpha(i) * PASTOUT(i,2)**beta(i))

! CES version

  alphaCES(i) = FEEDOUT(i,2) ** (1.0d0 - agrho(i)) * p_feed /  p_AnimPrds
  alphaCES(i) = alphaCES(i) ** (1.0d0 / agrho(i))
  
  betaCES(i) = PASTOUT(i,2) ** (1.0d0 - agrho(i)) * p_Pasture /  p_AnimPrds
  betaCES(i) = betaCES(i) ** (1.0d0 / agrho(i))

!   BREAK = 1.0d0

END DO 

! BREAK = 1.0d0

CALL Calout(avgyield1,avgyield2,fyield) ! Prints out calibration data

END SUBROUTINE Ag2Calibrate


! ==================================================================

SUBROUTINE ForestCal(fyield,forestreq)

! This subroutine replicates the calculations formerly done in the
! Forest2 worksheet ending the use of Excel Solver from the Input
! spreadsheet.
! This subroutine returns data into DemCoef & CALP and passes back fyield and
! forestreq to the calibration subroutine
! written 6/01 ktg

USE Ag2Global8

IMPLICIT NONE

! local variables: supply of forest prds, tech change improvement factor,
!  total 2035 Forest prds supply, total 2035 Forest prds demand, indust
!  products total demand, fuel prds total demand, future forest prds price,
!  future forest yield for 2035, total land required for managed forest
!  tolerance level of accuracy for bisection, how close the bisection is,
!  upper and lower limits for bisection
REAL(8) ForestSupply(NLP,4),techfactor, totsupply, totdem, inddem, &
		fueldem, testpr, fyield(NLP), forestreq(NLP), tolerance, &
		accuracy, upper, lower
		
! i is the region
INTEGER i

totsupply = 0.0_8

DO i=1,NLP

  ! Calculate 1990 calibration parameters - industrial wood & Fuel wood
  DemCoef(i,2)=Forestprds(i,2)/(popu(i,2)*GDPcap(i,2)**IE(i,2)*CALP(i,1)**PE(i,2))
  DemCoef(i,1)=Forestprds(i,3)/(popu(i,2)*GDPcap(i,2)**IE(i,1)*CALP(i,1)**PE(i,1))

  ! Calculate the tech change improvment
  techfactor = (1+KJIMP(2,i,3))**STEP

  ! Calculate Forest Products Supply
  ForestSupply(i,1) = Forestprds(i,1)
  ForestSupply(i,2) = ForestSupply(i,1)*techfactor
  ForestSupply(i,3) = ForestSupply(i,2)*techfactor
  ForestSupply(i,4) = ForestSupply(i,3)*techfactor

  ! Calculate 2035 Forest yield
  fyield(i) = TREEYIELD(i,2)*(1+KJIMP(2,i,3))**(STEP*3.0d0)

  ! Calculate vintage 1990 forest land required and extend array
  Forest(1,i,2)=STEP*ForestSupply(i,4)/fyield(i)
  Forest(2,i,2)=STEP*ForestSupply(i,4)/fyield(i)
  Forest(3,i,2)=STEP*ForestSupply(i,4)/fyield(i)

  ! Calculate 3x the 1990 vintage land requirement for tot land requirement
  forestreq(i) = Forest(1,i,2) * 3.0d0

  totsupply = totsupply + ForestSupply(i,4)

END DO

! Use a bisection algorithm to find the optimal 1990 fwd forest prds price

! initialize variables
testpr = CALP(1,1) ! set testprice to start out at 1990 forest prds price
DO i=1,NLP
  inddem = inddem+(DemCoef(i,2)*popu(i,5)*gdpcap(i,5)**IE(i,2)*testpr**PE(i,2))
  fueldem = fueldem+(DemCoef(i,1)*popu(i,5)*gdpcap(i,5)**IE(i,1)*testpr**PE(i,1))
END DO

totdem = inddem + fueldem
accuracy = totdem - totsupply

! Set how close we want the bisection to get us
tolerance = .000001_8  

! Set upper and lower bounds for bisection
IF (accuracy.GT.tolerance) THEN
  lower = testpr
  upper = testpr*5.0d0
ELSE IF (accuracy.LT.-tolerance) THEN
  upper = testpr
  lower = 0.0_8
END IF

! Perform bisection 
DO WHILE (accuracy.GT.tolerance .OR. accuracy.LT.-tolerance)
  
  IF (accuracy.GT.tolerance) THEN
    lower = testpr
	testpr = (upper+lower)/2.0_8
  ELSE IF (accuracy.LT.-tolerance) THEN
	upper = testpr
	testpr = (upper+lower)/2.0_8
  END IF 

  ! reinitialize industrial and fuel demands
  inddem = 0.0_8
  fueldem = 0.0_8

  DO i=1,NLP
	inddem = inddem+(DemCoef(i,2)*popu(i,5)*gdpcap(i,5)**IE(i,2)*testpr**PE(i,2))
	fueldem = fueldem+(DemCoef(i,1)*popu(i,5)*gdpcap(i,5)**IE(i,1)*testpr**PE(i,1))
  END DO
  
  totdem = inddem + fueldem  
  accuracy = totdem - totsupply

END DO

DO i=1, NLP
  CALP(i,2) = testpr
END DO

! If we are running the US-only model, calculate US fwd forest prds net exports
SELECT CASE(ICASE)

  CASE(2)

	i=1
	inddem = DemCoef(i,2)*popu(i,5)*gdpcap(i,5)**IE(i,2)*testpr**PE(i,2)
	fueldem = DemCoef(i,1)*popu(i,5)*gdpcap(i,5)**IE(i,1)*testpr**PE(i,1)
	totdem = inddem + fueldem

	NetExport(i,2) = ForestSupply(i,4)-totdem

  CASE DEFAULT

END SELECT

! BREAK = 1.0d0

END SUBROUTINE ForestCal


! ================================================================================

SUBROUTINE Calout(avgyield1,avgyield2,fyield)

USE Ag2Global8

IMPLICIT NONE

INTEGER L,YEAR

REAL(8) avgyield1(NLP,2),avgyield2(NLP,5),fyield(NLP)

OPEN(2,FILE='Calout.csv')

! Write units to output file for intrinsic yields

WRITE(2,1001) 
1001 FORMAT (',,(Gcal/ha),(Gcal/ha),(Gcal/ha),(Gcal/ha),(Gcal/ha),(Gcal/ha),&
&(Gcal/ha),(Gcal/ha),(Gcal/ha)')

! Write labels to output file

WRITE(2,1002)
1002 FORMAT('Region,Year,Unmanintriny,Forestintriny,Animintriny,Pastintriny,Biointriny,Foodintriny, &
&Coarseintriny,Oilintriny,Miscintriny')

! Write output

YEAR = 1990
DO L=1,NLP
    WRITE(2,1003) L,YEAR,RJ0(1,L),RJ0(2,L),RJ0(3,L),RJ0(4,L),RJ0(5,L),RJ0(6,L),RJ0(7,L), &
				  RJ0(8,L),RJ0(9,L)
END DO

DO L=1,2
	WRITE(2,*) " "
END DO

! Write units to output file for average yields

WRITE(2,1004) 
1004 FORMAT (',,(Gcal/ha),(Gcal/ha),(Gcal/ha),(Gcal/ha),(Gcal/ha),(Gcal/ha)')

! Write labels to output file

WRITE(2,1005)
1005 FORMAT('Region,Year,Forest_avgy,Food_avgy,Coarse_avgy,Oil_avgy,Misc_avgy,Bio_avgy')

! Write output

YEAR = 1990
DO L=1,NLP
    WRITE(2,1003) L,YEAR,avgyield1(L,1),avgyield2(L,2),avgyield2(L,3), &             
				  avgyield2(L,4),avgyield2(L,5),avgyield2(L,1)
END DO

DO L=1,2
	WRITE(2,*) " "
END DO

! Write units to output file for other calibration output

WRITE(2,1006) 
1006 FORMAT (',,($/c.m),($/Gcal),($/Gcal),($/Gcal),($/Gcal),($/Gcal),($/?), &
&($/Gcal),($/Gcal),($/c.m),($/Gcal),(c.m/ha)')

! Write labels to output file

WRITE(2,1007)
1007 FORMAT('Region,Year,vc_manfor,vc_bio,vc_foodgr,vc_coarsegr,vc_oilcr,vc_misccr,calp_unman,calp_past, &
&calp_anim,calp_ffor,calp_bio,2035fyield')

! Write output

YEAR = 1990
DO L=1,NLP
    WRITE(2,1003) L,YEAR,GJ(2,L,2),GJ(5,L,2),GJ(6,L,2),GJ(7,L,2),GJ(8,L,2), &
				  GJ(9,L,2),CALP(L,9),CALP(L,10),CALP(L,3),CALP(L,2),CALP(L,8),fyield(L)
END DO
1003 FORMAT(2(I4,','),50(F15.4,','))


CLOSE(2)

END SUBROUTINE Calout
