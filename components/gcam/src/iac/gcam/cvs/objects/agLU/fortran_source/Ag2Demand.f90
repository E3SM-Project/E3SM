SUBROUTINE Ag2Demand(Region,t)

USE Ag2Global8

IMPLICIT NONE

INTEGER Region,t,s,i,jwf,jwp,mkt,k

REAL(8) WoodFuel,WoodPrds,BeefMutton1,BeefMutton2,FoodGrain1, &
        FoodGrain2,CoarseGrain1,CoarseGrain2,OilCrop1,OilCrop2,MiscCrop1, &
        MiscCrop2,ProcCrop1,ProcCrop2,p_ProcCrop,p0_ProcCrop(20),TempSum,ProcCrop, &
		p_Pork,p0_Pork(20),Pork1,Pork2,Pork,p_poultry,p0_poultry(20),poultry1, &
		poultry2,poultry,ShareFeed(6),FeedCrop,beefdem

! Ordering of regional markets
!
!  1  Wood (current)
!  2  Wood (forward)
!  3  Pasture [Beef/Mutton]
!  4  Food Grains (wheat and rice)
!  5  Coarse Grains (all other cereals)
!  6  Oil Crops
!  7  Miscellaneous Crops (fruits, vegetabels, starchy roots, pulses)
!  8  Commercial Biomass
! [9  Processed Crops]
! [10 Pork]
! [11 Poultry]


i = Region

CALL BEEFPRICE(i,T,ShareFeed) ! Get price of beef prds	

jwf = 1       ! wood for fuel
jwp = 2       ! wood products

s = 3


! Present demand for forest products

WoodFuel = popu(i,t) * DemCoef(i,jwf) * gdpcap(i,t)**IE(i,jwf) * &
           price(1,i)**PE(i,jwf)

WoodPrds = popu(i,t) * DemCoef(i,jwp) * gdpcap(i,t)**IE(i,jwp) * &
           price(1,i)**PE(i,jwp)

qdem(1,i,t) = WoodFuel + WoodPrds
AgFD(1,1,i,t) = qdem(1,i,t)


! Forward demand for forest products

WoodFuel = popu(i,t+s) * DemCoef(i,jwf) * gdpcap(i,t+s)**IE(i,jwf) * &
           price(2,i)**PE(i,jwf)
WoodPrds = popu(i,t+s) * DemCoef(i,jwp) * gdpcap(i,t+s)**IE(i,jwp) * &
           price(2,i)**PE(i,jwp)
qdem(2,i,t) = WoodFuel + WoodPrds
AgFD(2,1,i,t) = qdem(2,i,t)

! For all sectors below, there are two demands.
! For example:
! BeefMutton1 = demand for food
! BeefMutton2 = demand for other non-food uses


! Demand for Beef and Mutton products with price feedback
! Now derive demand for pasture as well
mkt = 3

BeefMutton1 = popu(i,t) * kcal1(i,6,t) * (price(mkt,i)/CALP(i,mkt))**PE(i,mkt) &
           * (gdpcap(i,t)/gdpcap(i,2))**IE(i,mkt) * 365.0d0 / 10000000.0d0

BeefMutton2 = popu(i,t) * kcal2(i,6) * (price(mkt,i)/CALP(i,mkt))**PE(i,mkt) &
           * (gdpcap(i,t)/gdpcap(i,2))**IE(i,mkt) * 365.0d0 / 10000000.0d0

SaveDiet(mkt,i,t) = BeefMutton1 * 10000000.0d0 / (popu(i,t) * 365.0d0)

  ! Adjust for fixed trade in beefprds demand
beefdem = BeefMutton1 + BeefMutton2 + NetExport(i,8)


  ! Calculate gross beef production with self-consumption accounted for
beefdem = beefdem / (1.0d0-SelfCons(i,6))

AgZ(3,3,i,t) = beefdem*SelfCons(i,6)	! Store self consumption
AgFD(mkt,1,i,t) = BeefMutton1+BeefMutton2+(AgZ(3,3,i,t)) ! Store beef dem w/o trade

qdem(mkt,i,t) = beefdem*PASTOUT(i,T) ! Derive Pasture Demand


! Demand for crop #1 with price feedback
mkt = 4

FoodGrain1 = popu(i,t) * kcal1(i,1,t) * (price(mkt,i)/CALP(i,mkt))**PE(i,mkt) &
             * (gdpcap(i,t)/gdpcap(i,2))**IE(i,mkt) * 365.0d0 / 10000000.0d0

FoodGrain2 = popu(i,t) * kcal2(i,1) * (price(mkt,i)/CALP(i,mkt))**PE(i,mkt) &
             * (gdpcap(i,t)/gdpcap(i,2))**IE(i,mkt) * 365.0d0 / 10000000.0d0

SaveDiet(mkt,i,t) = FoodGrain1 * 10000000.0d0 / (popu(i,t) * 365.0d0)

qdem(mkt,i,t) = FoodGrain1 + FoodGrain2
AgFD(mkt,1,i,t) = qdem(mkt,i,t)


! Demand for crop #2 with price feedback
mkt = 5

CoarseGrain1 = popu(i,t) * kcal1(i,2,t) * (price(mkt,i)/CALP(i,mkt))**PE(i,mkt) &
               * (gdpcap(i,t)/gdpcap(i,2))**IE(i,mkt) * 365.0d0 / 10000000.0d0

CoarseGrain2 = popu(i,t) * kcal2(i,2) * (price(mkt,i)/CALP(i,mkt))**PE(i,mkt) &
               * (gdpcap(i,t)/gdpcap(i,2))**IE(i,mkt) * 365.0d0 / 10000000.0d0

SaveDiet(mkt,i,t) = CoarseGrain1 * 10000000.0d0 / (popu(i,t) * 365.0d0)

qdem(mkt,i,t) = CoarseGrain1 + CoarseGrain2
AgFD(mkt,1,i,t) = qdem(mkt,i,t)


! Demand for crop #3 with price feedback
mkt = 6

OilCrop1 = popu(i,t) * kcal1(i,3,t) * (price(mkt,i)/CALP(i,mkt))**PE(i,mkt) &
           * (gdpcap(i,t)/gdpcap(i,2))**IE(i,mkt) * 365.0d0 / 10000000.0d0

OilCrop2 = popu(i,t) * kcal2(i,3) * (price(mkt,i)/CALP(i,mkt))**PE(i,mkt) &
           * (gdpcap(i,t)/gdpcap(i,2))**IE(i,mkt) * 365.0d0 / 10000000.0d0

SaveDiet(mkt,i,t) = OilCrop1 * 10000000.0d0 / (popu(i,t) * 365.0d0)

qdem(mkt,i,t) = OilCrop1 + OilCrop2
AgFD(mkt,1,i,t) = qdem(mkt,i,t)


! Demand for crop #4 with price feedback
mkt = 7

MiscCrop1 = popu(i,t) * kcal1(i,4,t) * (price(mkt,i)/CALP(i,mkt))**PE(i,mkt) &
            * (gdpcap(i,t)/gdpcap(i,2))**IE(i,mkt) * 365.0d0 / 10000000.0d0

MiscCrop2 = popu(i,t) * kcal2(i,4) * (price(mkt,i)/CALP(i,mkt))**PE(i,mkt) &
            * (gdpcap(i,t)/gdpcap(i,2))**IE(i,mkt) * 365.0d0 / 10000000.0d0

SaveDiet(mkt,i,t) = MiscCrop1 * 10000000.0d0 / (popu(i,t) * 365.0d0)

qdem(mkt,i,t) = MiscCrop1 + MiscCrop2
AgFD(mkt,1,i,t) = qdem(mkt,i,t)


! Demand for processed crops with price feedback
! NEW -- elasticities read in. sjs. 11/01

mkt = 10

! Get a price for processed crops

TempSum = 0.0
p_ProcCrop = 0.0d0

do k=1,4
  p_ProcCrop = p_ProcCrop + price(3+k,i) * InputProcCrop(i,k)
  ! TempSum = TempSum + InputProcCrop(i,k)
end do
! p_ProcCrop = p_ProcCrop / TempSum


if (t .eq. 2)  p0_ProcCrop(i) = p_ProcCrop

ProcCrop1 = popu(i,t) * kcal1(i,5,t) * (p_ProcCrop/p0_ProcCrop(i))**PE(i,mkt) &
            * (gdpcap(i,t)/gdpcap(i,2))**IE(i,mkt) * 365.0d0 / 10000000.0d0

ProcCrop2 = popu(i,t) * kcal2(i,5) * (p_ProcCrop/p0_ProcCrop(i))**PE(i,mkt) &
            * (gdpcap(i,t)/gdpcap(i,2))**IE(i,mkt) * 365.0d0 / 10000000.0d0

SaveDiet(8,i,t) = ProcCrop1 * 10000000.0d0 / (popu(i,t) * 365.0d0)

! Include derived demands for processed crops, adjusting for net exports

ProcCrop = ProcCrop1 + ProcCrop2 + NetExport(i,7)
AgFD(9,1,i,t) = ProcCrop1 + ProcCrop2
AgFD(9,2,i,t) = NetExport(i,7)

qdem(4,i,t) = qdem(4,i,t) + ProcCrop * InputProcCrop(i,1)
qdem(5,i,t) = qdem(5,i,t) + ProcCrop * InputProcCrop(i,2)
qdem(6,i,t) = qdem(6,i,t) + ProcCrop * InputProcCrop(i,3)
qdem(7,i,t) = qdem(7,i,t) + ProcCrop * InputProcCrop(i,4)

AgZ(4,9,i,t) = ProcCrop * InputProcCrop(i,1)
AgZ(5,9,i,t) = ProcCrop * InputProcCrop(i,2)
AgZ(6,9,i,t) = ProcCrop * InputProcCrop(i,3)
AgZ(7,9,i,t) = ProcCrop * InputProcCrop(i,4)


! Get a price for PorkPrds

!TempSum = 0.0
p_Pork = 0.0d0

do k=1,4
  p_Pork = p_Pork + price(3+k,i) * InputAnimPrds(i,k+4)
  ! TempSum = TempSum + InputAnimPrds(i,k+4)
end do
! p_Pork = p_Pork / TempSum

if (t .eq. 2)  p0_Pork(i) = p_Pork

Pork1 = popu(i,t) * kcal1(i,7,t) * (p_Pork/p0_Pork(i))**PE(i,8) &
       * (gdpcap(i,t)/gdpcap(i,2))**IE(i,8) * 365.0d0 / 10000000.0d0

Pork2 = popu(i,t) * kcal2(i,7) * (p_Pork/p0_Pork(i))**PE(i,8) &
       * (gdpcap(i,t)/gdpcap(i,2))**IE(i,8) * 365.0d0 / 10000000.0d0

SaveDiet(9,i,t) = Pork1 * 10000000.0d0 / (popu(i,t) * 365.0d0)

Porkprod(i,T) = Pork1 + Pork2

! Include derived demands for pork, adjusting for net exports

Pork = Pork1 + Pork2 + NetExport(i,9)

qdem(4,i,t) = qdem(4,i,t) + Pork * InputAnimPrds(i,5)
qdem(5,i,t) = qdem(5,i,t) + Pork * InputAnimPrds(i,6)
qdem(6,i,t) = qdem(6,i,t) + Pork * InputAnimPrds(i,7)
qdem(7,i,t) = qdem(7,i,t) + Pork * InputAnimPrds(i,8)

AgFD(10,1,i,t) = Pork1 + Pork2
AgFD(10,2,i,t) = NetExport(i,9)

AgZ(4,10,i,t) = Pork * InputAnimPrds(i,5)
AgZ(5,10,i,t) = Pork * InputAnimPrds(i,6)
AgZ(6,10,i,t) = Pork * InputAnimPrds(i,7)
AgZ(7,10,i,t) = Pork * InputAnimPrds(i,8)


! Get a price for Poultry

!TempSum = 0.0
p_poultry = 0.0d0

do k=1,4
  p_poultry = p_poultry + price(3+k,i) * InputAnimPrds(i,k+8)
  ! TempSum = TempSum + InputAnimPrds(i,k+8)
end do
! p_poultry = p_poultry / TempSum

if (t .eq. 2)  p0_poultry(i) = p_poultry

poultry1 = popu(i,t) * kcal1(i,8,t) * (p_poultry/p0_poultry(i))**PE(i,9) &
         * (gdpcap(i,t)/gdpcap(i,2))**IE(i,9) * 365.0d0 / 10000000.0d0

poultry2 = popu(i,t) * kcal2(i,8) * (p_poultry/p0_poultry(i))**PE(i,9) &
         * (gdpcap(i,t)/gdpcap(i,2))**IE(i,9) * 365.0d0 / 10000000.0d0

SaveDiet(10,i,t) = poultry1 * 10000000.0d0 / (popu(i,t) * 365.0d0)

Poultryprod(i,T) = poultry1 + poultry2

! Include derived demands for poultry, adjusting for net exports

poultry = poultry1 + poultry2 + NetExport(i,10)

qdem(4,i,t) = qdem(4,i,t) + poultry * InputAnimPrds(i,9)
qdem(5,i,t) = qdem(5,i,t) + poultry * InputAnimPrds(i,10)
qdem(6,i,t) = qdem(6,i,t) + poultry * InputAnimPrds(i,11)
qdem(7,i,t) = qdem(7,i,t) + poultry * InputAnimPrds(i,12)

AgFD(11,1,i,t) = Poultry1 + Poultry2
AgFD(11,2,i,t) = NetExport(i,10)

AgZ(4,11,i,t) = Poultry * InputAnimPrds(i,9)
AgZ(5,11,i,t) = Poultry * InputAnimPrds(i,10)
AgZ(6,11,i,t) = Poultry * InputAnimPrds(i,11)
AgZ(7,11,i,t) = Poultry * InputAnimPrds(i,12)

AgPrice(9,i,t)  = p_ProcCrop
AgPrice(10,i,t) = p_pork
AgPrice(11,i,t) = p_poultry

! Adjust for crops fed to beef/mutton products
FeedCrop = beefdem * FEEDOUT(i,T)
DO k=1,NC(i)
  qdem(3+k,i,T) = qdem(3+k,i,T) + FeedCrop * ShareFeed(k)
  AgZ(3+k,3,i,t) = FeedCrop * ShareFeed(k)
END DO

END SUBROUTINE Ag2Demand


! ========================================================================

SUBROUTINE BEEFPRICE(i,T,ShareFeed)

! Subroutine to calculate the price of beefproducts using a pasture-feed
! substitution curve (either Cobb-Douglas or CES)

USE Ag2Global8

IMPLICIT NONE

INTEGER i, k, T

REAL(8) p_past,p_beef,ShareFeed(6),p_feed, TempSum, r


p_past = price(10,i)

! We are using a Cobb-Douglas production function for beef/mutton
! products to allow substitution between feed and pasture.
! Therefore, FEEDOUT(i) and PASTOUT(i) are both endogenous.

! First calculate shares of animal feed coming from each type of crop

TempSum = InputAnimPrds(i,1) + InputAnimPrds(i,2) + &
          InputAnimPrds(i,3) + InputAnimPrds(i,4)

ShareFeed(1) = InputAnimPrds(i,1)/TempSum
ShareFeed(2) = InputAnimPrds(i,2)/TempSum
ShareFeed(3) = InputAnimPrds(i,3)/TempSum
ShareFeed(4) = InputAnimPrds(i,4)/TempSum

p_feed = 0.0d0
do k=1,NC(i)
  p_feed = p_feed + ShareFeed(k) * price(3+k,i)
end do


! Use price equation to get p_beef

! Cobb-Douglas version
IF (agrho(i) .LT. .001d0) THEN

  p_beef = (1/constant(i)) * (p_feed/alpha(i))**alpha(i) * (p_past/beta(i))**beta(i)

  FEEDOUT(i,T) = p_beef * alpha(i) / p_feed

  PASTOUT(i,T) = p_beef * beta(i) / p_past

! CES version
ELSE

  r = agrho(i)/(agrho(i) - 1.0d0)

  p_beef = ( (p_feed/alphaCES(i))**r + (p_past/betaCES(i))**r )**(1.0d0/r)

  FEEDOUT(i,T) = alphaCES(i) ** agrho(i) * p_beef / p_feed
  FEEDOUT(i,T) = FEEDOUT(i,T) ** (1.0d0 / (1.0d0 - agrho(i)))

  PASTOUT(i,T) = betaCES(i) ** agrho(i) * p_beef / p_past
  PASTOUT(i,T) = PASTOUT(i,T) ** (1.0d0 / (1.0d0 - agrho(i)))

END IF

price(3,i) = p_beef

END SUBROUTINE BEEFPRICE