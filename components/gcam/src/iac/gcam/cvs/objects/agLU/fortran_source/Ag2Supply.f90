SUBROUTINE Ag2Supply(Region,t)

USE Ag2Global8

IMPLICIT NONE

INTEGER Region,t,s,i,ncrops,k

REAL(8) ModeCrop(10),ModeBio,ModePast,ModeTreeF,  &
        RateCrop(10),RateBio,RatePast,RateTreeF,RateUnMan,  &
        IntRate,Factor,sigma,sigma1,sigma2,sigma3,  &
        TempCrop(10),TempBio,ShareCrop(10),ShareBio,  &
        TempCropBio,TempPast,TempTreeF,TempUnMan,TempSum,RateAvg,RateCropBio,  &
        ShareCropBio,SharePast,ShareTreeF,ShareUnMan,  &
        YieldCrop(10),YieldBio,YieldPast,YieldTreeF, &
        SelfCons4,SelfCons5,SelfCons6,SelfCons7

REAL(8) p_woodt,p_woodf,p_crop(6),p_bio,p_past, &
        vc_crop(10),vc_bio,vc_past,vc_woodf, &
        RateMgdLand,TempMgdLand,ShareMgdLand,climate2(NZMAX,0:NLP,0:NMP)

i = Region
s = 3

ncrops = NC(i)

sigma = SIG(i)

sigma1 = sigma * SQRT(1.0d0 - CORRS(i,1))
sigma2 = sigma * SQRT(1.0d0 - CORRS(i,2))
sigma3 = sigma * SQRT(1.0d0 - CORRS(i,3))


! Set regional prices.

p_woodt    = price(1,i)
p_woodf    = price(2,i)
p_crop(1)  = price(4,i)
p_crop(2)  = price(5,i)
p_crop(3)  = price(6,i)
p_crop(4)  = price(7,i)
p_bio      = price(8,i)
p_past	   = price(10,i)

! Set variable cost parameters.

vc_woodf   = GJ(2,i,t)
vc_past    = GJ(3,i,t)
vc_bio     = GJ(5,i,t)

do k=1,ncrops
  vc_crop(k)  = GJ(5+k,i,t)
end do


! Adjust yield distributions for technical change and climate impacts.

climate2(:,i,t) = climate(:,i,t)
IF(CLIMTOGGLE(1) .EQ. 0) climate2(:,i,t) = 1.0d0	! Turn the climate impacts off

ModeTreeF   = RJ0(2,i) * KJ(2,i,t) * climate2(2,i,t)
ModePast    = RJ0(3,i) * KJ(3,i,t) * climate2(3,i,t)
ModeBio     = RJ0(5,i) * KJ(5,i,t) * climate2(5,i,t)

do k=1,ncrops
  ModeCrop(k) = RJ0(5+k,i) * KJ(5+k,i,t) * climate2(5+k,i,t)
end do


! Calculate rate of return per hectare for crop/bio nest.

RateBio = ModeBio * (p_bio - vc_bio)

do k=1,ncrops
  RateCrop(k) = ModeCrop(k) * (p_crop(k) - vc_crop(k))
end do

! Calculate Rate of Return of pasture
RatePast = ModePast * (p_past - vc_past)

! Calculate rate of return for for forests planted today.

IntRate = 0.02d0
Factor = IntRate / ((1.0d0 + IntRate)**45.0d0 - 1.0d0)

RateTreeF = Factor * ModeTreeF * (p_woodf - vc_woodf)

RateUnMan = CALP(i,9)


! Check for profit rates less than zero

RateTreeF   = max(RateTreeF,0.0001d0)
RatePast    = max(RatePast, 0.0001d0)
RateBio     = max(RateBio,0.000d0)

do k=1,ncrops
  RateCrop(k)  = max(RateCrop(k),0.0001d0)
end do

! Get land shares and average profit rate for crop/bio nest.

IF (T .EQ. 2) RateBio = 0.0_8 ! Initialize Biomass to zero base yr prod

TempBio = RateBio  ** (1/sigma3)

do k=1,ncrops
  TempCrop(k) = RateCrop(k) ** (1/sigma3)
end do

TempSum = TempBio
do k=1,ncrops
  TempSum = TempSum + TempCrop(k)
end do

RateCropBio = TempSum ** sigma3

ShareBio = TempBio  / TempSum
SaveLandShare(5,i,t) = ShareBio

do k=1,ncrops
  ShareCrop(k) = TempCrop(k) / TempSum
  SaveLandShare(5+k,i,t) = ShareCrop(k)
end do

! Get land shares and average profit rate for managed land

TempCropBio = RateCropBio ** (1/sigma2)
TempPast = RatePast ** (1/sigma2)
TempTreeF = RateTreeF ** (1/sigma2)

TempSum = TempCropBio + TempPast + TempTreeF
RateMgdLand = TempSum ** sigma2

ShareCropBio = TempCropBio / TempSum
SharePast    = TempPast    / TempSum
ShareTreeF   = TempTreeF   / TempSum

SaveLandShare(2,i,t) = ShareTreeF
SaveLandShare(3,i,t) = SharePast
SaveLandShare(4,i,t) = ShareCropBio

! Managed/Unmanaged land is the top nest

TempUnMan = RateUnMan ** (1/sigma1)
TempMgdLand = RateMgdLand ** (1/sigma1)
TempSum = TempUnMan + TempMgdLand
RateAvg = TempSum ** sigma1

ShareUnMan   = TempUnMan   / TempSum
ShareMgdLand = TempMgdLand / TempSum

SaveLandShare(1,i,t) = ShareUnMan


! Back out as-operated yield

IF (CLIMTOGGLE(1) .EQ. 0) climate2(:,i,t) = climate(:,i,t) ! Reset climate impacts
IF (CLIMTOGGLE(2) .EQ. 0) climate2(:,i,t) = 1.0d0	! Turn off climate impacts

YieldTreeF = (RateAvg / ((p_woodf - vc_woodf) * Factor))*climate2(2,i,T)
YieldTreeF = max(YieldTreeF,0.0d0)
YieldPast = (RateAvg / (p_past - vc_past))*climate2(3,i,T)
YieldPast = max(YieldPast, 0.0d0)
YieldBio = (RateAvg / (p_bio - vc_bio))*climate2(5,i,T)
YieldBio = max(YieldBio,0.000d0) ! Don't let yields go below zero
YieldBio = min(YieldBio,bioyieldcap(i)) ! Don't allow yields above unreasonable levels

do k=1,ncrops
  YieldCrop(k) = (RateAvg / (p_crop(k) - vc_crop(k)))*climate2(5+k,i,T)
  YieldCrop(k) = max(YieldCrop(k), 0.0d0) ! if the price is less than vc, no crops & no yield
  SaveYield(5+k,i,t) = YieldCrop(k)
end do

! Save tree yield for supply calculation when trees are cut
TreeYield(i,T+S) = YieldTreeF

! new amount of forest land
ForestN(i,T) = LAND(i,6) * ShareMgdLand * ShareTreeF - Forest(1,i,T) - Forest(2,i,T)
ForestN(i,T) = max(ForestN(i,T),0.0d0) ! Next vintage of forest land can't be less than zero

! Calculate supply for Forests, Fwd Forests and Pasture
qsup(1,i,T) = Forest(3,i,t) * TreeYield(i,T) / STEP
qsup(2,i,T) = ForestN(i,t) * YieldTreeF / STEP 
qsup(3,i,T) = LAND(i,6) * ShareMgdLand * SharePast * YieldPast / 10.0d0

do k=1,ncrops
  qsup(3+k,i,t) = LAND(i,6) * ShareMgdLand * ShareCropBio * ShareCrop(k) * YieldCrop(k) / 10.0d0
end do

! supply of biomass crops
qsup(8,i,t) = LAND(i,6) * ShareMgdLand * ShareCropBio * ShareBio * YieldBio / 10.0d0

AgP(1,i,t) = qsup(1,i,t)
AgP(2,i,t) = qsup(2,i,t)
AgP(3,i,t) = qsup(3,i,t)
AgP(4,i,t) = qsup(4,i,t)
AgP(5,i,t) = qsup(5,i,t)
AgP(6,i,t) = qsup(6,i,t)
AgP(7,i,t) = qsup(7,i,t)
AgP(8,i,t) = qsup(8,i,t)

AgP(12,i,t) = LAND(i,6) * ShareMgdLand * SharePast * YieldPast / 10.0d0
AgZ(12,3,i,t) = AgP(12,i,t)

! Adjust supply of crops for self consumption (beefprds done in demand.f90)

SelfCons4 = qsup(4,i,T) * SelfCons(i,1)
SelfCons5 = qsup(5,i,T) * SelfCons(i,2)
SelfCons6 = qsup(6,i,T) * SelfCons(i,3)
SelfCons7 = qsup(7,i,T) * SelfCons(i,4)

qsup(4,i,T) = qsup(4,i,T) - SelfCons4
qsup(5,i,T) = qsup(5,i,T) - SelfCons5
qsup(6,i,T) = qsup(6,i,T) - SelfCons6
qsup(7,i,T) = qsup(7,i,T) - SelfCons7

AgZ(4,4,i,t) = SelfCons4
AgZ(5,5,i,t) = SelfCons5
AgZ(6,6,i,t) = SelfCons6
AgZ(7,7,i,t) = SelfCons7


! Save selected output

SaveLand(1,i,T) = 0.0d0
do k=1,ncrops
  SaveCropLand(k,i,T) = LAND(i,6) * ShareMgdLand * ShareCropBio * ShareCrop(k) ! rs - 11/01
  SaveLand(1,i,T) = SaveLand(1,i,T) + SaveCropLand(k,i,T)
end do

SaveLand(2,i,T) = LAND(i,6) * ShareMgdLand * SharePast
SaveLand(3,i,T) = LAND(i,6) * ShareMgdLand * ShareTreeF
SaveLand(4,i,T) = LAND(i,6) * ShareMgdLand * ShareCropBio * ShareBio
SaveLand(5,i,T) = LAND(i,6) * ShareUnMan

SaveProf(i,T) = RateAvg

! Calculate net exports.

AgFD(1,2,i,t) = AgP(1,i,t) - AgFD(1,1,i,t) - sum(AgZ(1,:,i,t))
AgFD(2,2,i,t) = AgP(2,i,t) - AgFD(2,1,i,t) - sum(AgZ(2,:,i,t))
AgFD(3,2,i,t) = NetExport(i,8)
AgFD(4,2,i,t) = AgP(4,i,t) - AgFD(4,1,i,t) - sum(AgZ(4,:,i,t))
AgFD(5,2,i,t) = AgP(5,i,t) - AgFD(5,1,i,t) - sum(AgZ(5,:,i,t))
AgFD(6,2,i,t) = AgP(6,i,t) - AgFD(6,1,i,t) - sum(AgZ(6,:,i,t))
AgFD(7,2,i,t) = AgP(7,i,t) - AgFD(7,1,i,t) - sum(AgZ(7,:,i,t))
AgFD(8,2,i,t) = AgP(8,i,t) - AgFD(8,1,i,t) - sum(AgZ(8,:,i,t))

! Calculate output for other sectors.

AgP(9,i,t)  = AgFD(9,1,i,t)  + AgFD(9,2,i,t)  + sum(AgZ(9,:,i,t))
AgP(10,i,t) = AgFD(10,1,i,t) + AgFD(10,2,i,t) + sum(AgZ(10,:,i,t))
AgP(11,i,t) = AgFD(11,1,i,t) + AgFD(11,2,i,t) + sum(AgZ(11,:,i,t))

! Put variable cost into intermediate flows array, row 13, and view these
! payments as an import to the agricultural system.

AgZ(13,2,i,t)  = AgP(2,i,t)  * vc_woodf * Factor
AgZ(13,4,i,t)  = AgP(4,i,t)  * vc_crop(1)
AgZ(13,5,i,t)  = AgP(5,i,t)  * vc_crop(2)
AgZ(13,6,i,t)  = AgP(6,i,t)  * vc_crop(3)
AgZ(13,7,i,t)  = AgP(7,i,t)  * vc_crop(4)
AgZ(13,8,i,t)  = AgP(8,i,t)  * vc_bio
AgZ(13,12,i,t) = AgP(12,i,t) * vc_past

AgFD(13,2,i,t) = -sum(AgZ(13,:,i,t))

! Calculate value of output at base-year prices.

AgPrice(1,i,t) = p_woodt
AgPrice(2,i,t) = p_woodf * Factor
AgPrice(3,i,t) = price(3,i)
AgPrice(4,i,t) = p_crop(1)
AgPrice(5,i,t) = p_crop(2)
AgPrice(6,i,t) = p_crop(3)
AgPrice(7,i,t) = p_crop(4)
AgPrice(8,i,t) = p_bio
AgPrice(12,i,t) = p_past

END SUBROUTINE Ag2Supply
