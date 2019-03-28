SUBROUTINE AG2OUTPUT

! Procedure to write output to a text file.

USE Ag2Global8

IMPLICIT NONE

INTEGER L,M,YEAR,J,i

REAL(8) WorldTot(9,15)


OPEN(2,FILE='AgOut.csv')

! Write units to output file

WRITE(2,1001) 
1001 FORMAT (',,(10^12 $),(x000),(000 ha),(000 ha),(000 ha),(000 ha),(000 ha),(000 ha),&
&(TgC),($/c.m.),($/c.m.),($/c.m.),($/Gcal),($/Gcal),($/Gcal),($/Gcal),($/Gcal),($/Gcal),($/Gcal),,,&
&(kcal/cap),(kcal/cap),(kcal/cap),(kcal/cap),(kcal/cap),(kcal/cap),(kcal/cap),(kcal/cap),(TgC),(TgC),(TgC),(TgN),(TgN),(TgN), &
&(10^13cal),(10^13cal),(10^13cal),,(GCal/ha),(GCal/ha),(GCal/ha),(GCal/ha),&
&(000 ha),(000 ha),(000 ha),(000 ha)')


! (10^10 kcal),(10^10 kcal),&
! (10^10 kcal),(10^10 kcal),(000 c.m.),(10^10 kcal),(10^10 kcal),(10^10 kcal),&
! ,,&
! ($/Gcal),($/Gcal),($/Gcal),,,(kcal/person/day)')


! Write labels to output file

WRITE(2,1002)
1002 FORMAT('Region,Year,GDP,POP,&
&CropLand,PastureLand,ForestLand,BioLand,UnManLand,OthLand,LUCEmiss,&
&TreePrice,CalTFPrice,TreeFPrice,BeefPrice,p_FoodGr,p_CoarseGr,p_OilCrop,p_MiscCrop,&
&p_bio,p_past,FEEDOUT,PASTOUT,BeefPrds,FoodGr,CoarseGr,OilCrop,MiscCrop,ProcCrop,Pork,Poultry,&
&enterCH4,wasteCH4,riceCH4,fertilN2O,manwastN2O,unmanwastN2O,BioSupply,PastSupply,PastDem,Fwoodprd,&
&YieldCrop1,YieldCrop2,YieldCrop3,YieldCrop4,CropLand1,CropLand2,CropLand3,CropLand4')


! ForestProd,AnimProd,FieldCrops,OthCrops,BiomassProd,&
! ForestDem,AnimDem,FCropDem,OthCropDem,,&
! ,,,FCropPrice,OthCropPrice,BioPrice,&
! FEEDOUT,PASTOUT,FCrops,OthCrops,ProcCrop,AnimPrds,OthCrop,FCrop1,FCrop2,FCrop3,FCrop4,FCrop5,&
! UnMan,Forest,Pasture,CropBio,Bio,OthCrop,FCrop1,FCrop2,FCrop3,FCrop4,FCrop5')

WorldTot = 0.0d0
DO L=1, NLP
  WRITE(2,*)
  DO M=2,NMP
    YEAR = 1960 + M*15

    WRITE(2,1003) L,YEAR,GDP(L,M),popu(L,M),SaveLand(1,L,M),SaveLand(2,L,M),SaveLand(3,L,M), &
	              SaveLand(4,L,M),SaveLand(5,L,M),LAND(L,4),CarbEmiss(L,M),sprice(1,L,M),CALP(L,2),sprice(2,L,M),&
                  sprice(3,L,M),sprice(4,L,M),sprice(5,L,M),sprice(6,L,M),sprice(7,L,M), &
                  sprice(8,L,M),sprice(10,L,M),FEEDOUT(L,M),PASTOUT(L,M),SaveDiet(3,L,M), &
                  SaveDiet(4,L,M),SaveDiet(5,L,M),SaveDiet(6,L,M),SaveDiet(7,L,M),SaveDiet(8,L,M), &
				  SaveDiet(9,L,M),SaveDiet(10,L,M),Ag2CH4(1,L,M),Ag2CH4(2,L,M),Ag2CH4(3,L,M), &
				  Ag2N2O(1,L,M),Ag2N2O(2,L,M),Ag2N2O(3,L,M),qsup(8,L,M),qsup(3,L,M),qdem(3,L,M),qsup(2,L,M), &
                  SaveYield(6,L,M),SaveYield(7,L,M),SaveYield(8,L,M),SaveYield(9,L,M), &
                  SaveCropLand(1,L,M),SaveCropLand(2,L,M),SaveCropLand(3,L,M),SaveCropLand(4,L,M)
	
	
! qsup(1,L,M),qsup(3,L,M),qsup(4,L,M), &
! qsup(5,L,M),qsup(6,L,M),qdem(1,L,M),qdem(3,L,M),qdem(4,L,M), &
! qdem(5,L,M),, &
! ,,SaveProF(L,M), &
! ,sprice(5,L,M), &
! sprice(6,L,M),FEEDOUT(L,M),PASTOUT(L,M),SaveDiet(1,L,M),SaveDiet(2,L,M), &
! SaveDiet(3,L,M),SaveDiet(4,L,M),SaveYield(6,L,M),SaveYield(7,L,M), &
! SaveYield(8,L,M),SaveYield(9,L,M),SaveYield(10,L,M),SaveYield(11,L,M), &
! SaveLandShare(1,L,M),SaveLandShare(2,L,M),SaveLandShare(3,L,M), &
! SaveLandShare(4,L,M),SaveLandShare(5,L,M),SaveLandShare(6,L,M), &
! SaveLandShare(7,L,M),SaveLandShare(8,L,M),SaveLandShare(9,L,M), &
! SaveLandShare(10,L,M),SaveLandShare(11,L,M),SaveLandShare(12,L,M)

    WorldTot(M,1) = WorldTot(M,1) + GDP(L,M)
    WorldTot(M,2) = WorldTot(M,2) + popu(L,M)
    WorldTot(M,3) = WorldTot(M,3) + SaveLand(1,L,M)
    WorldTot(M,4) = WorldTot(M,4) + SaveLand(2,L,M)
    WorldTot(M,5) = WorldTot(M,5) + SaveLand(3,L,M)
    WorldTot(M,6) = WorldTot(M,6) + SaveLand(4,L,M)
    WorldTot(M,7) = WorldTot(M,7) + SaveLand(5,L,M)
	WorldTot(M,8) = WorldTot(M,8) + LAND(L,4)
    WorldTot(M,9) = WorldTot(M,9) + CarbEmiss(L,M)

  END DO
END DO

L=99
WRITE(2,*)
DO M=2,NMP
  YEAR = 1960 + M*15
  WRITE(2,1003) L,YEAR,(WorldTot(M,J),J=1,9)
END DO

1003 FORMAT(2(I4,','),50(F15.4,','))

CLOSE(2)


OPEN(3,FILE='AgBal.csv')

WRITE(3,*) 'Agricultural product balances for United States'

DO M=2,NMP
  WRITE(3,*)
  WRITE(3,*)
  WRITE(3,*)
  WRITE(3,*) 'M =',M
  WRITE(3,1004)
  WRITE(3,1005)
  DO i=1,13
    WRITE(3,1003) 1,i,(AgZ(i,j,1,M),j=1,13),(AgFD(i,j,1,M),j=1,2),AgP(i,1,M),AgPrice(i,1,M)
  END DO
END DO

1004 FORMAT(',,1,2,3,4,5,6,7,8,9,10,11,12,13')

1005 FORMAT('Region,Input,ForestPrds,FForest,Beef,&
&FoodGrains,CoarseGrains,OilCrops,OtherAg,Biomass,&
&ProcCrop,Pork,Poultry,Pasture,OthVarCost,Consumption,NetExports,Production,AgPrice')

CLOSE(3)

OPEN(4,FILE='AgEmiss.csv')

WRITE(4,*) 'Land Use Change CO2 Emissions (1975 numbers are missing residual soil emissions from prior time periods)'

DO M=1,NMP
  WRITE(4,*)
  WRITE(4,*)
  WRITE(4,*) 'M=',M
  WRITE(4,1006)
  DO L=1,NLP
	YEAR = 1960 + M*15
	WRITE(4,1003) L,YEAR,CarbEmiss(L,M),SoilEmiss(1,L,M),SoilEmiss(2,L,M), &
		SoilEmiss(3,L,M),SoilEmiss(4,L,M)
  END DO
END DO

1006 FORMAT('Region,Year,CarbEmiss,TotSoilEmiss,ThisPerSoilEmiss,NextPerSoilEmiss,FinPerSoilEmiss')

CLOSE(4)

END SUBROUTINE AG2OUTPUT
