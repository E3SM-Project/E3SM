SUBROUTINE Ag2GetData(AGFILE,INF)

! Procedure to read data and store in global arrays.
! AGFILE and INF are passed in to notify if an input file is missing

use Ag2Global8

implicit none

! Arguments
CHARACTER(len=72), intent( in ) :: AGFILE( MAXFILES )
integer, intent( in ) :: INF

! Variables
CHARACTER(len=80) :: MsgStr
CHARACTER(len=11) :: TABLEMRK
INTEGER IVARNUM,L,J,INDIC,NJ,MSTART,M,IVMAX

DATA TABLEMRK/'INPUT_TABLE'/

IVMAX = 3
NJ = 4

CALL Ag2NextTable(TABLEMRK,INDIC)

! If one of the input files can not be found
IF(INDIC.EQ.0) THEN
    WRITE(*,*) "Missing input file."
END IF

DO WHILE(INDIC.EQ.1)

  READ(1,*)
! Read the variable indicator label and select the right code.

  READ(1,'(I3)')IVARNUM 
      
  READ(1,*)
  READ(1,*)

  SELECT CASE(IVARNUM)

! Read the number of regions
!  CASE(601)
!     READ (1,*) NLP

! Read population data (thousands of people)
  CASE(506)
    DO L=1,NLP
      READ (1,*) IDUM,(ZLM(L,M),M=1,NMP+1) 
    END DO

! Read base-year GDP
!  CASE(507) 
!    DO L=1,NLP
!      READ (1,*) IDUM,GNPBL(L) 
!    END DO

! Read GDP Data
!  CASE(509)
!    DO L=1,NLP
!      READ (1,*) IDUM,(GDP(L,M),M=2,NMP)
!    END DO

! Read labor productivity growth rates - now GDP is read directly
!  CASE(508)
!    DO L=1,NLP
!      READ (1,*) IDUM,(PROLM(L,M),M=1,NMP)
!    END DO
         
! Technical change parameters
  CASE (602)
    READ(1,*) MSTART,L,(KJIMP(J,L,MSTART),J=1,NZ(L))
      DO WHILE (L.GT.0)
        DO M=MSTART,NMP+3
          DO J=1,NZ(L)
            KJIMP(J,L,M) = KJIMP(J,L,MSTART)
          END DO
        END DO
        READ(1,*) MSTART,L,(KJIMP(J,L,MSTART),J=1,NZ(L))
      END DO

    KJ = 1.0d0

    DO M=2,NMP+3
      DO L=1,NLP
        DO J=1,NZ(L)
          KJ(J,L,M) = KJ(J,L,M-1)*(1.0+KJIMP(J,L,M))**STEP
        END DO
      END DO
    END DO

! Calibration prices
  CASE (603)
    DO L=1,NLP
      READ(1,*) IDUM,(CALP(L,J),J=1,10)
    END DO
        
! Feed-output ratio
  CASE (606)
    DO L=1,NLP
      READ(1,*) FEEDOUT(L,2)
    END DO
                
! Demand function coefficients                  
  CASE (608)
    DO L=1,NLP
      READ(1,*) IDUM,(DemCoef(L,J),J=1,7)
    END DO

! Income elasticity of demand
  CASE (609)
    DO L=1,NLP
      READ(1,*) IDUM,(IE(L,J),J=1,10)
    END DO
             
! Price elasticity of demand                 
  CASE (610)
    DO L=1,NLP
      READ(1,*) IDUM,(PE(L,J),J=1,10)
    END DO

! Biomass Technical Parameters                 
  CASE (612)
    DO L=1,NLP
      READ(1,*) IDUM,GJ(5,L,2),Bioyield(L),GJperTON(L),RJ0(5,L),bioyieldcap(L)
    END DO

! Cost of transforming biomass to liquid fuel ($ per GJ)         
  CASE (613)
    READ(1,*) MSTART,(BioTranCost(L,MSTART), L=1,NLP) 
    DO WHILE (MSTART.GT.0)
      DO L=1,NLP
        DO M=MSTART,NMP
          BioTranCost(L,M) = BioTranCost(L,MSTART)
        END DO
      END DO
      READ(1,*) MSTART,(BioTranCost(L,MSTART), L=1,NLP)
    END DO

! Forest land in place in period 2 - now calculated in calibration.f90
!  CASE (614)
!    DO IVIN = 1,IVMAX
!      READ(1,*) (Forest(IVIN,L,2), L=1,NLP)
!    END DO
     
! Land carbon density parameters                
  CASE (615)
    DO L=1,NLP
      READ(1,*) IDUM,(CDensity(L,J),J=1,10)
    END DO
         
! Historical Land Use (for LUC emissions calculation)                
  CASE (616)
    DO L=1,NLP
      READ(1,*) IDUM,(HistLand(L,J),J=1,10)
    END DO
         
! 1990 land use carbon emissions
!  CASE (617)
!    DO L=1,NLP
!      READ(1,*) IDUM,CarbEmiss(L,2)
!    END DO
         
! Read demand parameters for crops
!  CASE (617)
!    DO L=1,NLP
!      READ(1,*) IDUM,(DemCrop(L,J),J=1,4)
!    END DO
         
! Read demand parameters for animal products
!  CASE (618)
!    DO L=1,NLP
!      READ(1,*) IDUM,(DemAnim(L,J),J=1,4)
!    END DO
                           
! Base-year tree yield
  CASE (621)
    DO L=1,NLP
      READ(1,*) TreeYield(L,2)
    END DO
         
! World oil price and world carbon price
  CASE (622)
    DO M=1,NMP
      READ(1,*) IDUM,OilPrice(M),CarbonPrice(M)
    END DO
   
! Base-year demand for food in kcal per person per day
  CASE (731)
    DO L=1,NLP
      READ(1,*) IDUM,(kcal1(L,J,2),J=1,10)
    END DO
  
! End-year demand for food in kcal per person per day
  CASE (736)
    DO L=1,NLP
      READ(1,*) IDUM,(kcal1(L,J,9),J=1,10)
    END DO

! Base-year demand - other uses of food products
  CASE (732)
    DO L=1,NLP
      READ(1,*) IDUM,(kcal2(L,J),J=1,10)
    END DO
          
! Net exports of food products in base year (10^10 kcal)
  CASE (733)
    DO L=1,NLP
      READ(1,*) IDUM,(NetExport(L,J),J=1,12)
    END DO
                   
! Calories of crop needed per calorie of processed crop
  CASE (634)
    DO L=1,NLP
      READ(1,*) IDUM,(InputProcCrop(L,J),J=1,4)
    END DO

! Calories of crop needed per calorie of animal product
  CASE (735)
    DO L=1,NLP
      READ(1,*) IDUM,(InputAnimPrds(L,J),J=1,12)
    END DO 
    
! Self-consumption fractions for crops and animal products
  CASE (637)
    DO L=1,NLP
      READ(1,*) IDUM,(SelfCons(L,J),J=1,6)
    END DO

! Yield multipliers for climate impacts
  CASE(623)

    READ(1,*) MSTART,L,(CLIMATE(J,L,MSTART),J=1,NZ(L))
    DO WHILE (L.GT.0)
      DO M=MSTART,NMP
        DO J=1,NZ(L)
          CLIMATE(J,L,M) = CLIMATE(J,L,MSTART)
        END DO
      END DO
      READ(1,*) MSTART,L,(CLIMATE(J,L,MSTART),J=1,NZ(L))
    END DO

! Climate Toggle Variable
  CASE (625)
    READ(1,*) IDUM,CLIMTOGGLE(1),CLIMTOGGLE(2)


! Number of crops, variable by region
  CASE (624)
    DO L=1,NLP
      READ(1,*) IDUM,NC(L)
	  NZ(L) = 5 + NC(L)
    END DO

! Parameters of probability distributions
  CASE (638)
    DO L=1,NLP
      READ(1,*) IDUM,SIG(L),(CORRS(L,J),J=1,3)
    END DO

! CES elasticity
  CASE (639)
    DO L=1,NLP
      READ(1,*) IDUM,agrho(L)
    END DO

! Area Harvested Data
  CASE (701)
	DO L=1,NLP
      READ(1,*) IDUM,(AreaHarv(L,J),J=1,6)
    END DO

! Crop Gross Prod Data
  CASE (702)
	DO L=1,NLP
      READ(1,*) IDUM,(CropSupply(L,J),J=1,4)
    END DO

! Land Area Data
  CASE (703)
	DO L=1,NLP
      READ(1,*) IDUM,(LAND(L,J),J=1,6)
    END DO

! Feed Output Data- used in calibration
  CASE (704)
	DO L=1,NLP
      READ(1,*) IDUM,(FeedOutCalc(L,J),J=1,3)
    END DO
      
! Animal Feed Data
  CASE (705)
	DO L=1,NLP
      READ(1,*) IDUM,(AnimFeed(L,J),J=1,5)
    END DO

! Forest Products Input Data
  CASE (706)
	DO L=1,NLP
      READ(1,*) IDUM,(Forestprds(L,J),J=1,3)
    END DO

! Animal Feed Requirements
  CASE (707)
	DO L=1,NLP
      READ(1,*) IDUM,(FeedReq(L,J),J=1,3)
    END DO

! 1990 Base Year CH4 and N2O Emissions
  CASE (301)
	DO L=1,NLP
      READ(1,*) IDUM,(BaseCH4N2O(L,J),J=1,6)
    END DO

! 1990 percentage of Rice in food grains
  CASE (302)
	DO L=1,NLP
      READ(1,*) IDUM, PCTRICE(L)
    END DO

! All other input blocks are ignored                                
  CASE DEFAULT

  END SELECT

  CALL Ag2NextTable(TABLEMRK,INDIC)
END DO

CLOSE(1)

! Interpolate consumption, in kcal per person per day, between
! base-year and end-year values for the expanded kcal

DO L=1,NLP
  DO J=1,10
    kcal1(L,J,3) = (6.0d0*kcal1(L,J,2) + 1.0d0*kcal1(L,J,9))/7.0d0
    kcal1(L,J,4) = (5.0d0*kcal1(L,J,2) + 2.0d0*kcal1(L,J,9))/7.0d0    
    kcal1(L,J,5) = (4.0d0*kcal1(L,J,2) + 3.0d0*kcal1(L,J,9))/7.0d0
    kcal1(L,J,6) = (3.0d0*kcal1(L,J,2) + 4.0d0*kcal1(L,J,9))/7.0d0
    kcal1(L,J,7) = (2.0d0*kcal1(L,J,2) + 5.0d0*kcal1(L,J,9))/7.0d0
    kcal1(L,J,8) = (1.0d0*kcal1(L,J,2) + 6.0d0*kcal1(L,J,9))/7.0d0
  END DO
END DO

END SUBROUTINE Ag2GetData

! =====================================================================

SUBROUTINE Ag2GDPpercap

USE Ag2Global8

IMPLICIT NONE


! Bring in the read-in population data into popu and calculate the gdp per capita
! This was previously done in AgLUmain

popu(:,1:9) = ZLM(:,2:10)
gdpcap(:,1:9) = 1.0d9 * gdp(:,1:9) / popu(:,1:9)

END SUBROUTINE Ag2GDPpercap


! =====================================================================

SUBROUTINE Ag2GetGDP

! This subroutine is no longer used because GDP is being read in directly
! Some lines are commented out because PROLM & GNPBL are commented out in Global8

! Procedure to calculate per-capita GDP from population and 
! labor productivity.

USE Ag2Global8

IMPLICIT NONE

INTEGER J

REAL(8) gdpindex(NLP,9)

gdpindex = 1.0
!gdp(:,1) = GNPBL / 1.0d6

DO J=2,NMP

!  gdpindex(:,J) = gdpindex(:,J-1) * (1.0 + PROLM(:,J))**15.0 *  &
!                  (ZLM(:,J) / ZLM(:,J-1))

!  gdp(:,J) = GNPBL * gdpindex(:,J) / 1.0d6

END DO

END SUBROUTINE Ag2GetGDP

! =====================================================================

SUBROUTINE Ag2Extend

! Procedure to extend population and gdp data three more time
! steps. This is needed to calculate demand for future
! forest products.

USE Ag2Global8

IMPLICIT NONE

INTEGER J 


DO J=NMP+1,NMP+3

  gdp(:,J) = gdp(:,J-1)
  gdpcap(:,J) = gdpcap(:,J-1)
  popu(:,J) = popu(:,J-1)

END DO

END SUBROUTINE Ag2Extend

! =====================================================================

SUBROUTINE Ag2NextTable(TABLEMRK,INDIC)

IMPLICIT NONE

CHARACTER*11 TABLEIN,TABLEMRK
INTEGER(4) ICNT,INDIC

READ(1,3000,END=10) TABLEIN

ICNT=0
DO WHILE (TABLEIN .NE. TABLEMRK .AND. ICNT .LT. 1000)
  ICNT = ICNT+1
  READ(1,3000,END=10) TABLEIN
END DO
INDIC=1
RETURN

! On end of file
10 CONTINUE
INDIC=0
RETURN

3000 FORMAT(A11)

END SUBROUTINE Ag2NextTable
