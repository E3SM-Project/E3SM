MODULE Ag2Global8

INTEGER,PARAMETER :: NLP = 14, NMP = 9, NUMAX = 100, NZMAX = 12, MAXFILES = 20, NLPMAX = 18, NNLPMax = 19
INTEGER, PARAMETER :: NINP = 12 ! Maximum number of inputs or market goods 
REAL(8),PARAMETER :: InterestR = 0.0200_8, STEP = 15.0_8
REAL(8),PARAMETER :: GJperGcal = 4.1868_8 ! Biomass energy to crop output(GJ per Gcal)
real(8), parameter :: CVRT90 = 2.212 ! GDP deflator from 1975$ to 1990$  
INTEGER IDUM, ICASE,NU,NZ(0:NLP),NC(NLP),CLIMTOGGLE(2)

REAL(8)                     &
YLM(NLPMax,0:NMP),           &
ZLM(NLP,NMP+1),             & ! Population (thousands of persons)
!GNPBL(NLP),                 & ! Base year GDP
!PROLM(NLP,NMP),             & ! Labor productivity growth rates

RJ0(NZMAX,0:NLP),           & ! Mode of yield distribution
KJIMP(NZMAX,0:NLP,0:NMP+3), & ! Annual rate of technical improvement
CALP(NLP,10),               & ! Calibration prices
GJ(NZMAX,0:NLP,0:NMP),      & ! Variable-cost parameter
FEEDOUT(NLP,NMP),           & ! Feed-output ratio
PASTOUT(NLP,NMP),           & ! Pasture-output ratio
DemCoef(NLP,7),             & ! Demand function coefficients
IE(NLP,10),                  & ! Income elasticity of demand. Increased to include proccrops. sjs - 11/01
PE(NLP,10),                  & ! Price elasticity of demand. Increased to include proccrops. sjs - 11/01
BioPrice(NLP,NMP),          & ! Exogenous price of biomass
BioTranCost(NLP,0:NMP),     & ! Cost per GJ to transform biomass to a liquid fuel
Forest(3,NLP,NMP+1),        & ! Forest land in place in period 2
CDensity(NLP,10),            & ! Land carbon density parameters
HistLand(NLP,10),			& ! Historical Land Use for 1960 and 1975
!DemCrop(NLP,4),             & ! Demand parameters for crops
!DemAnim(NLP,4),             & ! Demand parameters for animal products
!SupAnim1990(NLP),           & ! Base-year supply of animal products
!DemAnim1990(NLP),           & ! Base-year demand for animal products
TreeYield(NLP,NMP+3),       & ! Tree yield
OilPrice(NMP),              & ! World oil price
CarbonPrice(NMP),           & ! World carbon price

KJ(NZMAX,NLP,NMP+3),        & ! Technical change coefficients
GDP(NLP,NMP+3),             & ! Gross domestic product
GDPCAP(NLP,NMP+3),          & ! GDP per capita
popu(NLP,NMP+3),            & ! Population

GlobDem(NUMAX),             & !
GlobSup(NUMAX),             &
qdem(8,NLP,NMP),            &
qsup(8,NLP,NMP),            &
price(10,NLP),              &
WPrice(NUMAX),              &
ForestN(NLP,NMP),           & ! New forest land
SaveLand(6,NLP,NMP),        &
SoilEmiss(5,NLP,NMP),		& ! Soil Emissions to be shared out over time periods
CarbEmiss(NLP,NMP),         & ! Carbon emissions from change in land use
SaveProf(NLP,NMP),          &
sprice(10,NLP,NMP),         & ! Store prices for output file
SaveYield(12,NLP,NMP),      & ! Store yields for output file
SaveLandShare(12,NLP,NMP),  & ! Store logit shares
SaveCropLand(4,NLP,NMP),    & ! Store crop land for output.  rs -- 11/01

alpha(NLP),                 & ! Cobb-Douglas alpha
alphaCES(NLP),              & ! CES alpha
beta(NLP),                  & ! Cobb-Douglas beta
betaCES(NLP),               & ! CES beta
constant(NLP),              & ! CES constant
agrho(NLP),					& ! rho value (determined by elasticity)

kcal1(NLP, 10, NMP),        & ! food demand (kcal per person per day)
kcal2(NLP,10),              & ! other uses for food
NetExport(NLP,12),          & ! Net Exports (10^10 kcal)
InputProcCrop(NLP,4),       & ! Inverse efficiencies
InputAnimPrds(NLP,12),		& ! Inverse efficiencies broken up by anim prd
SelfCons(NLP,6),            & ! Self-consumption fractions
SIG(NLP),                   & ! sigma of probability distribution
CORRS(NLP,3),               & ! correlation coefficients
SaveDiet(11,NLP,NMP),       & ! Final consumption in kcal/day/person
CLIMATE(NZMAX,0:NLP,0:NMP), & ! Yield scale factors for climate impacts

AreaHarv(NLP,6),			& ! Crop Area Harvested Data
CropSupply(NLP,4),          & ! Crop Gross Production
LAND(NLP,6),				& ! Land Area Data
FeedOutCalc(NLP,3),			& ! Feed Output calculation data
AnimFeed(NLP,5),			& ! Animal Feed Data
Forestprds(NLP,3),			& ! Forest prds 1990 Supply and Demand
FeedReq(NLP,3),				& ! Animal Feed Requirements
Bioyield(NLP),				& ! Biomass avg observed yields
GJperTON(NLP),				& ! Gigajoules per Ton conversion factor
bioyieldcap(NLP),			& ! Upper limit of biomass yields

AG2CH4(4,NLP,NMP),			& ! Methane Emissions
AG2N2O(4,NLP,NMP),			& ! Nitrous oxide Emissions
BaseCH4N2O(NLP,6),			& ! 1990 CH4 & N20 Emissions
totanimprd(NLP,NMP),		& ! Total animal production
totcropprd(NLP,NMP),		& ! Total crop production
PCTRICE(NLP),				& ! 1990 % Rice in food grains
Porkprod(NLP,NMP),			& ! Total Pork production
Poultryprod(NLP,NMP),		& ! Total Poultry production
CCURVRDX(4,10,NLP,NMP),		& ! Cost curve for mitigation
AGCH4ACT(3,NNLPMax), AGN2OACT(3,NNLPMax), &              
							  ! Arrays for storing output
AgZ(14,14,0:NLP,0:NMP),     & ! Intermediate flows (physical units) for US
AgFD(14,2,0:NLP,0:NMP),     & ! Final demand (physical units) for US
AgP(14,0:NLP,0:NMP),        & ! Production by sector (physical units)
AgPrice(14,0:NLP,0:NMP)       ! Production by sector (at base-year prices)

integer :: JWood, JFWood, JBeef, JFoodGr, JCoarseGr, JOilCrops, JMisccrops, JBio, JPast ! Indices for Ag

! More indices
integer:: INOIL, INGAS, INCOAL, INBMASS, INCARB, INFOREST, INFFOREST, INFOODGR, INCOARSEGR, INOILCROPS, INMISCCROPS, INPAST
END MODULE Ag2Global8