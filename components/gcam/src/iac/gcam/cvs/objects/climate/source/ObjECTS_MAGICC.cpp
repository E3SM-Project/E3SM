#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <time.h>
#include <stdlib.h>


#include "climate/include/ObjECTS_MAGICC.h"

#define DEBUG_IO false

using namespace std;


// The climat() function is up here so as to encapsulate all these stinking variables;
// we're not going to allow any globals in the C++ code
void CLIMAT()
{
    ofstream outfile8; // need to do this here; see line F395 and F658
    openfile_write( &outfile8, "./mag_c.csv", DEBUG_IO );
    
    //F   1 ! MAGTAR.FOR
    //F   2 !
    //F   3 !   THE OUTPUTS PASSED TO THE Tcl CODE ARE:
    //F   4 !    (1) THE FIVE DISPLAY FILES, *.dis, FOR GLOBAL-MEAN TEMP AND MSL,
    //F   5 !        AND EMISSIONS, CONCENTRATIONS, AND RADIATIVE FORCING; AND
    //F   6 !    (2) THE FOUR SCENGEN DRIVER FILES lo/mid/hi/usrdrive.out.
    //F   7 !
    //F   8 ! Revision history:
    //F   9 ! 082009  *  Added fix to MAGOUT file that accounts for USER choice for output
    //F  10 !            reference years other than 1990
    //F  11 ! 062309  *  Added pre-industrial emissions point read-in for to BC-OC forcing. 
    //F  12 ! 022409  *  Added read-in and calc of BC+OC forcing directly from emissions
    //F  13 !            Enable by setting MAGICC input parameter IFOC = 3
    //F  14 ! 013009  *  Updated output arrays for ObjECTS MiniCAM
    //F  15 ! 091708  *  Changed file paths for use as stand-alone model in own directory sjs
    //F  16 
    //F  17 ! 012605  *  Added output of total Kyoto forcing to MiniCAM sjs
    //F  18 ! 030526  *  Revised to work with MiniCAM 5/2003 mrj
    //F  19 ! changes based on those by sjs and hmp in prior MAGICC version
    //F  20 
    //F  21 ! 030518  * sjs added extra halocarbon outputs to mag.csv
    //F  22 ! 030321  * sjs changed mag.out to mag.csv output
    //F  23 ! 03032?  * sjs added SAVE to DELTAQ to correct error
    //F  24 
    //F  25 ! 080625  * FOR NOUT=1, MSL TOTAL LESS TAR "OTHER" ADDED TO OUTPUT.
    //F  26 ! 080619  * ADDED OUTPUT FILE FOR CCSM INPUT, OUTPUT ONLY IF ICCSM=1
    //F  27 ! 080619  * FOR FORCINGS FROM 1765, CFC12 AND EFFECTIVE CFC11 CONCS
    //F  28 !            ADDED AS OUTPUT. THIS REQUIRED ADDING CFC12 CONCS AS AN
    //F  29 !            INPUT IN QHALOS.IN
    //F  30 ! 080617  * MAG.OUT OUTPUT REDUCED EVEN MORE FOR FULL VERSION. THIS
    //F  31 !            REDUCES COLLATION TIME MORE. (SEE LINE 1542.) THIS HAS THE
    //F  32 !            DESIRED EFFECT ON RUN  TIME, BUT IT REMOVES CRUCIAL OUTPUT
    //F  33 !            FROM THE REPORTS FILES. HENCE, BACK TO PREVIOUS VERSION.
    //F  34 !           ERROR IDENTIFIED IN ***drive OUTPUTS THAT DEFINE THE WEIGHTS
    //F  35 !            FOR THE AEROSOL PATTERNS IN SCENGEN. BECAUSE OF SMOOTHING,
    //F  36 !            THESE WERE NOT BEING PRODUCED FOR THE LAST THREE YEARS.
    //F  37 !            CORRECTED USING LINEAR EXTRAPOLATION.
    //F  38 !           THIS IS THE VERSION HANDED OVER TO SETH.
    //F  39 ! 080612  * VZERO IN GSIC MODEL CHANGED TO 18, 29, 44 CM. THESE ARE 1.2
    //F  40 !            TIMES THE AR4 NUMBERS (AS IN AR4) TO ACCOUNT FOR GREENLAND
    //F  41 !            AND ANTARCTICA GSICs. NOTE THAT THE CENTRAL NUMBER MUST BE
    //F  42 !            SPECIFIED IN MAGICE.CFG.
    //F  43 !           THE METHOD FOR CALCULATING SLT SEEMS TO BE WRONG. IT ATTEMPTS
    //F  44 !            TO USE QUADRATURE FOR THE ERRORS IN INDIVIDUAL COMPONENTS,
    //F  45 !            BUT DOES NOT DO THIS CONSISTENTLY. SO IT IS BETTER TO SUM
    //F  46 !            THE INDIVIDUAL COMPONENT VALUES. THIS GIVES EXTREMES THAT
    //F  47 !            ARE TOO HIGH OR TOO LOW. TO GET THE 5th AND 95th %ILES,
    //F  48 !            A GOOD ESTIMATE IS TO HALVE THE DEPARTURES OF THESE EXTREME
    //F  49 !            FROM THE BEST ESTIMATE.
    //F  50 ! 080611  * FOR FULL VERSION (ISCENGEN=1) REDUCE MAG.OUT OUTPUT TO GIVE
    //F  51 !            RESULTS ONLY FOR FULL SO2 EMISSIONS CASE (NESO2=1). THIS
    //F  52 !            IS TO REDUCE COLLATION TIME IN RUNNING MAGICC.
    //F  53 ! 080611  * HANDED OVER TO SETH
    //F  54 ! 080608  * FOR NOUT=1, SUM OF GREENLAND+ANTARCTICA NOW OUTPUTTED.
    //F  55 !           CO2 UNCERTAINTY VALUES FOR DN80S UPDATED TO USE NEW BEST
    //F  56 !            ESTIMATE OF 1.5 (BUT STILL +/- 0.7).
    //F  57 ! 080605  * NUMEROUS CORRECTIONS TO HANDLE QMN CORRECTLY AND KEYED TO 
    //F  58 !            MAGINV.FOR
    //F  59 !           CH4 AND N2O BALANCE METHOD GENERALIZED SO THAT YEAR 2000 DOES 
    //F  60 !            NOT HAVE TO BE D*(2).
    //F  61 ! 080531  * NO3 PLUS MINERAL AEROSOL FORCING (QMN) OUTPUT SEPARATELY
    //F  62 ! 080528  * BEST GUESS SENSITIVITY BACK TO 3.0
    //F  63 !           IDIS NOW BACK TO BEING SET IN MAGRUN.CFG
    //F  64 ! 080527  * NEW GSIC MODEL ADDED FROM MAGMSL.FOR, SING WIGLEY AND RAPER
    //F  65 !            (2005). NEW GSIC PARAMETERS PUT INTO A NEW CFG FILE
    //F  66 !            'MAGICE.CFG'.
    //F  67 ! 080520  * INTERIM VERSION WITH BEST GUESS DT2x = 2.6C TO ACCORD
    //F  68 !            WITH OLD MAG4.1 CFG FILES FOR DEFAULT.
    //F  69 ! 080517  * THREE NEW FORCINGS ADDED: DIRECT NITRATE AEROSOL FORCING
    //F  70 !            QNO3 (-0.1 IN 1990), MINERAL DUST QMIN (-0.1 IN 1990),
    //F  71 !            AND LAND ALBEDO QLAND (-0.2 IN 1990). THESE RAMP UP
    //F  72 !            LINEARLY TO 1990 AND ARE KEPT CONSTANT AFTER THIS.
    //F  73 !           QNO3 AND QMIN ARE SET IN SUBROUTINE SUPHATE AND ADDED TO
    //F  74 !            INDIRECT SO4 AEROSOL FORCING
    //F  75 !           QLAND IS SET IN SUBROUTINE DELTAQ AND PASSED TO MAIN
    //F  76 !           DEFAULT VALUES OF S90DIR, S90IND, S90BIO AND FOC90 RESET
    //F  77 !           BASE TROP OZONE FORCING CORRECTED TO GIVE 0.35 IN 2005 AS
    //F  78 !            IN AR4
    //F  79 !           CLIMATE SENSITIVITY 90% RANGE AND BEST ESTIMATE CHANGED
    //F  80 !            FROM 1.5(2.6)4.5 TO 1.5(3.0)6.0, IN ACCORD WITH AR4.
    //F  81 ! 080517  * LAST 4.1 VERSION BEFORE UPDATING TO VERSION 5.3. SEE THIS
    //F  82 !            VERSION FOR EARLIER HISTORY.
    //F  83 !
    //F  84 !-------------------------------------------------------------------------------
    //F  85 ! Written by Tom Wigley, Sarah Raper & Mike Salmon, Climatic Research Unit, UEA.
    //F  86 !-------------------------------------------------------------------------------
    //F  87 !
    //F  88 !      SUBROUTINE CLIMAT (IWrite, MAGICCCResults,MagEM)	
    //F  89 ! MiniCAM Header with inputs and outputs passed directly
    //F  90 
    //F  91       PROGRAM CLIMAT	
    //F  92 ! Expose subroutine climat to users of this DLL
    //F  93 !
    //F  94 !DEC$ ATTRIBUTES DLLEXPORT::climat
    //F  95       IMPLICIT REAL*4 (a-h,o-z), Integer (I-N)
    //F  96 
    //F  97 !
    //F  98 !   THIS IS THE CLIMATE MODEL MODULE.
    //F  99 !
    //F 100       parameter (iTp=740)
    // This is now defined in header file
    //F 101 !
    //F 102   	  REAL*4 MAGICCCResults(0:30,75) ! sjs pass data to ObjECTS
    float MAGICCCResults[ 22 ][ 75+1 ];
    //F 103 	  REAL*4 getForcing, getSLR
    //F 104   	  
    //F 105       INTEGER IY1(100),OVRWRITE
    int IY1[ 100+1 ], OVRWRITE;
    //F 106 !
    //F 107 ! sjs removed TEMUSER(iTp), QSO2SAVE(0:iTp+1),QDIRSAVE(0:iTp+1) from dimension statement since now in common block
    //F 108       DIMENSION FOS(100),DEF(100),DCH4(100),DN2O(100), &
    float FOS[ 100+1 ], DEF[ 100+1 ], DCH4[ 100+1 ], DN2O[ 100+ 1 ];
    //F 109       DNOX(100),DVOC(100),DCO(100),DSO2(100), &
    float DNOX[ 100+1 ], DVOC[ 100+1 ], DCO[ 100+1 ], DSO2[ 100+1 ]; 
    //F 110       DSO21(100),DSO22(100),DSO23(100),DCF4(100),DC2F6(100), &
    float DSO21[ 100+1 ], DSO22[ 100+1 ], DSO23[ 100+1 ], DCF4[ 100+1 ], DC2F6[ 100+1 ]; 
    //F 111       D125(100),D134A(100),D143A(100),D227(100),D245(100),DSF6(100), &
    float D125[ 100+1 ], D134A[ 100+1 ], D143A[ 100+1 ], D227[ 100+1 ], D245[ 100+1 ], DSF6[ 100+1 ];
    //F 112       DBC(100),DOC(100), &
    float DBC[ 100+1 ], DOC[ 100+1 ];
    //F 113       TEMLO(iTp),TEMMID(iTp),TEMHI(iTp),TEMNOSO2(iTp), &
    float TEMLO[ iTp+1 ], TEMMID[ iTp+1 ], TEMHI[ iTp+1 ], TEMNOSO2[ iTp+1 ];
    //F 114       SLUSER(iTp),SLLO(iTp),SLMID(iTp),SLHI(iTp), &
    float SLUSER[ iTp+1 ], SLLO[ iTp+1 ], SLMID[ iTp+1 ], SLHI[ iTp+1 ];
    //F 115       TALL(4,iTp-225),TGHG(4,iTp-225),TSO21(4,iTp-225), &
    float TALL[ 4+1 ][ iTp-225+1 ], TGHG[ 4+1 ][ iTp-225+1 ], TSO21[ 4+1 ][ iTp-225+1 ];
    //F 116       TSO22(4,iTp-225),TSO23(4,iTp-225),TREF(4),XSO21(4,iTp-225), &
    float TSO22[ 4+1 ][ iTp-225+1 ], TSO23[ 4+1 ][ iTp-225+1 ], TREF[ 4+1 ], XSO21[ 4+1 ][ iTp-225+1 ];
    //F 117       XSO22(4,iTp-225),XSO23(4,iTp-225),XGHG(4,iTp-225), &
    float XSO22[ 4+1 ][ iTp-225+1 ], XSO23[ 4+1 ][ iTp-225+1 ], XGHG[ 4+1 ][ iTp-225+1 ];
    //F 118       SCALER(197:iTp),SCALAR(197:iTp)
    float SCALER[ iTp+1 ], SCALAR[ iTp+1 ];
    //F 119 !
    //F 120       REAL*4 FBC1990, FOC1990, FSO2_dir1990,FSO2_ind1990, aBCUnitForcing, aOCUnitForcing, &
    //F 121              aBCBaseEmissions, aOCBaseEmissions !sjs
    //F 122       COMMON/BCOC/FBC1990, FOC1990, FSO2_dir1990,FSO2_ind1990, aBCUnitForcing, aOCUnitForcing, &
    //F 123              aBCBaseEmissions, aOCBaseEmissions !sjs
    BCOC_block BCOC;
    //F 124       DATA aBCUnitForcing/0.0/,aOCUnitForcing/0.0/
    //F 125 
    //F 126       DIMENSION EESS1(iTp),EESS2(iTp),EESS3(iTp),EESST(iTp),QTROZ(iTp), &
    float QTROZ[ iTp+1 ];
    //F 127       QSTROZ(iTp),C11EFF(iTp)
    float QSTROZ[ iTp+1 ], C11EFF[ iTp+1 ];
    //F 128 !
    //F 129       common /Limits/KEND
    Limits_block Limits;
    //F 130 !
    //F 131       COMMON/OZ/OZ00CH4,OZCH4,OZNOX,OZCO,OZVOC
    OZ_block OZ;
    //F 132 !
    //F 133       COMMON/CLIM/IC,IP,KC,DT,DZ,FK,HM,Q2X,QXX,PI,T,TE,TEND,W0,XK,XKLO, &
    //F 134       XKNS,XLAM,FL(2),FO(2),FLSUM,FOSUM,HEM(2),P(40),TEM(40),TO(2,40), &
    //F 135       AL,BL,CL,DTH,DTZ,DZ1,XLL,WWW,XXX,YYY,RHO,SPECHT,HTCONS,Y(4)
    CLIM_block CLIM;
    //F 136 !
    //F 137 !  NOTE THAT EMISSIONS START WITH J=226, =1990.
    //F 138 !
    //F 139       COMMON/CONCS/CH4(0:iTp),CN2O(0:iTp),ECH4(226:iTp+1), &
    //F 140       EN2O(226:iTp+1),ECO(226:iTp+1),COE(iTp+1),EVOC(226:iTp+1), &
    //F 141       ENOX(226:iTp+1),ESO2(0:iTp+1),ESO2SUM(226:iTp+1), &
    //F 142       ESO21(226:iTp+1),ESO22(226:iTp+1),ESO23(226:iTp+1), &
    //F 143       EBC(226:iTp+1), EOC(226:iTp+1) ! sjs- add BC-OC
    CONCS_block CONCS;
    CONCS.ECH4.init( "CONCs.ECH4", 226, iTp+1 );
    CONCS.EN2O.init( "CONCS.EN2O", 226, iTp+1 );
    CONCS.ECO.init( "CONCS.ECO", 226, iTp+1 );
    CONCS.COE.init( "CONCS.COE", 1, iTp+1 );
    CONCS.EVOC.init( "CONCS.EVOC", 226, iTp+1 );
    CONCS.ENOX.init( "CONCS.ENOX", 226, iTp+1 );
    CONCS.ESO2.init( "CONCS.ESO2", 0, iTp+1 );
    CONCS.ESO2SUM.init( "CONCS.ESO2SUM", 226, iTp+1 );
    CONCS.ESO21.init( "CONCS.ESO21", 226, iTp+1 );
    CONCS.ESO22.init( "CONCS.ESO22", 226, iTp+1 );
    CONCS.ESO23.init( "CONCS.ESO23", 226, iTp+1 );
    CONCS.EBC.init( "CONCS.EBC", 226, iTp+1 );
    CONCS.EOC.init( "CONCS.EOC", 226, iTp+1 );
    //F 144 !
    //F 145       COMMON/NEWCONCS/CF4(iTp),C2F6(iTp),C125(iTp),C134A(iTp), &
    //F 146       C143A(iTp),C227(iTp),C245(iTp),CSF6(iTp), &
    //F 147       ECF4(226:iTp+1),EC2F6(226:iTp+1),E125(226:iTp+1),E134A(226:iTp+1), &
    //F 148       E143A(226:iTp+1),E227(226:iTp+1),E245(226:iTp+1),ESF6(226:iTp+1)
    NEWCONCS_block NEWCONCS;
    NEWCONCS.ECF4.init( "NEWCONCS.ECF4", 226, iTp+1 );
    NEWCONCS.EC2F6.init( "NEWCONCS.EC2F6", 226, iTp+1 );
    NEWCONCS.E125.init( "NEWCONCS.E125", 226, iTp+1 );
    NEWCONCS.E134A.init( "NEWCONCS.E134A", 226, iTp+1 );
    NEWCONCS.E143A.init( "NEWCONCS.E143A", 226, iTp+1 );
    NEWCONCS.E227.init( "NEWCONCS.E227", 226, iTp+1 );
    NEWCONCS.E245.init( "NEWCONCS.E245", 226, iTp+1 );
    NEWCONCS.ESF6.init( "NEWCONCS.ESF6", 226, iTp+1 );
    //F 149 !
    //F 150       COMMON/COBS/COBS(0:236)
    COBS_block COBS;
    //F 151 !
    //F 152       COMMON/CARB/CCO2(4,224:iTp),EDGROSS(4,226:iTp),EF(226:iTp+1), &
    //F 153       REGROW(4,226:iTp),PL(4,226:iTp),HL(4,226:iTp),SOIL(4,226:iTp), &
    //F 154       TTT(226:iTp),ESUM(226:iTp),ETOT(4,226:iTp),EDNET90(4), &
    //F 155       FOC(4,226:iTp),CO2(0:iTp),CO2SAVE(0:iTp)
    CARB_block CARB;
    CARB.CCO2.init( "CARB.CCO2", 1, 4, 224, iTp );
    CARB.EDGROSS.init( "CARB.EDGROSS", 1, 4, 226, iTp );
    CARB.EF.init( "CARB.EF", 226, iTp+1 );
    CARB.REGROW.init( "CARB.REGROW", 1, 4, 226, iTp );
    CARB.PL.init( "CARB.PL", 1, 4, 226, iTp );
    CARB.HL.init( "CARB.HL", 1, 4, 226, iTp );
    CARB.SOIL.init( "CARB.SOIL", 1, 4, 226, iTp );
    CARB.TTT.init( "CARB.TTT", 226, iTp+1 );
    CARB.ESUM.init( "CARB.ESUM", 226, iTp+1 );
    CARB.ETOT.init( "CARB.ETOT", 1, 4, 226, iTp );
    CARB.FOC.init( "CARB.FOC", 1, 4, 226, iTp );
    //F 156 !
    //F 157       COMMON/TANDSL/TEQU(iTp),TGAV(iTp),TNHO(iTp), &
    //F 158       TSHO(iTp),TNHL(iTp),TSHL(iTp),TDEEP(iTp),TNHAV(iTp),TSHAV(iTp), &
    //F 159       TLAND(iTp),TOCEAN(iTp),TOCN(40),TOCNPREV(40), &
    //F 160       SIP,SGP,SAP,SLI(iTp),SLG(iTp),SLA(iTp),EX(0:iTp),SLT(iTp), &
    //F 161       QTOT(0:iTp),QGH(0:iTp),QOZ(0:iTp),QBIO(0:iTp),SLO(iTp), &
    //F 162       QSO2(0:iTp+1),QDIR(0:iTp+1),QLAND(0:iTp),QMN(0:iTp+1)
    TANDSL_block TANDSL;
    //F 163 !
    //F 164       COMMON/CAR/EL1,EL2,EL3,TINV0(5),TINV(4,5),A(3,5),AA(4,5), &
    //F 165       BCO2(4),BTGPP,BTRESP,BTHUM,GAMP,GPP0,RESP0,QA0,U0,C0,B340(4), &
    //F 166       PHI,RG,TAUP,TAUH,TAUS,THP,THS,THH0,THS0,THPL,G1,G2,G3,FACTOR, &
    //F 167       EL21,EL32,XX1,XX2,XX3,XX4,XX5,XX6,DEE1,DEE2,DEE3,DEE4,DEE5,DEE6, &
    //F 168       FL1,FL2,FL3,XL,GAMH,GAMS,QS0,BTSOIL,FERTTYPE,TOTEM,CONVTERP, &
    //F 169       R(4),CPART(4,5),DELMASS(4,226:iTp),ABFRAC(4,226:iTp)
    CAR_block CAR;
    CAR.TINV0.init( "CAR.TINV0", 1, 5 );
    CAR.TINV.init( "CAR.TINV", 1, 4, 1, 5 );
    CAR.A.init( "CAR.A", 1, 3, 1, 5 );
    CAR.AA.init( "CAR.AA", 1, 4, 1, 5 );
    CAR.DELMASS.init( "CAR.DELMASS", 1, 4, 226, iTp );
    CAR.ABFRAC.init( "CAR.ABFRAC", 1, 4, 226, iTp );
    //F 170 !
    //F 171       COMMON /METH1/emeth(226:iTp),imeth,ch4l(225:iTp),ch4b(225:iTp), &
    //F 172       ch4h(225:iTp),ef4(226:iTp),StratH2O,TCH4(iTp),iO3feed, &
    //F 173       ednet(226:iTp+1),DUSER,FUSER,CORRUSER,CORRMHI,CORRMMID,CORRMLO
    METH1_block METH1;
    METH1.emeth.init( "METH1.emeth", 226, iTp );
    METH1.ch4l.init( "METH1.ch4l", 225, iTp );
    METH1.ch4b.init( "METH1.ch4b", 225, iTp );
    METH1.ch4h.init( "METH1.ch4h", 225, iTp );
    METH1.ef4.init( "METH1.ef4", 226, iTp );
    METH1.ednet.init( "METH1.ednet", 226, iTp+1 );
    //F 174 !
    //F 175       COMMON /FORCE/qco2(0:iTp),qm(0:iTp),qn(0:iTp),QCFC(0:iTp), &
    //F 176       QMONT(0:iTp),QOTHER(0:iTp),QSTRATOZ(0:iTp),QCH4O3(0:iTp), &
    //F 177       CFC12(0:iTp), QCH4H2O(0:iTp),QBC(0:iTp),QOC(0:iTp)
    FORCE_block FORCE;
    //F 178 !
    //F 179       COMMON /METH2/LEVCH4,ch4bar90,QQQN2O
    METH2_block METH2;
    //F 180 !
    //F 181       COMMON /METH3/TCH4CON,TAUINIT,SCH4,DELSS,DELTAU, &
    //F 182       ANOX,ACO,AVOC,DELANOX,DELACO,DELAVOC,ICH4FEED
    METH3_block METH3;
    //F 183 !
    //F 184       COMMON /METH4/GAM,TAUOTHER,BBCH4,CM00
    METH4_block METH4;
    //F 185       common /TauNitr/TN2000,BBN2O,SN2O,CN00,NOFFSET
    TauNitr_block TauNitr;
    //F 186       common /Sulph/S90DIR,S90IND,S90BIO,enat,ES1990,ECO90,FOC90,IFOC
    Sulph_block Sulph;
    //F 187       COMMON /CO2READ/ICO2READ,XC(226:iTp),CO2SCALE,qtot86,LEVCO2
    CO2READ_block CO2READ;
    //F 188       COMMON /DSENS/IXLAM,XLAML,XLAMO,ADJUST
    DSENS_block DSENS;
    //F 189       COMMON /VARW/Z(40),W(2),DW(2),TO0(2),TP0(2),WNH(iTp),WSH(iTp), &
    //F 190       TW0NH,TW0SH,IVARW,KEYDW
    VARW_block VARW;
    //F 191       COMMON /QSPLIT/QNHO,QNHL,QSHO,QSHL,QGLOBE(0:iTp), &
    //F 192       QQNHO(0:iTp),QQNHL(0:iTp),QQSHO(0:iTp),QQSHL(0:iTp), &
    //F 193       QQQNHO(0:iTp),QQQNHL(0:iTp),QQQSHO(0:iTp),QQQSHL(0:iTp), &
    //F 194       EHistBC(iTp),EHistOC(iTp) ! Vars to store read-in BC-OC history.
    QSPLIT_block QSPLIT;
    //F 195 !
    //F 196       COMMON /ICE/T1990,G1990,SEN,SENG,SENA,ERRG,ERRA, &
    //F 197       DMG,DMA,SENI,SENP,SENS,DSENI,DSENP,DSENS,ICE,MODEL, &
    //F 198       NEWGSIC,IXG,VZERO,XG
    ICE_block ICE;
    //F 199 !
    //F 200       COMMON /AREAS/FNO,FNL,FSO,FSL
    AREAS_block AREAS;
    //F 201 !
    //F 202       COMMON /QADD/IQREAD,OrgIQREAD,JQFIRST,JQLAST,QEX(0:iTp),QEXNH(0:iTp), &
    //F 203       QEXSH(0:iTp),QEXNHO(0:iTp),QEXNHL(0:iTp),QEXSHO(0:iTp), &
    //F 204       QEXSHL(0:iTp),IOLDTZ
    QADD_block QADD;
    //F 205 !
    //F 206       COMMON /NSIM/NSIM,NCLIM,ISCENGEN,TEMEXP(2,40),IWNHOFF,IWSHOFF, &
    //F 207       WTHRESH
    NSIM_block NSIM;
    //F 208 !
    //F 209       COMMON /JSTART/JSTART,FOSSHIST(0:236),QKYMAG(0:iTp),IGHG, &
    //F 210       QCH4OZ,QFOC(0:iTp),ICO2CORR,TROZSENS
    JSTART_block JSTART;
    //F 211 !
    //F 212       COMMON /CORREN/CORREN1,CORREN2,CORREN3,CORREN4,CORREN
    CORREN_block CORREN;
    //F 213 !
    //F 214 ! sjs -- add storage for halocarbon variables
    //F 215       COMMON /HALOF/QCF4_ar(0:iTp),QC2F6_ar(0:iTp),qSF6_ar(0:iTp), &
    //F 216        Q125_ar(0:iTp),Q134A_ar(0:iTp), &
    //F 217        Q143A_ar(0:iTp),Q227_ar(0:iTp),Q245_ar(0:iTp)
    HALOF_block HALOF;
    //F 218 
    //F 219 !sjs -- parameters that can be modified directly from ObjECTS
    //F 220 ! Note need to add the appropriate model variables to a common block if they are not in one
    //F 221 ! already so that they can be set via subroutine
    //F 222 	  REAL*4 aNewClimSens, aNewBTsoil, aNewBTGPP,aNewBTHumus,aNewDUSER,aNewFUSER, aNewSO2dir1990, aNewSO2ind1990
    //F 223       COMMON/NEWPARAMS/aNewClimSens, aNewBTsoil, DT2XUSER,aNewBTGPP,aNewBTHumus,aNewDUSER,aNewFUSER, &
    //F 224       					aNewSO2dir1990, aNewSO2ind1990
    NEWPARAMS_block NEWPARAMS;
    //F 225       DATA aNewClimSens/-1.0/,aNewBTsoil/-1.0/,aNewBTHumus/-1.0/
    //F 226       DATA aNewDUSER/-1.0/,aNewFUSER/-1.0/,aNewBTGPP/-1.0/
    //F 227 
    //F 228 !Store temperature values
    //F 229       COMMON/STOREDVALS/ TEMUSER(iTp),QSO2SAVE(0:iTp+1),QDIRSAVE(0:iTp+1), &
    //F 230       KYRREF
    STOREDVALS_block STOREDVALS;
    //F 231 
    
    int NCOLS, IQFIRST, IQLAST;
    float EQUIVCO2, TORREF;

    // The MAGICC routines need access to these data structures, so set some global references
    setLocals( &CARB, &TANDSL, &CONCS, &NEWCONCS, 
                &STOREDVALS, &NEWPARAMS, &BCOC, 
                &METH1, &CAR, &FORCE, &JSTART,
                &QADD, &HALOF );
    
    // Code to mimic the climat() data statements (lines F3038-F3055)
    //F3038       DATA FL(1)/0.420/,FL(2)/0.210
    CLIM.FL[ 1 ] = 0.420; CLIM.FL[ 2 ] = 0.210;
    //F3040       DATA RHO/1.026/,SPECHT/0.9333/,HTCONS/4.1856/
    CLIM.RHO = 1.026; CLIM.SPECHT = 0.9333; CLIM.HTCONS = 4.1856;
    //F3045       DATA EL1/141./,EL2/565./,EL3/2500./
    CAR.EL1 = 141.0; CAR.EL2 = 565.0; CAR.EL3 = 2500.0;
    //F3047       DATA DEE1/0.25/,DEE2/0.5/,DEE3/1.0/,DEE4/2.0/ &
    //F3048       ,DEE5/4.0/,DEE6/8.0/
    CAR.DEE1 = 0.25; CAR.DEE2 = 0.5; CAR.DEE3 = 1.0; CAR.DEE4 = 2.0; 
    CAR.DEE5 = 4.0; CAR.DEE6 = 8.0;
    //F3052       DATA (TINV0(J),J=1,5)/0.0,0.0030303,0.0125,0.05,0.625/
    CAR.TINV0.setval( 0.0, 1 );
    CAR.TINV0.setval( 0.0030303, 2 );
    CAR.TINV0.setval( 0.0125, 3 );
    CAR.TINV0.setval( 0.05, 4 );
    CAR.TINV0.setval( 0.625, 5 );
  /*  CAR.TINV0[ 1 ] = 0.0; CAR.TINV0[ 2 ] = 0.0030303; CAR.TINV0[ 3 ] = 0.0125;
    CAR.TINV0[ 4 ] = 0.05; CAR.TINV0[ 5 ] = 0.625; */
    //F3053       DATA (A(1,J),J=1,5)/0.131,0.216,0.261,0.294,0.098/
    CAR.A.setval( 0.131, 1, 1 );
    CAR.A.setval( 0.216, 1, 2 );
    CAR.A.setval( 0.261, 1, 3 );
    CAR.A.setval( 0.294, 1, 4 );
    CAR.A.setval( 0.098, 1, 5 );
/*    CAR.A[ 1 ][ 1 ] = 0.131; CAR.A[ 1 ][ 2 ] = 0.216; CAR.A[ 1 ][ 3 ] = 0.261;
    CAR.A[ 1 ][ 4 ] = 0.294; CAR.A[ 1 ][ 5 ] = 0.098; */
    //F3054       DATA (A(2,J),J=1,5)/0.142,0.230,0.335,0.198,0.095/
    CAR.A.setval( 0.142, 2, 1 );
    CAR.A.setval( 0.230, 2, 2 );
    CAR.A.setval( 0.335, 2, 3 );
    CAR.A.setval( 0.198, 2, 4 );
    CAR.A.setval( 0.095, 2, 5 );
/*    CAR.A[ 2 ][ 1 ] = 0.142; CAR.A[ 2 ][ 2 ] = 0.230; CAR.A[ 2 ][ 3 ] = 0.335;
    CAR.A[ 2 ][ 4 ] = 0.198; CAR.A[ 2 ][ 5 ] = 0.095; */
    //F3055       DATA (A(3,J),J=1,5)/0.166,0.363,0.304,0.088,0.079/
    CAR.A.setval( 0.166, 3, 1 );
    CAR.A.setval( 0.363, 3, 2 );
    CAR.A.setval( 0.304, 3, 3 );
    CAR.A.setval( 0.088, 3, 4 );
    CAR.A.setval( 0.079, 3, 5 );
/*    CAR.A[ 3 ][ 1 ] = 0.166; CAR.A[ 3 ][ 2 ] = 0.363; CAR.A[ 3 ][ 3 ] = 0.304;
    CAR.A[ 3 ][ 4 ] = 0.088; CAR.A[ 3 ][ 5 ] = 0.079; */
    
    
    //F 232 
    //F 233 !  ********************************************************************
    //F 234 !
    //F 235       CHARACTER*4 LEVEL(4)
    string LEVEL[ 4+1 ];
    //F 236 !
    //F 237       CHARACTER*7 MODNAME(0:10)
    string MODNAME[ 10+1 ];
    //F 238 !
    //F 239       character*20 mnem
    string mnem;
    //F 240       character*3  month(12)
    string month[ 12+1 ];
    //F 241 !
    //F 242       DATA (LEVEL(I),I=1,4) /' LOW',' MID','HIGH','USER'/
    LEVEL[ 1 ] = " LOW"; LEVEL[ 2 ] = " MID"; LEVEL[ 3 ] = "HIGH"; LEVEL[ 4 ] = "USER";
    //F 243 !
    //F 244       DATA (MODNAME(I),I=0,10) /' MAGICC','   GFDL','  CSIRO', &
    //F 245       ' HadCM3',' HadCM2',' ECHAM4','    CSM','    PCM', &
    //F 246       'OAGCM 8','OAGCM 9','OAGCM 0'/
    MODNAME[ 0 ] = " MAGICC"; MODNAME[ 1 ] = "   GFDL"; MODNAME[ 2 ] = "  CSIRO";
    MODNAME[ 3 ] = " HadCM3"; MODNAME[ 4 ] = " HadCM2"; MODNAME[ 5 ] = " ECHAM4";
    MODNAME[ 6 ] = "    CSM"; MODNAME[ 7 ] = "    PCM"; MODNAME[ 8 ] = "OAGCM 8";
    MODNAME[ 9 ] = "OAGCM 9";MODNAME[ 10 ] = "OAGCM 0"; 
    //F 247 !
    //F 248       data (month(i),i=1,12) /'Jan','Feb','Mar','Apr','May','Jun', &
    //F 249                               'Jul','Aug','Sep','Oct','Nov','Dec'/
    month[ 1 ] = "Jan"; month[ 2 ] = "Feb"; month[ 3 ] = "Mar"; month[ 4 ] = "Apr";
    month[ 5 ] = "May"; month[ 6 ] = "Jun"; month[ 7 ] = "Jul"; month[ 8 ] = "Aug";
    month[ 9 ] = "Sep"; month[ 10 ] = "Oct"; month[ 11 ] = "Nov"; month[ 12 ] = "Dec";
    //F 250 !
    //F 251       JSTART=236
    JSTART.JSTART = 236;
    //F 252 !
    //F 253 !  READ CO2 CONCENTRATION AND FOSSIL EMISSIONS HISTORIES
    //F 254 !
    //F 255       lun = 42   ! spare logical unit no.
    //F 256       open(unit=lun,file='./magicc_files/CO2HIST.IN',status='OLD')
    ifstream infile;
    openfile_read( &infile, "./co2hist_c.in", DEBUG_IO );
    //F 257       DO ICO2=0,JSTART
    for( int ICO2=0; ICO2<=JSTART.JSTART; ICO2++ ) {
        //F 258       READ(LUN,4445)IIII,COBS(ICO2),FOSSHIST(ICO2)
        double yr; infile >> yr >> COBS.COBS[ ICO2 ] >> JSTART.FOSSHIST[ ICO2 ];
        //F 259       END DO
    }
    //F 260       CLOSE(lun)
    infile.close();
    //F 261 !
    //F 262 !  READ PARAMETERS FROM MAGUSER.CFG.
    //F 263 !
    //F 264       lun = 42   ! spare logical unit no.
    //F 265       open(unit=lun,file='./magicc_files/MAGUSER.CFG',status='OLD')
    openfile_read( &infile, "./maguser_c.cfg", DEBUG_IO );
    //F 266 !
    //F 267         READ(LUN,4240) LEVCO2
    CO2READ.LEVCO2 = read_and_discard( &infile, false );
    //F 268         READ(LUN,4240) IFEED
    const int IFEED = read_and_discard( &infile, false );
    //F 269         READ(LUN,4240) LEVSO4
    int LEVSO4 = read_and_discard( &infile, false );
    //F 270         READ(LUN,4241) DT2XUSER
    NEWPARAMS.DT2XUSER = read_and_discard( &infile, false );
    //F 271         READ(LUN,4241) YK
    float YK = read_and_discard( &infile, false );
    //F 272         READ(LUN,4240) IVARW
    VARW.IVARW = read_and_discard( &infile, false );
    //F 273         READ(LUN,4240) MODEL
    ICE.MODEL = read_and_discard( &infile, false );
    //F 274         READ(LUN,4240) KYRREF
    const int KYRREF = read_and_discard( &infile, false );
    //F 275         READ(LUN,4240) IY0
    const int IY0 = read_and_discard( &infile, false );
    //F 276         READ(LUN,4240) LASTYEAR
    int LASTYEAR = read_and_discard( &infile, false );
    //F 277         READ(LUN,4240) ITEMPRT
    const int ITEMPRT = read_and_discard( &infile, false );
    //F 278         READ(LUN,4240) ICE
    ICE.ICE = read_and_discard( &infile, false );
    //F 279         NONOFF=0
    const int NONOFF = 0;
    //F 280 !
    //F 281       close(lun)
    infile.close();
    //F 282 !
    //F 283       LASTMAX=1764+iTp
    const int LASTMAX = 1764 + iTp;
    //F 284       IF(LASTYEAR.GT.LASTMAX)LASTYEAR=LASTMAX
    if( LASTYEAR > LASTMAX ) LASTYEAR = LASTMAX;
    //F 285 !
    //F 286 !  ********************************************************************
    //F 287 !
    //F 288 !  READ PARAMETERS FROM MAGICE.CFG.
    //F 289 !
    //F 290       lun = 42   ! spare logical unit no.
    //F 291       open(unit=lun,file='./magicc_files/MAGICE.CFG',status='OLD')
    openfile_read( &infile, "./magice_c.cfg", DEBUG_IO );
    //F 292 !
    //F 293         READ(LUN,4240) NEWGSIC  ! SET = 1 TO USE NEW ALGORITHM
    ICE.NEWGSIC = read_and_discard( &infile, false );
    //F 294         READ(LUN,4241) VZERO
    ICE.VZERO = read_and_discard( &infile, false );
    //F 295         READ(LUN,4241) XG
    ICE.XG = read_and_discard( &infile, false );
    //F 296         READ(LUN,4240) IXG
    ICE.IXG = read_and_discard( &infile, false );
    //F 297         READ(LUN,4241) ASEN
    const float ASEN = read_and_discard( &infile, false );
    //F 298 !
    //F 299       CLOSE(lun)
    infile.close();
    //F 300 !
    //F 301 !  ********************************************************************
    //F 302 !
    //F 303 !  READ PARAMETERS FROM MAGGAS.CFG. ORDER OF PARAMETERS CHANGED
    //F 304 !   AND OTHER ITEMS ADDED IN SEPT 2000.
    //F 305 !
    //F 306       lun = 42   ! spare logical unit no.
    //F 307       open(unit=lun,file='./magicc_files/MAGGAS.CFG',status='OLD')
    openfile_read( &infile, "./maggas_c.cfg", DEBUG_IO );
    //F 308 !
    //F 309         READ(LUN,4240) OVRWRITE
    OVRWRITE = read_and_discard( &infile, false );
    //F 310         READ(LUN,4241) DUSER
    METH1.DUSER = read_and_discard( &infile, false );
    //F 311         READ(LUN,4241) FUSER
    METH1.FUSER = read_and_discard( &infile, false );
    //F 312         READ(LUN,4241) BTGPP
    CAR.BTGPP = read_and_discard( &infile, false );
    //F 313         READ(LUN,4241) BTRESP
    CAR.BTRESP = read_and_discard( &infile, false );
    //F 314         READ(LUN,4241) BTHUM
    CAR.BTHUM = read_and_discard( &infile, false );
    //F 315         READ(LUN,4241) BTSOIL
    CAR.BTSOIL = read_and_discard( &infile, false );
    //F 316         READ(LUN,4240) IMETH
    METH1.IMETH = read_and_discard( &infile, false );
    //F 317         READ(LUN,4240) LEVCH4
    METH2.LEVCH4 = read_and_discard( &infile, false );
    //F 318         READ(LUN,4241) TCH4CON
    METH3.TCH4CON = read_and_discard( &infile, false );
    //F 319         READ(LUN,4241) TAUINIT
    METH3.TAUINIT = read_and_discard( &infile, false );
    //F 320         READ(LUN,4241) DELTAU
    METH3.DELTAU = read_and_discard( &infile, false );
    //F 321         READ(LUN,4241) TAUSOIL
    const float TAUSOIL = read_and_discard( &infile, false );
    //F 322         READ(LUN,4241) TAUSTRAT
    const float TAUSTRAT = read_and_discard( &infile, false );
    //F 323         READ(LUN,4241) CMBAL
    const float CMBAL = read_and_discard( &infile, false );
    //F 324         READ(LUN,4241) DCMBAL
    const float DCMBAL = read_and_discard( &infile, false );
    //F 325         READ(LUN,4241) TN2000
    TauNitr.TN2000 = read_and_discard( &infile, false );
    //F 326         READ(LUN,4241) CNBAL
    const float CNBAL = read_and_discard( &infile, false );
    //F 327         READ(LUN,4241) DCNBAL
    const float DCNBAL = read_and_discard( &infile, false );
    //F 328         READ(LUN,4241) QQQN2O
    METH2.QQQN2O = read_and_discard( &infile, false );
    //F 329         READ(LUN,4241) S90Duser
    float S90Duser = read_and_discard( &infile, false );
    //F 330         READ(LUN,4241) S90Iuser
    float S90Iuser = read_and_discard( &infile, false );
    //F 331         READ(LUN,4241) S90Buser
    float S90Buser = read_and_discard( &infile, false );
    //F 332         READ(LUN,4241) FOC90usr
    float FOC90usr = read_and_discard( &infile, false );
    //F 333         READ(LUN,4240) IFOC
    Sulph.IFOC = read_and_discard( &infile, false );
    //F 334         READ(LUN,4241) TROZSENS
    JSTART.TROZSENS = read_and_discard( &infile, false );
    //F 335         READ(LUN,4240) IO3FEED
    METH1.IO3FEED = read_and_discard( &infile, false );
    //F 336         READ(LUN,4240) IHALOETC
    const int IHALOETC = read_and_discard( &infile, false );
    //F 337         READ(LUN,4241) OZ00CH4
    OZ.OZ00CH4 = read_and_discard( &infile, false );
    //F 338         READ(LUN,4241) OZCH4
    OZ.OZCH4 = read_and_discard( &infile, false );
    //F 339         READ(LUN,4241) OZNOX
    OZ.OZNOX = read_and_discard( &infile, false );
    //F 340         READ(LUN,4241) OZCO
    OZ.OZCO = read_and_discard( &infile, false );
    //F 341         READ(LUN,4241) OZVOC
    OZ.OZVOC = read_and_discard( &infile, false );
    //F 342         READ(LUN,4240) ICH4FEED
    //F 343         IOLDOZ=0
    // This appears to be bug in original Fortran code; variables are being read out of order
    // We'll read in correct order
    const int IOLDOZ = read_and_discard( &infile, false ); 
    //F         READ(LUN,4240) ICH4FEED
    METH3.ICH4FEED = read_and_discard( &infile, false );
    //F 344 !
    //F 345       close(lun)
    infile.close();
    //F 346 
    //! Initiailize internal BC-OC vars
    //aBCUnitForcing = 0
    //aOCUnitForcing = 0
    //aNewSO2dir1990 = 0
    //aNewSO2ind1990 = 0
    BCOC.aBCUnitForcing = BCOC.aOCUnitForcing = NEWPARAMS.aNewSO2dir1990 = NEWPARAMS.aNewSO2ind1990 = 0.0;
	BCOC.FSO2_dir1990 = BCOC.FSO2_ind1990 = BCOC.FBC1990 = 0.0;
    
    //F 347 !
    //F 348 !   Call overrite subroutine after each file that may have parameters to overwrite
    //F 349       call overrideParameters( )	! sjs
    overrideParameters( &NEWPARAMS, &CAR, &METH1, &BCOC );
    //F 350 
    //F 351        IF ( FSO2_dir1990 .LT. 0) THEN
    if( BCOC.FSO2_dir1990 < 0 ) {
        //F 352           S90Duser = FSO2_dir1990
        S90Duser = BCOC.FSO2_dir1990;
        //F 353 	      S90Iuser = FSO2_ind1990
        S90Iuser = BCOC.FSO2_ind1990;
        //F 354 	   END IF
    }
    //F 355 
    //F 356        IF ( FBC1990 .NE. 0) THEN
    if( BCOC.FBC1990 != 0 ) {
        //F 357           FOC90usr = 0
        //F 358 	      S90Buser = 0
        FOC90usr = S90Buser = 0;
        //F 359 	   END IF
    }
    //F 360 
    //F 361 !
    //F 362       IF(OVRWRITE.EQ.1)THEN
    if( OVRWRITE == 1 ) {
        //F 363         LEVCO2=4
        //F 364         LEVSO4=4
        CO2READ.LEVCO2 = LEVSO4 = 4;
        //F 365       ENDIF
    }
    //F 366 !
    //F 367 !  NEW METHOD FOR CALCULATING CO2 OUTPUT PDF FROM INVERSE VERSION
    //F 368 !   OF MAGICC. XB=0.074 gives DUSER=1.1, BTGPP=0.015
    //F 369 !
    //F 370       XB=BTSOIL
    //UNUSED const float XB = CAR.BTSOIL;
    //F 371 !      DUSER=1.8*EXP(-6.6551*XB)
    //F 372 !      BTGPP=0.03*EXP(-9.3669*XB)
    //F 373 !
    //F 374       IF(IOLDOZ.EQ.1)THEN
    if( IOLDOZ == 1 ) {
        //F 375          OZ00CH4 = 0.168
        OZ.OZ00CH4 = 0.168;
        //F 376          OZCH4   = 6.2048
        OZ.OZCH4 = 6.2048;
        //F 377          OZNOX   = 0.17
        OZ.OZNOX = 0.17;
        //F 378          OZCO    = 0.0014
        OZ.OZCO = 0.0014;
        //F 379          OZVOC   = 0.0042
        OZ.OZVOC = 0.0042;
        //F 380       ENDIF
    }
    //F 381 !
    //F 382       TAUOTHER=1.0/(1.0/TAUSOIL+1.0/TAUSTRAT)
    METH4.TAUOTHER=1.0/(1.0/TAUSOIL+1.0/TAUSTRAT);
    //F 383 !
    //F 384       IF(IFEED.EQ.0)THEN
    if( IFEED == 0 ) {
        //F 385         BTGPP  = 0.0
        //F 386         BTRESP = 0.0
        //F 387         BTHUM  = 0.0
        //F 388         BTSOIL = 0.0
        CAR.BTGPP = CAR.BTRESP = CAR.BTHUM = CAR.BTSOIL = 0.0;
        //F 389       ENDIF
    }
    //F 390 !
    //F 391 !  TRAP IN CASE LEV* MIS-SPECIFIED OUTSIDE PERMISSIBLE RANGE. IF SO, 
    //F 392 !   RE-SET TO BEST GUESS CASE.
    //F 393 !
    //F 394       IF((LEVCO2.GT.4).OR.(LEVCH4.GT.4).OR.(LEVSO4.GT.4))THEN
    if( ( CO2READ.LEVCO2 > 4 ) || ( METH2.LEVCH4 > 4 ) || (LEVSO4 > 4 ) ) {
        //F 395         WRITE(8,115)
        outfile8 << "LEVCO2, LEVCH4 AND/OR LEVSO4 WRONGLY SET > 4 : RESET AT 2" << endl;
        //F 396       ENDIF
    }
    //F 397 !
    //F 398       IF(LEVCO2.GT.4)LEVCO2=2
    if( CO2READ.LEVCO2 > 4 ) CO2READ.LEVCO2 = 2;
    //F 399       IF(LEVCH4.GT.4)LEVCH4=2
    if( METH2.LEVCH4 > 4 ) METH2.LEVCH4 = 2;
    //F 400       IF(LEVSO4.GT.4)LEVSO4=2
    if( LEVSO4 > 4 ) LEVSO4 = 2;
    //F 401 !
    //F 402 !  READ PARAMETERS FROM MAGMOD.CFG
    //F 403 !
    //F 404       lun = 42   ! spare logical unit no.
    //F 405       open(unit=lun,file='./magicc_files/MAGMOD.CFG',status='OLD')
    openfile_read( &infile, "./magmod_c.cfg", DEBUG_IO );
    //F 406 !
    //F 407         READ(LUN,4241) ADJUST
    DSENS.ADJUST = read_and_discard( &infile, false );
    //F 408         READ(LUN,4241) CO2DELQ
    float CO2DELQ = read_and_discard( &infile, false );
    //F 409         READ(LUN,4241) RLO
    float RLO = read_and_discard( &infile, false );
    //F 410         READ(LUN,4241) HM
    CLIM.HM = read_and_discard( &infile, false );
    //F 411         READ(LUN,4241) W0
    CLIM.W0 = read_and_discard( &infile, false );
    //F 412         READ(LUN,4241) PI
    CLIM.PI = read_and_discard( &infile, false );
    //F 413         READ(LUN,4241) TW0NH
    VARW.TW0NH = read_and_discard( &infile, false );
    //F 414         READ(LUN,4241) TW0SH
    VARW.TW0SH = read_and_discard( &infile, false );
    //F 415         READ(LUN,4241) XKLO
    CLIM.XKLO = read_and_discard( &infile, false );
    //F 416         READ(LUN,4241) XKNS
    CLIM.XKNS = read_and_discard( &infile, false );
    //F 417         READ(LUN,4240) ICO2READ
    CO2READ.ICO2READ = read_and_discard( &infile, false );
    //F 418         READ(LUN,4241) CO2SCALE
    CO2READ.CO2SCALE = read_and_discard( &infile, false );
    //F 419         READ(LUN,4240) IQREAD
    QADD.IQREAD = read_and_discard( &infile, false );
    //F 420         READ(LUN,4241) QOFFSET
    const float QOFFSET = read_and_discard( &infile, false );
    //F 421         READ(LUN,4241) QFACTOR
    const float QFACTOR = read_and_discard( &infile, false );
    //F 422 !
    //F 423       IF(NONOFF.EQ.1)ICO2READ=5
    if( NONOFF == 1 ) CO2READ.ICO2READ = 5;
    //F 424       IF(RLO.GT.2.0)THEN
    if( RLO > 2.0 ) {
        //F 425         RLO=1.05+0.6228*EXP(-0.339*DT2XUSER)
        RLO = 1.05 + 0.6228 * exp( float( -0.339 * NEWPARAMS.DT2XUSER ) );
        //F 426       ENDIF
    }
    //F 427 !
    //F 428       close(lun)
    infile.close();
    //F 429 !
    //F 430 !   Call overrite subroutine after each file that may have parameters to overwrite
    //F 431       call overrideParameters( )	! sjs
    overrideParameters( &NEWPARAMS, &CAR, &METH1, &BCOC );
    //F 432 !
    //F 433       IF(IVARW.EQ.2)THEN
    if( VARW.IVARW == 2 ) {
        //F 434         TW0SH=TW0NH
        VARW.TW0SH = VARW.TW0NH;
        //F 435         TW0=TW0NH
        //UNUSED const float TW0 = VARW.TW0NH;
        //F 436       ENDIF
    }
    //F 437 !
    //F 438 !  'DEFAULT' (= AOGCM MEAN) ICE MELT PARAMETERS FOR MAGICC. NOTE
    //F 439 !    THAT T1990 IS NOT USED FOR THIS CASE.
    //F 440 !
    //F 441       T1990    =  0.607
    ICE.T1990 = 0.607;
    //F 442       G1990    =  2.14     ! CM
    ICE.G1990 = 2.14; // CM
    //F 443       SEN      =  0.0625   ! CM/YR-degC (GSIC)
    ICE.SEN = 0.0625; // CM/YR-degC (GSIC)
    //F 444       SENG     =  0.0110   ! CM/YR-degC (GREENLAND)
    ICE.SENG = 0.0110; // CM/YR-degC (GREENLAND)
    //F 445       SENA     = -0.0341   ! CM/YR-degC (ANTARCTICA)
    ICE.SENA = -0.0341; // CM/YR-degC (ANTARCTICA)
    //F 446       ERRG     =  1.896    ! GREENLAND
    ICE.ERRG = 1.896; // GREENLAND
    //F 447       ERRA     =  1.242    ! GREENLAND
    ICE.ERRA = 1.242; // GREENLAND
    //F 448 !
    //F 449 !  ICE PARAMETERS COMMON TO ALL CASES
    //F 450 !
    //F 451       DMG      =  0.005    ! CM/YR-degC
    ICE.DMG = 0.005; // CM/YR-degC
    //F 452       DMA      =  0.008    ! CM/YR-degC
    ICE.DMA = 0.08; // CM/YR-degC
    //F 453       SENI     =  0.025    ! CM/YR
    ICE.SENI = 0.025; // CM/YR
    //F 454       DSENI    =  0.025    ! CM/YR
    ICE.DSENI = 0.025; // CM/YR
    //F 455       SENP     =  0.01136  ! CM/YR
    ICE.SENP = 0.01136; // CM/YR
    //F 456       DSENP    =  0.01136  ! CM/YR
    ICE.DSENP = 0.01136; // CM/YR
    //F 457       SENS     =  0.0025   ! CM/YR
    ICE.SENS = 0.0025; // CM/YR
    //F 458       DSENS    =  0.0025   ! CM/YR
    ICE.DSENS = 0.0025; // CM/YR
    //F 459 !
    //F 460 !  IF MODEL.NE.0, THEN SELECT MODEL FROM 'LIBRARY' GIVEN BELOW.
    //F 461 !
    //F 462       IF(MODEL.EQ.1)THEN
    if( ICE.MODEL == 1 ) {
        //F 463         CO2DELQ  =  5.352
        //F 464         DT2XUSER =  4.20
        //F 465         TW0NH    =  8.00
        //F 466         YK       =  2.30
        //F 467         RLO      =  1.20
        //F 468         XKLO     =  1.00
        //F 469         T1990    =  0.635
        //F 470         G1990    =  1.5      ! CM
        //F 471         SEN      =  0.0576   ! CM/YR-degC (GSIC)
        //F 472         SENG     =  0.0121   ! CM/YR-degC (GREENLAND)
        //F 473         SENA     = -0.0177   ! CM/YR-degC (ANTARCTICA)
        //F 474         ERRG     =  1.879    ! GREENLAND
        //F 475         ERRA     =  0.799    ! GREENLAND
        CO2DELQ             = 5.352;
        NEWPARAMS.DT2XUSER  = 4.20;
        VARW.TW0NH          = 8.00;
        YK                  = 2.30;
        RLO                 = 1.20;
        CLIM.XKLO           = 1.00;
        ICE.T1990           = 0.635;
        ICE.G1990           = 1.5;      // CM
        ICE.SEN             = 0.0576;   // CM/YR-degC (GSIC)
        ICE.SENG            = 0.0121;   // CM/YR-degC (GREENLAND)
        ICE.SENA            = -0.0177;  // CM/YR-degC (ANTARCTICA)
        ICE.ERRG            = 1.879;    // GREENLAND
        ICE.ERRA            = 0.799;    // ANTARCTICA
        //F 476       ENDIF
    }
    //F 477 !
    //F 478       IF(MODEL.EQ.2)THEN
    if( ICE.MODEL == 2 ) {
        //F 479         CO2DELQ  =  4.977
        //F 480         DT2XUSER =  3.70
        //F 481         TW0NH    =  5.00
        //F 482         YK       =  1.60
        //F 483         RLO      =  1.20
        //F 484         XKLO     =  1.00
        //F 485         T1990    =  0.593
        //F 486         G1990    =  2.2      ! CM
        //F 487         SEN      =  0.0733   ! CM/YR-degC (GSIC)
        //F 488         SENG     =  0.0157   ! CM/YR-degC (GREENLAND)
        //F 489         SENA     = -0.0373   ! CM/YR-degC (ANTARCTICA)
        //F 490         ERRG     =  2.042    ! GREENLAND
        //F 491         ERRA     =  1.120    ! GREENLAND
        CO2DELQ             =  4.977;
        NEWPARAMS.DT2XUSER  =  3.70;
        VARW.TW0NH          =  5.00;
        YK                  =  1.60;
        RLO                 =  1.20;
        CLIM.XKLO           =  1.00;
        ICE.T1990           =  0.593;
        ICE.G1990           =  2.2;          // CM
        ICE.SEN             =  0.0733;       // CM/YR-degC (GSIC)
        ICE.SENG            =  0.0157;       // CM/YR-degC (GREENLAND)
        ICE.SENA            = -0.0373;       // CM/YR-degC (ANTARCTICA)
        ICE.ERRG            =  2.042;        // GREENLAND
        ICE.ERRA            =  1.120;        // ANTARCTICA
        //F 492       ENDIF
    }
    //F 493 !
    //F 494       IF(MODEL.EQ.3)THEN
    if( ICE.MODEL == 3 ) {
        //F 495         CO2DELQ  =  5.396
        //F 496         DT2XUSER =  3.00
        //F 497         TW0NH    = 25.00
        //F 498         YK       =  1.90
        //F 499         RLO      =  1.40
        //F 500         XKLO     =  0.50
        //F 501         T1990    =  0.562
        //F 502         G1990    =  2.1      ! CM
        //F 503         SEN      =  0.0622   ! CM/YR-degC (GSIC)
        //F 504         SENG     =  0.0085   ! CM/YR-degC (GREENLAND)
        //F 505         SENA     = -0.0354   ! CM/YR-degC (ANTARCTICA)
        //F 506         ERRG     =  1.443    ! GREENLAND
        //F 507         ERRA     =  1.288    ! GREENLAND
        CO2DELQ             = 5.396;
        NEWPARAMS.DT2XUSER  = 3.00;
        VARW.TW0NH          = 25.00;
        YK                  = 1.90;
        RLO                 = 1.40;
        CLIM.XKLO           = 0.50;
        ICE.T1990           = 0.562;
        ICE.G1990           = 2.1;      // CM
        ICE.SEN             = 0.0622;   // CM/YR-degC (GSIC)
        ICE.SENG            = 0.0085;   // CM/YR-degC (GREENLAND)
        ICE.SENA            = -0.0354;  // CM/YR-degC (ANTARCTICA)
        ICE.ERRG            = 1.443;    // GREENLAND
        ICE.ERRA            = 1.288;    // ANTARCTICA
        //F 508       ENDIF
    }
    //F 509 !
    //F 510       IF(MODEL.EQ.4)THEN
    if( ICE.MODEL == 4 ) {
        //F 511         CO2DELQ  =  5.006
        //F 512         DT2XUSER =  2.50
        //F 513         TW0NH    = 12.00
        //F 514         YK       =  1.70
        //F 515         RLO      =  1.40
        //F 516         XKLO     =  0.50
        //F 517         T1990    =  0.603
        //F 518         G1990    =  2.7      ! CM
        //F 519         SEN      =  0.0613   ! CM/YR-degC (GSIC)
        //F 520         SENG     =  0.0096   ! CM/YR-degC (GREENLAND)
        //F 521         SENA     = -0.0214   ! CM/YR-degC (ANTARCTICA)
        //F 522         ERRG     =  1.441    ! GREENLAND
        //F 523         ERRA     =  1.239    ! GREENLAND
        CO2DELQ             = 5.006;
        NEWPARAMS.DT2XUSER  = 2.50;
        VARW.TW0NH          = 12.00;
        YK                  = 1.70;
        RLO                 = 1.40;
        CLIM.XKLO           = 0.50;
        ICE.T1990           = 0.603;
        ICE.G1990           = 2.7;      // CM
        ICE.SEN             = 0.0613;   // CM/YR-degC (GSIC)
        ICE.SENG            = 0.0096;   // CM/YR-degC (GREENLAND)
        ICE.SENA            = -0.0214;  // CM/YR-degC (ANTARCTICA)
        ICE.ERRG            = 1.441;    // GREENLAND
        ICE.ERRA            = 1.239;    // ANTARCTICA
        //F 524       ENDIF
    }
    //F 525 !
    //F 526       IF(MODEL.EQ.5)THEN
    if( ICE.MODEL == 5 ) {
        //F 527         CO2DELQ  =  5.482
        //F 528         DT2XUSER =  2.60
        //F 529         TW0NH    = 20.00
        //F 530         YK       =  9.00
        //F 531         RLO      =  1.40
        //F 532         XKLO     =  0.50
        //F 533         T1990    =  0.780
        //F 534         G1990    =  2.7      ! CM
        //F 535         SEN      =  0.0637   ! CM/YR-degC (GSIC)
        //F 536         SENG     =  0.0029   ! CM/YR-degC (GREENLAND)
        //F 537         SENA     = -0.0478   ! CM/YR-degC (ANTARCTICA)
        //F 538         ERRG     =  1.153    ! GREENLAND
        //F 539         ERRA     =  1.484    ! GREENLAND
        CO2DELQ             = 5.482;
        NEWPARAMS.DT2XUSER  = 2.60;
        VARW.TW0NH          = 20.00;
        YK                  = 9.00;
        RLO                 = 1.40;
        CLIM.XKLO           = 0.50;
        ICE.T1990           = 0.780;
        ICE.G1990           = 2.7;      // CM
        ICE.SEN             = 0.0637;   // CM/YR-degC (GSIC)
        ICE.SENG            = 0.0029;   // CM/YR-degC (GREENLAND)
        ICE.SENA            = -0.0478;  // CM/YR-degC (ANTARCTICA)
        ICE.ERRG            = 1.153;    // GREENLAND
        ICE.ERRA            = 1.484;    // ANTARCTICA
        //F 540       ENDIF
    }
    //F 541 !
    //F 542       IF(MODEL.EQ.6)THEN
    if( ICE.MODEL == 6 ) {
        //F 543         CO2DELQ  =  5.194
        //F 544         DT2XUSER =  1.90
        //F 545         TW0NH    = 1000.00
        //F 546         YK       =  2.30
        //F 547         RLO      =  1.40
        //F 548         XKLO     =  0.50
        //F 549         T1990    =  0.567
        //F 550         G1990    =  2.1      ! CM
        //F 551         SEN      =  0.0608   ! CM/YR-degC (GSIC)
        //F 552         SENG     =  0.0146   ! CM/YR-degC (GREENLAND)
        //F 553         SENA     = -0.0305   ! CM/YR-degC (ANTARCTICA)
        //F 554         ERRG     =  3.147    ! GREENLAND
        //F 555         ERRA     =  1.143    ! GREENLAND
        CO2DELQ             = 5.194;
        NEWPARAMS.DT2XUSER  = 1.90;
        VARW.TW0NH          = 1000.00;
        YK                  = 2.30;
        RLO                 = 1.40;
        CLIM.XKLO           = 0.50;
        ICE.T1990           = 0.567;
        ICE.G1990           = 2.1;      // CM
        ICE.SEN             = 0.0608;   // CM/YR-degC (GSIC)
        ICE.SENG            = 0.0146;   // CM/YR-degC (GREENLAND)
        ICE.SENA            = -0.0305;  // CM/YR-degC (ANTARCTICA)
        ICE.ERRG            = 3.147;    // GREENLAND
        ICE.ERRA            = 1.143;    // ANTARCTICA
        //F 556       ENDIF
    }
    //F 557 !
    //F 558       IF(MODEL.EQ.7)THEN
    if( ICE.MODEL == 7 ) {
        //F 559         CO2DELQ  =  5.194
        //F 560         DT2XUSER =  1.70
        //F 561         TW0NH    = 14.00
        //F 562         YK       =  2.30
        //F 563         RLO      =  1.40
        //F 564         XKLO     =  0.50
        //F 565         T1990    =  0.510
        //F 566         G1990    =  1.7      ! CM
        //F 567         SEN      =  0.0587   ! CM/YR-degC (GSIC)
        //F 568         SENG     =  0.0136   ! CM/YR-degC (GREENLAND)
        //F 569         SENA     = -0.0484   ! CM/YR-degC (ANTARCTICA)
        //F 570         ERRG     =  2.165    ! GREENLAND
        //F 571         ERRA     =  1.618    ! GREENLAND
        CO2DELQ             = 5.194;
        NEWPARAMS.DT2XUSER  = 1.70;
        VARW.TW0NH          = 14.00;
        YK                  = 2.30;
        RLO                 = 1.40;
        CLIM.XKLO           = 0.50;
        ICE.T1990           = 0.510;
        ICE.G1990           = 1.7;      // CM
        ICE.SEN             = 0.0587;   // CM/YR-degC (GSIC)
        ICE.SENG            = 0.0136;   // CM/YR-degC (GREENLAND)
        ICE.SENA            = -0.0484;  // CM/YR-degC (ANTARCTICA)
        ICE.ERRG            = 2.165;    // GREENLAND
        ICE.ERRA            = 1.618;    // ANTARCTICA
        //F 572       ENDIF
    }
    //F 573 !
    //F 574 !  OVERWRITE GSIC SENSITIVITY
    //F 575 !
    //F 576       SEN=ASEN*SEN
    ICE.SEN = ASEN * ICE.SEN;
    //F 577 !
    //F 578       IF(MODEL.NE.0)THEN
    if( ICE.MODEL != 0 ) {
        //F 579         TW0SH=TW0NH
        VARW.TW0SH = VARW.TW0NH;
        //F 580         XKNS=XKLO
        CLIM.XKNS = CLIM.XKLO;
        //F 581       ENDIF
    }
    //F 582 !
    //F 583 !  ********************************************************************
    //F 584 !
    //F 585 !  READ PARAMETERS FROM MAGRUN.CFG
    //F 586 !  NOTE (051308): IDIS NOW DEFINED IN MAGUSER.CFG TO BE THE SAME
    //F 587 !   AS ITEMPRT. IDIS BELOW RELABELLED AS JDIS TO AVOID OVERWRITING.
    //F 588 !  052808: IDIS BACK TO BEING SPECIFIED HERE
    //F 589 !
    //F 590       lun = 42   ! spare logical unit no.
    //F 591       open(unit=lun,file='./magicc_files/MAGRUN.CFG',status='OLD')
    openfile_read( &infile, "./magrun_c.cfg", DEBUG_IO );
    //F 592 !
    //F 593         READ(LUN,4240) ISCENGEN
    NSIM.ISCENGEN = read_and_discard( &infile, false );
    //F 594         READ(LUN,4240) IDIS
    const int IDIS = read_and_discard( &infile, false );
    //F 595 !        READ(LUN,4240) JDIS
    //F 596         READ(LUN,4240) KSTART
    const int KSTART = read_and_discard( &infile, false );
    //F 597         READ(LUN,4240) ICO2CORR
    JSTART.ICO2CORR = read_and_discard( &infile, false );
    //F 598         READ(LUN,4240) KEYDW
    VARW.KEYDW = read_and_discard( &infile, false );
    //F 599         READ(LUN,4240) IMAGTAR
    /* //UNUSED const int IMAGTAR = */ read_and_discard( &infile, false );
    //F 600         READ(LUN,4240) IPEAK
    /* //UNUSED const int IPEAK = */ read_and_discard( &infile, false );
    //F 601         READ(LUN,4241) DPEAK
    const float DPEAK = read_and_discard( &infile, false );
    //F 602         READ(LUN,4241) D2400
    /* //UNUSED const float D2400 = */ read_and_discard( &infile, false );
    //F 603 !
    //F 604       close(lun)
    infile.close();
    //F 605 !
    //F 606 !  ********************************************************************
    //F 607 !
    //F 608 !  READ PARAMETERS FROM MAGXTRA.CFG
    //F 609 !
    //F 610       lun = 42   ! spare logical unit no.
    //F 611       open(unit=lun,file='./magicc_files/MAGXTRA.CFG',status='OLD')
    openfile_read( &infile, "./magxtra_c.cfg", DEBUG_IO );
    //F 612 !
    //F 613         READ(LUN,4240) IOLDTZ
    QADD.IOLDTZ = read_and_discard( &infile, false );
    //F 614         READ(LUN,4241) DT
    CLIM.DT = read_and_discard( &infile, false );
    //F 615         READ(LUN,4240) NOUT
    const int NOUT = read_and_discard( &infile, false );
    //F 616         READ(LUN,4240) IEMPRT
    const int IEMPRT = read_and_discard( &infile, false );
    //F 617         READ(LUN,4240) ICO2PRT
    const int ICO2PRT = read_and_discard( &infile, false );
    //F 618         READ(LUN,4240) ICONCPRT
    const int ICONCPRT = read_and_discard( &infile, false );
    //F 619         READ(LUN,4240) IQGASPRT
    const int IQGASPRT = read_and_discard( &infile, false );
    //F 620         READ(LUN,4241) TOFFSET
    const float TOFFSET = read_and_discard( &infile, false );
    //F 621         READ(LUN,4241) STRATH2O
    METH1.STRATH2O = read_and_discard( &infile, false );
    //F 622         READ(LUN,4241) ENAT
    Sulph.ENAT = read_and_discard( &infile, false );
    //F 623         READ(LUN,4240) IGHG
    JSTART.IGHG = read_and_discard( &infile, false );
    //F 624         READ(LUN,4241) BBCH4
    METH4.BBCH4 = read_and_discard( &infile, false );
    //F 625         READ(LUN,4241) SCH4
    METH3.SCH4 = read_and_discard( &infile, false );
    //F 626         READ(LUN,4241) DELSS
    METH3.DELSS = read_and_discard( &infile, false );
    //F 627         READ(LUN,4241) GAM
    METH4.GAM = read_and_discard( &infile, false );
    //F 628         READ(LUN,4241) ANOX
    METH3.ANOX = read_and_discard( &infile, false );
    //F 629         READ(LUN,4241) DELANOX
    METH3.DELANOX = read_and_discard( &infile, false );
    //F 630         READ(LUN,4241) ACO
    METH3.ACO = read_and_discard( &infile, false );
    //F 631         READ(LUN,4241) DELACO
    METH3.DELACO = read_and_discard( &infile, false );
    //F 632         READ(LUN,4241) AVOC
    METH3.AVOC = read_and_discard( &infile, false );
    //F 633         READ(LUN,4241) DELAVOC
    METH3.DELAVOC = read_and_discard( &infile, false );
    //F 634         READ(LUN,4240) NOOTHER
    const int NOOTHER = read_and_discard( &infile, false );
    //F 635         READ(LUN,4241) BBN2O
    TauNitr.BBN2O = read_and_discard( &infile, false );
    //F 636         READ(LUN,4241) SN2O
    TauNitr.SN2O = read_and_discard( &infile, false );
    //F 637         READ(LUN,4240) NOFFSET
    TauNitr.NOFFSET = read_and_discard( &infile, false );
    //F 638         READ(LUN,4241) WTHRESH
    NSIM.WTHRESH = read_and_discard( &infile, false );
    //F 639         READ(LUN,4241) ES1990
    Sulph.ES1990 = read_and_discard( &infile, false );
    //F 640         READ(LUN,4240) ICCSM
    const int ICCSM = read_and_discard( &infile, false );
    //F 641         IYRQALL=1990
    const int IYRQALL = 1990;
    //F 642 !
    //F 643       close(lun)
    infile.close();
    //F 644 !
    //F 645 !   Call overrite subroutine after each file that may have parameters to overwrite
    //F 646       call overrideParameters( ) !sjs
    overrideParameters( &NEWPARAMS, &CAR, &METH1, &BCOC );
    //F 647 
    //F 648       IF(NOOTHER.EQ.1)THEN
    //F 649         ANOX=0.0
    //F 650         ACO=0.0
    //F 651         AVOC=0.0
    //F 652       ENDIF
    if( NOOTHER == 1 ) METH3.ANOX = METH3.ACO = METH3.AVOC = 0.0;
    //F 653 !
    //F 654 !  ************************************************************
    //F 655 !
    //F 656 !  OPEN MAIN OUTPUT FILE (MAG.OUT).
    //F 657 !
    //F 658    OPEN(UNIT=8,file='./outputs/MAG.CSV',STATUS='UNKNOWN')
    //F 659 ! sjs changed to .csv
    // Handled at beginning of climat()
    //F 660       OPEN(UNIT=88,file='./outputs/CCSM.TXT', STATUS='UNKNOWN')
    ofstream outfile88;
    openfile_write( &outfile88, "./ccsm_c.txt", DEBUG_IO );
    //F 661 !
    //F 662 !  INTERIM CORRECTION TO AVOID CRASH IF S90IND SET TO ZERO IN
    //F 663 !   MAGUSER.CFG
    //F 664 !
    //F 665       IF(S90IND.EQ.0.0)S90IND=-0.0001
    if( Sulph.S90IND == 0.0 ) Sulph.S90IND = -0.0001;
    //F 666 !
    //F 667       XK=YK*3155.76
    CLIM.XK = YK * 3155.76;
    //F 668 !
    //F 669 !  CO2 AND CH4 GAS CYCLE PARAMETERS ARE SELECTED IN DELTAQ.
    //F 670 !   SO4 AEROSOL LEVEL IS SET HERE.
    //F 671 !
    //F 672       IF(LEVSO4.EQ.1)THEN      ! LOW
    if( LEVSO4 == 1 ) {
        //F 673         S90DIR   = -0.2
        //F 674         S90IND   = -0.3
        //F 675         S90BIO   = -0.095
        //F 676         FOC90    =  0.094
        Sulph.S90DIR    = -0.2;
        Sulph.S90IND    = -0.3;
        Sulph.S90BIO    = -0.095;
        Sulph.FOC90     = 0.094;
        //F 677       ENDIF
    }
    //F 678       IF(LEVSO4.EQ.2)THEN      ! MID
    if( LEVSO4 == 2 ) {
        //F 679         S90DIR   = -0.4
        //F 680         S90IND   = -0.7
        //F 681         S90BIO   =  0.025
        //F 682         FOC90    =  0.244
        //F 683       ENDIF
        Sulph.S90DIR    = -0.4;
        Sulph.S90IND    = -0.7;
        Sulph.S90BIO    = 0.025;
        Sulph.FOC90     = 0.244;
    }
    //F 684       IF(LEVSO4.EQ.3)THEN      ! HIGH
    if( LEVSO4 == 3 ) {
        //F 685         S90DIR   = -0.6
        //F 686         S90IND   = -1.1
        //F 687         S90BIO   =  0.145
        //F 688         FOC90    =  0.394
        //F 689       ENDIF
        Sulph.S90DIR    = -0.6;
        Sulph.S90IND    = -1.1;
        Sulph.S90BIO    = -0.025;
        Sulph.FOC90     = 0.244;
    }
    //F 690       IF(LEVSO4.EQ.4)THEN      ! user
    if( LEVSO4 == 4 ) {
        //F 691         S90DIR   = S90Duser
        //F 692         S90IND   = S90Iuser
        //F 693         S90BIO   = S90Buser
        //F 694         FOC90    = FOC90usr
        //F 695       ENDIF
        Sulph.S90DIR    = S90Duser;
        Sulph.S90IND    = S90Iuser;
        Sulph.S90BIO    = S90Buser;
        Sulph.FOC90     = FOC90usr;
    }
    //F 696 !
    //F 697 !  READ IN HALOCARBON FORCINGS. FORCINGS ARE ZERO TO AND INCLUDING
    //F 698 !   1930 (1930-1940 FORCINGS ARE LINEARLY INTERPOLATED FROM ZERO IN
    //F 699 !   193O TO THE 1940 VALUE CALCULATED OFF LINE IN HALOSRES.FOR).
    //F 700 !  INPUT FILE ENDS IN IHALO1.
    //F 701 !  FORCING AFTER 2100 DETERMINED BY IHALOETC (EMS TO ZERO BY 2200 IF
    //F 702 !   IHALOETC=1, CONST AFTER 2100 IF IHALOETC=2).
    //F 703 !  CONST FORCING ASSUMED AFTER IHALO1.
    //F 704 !  FORCINGS ARE END OF YEAR VALUES.
    //F 705 !  QHALOS.IN FORCINGS ARE BROKEN DOWN INTO MONTREAL GASES, STRAT
    //F 706 !   OZONE, MAGICC KYOTO GASES (I.E., THE 8 MAJOR CONTRIBUTORS),
    //F 707 !   AND OTHER GASES (TWO CASES FOR LAST TWO DEPENDING ON IHALOETC).
    //F 708 !  FOR 1991+, QHALOS.IN ONLY GIVES FORCINGS FOR MONTREAL
    //F 709 !   GASES, OTHER GASES AND STRAT OZONE BECAUSE THE CODE
    //F 710 !   CALCULATES THE MAGICC-KYOTO COMPONENT.
    //F 711 !
    //F 712 !  FIRST INITIALIZE QCFC, QMONT, ETC. ARRAYS WITH ZEROES.
    //F 713 !
    //F 714       DO JCFC=0,LASTYEAR-1764
    for( int JCFC=0; JCFC<=(LASTYEAR-1764); JCFC++) {
        //F 715         QCFC(JCFC)=0.0
        FORCE.QCFC[ JCFC ] = 0.0;
        //F 716         QMONT(JCFC)=0.0
        FORCE.QMONT[ JCFC ] = 0.0;
        //F 717         QOTHER(JCFC)=0.0
        FORCE.QOTHER[ JCFC ] = 0.0;
        //F 718         QSTRATOZ(JCFC)=0.0
        FORCE.QSTRATOZ[ JCFC ] = 0.0;
        //F 719         QKYMAG(JCFC)=0.0
        JSTART.QKYMAG[ JCFC ] = 0.0;
        //F 720         CFC12(JCFC)=0.0
        FORCE.CFC12[ JCFC ] = 0.0;
        //F 721       END DO
    }
    //F 722 !
    //F 723       lun = 42   ! spare logical unit no.
    //F 724       open(unit=lun,file='./magicc_files/QHALOS.IN',status='OLD')
    openfile_read( &infile, "./qhalos_c.in", DEBUG_IO );
    //F 725 !
    //F 726       READ(LUN,4446)IHALO1
    int IHALO1 = read_and_discard( &infile, false );
    //F 727       IF(IHALO1.GT.LASTYEAR)IHALO1=LASTYEAR
    if( IHALO1 > LASTYEAR ) IHALO1 = LASTYEAR;
    //F 728 !
    //F 729       READ(LUN,*)
    read_and_discard( &infile, false );
    //F 730       READ(LUN,*)
    read_and_discard( &infile, false );
    //F 731       READ(LUN,*)
    read_and_discard( &infile, false );
    //F 732 !
    //F 733       L0=1930-1765
    const int L0 = 1930-1765;
    //F 734       LAST=IHALO1-1764
    const int LAST = IHALO1 - 1764;
    //F 735 !
    //F 736 !  READ QCFC INPUT ARRAY
    //F 737 !
    //F 738       DO JCFC=L0+1,LAST
    for( int JCFC=L0+1; JCFC<=LAST; JCFC++ ) {
        int IYCFC;
        float JUNK;
        //F 739       IF(IHALOETC.EQ.2)THEN
        if( IHALOETC == 2 ) {
            
            //F 740         READ(LUN,4448)IYCFC,QMONT(JCFC),QSTRATOZ(JCFC),QKYMAG(JCFC), &
            //F 741         QOTHER(JCFC),CFC12(JCFC)
            infile >> IYCFC >> FORCE.QMONT[ JCFC ] >> FORCE.QSTRATOZ[ JCFC ]
            >> JUNK >> JUNK >> JSTART.QKYMAG[ JCFC ] >> FORCE.QOTHER[ JCFC ]
            >> FORCE.CFC12[ JCFC ]; 
            // note skipping two fields--per F3010  4448 FORMAT(1X,I5,2F10.0,20X,3F10.0)
            //F 742       ELSE
        } else {
            //F 743         READ(LUN,4447)IYCFC,QMONT(JCFC),QSTRATOZ(JCFC),QKYMAG(JCFC), &
            //F 744         QOTHER(JCFC),CFC12(JCFC)
            infile >> IYCFC >> FORCE.QMONT[ JCFC ] >> FORCE.QSTRATOZ[ JCFC ]
            >> JSTART.QKYMAG[ JCFC ] >> FORCE.QOTHER[ JCFC ] >> JUNK 
            >> JUNK >> FORCE.CFC12[ JCFC ];
            // note skipping two fields--per F3009  4447 FORMAT(1X,I5,4F10.0,20X,F10.0)
            //F 745       ENDIF
        }
        //F 746       END DO
    }
    //F 747 !
    //F 748       CLOSE(lun)
    infile.close();
    //F 749 !
    //F 750 !  TAU FOR CH4 SOIL SINK CHANGED TO ACCORD WITH IPCC94 (160 yr).
    //F 751 !  SPECIFICATION OF TauSoil MOVED TO MAGEXTRA.CFG ON 1/10/97.
    //F 752 !
    //F 753       IF(ISCENGEN.EQ.1)THEN
    //F 754         NUMSIMS=20
    //F 755       ELSE
    //F 756         NUMSIMS=4
    //F 757       ENDIF
    int NUMSIMS = ( NSIM.ISCENGEN == 1 ) ? 20 : 4;
    //F 758 !
    //F 759 !  SWITCH TO PRODUCE ONLY 1 SIMULATIONS.
    //F 760 !
    //F 761       IF(ISCENGEN.EQ.9)NUMSIMS=1
    if( NSIM.ISCENGEN == 9 ) NUMSIMS = 1;
    //F 762 !
    //F 763       IXLAM=1
    DSENS.IXLAM = 1;
    //F 764       IF(RLO.EQ.1.0) IXLAM=0
    if( RLO == 1 ) DSENS.IXLAM = 0;
    //F 765 !
    //F 766       IF(ICO2READ.GE.1)IMETH=0
    if( CO2READ.ICO2READ > 1 ) METH1.IMETH = 0;
    //F 767 !
    //F 768       IF(DT.GT.1.0)DT=1.0
    if( CLIM.DT > 1.0 ) CLIM.DT = 1.0;
    //F 769 !
    //F 770       QXX=CO2DELQ
    CLIM.QXX = CO2DELQ;
    //F 771       Q2X=QXX*ALOG(2.)
    CLIM.Q2X = CLIM.QXX * log( float( 2.0 ) );
    //F 772 !
    //F 773 !   FO(I) AND FL(I) ARE N.H. AND S.H. OCEAN AND LAND FRACTIONS.
    //F 774 !
    //F 775       FO(1)=1.0-FL(1)
    CLIM.FO[ 1 ] = 1.0 - CLIM.FL[ 1 ];
    //F 776       FO(2)=1.0-FL(2)
    CLIM.FO[ 2 ] = 1.0 - CLIM.FL[ 2 ];
    //F 777 !
    //F 778       FK=RHO*SPECHT*HTCONS/31.5576
    CLIM.FK = CLIM.RHO * CLIM.SPECHT * CLIM.HTCONS / 31.5576;
    //F 779 !
    //F 780 !  TRAP TO CATCH AND OVERWRITE UNREALISTIC D80SIN
    //F 781 !
    //F 782       IF(DUSER.LT.-0.5)THEN
    if( METH1.DUSER < -0.5 ) {
        //F 783         DOLD=DUSER
        float DOLD = METH1.DUSER;
        //F 784         DUSER=-0.5
        METH1.DUSER = -0.5;
        //F 785         WRITE(8,808)DOLD,DUSER
        cout << "*** ERROR : D80SIN SET TOO LOW AT " << DOLD << " : RESET AT " 
        << METH1.DUSER << " ***" << endl;
        //F 786       ENDIF
    }
    //F 787       IF(DUSER.GT.3.0)THEN
    if( METH1.DUSER > 3.0 ) {
        //F 788         DOLD=DUSER
        float DOLD = METH1.DUSER;
        //F 789         DUSER=3.0
        METH1.DUSER = 3.0;
        //F 790         WRITE(8,809)DOLD,DUSER
        cout << "*** ERROR : D80SIN SET TOO HIGH AT " << DOLD << " : RESET AT " 
        << METH1.DUSER << " ***" << endl;
        //F 791       ENDIF
    }
    //F 792 !
    //F 793 !  READ CO2 CONCS DIRECTLY FROM CO2INPUT.DAT IF ICO2READ=1,2,3,4
    //F 794 !
    //F 795       IF(ICO2READ.GE.1.AND.ICO2READ.LE.4)THEN
    if( CO2READ.ICO2READ >= 1 && CO2READ.ICO2READ <= 4 ) {
        //F 796         lun = 42   ! spare logical unit no.
        //F 797         open(unit=lun,file='./magicc_files/Co2input.dat',status='OLD')
        openfile_read( &infile, "./Co2input_c.dat", DEBUG_IO );
        //F 798 !
        //F 799 !  CO2INPUT.DAT MUST HAVE FIRST YEAR = 1990 AND MUST HAVE ANNUAL END
        //F 800 !   OF YEAR VALUES. FIRST LINE OF FILE GIVES LAST YEAR OF ARRAY.
        //F 801 !
        //F 802         READ(lun,900)LCO2
        const int LCO2 = read_and_discard( &infile, false );
        //F 803         ILCO2=LCO2-1764
        const int ILCO2 = LCO2 - 1764;
        //F 804         DO JCO2=226,ILCO2
        for( int JCO2=226; JCO2<=ILCO2; JCO2++ ) {
            //F 805         READ(lun,902)JYEAR,XC(JCO2)
            int JYEAR; infile >> JYEAR >> CO2READ.XC[ JCO2 ];
            //F 806         END DO
        }
        //F 807 !
        //F 808 !  IF LAST YEAR OF INPUT CO2 DATA LESS THAN LASTYEAR FILL OUT
        //F 809 !   ARRAY WITH CONSTANT CO2
        //F 810 !
        //F 811         IF(LASTYEAR.GT.LCO2)THEN
        if( LASTYEAR > LCO2 ) {
            //F 812           DO JCO2=ILCO2+1,LASTYEAR-1764
            for( int JCO2=ILCO2+1; JCO2<=LASTYEAR-1764; JCO2++ )
                //F 813           XC(JCO2)=XC(ILCO2)
                CO2READ.XC[ JCO2 ] = CO2READ.XC[ ILCO2 ];
            //F 814           END DO
            //F 815         ENDIF
        }
        //F 816         close(lun)
        infile.close();
        //F 817       ENDIF
    }
    //F 818 !
    //F 819 !  ************************************************************
    //F 820 !
    //F 821 !  READ EXTRA FORCING IF IQREAD=1 OR 2. IF IQREAD=1, FORCING
    //F 822 !   IN QEXTRA.IN IS ADDED TO ANTHROPOGENIC FORCING. IF IQREAD=2,
    //F 823 !   QEXTRA.IN FORCING IS USED ALONE. QEXTRA.IN HAS A FLAG (NCOLS)
    //F 824 !   TO TELL WHETHER THE DATA ARE GLOBAL (ONE Q COLUMN), HEMISPHERIC
    //F 825 !   (TWO Q COLUMNS, NH THEN SH) OR FOR ALL BOXES (FOUR Q COLUMNS,
    //F 826 !   IN ORDER NHO, NHL, SHO, SHL)
    //F 827 !
    //F 828       IF(IQREAD.GE.1)THEN
    if( QADD.IQREAD >= 1 ) {
        //F 829         lun = 42   ! spare logical unit no.
        //F 830         open(unit=lun,file='./magicc_files/qextra.in',status='OLD')
        openfile_read( &infile, "./qextra_c.in", DEBUG_IO );
        //F 831 !
        //F 832         READ(LUN,900)NCOLS
        //F 833         READ(lun,901)IQFIRST,IQLAST
        NCOLS = read_and_discard( &infile, false );
        infile >> IQFIRST >> IQLAST;
        //F 834         JQFIRST=IQFIRST-1764
        QADD.JQFIRST = IQFIRST - 1764;
        
        //F 835 !
        //F 836 !  TRAP IN CASE FIRST YEAR IS BEFORE 1765
        //F 837 !
        //F 838         IF(JQFIRST.LT.1)THEN
        if( QADD.JQFIRST < 1 ) {
            //F 839           DO JQ=JQFIRST,0
            int JYEAR;
            float QQQGL, QQQNH, QQQSH, QQQNHO, QQQNHL, QQQSHO, QQQSHL;
            for( int JQ=QADD.JQFIRST; JQ<=0; JQ++ ) {
                //F 840           IF(NCOLS.EQ.1)READ(lun,902)JYEAR,QQQGL
                if( NCOLS == 1 ) infile >> JYEAR >> QQQGL;
                //F 841           IF(NCOLS.EQ.2)READ(lun,903)JYEAR,QQQNH,QQQSH
                if( NCOLS == 2 ) infile >> JYEAR >> QQQNH >> QQQSH;
                //F 842           IF(NCOLS.EQ.4)READ(lun,904)JYEAR,QQQNHO,QQQNHL,QQQSHO,QQQSHL
                if( NCOLS == 4) infile >> JYEAR >> QQQNHO >> QQQNHL >> QQQSHO >> QQQSHL;
                //F 843           END DO
            }
            //F 844           JQFIRST=1
            QADD.JQFIRST = 1;
            //F 845           IQFIRST=1765
            IQFIRST = 1765;
            //F 846         ENDIF
        }
        //F 847 !
        //F 848         JQLAST=IQLAST-1764
        QADD.JQLAST = IQLAST - 1764;
        //F 849         DO JQ=JQFIRST,JQLAST
        for( int JQ=QADD.JQFIRST; JQ<=QADD.JQLAST; JQ++ ) {
            //F 850           IF(NCOLS.EQ.1)THEN
            int JYEAR;
            if( NCOLS == 1 ) {
                //F 851             READ(lun,902)JYEAR,QEX(JQ)
                infile >> JYEAR >> QADD.QEX[ JQ ];
                //F 852             QEXNHO(JQ)=(QEX(JQ)-QOFFSET)*QFACTOR
                //F 853             QEXNHL(JQ)=(QEX(JQ)-QOFFSET)*QFACTOR
                //F 854             QEXSHO(JQ)=(QEX(JQ)-QOFFSET)*QFACTOR
                //F 855             QEXSHL(JQ)=(QEX(JQ)-QOFFSET)*QFACTOR
                QADD.QEXNHO[ JQ ] = QADD.QEXNHL[ JQ ] = QADD.QEXSHO[ JQ ] = 
                QADD.QEXSHL[ JQ ] = ( QADD.QEX[ JQ ]-QOFFSET )*QFACTOR;
                //F 856           ENDIF
            }
            //F 857           IF(NCOLS.EQ.2)THEN
            if( NCOLS == 2 ) {
                //F 858             READ(lun,903)JYEAR,QEXNH(JQ),QEXSH(JQ)
                infile >> JYEAR >> QADD.QEXNH[ JQ ] >> QADD.QEXSH[ JQ ];
                //F 859             QEXNHO(JQ)=(QEXNH(JQ)-QOFFSET)*QFACTOR
                //F 860             QEXNHL(JQ)=(QEXNH(JQ)-QOFFSET)*QFACTOR
                QADD.QEXNHO[ JQ ] = QADD.QEXNHL[ JQ ] = ( QADD.QEXNH[ JQ ]-QOFFSET )*QFACTOR;
                //F 861             QEXSHO(JQ)=(QEXSH(JQ)-QOFFSET)*QFACTOR
                //F 862             QEXSHL(JQ)=(QEXSH(JQ)-QOFFSET)*QFACTOR
                QADD.QEXSHO[ JQ ] = QADD.QEXSHL[ JQ ] = ( QADD.QEXSH[ JQ ]-QOFFSET )*QFACTOR;
                //F 863           ENDIF
            }
            //F 864           IF(NCOLS.EQ.4)THEN
            if( NCOLS == 4 ) {
                //F 865             READ(lun,904)JYEAR,QEXNHO(JQ),QEXNHL(JQ),QEXSHO(JQ), &
                //F 866             QEXSHL(JQ)
                infile >> JYEAR >> QADD.QEXNHO[ JQ ] >> QADD.QEXNHL[ JQ ] 
                >> QADD.QEXSHO[ JQ ] >> QADD.QEXSHL[ JQ ];
                //F 867             QEXNHO(JQ)=(QEXNHO(JQ)-QOFFSET)*QFACTOR
                QADD.QEXNHO[ JQ ] = ( QADD.QEXNHO[ JQ ]-QOFFSET )*QFACTOR;
                //F 868             QEXNHL(JQ)=(QEXNHL(JQ)-QOFFSET)*QFACTOR
                QADD.QEXNHL[ JQ ] = ( QADD.QEXNHL[ JQ ]-QOFFSET )*QFACTOR;
                //F 869             QEXSHO(JQ)=(QEXSHO(JQ)-QOFFSET)*QFACTOR
                QADD.QEXSHO[ JQ ] = ( QADD.QEXSHO[ JQ ]-QOFFSET )*QFACTOR;
                //F 870             QEXSHL(JQ)=(QEXSHL(JQ)-QOFFSET)*QFACTOR
                QADD.QEXSHL[ JQ ] = ( QADD.QEXSHL[ JQ ]-QOFFSET )*QFACTOR;
                //F 871           ENDIF
            }
            //F 872         END DO
        }
        //F 873         IF(NCOLS.EQ.1.OR.NCOLS.EQ.2)THEN
        if( NCOLS == 1 || NCOLS == 2 ) {
            //F 874           QEXNH(JQFIRST-1)=QEXNH(JQFIRST)
            QADD.QEXNH[ QADD.JQFIRST-1 ] = QADD.QEXNH[ QADD.JQFIRST ];
            //F 875           QEXSH(JQFIRST-1)=QEXSH(JQFIRST)
            QADD.QEXSH[ QADD.JQFIRST-1 ] = QADD.QEXSH[ QADD.JQFIRST ];
            //F 876         ENDIF
        }
        //F 877         IF(NCOLS.EQ.4)THEN
        if( NCOLS == 4 ) {
            //F 878           QEXNHO(JQFIRST-1)=QEXNHO(JQFIRST)
            QADD.QEXNHO[ QADD.JQFIRST-1 ] = QADD.QEXNHO[ QADD.JQFIRST ];
            //F 879           QEXNHL(JQFIRST-1)=QEXNHL(JQFIRST)
            QADD.QEXNHL[ QADD.JQFIRST-1 ] = QADD.QEXNHL[ QADD.JQFIRST ];
            //F 880           QEXSHO(JQFIRST-1)=QEXSHO(JQFIRST)
            QADD.QEXSHO[ QADD.JQFIRST-1 ] = QADD.QEXSHO[ QADD.JQFIRST ];
            //F 881           QEXSHL(JQFIRST-1)=QEXSHL(JQFIRST)
            QADD.QEXSHL[ QADD.JQFIRST-1 ] = QADD.QEXSHL[ QADD.JQFIRST ];
            //F 882         ENDIF
        }
        //F 883         close(lun)
        infile.close();
        //F 884       ELSE
    } else {
        //F 885         JQLAST=2100-1764
        QADD.JQLAST = 2100 - 1764;
        //F 886         DO JQ=1,JQLAST ! sjs
        for( int JQ=1; JQ<=QADD.JQLAST; JQ++ ) {
            //F 887             QEXNHO(JQ)=0.0
            //F 888             QEXNHL(JQ)=0.0
            //F 889             QEXSHO(JQ)=0.0
            //F 890             QEXSHL(JQ)=0.0
            QADD.QEXNHO[ JQ ] = QADD.QEXNHL[ JQ ] = QADD.QEXSHO[ JQ ] = QADD.QEXSHL[ JQ ] = 0.0;
            //F 891           END DO
        }
        //F 892       ! If QEXTRA was not read-in initialize arrays with zeros in case used for BCOC forcing
        //F 893       ENDIF
    }
    //F 894 
    //F 895 
    //F 896 !
    //F 897 !  ************************************************************
    //F 898 !
    //F 899 !  Read in historical BC and OC emissions and translate into
    //F 900 !  radiative forcing. Add to QXTRA forcing and set QEXTRA to be
    //F 901 !  true if not already. 
    //F 902 !
    //F 903 !  Emissions are read in as global values. sjs
    //F 904 
    //F 905 ! Store original value
    //F 906       OrgIQREAD = IQREAD
    QADD.OrgIQREAD = QADD.IQREAD;
    //F 907 
    //F 908 ! Read in BC emissions if IFOC is set to 3  
    //F 909 ! and if IQREAD is NE 2 (which means only qextra should be used)
    //F 910      IF( IFOC.GE.3 .AND. IQREAD.NE.2 )THEN
    if( Sulph.IFOC >= 3 && QADD.IQREAD != 2 ) {
        //F 911 
        //F 912         lun = 42   ! spare logical unit no.
        //F 913         open(unit=lun,file='./BCOCHist.csv',status='OLD')
        openfile_read( &infile, "./BCOCHist_c.csv", DEBUG_IO );  //FIX location
        //F 914 !
        //F 915         READ(LUN,*)QtempBCUnitForcing, aBCBaseEmissions
        float QtempBCUnitForcing, QtempOCUnitForcing;
        QtempBCUnitForcing = read_csv_value( &infile, false );
        BCOC.aBCBaseEmissions = read_and_discard( &infile, false );
        //F 916         READ(LUN,*)QtempOCUnitForcing, aOCBaseEmissions
        QtempOCUnitForcing = read_csv_value( &infile, false );
        BCOC.aOCBaseEmissions = read_and_discard( &infile, false );
        //F 917 
        //F 918 ! Convert to W/m^2 per Gg
        //F 919         QtempBCUnitForcing = QtempBCUnitForcing / 1000.
        QtempBCUnitForcing /= 1000.0;
        //F 920         QtempOCUnitForcing = QtempOCUnitForcing / 1000.
        QtempOCUnitForcing /= 1000.0;
        //F 921         aOCUnitForcing = aOCUnitForcing / 1000.
        BCOC.aOCUnitForcing /= 1000.0;
        //F 922         aBCUnitForcing = aBCUnitForcing / 1000.
        BCOC.aBCUnitForcing /= 1000.0;
        //F 923         
        //F 924 ! Use default read-in values if have not been otherwise set.
        //F 925         IF(aBCUnitForcing.EQ.0)THEN
        //F 926            aBCUnitForcing   = QtempBCUnitForcing
        //F 927         ENDIF
        if( BCOC.aBCUnitForcing == 0 ) BCOC.aBCUnitForcing = QtempBCUnitForcing;
        //F 928         
        //F 929         IF(aOCUnitForcing.EQ.0)THEN
        //F 930            aOCUnitForcing   = QtempOCUnitForcing
        //F 931         ENDIF
        if( BCOC.aOCUnitForcing == 0 ) BCOC.aOCUnitForcing = QtempOCUnitForcing;
        //F 932 
        //F 933         READ(lun,*)IQFIRST,IQLAST
        int IQFIRST = read_csv_value( &infile, false );
        const int IQLAST = read_and_discard( &infile, false );
        //F 934         JQFIRST=IQFIRST-1764
        QADD.JQFIRST = IQFIRST - 1764;
        //F 935 !
        //F 936 !  TRAP IN CASE FIRST YEAR IS BEFORE 1765 with dummy read
        //F 937 !
        //F 938         IF(JQFIRST.LT.1)THEN
        if( QADD.JQFIRST < 1 ) {
            //F 939           DO JQ=JQFIRST,0
            for( int JQ=QADD.JQFIRST; JQ<=0; JQ++ ) {
                //F 940              READ(lun,902)JYEAR,QQQGL
                read_and_discard( &infile, false );
                //F 941           END DO
            }
            //F 942           JQFIRST=1
            QADD.JQFIRST = 1;
            //F 943           IQFIRST=1765
            IQFIRST = 1765;
            //F 944         ENDIF
        }
        //F 945 !
        //F 946         JQLAST=IQLAST-1764
        QADD.JQLAST = IQLAST - 1764;
        //F 947         DO JQ=JQFIRST,JQLAST
        for( int JQ=QADD.JQFIRST; JQ<=QADD.JQLAST; JQ++ ) {
            //F 948             READ(lun,*)JYEAR,EHistBC(JQ),EHistOC(JQ)
            /* int JYEAR = */ read_csv_value( &infile, false );
            QSPLIT.EHistBC[ JQ ] = read_csv_value( &infile, false );
            QSPLIT.EHistOC[ JQ ] = read_and_discard( &infile, false );
            //F 949         END DO
        }
        //F 950         
        //F 951         close(lun)
        infile.close();
        //F 952         
        //F 953         ! Flag to use QExtra forcing
        //F 954         IQREAD = 1
        QADD.IQREAD = 1;
        //F 955         
        //F 956         !Remove forcing from MAGICC internal calc if BCOC is read in
        //F 957         S90BIO   = 0
        //F 958         FOC90    = 0
        Sulph.S90BIO = Sulph.FOC90 = 0;
        //F 959 
        //F 960       ENDIF
    } // if
    //F 961 !
    //F 962 !  ******************************************************************
    //F 963 !
    //F 964 !  Read in gas emissions from GAS.EMK
    //F 965 !
    //F 966       lun = 42   ! spare logical unit no.
    //F 967 !
    //F 968       open(unit=lun,file='GAS.EMK',status='OLD')
    openfile_read( &infile, "gas.emk", DEBUG_IO );
    //F 969 !
    //F 970 !  READ HEADER AND NUMBER OR ROWS OF EMISIONS DATA FROM GAS.EMK
    //F 971 !
    //F 972       read(lun,4243)  NVAL
    const int NVAL = read_and_discard( &infile, DEBUG_IO );
    //F 973       read(lun,'(a)') mnem
    getline( infile, mnem );
    //F 974       read(lun,*) !   skip description
    skipline( &infile, DEBUG_IO );
    //F 975       read(lun,*) !   skip column headings
    skipline( &infile, DEBUG_IO );
    //F 976       read(lun,*) !   skip units
    skipline( &infile, DEBUG_IO );
    //F 977 !
    //F 978 !  READ INPUT EMISSIONS DATA FROM GAS.EMK
    //F 979 !  SO2 EMISSIONS (BY REGION) MUST BE INPUT AS CHANGES FROM 1990.
    //F 980 !
    //F 981 	iReadNative = 0
    const int iReadNative = 0;
    int ICORR = 0;
    //F 982 
    //F 983 	! Maintain code to read in original magicc input file format
    //F 984     do i=1,NVAL
    for( int i=1; i<=NVAL; i++) {
        //F 985 	  if ( iReadNative .EQ. 1 )THEN
        if( iReadNative == 1 ) {
            //F 986         read(lun,4242) IY1(I),FOS(I),DEF(I),DCH4(I),DN2O(I), &
            IY1[ i ] = read_csv_value( &infile, DEBUG_IO );
            FOS[ i ] = read_csv_value( &infile, DEBUG_IO );
            DEF[ i ] = read_csv_value( &infile, DEBUG_IO );
            DCH4[ i ] = read_csv_value( &infile, DEBUG_IO );
            DN2O[ i ] = read_csv_value( &infile, DEBUG_IO );
            //F 987         DNOX(I),DVOC(I),DCO(I), &
            DNOX[ i ] = read_csv_value( &infile, DEBUG_IO );
            DVOC[ i ] = read_csv_value( &infile, DEBUG_IO );
            DCO[ i ] = read_csv_value( &infile, DEBUG_IO );
            //F 988         DSO21(I),DSO22(I),DSO23(I),DCF4(I),DC2F6(I),D125(I), &
            DSO21[ i ] = read_csv_value( &infile, DEBUG_IO );
            DSO22[ i ] = read_csv_value( &infile, DEBUG_IO );
            DSO23[ i ] = read_csv_value( &infile, DEBUG_IO );
            DCF4[ i ] = read_csv_value( &infile, DEBUG_IO );
            DC2F6[ i ] = read_csv_value( &infile, DEBUG_IO );
            D125[ i ] = read_csv_value( &infile, DEBUG_IO );
            //F 989         D134A(I),D143A(I),D227(I),D245(I),DSF6(I)
            D134A[ i ] = read_csv_value( &infile, DEBUG_IO );
            D143A[ i ] = read_csv_value( &infile, DEBUG_IO );
            D227[ i ] = read_csv_value( &infile, DEBUG_IO );
            D245[ i ] = read_csv_value( &infile, DEBUG_IO );
            DSF6[ i ] = read_and_discard( &infile, DEBUG_IO );
            //F 990   	 END IF
        }
        //F 991 
        //F 992 ! For objects, read in our csv format.
        //F 993 	 IF ( iReadNative .EQ. 0 )THEN
        if( iReadNative == 0 ) {
            //F 994         read(lun,*) IY1(I),FOS(I),DEF(I),DCH4(I),DN2O(I), &
            IY1[ i ] = read_csv_value( &infile, DEBUG_IO );
            FOS[ i ] = read_csv_value( &infile, DEBUG_IO );
            DEF[ i ] = read_csv_value( &infile, DEBUG_IO );
            DCH4[ i ] = read_csv_value( &infile, DEBUG_IO );
            DN2O[ i ] = read_csv_value( &infile, DEBUG_IO );
            //F 995        DSO21(I),DSO22(I),DSO23(I),DCF4(I),DC2F6(I),D125(I), &
            DSO21[ i ] = read_csv_value( &infile, DEBUG_IO );
            DSO22[ i ] = read_csv_value( &infile, DEBUG_IO );
            DSO23[ i ] = read_csv_value( &infile, DEBUG_IO );
            DCF4[ i ] = read_csv_value( &infile, DEBUG_IO );
            DC2F6[ i ] = read_csv_value( &infile, DEBUG_IO );
            D125[ i ] = read_csv_value( &infile, DEBUG_IO );
            //F 996        D134A(I),D143A(I),D227(I),D245(I),DSF6(I), &
            D134A[ i ] = read_csv_value( &infile, DEBUG_IO );
            D143A[ i ] = read_csv_value( &infile, DEBUG_IO );
            D227[ i ] = read_csv_value( &infile, DEBUG_IO );
            D245[ i ] = read_csv_value( &infile, DEBUG_IO );
            DSF6[ i ] = read_csv_value( &infile, DEBUG_IO );
            //F 997        DNOX(I),DVOC(I),DCO(I), DBC(I), DOC(I)  ! Change to match order of writeout -- this is different than magicc default - sjs
            DNOX[ i ] = read_csv_value( &infile, DEBUG_IO );
            DVOC[ i ] = read_csv_value( &infile, DEBUG_IO );
            DCO[ i ] = read_csv_value( &infile, DEBUG_IO );
            DBC[ i ] = read_csv_value( &infile, DEBUG_IO );
            DOC[ i ] = read_and_discard( &infile, DEBUG_IO );
            //F 998   	 END IF
        }
        //F 999  
        //F1000 	IF(i.eq.1) THEN !collect our 1990 values since MAGICC wants differences from 1990
        float DSO211990, DSO221990, DSO231990;
        if( i == 1 ) {
            //F1001 	  ES1990 = DSO21(1) + DSO22(1) + DSO23(1) !global 1990 emissions
            Sulph.ES1990 = DSO21[ 1 ] + DSO22[ 1 ] + DSO23[ 1 ];
            //F1002 	  DSO211990 = DSO21(1) !regional 1990 emissions
            DSO211990 = DSO21[ 1 ];
            //F1003 	  DSO221990 = DSO22(1)
            DSO221990 = DSO22[ 1 ];
            //F1004 	  DSO231990 = DSO23(1)
            DSO231990 = DSO23[ 1 ];
            //F1005 	END IF
        }
        //F1006 
        //F1007 !
        //F1008         DERROR=DPEAK
        float DERROR = DPEAK;
        //F1009         FOS(I)=FOS(I)-DERROR
        FOS[ i ] -= DERROR;
        //F1010         IF(IY1(I).EQ.2000)ICORR=I
        if( IY1[ i ] == 2000 ) ICORR = i;
        //F1011 !
        //F1012 !  ADJUST SO2 EMISSIONS INPUT
        //F1013 !
        //F1014     IF ( iReadNative .EQ. 1 )THEN ! Original MAGICC code
        if( iReadNative == 1 ) {
            //F1015         DSO21(I)= DSO21(I)+ES1990
            DSO21[ i ] += Sulph.ES1990;
            //F1016         DSO22(I)= DSO22(I)+ES1990
            DSO22[ i ] += Sulph.ES1990;
            //F1017         DSO23(I)= DSO23(I)+ES1990
            DSO23[ i ] += Sulph.ES1990;
            //F1018         DSO2(I) = DSO21(I)+DSO22(I)+DSO23(I)-2.0*ES1990
            DSO2[ i ] = DSO21[ i ] + DSO22[ i ] + DSO23[ i ] - 2.0 * Sulph.ES1990;
            //F1019     END IF
        }
        //F1020 !
        //F1021     IF ( iReadNative .EQ. 0 )THEN 
        if( iReadNative == 0 ) {
            //F1022 		DSO21(I) = DSO21(I) - DSO211990 + ES1990 !where ES1990 is the MAGICC global 1990 emissions
            DSO21[ i ] = DSO21[ i ] - DSO211990 + Sulph.ES1990;
            //F1023 		DSO22(I) = DSO22(I) - DSO221990 + ES1990 !note it is added in to all three regions...
            DSO22[ i ] = DSO22[ i ] - DSO221990 + Sulph.ES1990;
            //F1024 		DSO23(I) = DSO23(I) - DSO231990 + ES1990
            DSO23[ i ] = DSO23[ i ] - DSO231990 + Sulph.ES1990;
            //F1025 	
            //F1026 		DSO2(I) = DSO21(I)+DSO22(I)+DSO23(I)-2.0*ES1990 !but here correct for it in the global number
            DSO2[ i ] = DSO21[ i ] + DSO22[ i ] + DSO23[ i ] - 2.0 * Sulph.ES1990;
            //F1027      END IF
        }
        //F1028      
        //F1029      END DO
    } // for
    //F1030      close(lun)
    infile.close();
    // Error checking if year 2000 is not present which is assumed by MAGICC
    if( !ICORR ) {
        cout << "Year 2000 missing from gas.emk." << endl;
        exit( 1 );
    }
    //F1031 !
    //F1032 !  TRAP TO CATCH INCONSISTENCY BETWEEN LASTYEAR FROM CFG FILE
    //F1033 !   AND LAST YEAR OF EMISSIONS INPUT OR MAX ARRAY SIZE
    //F1034 !
    //F1035       IF(ICO2READ.EQ.0)THEN
    if( CO2READ.ICO2READ == 0 ) {
        //F1036         IF((LASTYEAR-1764).GT.iTp)LASTYEAR=iTp+1764
        if( LASTYEAR-1764 > iTp ) LASTYEAR = iTp + 1764;
        //F1037         IF(LASTYEAR.GT.IY1(NVAL)) LASTYEAR=IY1(NVAL)
        if( LASTYEAR > IY1[ NVAL ] ) LASTYEAR = IY1[ NVAL ];
        //F1038       ENDIF
    }
    //F1039       IYEND = LASTYEAR
    const int IYEND = LASTYEAR;
    //F1040       KEND  = IYEND-1764
    Limits.KEND = IYEND - 1764;
    //F1041       KREF  = KYRREF-1764
    const int KREF = KYRREF - 1764;
    //F1042       TEND  = FLOAT(KEND-1)
    CLIM.TEND = Limits.KEND - 1;
    //F1043 !
    //F1044 !  Offset IY1 entries from 1990 : i.e., make IY1(1)=0,
    //F1045 !   IY1(2)=IY1(2)-1990, etc.
    //F1046 !
    //F1047       do i=1,NVAL
    for( int i=1; i<=NVAL; i++)
        //F1048         IY1(i) = IY1(i) - 1990
        IY1[ i ] -= 1990;
    //F1049       end do
    //F1050 !
    //F1051 ! INITIAL (1990) METHANE VALUES: EMISS (LO, MID, HI OR CON) IS THE
    //F1052 !  'CORRECT' 1990 VALUE CALCULATED TO BE CONSISTENT WITH THE
    //F1053 !  CORRESPONDING VALUE OF THE 1990 CH4 TAU. IN GENERAL, EMISS
    //F1054 !  WILL BE INCONSISTENT WITH THE 1990 INPUT VALUE. THIS IS CORRECTED
    //F1055 !  BY OFFSETTING ALL INPUT VALUES BY THE 1990 'ERROR'. SINCE
    //F1056 !  THIS ERROR DEPENDS ON TAU, DIFFERENT OFFSETS MUST BE CALCULATED FOR
    //F1057 !  EACH 1990 TAU VALUE.
    //F1058 ! Note that C and dC/dt must be for mid 1990. VALUES CORRECTED TO
    //F1059 !  AGREE WITH IPCC SAR ON 1/10/98. HISTORY CORRECTED TOO.
    //F1060 !
    //F1061 ! CHANGED TO BALANCE BUDGET IN 2000 FOR TAR (SEPT 2000). BALANCE
    //F1062 !  NUMBERS NOW SPECIFIED IN MAGXTRA.CFG
    //F1063 !
    //F1064       TTLO    = TAUINIT-DELTAU
    const float TTLO = METH3.TAUINIT - METH3.DELTAU;
    //F1065       TTMID   = TAUINIT
    const float TTMID = METH3.TAUINIT;
    //F1066       TTHI    = TAUINIT+DELTAU
    const float TTHI = METH3.TAUINIT + METH3.DELTAU;
    //F1067       TTCON   = TCH4CON
    const float TTCON = METH3.TCH4CON;
    //F1068 !
    //F1069       CBAL=CMBAL
    float CBAL = CMBAL;
    //F1070       DCDT=DCMBAL
    float DCDT = DCMBAL;
    //F1071 !
    //F1072       EMISSLO  = BBCH4*(DCDT +CBAL/TTLO  +CBAL/TAUOTHER)
    float EMISSLO = METH4.BBCH4 * ( DCDT + CBAL/TTLO + CBAL/METH4.TAUOTHER );
    //F1073       EMISSMID = BBCH4*(DCDT +CBAL/TTMID +CBAL/TAUOTHER)
    float EMISSMID = METH4.BBCH4 * ( DCDT + CBAL/TTMID + CBAL/METH4.TAUOTHER );
    //F1074       EMISSHI  = BBCH4*(DCDT +CBAL/TTHI  +CBAL/TAUOTHER)
    float EMISSHI = METH4.BBCH4 * ( DCDT + CBAL/TTHI + CBAL/METH4.TAUOTHER );
    //F1075       EMISSCON = BBCH4*(DCDT +CBAL/TTCON +CBAL/TAUOTHER)
    float EMISSCON = METH4.BBCH4 * ( DCDT + CBAL/TTCON + CBAL/METH4.TAUOTHER );
    //F1076 !
    //F1077 ! INITIAL (2000) N2O VALUE: FOLLOWS CH4 CASE, BUT THERE IS ONLY ONE
    //F1078 !  CORRECTION FACTOR (emissN). THIS IS CALCULATED to be consistent
    //F1079 !  with TN2000. Note that C and dC/dt must be for mid 2000.
    //F1080 !
    //F1081       CBAL=CNBAL
    CBAL = CNBAL;
    //F1082       DCDT=DCNBAL
    DCDT = DCNBAL;
    //F1083 !
    //F1084       EMISSN   = BBN2O*(DCDT +CBAL/TN2000)
    float EMISSN = TauNitr.BBN2O * ( DCDT + CBAL/TauNitr.TN2000 );
    //F1085 !
    //F1086 !  ADD (OR SUBTRACT) CONSTANT TO ALL CH4 AND N2O EMISSIONS TO GIVE
    //F1087 !   1990 VALUE CONSISTENT WITH LIFETIME, CONC AND DC/DT.
    //F1088 !  FOR CH4, ONLY THE CORRECTION FOR THE USER-SPECIFIED 1990 LIFETIME
    //F1089 !   IS APPLIED (GIVEN BY THE CHOICE OF LEVCH4).
    //F1090 !
    //F1091 !  SPECIFY USER LIFETIME. NOTE, THIS IS RESPECIFIED IN DELTAQ.
    //F1092 !
    //F1093       IF(LEVCH4.EQ.1) TTUSER = TTLO
    float TTUSER;
    if( METH2.LEVCH4 == 1 ) TTUSER = TTLO;
    //F1094       IF(LEVCH4.EQ.2) TTUSER = TTMID
    if( METH2.LEVCH4 == 2 ) TTUSER = TTMID;
    //F1095       IF(LEVCH4.EQ.3) TTUSER = TTHI
    if( METH2.LEVCH4 == 3 ) TTUSER = TTHI;
    //F1096       IF(LEVCH4.EQ.4) TTUSER = TTCON
    if( METH2.LEVCH4 == 4 ) TTUSER = TTCON;
    //F1097 !
    //F1098 !  NOTE : D*(2) MUST BE 2000 VALUE
    //F1099 !  THIS IS I=ICORR
    //F1100 !
    //F1101 
    //F1102       CORRMLO  = EMISSLO  - DCH4(ICORR)
    METH1.CORRMLO = EMISSLO - DCH4[ ICORR ];
    //F1103       CORRMMID = EMISSMID - DCH4(ICORR)
    METH1.CORRMMID = EMISSMID - DCH4[ ICORR ];
    //F1104       CORRMHI  = EMISSHI  - DCH4(ICORR)
    METH1.CORRMHI = EMISSHI - DCH4[ ICORR ];
    //F1105       CORRMCON = EMISSCON - DCH4(ICORR)
    float CORRMCON = EMISSCON - DCH4[ ICORR ];
    //F1106 !
    //F1107       CORRN2O  = EMISSN   - DN2O(ICORR)
    float CORRN2O = EMISSN - DN2O[ ICORR ];
    //F1108 !
    //F1109       IF(LEVCH4.EQ.1) CORRUSER = CORRMLO
    if( METH2.LEVCH4 == 1 ) METH1.CORRUSER = METH1.CORRMLO;
    //F1110       IF(LEVCH4.EQ.2) CORRUSER = CORRMMID
    if( METH2.LEVCH4 == 2 ) METH1.CORRUSER = METH1.CORRMMID;
    //F1111       IF(LEVCH4.EQ.3) CORRUSER = CORRMHI
    if( METH2.LEVCH4 == 3 ) METH1.CORRUSER = METH1.CORRMHI;
    //F1112       IF(LEVCH4.EQ.4) CORRUSER = CORRMCON
    if( METH2.LEVCH4 == 4 ) METH1.CORRUSER = CORRMCON;
    //F1113 !
    //F1114       do i=1,NVAL
    for( int i=1; i<=NVAL; i++ ) {
        //F1115         DCH4(I)=DCH4(I)+CORRUSER
        DCH4[ i ] += METH1.CORRUSER;
        //F1116         DN2O(I)=DN2O(I)+CORRN2O
        DN2O[ i ] += CORRN2O;
        //F1117       end do
    }
    //F1118 !
    //F1119 !  ***************************************************************
    //F1120 !
    //F1121       call interp(NVAL,226,IY1,fos,ef)
    interp( NVAL, 226, IY1, FOS, &CARB.EF, Limits.KEND );
    //F1122       call interp(NVAL,226,IY1,def,ednet)
    interp( NVAL, 226, IY1, DEF, &METH1.ednet, Limits.KEND );
    //F1123       call interp(NVAL,226,IY1,DCH4,ECH4)
    interp( NVAL, 226, IY1, DCH4, &CONCS.ECH4, Limits.KEND );
    //F1124       call interp(NVAL,226,IY1,DN2O,EN2O)
    interp( NVAL, 226, IY1, DN2O, &CONCS.EN2O, Limits.KEND );
    //F1125       call interp(NVAL,226,IY1,DNOX,ENOX)
    interp( NVAL, 226, IY1, DNOX, &CONCS.ENOX, Limits.KEND );
    //F1126       call interp(NVAL,226,IY1,DVOC,EVOC)
    interp( NVAL, 226, IY1, DVOC, &CONCS.EVOC, Limits.KEND );
    //F1127       call interp(NVAL,226,IY1,DCO,ECO)
    interp( NVAL, 226, IY1, DCO, &CONCS.ECO, Limits.KEND );
    //F1128       ECO90=ECO(226)
    Sulph.ECO90 = CONCS.ECO.getval( 226 );
    //F1129 !
    //F1130 !  NOTE, IF ESO2 WERE BEING INTERPOLATED, WOULD HAVE TO HAVE ESO2(226)
    //F1131 !   AS LAST ARGUMENT BECAUSE OF MISMATCH OF ESO2 AND Y ARRAYS IN MAIN
    //F1132 !   AND SUBROUTINE INTERP. THIS IS AVOIDED BY USING ESO2SUM.
    //F1133 !
    //F1134       call interp(NVAL,226,IY1,dSO2,ESO2SUM)
    interp( NVAL, 226, IY1, DSO2, &CONCS.ESO2SUM, Limits.KEND );
    //F1135       call interp(NVAL,226,IY1,dSO21,ESO21)
    interp( NVAL, 226, IY1, DSO21, &CONCS.ESO21, Limits.KEND );
    //F1136       call interp(NVAL,226,IY1,dSO22,ESO22)
    interp( NVAL, 226, IY1, DSO22, &CONCS.ESO22, Limits.KEND );
    //F1137       call interp(NVAL,226,IY1,dSO23,ESO23)
    interp( NVAL, 226, IY1, DSO23, &CONCS.ESO23, Limits.KEND );
    //F1138       call interp(NVAL,226,IY1,DCF4 ,ECF4 )
    interp( NVAL, 226, IY1, DCF4, &NEWCONCS.ECF4, Limits.KEND );
    //F1139       call interp(NVAL,226,IY1,DC2F6,EC2F6)
    interp( NVAL, 226, IY1, DC2F6, &NEWCONCS.EC2F6, Limits.KEND );
    //F1140       call interp(NVAL,226,IY1,D125 ,E125 )
    interp( NVAL, 226, IY1, D125, &NEWCONCS.E125, Limits.KEND );
    //F1141       call interp(NVAL,226,IY1,D134A,E134A)
    interp( NVAL, 226, IY1, D134A, &NEWCONCS.E134A, Limits.KEND );
    //F1142       call interp(NVAL,226,IY1,D143A,E143A)
    interp( NVAL, 226, IY1, D143A, &NEWCONCS.E143A, Limits.KEND );
    //F1143       call interp(NVAL,226,IY1,D227 ,E227 )
    interp( NVAL, 226, IY1, D227, &NEWCONCS.E227, Limits.KEND );
    //F1144       call interp(NVAL,226,IY1,D245 ,E245 )
    interp( NVAL, 226, IY1, D245, &NEWCONCS.E245, Limits.KEND );
    //F1145       call interp(NVAL,226,IY1,DSF6 ,ESF6 )
    interp( NVAL, 226, IY1, DSF6, &NEWCONCS.ESF6, Limits.KEND );
    //F1146       call interp(NVAL,226,IY1,DBC  ,EBC  )
    interp( NVAL, 226, IY1, DBC, &CONCS.EBC, Limits.KEND );
    //F1147       call interp(NVAL,226,IY1,DOC  ,EOC  )
    interp( NVAL, 226, IY1, DOC, &CONCS.EOC, Limits.KEND );
    //F1148 !
    //F1149 !  SET ESO2 ARRAY
    //F1150 !
    //F1151       DO KE=226,KEND
    for( int KE=226; KE<=Limits.KEND; KE++)
        //F1152         ESO2(KE)=ESO2SUM(KE)
        CONCS.ESO2.setval( CONCS.ESO2SUM.getval( KE ), KE );
    //F1153       END DO
    //F1154 !
    //F1155 !  FIRST PRINT OUTS TO MAG.OUT
    //F1156 !  PRINT OUT DATE HEADER
    //F1157 !
    //F1158 !      call getdat(myr,imon,iday)  ! sjs - comment out since not available on all platforms
    time_t rawtime;
    time ( &rawtime );
    struct tm * timeinfo = localtime ( &rawtime );    
    //F1159 	myr = 0
    //F1160 	imon = 0
    //F1161 	iday = 0
    //F1162 
    //F1163       write(8,87) mnem,iday,month(imon),myr
    //F1164       write(88,87) mnem,iday,month(imon),myr
    //F1165   87  format(' Emissions profile: ',a20,20x,' Date: ',i2,1x,a3,1x,i4,/)
    outfile8 << " Emissions profile: " << mnem << "  Date: " << asctime( timeinfo ) << endl;
    outfile88 << " Emissions profile: " << mnem << "  Date: " << asctime( timeinfo ) << endl;
    //F1166 !
    //F1167 !  PRINT OUT CO2, CH4 AND SO4 AEROSOL CHOICES (IN WORDS)
    //F1168 !
    //F1169         WRITE(8,110) LEVEL(LEVCO2)
    outfile8 << LEVEL[ CO2READ.LEVCO2 ] <<  " CONCENTRATION PROJECTION FOR CO2" << endl;
    //F1170         IF(IFEED.EQ.0)WRITE(8,1100)
    if( IFEED == 0 ) outfile8 << "(CO2-CLIMATE FEEDBACK NOT INCLUDED)" << endl;
    //F1171         IF(IFEED.EQ.1)WRITE(8,1101)
    if( IFEED == 1 ) outfile8 << "(CO2-CLIMATE FEEDBACK INCLUDED)" << endl;
    //F1172         IF(LEVCO2.EQ.4)WRITE(8,111) DUSER,FUSER
    if( CO2READ.LEVCO2 == 4 ) outfile8 << "Dn(1980s) =" << METH1.DUSER << ": Foc(1980s) =" << METH1.FUSER << endl;
    //F1173         IF(LEVCH4.LE.3) WRITE (8,112) LEVEL(LEVCH4)
    if( METH2.LEVCH4 <= 3 ) outfile8 << LEVEL[ METH2.LEVCH4 ] << " CONCENTRATION PROJECTION FOR CH4" << endl;
    //F1174         IF(LEVCH4.EQ.4) WRITE (8,113) TCH4CON
    if( METH2.LEVCH4 == 4 ) outfile8 << "CH4 CONCS USE CONSTANT LIFETIME OF " << METH3.TCH4CON << endl;
    //F1175         WRITE(8,114) LEVEL(LEVSO4)
    outfile8 << LEVEL[ LEVSO4 ] << " 1990 FORCINGS FOR SO4 AEROSOL" << endl;
    //F1176 !
    //F1177         WRITE(88,110) LEVEL(LEVCO2)
    outfile88 << LEVEL[ CO2READ.LEVCO2 ] <<  " CONCENTRATION PROJECTION FOR CO2" << endl;
    //F1178         IF(IFEED.EQ.0)WRITE(88,1100)
    if( IFEED == 0 ) outfile88 << "(CO2-CLIMATE FEEDBACK NOT INCLUDED)" << endl;
    //F1179         IF(IFEED.EQ.1)WRITE(88,1101)
    if( IFEED == 1 ) outfile88 << "(CO2-CLIMATE FEEDBACK INCLUDED)" << endl;
    //F1180         IF(LEVCO2.EQ.4)WRITE(88,111) DUSER,FUSER
    if( CO2READ.LEVCO2 == 4 ) outfile88 << "Dn(1980s) =" << METH1.DUSER << ": Foc(1980s) =" << METH1.FUSER << endl;
    //F1181         IF(LEVCH4.LE.3) WRITE (88,112) LEVEL(LEVCH4)
    if( METH2.LEVCH4 <= 3 ) outfile88 << LEVEL[ METH2.LEVCH4 ] << " CONCENTRATION PROJECTION FOR CH4" << endl;
    //F1182         IF(LEVCH4.EQ.4) WRITE (88,113) TCH4CON
    if( METH2.LEVCH4 == 4 ) outfile88 << "CH4 CONCS USE CONSTANT LIFETIME OF " << METH3.TCH4CON << endl;
    //F1183 !
    //F1184 !  PRINT OUT HALOCARBON CHOICES
    //F1185 !
    //F1186         IF(IO3FEED.EQ.0) WRITE(8,117)
    if( METH1.IO3FEED == 0 ) outfile8 << "STRAT OZONE DEPLETION FEEDBACK OMITTED" << endl;
    //F1187         IF(IO3FEED.NE.0) WRITE(8,1171)
    if( METH1.IO3FEED != 0 ) outfile8 << "STRAT OZONE DEPLETION FEEDBACK INCLUDED" << endl;
    //F1188         IF(IHALOETC.EQ.2) WRITE(8,1181)
    if( IHALOETC == 2 ) outfile8 << "FOR HALOCARBONS NOT IN GAS.EMK, EMS CONSTANT AFTER 2100" << endl;
    //F1189         IF(IHALOETC.NE.2) WRITE(8,118)
    if( IHALOETC != 2 ) outfile8 << "FOR HALOCARBONS NOT IN GAS.EMK, EMS DROP TO ZERO OVER 2100-2200" << endl;
    //F1190 !
    //F1191 !  PRINT OUT CLIMATE MODEL SELECTED
    //F1192 !
    //F1193         WRITE(8,116) MODNAME(MODEL)
    outfile8 << endl << "CLIMATE MODEL SELECTED = " << MODNAME[ ICE.MODEL ] << endl;
    //F1194 !
    //F1195 !  PRINT OUT ICE MELT SELECTED (LOW, MID OR HIGH)
    //F1196 !
    //F1197         IF(ICE.EQ.1)WRITE(8,1161)
    if( ICE.ICE == 1 ) outfile8 << "USER ICE MELT = LOW" << endl;
    //F1198         IF(ICE.EQ.2)WRITE(8,1162)
    if( ICE.ICE == 2 ) outfile8 << "USER ICE MELT = MID" << endl;
    //F1199         IF(ICE.EQ.3)WRITE(8,1163)
    if( ICE.ICE == 3 ) outfile8 << "USER ICE MELT = HIGH" << endl;
    //F1200 !
    //F1201       NCLIM=1
    NSIM.NCLIM = 1;
    //F1202       CALL INIT
    init( &Limits, &CLIM, &CONCS, &TANDSL, &FORCE, &Sulph, &VARW, &ICE, &AREAS, &NSIM,
         &OZ, &NEWCONCS, &CARB, &CAR, &METH1, &METH2, &METH3, &METH4, &CO2READ, &JSTART,
         &CORREN, &HALOF, &COBS, &TauNitr );
    //F1203 !
    //F1204 !  LINEARLY EXTRAPOLATE LAST ESO2 VALUES FOR ONE YEAR
    //F1205 !
    //F1206       ESO2SUM(KEND+1) = 2.*ESO2SUM(KEND)-ESO2SUM(KEND-1)
    CONCS.ESO2SUM.setval(  2.0 * CONCS.ESO2SUM.getval( Limits.KEND ) - CONCS.ESO2SUM.getval( Limits.KEND-1 ), Limits.KEND+1 );
    //F1207       ESO21(KEND+1)   = 2.*ESO21(KEND)-ESO21(KEND-1)
    CONCS.ESO21.setval(  2.0 * CONCS.ESO21.getval( Limits.KEND ) - CONCS.ESO21.getval( Limits.KEND-1 ), Limits.KEND+1 );
    //F1208       ESO22(KEND+1)   = 2.*ESO22(KEND)-ESO22(KEND-1)
    CONCS.ESO22.setval(  2.0 * CONCS.ESO22.getval( Limits.KEND ) - CONCS.ESO22.getval( Limits.KEND-1 ), Limits.KEND+1 );
    //F1209       ESO23(KEND+1)   = 2.*ESO23(KEND)-ESO23(KEND-1)
    CONCS.ESO23.setval(  2.0 * CONCS.ESO23.getval( Limits.KEND ) - CONCS.ESO23.getval( Limits.KEND-1 ), Limits.KEND+1 );
    //F1210       ESO2(KEND+1)    = ESO2SUM(KEND+1)
    CONCS.ESO2.setval(  CONCS.ESO2SUM.getval( Limits.KEND+1 ), Limits.KEND+1 );
    //F1211 !
    //F1212 !  DEFINE ECO FOR J=1,225
    //F1213 !
    //F1214       DO KC=1,225
    for( int KC=1; KC<=225; KC++ ) {
        //F1215         COE(KC)=ECO90*KC/226
        CONCS.COE.setval( Sulph.ECO90 * KC / 226, KC );
        //F1216       END DO
    }
    //F1217       DO KC=226,KEND
    for( int KC=226; KC<=Limits.KEND; KC++ ) {
        //F1218         COE(KC)=ECO(KC)
        CONCS.COE.setval( CONCS.ECO.getval( KC ), KC );
        //F1219       END DO
    }
    //F1220 !
    //F1221 !  WRITE OUT HEADER INFORMATION FOR MAG.OUT
    //F1222 !
    //F1223 !  SCALING FACTOR FOR CO2 FORCING :
    //F1224 !   (QTOT-Q1990)=(QCO2-QCO2.1990)*CO2SCALE
    //F1225 !
    //F1226       SCAL=100.*(CO2SCALE-1.)
    const float SCAL = 100.0 * ( CO2READ.CO2SCALE - 1.0 );
    //F1227       IF(ICO2READ.EQ.1)WRITE(8,871)SCAL
    if( CO2READ.ICO2READ == 1 ) outfile8 << "CO2 CONC INPUT OVERWRITES CO2 EMISSIONS :"
        << " POST-1990 CO2 FORCING SCALED UP BY" << SCAL << endl;
    //F1228       IF(ICO2READ.EQ.2)WRITE(8,872)
    if( CO2READ.ICO2READ == 2 ) outfile8 << "CO2 CONC INPUT OVERWRITES CO2 EMISSIONS :"
        << " OTHER GAS EMISSIONS AS SPECIFIED IN GAS.EMK" << endl;
    //F1229       IF(ICO2READ.EQ.3)WRITE(8,873)
    if( CO2READ.ICO2READ == 3 ) outfile8 << "CO2 CONC INPUT OVERWRITES CO2 EMISSIONS :"
        << " SO2 EMISSIONS AS SPECIFIED IN GAS.EMK" << endl;
    //F1230 !
    //F1231       IF(IQREAD.EQ.0)WRITE(8,756)
    if( QADD.IQREAD == 0 ) outfile8 << "NO EXTRA FORCING ADDED" << endl;
    //F1232       IF(OrgIQREAD.GE.1)THEN
    if( QADD.OrgIQREAD >= 1 ) {
        //F1233         IF(NCOLS.EQ.1)WRITE(8,757)IQFIRST,IQLAST
        if( NCOLS == 1 ) outfile8 << "EXTRA GLOBAL MEAN FORCING ADDED FROM QEXTRA.IN OVER " 
            << IQFIRST << " TO " << IQLAST << " INCLUSIVE" << endl;
        //F1234         IF(NCOLS.EQ.2)WRITE(8,758)IQFIRST,IQLAST
        if( NCOLS == 2 ) outfile8 << "EXTRA HEMISPHERIC FORCINGS ADDED FROM QEXTRA.IN OVER " 
            << IQFIRST << " TO " << IQLAST << " INCLUSIVE" << endl;
        //F1235         IF(NCOLS.EQ.4)WRITE(8,759)IQFIRST,IQLAST
        if( NCOLS == 4 ) outfile8 << "EXTRA NHO, NHL, SHO, SHL FORCINGS ADDED FROM QEXTRA.IN OVER " 
            << IQFIRST << " TO " << IQLAST << " INCLUSIVE" << endl;
        //F1236         IF(QOFFSET.NE.0.0)WRITE(8,760)QOFFSET
        if( QOFFSET != 0 ) outfile8 << QOFFSET << " W/m**2 SUBTRACTED FROM ALL VALUES" << endl;
        //F1237         IF(QFACTOR.NE.1.0)WRITE(8,761)QFACTOR
        if( QFACTOR != 1 ) outfile8 << "FORCING SCALED BY " << QFACTOR << " AFTER OFFSET" << endl;
        //F1238       ENDIF
    }
    //F1239       IF(IQREAD.EQ.2)WRITE(8,762)
    if( QADD.IQREAD == 2 ) outfile8 << "QEXTRA FORCING USED ALONE" << endl;
    //F1240       
    //F1241       IF ( IFOC.GE.3 ) THEN	! sjs
    if( Sulph.IFOC >= 3 ) {
        //F1242          WRITE(8,*) "ObjECTS Custom BC-OC forcing Used."
        outfile8 << "ObjECTS Custom BC-OC forcing Used." << endl;
        //F1243       ELSE
    } else {
        //F1244          WRITE(8,*) "MAGICC Internal BC-OC forcing Used."
        outfile8 << "MAGICC Internal BC-OC forcing Used." << endl;
        //F1245       ENDIF
    }
    //F1246       
    //F1247 !
    //F1248       WRITE (8,10) Q2X
    outfile8 << endl << "CO2-DOUBLING FORCING IN W/M**2 = " << CLIM.Q2X << endl;
    //F1249 !      WRITE (8,11) FO(1),FO(2),FL(1),FL(2)
    /* outfile8 << "FNHOC= " << CLIM.FO[ 1 ] << " * FSHOC= " << CLIM.FO[ 2 ] 
     << " * FNHLAND= " << CLIM.FL[ 1 ] << " * FSHLAND= " << CLIM.FL[ 2 ] << endl; */
    //F1250 !
    //F1251       IF(S90DIR.EQ.0.0) write(8,*) 'DIRECT AEROSOL FORCING IGNORED'
    if( Sulph.S90DIR == 0.0 ) outfile8 << "DIRECT AEROSOL FORCING IGNORED" << endl;
    //F1252       IF(ABS(S90DIR).GT.0.0) write(8,60)S90DIR
    if( fabs( Sulph.S90DIR ) > 0.0 ) outfile8 << "1990 DIRECT AEROSOL FORCING          =" << Sulph.S90DIR << "W/m**2" << endl;
    //F1253       IF(S90IND.EQ.0.0) write(8,*) 'INDIRECT AEROSOL FORCING IGNORED'
    if( Sulph.S90IND == 0.0 ) outfile8 << "INDIRECT AEROSOL FORCING IGNORED" << endl;
    //F1254       IF(ABS(S90IND).GT.0.0) write(8,61)S90IND
    if( fabs( Sulph.S90IND ) > 0.0 ) outfile8 << "1990 INDIRECT AEROSOL FORCING        =" << Sulph.S90IND << "W/m**2" << endl;
    //F1255       IF(S90BIO.EQ.0.0) write(8,*) 'BIOMASS AEROSOL FORCING IGNORED'
    if( Sulph.S90BIO == 0.0 ) outfile8 << "BIOMASS AEROSOL FORCING IGNORED" << endl;
    //F1256       IF(ABS(S90BIO).GT.0.0) write(8,62)S90BIO
    if( fabs( Sulph.S90BIO ) > 0.0 ) outfile8 << "1990 BIOMASS AEROSOL FORCING         =" << Sulph.S90BIO << "W/m**2" << endl;
    //F1257       IF(FOC90.EQ.0.0) write(8,*) 'FOSSIL OB+BC AEROSOL FORCING', &
    //F1258       ' IGNORED'
    if( Sulph.FOC90 == 0.0 ) outfile8 << "FOSSIL OB+BC AEROSOL FORCING IGNORED" << endl;
    //F1259       IF(ABS(FOC90).GT.0.0) write(8,63)FOC90
    if( fabs( Sulph.FOC90 ) > 0.0 ) outfile8 << "1990 FOSSIL ORG C + BLACK C FORCING  =" << Sulph.FOC90 << "W/m**2" << endl;
    //F1260       write(8,53)STRATH2O
    outfile8 << "STRAT H2O FROM CH4 : DELQH2O/DELQCH4 = " << METH1.STRATH2O << endl;
    //F1261 !
    //F1262 !  ****************************************************************
    //F1263 !
    //F1264 !  Run model NUMSIMS times for different values of DT2X.
    //F1265 !    The input parameter ICEOPT determines what ice melt parameter
    //F1266 !    values are used in each case.
    //F1267 !
    //F1268 !  THE KEY FOR NSIM IS AS FOLLOWS (CASES 17-20 ADDED FEB 7, 1998) ...
    //F1269 !   NOTE : IF ISCENGEN=9, ONLY NSIM=1 IS RUN, BUT NCLIM IS SET TO 4.
    //F1270 !
    //F1271 !  NSIM  CLIM MODEL  EMISSIONS                     NESO2  NCLIM
    //F1272 !     1     LOW         ALL                          1      1
    //F1273 !     2     MID         ALL                          1      2
    //F1274 !     3    HIGH         ALL                          1      3
    //F1275 !     4    USER         ALL                          1      4
    //F1276 !     5     LOW         ESO2 = CONST AFTER 1990      2      1
    //F1277 !     6     MID         ESO2 = CONST AFTER 1990      2      2
    //F1278 !     7    HIGH         ESO2 = CONST AFTER 1990      2      3
    //F1279 !     8    USER         ESO2 = CONST AFTER 1990      2      4
    //F1280 !     9     LOW         ESO2 = ESO2(REGION 1)        3      1
    //F1281 !    10     MID         ESO2 = ESO2(REGION 1)        3      2
    //F1282 !    11    HIGH         ESO2 = ESO2(REGION 1)        3      3
    //F1283 !    12    USER         ESO2 = ESO2(REGION 1)        3      4
    //F1284 !    13     LOW         ESO2 = ESO2(REGION 2)        4      1
    //F1285 !    14     MID         ESO2 = ESO2(REGION 2)        4      2
    //F1286 !    15    HIGH         ESO2 = ESO2(REGION 2)        4      3
    //F1287 !    16    USER         ESO2 = ESO2(REGION 2)        4      4
    //F1288 !    17     LOW         ESO2 = ESO2(REGION 3)        5      1
    //F1289 !    18     MID         ESO2 = ESO2(REGION 3)        5      2
    //F1290 !    19    HIGH         ESO2 = ESO2(REGION 3)        5      3
    //F1291 !    20    USER         ESO2 = ESO2(REGION 3)        5      4
    //F1292 !
    //F1293 !  ALTERNATIVE WAY TO DO REGIONAL BREAKDOWN OF AEROSOL EFFECTS
    //F1294 !
    //F1295 !     9     LOW         ESO2 = ESO2(REG 1+2)         3      1
    //F1296 !    10     MID         ESO2 = ESO2(REG 1+2)         3      2
    //F1297 !    11    HIGH         ESO2 = ESO2(REG 1+2)         3      3
    //F1298 !    12    USER         ESO2 = ESO2(REG 1+2)         3      4
    //F1299 !    13     LOW         ESO2 = ESO2(REG 2+3)         4      1
    //F1300 !    14     MID         ESO2 = ESO2(REG 2+3)         4      2
    //F1301 !    15    HIGH         ESO2 = ESO2(REG 2+3)         4      3
    //F1302 !    16    USER         ESO2 = ESO2(REG 2+3)         4      4
    //F1303 !    17     LOW         ESO2 = ESO2(REG 3+1)         5      1
    //F1304 !    18     MID         ESO2 = ESO2(REG 3+1)         5      2
    //F1305 !    19    HIGH         ESO2 = ESO2(REG 3+1)         5      3
    //F1306 !    20    USER         ESO2 = ESO2(REG 3+1)         5      4
    //F1307 !
    //F1308 !  NOTE : NSIM=5-20 ONLY USED IF ISCENGEN=1 (I.E., USER PLANS TO
    //F1309 !    GO INTO SCENGEN AFTER MAGICC).
    //F1310 !
    //F1311       DO 1 NSIM=1,NUMSIMS
    for( NSIM.NSIM=1; NSIM.NSIM<=NUMSIMS; NSIM.NSIM++ ) {
        //F1312       NESO2=1+INT((NSIM-0.1)/4.0)
        int NESO2 = 1 + int( ( NSIM.NSIM-0.1 ) / 4.0 );
        //F1313       NCLIM=NSIM
        NSIM.NCLIM = NSIM.NSIM;
        //F1314 !
        //F1315 !  RE-SET NCLIM FOR SULPHATE PATTERN WEIGHT CASES (NSIM.GE.5).
        //F1316 !
        //F1317       IF(NSIM.GE.5)NCLIM=NSIM-4
        if( NSIM.NSIM >= 5 ) NSIM.NCLIM = NSIM.NSIM - 4;
        //F1318       IF(NSIM.GE.9)NCLIM=NSIM-8
        if( NSIM.NSIM >= 9 ) NSIM.NCLIM = NSIM.NSIM - 8;
        //F1319       IF(NSIM.GE.13)NCLIM=NSIM-12
        if( NSIM.NSIM >= 13 ) NSIM.NCLIM = NSIM.NSIM - 12;
        //F1320       IF(NSIM.GE.17)NCLIM=NSIM-16
        if( NSIM.NSIM >= 17 ) NSIM.NCLIM = NSIM.NSIM - 16;
        //F1321 !
        //F1322 !  RE-SET NCLIM=4 (USER CASE) IF ONLY ONE SIMULATION (ISCENGEN=9).
        //F1323 !
        //F1324       IF(ISCENGEN.EQ.9)NCLIM=4
        if( NSIM.ISCENGEN == 9 ) NSIM.NCLIM = 4;
        //F1325 !
        //F1326 !  ****************************************************************
        //F1327 !
        //F1328       IF(NESO2.EQ.1)THEN
        if( NESO2 == 1 ) {
            //F1329         DO KE=226,KEND+1
            for( int KE=226; KE<=Limits.KEND+1; KE++ )
                //F1330         ESO2(KE)=ESO2SUM(KE)
                CONCS.ESO2.setval( CONCS.ESO2SUM.getval( KE ), KE);
            //F1331         END DO
            //F1332       ENDIF
        }
        //F1333 !
        //F1334       IF(NESO2.EQ.2)THEN
        if( NESO2 == 2 ) {
            //F1335         DO KE=226,KEND+1
            for( int KE=226; KE<=Limits.KEND+1; KE++ )
                //F1336         ESO2(KE)=ES1990
                CONCS.ESO2.setval( Sulph.ES1990, KE );
            //F1337         END DO
            //F1338       ENDIF
        }
        //F1339 !
        //F1340       IF(NESO2.EQ.3)THEN
        if( NESO2 == 3 ) {
            //F1341         DO KE=226,KEND+1
            for( int KE=226; KE<=Limits.KEND+1; KE++ )
                //F1342         ESO2(KE)=ESO21(KE)
                CONCS.ESO2.setval( CONCS.ESO21.getval( KE ), KE );
            //F1343         END DO
            //F1344       ENDIF
        }
        //F1345 !
        //F1346       IF(NESO2.EQ.4)THEN
        if( NESO2 == 4 ) {
            //F1347         DO KE=226,KEND+1
            for( int KE=226; KE<=Limits.KEND+1; KE++ )
                //F1348         ESO2(KE)=ESO22(KE)
                CONCS.ESO2.setval( CONCS.ESO22.getval( KE ), KE );
            //F1349         END DO
            //F1350       ENDIF
        }
        //F1351 !
        //F1352       IF(NESO2.EQ.5)THEN
        if( NESO2 == 5 ) {
            //F1353         DO KE=226,KEND+1
            for( int KE=226; KE<=Limits.KEND+1; KE++ )
                //F1354         ESO2(KE)=ESO23(KE)
                CONCS.ESO2.setval( CONCS.ESO23.getval( KE ), KE );
            //F1355         END DO
            //F1356       ENDIF
        }
        //F1357 !
        //F1358 !  SET CLIMATE SENSITIVITY
        //F1359 !
        //F1360       IF(NCLIM.EQ.1)TE=1.5        ! LOW
        if( NSIM.NCLIM == 1 ) CLIM.TE = 1.5;
        //F1361       IF(NCLIM.EQ.2)TE=3.0        ! MID
        if( NSIM.NCLIM == 2 ) CLIM.TE = 3.0;
        //F1362       IF(NCLIM.EQ.3)TE=6.0        ! HIGH
        if( NSIM.NCLIM == 3 ) CLIM.TE = 6.0;
        //F1363       IF(NCLIM.EQ.4)TE=DT2XUSER   ! USER
        if( NSIM.NCLIM == 4 ) CLIM.TE = NEWPARAMS.DT2XUSER;
        //F1364 !
        //F1365       IF(IXLAM.EQ.1)THEN
        if( DSENS.IXLAM == 1 )
            //F1366         CALL LAMCALC(Q2X,FL(1),FL(2),XKLO,XKNS,TE,RLO,XLAMO,XLAML)
            lamcalc( CLIM.Q2X, CLIM.FL[ 1 ], CLIM.FL[ 2 ], CLIM.XKLO, CLIM.XKNS, CLIM.TE, RLO, &DSENS.XLAMO, &DSENS.XLAML );
        //F1367       ENDIF
        //F1368 !
        //F1369       XLAM=Q2X/TE
        CLIM.XLAM - CLIM.Q2X / CLIM.TE;
        //F1370 !
        //F1371       CALL INIT
        init( &Limits, &CLIM, &CONCS, &TANDSL, &FORCE, 
             &Sulph, &VARW, &ICE, &AREAS, &NSIM,
             &OZ, &NEWCONCS, &CARB, &CAR, &METH1,
             &METH2, &METH3, &METH4, &CO2READ, &JSTART,
             &CORREN, &HALOF, &COBS, &TauNitr );     
        //F1372 !
        //F1373       IF(NESO2.EQ.1)THEN
        if( NESO2 == 1 ) {
            //F1374         WRITE (8,179)
            outfile8 << endl << "*************************************************" << endl << endl;
            //F1375         WRITE (8,176) NSIM,TE
            outfile8 << "NSIM = " << NSIM.NSIM << ": DELT(2XCO2) = " << CLIM.TE << "DEGC" << endl;
            //F1376         WRITE (88,1761) TE
            outfile88 << "DELT(2XCO2) = " << CLIM.TE << " DEGC" << endl;
            //F1377         IF(NESO2.EQ.1)WRITE(8,186)
            if( NESO2 == 1 ) outfile8 << "FULL GLOBAL SO2 EMISSIONS" << endl;
            //F1378         IF(NESO2.EQ.2)WRITE(8,187)
            if( NESO2 == 2 ) outfile8 << "SO2 EMISSIONS CONSTANT AFTER 1990" << endl;
            //F1379         IF(NESO2.EQ.3)WRITE(8,188)
            if( NESO2 == 3 ) outfile8 << "REGION 1 SO2 EMISSIONS" << endl;
            //F1380         IF(NESO2.EQ.4)WRITE(8,189)
            if( NESO2 == 4 ) outfile8 << "REGION 2 SO2 EMISSIONS" << endl;
            //F1381         IF(NESO2.EQ.5)WRITE(8,190)
            if( NESO2 == 5 ) outfile8 << "REGION 3 SO2 EMISSIONS" << endl;
            //F1382         WRITE(8,1220)IVARW
            outfile8 << endl << "IVARW SET AT " << VARW.IVARW << endl;
            //F1383         IF(IVARW.EQ.0)THEN
            if( VARW.IVARW == 0 ) {
                //F1384           WRITE (8,122)
                outfile8 << "CONSTANT W CASE" << endl;
                //F1385         ENDIF
            }
            //F1386         IF(IVARW.EQ.1)THEN
            if( VARW.IVARW == 1 ) {
                //F1387           WRITE (8,123) TW0NH
                outfile8 << "VARIABLE W : NH W = ZERO WHEN TEMPERATURE = " << VARW.TW0NH << " degC" << endl;
                //F1388           WRITE (8,124) TW0SH
                outfile8 << "VARIABLE W : SH W = ZERO WHEN TEMPERATURE = " << VARW.TW0SH << " degC" << endl;
                //F1389           IF(KEYDW.EQ.1)WRITE (8,1231)
                if( VARW.KEYDW == 1 ) outfile8 << "FULL W SCALED WITH GLOBAL-MEAN TEMPERATURE" << endl;
                //F1390           IF(KEYDW.EQ.2)WRITE (8,1232)
                if( VARW.KEYDW == 2 ) outfile8 << "FULL W SCALED WITH GLOBAL-MEAN OCEAN TEMPERATURE" << endl;
                //F1391           IF(KEYDW.EQ.3)WRITE (8,1233)
                if( VARW.KEYDW == 3 ) outfile8 << "FULL W SCALED WITH HEMISPHERIC-MEAN OCEAN TEMPERATURE" << endl;
                //F1392           IF(KEYDW.EQ.4)WRITE (8,1234)
                if( VARW.KEYDW == 4 ) outfile8 << "ACTIVE W SCALED WITH GLOBAL-MEAN TEMPERATURE" << endl;
                //F1393           IF(KEYDW.EQ.5)WRITE (8,1235)
                if( VARW.KEYDW == 5 ) outfile8 << "ACTIVE W SCALED WITH GLOBAL-MEAN OCEAN TEMPERATURE" << endl;
                //F1394         ENDIF
            }
            //F1395         IF(IVARW.EQ.2)THEN
            if( VARW.IVARW == 2 ) {
                //F1396           WRITE (8,126) WTHRESH
                outfile8 << "PERMANENT THC SHUTDOWN AT W = " << NSIM.WTHRESH << " M/YR" << endl;
                //F1397           WRITE (8,127) TW0NH
                outfile8 << "W = ZERO WHEN TEMPERATURE = " << VARW.TW0NH << "degC" << endl;
                //F1398           IF(KEYDW.EQ.1)WRITE (8,1231)
                if( VARW.KEYDW == 1 ) outfile8 << "FULL W SCALED WITH GLOBAL-MEAN TEMPERATURE" << endl;
                //F1399           IF(KEYDW.EQ.2)WRITE (8,1232)
                if( VARW.KEYDW == 2 ) outfile8 << "FULL W SCALED WITH GLOBAL-MEAN OCEAN TEMPERATURE" << endl;
                //F1400           IF(KEYDW.EQ.3)WRITE (8,1233)
                if( VARW.KEYDW == 3 ) outfile8 << "FULL W SCALED WITH HEMISPHERIC-MEAN OCEAN TEMPERATURE" << endl;
                //F1401           IF(KEYDW.EQ.4)WRITE (8,1234)
                if( VARW.KEYDW == 4 ) outfile8 << "ACTIVE W SCALED WITH GLOBAL-MEAN TEMPERATURE" << endl;
                //F1402           IF(KEYDW.EQ.5)WRITE (8,1235)
                if( VARW.KEYDW == 5 ) outfile8 << "ACTIVE W SCALED WITH GLOBAL-MEAN OCEAN TEMPERATURE" << endl;
                //F1403         ENDIF
            }
            //F1404         IF(IVARW.EQ.3)THEN
            if( VARW.IVARW == 3 ) {
                //F1405           WRITE (8,125)
                outfile8 << "VARIABLE W : NH AND SH W(t) SPECIFIED IN WINPUT.IN" << endl;
                //F1406         ENDIF
            }
            //F1407 !
            //F1408         WRITE (8,12) XKNS,XKLO
            outfile8 << endl << "XKNS = " << CLIM.XKNS << " : XKLO = " << CLIM.XKLO << endl;
            //F1409         WRITE (8,120) HM,YK
            outfile8 << "HM = " << CLIM.HM << "M : XK = " << /*CLIM.*/YK << " CM**2/SEC" << endl;  //FIX this is a MAGICC bug
            //F1410         WRITE (8,121) PI,W0
            outfile8 << "PI = " << CLIM.PI << " : INITIAL W= " << CLIM.W0 << " M/YR" << endl << endl << endl;
            //F1411 !
            //F1412         IF(IXLAM.EQ.1)THEN
            if( DSENS.IXLAM == 1 ) {
                //F1413           WRITE(8,914) RLO,XLAML,XLAMO
                outfile8 << "DIFF/L SENSITIVITY CASE : RLO = " << RLO << " : XLAML = " << DSENS.XLAML << " : XLAMO = " << DSENS.XLAMO << endl;
                //F1414           IF(XLAML.LT.0.0)WRITE(8,916)
                if( DSENS.XLAML < 0.0 )
                    outfile8 << "  **  WARNING, XLAML<0.0 : USE SMALLER XKLO **" << endl;
                //F1415         ELSE
                //F1416           WRITE(8,915) XLAM
            }
            else 
                outfile8 << "GLOBAL SENSITIVITY CASE : INITIAL XLAM =" << CLIM.XLAM << endl;
            //F1417         ENDIF
            //F1418       ENDIF
        }
        //F1419 !
        //F1420 !  ***********************************************************
        //F1421 !
        //F1422       CALL RUNMOD
        runmod( &Limits, &CLIM, &CONCS, &CARB, &TANDSL,
               &CO2READ, &Sulph, &DSENS, &VARW, &QSPLIT,
               &AREAS, &QADD, &BCOC, &FORCE, &NSIM,
               &OZ, &NEWCONCS, &CAR, &METH1, &METH2, &METH3, &METH4, &TauNitr,
               &JSTART, &CORREN, &HALOF, &COBS, &ICE, &outfile8 );
        //F1423 !
        //F1424 !  EXTRA CALL TO RUNMOD TO GET FINAL FORCING VALUES FOR K=KEND
        //F1425 !   WHEN DT=1.0
        //F1426 !
        //F1427 !      IF((K.EQ.KEND).AND.(DT.EQ.1.0))CALL RUNMOD
        //F1428 !
        //F1429 !  SAVE SULPHATE AEROSOL FORCINGS IN FIRST PASS THROUGH OF NSIM
        //F1430 !   LOOP, WHEN TOTAL SO2 EMISSIONS ARE BEING USED.
        //F1431 !
        //F1432       IF(NSIM.EQ.1)THEN
        if( NSIM.NSIM == 1 ) {
            //F1433         DO K=1,KEND
            for( int K=1; K<=Limits.KEND; K++ ) {
                //F1434           QSO2SAVE(K)=QSO2(K)
                STOREDVALS.QSO2SAVE[ K ] = TANDSL.QSO2[ K ];
                //F1435           QDIRSAVE(K)=QDIR(K)
                STOREDVALS.QDIRSAVE[ K ] = TANDSL.QDIR[ K ];
                //F1436         END DO
            }
            //F1437       ENDIF
        }
        //F1438 !
        //F1439 !  PRINT OUT RESULTS
        //F1440 !
        //F1441       DT1    = TGAV(226)-TGAV(116)
        const float DT1 = TANDSL.TGAV[ 226 ] - TANDSL.TGAV[ 116 ];
        //F1442       DMSL1  = SLT(226)-SLT(116)
        const float DMSL1 = TANDSL.SLT[ 226 ] - TANDSL.SLT[ 116 ];
        //F1443       DTNH1  = TNHAV(226)-TNHAV(116)
        const float DTNH1 = TANDSL.TNHAV[ 226 ] - TANDSL.TNHAV[ 116 ];
        //F1444       DTSH1  = TSHAV(226)-TSHAV(116)
        const float DTSH1 = TANDSL.TSHAV[ 226 ] - TANDSL.TSHAV[ 116 ];
        //F1445       DTLAND = TLAND(226)-TLAND(116)
        const float DTLAND = TANDSL.TLAND[ 226 ] - TANDSL.TLAND[ 116 ];
        //F1446       DTOCEAN= TOCEAN(226)-TOCEAN(116)
        const float DTOCEAN = TANDSL.TOCEAN[ 226 ] - TANDSL.TOCEAN[ 116 ];
        //F1447       DTNHO  = TNHO(226)-TNHO(116)
        const float DTNHO = TANDSL.TNHO[ 226 ] - TANDSL.TNHO[ 116 ];
        //F1448       DTSHO  = TSHO(226)-TSHO(116)
        const float DTSHO = TANDSL.TSHO[ 226 ] - TANDSL.TSHO[ 116 ];
        //F1449       DTNHL  = TNHL(226)-TNHL(116)
        const float DTNHL = TANDSL.TNHL[ 226 ] - TANDSL.TNHL[ 116 ];
        //F1450       DTSHL  = TSHL(226)-TSHL(116)
        const float DTSHL = TANDSL.TSHL[ 226 ] - TANDSL.TSHL[ 116 ];
        //F1451 !
        //F1452       IF(NESO2.EQ.1)THEN
        if( NESO2 == 1 ) {
            //F1453         WRITE (8,140) DT1,DMSL1
            outfile8 << endl << "1880-1990 CHANGES : GLOBAL DTEMP = " << DT1 << " :   DMSL = " << DMSL1 << endl;
            //F1454         WRITE (8,141) DTNHL,DTNHO,DTSHL,DTSHO
            outfile8 << "          DTNHL = " << DTNHL << " : DTNHO = " << DTNHO << " :  DTSHL = " << DTSHL << " :   DTSHO = " << DTSHO << endl;
            //F1455         WRITE (8,142) DTNH1,DTSH1,DTLAND,DTOCEAN
            outfile8 << "          DTNH = " << DTNH1 << " :  DTSH = " << DTSH1 << " : DTLAND = " << DTLAND << " : DTOCEAN =" << DTOCEAN << endl;
            //F1456         WRITE (8,15) KYRREF
            outfile8 << endl << "** TEMPERATURE AND SEA LEVEL CHANGES FROM " << KYRREF << endl;
            //F1457         WRITE (8,16)
            outfile8 << "     (FIRST LINE GIVES 1765-1990 CHANGES : ALL VALUES ARE MID-YEAR TO MID-YEAR)" << endl << endl;
            //F1458 !
            //F1459         IF(IVARW.GE.1)THEN
            if( VARW.IVARW >= 1 ) {
                //F1460           WRITE(8,178)TE
                outfile8 << "DT2X = " << CLIM.TE << " : VARIABLE W" << endl << endl;
                //F1461         ELSE
            } else {
                //F1462           WRITE(8,177)TE
                outfile8 << "DT2X = " << CLIM.TE << " : CONSTANT W" << endl << endl;
                //F1463         ENDIF
            }
            //F1464 !
            //F1465         IF(NCLIM.EQ.1)WRITE(8,161)
            if( NSIM.NCLIM == 1 ) outfile8 << "LOW CLIMATE AND SEA LEVEL MODEL PARAMETERS" << endl;
            //F1466         IF(NCLIM.EQ.2)WRITE(8,162)
            if( NSIM.NCLIM == 2 ) outfile8 << "MID CLIMATE AND SEA LEVEL MODEL PARAMETERS" << endl;
            //F1467         IF(NCLIM.EQ.3)WRITE(8,163)
            if( NSIM.NCLIM == 3 ) outfile8 << "HIGH CLIMATE AND SEA LEVEL MODEL PARAMETERS" << endl;
            //F1468         IF(NCLIM.EQ.4)WRITE(8,164)
            if( NSIM.NCLIM == 4 ) outfile8 << "USER CLIMATE AND SEA LEVEL MODEL PARAMETERS" << endl;
            //F1469 !
            //F1470         IF(NOUT.EQ.1)WRITE(8,171)
            if( NOUT == 1 ) outfile8 << " YEAR,DELTAQ, TEQU, TEMP, EXPN, GLAC, GREENL,ANTAR,Z-XTRA,MSLTOT, TNH, TSH, WNH, WSH,YEAR, GR+ANT,ZTOT-ZXTRA" << endl;
            //F1471         IF(NOUT.EQ.2)WRITE(8,172)
            if( NOUT == 2 ) outfile8 << " YEAR,DELTAQ, TEQU, TEMP, EXPN, GLAC, GREENL,ANTAR, Z-XTRA, MSLTOT, TLAND,  TOCN, TL/TO,   WNH,   WSH, YEAR" << endl;
            //F1472         IF(NOUT.EQ.3)WRITE(8,173)
            if( NOUT == 3 ) outfile8 << "YEAR,DELTAQ, TEQU, TEMP, EXPN, GLAC, GREENL,  ANTAR, Z-XTRA, MSLTOT, TEQ-T, TDEEP,   WNH,   WSH,  YEAR" << endl;
            //F1473         IF(NOUT.EQ.4)WRITE(8,174)
            if( NOUT == 4 ) outfile8 << " YEAR,EQVCO2, TEQU, TEMP, EXPN, GLAC, GREENL,ANTAR.Z-XTRA.MSLTOT.TEQ-T.TDEEP, WNH, WSH,YEAR" << endl;
            //F1474         IF(NOUT.EQ.5)THEN
            //F1475           WRITE(8,175)
            if( NOUT == 5 ) outfile8 << " YEAR,DELTAQ,   TEMP,TL/TO,MSLTOT, EXPN, GLAC, GREENL,ANTAR Z-XTRA, WNH,YEAR" << endl;
            //F1476         ENDIF
            //F1477 !
            //F1478 !  PRINTOUT OPTIONS
            //F1479 !
            //F1480         IF(NOUT.EQ.1)THEN
            float COL9, COL10, COL11;
            if( NOUT == 1 ) {
                //F1481           COL9   = TNHAV(226)
                COL9 = TANDSL.TNHAV[ 226 ];
                //F1482           COL10  = TSHAV(226)
                COL10 = TANDSL.TSHAV[ 226 ];
                //F1483         ENDIF
            }
            //F1484         IF(NOUT.EQ.2.OR.NOUT.EQ.5)THEN
            if( NOUT == 2 || NOUT == 5 ) {
                //F1485           COL9   = TLAND(226)
                COL9 = TANDSL.TLAND[ 226 ];
                //F1486           COL10  = TOCEAN(226)
                COL10 = TANDSL.TOCEAN[ 226 ];
                //F1487           IF(COL10.NE.0.0)THEN
                if( COL10 != 0.0 ) {
                    //F1488             COL11= COL9/COL10
                    COL11 = COL9 / COL10;
                    //F1489           ELSE
                } else {
                    //F1490             COL11= 9.999
                    COL11 = 9.999;
                    //F1491           ENDIF
                }
                //F1492         ENDIF
            }
            //F1493         IF(NOUT.EQ.3.OR.NOUT.EQ.4)THEN
            if( NOUT == 3 || NOUT == 4 ) {
                //F1494           COL9   = TEQU(226)-TGAV(226)
                COL9 = TANDSL.TEQU[ 226 ] - TANDSL.TGAV[ 226 ];
                //F1495           COL10  = TDEEP(226)
                COL10 = TANDSL.TDEEP[ 226 ];
                //F1496         ENDIF
            }
            //F1497 !
            //F1498         TEOUT=TEQU(226)
            const float TEOUT = TANDSL.TEQU[ 226 ];
            //F1499 !
            //F1500 !  THE METHOD FOR CALCULATING SLT SEEMS TO BE WRONG. IT ATTEMPTS TO
            //F1501 !   USE QUADRATURE FOR THE ERRORS IN INDIVIDUAL COMPONENTS, BUT DOES
            //F1502 !   NOT DO THIS CONSISTENTLY. SO IT IS BETTER TO SUM THE INDIVIDUAL
            //F1503 !   VALUES.
            //F1504 !
            //F1505         SLT(226)=EX(226)+SLI(226)+SLG(226)+SLA(226)+SLO(226)   
            TANDSL.SLT[ 226 ] = TANDSL.EX[ 226 ] + TANDSL.SLI[ 226 ] + TANDSL.SLG[ 226 ] + TANDSL.SLA[ 226 ] + TANDSL.SLO[ 226 ];
            //F1506         IF(NOUT.EQ.1)THEN
            if( NOUT == 1 ) {
                //F1507           SLICE=SLA(226)+SLG(226)
                const float SLICE = TANDSL.SLA[ 226 ] + TANDSL.SLG[ 226 ];
                //F1508           SLRAW=SLT(226)-SLO(226)
                const float SLRAW = TANDSL.SLT[ 226 ] - TANDSL.SLO[ 226 ];
                //F1509           WRITE(8,181)QGLOBE(226),TEOUT,TGAV(226),EX(226),SLI(226), &
                //F1510           SLG(226),SLA(226),SLO(226),SLT(226),COL9,COL10,WNH(226), &
                //F1511           WSH(226),SLICE,SLRAW
                outfile8 << "TO1990," << QSPLIT.QGLOBE[ 226 ] << "," << TEOUT << "," << TANDSL.TGAV[ 226 ] << ","
                << TANDSL.EX[ 226 ] << "," << TANDSL.SLI[ 226 ] << "," << TANDSL.SLG[ 226 ] << "," 
                << TANDSL.SLA[ 226 ] << "," << TANDSL.SLO[ 226 ] << "," << TANDSL.SLT[ 226 ] << "," << COL9 << "," 
                << COL10 << "," << VARW.WNH[ 226 ] << "," << VARW.WSH[ 226 ] << "," << SLICE << "," << SLRAW << endl;
                //F1512         ENDIF
            }
            //F1513         IF(NOUT.EQ.2) THEN
            if( NOUT == 2 ) {
                //F1514           WRITE(8,182)QGLOBE(226),TEOUT,TGAV(226),EX(226),SLI(226), &
                //F1515           SLG(226),SLA(226),SLO(226),SLT(226),COL9,COL10,COL11,WNH(226), &
                //F1516           WSH(226)
                outfile8 << "TO1990," << QSPLIT.QGLOBE[ 226 ] << "," << TEOUT << "," << TANDSL.TGAV[ 226 ] << ","
                << TANDSL.EX[ 226 ] << "," << TANDSL.SLI[ 226 ] << "," << TANDSL.SLG[ 226 ] << "," 
                << TANDSL.SLA[ 226 ] << "," << TANDSL.SLO[ 226 ] << "," << TANDSL.SLT[ 226 ] << "," << COL9 << "," 
                << COL10 << "," << COL11 << "," << VARW.WNH[ 226 ] << "," << VARW.WSH[ 226 ] << endl;
                //F1517         ENDIF
            }
            //F1518         IF(NOUT.EQ.3)THEN
            if( NOUT == 3 ) {
                //F1519           WRITE(8,183)QGLOBE(226),TEOUT,TGAV(226),EX(226),SLI(226), &
                //F1520           SLG(226),SLA(226),SLO(226),SLT(226),COL9,COL10,WNH(226), &
                //F1521           WSH(226)
                outfile8 << "TO1990," << QSPLIT.QGLOBE[ 226 ] << "," << TEOUT << "," << TANDSL.TGAV[ 226 ] << ","
                << TANDSL.EX[ 226 ] << "," << TANDSL.SLI[ 226 ] << "," << TANDSL.SLG[ 226 ] << "," 
                << TANDSL.SLA[ 226 ] << "," << TANDSL.SLO[ 226 ] << "," << TANDSL.SLT[ 226 ] << "," << COL9 << "," 
                << COL10 << "," << VARW.WNH[ 226 ] << "," << VARW.WSH[ 226 ] << endl;
                //F1522         ENDIF
            }
            //F1523         IF(NOUT.EQ.4)THEN
            if( NOUT == 4 ) {
                //F1524 !
                //F1525 !  CONVERT W/M**2 TO EQUIV CO2 RELATIVE TO END-1765 CO2 CONC
                //F1526 !
                //F1527           EQUIVCO2=COBS(1)*EXP(QGLOBE(226)/QXX)
                EQUIVCO2 = COBS.COBS[ 1 ] * exp( float( QSPLIT.QGLOBE[ 226 ] / CLIM.QXX ) );
                //F1528           WRITE(8,184)EQUIVCO2,TEOUT,TGAV(226),EX(226),SLI(226), &
                //F1529           SLG(226),SLA(226),SLO(226),SLT(226),COL9,COL10,WNH(226), &
                //F1530           WSH(226)
                outfile8 << "TO1990," << EQUIVCO2 << "," << TEOUT << "," << TANDSL.TGAV[ 226 ] << ","
                << TANDSL.EX[ 226 ] << "," << TANDSL.SLI[ 226 ] << "," << TANDSL.SLG[ 226 ] << "," 
                << TANDSL.SLA[ 226 ] << "," << TANDSL.SLO[ 226 ] << "," << TANDSL.SLT[ 226 ] << "," << COL9 << "," 
                << COL10 << "," << VARW.WNH[ 226 ] << "," << VARW.WSH[ 226 ] << endl;
                //F1531         ENDIF
            }
            //F1532         IF(NOUT.EQ.5)THEN
            if( NOUT == 5 ) {
                //F1533           WRITE(8,185)QGLOBE(226),TGAV(226),COL11,SLT(226),EX(226), &
                //F1534           SLI(226),SLG(226),SLA(226),SLO(226),WNH(226)
                outfile8 << "TO1990," << QSPLIT.QGLOBE[ 226 ] << "," << TANDSL.TGAV[ 226 ] << "," << COL11 << ","
                << TANDSL.SLT[ 226 ] << "," << TANDSL.SLG[ 226 ] << "," << TANDSL.SLA[ 226 ] << "," << TANDSL.SLO[ 226 ] 
                << VARW.WNH[ 226 ] << endl;
                //F1535         ENDIF
            }
            //F1536 !
            //F1537 !  ****************************************************************
            //F1538 !
            //F1539 !  MAIN PRINT OUT LOOP
            //F1540 !
            //F1541         TE1 = 0.
            //F1542         TT1 = 0.
            //F1543         TN1 = 0.
            //F1544         TS1 = 0.
            //UNUSED float TE1 = 0.0, TT1 = 0.0, TN1 = 0.0, TS1 = 0.0;
            //F1545 !
            //F1546         QR    = QGLOBE(KREF)
            const float QR = QSPLIT.QGLOBE[ KREF ];
            
            //F1547         XPRT  = FLOAT(ITEMPRT)
            const float XPRT = float( ITEMPRT );
            //F1548         NPRT  = INT(225./XPRT +0.01)
            const int NPRT = int( 225.0 / XPRT + 0.01 );
            //F1549         MPRT  = NPRT*ITEMPRT
            const int MPRT = NPRT * ITEMPRT;
            //F1550         KYEAR0= 1990-MPRT
            const int KYEAR0 = 1990 - MPRT;
            //F1551 !
            //F1552         TREFSUM=0.0
            float TREFSUM = 0.0;
            //F1553 !
            //F1554         DO 987 K=1,KEND
            for( int K=1; K<=Limits.KEND; K++ ) {
                //F1555 !
                //F1556           KYEAR=1764+K
                int KYEAR = 1764 + K;
                //F1557           QK = QGLOBE(K)
                float QK = QSPLIT.QGLOBE[ K ]; 
                //F1558           Q1 = QK-QR
                float Q1 = QK - QR;
                //F1559           TEQUIL = TEQU(K)
                float TEQUIL = TANDSL.TEQU[ K ];
                //F1560           TEQUIL0= TEQU(KREF)
                float TEQUIL0 = TANDSL.TEQU[ KREF ];
                //F1561 !
                //F1562           TE1 = TEQUIL-TEQUIL0
                float TE1 = TEQUIL - TEQUIL0;
                //F1563           TT1 = TGAV(K)-TGAV(KREF)
                float TT1 = TANDSL.TGAV[ K ] - TANDSL.TGAV[ KREF ];
                //F1564           TN1 = TNHAV(K)-TNHAV(KREF)
                float TN1 = TANDSL.TNHAV[ K ] - TANDSL.TNHAV[ KREF ];
                //F1565           TS1 = TSHAV(K)-TSHAV(KREF)
                float TS1 = TANDSL.TSHAV[ K ] - TANDSL.TSHAV[ KREF ];
                //F1566 !
                //F1567 !  CALCULATE 1981-2000 MEAN TEMPERATURE AS REFERENCE LEVEL FOR
                //F1568 !   CALCULATION OF INPUT INTO SCENGEN DRIVER FILES.
                //F1569 !  NOTE : TREF DEPENDS ON CLIMATE MODEL PARAMS (I.E., ON NCLIM)
                //F1570 !
                //F1571           IF(K.GE.217.AND.K.LE.236)TREFSUM=TREFSUM+TGAV(K)
                if( K >= 217 && K <= 236 ) TREFSUM += TANDSL.TGAV[ K ];
                //F1572           IF(K.EQ.236)TREF(NCLIM)=TREFSUM/20.
                if( K == 236 ) TREF[ NSIM.NCLIM ] = TREFSUM / 20.0;
                //F1573 !
                //F1574 !  ******************************************************
                //F1575 !
                //F1576 !  PRINTOUT OPTIONS
                //F1577 !
                //F1578           IF(NOUT.EQ.1)THEN
                if( NOUT == 1 ) {
                    //F1579             COL9   = TN1
                    COL9 = TN1;
                    //F1580             COL10  = TS1
                    COL10 = TS1;
                    //F1581           ENDIF
                }
                //F1582 !
                //F1583           IF(NOUT.EQ.2.OR.NOUT.EQ.5)THEN
                if( NOUT == 2 || NOUT == 5 ) {
                    //F1584             COL9   = TLAND(K)-TLAND(KREF)
                    COL9 = TANDSL.TLAND[ K ] - TANDSL.TLAND[ KREF ];
                    //F1585             COL10  = TOCEAN(K)-TOCEAN(KREF)
                    COL10 = TANDSL.TOCEAN[ K ] - TANDSL.TOCEAN[ KREF ];
                    //F1586             IF(TOCEAN(K).NE.0.0)THEN
                    if( TANDSL.TOCEAN[ K ] != 0.0 ) {
                        //F1587               COL11= TLAND(K)/TOCEAN(K)
                        COL11 = TANDSL.TLAND[ K ] / TANDSL.TOCEAN[ K ];
                        //F1588             ELSE
                    } else {
                        //F1589               COL11= 9.999
                        COL11 = 9.999;
                        //F1590             ENDIF
                    }
                    //F1591           ENDIF
                }
                //F1592 !
                //F1593           IF(NOUT.EQ.3.OR.NOUT.EQ.4)THEN
                if( NOUT == 3 || NOUT == 4 ) {
                    //F1594             COL9   = TE1-TT1
                    COL9 = TE1 - TT1;
                    //F1595             COL10  = TDEEP(K)-TDEEP(KREF)
                    COL10 = TANDSL.TDEEP[ K ] - TANDSL.TDEEP[ KREF ];
                    //F1596           ENDIF
                }
                //F1597 !
                //F1598           EX1=EX(K) -EX(KREF)
                float EX1 = TANDSL.EX[ K ] - TANDSL.EX[ KREF ];
                //F1599           SI1=SLI(K)-SLI(KREF)
                float SI1 = TANDSL.SLI[ K ] - TANDSL.SLI[ KREF ];
                //F1600           SG1=SLG(K)-SLG(KREF)
                float SG1 = TANDSL.SLG[ K ] - TANDSL.SLG[ KREF ];
                //F1601           SA1=SLA(K)-SLA(KREF)
                float SA1 = TANDSL.SLA[ K ] - TANDSL.SLA[ KREF ];
                //F1602           ST1=SLT(K)-SLT(KREF)
                float ST1 = TANDSL.SLT[ K ] - TANDSL.SLT[ KREF ];
                //F1603           SO1=SLO(K)-SLO(KREF)
                float SO1 = TANDSL.SLO[ K ] - TANDSL.SLO[ KREF ];
                //F1604           ST1=EX1+SI1+SG1+SA1+SO1 
                ST1 = EX1 + SI1 + SG1 + SA1 + SO1 ;
                //F1605 !
                //F1606 !  PUT TEMPERATURE AND SEA LEVEL RESULTS FOR FULL GLOBAL FORCING
                //F1607 !   INTO DISPLAY OUTPUT FILES
                //F1608 !
                //F1609           IF(ISCENGEN.NE.9)THEN
                if( NSIM.ISCENGEN != 9 ) {
                    //F1610             IF(NSIM.EQ.1)THEN
                    if( NSIM.NSIM == 1 ) {
                        //F1611               TEMLO(K)  = TT1
                        TEMLO[ K ] = TT1;
                        
                        //F1612               SLLO(K)   = ST1
                        SLLO[ K ] = ST1;
                        //F1613             ENDIF
                    }
                    //F1614           ENDIF
                }
                //F1615 !
                //F1616           IF(NSIM.EQ.2)THEN
                if( NSIM.NSIM == 2 ) {
                    //F1617             TEMMID(K) = TT1
                    TEMMID[ K ] = TT1;
                    //F1618             SLMID(K)  = ST1
                    SLMID[ K ] = ST1;
                    //F1619           ENDIF
                }
                //F1620 !
                //F1621           IF(NSIM.EQ.3)THEN
                if( NSIM.NSIM == 3 ) {
                    //F1622             TEMHI(K)  = TT1
                    TEMHI[ K ] = TT1;
                    //F1623             SLHI(K)   = ST1
                    SLHI[ K ] = ST1;
                    //F1624           ENDIF
                }
                //F1625 !
                //F1626           IF((ISCENGEN.EQ.9).OR.(NSIM.EQ.4))THEN
                if( NSIM.ISCENGEN == 9 || NSIM.NSIM == 4 ) {
                    //F1627             TEMUSER(K)= TT1
                    STOREDVALS.TEMUSER[ K ] = TT1;
                    //F1628             SLUSER(K) = ST1
                    SLUSER[ K ] = ST1;
                    //F1629           ENDIF
                }
                //F1630 !
                //F1631 !  RESULTS FOR ESO2 CONST AFTER 1990 STORED ONLY FOR MID CLIMATE CASE.
                //F1632 !   ZERO VALUES STORED IF ISCENGEN=0 OR =9
                //F1633 !
                //F1634           IF(ISCENGEN.EQ.0.OR.ISCENGEN.EQ.9)THEN
                if( NSIM.ISCENGEN == 0 || NSIM.ISCENGEN == 9 )
                    //F1635             TEMNOSO2(K)= 0.0
                    TEMNOSO2[ K ] = 0.0;
                //F1636           ENDIF
                //F1637 !
                //F1638           IF(NSIM.EQ.6)THEN
                if( NSIM.NSIM == 6 )
                    //F1639             TEMNOSO2(K)= TT1
                    TEMNOSO2[ K ] = TT1;
                //F1640           ENDIF
                //F1641 !
                //F1642 !  PRINT OUT FLAG IS KKKK=1
                //F1643 !
                //F1644           KKKK=0
                int KKK = 0;
                //F1645 !
                //F1646 !  ALWAYS PRINT OUT 1765, IY0 AND 1990 VALUES
                //F1647 !
                //F1648           IF(KYEAR.EQ.1764.OR.KYEAR.EQ.IY0.OR.KYEAR.EQ.1990)KKKK=1
                if( KYEAR == 1764 || KYEAR == IY0 || KYEAR == 1990 ) KKK = 1;
                //F1649 !
                //F1650           IF(KYEAR.GE.IY0)THEN
                if( KYEAR >= IY0 ) {
                    //F1651             PRIN=(KYEAR-KYEAR0+0.01)/XPRT
                    float PRIN = ( KYEAR - KYEAR0 + 0.01 ) / XPRT;
                    //F1652             BIT=PRIN-INT(PRIN)
                    float BIT = PRIN - int( PRIN );
                    //F1653             IF(PRIN.GT.0.0.AND.BIT.LT.0.02)KKKK=1
                    if( PRIN > 0.0 && BIT < 0.02 ) KKK = 1;
                    //F1654             IF(KKKK.EQ.1)THEN
                    if( KKK == 1 ) {
                        //F1655 !
                        //F1656 !  ADD CONSTANT TO ALL TEMPS FOR IPCC DETEX TIME FIGURE
                        //F1657 !
                        //F1658               TT1=TT1+TOFFSET
                        TT1 += TOFFSET;
                        //F1659 !
                        //F1660               IF(NOUT.EQ.1)THEN
                        if( NOUT == 1 ) {
                            //F1661                 SLICE1=SG1+SA1
                            float SLICE1 = SG1 + SA1;
                            //F1662                 SLRAW1=ST1-SO1
                            float SLRAW1 = ST1 - SO1;
                            //F1663                 WRITE (8,191) KYEAR,Q1,TE1,TT1,EX1,SI1,SG1,SA1,SO1,ST1, &
                            //F1664                  COL9,COL10,WNH(K),WSH(K),KYEAR,SLICE1,SLRAW1
                            outfile8 << KYEAR << "," << Q1 << "," << TE1 << "," << TT1 << "," 
                            << EX1 << "," << SI1 << "," << SG1 << "," << SA1 << "," << SO1 << ","
                            << ST1 << "," << COL9 << "," << COL10 << "," << VARW.WNH[ K ] << ","
                            << VARW.WSH[ K ] << ","<< KYEAR << "," << SLICE1 << "," << SLRAW1 << endl;
                            //F1665               ENDIF
                        }
                        //F1666 !
                        //F1667               IF(NOUT.EQ.2)THEN
                        if( NOUT == 2 ) {
                            //F1668                 WRITE (8,192) KYEAR,Q1,TE1,TT1,EX1,SI1,SG1,SA1,SO1,ST1, &
                            //F1669                  COL9,COL10,COL11,WNH(K),WSH(K),KYEAR
                            outfile8 << KYEAR << "," << Q1 << "," << TE1 << "," << TT1 << "," 
                            << EX1 << "," << SI1 << "," << SG1 << "," << SA1 << "," << SO1 << ","
                            << ST1 << "," << COL9 << "," << COL10 << "," << COL11 << "," << VARW.WNH[ K ] << ","
                            << VARW.WSH[ K ] << KYEAR << endl;                            
                            //F1670               ENDIF
                        }
                        //F1671 !
                        //F1672               IF(NOUT.EQ.3)THEN
                        if( NOUT == 3 ) {
                            //F1673                 WRITE (8,193) KYEAR,Q1,TE1,TT1,EX1,SI1,SG1,SA1,SO1,ST1, &
                            //F1674                  COL9,COL10,WNH(K),WSH(K),KYEAR
                            outfile8 << KYEAR << "," << Q1 << "," << TE1 << "," << TT1 << "," 
                            << EX1 << "," << SI1 << "," << SG1 << "," << SA1 << "," << SO1 << ","
                            << ST1 << "," << COL9 << "," << COL10 << "," << VARW.WNH[ K ] << ","
                            << VARW.WSH[ K ] << KYEAR << endl;                            
                            //F1675               ENDIF
                        }
                        //F1676 !
                        //F1677               IF(NOUT.EQ.4)THEN
                        if( NOUT == 4 ) {
                            //F1678                 EQUIVCO2=COBS(1)*EXP(QK/QXX)
                            EQUIVCO2 = COBS.COBS[ 1 ] * exp( float( QK / CLIM.QXX ) );
                            //F1679                 WRITE (8,194) KYEAR,EQUIVCO2,TE1,TT1,EX1,SI1,SG1,SA1, &
                            //F1680                  SO1,ST1,COL9,COL10,WNH(K),WSH(K),KYEAR
                            outfile8 << KYEAR << "," << EQUIVCO2 << "," << TE1 << "," << TT1 << "," 
                            << EX1 << "," << SI1 << "," << SG1 << "," << SA1 << "," << SO1 << ","
                            << ST1 << "," << COL9 << "," << COL10 << "," << VARW.WNH[ K ] << ","
                            << VARW.WSH[ K ] << KYEAR << endl;                            
                            //F1681               ENDIF
                        }
                        //F1682 !
                        //F1683               IF(NOUT.EQ.5)THEN
                        if( NOUT == 5 ) {
                            //F1684                 WRITE(8,195) KYEAR,Q1,TT1,COL11,ST1,EX1,SI1,SG1,SA1,SO1, &
                            //F1685                  WNH(K),KYEAR
                            outfile8 << KYEAR << "," << Q1 << "," << TT1 << "," << COL11 << ","
                            << ST1 << "," << EX1 << "," << SI1 << "," << SG1 << "," << SA1 << "," << SO1 << ","
                            << VARW.WNH[ K ] << "," << KYEAR << endl;                            
                            //F1686               ENDIF
                        }
                        //F1687 !
                        //F1688             ENDIF
                    }
                    //F1689           ENDIF
                }
                //F1690  987  CONTINUE
            } // for
            //F1691 !
            //F1692           IF(NOUT.EQ.1)WRITE(8,171)
            if( NOUT == 1 ) outfile8 << " YEAR,DELTAQ, TEQU, TEMP, EXPN, GLAC, GREENL,ANTAR,Z-XTRA,MSLTOT, TNH, TSH, WNH, WSH,YEAR, GR+ANT,ZTOT-ZXTRA" << endl;
            //F1693           IF(NOUT.EQ.2)WRITE(8,172)
            if( NOUT == 2 ) outfile8 << " YEAR,DELTAQ, TEQU, TEMP, EXPN, GLAC, GREENL,ANTAR, Z-XTRA, MSLTOT, TLAND,  TOCN, TL/TO,   WNH,   WSH,  YEAR" << endl;
            //F1694           IF(NOUT.EQ.3)WRITE(8,173)
            if( NOUT == 3 ) outfile8 << " YEAR,DELTAQ, TEQU, TEMP, EXPN, GLAC, GREENL,  ANTAR, Z-XTRA, MSLTOT, TEQ-T, TDEEP,   WNH,   WSH,  YEAR," << endl;
            //F1695           IF(NOUT.EQ.4)WRITE(8,174)
            if( NOUT == 4 ) outfile8 << " YEAR,EQVCO2, TEQU, TEMP, EXPN, GLAC, GREENL,ANTAR.Z-XTRA.MSLTOT.TEQ-T.TDEEP, WNH, WSH,YEAR," << endl;
            //F1696           IF(NOUT.EQ.5)WRITE(8,175)
            if( NOUT == 5 ) outfile8 << " YEAR,DELTAQ,   TEMP,TL/TO,MSLTOT, EXPN, GLAC, GREENL,ANTAR Z-XTRA, WNH,YEAR" << endl;
            //F1697           WRITE(8,30)
            outfile8 << endl;
            //F1698       ENDIF
        }
        //F1699 !
        //F1700 !  **************************************************************
        //F1701 !
        //F1702 !  DEFINE TEMPERATURE ARRAYS FOR WRITING TO SCENGEN DRIVER FILES.
        //F1703 !   ARRAY SUBSCRIPT NCLIM=1,2,3,4 CORRESPONDS TO LO, MID, HIGH
        //F1704 !   AND USER CLIMATE MODEL PARAMETER SETS.
        //F1705 !  NOTE THAT THESE ARRAYS START WITH KSG=1 IN 1990.
        //F1706 !
        //F1707 !  TSO21,2,3 ARE THE RAW TEMPERATURES IN RESPONSE TO REGIONAL
        //F1708 !   FORCING. THERE ARE AT LEAST TWO WAYS TO CALCULATE THESE.
        //F1709 !  ORIGINALLY (METHOD-1 = SAR VERSION OF SCENGEN) I COMPARED
        //F1710 !   GHG-ALONE RESULTS WITH (GHG) + (REGi EMISSIONS) RESULTS; THE
        //F1711 !   DIFFERENCE GIVING THE RESPONSE TO REGi EMISSIONS.
        //F1712 !  METHOD 2 IS TO COMPARE THE RESULTS FOR 'ALL' EMISSIONS WITH 
        //F1713 !   THOSE FOR .....
        //F1714 !   (ALL) MINUS (REGi EMISSIONS), WHICH EQUALS
        //F1715 !   (GHG)  + (REGj EMISSIONS) + (REGk EMISSIONS)
        //F1716 !   WHERE j AMD k DIFFER FROM i.
        //F1717 !  IN BOTH CASES THERE ARE INTERNAL INCONSISTENCIES BECAUSE OF THE
        //F1718 !   NONLINEAR RELATIONSHIP BETWEEN ESO2 AND INDIRECT AEROSOL FORCING.
        //F1719 !   IN OTHER WORDS, IN GENERAL, TALL MINUS TGHG WILL NOT EQUAL
        //F1720 !   TSO21+TSO22+TSO23.
        //F1721 !  TO CORRECT FOR THIS, I SCALE TSO2i (GIVING XSO2i) BY
        //F1722 !   (DIFF)/(SUM TSO2i) WHERE DIFF = TALL-TGHG.
        //F1723 !  THIS CORRECTION CAN BE A LITTLE ODD AT TIMES. FOR INSTANCE,
        //F1724 !   SUM TSO2i MAY CHANGE SIGN AT A DIFFERENT TIME FROM DIFF, LEADING
        //F1725 !   TO 'WILD' FLUCTUATIONS IN THE RATIO.
        //F1726 !
        //F1727           DO K=197,KEND
        float TALLREF, TGHGREF, TSO21REF, TSO22REF, TSO23REF;
        for( int K=197; K<=Limits.KEND; K++ ) {
            //F1728         KSG=K-196
            int KSG = K - 196;
            //F1729         KKYR=K+1764
            int KKYR = K + 1764;
            //F1730 !
            //F1731         IF(NESO2.EQ.1)THEN
            if( NESO2 == 1 )
                //F1732           TALL(NCLIM,KSG)=TGAV(K)
                TALL[ NSIM.NCLIM ][ KSG ] = TANDSL.TGAV[ K ];
            
            //F1733         ENDIF
            //F1734 !
            //F1735         IF(NESO2.EQ.2)THEN
            if( NESO2 == 2 )
                //F1736           TGHG(NCLIM,KSG)=TGAV(K)
                TGHG[ NSIM.NCLIM ][ KSG ] = TANDSL.TGAV[ K ];
            //F1737         ENDIF
            //F1738 !
            //F1739         IF(NESO2.EQ.3)TSO21(NCLIM,KSG)=TGAV(K)-TGHG(NCLIM,KSG)
            if( NESO2 == 3 ) TSO21[ NSIM.NCLIM ][ KSG ] = TANDSL.TGAV[ K ] - TGHG[ NSIM.NCLIM ][ KSG ];
            //F1740         IF(NESO2.EQ.4)TSO22(NCLIM,KSG)=TGAV(K)-TGHG(NCLIM,KSG)
            if( NESO2 == 4 ) TSO22[ NSIM.NCLIM ][ KSG ] = TANDSL.TGAV[ K ] - TGHG[ NSIM.NCLIM ][ KSG ];
            //F1741         IF(NESO2.EQ.5)TSO23(NCLIM,KSG)=TGAV(K)-TGHG(NCLIM,KSG)
            if( NESO2 == 5 ) TSO23[ NSIM.NCLIM ][ KSG ] = TANDSL.TGAV[ K ] - TGHG[ NSIM.NCLIM ][ KSG ];
            //F1742 !
            //F1743         IF(KKYR.EQ.1990)THEN
            if( KKYR == 1990 ) {
                //F1744           TALLREF =TALL(NCLIM,KSG)
                TALLREF = TALL[ NSIM.NCLIM ][ KSG ];
                //F1745           TGHGREF =TGHG(NCLIM,KSG)
                TGHGREF = TGHG[ NSIM.NCLIM ][ KSG ];
                //F1746           TSO21REF=TSO21(NCLIM,KSG)
                TSO21REF = TSO21[ NSIM.NCLIM ][ KSG ];
                //F1747           TSO22REF=TSO22(NCLIM,KSG)
                TSO22REF = TSO22[ NSIM.NCLIM ][ KSG ];
                //F1748           TSO23REF=TSO23(NCLIM,KSG)
                TSO23REF = TSO23[ NSIM.NCLIM ][ KSG ];
                //F1749         ENDIF
            }
            //F1750       END DO
        }
        //F1751 !
        //F1752       DO K=197,KEND
        for( int K=197; K<=Limits.KEND; K++ ) {
            //F1753         KSG=K-196
            int KSG = K - 196;
            //F1754         TALL(NCLIM,KSG)=TALL(NCLIM,KSG)-TALLREF
            TALL[ NSIM.NCLIM ][ KSG ] -= TALLREF;
            //F1755         TGHG(NCLIM,KSG)=TGHG(NCLIM,KSG)-TGHGREF
            TGHG[ NSIM.NCLIM ][ KSG ] -= TGHGREF;
            //F1756         TSO21(NCLIM,KSG)=TSO21(NCLIM,KSG)-TSO21REF
            TSO21[ NSIM.NCLIM ][ KSG ] -= TSO21REF;
            //F1757         TSO22(NCLIM,KSG)=TSO22(NCLIM,KSG)-TSO22REF
            TSO22[ NSIM.NCLIM ][ KSG ] -= TSO22REF;
            //F1758         TSO23(NCLIM,KSG)=TSO23(NCLIM,KSG)-TSO23REF
            TSO23[ NSIM.NCLIM ][ KSG ] -= TSO23REF;
            //F1759       END DO
        }
        //F1760 !
        //F1761 !  **************************************************************
        //F1762 !
        
        //MOVED from F1954-F1979, because of C scoping rules and b/c this was poorly placed
        //F1954 !  FIRST CALCULATE REFERENCE VALUES (MID IYRQALL)
        //F1955 !
        //F1956          M00=IYRQALL-1764
        const int M00 = IYRQALL - 1764;
        //F1957          M01=M00-1
        const int M01 = M00 - 1;
        //F1958 !
        //F1959          QQQCO2R  = (qco2(M00)     +qco2(M01))     /2.
        const float QQQCO2R = ( FORCE.QCO2[ M00 ] + FORCE.QCO2[ M01 ] ) / 2.0;
        //F1960          QQQMR    = (QM(M00)       +QM(M01))       /2.
        const float QQQMR = ( FORCE.QM[ M00 ] + FORCE.QM[ M01 ] ) / 2.0;
        //F1961          QQQNR    = (QN(M00)       +QN(M01))       /2.
        const float QQQNR = ( FORCE.QN[ M00 ] + FORCE.QN[ M01 ] ) / 2.0;
        //F1962          QQQCFCR  = (QCFC(M00)     +QCFC(M01))     /2.
        const float QQQCFCR = ( FORCE.QCFC[ M00 ] + FORCE.QCFC[ M01 ] ) / 2.0;
        //F1963          QQQSO2R  = (QSO2SAVE(M00) +QSO2SAVE(M01)) /2.
        const float QQQSO2R = ( STOREDVALS.QSO2SAVE[ M00 ] + STOREDVALS.QSO2SAVE[ M01 ] ) / 2.0;
        //F1964          QQQDIRR  = (QDIRSAVE(M00) +QDIRSAVE(M01)) /2.
        const float QQQDIRR = ( STOREDVALS.QDIRSAVE[ M00 ] + STOREDVALS.QDIRSAVE[ M01 ] ) / 2.0;
        //F1965          QQQFOCR  = (QFOC(M00)     +QFOC(M01))     /2.
        const float QQQFOCR = ( JSTART.QFOC[ M00 ] + JSTART.QFOC[ M01 ] ) / 2.0;
        //F1966          QQQMNR   = (QMN(M00)      +QMN(M01))      /2.
        const float QQQMNR = ( TANDSL.QMN[ M00 ] + TANDSL.QMN[ M01 ] ) / 2.0;
        //F1967 !
        //F1968 ! NOTE SPECIAL CASE FOR QOZ BECAUSE OF NONLINEAR CHANGE OVER 1990
        //F1969 !
        //F1970          QQQOZR=(QOZ(M00)+QOZ(M01))/2.
        float QQQOZR = ( TANDSL.QOZ[ M00 ] + TANDSL.QOZ[ M01 ] ) / 2.0;
        //F1971          IF(IYRQALL.EQ.1990)QQQOZR=QOZ(M00)
        if( IYRQALL == 1990 ) QQQOZR = TANDSL.QOZ[ M00 ];
        //F1972 !
        //F1973          QQQLANDR = (QLAND(M00)    +QLAND(M01))     /2.
        const float QQQLANDR = ( TANDSL.QLAND[ M00 ] + TANDSL.QLAND[ M01 ] ) / 2.0;
        //F1974          QQQBIOR  = (QBIO(M00)     +QBIO(M01))     /2.
        const float QQQBIOR = ( TANDSL.QBIO[ M00 ] + TANDSL.QBIO[ M01 ] ) / 2.0;
        //F1975          QQQMO3R  = (QCH4O3(M00)   +QCH4O3(M01))   /2.
        const float QQQMO3R = ( FORCE.QCH4O3[ M00 ] + FORCE.QCH4O3[ M01 ] ) / 2.0;
        //F1976          QQRSTROZ = (QSTRATOZ(M00) +QSTRATOZ(M01)) /2.
        const float QQRSTROZ = ( FORCE.QSTRATOZ[ M00 ] + FORCE.QSTRATOZ[ M01 ] ) / 2.0;
        //F1977          QQRKYMAG = (QKYMAG(M00)   +QKYMAG(M01))   /2.
        const float QQRKYMAG = ( JSTART.QKYMAG[ M00 ] + JSTART.QKYMAG[ M01 ] ) / 2.0;
        //F1978          QQRMONT  = (QMONT(M00)    +QMONT(M01))    /2.
        const float QQRMONT = ( FORCE.QMONT[ M00 ] + FORCE.QMONT[ M01 ] ) / 2.0;
        //F1979          QQROTHER = (QOTHER(M00)   +QOTHER(M01))   /2.
        const float QQROTHER = ( FORCE.QOTHER[ M00 ] + FORCE.QOTHER[ M01 ] ) / 2.0;
        
        // Store 1990 BC and OC values for later use.
        const float QQQBCR = ( FORCE.QBC[ M00 ] + FORCE.QBC[ M01 ] ) / 2.0;
        const float QQQOCR = ( FORCE.QOC[ M00 ] + FORCE.QOC[ M01 ] ) / 2.0;
        
        //F1763 !  PRINT OUT EMISSIONS, CONCS AND FORCING DETAILS
        //F1764 !
        //F1765 !  PRINT OUT INPUT EMISSIONS
        //F1766 !
        //F1767         IF(NESO2.EQ.1)THEN
        if( NESO2 == 1 ) {
            //F1768           WRITE (8,30)
            //F1769           WRITE (8,31)
            //F1770           WRITE (8,30)
            outfile8 << endl << "****************************************************" << endl;
            //F1771           WRITE (8,23)
            outfile8 << endl << "** INPUT EMISSIONS **" << endl;
            //F1772           WRITE (8,231)
            outfile8 << "BALANCED EMISSIONS FOR CH4 & N2O : SO2 EMISSIONS RELATIVE TO 1990" << endl;
            //F1773           WRITE (8,21)
            outfile8 << "YEAR,EFOSS,NETDEF,CH4,N2O,NOX,VOC,CO,SO2REG1,SO2REG2,SO2REG3,CF4,C2F6,HFC125,HFC134A,HFC143A,HFC227ea,HFC245ca,SF6,ESO2TOT,YEAR" << endl;
            //F1774 !
            //F1775 !  PRINTOUT INTERVAL IS DET BY VALUE OF IEMPRT
            //F1776 !
            //F1777           DO K=226,KEND,IEMPRT
            for( int K=226; K<=Limits.KEND; K += IEMPRT ) {
                //F1778             IYEAR=1764+K
                int IYEAR = 1764 + K;
                //F1779             ES1=ESO21(K)-ES1990
                float ES1 = CONCS.ESO22.getval( K ) - Sulph.ES1990;
                //F1780             ES2=ESO22(K)-ES1990
                int ES2 = CONCS.ESO22.getval( K ) - Sulph.ES1990;
                //F1781             ES3=ESO23(K)-ES1990
                int ES3 = CONCS.ESO23.getval( K ) - Sulph.ES1990;
                //F1782             EST=ES1+ES2+ES3
                float EST = ES1 + ES2 + ES3;
                //F1783             WRITE (8,222) IYEAR,EF(K),EDNET(K),ECH4(K),EN2O(K),ENOX(K), &
                //F1784              EVOC(K),ECO(K),ES1,ES2,ES3,ECF4(K),EC2F6(K),E125(K), &
                //F1785              E134A(K),E143A(K),E227(K),E245(K),ESF6(K),EST,EBC(K),EOC(K),IYEAR
                outfile8 << IYEAR << "," << CARB.EF.getval( K ) << "," << METH1.ednet.getval( K ) << "," << CONCS.ECH4.getval( K ) << "," 
                << CONCS.EN2O.getval( K ) << "," << CONCS.ENOX.getval( K ) << "," << CONCS.EVOC.getval( K ) << "," << CONCS.ECO.getval( K ) << "," 
                << ES1 << "," << ES2 << "," << ES3 << "," << NEWCONCS.ECF4.getval( K ) << "," << NEWCONCS.EC2F6.getval( K ) << "," 
                << NEWCONCS.E125.getval( K ) << "," << NEWCONCS.E134A.getval( K ) << "," << NEWCONCS.E143A.getval( K ) << "," 
                << NEWCONCS.E227.getval( K ) << "," << NEWCONCS.E245.getval( K ) << "," << NEWCONCS.ESF6.getval( K ) << "," << EST << "," 
                << CONCS.EBC.getval( K ) << "," << CONCS.EOC.getval( K ) << "," << IYEAR << endl;
                //F1786           END DO
            }
            //F1787 !
            //F1788           WRITE (8,21)
            outfile8 << "YEAR,EFOSS,NETDEF,CH4,N2O,NOX,VOC,CO,SO2REG1,SO2REG2,SO2REG3,CF4,C2F6,HFC125,HFC134A,HFC143A,HFC227ea,HFC245ca,SF6,ESO2TOT,YEAR" << endl;
            //F1789           WRITE (8,30)
            //F1790           WRITE (8,31)
            //F1791           WRITE (8,30)
            outfile8 << endl << "****************************************************" << endl;
            //F1792 !
            //F1793 !  **************************************************************
            //F1794 !
            //F1795 !  PRINT OUT USER CARBON CYCLE DETAILS
            //F1796 !
            //F1797           WRITE(8,24)
            outfile8 << endl << "** CARBON CYCLE DETAILS **" << endl;
            //F1798           WRITE(8,241)LEVCO2
            outfile8 << "CONCENTRATIONS ARE UNCORRECTED MODEL OUTPUT : LEVCO2 =" << CO2READ.LEVCO2 << endl;
            //F1799           WRITE(8,800)R(1)
            outfile8 << "LOW CONC CASE  : NETDEF(80s) = 1.80GtC/yr : GIFFORD FERTILIZATION FACTOR =" << CAR.R[ 1 ] << endl;
            //F1800           WRITE(8,801)R(2)
            outfile8 << "MID CONC CASE  : NETDEF(80s) = 1.10GtC/yr : GIFFORD FERTILIZATION FACTOR =" << CAR.R[ 2 ] << endl;
            //F1801           WRITE(8,802)R(3)
            outfile8 << "HIGH CONC CASE  : NETDEF(80s) = 0.40GtC/yr : GIFFORD FERTILIZATION FACTOR =" << CAR.R[ 3 ] << endl;
            //F1802           WRITE(8,803)DUSER,R(4)
            outfile8 << "USER CONC CASE : NETDEF(80s) =" << METH1.DUSER << "GtC/yr : GIFFORD FERTILIZATION FACTOR =" << CAR.R[ 4 ] << endl;
            //F1803           WRITE(8,804)
            outfile8 << endl << "ALL CASES USE 1980s MEAN OCEAN FLUX OF 2.0GtC/yr" << endl;
            //F1804           WRITE(8,805)
            outfile8 << "DETAILED CARBON CYCLE OUTPUT IS FOR LEVCO2 CASE ONLY" << endl;
            //F1805 !
            //F1806           if(iMeth.eq.1) then
            if( METH1.IMETH == 1 ) {
                //F1807             write(8,806)
                outfile8 << "METHANE OXIDATION TERM INCLUDED IN EMISSIONS" << endl;
                //F1808           else
            } else {
                //F1809             write(8,807)
                outfile8 << "METHANE OXIDATION TERM NOT INCLUDED IN EMISSIONS" << endl;
                //F1810           endif
            }
            //F1811 !
            //F1812           WRITE(8,8071)
            outfile8 << "NOTE: CORRECTION TO MATCH OBSERVED IN 2000 NOT APPLIED IN THIS SECTION" << endl << endl;
            //F1813           MID=0
            const int MID = 0;
            //F1814           IF(MID.NE.1)WRITE(8,810)
            if( MID != 1 ) outfile8 << "ENDYEAR" << endl;
            //F1815           IF(MID.EQ.1)WRITE(8,811)
            if( MID == 1 ) outfile8 << "MIDYEAR" << endl;
            //F1816           WRITE(8,812)
            outfile8 << "YEAR, ETOTAL,  EFOSS, CH4OXN,   NETD, GROSSD,  OFLUX, ABFRAC, PLANT C, HLITT,    SOIL,    CONC,  DEL-M,  YEAR" << endl;
            //F1817 !
            //F1818 !  PRINTOUT INTERVAL IS DET BY VALUE OF ICO2PRT. NOTE THAT CARBON
            //F1819 !   CYCLE MODEL RESULTS GIVE RAW (UNCORRECTED) CO2 CONC OUTPUT.
            //F1820 !
            //F1821           DO K=226,KEND,ICO2PRT
            for( int K=226; K<=Limits.KEND; K += ICO2PRT ) {
                //F1822             IYEAR=1764+K
                int IYEAR = 1764 + K;
                //F1823             CONCOUT=CCO2(LEVCO2,K)
                float CONCOUT = CARB.CCO2.getval( CO2READ.LEVCO2, K );
                //F1824             IF(MID.EQ.1)CONCOUT=(CCO2(LEVCO2,K-1)+CCO2(LEVCO2,K))/2.
                if( MID == 1 ) CONCOUT = ( CARB.CCO2.getval( CO2READ.LEVCO2, K-1 ) + CARB.CCO2.getval( CO2READ.LEVCO2, K ) ) / 2.0;
                //F1825 !
                //F1826             IF(IMETH.EQ.0)THEN
                float TOTE;
                if( METH1.IMETH == 0 )
                    //F1827               TOTE=EF(K)+EDNET(K)
                    TOTE = CARB.EF.getval( K ) + METH1.ednet.getval( K );
                //F1828             ELSE
                else
                    //F1829               TOTE=EF(K)+EDNET(K)+EMETH(K)
                    TOTE = CARB.EF.getval( K ) + METH1.ednet.getval( K ) + METH1.emeth.getval( K );
                //F1830             ENDIF
                //F1831 !
                //F1832             IF(TOTE.EQ.0.0)THEN
                if( TOTE == 0.0 ) {
                    //F1833               IF(DELMASS(4,K).EQ.0.0)THEN
                    float ABX;
                    if( CAR.DELMASS.getval( 4, K ) == 0.0 )
                        //F1834                 ABX=1.0
                        ABX = 1.0;
                    //F1835               ELSE
                    else
                        //F1836                 ABX=DELMASS(4,K)/ABS(DELMASS(4,K))
                        ABX = CAR.DELMASS.getval( 4, K ) / fabs( CAR.DELMASS.getval( 4, K ) );
                    //F1837               ENDIF
                    //F1838               ABFRAC(4,K)=ABX*9.999
                    CAR.ABFRAC.setval( ABX * 9.999, 4, K );
                    //F1839             ELSE
                } else {
                    //F1840               ABFRAC(4,K)=DELMASS(4,K)/TOTE
                    CAR.ABFRAC.setval( CAR.DELMASS.getval( 4, K) / TOTE, 4, K );
                    //F1841             ENDIF
                }
                //F1842 !
                //F1843             IF(ABFRAC(4,K).GT.9.999)ABFRAC(4,K)=9.999
                if( CAR.ABFRAC.getval( 4, K ) > 9.999 ) CAR.ABFRAC.setval( 9.999, 4, K );
                //F1844             IF(ABFRAC(4,K).LT.-9.999)ABFRAC(4,K)=-9.999
                if( CAR.ABFRAC.getval( 4, K ) < -9.999 ) CAR.ABFRAC.setval( -9.999, 4, K );
                //F1845 !
                //F1846             ECH4OX=EMETH(K)
                float ECH4OX = METH1.emeth.getval( K );
                //F1847             IF(IMETH.EQ.0)ECH4OX=0.0
                if( METH1.IMETH == 0 ) ECH4OX = 0.0;
                //F1848             WRITE(8,813)IYEAR,TOTE,EF(K),ECH4OX,EDNET(K),EDGROSS(4,K), &
                //F1849              FOC(4,K),ABFRAC(4,K),PL(4,K),HL(4,K),SOIL(4,K),CONCOUT, &
                //F1850              DELMASS(4,K),IYEAR
                outfile8 << IYEAR << "," << TOTE << "," << CARB.EF.getval( K ) << "," << ECH4OX << "," << METH1.ednet.getval( K ) << "," 
                << CARB.EDGROSS.getval( 4, K ) << "," <<   CARB.FOC.getval( 4, K ) << "," << CAR.ABFRAC.getval( 4, K ) << "," << CARB.PL.getval( 4, K ) << "," 
                << CARB.HL.getval( 4, K ) << "," << CARB.SOIL.getval( 4, K ) << "," << CONCOUT << "," << CAR.DELMASS.getval( 4, K ) << "," << IYEAR << endl;
                //F1851 !
                //F1852           END DO
            }
            //F1853           WRITE(8,812)
            outfile8 << "YEAR, ETOTAL,  EFOSS, CH4OXN,   NETD, GROSSD,  OFLUX, ABFRAC, PLANT C, HLITT,    SOIL,    CONC,  DEL-M,  YEAR" << endl;
            //F1854 !
            //F1855 !  **************************************************************
            //F1856 !
            //F1857 !  PRINT OUT CONCENTRATIONS
            //F1858 !
            //F1859           WRITE (8,30)
            //F1860           WRITE (8,31)
            //F1861           WRITE (8,30)
            outfile8 << endl << "****************************************************" << endl << endl;
            //F1862           IF(KSTART.EQ.0)WRITE (8,20)
            if( KSTART == 0 ) outfile8 << "*** CONCENTRATIONS (CO2,PPM : CH4,N2O,PPB) ***" << endl << "*** MIDYEAR VALUES ***" << endl;
            //F1863           IF(KSTART.EQ.1)WRITE (8,202)
            if( KSTART == 1 ) outfile8 << "*** CONCENTRATIONS (CO2,PPM : CH4,N2O,PPB) ***" << endl << "*** START OF YEAR VALUES FOR YR.GE.1990 IN COLS 2,3,4 ***" << endl;
            //F1864           IF(KSTART.EQ.2)WRITE (8,203)
            if( KSTART == 2 ) outfile8 << "*** CONCENTRATIONS (CO2,PPM : CH4,N2O,PPB) ***" << endl << "*** END OF YEAR VALUES FOR YR.GE.1990 IN COLS 2,3,4 ***" << endl;
            //F1865           WRITE (8,201)
            outfile8 << "<- USER MODEL CONCS -><------ CH4 & CO2 MID CONCS & RANGES ------->" << endl;
            //F1866           WRITE (8,210)
            outfile8 << "YEAR      CO2     CH4    N2O   CH4LO  CH4MID   CH4HI   CO2LO  CO2MID   CO2HI  YEAR TAUCH4" << endl;
            //F1867 !
            //F1868 !  PRINTOUT INTERVAL IS DET BY VALUE OF ICONCPRT
            //F1869 !
            //F1870          CO2(0) =CO2(1)
            CARB.CO2[ 0 ] = CARB.CO2[ 1 ];
            //F1871          CH4(0) =CH4(1)-0.4
            CONCS.CH4[ 0 ] = CONCS.CH4[ 1 ] - 0.4;
            //F1872          CN2O(0)=CN2O(1)
            CONCS.CN2O[ 0 ] = CONCS.CN2O[ 1 ];
            //F1873 !
            //F1874          DO K=1,KEND,ICONCPRT
            for( int K=1; K<=Limits.KEND; K += ICONCPRT ) {
                //F1875            IYEAR=1764+K
                int IYEAR = 1764 + K;
                //F1876 !
                //F1877 !  CONVERT END OF YEAR TO MIDYEAR CONCS
                //F1878 !
                //F1879            CO2MID =(CO2(K)+CO2(K-1))/2.
                float CO2MID = ( CARB.CO2[ K ] + CARB.CO2[ K-1 ] ) / 2.0;
                //F1880            CH4MID =(CH4(K)+CH4(K-1))/2.
                float CH4MID = ( CONCS.CH4[ K ] + CONCS.CH4[ K-1 ] ) / 2.0;
                //F1881            CN2OMID=(CN2O(K)+CN2O(K-1))/2.
                float CN2OMID = ( CONCS.CN2O[ K ] + CONCS.CN2O[ K-1 ] ) / 2.0;
                //F1882 !
                //F1883            IF(K.GE.226)THEN
                if( K >= 226 ) {
                    //F1884 !
                    //F1885              CH4LMID=(CH4L(K)+CH4L(K-1))/2.
                    float CH4LMID = ( METH1.ch4l.getval( K ) + METH1.ch4l.getval( K-1 ) ) / 2.0;
                    //F1886              CH4BMID=(CH4B(K)+CH4B(K-1))/2.
                    float CH4BMID = ( METH1.ch4b.getval( K ) + METH1.ch4b.getval( K-1 ) ) / 2.0;
                    //F1887              CH4HMID=(CH4H(K)+CH4H(K-1))/2.
                    float CH4HMID = ( METH1.ch4h.getval( K ) + METH1.ch4h.getval( K-1 ) ) / 2.0;
                    //F1888 !
                    //F1889              CO2LMID=(CCO2(1,K)+CCO2(1,K-1))/2.
                    float CO2LMID = ( CARB.CCO2.getval( 1, K ) + CARB.CCO2.getval( 1, K-1 ) ) / 2.0;
                    //F1890              CO2BMID=(CCO2(2,K)+CCO2(2,K-1))/2.
                    float CO2BMID = ( CARB.CCO2.getval( 2, K ) + CARB.CCO2.getval( 2, K-1 ) ) / 2.0;
                    //F1891              CO2HMID=(CCO2(3,K)+CCO2(3,K-1))/2.
                    float CO2HMID = ( CARB.CCO2.getval( 3, K ) + CARB.CCO2.getval( 3, K-1 ) ) / 2.0;
                    //F1892 !
                    //F1893 !  ADD CORRECTIONS TO LO, MID, HI CO2 TO FIT OBSERVED DATA IN 2000
                    //F1894 !
                    //F1895              IF(ICO2CORR.EQ.1)THEN
                    if( JSTART.ICO2CORR == 1 ) {
                        //F1896                CO2LMID=CO2LMID+CORREN1
                        CO2LMID += CORREN.CORREN1;
                        //F1897                CO2BMID=CO2BMID+CORREN2
                        CO2BMID += CORREN.CORREN2;
                        //F1898                CO2HMID=CO2HMID+CORREN3
                        CO2HMID += CORREN.CORREN3;
                        //F1899              ENDIF
                    }
                    //F1900 !
                    //F1901              IF(K.LE.236)THEN
                    if( K <= 236 ) {
                        //F1902                CO2LMID=CO2MID
                        //F1903                CO2BMID=CO2MID
                        //F1904                CO2HMID=CO2MID
                        //F1905              ENDIF
                        CO2LMID = CO2BMID = CO2HMID = CO2MID;
                    }
                    //F1906 !
                    //F1907 !  DEFINE LOW, MID AND HIGH CH4 VALUES OVER 1991 TO JSTART YEAR
                    //F1908 !
                    //F1909              IF(K.LT.236)THEN
                    if( K < 236 ) {
                        //F1910                CH4LMID=CH4MID
                        //F1911                CH4BMID=CH4MID
                        //F1912                CH4HMID=CH4MID
                        CH4LMID = CH4BMID = CH4HMID = CH4MID;
                        //F1913              ENDIF
                    }
                    //F1914 !
                    //F1915 !  SPECIFY METHANE LIFETIME OUTPUT
                    //F1916 !
                    //F1917              IF(K.LE.236)THEN
                    float TOR;
                    if( K <= 236 ) {
                        //F1918                TOR=TTUSER
                        TOR = TTUSER;
                        //F1919              ELSE
                    } else {
                        //F1920                TOR=TCH4(K)
                        TOR = METH1.TCH4[ K ];
                        //F1921              ENDIF
                    }
                    //F1922              IF(K.EQ.236)TORREF=TOR
                    if( K == 236 ) TORREF = TOR;
                    //F1923 !
                    //F1924 !  IF KSTART SPECIFIED IN MAGRUN.CFG AS 1, OVERWRITE MIDYEAR
                    //F1925 !   CONCENTRATION VALUES WITH START OR END YEAR VALUES (START IS
                    //F1926 !   WHAT TAR USES, AT LEAST FOR NON-CO2 GASES).
                    //F1927 !
                    //F1928              IF(KSTART.EQ.1)THEN
                    if( KSTART == 1 ) {
                        //F1929                CO2MID=CO2(K-1)
                        CO2MID = CARB.CO2[ K-1 ];
                        //F1930                CH4MID=CH4(K-1)
                        CH4MID = CONCS.CH4[ K-1 ];
                        //F1931                CN2OMID=CN2O(K-1)
                        CN2OMID = CONCS.CN2O[ K-1 ];
                        //F1932              ENDIF
                    }
                    //F1933 !
                    //F1934              IF(KSTART.EQ.2)THEN
                    if( KSTART == 2 ) {
                        //F1935                CO2MID=CO2(K)
                        CO2MID = CARB.CO2[ K ];
                        //F1936                CH4MID=CH4(K)
                        CH4MID = CONCS.CH4[ K ];
                        //F1937                CN2OMID=CN2O(K)
                        CN2OMID = CONCS.CN2O[ K ];
                        //F1938              ENDIF
                    }
                    //F1939 !
                    //F1940              WRITE (8,220) IYEAR,CO2MID,CH4MID,CN2OMID, &
                    //F1941               ch4lMID,ch4bMID,ch4hMID, &
                    //F1942               CO2LMID,CO2BMID,CO2HMID,IYEAR,TOR
                    outfile8 << IYEAR << "," << CO2MID << "," << CH4MID << "," << CN2OMID << ","
                    << CH4LMID << "," << CH4BMID << "," << CH4HMID << ","
                    << CO2LMID << "," << CO2BMID << "," << CO2HMID << "," << IYEAR << "," << TOR << endl;
                    //F1943            ELSE
                    //F1944              WRITE (8,221) IYEAR,CO2MID,CH4MID,CN2OMID,IYEAR
                }
                else outfile8 << IYEAR << "," << CO2MID << "," << CH4MID << "," << CN2OMID << ",,,,,," << IYEAR << endl;
                    //F1945            ENDIF
                //F1946          END DO
            }
            //F1947 !
            //F1948          WRITE (8,210)
            outfile8 << "YEAR      CO2     CH4    N2O   CH4LO  CH4MID   CH4HI   CO2LO  CO2MID   CO2HI  YEAR TAUCH4" << endl;
            //F1949          WRITE (8,201)
            outfile8 << "<- USER MODEL CONCS -><------ CH4 & CO2 MID CONCS & RANGES ------->" << endl;
            //F1950 !
            //F1951 !  **************************************************************
            //F1952 !
            //F1953 !  PRINT OUT TABLES OF DELTA-Q FROM IYRQALL AND 1765 TO MAG.OUT.
            //MOVED lines F1954-F1979 moved to after F1762 -- see note there
            //F1980 !
            //F1981 !   PRINT OUT DELTA-Q FROM MID 1990 TO MAG.OUT
            //F1982 !
            //F1983          write(8,30)
            //F1984          write(8,31)
            //F1985          WRITE(8,30)
            outfile8 << endl << "****************************************************" << endl << endl;
            //F1986          write(8,55)IYRQALL
            outfile8 << "** GAS BY GAS DELTA-Q FROM " << IYRQALL << " : MIDYEAR VALUES **" << endl;
            //F1987          write(8,56)
            outfile8 << "(BASED ON MID-YEAR FORCING VALUES : QEXTRA NOT INCLUDED)" << endl;
            //F1988          write(8,561)
            outfile8 << "CH4tot INCLUDES STRATH2O : TROPO3 INCLUDES CH4 COMPONENT" << endl;
            //F1989          write(8,5611)
            outfile8 << "QAERMN IS THE SUM OF NITRATE AND MINERAL DUST AEROSOL FORCING" << endl;
            //F1990          IF(IO3FEED.EQ.1)write(8,562)
            outfile8 << "HALOtot INCLUDES STRAT O3" << endl;
            //F1991          IF(IO3FEED.EQ.0)write(8,563)
            outfile8 << "HALOtot (AND QTOTAL) DOES NOT INCLUDE STRAT O3" << endl;
            //F1992          write(8,57)
            outfile8 << "YEAR,CO2,CH4tot,N2O, HALOtot,TROPOZ,SO4DIR,SO4IND,BIOAER,FOC+FBC,QAERMN,QLAND, TOTAL, YEAR,CH4-O3, STRATO3, MONTDIR,QKYOTO,BC,OC" << endl;
            //F1993 !
            //F1994 !  PRINTOUT INTERVAL IS DET BY VALUE OF IQGASPRT
            //F1995 !
            //F1996          DO K=1990,IYEND,IQGASPRT
            for( int K=1990; K<=IYEND; K += IQGASPRT ) {
                //F1997            IYR = K-1990+226
                int IYR = K - 1990 + 226;
                //F1998            IYRP=IYR-1
                int IYRP = IYR - 1;
                //F1999 ! 
                //F2000            DELQCO2 = (QCO2(IYR)+QCO2(IYRP))/2.-QQQCO2R
                float DELQCO2 = ( FORCE.QCO2[ IYR ] + FORCE.QCO2[ IYRP ] ) / 2.0 - QQQCO2R;
                //F2001            DELQM   = (QM(IYR)+QM(IYRP))/2.    -QQQMR
                float DELQM = ( FORCE.QM[ IYR ] + FORCE.QM[ IYRP ] ) / 2.0 - QQQMR;
                //F2002            DELQN   = (QN(IYR)+QN(IYRP))/2.    -QQQNR
                float DELQN = ( FORCE.QN[ IYR ] + FORCE.QN[ IYRP ] ) / 2.0 - QQQNR;
                //F2003            DELQCFC = (QCFC(IYR)+QCFC(IYRP))/2.-QQQCFCR
                float DELQCFC = ( FORCE.QCFC[ IYR ] + FORCE.QCFC[ IYRP ] ) / 2.0 - QQQCFCR;
                //F2004 !
                //F2005 !  NOTE : DELQSO2 AND DELQDIR BOTH INCLUDE QFOC
                //F2006 !
                //F2007            DELQSO2 = (QSO2SAVE(IYR)+QSO2SAVE(IYRP))/2.-QQQSO2R
                float DELQSO2 = ( STOREDVALS.QSO2SAVE[ IYR ] + STOREDVALS.QSO2SAVE[ IYRP ] ) / 2.0 - QQQSO2R;
                //F2008            DELQDIR = (QDIRSAVE(IYR)+QDIRSAVE(IYRP))/2.-QQQDIRR
                float DELQDIR = ( STOREDVALS.QDIRSAVE[ IYR ] + STOREDVALS.QDIRSAVE[ IYRP ] ) / 2.0 - QQQDIRR;
                //F2009            DELQIND = DELQSO2-DELQDIR
                float DELQIND = DELQSO2 - DELQDIR;
                //F2010            DELQFOC = (QFOC(IYR)+QFOC(IYRP))/2.-QQQFOCR
                float DELQFOC = ( JSTART.QFOC[ IYR ] + JSTART.QFOC[ IYRP ] ) / 2.0 - QQQFOCR;
                //F2011            DELQD   = DELQDIR-DELQFOC
                float DELQD = DELQDIR - DELQFOC;
                //F2012            DELQMN  = (QMN(IYR)+QMN(IYRP))/2.-QQQMNR
                float DELQMN = ( TANDSL.QMN[ IYR ] + TANDSL.QMN[ IYRP ] ) / 2.0 - QQQMNR;
                //F2013 !
                //F2014 ! NOTE SPECIAL CASE FOR QOZ BECAUSE OF NONLINEAR CHANGE OVER 1990
                //F2015 !
                //F2016            IF(IYR.EQ.226)THEN
                float QOZMID;
                if( IYR == 226 ) {
                    //F2017              QOZMID= QOZ(IYR)
                    QOZMID = TANDSL.QOZ[ IYR ];
                    //F2018            ELSE
                } else {
                    //F2019              QOZMID= (QOZ(IYR)+QOZ(IYRP))/2.
                    QOZMID = ( TANDSL.QOZ[ IYR ] + TANDSL.QOZ[ IYRP ] ) / 2.0;
                    //F2020            ENDIF
                }
                //F2021            DELQOZ  = QOZMID-QQQOZR
                float DELQOZ = QOZMID - QQQOZR;
                //F2022 !
                //F2023            DELQLAND= (QLAND(IYR)+QLAND(IYRP))/2.-QQQLANDR
                float DELQLAND = ( TANDSL.QLAND[ IYR ] + TANDSL.QLAND[ IYRP ] ) / 2.0 - QQQLANDR;
                //F2024            DELQBIO = (QBIO(IYR)+QBIO(IYRP))/2.-QQQBIOR
                float DELQBIO = ( TANDSL.QBIO[ IYR ] + TANDSL.QBIO[ IYRP ] ) / 2.0 - QQQBIOR;
                //F2025            DELQTOT = DELQCO2+DELQM+DELQN+DELQCFC+DELQSO2+DELQBIO &
                //F2026             +DELQOZ+DELQLAND+DELQMN
                float DELQTOT = DELQCO2 + DELQM + DELQN + DELQCFC + DELQSO2 + DELQBIO + DELQOZ + DELQLAND + DELQMN;
                //F2027 !
                //F2028            DQCH4O3 = (QCH4O3(IYR)+QCH4O3(IYRP))/2.-QQQMO3R
                float DQCH4O3 = ( FORCE.QCH4O3[ IYR ] + FORCE.QCH4O3[ IYRP ] ) / 2.0 - QQQMO3R;
                //F2029            DELQM   = DELQM-DQCH4O3
                DELQM -= DQCH4O3;
                //F2030            DELQOZ  = DELQOZ+DQCH4O3
                DELQOZ += DQCH4O3;
                //F2031 !
                //F2032            DELSTROZ= (QSTRATOZ(IYR)+QSTRATOZ(IYRP))/2.-QQRSTROZ
                float DELSTROZ = ( FORCE.QSTRATOZ[ IYR ] + FORCE.QSTRATOZ[ IYRP ] ) / 2.0 - QQRSTROZ;
                //F2033            IF(IO3FEED.EQ.0)DELSTROZ=0.0
                //F2034 !
                //F2035            DELKYMAG = (QKYMAG(IYR) +QKYMAG(IYRP))  /2.-QQRKYMAG
                float DELKYMAG = ( JSTART.QKYMAG[ IYR ] + JSTART.QKYMAG[ IYRP ] ) / 2.0 - QQRKYMAG;
                //F2036            DELMONT  = (QMONT(IYR)  +QMONT(IYRP))   /2.-QQRMONT
                float DELMONT = ( FORCE.QMONT[ IYR ] + FORCE.QMONT[ IYRP ] ) / 2.0 - QQRMONT;
                //F2037            DELOTHER = (QOTHER(IYR) +QOTHER(IYRP))  /2.-QQROTHER
                float DELOTHER = ( FORCE.QOTHER[ IYR ] + FORCE.QOTHER[ IYRP ] ) / 2.0 - QQROTHER;
                //F2038            DELKYOTO = DELKYMAG+DELOTHER
                float DELKYOTO = DELKYMAG + DELOTHER;
                
                // Adding BC OC to this table.
                float DELQBC = ( FORCE.QBC[ IYR ] + FORCE.QBC[ IYRP ] ) / 2.0 - QQQBCR;
                float DELQOC = ( FORCE.QOC[ IYR ] + FORCE.QOC[ IYRP ] ) / 2.0 - QQQOCR;
                // Add BC and OC to total.
                DELQTOT += DELQBC + DELQOC;
                //F2039 !
                //F2040            WRITE(8,571)K,DELQCO2,DELQM,DELQN,DELQCFC,DELQOZ, &
                //F2041             DELQD,DELQIND,DELQBIO,DELQFOC,DELQMN,DELQLAND,DELQTOT, &
                //F2042             K,DQCH4O3,DELSTROZ,DELMONT,DELKYOTO
                outfile8 << K << "," << DELQCO2 << "," << DELQM << "," << DELQN << "," << DELQCFC << "," << DELQOZ << "," 
                << DELQD << "," << DELQIND << "," << DELQBIO << "," << DELQFOC << "," << DELQMN << "," << DELQLAND << "," 
                << DELQTOT << "," << K << "," << DQCH4O3 << "," << DELSTROZ << "," << DELMONT << "," << DELKYOTO << ","
                << DELQBC << "," << DELQOC << endl;
                //F2043          END DO
            }
            //F2044 !
            //F2045 !  ************************************************************
            //F2046 !
            //F2047 !  SAVE ARRAYS FOR TROPOZ, STRATOZ, CFC12, C11EFF
            //F2048 !
            //F2049         ALF11=0.000250
            float ALF11 = 0.000250;
            //F2050         ALF12=0.000320
            float ALF12 = 0.000320;
            //F2051       DO KK=1,IYEND-1764
            for( int KK=1; KK<=IYEND-1764; KK++ ) {
                //F2052         KKP=1
                int KKP = 1;
                //F2053         IF(KK.GT.1)KKP=KK-1
                if( KK > 1 ) KKP = KK - 1;
                //F2054         QQQCH4O3=(QCH4O3(KK)+QCH4O3(KKP))/2.
                float QQQCH4O3 = ( FORCE.QCH4O3[ KK ] + FORCE.QCH4O3[ KKP ] ) / 2.0;
                //F2055         QTROZ(KK)=QQQCH4O3+(QOZ(KK)+QOZ(KKP))/2.
                QTROZ[ KK ] = QQQCH4O3 + ( TANDSL.QOZ[ KK ] + TANDSL.QOZ[ KKP ] ) / 2.0;
                //F2056         QSTROZ(KK)=(QSTRATOZ(KK)+QSTRATOZ(KKP))/2.
                QSTROZ[ KK ] = ( FORCE.QSTRATOZ[ KK ] + FORCE.QSTRATOZ[ KKP ] ) / 2.0;
                //F2057         QQQCFC=(QCFC(KK)+QCFC(KKP))/2.
                float QQQCFC = ( FORCE.QCFC[ KK ] + FORCE.QCFC[ KKP ] ) / 2.0;
                //F2058         C11EFF(KK)=(QQQCFC-QSTROZ(KK)-CFC12(Kk)*ALF12)/ALF11
                C11EFF[ KK ] = ( QQQCFC - QSTROZ[ KK ] - FORCE.CFC12[ KK ] * ALF12 ) / ALF11;
                //F2059       END DO
            }
            //F2060 !
            //F2061 !  ************************************************************
            //F2062 !
            //F2063 !  NOW PRINT OUT FORCING CHANGES FROM MID 1765.
            //F2064 !
            //F2065          WRITE(8,57)
            outfile8 << "YEAR,CO2,CH4tot,N2O, HALOtot,TROPOZ,SO4DIR,SO4IND,BIOAER,FOC+FBC,QAERMN,QLAND, TOTAL, YEAR,CH4-O3, STRATO3, MONTDIR,QKYOTO,BC,OC" << endl;
            //F2066          write(8,30)
            //F2067          write(8,31)
            //F2068          WRITE(8,30)
            outfile8 << endl << "****************************************************" << endl << endl;
            //F2069          write(8,58)
            outfile8 << "** GAS BY GAS DELTA-Q FROM 1765 : MIDYEAR VALUES **" << endl;
            //F2070          write(8,56)
            outfile8 << "(BASED ON MID-YEAR FORCING VALUES : QEXTRA NOT INCLUDED)" << endl;
            //F2071          write(8,561)
            outfile8 << "CH4tot INCLUDES STRATH2O : TROPO3 INCLUDES CH4 COMPONENT" << endl;
            //F2072          write(8,5611)
            outfile8 << "QAERMN IS THE SUM OF NITRATE AND MINERAL DUST AEROSOL FORCING" << endl;
            //F2073          IF(IO3FEED.EQ.1)write(8,562)
            if( METH1.IO3FEED == 1 ) outfile8 << "HALOtot INCLUDES STRAT O3" << endl;
            //F2074          IF(IO3FEED.EQ.0)write(8,563)
            if( METH1.IO3FEED == 0 ) outfile8 << "HALOtot (AND QTOTAL) DOES NOT INCLUDE STRAT O3" << endl;
            //F2075          write(8,573)
            outfile8 << "YEAR,CO2,CH4tot,N2O, HALOtot,TROPOZ,SO4DIR,SO4IND,BIOAER,FOC+FBC,QAERMN,QLAND, TOTAL, YEAR,CH4-O3, STRATO3, MONTDIR,QKYOTO,BC,OC,QEXTRA" << endl;
            //F2076 !
            //F2077        DO K=1770,IYEND,IQGASPRT
            for( int K=1770; K<=IYEND; K+=IQGASPRT ) {
                //F2078          IYR = K-1990+226
                int IYR = K - 1990 + 226;
                //F2079          IYRP=IYR-1
                int IYRP = IYR-1;
                //F2080 !
                //F2081          QQQSO2 = 0.0
                float QQQSO2 = 0.0;
                //F2082          QQQDIR = 0.0
                float QQQDIR = 0.0;
                //F2083          IF(K.GT.1860)THEN
                if( K > 1860 ) {
                    //F2084            QQQSO2 = (QSO2SAVE(IYR)+QSO2SAVE(IYRP))/2.
                    QQQSO2 = ( STOREDVALS.QSO2SAVE[ IYR ] + STOREDVALS.QSO2SAVE[ IYRP ] ) / 2.0;
                    //F2085            QQQDIR = (QDIRSAVE(IYR)+QDIRSAVE(IYRP))/2.
                    QQQDIR = ( STOREDVALS.QDIRSAVE[ IYR ] + STOREDVALS.QDIRSAVE[ IYRP ] ) / 2.0;
                    //F2086          ENDIF
                }
                //F2087          QQQIND = QQQSO2-QQQDIR
                float QQQIND = QQQSO2 - QQQDIR;
                //F2088 ! 
                //F2089          QQQCO2 = (QCO2(IYR)+QCO2(IYRP))/2.
                float QQQCO2 = ( FORCE.QCO2[ IYR ] + FORCE.QCO2[ IYRP ] ) / 2.0;
                //F2090          QQQM   = (QM(IYR)+QM(IYRP))/2.
                float QQQM = ( FORCE.QM[ IYR ] + FORCE.QM[ IYRP ] ) / 2.0;
                //F2091          QQQN   = (QN(IYR)+QN(IYRP))/2.
                float QQQN = ( FORCE.QN[ IYR ] + FORCE.QN[ IYRP ] ) / 2.0;
                //F2092          QQQCFC = (QCFC(IYR)+QCFC(IYRP))/2.
                float QQQCFC = ( FORCE.QCFC[ IYR ] + FORCE.QCFC[ IYRP ] ) / 2.0;
                //F2093          QQQOZ  = (QOZ(IYR)+QOZ(IYRP))/2.
                float QQQOZ = ( TANDSL.QOZ[ IYR ] + TANDSL.QOZ[ IYRP ] ) / 2.0;
                //F2094          QQQFOC = (QFOC(IYR)+QFOC(IYRP))/2.
                float QQQFOC = ( JSTART.QFOC[ IYR ] + JSTART.QFOC[ IYRP ] ) / 2.0;
                //F2095          QQQMN  = (QMN(IYR)+QMN(IYRP))/2.
                float QQQMN = ( TANDSL.QMN[ IYR ] + TANDSL.QMN[ IYRP ] ) / 2.0;
                //F2096 !
                //F2097 ! NOTE SPECIAL CASE FOR QOZ BECAUSE OF NONLINEAR CHANGE OVER 1990
                //F2098 !
                //F2099          IF(IYR.EQ.226)QQQOZ=QOZ(IYR)
                if( IYR == 226 ) QQQOZ = TANDSL.QOZ[ IYR ];
                //F2100 !
                //F2101          QQQLAND= (QLAND(IYR)+QLAND(IYRP))/2.
                float QQQLAND = ( TANDSL.QLAND[ IYR ] + TANDSL.QLAND[ IYRP ] ) / 2.0;
                //F2102          QQQBIO = (QBIO(IYR)+QBIO(IYRP))/2.
                float QQQBIO = ( TANDSL.QBIO[ IYR ] + TANDSL.QBIO[ IYRP ] ) / 2.0;
                //F2103          QQQTOT = QQQCO2+QQQM+QQQN+QQQCFC+QQQSO2+QQQBIO+QQQOZ+QQQLAND &
                //F2104           +QQQMN
                float QQQTOT = QQQCO2 + QQQM + QQQN + QQQCFC + QQQSO2 + QQQBIO + QQQOZ + QQQLAND + QQQMN;
                //F2105 !
                //F2106          QQCH4O3= (QCH4O3(IYR)+QCH4O3(IYRP))/2.
                float QQCH4O3 = ( FORCE.QCH4O3[ IYR ] + FORCE.QCH4O3[ IYRP ] ) / 2.0;
                //F2107          QQQM   = QQQM-QQCH4O3
                QQQM = QQQM - QQCH4O3;
                //F2108          QQQOZ  = QQQOZ+QQCH4O3
                QQQOZ = QQQOZ + QQCH4O3;
                //F2109          QQQD   = QQQDIR-QQQFOC
                float QQQD = QQQDIR - QQQFOC;
                //F2110 !
                //F2111          QQQSTROZ= (QSTRATOZ(IYR)+QSTRATOZ(IYRP))/2.
                float QQQSTROZ = ( FORCE.QSTRATOZ[ IYR ] + FORCE.QSTRATOZ[ IYRP ] ) / 2.0;
                //F2112          IF(IO3FEED.EQ.0)QQQSTROZ=0.0
                //F2113 !
                //F2114          QQQKYMAG = (QKYMAG(IYR)+QKYMAG(IYRP))/2.
                float QQQKYMAG = ( JSTART.QKYMAG[ IYR ] + JSTART.QKYMAG[ IYRP ] ) / 2.0;
                //F2115          QQQMONT  = (QMONT(IYR) +QMONT(IYRP)) /2.
                float QQQMONT = ( FORCE.QMONT[ IYR ] + FORCE.QMONT[ IYRP ] ) / 2.0;
                //F2116          QQQOTHER = (QOTHER(IYR)+QOTHER(IYRP))/2.
                float QQQOTHER = ( FORCE.QOTHER[ IYR ] + FORCE.QOTHER[ IYRP ] ) / 2.0;
                //F2117          QQQKYOTO = QQQKYMAG+QQQOTHER
                float QQQKYOTO = QQQKYMAG + QQQOTHER;
                //F2118 !
//REVERT                //F2119 ! Add BC, OC, and QExtra forcing to output. Note this is not included in foring total since is in QExtra instead
                //F2120          QQQEXTRA = ( QEXNH(IYR)+QEXSH(IYR)+QEXNHO(IYR)+QEXNHL(IYR) + &
                //F2121                       QEXNH(IYRP)+QEXSH(IYRP)+QEXNHO(IYRP)+QEXNHL(IYRP) )/2.
                float QQQEXTRA = ( QADD.QEXNH[ IYR ] + QADD.QEXSH[ IYR ] + QADD.QEXNHO[ IYR ] + QADD.QEXNHL[ IYR ] 
                                  + QADD.QEXNH[ IYRP ] + QADD.QEXSH[ IYRP ] + QADD.QEXNHO[ IYRP ] + QADD.QEXNHL[ IYRP ] ) / 2.0;
                //F2122          QQQBC = ( QBC(IYR) + QBC(IYRP) )/2.
                float QQQBC = ( FORCE.QBC[ IYR ] + FORCE.QBC[ IYRP ] ) / 2.0;
                //F2123          QQQOC = ( QOC(IYR) + QOC(IYRP) )/2.
                float QQQOC = ( FORCE.QOC[ IYR ] + FORCE.QOC[ IYRP ] ) / 2.0;
                // Add BC and OC to total forcing.
                QQQTOT += QQQOC + QQQBC;
                //F2124 
                //F2125          WRITE(8,571)K,QQQCO2,QQQM,QQQN,QQQCFC,QQQOZ,QQQD,QQQIND, &
                //F2126          QQQBIO,QQQFOC,QQQMN,QQQLAND,QQQTOT,K,QQCH4O3,QQQSTROZ,QQQMONT, &
                //F2127          QQQKYOTO, QQQBC, QQQOC, QQQEXTRA
                outfile8 << K << "," << QQQCO2 << "," << QQQM << "," << QQQN << "," << QQQCFC << "," << QQQOZ << "," << QQQD << "," 
                << QQQIND << "," << QQQBIO << "," << QQQFOC << "," << QQQMN << "," << QQQLAND << "," << QQQTOT << "," << K << "," 
                << QQCH4O3 << "," << QQQSTROZ << "," << QQQMONT << "," << QQQKYOTO << "," << QQQBC << "," << QQQOC << "," << QQQEXTRA << endl;
                //F2128        END DO
            }
            //F2129 !
            //F2130        WRITE(8,573)
            outfile8 << "YEAR,CO2,CH4tot,N2O, HALOtot,TROPOZ,SO4DIR,SO4IND,BIOAER,FOC+FBC,QAERMN,QLAND, TOTAL, YEAR,CH4-O3, STRATO3, MONTDIR,QKYOTO,BC,OC,QEXTRA" << endl;
            //F2131 
            //F2132 
            //F2133 ! *******************************************************************************************
            //F2134 ! sjs -- write out halocarbon forcings separately as well
            //F2135 ! *******************************************************************************************
            //F2136 
            //F2137 ! Calculate 1990 halocarbon forcings so these can be output
            //F2138 	QCF4_ar(226) = CF4(226) * ( QCF4_ar(227)/CF4(227) )
            HALOF.QCF4_ar[ 226 ] = NEWCONCS.CF4[ 226 ] * ( HALOF.QCF4_ar[ 227 ] / NEWCONCS.CF4[ 227 ] );
            //F2139 	QC2F6_ar(226) = C2F6(226) * ( QC2F6_ar(227)/C2F6(227) )
            HALOF.QC2F6_ar[ 226 ] = NEWCONCS.C2F6[ 226 ] * ( HALOF.QC2F6_ar[ 227 ] / NEWCONCS.C2F6[ 227 ] );
            //F2140 	qSF6_ar(226) = CSF6(226) * ( qSF6_ar(227)/CSF6(227) )
            HALOF.qSF6_ar[ 226 ] = NEWCONCS.CSF6[ 226 ] * ( HALOF.qSF6_ar[ 227 ] / NEWCONCS.CSF6[ 227 ] );
            //F2141 
            //F2142 ! Approximate to same 1989 forcing
            //F2143 	QCF4_ar(225) = QCF4_ar(226)
            HALOF.QCF4_ar[ 225 ] = HALOF.QCF4_ar[ 226 ];
            //F2144 	QC2F6_ar(225) = QC2F6_ar(226)
            HALOF.QC2F6_ar[ 225 ] = HALOF.QC2F6_ar[ 226 ];
            //F2145 	qSF6_ar(225) = qSF6_ar(226)
            HALOF.qSF6_ar[ 225 ] = HALOF.qSF6_ar[ 226 ];
            //F2146 
            //F2147 !
            //F2148 !
            //F2149 
            //F2150        write(8,30)
            //F2151         write(8,31)
            //F2152         WRITE(8,30)
            outfile8 << endl << "****************************************************" << endl << endl;
            //F2153         WRITE(8,*) "Halocarbon Emissions"
            outfile8 << "Halocarbon Emissions" << endl;
            //F2154         WRITE(8,588)
            outfile8 << "YEAR,HFC245,HFC134A,HFC125,HFC227,HFC143A,CF4,C2F6,SF6," << endl;
            //F2155 
            //F2156         DO K=1990,IYEND,IQGASPRT
            for( int K=1990; K<=IYEND; K+=IQGASPRT ) {
                //F2157           IYR = K-1990+226
                int IYR = K - 1990 + 226;
                //F2158           IYRP=IYR-1
                //UNUSED int IYRP = IYR - 1;
                //F2159 !
                //F2160         WRITE(8,589)K,E245(IYR),E134A(IYR),E125(IYR),E227(IYR), &
                //F2161        E143A(IYR),ECF4(IYR),EC2F6(IYR),ESF6(IYR)
                outfile8 << K << "," << NEWCONCS.E245.getval( IYR ) << "," << NEWCONCS.E134A.getval( IYR ) << "," 
                << NEWCONCS.E125.getval( IYR ) << "," << NEWCONCS.E227.getval( IYR ) << "," << NEWCONCS.E143A.getval( IYR ) << "," 
                << NEWCONCS.ECF4.getval( IYR ) << "," << NEWCONCS.EC2F6.getval( IYR ) << "," << NEWCONCS.ESF6.getval( IYR ) << endl;
                //F2162          
                //F2163       END DO
            }
            //F2164 
            //F2165        write(8,30)
            //F2166         write(8,31)
            //F2167         WRITE(8,30)
            outfile8 << endl << "****************************************************" << endl << endl;
            //F2168         WRITE(8,*) "Halocarbon Concentrations"
            outfile8 << "Halocarbon Concentrations" << endl;
            //F2169         WRITE(8,588)
            outfile8 << "YEAR,HFC245,HFC134A,HFC125,HFC227,HFC143A,CF4,C2F6,SF6," << endl;
            //F2170 
            //F2171         DO K=1990,IYEND,IQGASPRT
            for( int K=1990; K<=IYEND; K+=IQGASPRT ) {
                //F2172           IYR = K-1990+226
                int IYR = K - 1990 + 226;
                //F2173           IYRP=IYR-1
                //UNUSED int IYRP = IYR - 1;
                //F2174 
                //F2175         WRITE(8,589)K,C245(IYR),C134A(IYR),C125(IYR),C227(IYR), &
                //F2176        C143A(IYR),CF4(IYR),C2F6(IYR),CSF6(IYR)
                outfile8 << K << "," << NEWCONCS.C245[ IYR ] << "," << NEWCONCS.C134A[ IYR ] << "," 
                << NEWCONCS.C125[ IYR ] << "," << NEWCONCS.C227[ IYR ] << "," << NEWCONCS.C143A[ IYR ] << "," 
                << NEWCONCS.CF4[ IYR ] << "," << NEWCONCS.C2F6[ IYR ] << "," << NEWCONCS.CSF6[ IYR ] << endl;
                //F2177           
                //F2178       END DO
            }
            //F2179 
            //F2180        write(8,30)
            //F2181         write(8,31)
            //F2182         WRITE(8,30)
            outfile8 << endl << "****************************************************" << endl << endl;
            //F2183         WRITE(8,*) "Halocarbon Forcing"
            outfile8 << "Halocarbon Forcing" << endl;
            //F2184         write(8,590)
            outfile8 << "YEAR,HFC245,HFC134A,HFC125,HFC227,HFC143A,CF4,C2F6,SF6,Qother,QMont,QStratOz,HaloTot,KyotoTot" << endl;
            //F2185 
            //F2186         DO K=1990,IYEND,IQGASPRT
            for( int K=1990; K<=IYEND; K+=IQGASPRT ) {
                //F2187           IYR = K-1990+226
                int IYR = K - 1990 + 226;
                //F2188           IYRP=IYR-1
                //UNUSED int IYRP = IYR - 1;
                //F2189 
                //F2190         WRITE(8,589)K, &
                //F2191        Q245_ar(IYR),Q134A_ar(IYR),Q125_ar(IYR), &
                //F2192        Q227_ar(IYR),Q143A_ar(IYR),QCF4_ar(IYR), &
                //F2193        QC2F6_ar(IYR),qSF6_ar(IYR), QOTHER(IYR), &
                //F2194        QMONT(IYR),QSTRATOZ(IYR),QMONT(IYR)+ &
                //F2195        QKYMAG(IYR)+QOTHER(IYR)+QSTRATOZ(IYR), &
                //F2196        QKYMAG(IYR)+QOTHER(IYR)
                outfile8 << HALOF.Q245_ar[ IYR ] << "," << HALOF.Q134A_ar[ IYR ] << "," << HALOF.Q125_ar[ IYR ] << "," 
                << HALOF.Q227_ar[ IYR ] << "," << HALOF.Q143A_ar[ IYR ] << "," << HALOF.QCF4_ar[ IYR ] << "," << HALOF.QC2F6_ar[ IYR ] << "," 
                << HALOF.qSF6_ar[ IYR ] << "," << FORCE.QOTHER[ IYR ] << "," << FORCE.QMONT[ IYR ] << "," << FORCE.QSTRATOZ[ IYR ] << "," 
                << FORCE.QMONT[ IYR ] + JSTART.QKYMAG[ IYR ] + FORCE.QOTHER[ IYR ] + FORCE.QSTRATOZ[ IYR ] << "," 
                << JSTART.QKYMAG[ IYR ] + FORCE.QOTHER[ IYR ] << endl;
                //F2197 
                //F2198 
                //F2199       END DO
            }
            //F2200       write(8,590)
            outfile8 << "YEAR,HFC245,HFC134A,HFC125,HFC227,HFC143A,CF4,C2F6,SF6,Qother,QMont,QStratOz,HaloTot,KyotoTot" << endl;
            //F2201         WRITE(8,30)
            //F2202         WRITE(8,30)
            outfile8 << endl << endl;
            //F2203 
            //F2204  588  FORMAT (1X,'YEAR,HFC245,HFC134A,HFC125,HFC227,HFC143A,','CF4,C2F6,SF6,')
            //F2205  590  FORMAT (1X,'YEAR,HFC245,HFC134A,HFC125,HFC227,HFC143A,','CF4,C2F6,SF6,Qother,QMont,QStratOz,HaloTot,','KyotoTot')
            //F2206  589  FORMAT (1X,I5,',',15(e18.10,',')) 
            //F2207 		
            //F2208 
            //F2209 ! *******************************************************************************************
            //F2210 ! End halocarbon write
            //F2211 ! *******************************************************************************************
            //F2212 
            //F2213       ENDIF
        }
        //F2214 
        //F2215 !*******************************************************************
        //F2216 !*******************************************************************
        //F2217 !	minicam csv output mrj 4/26/00
        //F2218 !   revised for TAR vsn 5/03 mrj
        //F2219 
        //F2220 
        //F2221 !**** following code is from mag.csv (formerly mag.out) computations
        //F2222         QQQCO2R  = (qco2(226)     +qco2(225))     /2.
        //UNNECESSARY const float QQQCO2R = ( FORCE.QCO2[ 226 ] + FORCE.QCO2[ 225 ] ) / 2.0;
        //F2223         QQQMR    = (QM(226)       +QM(225))       /2.
        //UNNECESSARY const float QQQMR = ( FORCE.QM[ 226 ] + FORCE.QM[ 225 ] ) / 2.0;
        //F2224         QQQNR    = (QN(226)       +QN(225))       /2.
        //UNNECESSARY const float QQQNR = ( FORCE.QN[ 226 ] + FORCE.QN[ 225 ] ) / 2.0;
        //F2225         QQQCFCR  = (QCFC(226)     +QCFC(225))     /2.
        //UNNECESSARY const float QQQCFCR = ( FORCE.QCFC[ 226 ] + FORCE.QCFC[ 225 ] ) / 2.0;
        //F2226         QQQSO2R  = (QSO2SAVE(226) +QSO2SAVE(225)) /2.
        //UNNECESSARY const float QQQSO2R = ( STOREDVALS.QSO2SAVE[ 226 ] + STOREDVALS.QSO2SAVE[ 225 ] ) / 2.0;
        //F2227         QQQDIRR  = (QDIRSAVE(226) +QDIRSAVE(225)) /2.
        //UNNECESSARY const float QQQDIRR = ( STOREDVALS.QDIRSAVE[ 226 ] + STOREDVALS.QDIRSAVE[ 225 ] ) / 2.0;
        //F2228         QQQFOCR  = (QFOC(226)     +QFOC(225))     /2.
        //UNNECESSARY const float QQQFOCR = ( JSTART.QFOC[ 226 ] + JSTART.QFOC[ 225 ] ) / 2.0;
        //F2229         QQQOZR   = QOZ(226)
        //UNNECESSARY const float QQQOZR = TANDSL.QOZ[ 226 ];
        //F2230         QQQBIOR  = (QBIO(226)     +QBIO(225))     /2.
        //UNNECESSARY const float QQQBIOR = ( TANDSL.QBIO[ 226 ] + TANDSL.QBIO[ 225 ] ) / 2.0;
        //F2231         QQQMO3R  = (QCH4O3(226)   +QCH4O3(225))   /2.
        //UNNECESSARY const float QQQMO3R = ( FORCE.QCH4O3[ 226 ] + FORCE.QCH4O3[ 225 ] ) / 2.0;
        //F2232         QQRSTROZ = (QSTRATOZ(226) +QSTRATOZ(225)) /2.
        //UNNECESSARY const float QQRSTROZ = ( FORCE.QSTRATOZ[ 226 ] + FORCE.QSTRATOZ[ 225 ] ) / 2.0;
        //F2233         QQRKYMAG = (QKYMAG(226)   +QKYMAG(225))   /2.
        //UNNECESSARY const float QQRKYMAG = ( JSTART.QKYMAG[ 226 ] + JSTART.QKYMAG[ 225 ] ) / 2.0;
        //F2234         QQRMONT  = (QMONT(226)    +QMONT(225))    /2.
        //UNNECESSARY const float QQRMONT = ( FORCE.QMONT[ 226 ] + FORCE.QMONT[ 225 ] ) / 2.0;
        //F2235         QQROTHER = (QOTHER(226)   +QOTHER(225))   /2.
        //UNNECESSARY const float QQROTHER = ( FORCE.QOTHER[ 226 ] + FORCE.QOTHER[ 225 ] ) / 2.0;
        //F2236 !*** end code block
        //F2237 
        //F2238 	OPEN (UNIT=9, file='./outputs/MAGOUT.CSV')

        // GetForcing now relies on globals, and these need to be set
        setGlobals( &CARB, &TANDSL, &CONCS, &NEWCONCS, 
                   &STOREDVALS, &NEWPARAMS, &BCOC, 
                   &METH1, &CAR, &FORCE, &JSTART,
                   &QADD, &HALOF );

        ofstream outfile9;
        openfile_write( &outfile9, "./magout_c.csv", DEBUG_IO ); //FIX filename
        //F2239 
        //F2240   100 FORMAT(I5,1H,,27(F15.5,1H,))
        //F2241 
        //F2242   101 FORMAT('Year,Temp,CO2Conc,CH4Conc,N2OConc,', &
        //F2243       'FcCO2,FcCH4,FcN2O,FcHALOS,FcTROPO3,', &
        //F2244       'FcSO4DIR,FcSO4IND,FcBIOAER,', &
        //F2245       'FcTOTAL,FOSSCO2,NETDEFOR,CH4Em,N2OEm,SO2-REG1,SO2-REG2,SO2-REG3,', &
        //F2246       'SeaLevel,FcKyoto,FcHFC,FcCFC+SF6,FcCH4H2O, FcBC, FcOC')
        //F2247 
        //F2248 	WRITE(9,101)  !header row
        outfile9 << "Year,CO2 Concentration,CH4 Concentration,N2O Concentration,Total RF,Kyoto Gas RF,";
        outfile9 << "CO2 RF,CH4 RF,N2O RF,Direct SO2 RF,Indirect SO2 RF,Fossil Fuel and Industrial CO2 Emissions,";
        outfile9 << "Land Use Change CO2 Emissions,CH4 Emissions,N2O Emissions,SO2 Emissions (Reg.1),";
        outfile9 << "SO2 Emissions (Reg.2),SO2 Emissions (Reg.3),Global Mean Temperature Rise,Sea Level Rise,";
        outfile9 << "Fossil BC/OC RF,Bio Burn RF" << endl;
        //F2249 
        //F2250         IIPRT=5	! sjs -- changed to 5 year interval in order to save more data points
        const int IIPRT = 5;
        //F2251         DO K=1990,IYEND,IIPRT
        for( int K=1990; K<=IYEND; K+=IIPRT ) {
            //F2252           IYR = K-1990+226
            int IYR = K - 1990 + 226;
            //F2253           IYRP=IYR-1
            //UNUSED int IYRP = IYR - 1;
            int yrindex = ( K - 1990 ) / IIPRT + 1;
            //F2254 
            //F2255 ! code to pass these items to MiniCAM in Results array
            //F2256 
            //F2257 	 MAGICCCResults(0,(K-1990)/IIPRT+1) = Float(K)
            MAGICCCResults[ 0 ][ yrindex ] = float( K );
            // Concentrations
            //F2259 	 MAGICCCResults(2,(K-1990)/IIPRT+1) = CO2(IYR)
            MAGICCCResults[ 1 ][ yrindex ] = CARB.CO2[ IYR ];
            //F2260 	 MAGICCCResults(3,(K-1990)/IIPRT+1) = CH4(IYR)
            MAGICCCResults[ 2 ][ yrindex ] = CONCS.CH4[ IYR ];
            //F2261 	 MAGICCCResults(4,(K-1990)/IIPRT+1) = CN2O(IYR)
            MAGICCCResults[ 3 ][ yrindex ] = CONCS.CN2O[ IYR ];

            // RADIATIVE FORCING
            //F2273 	 MAGICCCResults(13,(K-1990)/IIPRT+1) = GETFORCING( 0, K ) ! Total antro forcing
            MAGICCCResults[ 4 ][ yrindex ] = GETFORCING( 0, K );
            //F2282 	 MAGICCCResults(22,(K-1990)/IIPRT+1) = & !Kyoto Forcing
            //F2283 	    GETFORCING( 1, K ) + GETFORCING( 2, K )  + GETFORCING( 3, K ) + & ! CO2, CH4, and N2O
            //F2284 	    GETFORCING( 4, K ) + GETFORCING( 9, K ) + GETFORCING( 10, K ) + &! Long-lived F-gases
            //F2285 	    GETFORCING( 5, K ) + GETFORCING( 6, K ) + GETFORCING( 7, K ) + &
            //F2286 	    GETFORCING( 8, K ) + GETFORCING( 11, K ) + GETFORCING( 12, K ) ! Shorter-lived F-gases
            MAGICCCResults[ 5 ][ yrindex ] = GETFORCING( 1, K ) + GETFORCING( 2, K ) + GETFORCING( 3, K ) +
                GETFORCING( 4, K ) + GETFORCING( 9, K ) + GETFORCING( 10, K ) +
                GETFORCING( 5, K ) + GETFORCING( 6, K ) + GETFORCING( 7, K ) +
                GETFORCING( 8, K ) + GETFORCING( 11, K ) + GETFORCING( 12, K );
            //F2262 	 MAGICCCResults(5,(K-1990)/IIPRT+1) = GETFORCING( 1, K ) ! CO2
            MAGICCCResults[ 6 ][ yrindex ] = GETFORCING( 1, K );
            //F2263 	 MAGICCCResults(6,(K-1990)/IIPRT+1) = GETFORCING( 2, K ) ! CH4 (no indirect components)
            MAGICCCResults[ 7 ][ yrindex ] = GETFORCING( 2, K );
            //F2264 	 MAGICCCResults(7,(K-1990)/IIPRT+1) = GETFORCING( 3, K ) ! N2O
            MAGICCCResults[ 8 ][ yrindex ] = GETFORCING( 3, K );
            //F2270 	 MAGICCCResults(10,(K-1990)/IIPRT+1) = GETFORCING( 14, K ) ! SO2 direct only
            MAGICCCResults[ 9 ][ yrindex ] = GETFORCING( 14, K );
            //F2271 	 MAGICCCResults(11,(K-1990)/IIPRT+1) = GETFORCING( 13, K ) - GETFORCING( 14, K ) ! indirect only
            MAGICCCResults[ 10 ][ yrindex ] = GETFORCING( 13, K ) - GETFORCING( 14, K );

            // EMISSIONS
            //F2274 	 MAGICCCResults(14,(K-1990)/IIPRT+1) = EF(IYR)
            MAGICCCResults[ 11 ][ yrindex ] = CARB.EF.getval( IYR );
            //F2275 	 MAGICCCResults(15,(K-1990)/IIPRT+1) = EDNET(IYR)
            MAGICCCResults[ 12 ][ yrindex ] = METH1.ednet.getval( IYR );
            //F2276 	 MAGICCCResults(16,(K-1990)/IIPRT+1) = ECH4(IYR)
            MAGICCCResults[ 13 ][ yrindex ] = CONCS.ECH4.getval( IYR );
            //F2277 	 MAGICCCResults(17,(K-1990)/IIPRT+1) = EN2O(IYR)
            MAGICCCResults[ 14 ][ yrindex ] = CONCS.EN2O.getval( IYR );
            //F2278 	 MAGICCCResults(18,(K-1990)/IIPRT+1) = ESO21(IYR)-ES1990
            MAGICCCResults[ 15 ][ yrindex ] = CONCS.ESO21.getval( IYR ) - Sulph.ES1990;
            //F2279 	 MAGICCCResults(19,(K-1990)/IIPRT+1) = ESO22(IYR)-ES1990
            MAGICCCResults[ 16 ][ yrindex ] = CONCS.ESO22.getval( IYR ) - Sulph.ES1990;
            //F2280 	 MAGICCCResults(20,(K-1990)/IIPRT+1) = ESO23(IYR)-ES1990
            MAGICCCResults[ 17 ][ yrindex ] = CONCS.ESO23.getval( IYR ) - Sulph.ES1990;

            // TEMPERATURE AND SEA LEVEL RISE
            //F2258 	 MAGICCCResults(1,(K-1990)/IIPRT+1) = TEMUSER(IYR)+TGAV(226)
            MAGICCCResults[ 18 ][ yrindex ] = STOREDVALS.TEMUSER[ IYR ] + TANDSL.TGAV[ 226 ];
            //F2281 	 MAGICCCResults(21,(K-1990)/IIPRT+1) = getSLR( IYR ) ! getSLR is external fn with acutal year as argument
            MAGICCCResults[ 19 ][ yrindex ] = getSLR( K );

            // BC/OC FORCING
            //F2293 	 MAGICCCResults(26,(K-1990)/IIPRT+1) = GETFORCING( 24, K )	! BC forcing 
         //   MAGICCCResults[ 20 ][ yrindex ] = GETFORCING( 24, K );
            //F2294 	 MAGICCCResults(27,(K-1990)/IIPRT+1) = GETFORCING( 25, K )	! OC forcing 
         //   MAGICCCResults[ 21 ][ yrindex ] = GETFORCING( 25, K );
            // Fossil BC/OC Forcing
            MAGICCCResults[ 20 ][ yrindex ] = GETFORCING( 28, K );
            // Biomass Burning Aerosol Forcing
            MAGICCCResults[ 21 ][ yrindex ] = GETFORCING( 20, K );
            
            //F2295 
            //F2296 ! now we can write stuff out
            //F2297 
            //F2298 	   WRITE (9,100) K,MAGICCCResults(1:25,(K-1990)/IIPRT+1)
            outfile9 << K << ", ";
            for( int i=1; i<=21; i++ )
                outfile9 << MAGICCCResults[ i ][ yrindex ] << ",";
            outfile9 << endl;
            //F2299 
            //F2300 	END DO
        }
        //F2301 
        //F2302 	WRITE(9,*)
        outfile9 << endl;
        //F2303 !     ******* END MINICAM OUTPUT ***
        //F2304 	CLOSE (9)
        outfile9.close();
        //F2305 !
        //F2306 !
        //F2307 !  ***************************************************************
        //F2308 !
        //F2309 !  WRITE CONCENTRATIONS CONSISTENT WITH USER CLIMATE MODEL AND WITH
        //F2310 !   THE CORRESPONDING SET OF RESULTS IN THE MAG.OUT CONCS
        //F2311 !   DISPLAY FILE.
        //F2312 !
        //F2313       IF(ISCENGEN.EQ.9.OR.NSIM.EQ.4)THEN
        if( NSIM.ISCENGEN == 9 || NSIM.NSIM == 4 ) {
            //F2314 !
            //F2315       open(unit=9,file='./outputs/concs.dis',status='UNKNOWN')
            openfile_write( &outfile9, "./concs_c.dis", DEBUG_IO );
            //F2316 !
            //F2317         WRITE (9,211)
            outfile9 << "YEAR CO2USER   CO2LO  CO2MID   CO2HI CH4USER   CH4LO  CH4MID   CH4HI     N2O MIDTAUCH4" << endl;
            //F2318 !
            //F2319 !  PRINTOUT INTERVAL FOR DISPLAY PURPOSES SET BY IDIS IN MAGEXTRA.CFG
            //F2320 !
            //F2321         DO K=1,KEND,IDIS
            for( int K=1; K<=Limits.KEND; K += IDIS ) {
                float TOR;
                //F2322 !
                //F2323 !  CONVERT END OF YEAR TO MIDYEAR CONCS
                //F2324 !
                //F2325           CO2MID =(CO2(K)+CO2(K-1))/2.
                float CO2MID = ( CARB.CO2[ K ] + CARB.CO2[ K-1 ] ) / 2.0;
                //F2326           CH4MID =(CH4(K)+CH4(K-1))/2.
                float CH4MID = ( CONCS.CH4[ K ] + CONCS.CH4[ K-1 ] ) / 2.0;
                //F2327           CN2OMID=(CN2O(K)+CN2O(K-1))/2.
                float CN2OMID = ( CONCS.CN2O[ K ] + CONCS.CN2O[ K-1 ] ) / 2.0;
                //F2328 !
                //F2329           IYEAR=1764+K
                int IYEAR = 1764 + K;
                //F2330 !
                //F2331           IF(K.GE.226)THEN
                float CH4LMID, CH4BMID, CH4HMID, CO2LMID, CO2BMID, CO2HMID;
                if( K >= 226 ) {
                    //F2332 !
                    //F2333             CH4LMID=(CH4L(K)+CH4L(K-1))/2.
                    CH4LMID = ( METH1.ch4l.getval( K ) + METH1.ch4l.getval( K-1 ) ) / 2.0;
                    //F2334             CH4BMID=(CH4B(K)+CH4B(K-1))/2.
                    CH4BMID = ( METH1.ch4b.getval( K ) + METH1.ch4b.getval( K-1 ) ) / 2.0;
                    //F2335             CH4HMID=(CH4H(K)+CH4H(K-1))/2.
                    CH4HMID = ( METH1.ch4h.getval( K ) + METH1.ch4h.getval( K-1 ) ) / 2.0;
                    //F2336 !
                    //F2337             CO2LMID=(CCO2(1,K)+CCO2(1,K-1))/2.
                    CO2LMID = ( CARB.CCO2.getval( 1, K ) + CARB.CCO2.getval( 1, K-1 ) ) / 2.0;
                    //F2338             CO2BMID=(CCO2(2,K)+CCO2(2,K-1))/2.
                    CO2BMID = ( CARB.CCO2.getval( 2, K ) + CARB.CCO2.getval( 2, K-1 ) ) / 2.0;
                    //F2339             CO2HMID=(CCO2(3,K)+CCO2(3,K-1))/2.
                    CO2HMID = ( CARB.CCO2.getval( 3, K ) + CARB.CCO2.getval( 3, K-1 ) ) / 2.0;
                    //F2340 !
                    //F2341 !  ADD CORRECTIONS TO LO, MID, HI CO2 TO FIT OBSERVED DATA IN 2000
                    //F2342 !
                    //F2343             IF(ICO2CORR.EQ.1)THEN
                    if( JSTART.ICO2CORR == 1 ) {
                        //F2344               CO2LMID=CO2LMID+CORREN1
                        CO2LMID += CORREN.CORREN1;
                        //F2345               CO2BMID=CO2BMID+CORREN2
                        CO2BMID += CORREN.CORREN2;
                        //F2346               CO2HMID=CO2HMID+CORREN3
                        CO2HMID += CORREN.CORREN3;
                        //F2347             ENDIF
                    }
                    //F2348 !
                    //F2349             IF(K.LE.236)THEN
                    if( K <= 236 ) {
                        //F2350               CO2LMID=CO2MID
                        //F2351               CO2BMID=CO2MID
                        //F2352               CO2HMID=CO2MID
                        CO2LMID = CO2BMID = CO2HMID = CO2MID;
                        //F2353             ENDIF
                    }
                    //F2354 !
                    //F2355 !  DEFINE LOW, MID AND HIGH CH4 VALUES OVER 1991 TO JSTART YEAR
                    //F2356 !
                    //F2357             IF(K.LT.236)THEN
                    if( K < 236 ) {
                        //F2358               CH4LMID=CH4MID
                        //F2359               CH4BMID=CH4MID
                        //F2360               CH4HMID=CH4MID
                        CH4LMID = CH4BMID = CH4HMID = CH4MID;
                        //F2361             ENDIF
                    }
                    //F2362 !
                    //F2363 !  SPECIFY METHANE LIFETIME OUTPUT (CENTRAL VALUE)
                    //F2364 !
                    //F2365             IF(K.LE.236)THEN
                    if( K <= 236 ) {
                        //F2366               TOR=TORREF
                        TOR = TORREF;
                        //F2367             ELSE
                    } else {
                        //F2368               TOR=TCH4(K)
                        TOR = METH1.TCH4[ K ];
                        //F2369             ENDIF
                    }
                    //F2370 !
                    //F2371           ENDIF
                }
                //F2372 !
                //F2373           IF(IYEAR.LT.1990) &
                if( IYEAR < 1990 )
                    //F2374             WRITE (9,223) IYEAR,CO2MID,CO2MID,CO2MID,CO2MID, &
                    //F2375             CH4MID,CH4MID,CH4MID,CH4MID,CN2OMID
                    outfile9 << IYEAR << "," << CO2MID << "," << CO2MID << "," << CO2MID << "," << CO2MID << "," << CH4MID << "," << CH4MID << "," << CH4MID << "," << CH4MID << "," << CN2OMID << endl;
                
                //F2376           IF(IYEAR.GE.1990)THEN
                if( IYEAR >= 1990 ) {
                    //F2377             IF(IYEAR.LE.2000)THEN
                    if( IYEAR <= 2000 ) {
                        //F2378               CO2LMID=CO2MID
                        //F2379               CO2BMID=CO2MID
                        //F2380               CO2HMID=CO2MID
                        CO2LMID = CO2BMID = CO2HMID = CO2MID;
                        //F2381             ENDIF
                    }
                    //F2382             WRITE (9,224) IYEAR,CO2MID,CO2LMID,CO2BMID,CO2HMID, &
                    //F2383             CH4MID,CH4LMID,CH4BMID,CH4HMID,CN2OMID,TOR
                    outfile9 << IYEAR << "," << CO2MID << "," << CO2LMID << "," << CO2BMID << "," << CO2HMID << "," << CH4MID << "," << CH4LMID << "," << CH4BMID << "," << CH4HMID << "," << CN2OMID << "," << TOR << endl;
                    //F2384           ENDIF
                }
                //F2385 !
                //F2386         END DO
            }
            //F2387 !
            //F2388         WRITE (9,211)
            outfile9 << "YEAR CO2USER   CO2LO  CO2MID   CO2HI CH4USER   CH4LO  CH4MID   CH4HI     N2O MIDTAUCH4" << endl;
            //F2389 !
            //F2390       CLOSE(9)
            outfile9.close();
            //F2391 !
            //F2392 !  ************************************************************
            //F2393 !
            //F2394 !  WRITE FORCING CHANGES FROM MID-1990 TO MAG DISPLAY FILE
            //F2395 !
            //F2396       open(unit=9,file='./outputs/forcings.dis',status='UNKNOWN')
            openfile_write( &outfile9, "./forcings_c.dis", DEBUG_IO );
            //F2397 !
            //F2398         WRITE (9,57)
            outfile9 << "YEAR,CO2,CH4tot,N2O, HALOtot,TROPOZ,SO4DIR,SO4IND,BIOAER,FOC+FBC,QAERMN,QLAND, TOTAL, YEAR,CH4-O3," << endl;
            //F2399 !
            //F2400 !  PRINTOUT INTERVAL FOR DISPLAY PURPOSES SET BY IDIS IN MAGEXTRA.CFG
            //F2401 !
            //F2402         DO K=1990,IYEND,IDIS
            for( int K=1990; K<=IYEND; K += IDIS ) {
                //F2403           IYR = K-1990+226
                int IYR = K - 1990 + 226;
                //F2404           IYRP=IYR-1
                int IYRP = IYR - 1;
                //F2405 !
                //F2406           DELQCO2 = (QCO2(IYR)+QCO2(IYRP))/2.-QQQCO2R
                float DELQCO2 = ( FORCE.QCO2[ IYR ] + FORCE.QCO2[ IYRP ] ) / 2.0 - QQQCO2R;
                //F2407           DELQM   = (QM(IYR)+QM(IYRP))/2.    -QQQMR
                float DELQM = ( FORCE.QM[ IYR ] + FORCE.QM[ IYRP ] ) / 2.0 - QQQMR;
                //F2408           DELQN   = (QN(IYR)+QN(IYRP))/2.    -QQQNR
                float DELQN = ( FORCE.QN[ IYR ] + FORCE.QN[ IYRP ] ) / 2.0 - QQQNR;
                //F2409           DELQCFC = (QCFC(IYR)+QCFC(IYRP))/2.-QQQCFCR
                float DELQCFC = ( FORCE.QCFC[ IYR ] + FORCE.QCFC[ IYRP ] ) / 2.0 - QQQCFCR;
                //F2410 !
                //F2411 !  NOTE : DELQSO2 INCLUDES QFOC
                //F2412 !
                //F2413           DELQSO2 = (QSO2SAVE(IYR)+QSO2SAVE(IYRP))/2.-QQQSO2R
                float DELQSO2 = ( STOREDVALS.QSO2SAVE[ IYR ] + STOREDVALS.QSO2SAVE[ IYRP ] ) / 2.0 - QQQSO2R;
                //F2414           DELQDIR = (QDIRSAVE(IYR)+QDIRSAVE(IYRP))/2.-QQQDIRR
                float DELQDIR = ( STOREDVALS.QDIRSAVE[ IYR ] + STOREDVALS.QDIRSAVE[ IYRP ] ) / 2.0 - QQQDIRR;
                //F2415           DELQIND = DELQSO2-DELQDIR
                float DELQIND = DELQSO2 - DELQDIR;
                //F2416           DELQFOC = (QFOC(IYR)+QFOC(IYRP))/2.-QQQFOCR
                float DELQFOC = ( JSTART.QFOC[ IYR ] + JSTART.QFOC[ IYRP ] ) / 2.0 - QQQFOCR;
                //F2417           DELQD   = DELQDIR-DELQFOC
                float DELQD = DELQDIR - DELQFOC;
                //F2418 !
                //F2419 ! NOTE SPECIAL CASE FOR QOZ BECAUSE OF NONLINEAR CHANGE OVER 1990
                //F2420 !
                //F2421           IF(IYR.EQ.226)THEN
                float QOZMID;
                if( IYR == 226 )
                    //F2422             QOZMID= QOZ(IYR)
                    QOZMID = TANDSL.QOZ[ IYR ];
                //F2423           ELSE
                else
                    //F2424             QOZMID= (QOZ(IYR)+QOZ(IYRP))/2.
                    QOZMID = ( TANDSL.QOZ[ IYR ] + TANDSL.QOZ[ IYRP ] ) / 2.0;
                //F2425           ENDIF
                //F2426           DELQOZ  = QOZMID-QQQOZR
                float DELQOZ = QOZMID - QQQOZR;
                //F2427 !
                //F2428           DELQMN  = (QMN(IYR)+QMN(IYRP))/2.-QQQMNR
                float DELQMN = ( TANDSL.QMN[ IYR ] + TANDSL.QMN[ IYRP ] ) / 2.0 - QQQMNR;
                //F2429           DELQLAND= (QLAND(IYR)+QLAND(IYRP))/2.-QQQLANDR
                float DELQLAND = ( TANDSL.QLAND[ IYR ] + TANDSL.QLAND[ IYRP ] ) / 2.0 - QQQLANDR;
                //F2430           DELQBIO = (QBIO(IYR)+QBIO(IYRP))/2.-QQQBIOR
                float DELQBIO = ( TANDSL.QBIO[ IYR ] + TANDSL.QBIO[ IYRP ] ) / 2.0 - QQQBIOR;
                //F2431           DELQTOT = DELQCO2+DELQM+DELQN+DELQCFC+DELQSO2+DELQBIO &
                //F2432           +DELQOZ+DELQLAND+DELQMN
                float DELQTOT = DELQCO2 + DELQM + DELQN + DELQCFC + DELQSO2 + DELQBIO + DELQOZ + DELQLAND + DELQMN;
                //F2433 !
                //F2434           DQCH4O3 = (QCH4O3(IYR)+QCH4O3(IYRP))/2.-QQQMO3R
                float DQCH4O3 = ( FORCE.QCH4O3[ IYR ] + FORCE.QCH4O3[ IYRP ] ) / 2.0 - QQQMO3R;
                //F2435           DELQM   = DELQM-DQCH4O3
                DELQM -= DQCH4O3;
                //F2436           DELQOZ  = DELQOZ+DQCH4O3
                DELQOZ += DQCH4O3;
                //F2437           DELSTROZ= (QSTRATOZ(IYR)+QSTRATOZ(IYRP))/2.-QQRSTROZ
                float DELSTROZ = ( FORCE.QSTRATOZ[ IYR ] + FORCE.QSTRATOZ[ IYRP ] ) / 2.0 - QQRSTROZ;
                //F2438           IF(IO3FEED.EQ.0)DELSTROZ=0.0
                if( METH1.IO3FEED == 0 ) DELSTROZ = 0.0;
                //F2439 !
                //F2440           DELKYMAG = (QKYMAG(IYR) +QKYMAG(IYRP))  /2.-QQRKYMAG
                float DELKYMAG = ( JSTART.QKYMAG[ IYR ] + JSTART.QKYMAG[ IYRP ] ) / 2.0 - QQRKYMAG;
                //F2441           DELMONT  = (QMONT(IYR)  +QMONT(IYRP))   /2.-QQRMONT
                float DELMONT = ( FORCE.QMONT[ IYR ] + FORCE.QMONT[ IYRP ] ) / 2.0 - QQRMONT;
                //F2442           DELOTHER = (QOTHER(IYR) +QOTHER(IYRP))  /2.-QQROTHER
                float DELOTHER = ( FORCE.QOTHER[ IYR ] + FORCE.QOTHER[ IYRP ] ) / 2.0 - QQROTHER;
                //F2443           DELKYOTO = DELKYMAG+DELOTHER
                float DELKYOTO = DELKYMAG + DELOTHER;
                //F2444 !
                //F2445           WRITE(9,571)K,DELQCO2,DELQM,DELQN,DELQCFC,DELQOZ, &
                //F2446           DELQD,DELQIND,DELQBIO,DELQFOC,DELQMN,DELQLAND,DELQTOT,K, &
                //F2447           DQCH4O3,DELSTROZ,DELMONT,DELKYOTO
                outfile9 << K << "," << DELQCO2 << "," << DELQM << "," << DELQN << "," << DELQCFC << "," << DELQOZ << "," 
                << DELQD << "," << DELQIND << "," << DELQBIO << "," << DELQFOC << "," << DELQMN << "," << DELQLAND << "," 
                << DELQTOT << "," << K << "," << DQCH4O3 << "," << DELSTROZ << "," << DELMONT << "," << DELKYOTO << endl;
                //F2448         END DO
            }
            //F2449 !
            //F2450         WRITE (9,57)
            outfile9 << "YEAR,CO2,CH4tot,N2O, HALOtot,TROPOZ,SO4DIR,SO4IND,BIOAER,FOC+FBC,QAERMN,QLAND, TOTAL, YEAR,CH4-O3," << endl;
            //F2451 !
            //F2452       CLOSE(9)
            outfile9.close();
            //F2453 !
            //F2454       ENDIF
        }
        //F2455 !
        //F2456 !  **************************************************************
        //F2457 !
        //F2458 !  end of NSIM loop
        //F2459 !
        //F2460   1   CONTINUE
    }
    //F2461 !
    //F2462 !  **************************************************************
    //F2463 !
    //F2464 !  SCALE TEMPERATURES TO GET DRIVER TEMPERATURES, XSO2i. THIS IS IN
    //F2465 !   A LONG IF STATEMENT, BRACKETTED BY THE TRIPLE '******' LINES.
    //F2466 !
    //F2467       IF(ISCENGEN.EQ.1)THEN
    if( NSIM.ISCENGEN == 1 ) {
        //F2468 !
        //F2469         DO 565 NCLIM=1,4
        for( NSIM.NCLIM=1; NSIM.NCLIM<=4; NSIM.NCLIM++ ) {
            //F2470 !
            //F2471           DO K=197,KEND
            for( int K=197; K<=Limits.KEND; K++ ) {
                //F2472             KSG=K-196
                int KSG = K - 196;
                //F2473             iiiyr=ksg+1960
                //UNUSED int IIIYR = KSG + 1960;
                //F2474             XTSUM=TALL(NCLIM,KSG)-TGHG(NCLIM,KSG)
                float XTSUM = TALL[ NSIM.NCLIM ][ KSG ] - TGHG[ NSIM.NCLIM ][ KSG ];
                //F2475             YTSUM=TSO21(NCLIM,KSG)+TSO22(NCLIM,KSG)+TSO23(NCLIM,KSG)
                float YTSUM = TSO21[ NSIM.NCLIM ][ KSG ] + TSO22[ NSIM.NCLIM ][ KSG ] + TSO23[ NSIM.NCLIM ][ KSG ];
                //F2476 !
                //F2477             SCALER(K)=1.0
                SCALER[ K ] = 1.0;
                //F2478             IF(YTSUM.NE.0.0)THEN
                if( YTSUM != 0.0 ) {
                    //F2479               SCALER(K)=XTSUM/YTSUM
                    SCALER[ K ] = XTSUM / YTSUM;
                    //F2480               IF(SCALER(K).GT.2.0)SCALER(K)=2.0
                    if( SCALER[ K ] > 2.0 ) SCALER[ K ] = 2.0;
                    //F2481               IF(SCALER(K).LT.0.0)SCALER(K)=0.0
                    if( SCALER[ K ] < 0.0 ) SCALER[ K ] = 0.0;
                    //F2482             ENDIF
                }
                //F2483             XGHG(NCLIM,KSG) =XGHG(NCLIM,KSG)*SCALER(K)
                XGHG[ NSIM.NCLIM ][ KSG ] *= SCALER[ K ];
                //F2484 !
                //F2485           END DO
            }
            //F2486 !
            //F2487 !  SMOOTHED VERSION OF METHOD 1
            //F2488 !
            //F2489           ISMOOTH=1
            const int ISMOOTH = 1;
            //F2490           IF(ISMOOTH.EQ.1)THEN
            if( ISMOOTH == 1 ) {
                //F2491             DO K=200,KEND-3
                for( int K=200; K<=Limits.KEND-3; K++ ) {
                    //F2492               SS1=SCALER(K-3)
                    float SS1 = SCALER[ K-3 ];
                    //F2493               SS2=SCALER(K-2)
                    float SS2 = SCALER[ K-2 ];
                    //F2494               SS3=SCALER(K-1)
                    float SS3 = SCALER[ K-1 ];
                    //F2495               SS4=SCALER(K-0)
                    float SS4 = SCALER[ K-0 ];
                    //F2496               SS5=SCALER(K+1)
                    float SS5 = SCALER[ K+1 ];
                    //F2497               SS6=SCALER(K+2)
                    float SS6 = SCALER[ K+2 ];
                    //F2498               SS7=SCALER(K+3)
                    float SS7 = SCALER[ K+3 ];
                    //F2499 !              SCALAR(K)=SS1+SS7+6.0*(SS2+SS6)+15.0*(SS3+SS5)+20.0*SS4
                    //F2500 !              SCALAR(K)=SCALAR(K)/64.0
                    //F2501               SCALAR(K)=SS1+SS7+SS2+SS6+SS3+SS5+SS4
                    SCALAR[ K ] = SS1 + SS7 + SS2 + SS6 + SS3 + SS5 + SS4;
                    //F2502               SCALAR(K)=SCALAR(K)/7.0
                    SCALAR[ K ] /= 7.0;
                    //F2503               KSG=K-196
                    int KSG = K - 196;
                    //F2504                 XSO21(NCLIM,KSG)=TSO21(NCLIM,KSG)*SCALAR(K)
                    XSO21[ NSIM.NCLIM ][ KSG ] = TSO21[ NSIM.NCLIM ][ KSG ] * SCALAR[ K ];
                    //F2505                 XSO22(NCLIM,KSG)=TSO22(NCLIM,KSG)*SCALAR(K)
                    XSO22[ NSIM.NCLIM ][ KSG ] = TSO22[ NSIM.NCLIM ][ KSG ] * SCALAR[ K ];
                    //F2506                 XSO23(NCLIM,KSG)=TSO23(NCLIM,KSG)*SCALAR(K)
                    XSO23[ NSIM.NCLIM ][ KSG ] = TSO23[ NSIM.NCLIM ][ KSG ] * SCALAR[ K ];
                    //F2507             END DO
                }
                //F2508 !
                //F2509             DO K=197,199
                for( int K=197; K<=199; K++ )
                    //F2510               SCALAR(K)=SCALER(K)
                    SCALAR[ K ] = SCALER[ K ];
                //F2511             END DO
                //F2512 !
                //F2513             DO K=KEND-2,KEND
                for( int K=Limits.KEND-2; K<=Limits.KEND; K++ )
                    //F2514               SCALAR(K)=SCALER(K)
                    SCALAR[ K ] = SCALER[ K ];
                //F2515             END DO
                //F2516           ENDIF
            }
            //F2517 !
            //F2518  565    CONTINUE
        }
        //F2519 !
        //F2520 !  **************************************************************
        //F2521 !
        //F2522 !  WRITE TEMPERATURES TO OLD AND NEW SCENGEN DRIVER FILES.
        //F2523 !   UNSCALED TEMPS GO TO OLD FILES, SCALED TEMPS TO NEW FILES.
        //F2524 !
        //F2525         OPEN(UNIT=10,file='./outputs/lodrive.raw' ,STATUS='UNKNOWN')
        ofstream outfile10, outfile11, outfile12, outfile13, outfile14, outfile15, outfile16, outfile17;
        openfile_write( &outfile10, "./lodrive.raw", DEBUG_IO );
        //F2526         OPEN(UNIT=11,file='./outputs/middrive.raw',STATUS='UNKNOWN')
        openfile_write( &outfile11, "./middrive.raw", DEBUG_IO );
        //F2527         OPEN(UNIT=12,file='./outputs/hidrive.raw' ,STATUS='UNKNOWN')
        openfile_write( &outfile12, "./hidrive.raw", DEBUG_IO );
        //F2528         OPEN(UNIT=13,file='./outputs/usrdrive.raw',STATUS='UNKNOWN')
        openfile_write( &outfile13, "./usrdrive.raw", DEBUG_IO );
        //F2529 !
        //F2530         OPEN(UNIT=14,file='./outputs/lodrive.out' ,STATUS='UNKNOWN')
        openfile_write( &outfile14, "./lodrive.out", DEBUG_IO );
        //F2531         OPEN(UNIT=15,file='./outputs/middrive.out',STATUS='UNKNOWN')
        openfile_write( &outfile15, "./middrive.out", DEBUG_IO );
        //F2532         OPEN(UNIT=16,file='./outputs/hidrive.out' ,STATUS='UNKNOWN')
        openfile_write( &outfile16, "./hidrive.out", DEBUG_IO );
        //F2533         OPEN(UNIT=17,file='./outputs/usrdrive.out',STATUS='UNKNOWN')
        openfile_write( &outfile17, "./usrdrive.out", DEBUG_IO );
        //F2534 !
        //F2535         DO NCLIM=1,4
        for( NSIM.NCLIM=1; NSIM.NCLIM<=4; NSIM.NCLIM++ ) {
            //F2536         KSGL=KEND-196
            const int KSGL = Limits.KEND - 196;
            //F2537         DO KSG=1,KSGL
            for( int KSG=1; KSG<=KSGL; KSG++ ) {
                //F2538         KYY=KSG+1960
                int KYY = KSG + 1960;
                //F2539 !
                //F2540 !  BECAUSE OF SMOOTHING, TSOij AND XSOij ARE NOT DEFINED FOR THE
                //F2541 !   LAST 3 YEARS OF THE RUN. WE THEREFORE DEFINE THEN BY LINEAR
                //F2542 !   EXTRAPOLATION.
                //F2543 !
                //F2544         IF(KSG.EQ.KSGL-3)THEN
                float DT21, DT22, DT23, DX21, DX22, DX23;
                if( KSG == KSGL-3 ) {
                    //F2545           DT21=(TSO21(NCLIM,KSG)-TSO21(NCLIM,KSG-3))/3.0
                    DT21 = ( TSO21[ NSIM.NCLIM ][ KSG ] - TSO21[ NSIM.NCLIM ][ KSG-3 ] ) / 3.0;
                    //F2546           DT22=(TSO22(NCLIM,KSG)-TSO22(NCLIM,KSG-3))/3.0
                    DT22 = ( TSO22[ NSIM.NCLIM ][ KSG ] - TSO22[ NSIM.NCLIM ][ KSG-3 ] ) / 3.0;
                    //F2547           DT23=(TSO23(NCLIM,KSG)-TSO23(NCLIM,KSG-3))/3.0
                    DT23 = ( TSO23[ NSIM.NCLIM ][ KSG ] - TSO23[ NSIM.NCLIM ][ KSG-3 ] ) / 3.0;
                    //F2548           DX21=(XSO21(NCLIM,KSG)-XSO21(NCLIM,KSG-3))/3.0
                    DX21 = ( XSO21[ NSIM.NCLIM ][ KSG ] - XSO21[ NSIM.NCLIM ][ KSG-3 ] ) / 3.0;
                    //F2549           DX22=(XSO22(NCLIM,KSG)-XSO22(NCLIM,KSG-3))/3.0
                    DX22 = ( XSO22[ NSIM.NCLIM ][ KSG ] - XSO22[ NSIM.NCLIM ][ KSG-3 ] ) / 3.0;
                    //F2550           DX23=(XSO23(NCLIM,KSG)-XSO23(NCLIM,KSG-3))/3.0
                    DX23 = ( XSO23[ NSIM.NCLIM ][ KSG ] - XSO23[ NSIM.NCLIM ][ KSG-3 ] ) / 3.0;
                    //F2551         ENDIF
                }
                //F2552         IF(KSG.GT.KSGL-3)THEN
                if( KSG > KSGL-3 ) {
                    //F2553           KKK=KSG-(KSGL-3)
                    int KKK = KSG - (KSGL - 3 );
                    //F2554           TSO21(NCLIM,KSG)=TSO21(NCLIM,KSGL-3)+KKK*DT21
                    TSO21[ NSIM.NCLIM ][ KSG ] = TSO21[ NSIM.NCLIM ][ KSGL-3 ] + KKK * DT21;
                    //F2555           TSO22(NCLIM,KSG)=TSO22(NCLIM,KSGL-3)+KKK*DT22
                    TSO22[ NSIM.NCLIM ][ KSG ] = TSO22[ NSIM.NCLIM ][ KSGL-3 ] + KKK * DT22;
                    //F2556           TSO23(NCLIM,KSG)=TSO23(NCLIM,KSGL-3)+KKK*DT23
                    TSO21[ NSIM.NCLIM ][ KSG ] = TSO21[ NSIM.NCLIM ][ KSGL-3 ] + KKK * DT23;
                    //F2557           XSO21(NCLIM,KSG)=XSO21(NCLIM,KSGL-3)+KKK*DX21
                    XSO21[ NSIM.NCLIM ][ KSG ] = XSO21[ NSIM.NCLIM ][ KSGL-3 ] + KKK * DX21;
                    //F2558           XSO22(NCLIM,KSG)=XSO22(NCLIM,KSGL-3)+KKK*DX22
                    XSO22[ NSIM.NCLIM ][ KSG ] = XSO22[ NSIM.NCLIM ][ KSGL-3 ] + KKK * DX22;
                    //F2559           XSO23(NCLIM,KSG)=XSO23(NCLIM,KSGL-3)+KKK*DX23
                    XSO23[ NSIM.NCLIM ][ KSG ] = XSO23[ NSIM.NCLIM ][ KSGL-3 ] + KKK * DX23;
                    //F2560         ENDIF  
                }
                //F2561         IF(NCLIM.EQ.1)THEN
                if( NSIM.NCLIM == 1 ) {
                    //F2562           IF(KSG.EQ.1)THEN
                    if( KSG == 1 ) {
                        //F2563             WRITE(10,937)mnem
                        outfile10 << "PROFILE: " << mnem << endl;
                        //F2564             WRITE(10,930)
                        outfile10 << "LOW CLIMATE MODEL PARAMETERS" << endl;
                        //F2565             WRITE(10,934)
                        outfile10 << "YEAR       GHG     ESO21     ESO22     ESO23       ALL" << endl;
                        //F2566             WRITE(10,936)TREF(NCLIM)
                        outfile10 << "REF T " << TREF[ NSIM.NCLIM ] << endl;
                        //F2567             WRITE(14,937)mnem
                        outfile14 << "PROFILE: " << mnem << endl;
                        //F2568             WRITE(14,930)
                        outfile14 << "LOW CLIMATE MODEL PARAMETERS" << endl;
                        //F2569             WRITE(14,934)
                        outfile14 << "YEAR       GHG     ESO21     ESO22     ESO23       ALL" << endl;
                        //F2570             WRITE(14,936)TREF(NCLIM)
                        outfile14 << "REF T " << TREF[ NSIM.NCLIM ] << endl;
                        //F2571           ENDIF
                    }
                    //F2572           WRITE(10,935)KYY,TGHG(NCLIM,KSG),TSO21(NCLIM,KSG), &
                    //F2573           TSO22(NCLIM,KSG),TSO23(NCLIM,KSG),TALL(NCLIM,KSG)
                    outfile10 << KYY << "," << TGHG[ NSIM.NCLIM ][ KSG ] << "," << TSO21[ NSIM.NCLIM ][ KSG ] << "," << TSO22[ NSIM.NCLIM ][ KSG ] << "," 
                    << TSO23[ NSIM.NCLIM ][ KSG ] << "," << TALL[ NSIM.NCLIM ][ KSG ] << endl;
                    //F2574           WRITE(14,935)KYY,TGHG(NCLIM,KSG),XSO21(NCLIM,KSG), &
                    //F2575           XSO22(NCLIM,KSG),XSO23(NCLIM,KSG),TALL(NCLIM,KSG)
                    outfile14 << KYY << "," << TGHG[ NSIM.NCLIM ][ KSG ] << "," << TSO21[ NSIM.NCLIM ][ KSG ] << "," << TSO22[ NSIM.NCLIM ][ KSG ] << "," 
                    << TSO23[ NSIM.NCLIM ][ KSG ] << "," << TALL[ NSIM.NCLIM ][ KSG ] << endl;
                    //F2576         ENDIF
                }
                //F2577 !
                //F2578         IF(NCLIM.EQ.2)THEN
                if( NSIM.NCLIM == 2 ) {
                    //F2579           IF(KSG.EQ.1)THEN
                    if( KSG == 1 ) {
                        //F2580             WRITE(11,937)mnem
                        outfile11 << "PROFILE: " << mnem << endl;
                        //F2581             WRITE(11,931)
                        outfile11 << "MID CLIMATE MODEL PARAMETERS" << endl;
                        //F2582             WRITE(11,934)
                        outfile11 << "YEAR       GHG     ESO21     ESO22     ESO23       ALL" << endl;
                        //F2583             WRITE(11,936)TREF(NCLIM)
                        outfile11 << "REF T " << TREF[ NSIM.NCLIM ] << endl;
                        //F2584             WRITE(15,937)mnem
                        outfile15 << "PROFILE: " << mnem << endl;
                        //F2585             WRITE(15,931)
                        outfile15 << "MID CLIMATE MODEL PARAMETERS" << endl;
                        //F2586             WRITE(15,934)
                        outfile15 << "YEAR       GHG     ESO21     ESO22     ESO23       ALL" << endl;
                        //F2587             WRITE(15,936)TREF(NCLIM)
                        outfile15 << "REF T " << TREF[ NSIM.NCLIM ] << endl;
                        //F2588           ENDIF
                    }
                    //F2589           WRITE(11,935)KYY,TGHG(NCLIM,KSG),TSO21(NCLIM,KSG), &
                    //F2590           TSO22(NCLIM,KSG),TSO23(NCLIM,KSG),TALL(NCLIM,KSG)
                    outfile11 << KYY << "," << TGHG[ NSIM.NCLIM ][ KSG ] << "," << TSO21[ NSIM.NCLIM ][ KSG ] << "," << TSO22[ NSIM.NCLIM ][ KSG ] << "," 
                    << TSO23[ NSIM.NCLIM ][ KSG ] << "," << TALL[ NSIM.NCLIM ][ KSG ] << endl;
                    //F2591           WRITE(15,935)KYY,TGHG(NCLIM,KSG),XSO21(NCLIM,KSG), &
                    //F2592           XSO22(NCLIM,KSG),XSO23(NCLIM,KSG),TALL(NCLIM,KSG)
                    outfile15 << KYY << "," << TGHG[ NSIM.NCLIM ][ KSG ] << "," << TSO21[ NSIM.NCLIM ][ KSG ] << "," << TSO22[ NSIM.NCLIM ][ KSG ] << "," 
                    << TSO23[ NSIM.NCLIM ][ KSG ] << "," << TALL[ NSIM.NCLIM ][ KSG ] << endl;
                    //F2593         ENDIF
                }
                //F2594 !
                //F2595         IF(NCLIM.EQ.3)THEN
                if( NSIM.NCLIM == 3 ) {
                    //F2596           IF(KSG.EQ.1)THEN
                    if( KSG == 1 ) {
                        //F2597             WRITE(12,937)mnem
                        outfile12 << "PROFILE: " << mnem << endl;
                        //F2598             WRITE(12,932)
                        outfile12 << "HIGH CLIMATE MODEL PARAMETERS" << endl;
                        //F2599             WRITE(12,934)
                        outfile12 << "YEAR       GHG     ESO21     ESO22     ESO23       ALL" << endl;
                        //F2600             WRITE(12,936)TREF(NCLIM)
                        outfile12 << "REF T " << TREF[ NSIM.NCLIM ] << endl;
                        //F2601             WRITE(16,937)mnem
                        outfile16 << "PROFILE: " << mnem << endl;
                        //F2602             WRITE(16,932)
                        outfile16 << "HIGH CLIMATE MODEL PARAMETERS" << endl;
                        //F2603             WRITE(16,934)
                        outfile16 << "YEAR       GHG     ESO21     ESO22     ESO23       ALL" << endl;
                        //F2604             WRITE(16,936)TREF(NCLIM)
                        outfile16 << "REF T " << TREF[ NSIM.NCLIM ] << endl;
                        //F2605           ENDIF
                    }
                    //F2606           WRITE(12,935)KYY,TGHG(NCLIM,KSG),TSO21(NCLIM,KSG), &
                    //F2607           TSO22(NCLIM,KSG),TSO23(NCLIM,KSG),TALL(NCLIM,KSG)
                    outfile12 << KYY << "," << TGHG[ NSIM.NCLIM ][ KSG ] << "," << TSO21[ NSIM.NCLIM ][ KSG ] << "," << TSO22[ NSIM.NCLIM ][ KSG ] << "," 
                    << TSO23[ NSIM.NCLIM ][ KSG ] << "," << TALL[ NSIM.NCLIM ][ KSG ] << endl;
                    //F2608           WRITE(16,935)KYY,TGHG(NCLIM,KSG),XSO21(NCLIM,KSG), &
                    //F2609           XSO22(NCLIM,KSG),XSO23(NCLIM,KSG),TALL(NCLIM,KSG)
                    outfile16 << KYY << "," << TGHG[ NSIM.NCLIM ][ KSG ] << "," << TSO21[ NSIM.NCLIM ][ KSG ] << "," << TSO22[ NSIM.NCLIM ][ KSG ] << "," 
                    << TSO23[ NSIM.NCLIM ][ KSG ] << "," << TALL[ NSIM.NCLIM ][ KSG ] << endl;
                    //F2610         ENDIF
                }
                //F2611 !
                //F2612         IF(NCLIM.EQ.4)THEN
                if( NSIM.NCLIM == 4 ) {
                    //F2613           IF(KSG.EQ.1)THEN
                    if( KSG == 1 ) {
                        //F2614             WRITE(13,937)mnem
                        outfile13 << "PROFILE: " << mnem << endl;
                        //F2615             WRITE(13,933)
                        outfile13 << "USER CLIMATE MODEL PARAMETERS" << endl;
                        //F2616             WRITE(13,934)
                        outfile13 << "YEAR       GHG     ESO21     ESO22     ESO23       ALL" << endl;
                        //F2617             WRITE(13,936)TREF(NCLIM)
                        outfile13 << "REF T " << TREF[ NSIM.NCLIM ] << endl;
                        //F2618             WRITE(17,937)mnem
                        outfile17 << "PROFILE: " << mnem << endl;
                        //F2619             WRITE(17,933)
                        outfile17 << "USER CLIMATE MODEL PARAMETERS" << endl;
                        //F2620             WRITE(17,934)
                        outfile17 << "YEAR       GHG     ESO21     ESO22     ESO23       ALL" << endl;
                        //F2621             WRITE(17,936)TREF(NCLIM)
                        outfile17 << "REF T " << TREF[ NSIM.NCLIM ] << endl;
                        //F2622           ENDIF
                    }
                    //F2623           WRITE(13,935)KYY,TGHG(NCLIM,KSG),TSO21(NCLIM,KSG), &
                    //F2624           TSO22(NSIM.NCLIM,KSG),TSO23(NSIM.NCLIM,KSG),TALL(NSIM.NCLIM,KSG)
                    outfile13 << KYY << "," << TGHG[ NSIM.NCLIM ][ KSG ] << "," << TSO21[ NSIM.NCLIM ][ KSG ] << "," << TSO22[ NSIM.NCLIM ][ KSG ] << "," 
                    << TSO23[ NSIM.NCLIM ][ KSG ] << "," << TALL[ NSIM.NCLIM ][ KSG ] << endl;
                    //F2625           WRITE(17,935)KYY,TGHG(NSIM.NCLIM,KSG),XSO21(NSIM.NCLIM,KSG), &
                    //F2626           XSO22(NSIM.NCLIM,KSG),XSO23(NSIM.NCLIM,KSG),TALL(NSIM.NCLIM,KSG)
                    outfile17 << KYY << "," << TGHG[ NSIM.NCLIM ][ KSG ] << "," << TSO21[ NSIM.NCLIM ][ KSG ] << "," << TSO22[ NSIM.NCLIM ][ KSG ] << "," 
                    << TSO23[ NSIM.NCLIM ][ KSG ] << "," << TALL[ NSIM.NCLIM ][ KSG ] << endl;
                    //F2627         ENDIF
                }
                //F2628 !
                //F2629         END DO
            }
            //F2630         END DO
        }
        //F2631 !
        //F2632         CLOSE(10)
        outfile10.close();
        //F2633         CLOSE(11)
        outfile11.close();
        //F2634         CLOSE(12)
        outfile12.close();
        //F2635         CLOSE(13)
        outfile13.close();
        //F2636         CLOSE(14)
        outfile14.close();
        //F2637         CLOSE(15)
        outfile15.close();
        //F2638         CLOSE(16)
        outfile16.close();
        //F2639         CLOSE(17)
        outfile17.close();
        //F2640 !
        //F2641       ENDIF
    }
    //F2642 !
    //F2643 !  **************************************************************
    //F2644 !
    //F2645 !  WRITE TEMPERATURES TO MAG DISPLAY FILE. NOTE THAT APPROPRIATE
    //F2646 !   TEMP (AND SEA LEVEL) DATA ARE SAVED AT THE RIGHT POINTS IN THE
    //F2647 !   NSIM LOOP, SO THESE DISPLAY FILES CAN BE PRODUCED *OUTSIDE*
    //F2648 !   THE NSIM LOOP.
    //F2649 !
    //F2650       open(unit=9,file='./outputs/temps.dis',status='UNKNOWN')
    ofstream outfile9;
    openfile_write( &outfile9, "./temps_c.dis", DEBUG_IO );
    //F2651 !
    //F2652         WRITE (9,213)
    outfile9 << "YEAR  TEMUSER    TEMLO   TEMMID    TEMHI TEMNOSO2" << endl;
    //F2653 !
    //F2654 !  PRINTOUT INTERVAL FOR DISPLAY PURPOSES SET BY IDIS IN MAGEXTRA.CFG
    //F2655 !
    //F2656         DO K=1,KEND,IDIS
    for( int K=1; K<=Limits.KEND; K+=IDIS ) {
        //F2657           IYEAR=1764+K
        int IYEAR = 1764 + K;
        //F2658           WRITE (9,226) IYEAR,TEMUSER(K),TEMLO(K),TEMMID(K), &
        //F2659           TEMHI(K),TEMNOSO2(K)
        outfile9 << IYEAR << "," << STOREDVALS.TEMUSER[ K ] << "," << TEMLO[ K ] << "," << TEMMID[ K ] << "," << TEMHI[ K ] << "," << TEMNOSO2[ K ] << endl;
        //F2660         END DO
    }
    //F2661 !
    //F2662         WRITE (9,213)
    outfile9 << "YEAR  TEMUSER    TEMLO   TEMMID    TEMHI TEMNOSO2" << endl;
    //F2663 !
    //F2664       CLOSE(9)
    outfile9.close();
    //F2665 !
    //F2666 !  **************************************************************
    //F2667 !
    //F2668 !  WRITE SEALEVEL CHANGES TO MAG DISPLAY FILE
    //F2669 !
    //F2670       open(unit=9,file='./outputs/sealev.dis',status='UNKNOWN')
    openfile_write( &outfile9, "./sealev_c.dis", DEBUG_IO );
    //F2671 !
    //F2672         WRITE (9,214)
    outfile9 << "YEAR  MSLUSER    MSLLO   MSLMID    MSLHI" << endl;
    //F2673 !
    //F2674 !  PRINTOUT INTERVAL FOR DISPLAY PURPOSES SET BY IDIS IN MAGEXTRA.CFG
    //F2675 !
    //F2676         DO K=1,KEND,IDIS
    for( int K=1; K<=Limits.KEND; K+=IDIS ) {
        //F2677           IYEAR=1764+K
        int IYEAR = 1764 + K;
        //F2678             WRITE (9,227) IYEAR,SLUSER(K),SLLO(K),SLMID(K), &
        //F2679             SLHI(K)
        outfile9 << IYEAR << "," << SLUSER[ K ] << "," << SLLO[ K ] << "," << SLMID[ K ] << "," << SLHI[ K ] << endl;
        //F2680         END DO
    }
    //F2681 !
    //F2682         WRITE (9,214)
    outfile9 << "YEAR  MSLUSER    MSLLO   MSLMID    MSLHI" << endl;
    //F2683 !
    //F2684       CLOSE(9)
    outfile9.close();
    //F2685 !
    //F2686 !  **************************************************************
    //F2687 !
    //F2688 !  WRITE EMISSIONS TO MAG DISPLAY FILE
    //F2689 !
    //F2690       open(unit=9,file='./outputs/emiss.dis',status='UNKNOWN')
    openfile_write( &outfile9, "./emiss_c.dis", DEBUG_IO );
    //F2691 !
    //F2692         WRITE (9,212)
    outfile9 << "YEAR  FOSSCO2 NETDEFOR      CH4      N2O SO2-REG1 SO2-REG2 SO2-REG3   SO2-GL" << endl;
    //F2693 !
    //F2694 !  PRINTOUT INTERVAL FOR DISPLAY PURPOSES SET BY IDIS IN MAGEXTRA.CFG
    //F2695 !   NOTE THAT ESO2(K) IS OVERWRITTEN IN LAST LOOP OF CLIMATE MODEL
    //F2696 !   SIMULATIONS, BUT ESO2i(K) REMAIN AS INPUTTED.
    //F2697 !  FOR ESO2 DISPLAY, NEED TO SUBTRACT THE 1990 TOTAL VALUE FROM EACH
    //F2698 !   REGION THEN ADD THE REGIONAL 1990 VALUE. THE REGIONAL VALUES FROM
    //F2699 !   BEFORE WERE 37,28,10 TgSO4, SUMMING TO 75 TgSO4.
    //F2700 !
    //F2701         DO K=1,225
    float EESS1[ iTp+1 ], EESS2[ iTp+1 ], EESS3[ iTp+1 ], EESST[ iTp+1 ];
    for( int K=1; K<=225; K++ ) {
        //F2702           EESS1(K)=0.0
        //F2703           EESS2(K)=0.0
        //F2704           EESS3(K)=0.0
        //F2705           EESST(K)=0.0
        EESS1[ K ] = EESS2[ K ] = EESS3[ K ] = EESST[ K ] = 0.0;
        //F2706         END DO
    }
    //F2707 ! 
    //F2708         DO K=226,KEND,IDIS
    for( int K=226; K<=Limits.KEND; K+=IDIS ) {
        //F2709           IYEAR=1764+K
        int IYEAR = 1764 + K;
        //F2710           EESS1(K)=ESO21(K)-ES1990*(1.0-37.0/75.0)
        EESS1[ K ] = CONCS.ESO21.getval( K ) - Sulph.ES1990 * ( 1.0 - 37.0 / 75.0 );
        //F2711           EESS2(K)=ESO22(K)-ES1990*(1.0-28.0/75.0)
        EESS2[ K ] = CONCS.ESO22.getval( K ) - Sulph.ES1990 * ( 1.0 - 37.0 / 75.0 );
        //F2712           EESS3(K)=ESO23(K)-ES1990*(1.0-10.0/75.0)
        EESS3[ K ] = CONCS.ESO23.getval( K ) - Sulph.ES1990 * ( 1.0 - 37.0 / 75.0 );
        //F2713           EESST(K)=EESS1(K)+EESS2(K)+EESS3(K)
        EESST[ K ] = EESS1[ K ] + EESS2[ K ] + EESS3[ K ];
        //F2714           WRITE (9,225) IYEAR,EF(K),EDNET(K),ECH4(K),EN2O(K), &
        //F2715            EESS1(K),EESS2(K),EESS3(K),EESST(K)
        outfile9 << IYEAR << "," << CARB.EF.getval( K ) << "," << METH1.ednet.getval( K ) << "," << CONCS.ECH4.getval( K ) << "," << CONCS.EN2O.getval( K ) << "," 
        << EESS1[ K ] << "," << EESS2[ K ] << "," << EESS3[ K ] << "," << EESST[ K ] << endl;
        //F2716         END DO
    }
    //F2717 !
    //F2718         WRITE (9,212)
    outfile9 << "YEAR  FOSSCO2 NETDEFOR      CH4      N2O SO2-REG1 SO2-REG2 SO2-REG3   SO2-GL" << endl;
    //F2719 !
    //F2720       CLOSE(9)
    outfile9.close();
    //F2721 !
    //F2722 !  ************************************************************
    //F2723 !
    //F2724 !  OPEN NEW OUTPUT FILE (FRACLEFT.OUT).
    //F2725 !
    //F2726       OPEN(UNIT=888,file='./outputs/FRACLEFT.OUT',STATUS='UNKNOWN')
    ofstream outfile888;
    openfile_write( &outfile888, "./fracleft_c.out", DEBUG_IO );
    //F2727 !
    //F2728 !  FRACTION OF CO2 REMAINING IN ATMOSPHERE
    //F2729 !
    //F2730       WRITE(888,887)
    outfile888 << "  YEAR    ATMASS    CUMEMS  FRACLEFT" << endl;
    //F2731       EMTOT=300.0
    float EMTOT = 300.0;
    //F2732 !
    //F2733       DO K=226,KEND
    for( int K=226; K<=Limits.KEND; K++) {
        //F2734         EMTOT=EMTOT+EF(K)+EDNET(K)
        EMTOT += CARB.EF.getval( K ) + METH1.ednet.getval( K );
        //F2735         ATBIT=2.123*(CO2(K)-278.)
        float ATBIT = 2.123 * ( CARB.CO2[ K ] - 278.0 );
        //F2736         FRACLEFT=ATBIT/EMTOT
        float FRACLEFT = ATBIT / EMTOT;
        //F2737         IIYY=K+1764
        int IIYY = K + 1764;
        //F2738         WRITE(888,889)IIYY,ATBIT,EMTOT,FRACLEFT
        outfile888 << IIYY << "," << ATBIT << "," << EMTOT << "," << FRACLEFT << endl;
        //F2739       END DO
    }
    outfile888.close();
    //F2740  887  FORMAT(/1X,'  YEAR    ATMASS    CUMEMS  FRACLEFT')
    //F2741  889  FORMAT(1X,I6,2F10.3,F10.5)
    //F2742 !
    //F2743 !  **************************************************************
    //F2744 !
    //F2745 !  WRITE DATA TO CCSM FILE
    //F2746 !
    //F2747       IF(ICCSM.EQ.1)THEN
    if( ICCSM == 1 ) {
        //F2748         WRITE(88,883)
        outfile88 << "*** MIDYEAR CONCENTRATIONS : ESO2, BC & OC SET TO ZERO BEFORE 1990 ***" << endl;
        //F2749         WRITE(88,882)
        outfile88 << "YEAR       CO2       CH4       N2O     CFC12    C11EFF   QTROPOZ  QSTRATOZ      ESO2        BC        OC     YEAR" << endl;
        //F2750         DO K=1,IYEND-1764
        for( int K=1; K<=IYEND-1764; K++ ) {
            //F2751           KKYR=K+1764
            int KKYR = K + 1764;
            //F2752 !
            //F2753 !  CONVERT END OF YEAR TO MIDYEAR CONCS
            //F2754 !
            //F2755           CO2MID =(CO2(K)+CO2(K-1))/2.
            float CO2MID = ( CARB.CO2[ K ] + CARB.CO2[ K-1 ] ) / 2.0;
            //F2756           CH4MID =(CH4(K)+CH4(K-1))/2.
            float CH4MID = ( CONCS.CH4[ K ] + CONCS.CH4[ K-1 ] ) / 2.0;
            //F2757           CN2OMID=(CN2O(K)+CN2O(K-1))/2.
            float CN2OMID = ( CONCS.CN2O[ K ] + CONCS.CN2O[ K-1 ] ) / 2.0;
            //F2758           IF(K.LE.225)THEN
            float XBC, XOC;
            if( K <= 225 )
                //F2759             XBC=0.0
                //F2760             XOC=0.0
                XBC = XOC = 0.0;
            //F2761           ELSE
            else {
                //F2762             XBC=EBC(K)
                XBC = CONCS.EBC.getval( K );
                //F2763             XOC=EOC(K)
                XOC = CONCS.EOC.getval( K );
                //F2764           ENDIF
            }
            //F2765           WRITE(88,881)KKYR,CO2MID,CH4MID,CN2OMID,CFC12(K),C11EFF(K), &
            //F2766            QTROZ(K),QSTROZ(K),EESST(K),XBC,XOC,KKYR
            outfile88 << KKYR << "," << CO2MID << "," << CH4MID << "," << CN2OMID << "," 
            << FORCE.CFC12[ K ] << "," << C11EFF[ K ] << "," << QTROZ[ K ] << "," << QSTROZ[ K ] << "," 
            << EESST[ K ] << "," << XBC << "," << XOC << "," << KKYR << endl;
            //F2767         END DO
        }
        //F2768       ENDIF
    }
    outfile88.close();
    //F2769 !
    //F2770  881  FORMAT(1X,I5,5F10.3,2F10.4,3F10.3,I9)
    //F2771  882  FORMAT(/2X,'YEAR       CO2       CH4       N2O     CFC12', &
    //F2772       '    C11EFF   QTROPOZ  QSTRATOZ      ESO2', &
    //F2773       '        BC        OC     YEAR')
    //F2774  883  FORMAT(/1X,'*** MIDYEAR CONCENTRATIONS : ESO2, BC & OC SET', &
    //F2775       ' TO ZERO BEFORE 1990 ***')
    //F2776 !                  
    //F2777 !  **************************************************************
    //F2778 !  **************************************************************
    //F2779 !
    //F2780 !  FORMAT STATEMENTS
    //F2781 !
    //F2782  10   FORMAT (/1X,'CO2-DOUBLING FORCING IN W/M**2 =',F6.3)
    //F2783  11   FORMAT (1X,'FNHOC= ',F4.2,' * FSHOC= ',F4.2,' * FNHLAND= ',F4.2, &
    //F2784       ' * FSHLAND= ',F4.2)
    //F2785  110  FORMAT (1X,A4,' CONCENTRATION PROJECTION FOR CO2')
    //F2786  1100 FORMAT (3X,'(CO2-CLIMATE FEEDBACK NOT INCLUDED)')
    //F2787  1101 FORMAT (3X,'(CO2-CLIMATE FEEDBACK INCLUDED)')
    //F2788  111  FORMAT (1X,'Dn(1980s) =',F6.3,' : Foc(1980s) =',F6.3)
    //F2789  112  FORMAT (1X,A4,' CONCENTRATION PROJECTION FOR CH4')
    //F2790  113  FORMAT (1X,'CH4 CONCS USE CONSTANT LIFETIME OF',F7.3,'YEARS')
    //F2791  114  FORMAT (1X,A4,' 1990 FORCINGS FOR SO4 AEROSOL')
    //F2792  115  FORMAT (1X,'LEVCO2, LEVCH4 AND/OR LEVSO4 WRONGLY SET > 4 :', &
    //F2793       ' RESET AT 2')
    //F2794  117  FORMAT (2X,'STRAT OZONE DEPLETION FEEDBACK OMITTED')
    //F2795  1171 FORMAT (2X,'STRAT OZONE DEPLETION FEEDBACK INCLUDED')
    //F2796  118  FORMAT (2X,'FOR HALOCARBONS NOT IN GAS.EMK,', &
    //F2797       ' EMS DROP TO ZERO OVER 2100-2200')
    //F2798  1181 FORMAT (2X,'FOR HALOCARBONS NOT IN GAS.EMK,', &
    //F2799       ' EMS CONSTANT AFTER 2100')
    //F2800  116  FORMAT (/1X,'CLIMATE MODEL SELECTED = ',A7)
    //F2801  1161 FORMAT (1X,'USER ICE MELT = LOW')
    //F2802  1162 FORMAT (1X,'USER ICE MELT = MID')
    //F2803  1163 FORMAT (1X,'USER ICE MELT = HIGH')
    //F2804  1164 FORMAT (1X,'TAR GSIC SENSITIVITY =',F8.4,' CM/YR-DEGC : VZERO =', &
    //F2805       F5.1,'CM')
    //F2806  12   FORMAT (1X,'XKNS=',F4.1,' : XKLO=',F4.1)
    //F2807  120  FORMAT (1X,'HM=',F5.1,'M : XK=',F6.4,'CM**2/SEC')
    //F2808  121  FORMAT (1X,'PI=',F6.4,' : INITIAL W=',F5.2,'M/YR',/)
    //F2809  122  FORMAT (1X,'CONSTANT W CASE')
    //F2810  1220 FORMAT (1X,'IVARW SET AT',I2)
    //F2811  123  FORMAT (1X,'VARIABLE W : NH W = ZERO WHEN TEMPERATURE =',F6.2, &
    //F2812       'degC')
    //F2813  1231 FORMAT (1X,'FULL W SCALED WITH GLOBAL-MEAN TEMPERATURE',/)
    //F2814  1232 FORMAT (1X,'FULL W SCALED WITH GLOBAL-MEAN OCEAN TEMPERATURE',/)
    //F2815  1233 FORMAT (1X,'FULL W SCALED WITH HEMISPHERIC-MEAN OCEAN', &
    //F2816       'TEMPERATURE',/)
    //F2817  1234 FORMAT (1X,'ACTIVE W SCALED WITH GLOBAL-MEAN TEMPERATURE',/)
    //F2818  1235 FORMAT (1X,'ACTIVE W SCALED WITH GLOBAL-MEAN OCEAN TEMPERATURE',/)
    //F2819  124  FORMAT (1X,'VARIABLE W : SH W = ZERO WHEN TEMPERATURE =',F6.2, &
    //F2820       'degC')
    //F2821  125  FORMAT (1X,'VARIABLE W : NH AND SH W(t) SPECIFIED IN WINPUT.IN')
    //F2822  126  FORMAT (1X,'PERMANENT THC SHUTDOWN AT W =',F6.2,'M/YR')
    //F2823  127  FORMAT (1X,'W = ZERO WHEN TEMPERATURE =',F6.2, &
    //F2824       'degC')
    //F2825  140  FORMAT (/1X,'1880-1990 CHANGES : GLOBAL DTEMP =',F7.3, &
    //F2826       ' :   DMSL =',F7.3)
    //F2827  141  FORMAT (1X,'          DTNHL =',F7.3,' : DTNHO =',F7.3, &
    //F2828       ' :  DTSHL =',f7.3,' :   DTSHO =',f7.3)
    //F2829  142  FORMAT (1X,'           DTNH =',F7.3,' :  DTSH =',F7.3, &
    //F2830       ' : DTLAND =',f7.3,' : DTOCEAN =',f7.3)
    //F2831  15   FORMAT (/1X,'** TEMPERATURE AND SEA LEVEL CHANGES FROM',I5,' **')
    //F2832  16   FORMAT (1X,'     (FIRST LINE GIVES 1765-1990 CHANGES : ', &
    //F2833                        'ALL VALUES ARE MID-YEAR TO MID-YEAR)')
    //F2834  161  FORMAT (/1X,'LOW CLIMATE AND SEA LEVEL MODEL PARAMETERS')
    //F2835  162  FORMAT (/1X,'MID CLIMATE AND SEA LEVEL MODEL PARAMETERS')
    //F2836  163  FORMAT (/1X,'HIGH CLIMATE AND SEA LEVEL MODEL PARAMETERS')
    //F2837  164  FORMAT (/1X,'USER CLIMATE AND SEA LEVEL MODEL PARAMETERS')
    //F2838  171  FORMAT (1X,' YEAR,DELTAQ, TEQU, TEMP, EXPN, GLAC,', &
    //F2839       ' GREENL,ANTAR,Z-XTRA,MSLTOT, TNH, TSH, WNH, WSH,YEAR,', &
    //F2840       '  GR+ANT,ZTOT-ZXTRA')
    //F2841  172  FORMAT (1X,' YEAR,DELTAQ, TEQU, TEMP, EXPN, GLAC,', &
    //F2842       ' GREENL,ANTAR, Z-XTRA, MSLTOT, TLAND,  TOCN, TL/TO,   WNH,   WSH,', &
    //F2843       '  YEAR')
    //F2844  173  FORMAT (1X,' YEAR,DELTAQ, TEQU, TEMP, EXPN, GLAC,', &
    //F2845       ' GREENL,  ANTAR, Z-XTRA, MSLTOT, TEQ-T, TDEEP,   WNH,   WSH,  YEAR,')
    //F2846  174  FORMAT (1X,' YEAR,EQVCO2, TEQU, TEMP, EXPN, GLAC,', &
    //F2847       ' GREENL,ANTAR.Z-XTRA.MSLTOT.TEQ-T.TDEEP, WNH, WSH,YEAR,')
    //F2848  175  FORMAT (1X,' YEAR,DELTAQ,   TEMP,TL/TO,MSLTOT, EXPN, GLAC,', &
    //F2849       ' GREENL,ANTAR Z-XTRA, WNH,YEAR')
    //F2850  176  FORMAT (/1X,'NSIM =',I3,' : DELT(2XCO2) =',F6.3,'DEGC')
    //F2851 1761  FORMAT (/1X,' DELT(2XCO2) =',F6.3,'DEGC')
    //F2852  177  FORMAT (/1X,'DT2X =',F5.2,' : CONSTANT W')
    //F2853  178  FORMAT (/1X,'DT2X =',F5.2,' : VARIABLE W')
    //F2854  179  FORMAT (/1X,'*************************************************')
    //F2855  181  FORMAT ('TO1990, ',e18.10,20(',',e18.10))
    //F2856  182  FORMAT ('TO1990, ',e18.10,20(',',e18.10))
    //F2857  183  FORMAT ('TO1990, ',e18.10,20(',',e18.10))
    //F2858  184  FORMAT ('TO1990, ',e18.10,20(',',e18.10))
    //F2859  185  FORMAT ('TO1990, ',e18.10,20(',',e18.10))
    //F2860  186  FORMAT (1X,'FULL GLOBAL SO2 EMISSIONS',/)
    //F2861  187  FORMAT (1X,'SO2 EMISSIONS CONSTANT AFTER 1990',/)
    //F2862  188  FORMAT (1X,'REGION 1 SO2 EMISSIONS',/)
    //F2863  189  FORMAT (1X,'REGION 2 SO2 EMISSIONS',/)
    //F2864  190  FORMAT (1X,'REGION 3 SO2 EMISSIONS',/)
    //F2865  191  FORMAT (1X,I5,',',e18.10,12(',',e18.10),',',I6,2(',',e18.10))
    //F2866  192  FORMAT (1X,I5,',',e18.10,12(',',e18.10),',',I6)
    //F2867  193  FORMAT (1X,I5,',',e18.10,11(',',e18.10),',',I6)
    //F2868  194  FORMAT (1X,I5,',',e18.10,11(',',e18.10),',',I6)
    //F2869  195  FORMAT (1X,I5,',',e18.10,8(',',e18.10),',',I6)
    //F2870 !
    //F2871  20   FORMAT (1X,'*** CONCENTRATIONS (CO2,PPM : CH4,N2O,PPB) ***', &
    //F2872       /1X,'*** MIDYEAR VALUES ***')
    //F2873  202  FORMAT (1X,'*** CONCENTRATIONS (CO2,PPM : CH4,N2O,PPB) ***', &
    //F2874       /1X,'*** START OF YEAR VALUES FOR YR.GE.1990 IN COLS 2,3,4 ***')
    //F2875  203  FORMAT (1X,'*** CONCENTRATIONS (CO2,PPM : CH4,N2O,PPB) ***', &
    //F2876       /1X,'*** END OF YEAR VALUES FOR YR.GE.1990 IN COLS 2,3,4 ***')
    //F2877  201  FORMAT (5X,'<- USER MODEL CONCS ->', &
    //F2878       '<------ CH4 & CO2 MID CONCS & RANGES ------->')
    //F2879  21   FORMAT (1X,'YEAR,EFOSS,NETDEF,CH4,N2O,NOX,VOC,CO,SO2REG1,SO2REG2,SO2REG3,',&
    //F2880                  'CF4,C2F6,HFC125,HFC134A,HFC143A,HFC227ea,HFC245ca,SF6,ESO2TOT,YEAR')
    //F2881  210  FORMAT (1X,'YEAR      CO2     CH4    N2O', &
    //F2882       '   CH4LO  CH4MID   CH4HI   CO2LO  CO2MID   CO2HI  YEAR', &
    //F2883       ' TAUCH4')
    //F2884  211  FORMAT (1X,'YEAR CO2USER   CO2LO  CO2MID   CO2HI', &
    //F2885                  ' CH4USER   CH4LO  CH4MID   CH4HI     N2O', &
    //F2886       ' MIDTAUCH4')
    //F2887  212  FORMAT (1X,'YEAR  FOSSCO2 NETDEFOR      CH4      N2O', &
    //F2888       ' SO2-REG1 SO2-REG2 SO2-REG3   SO2-GL')
    //F2889  213  FORMAT (1X,'YEAR  TEMUSER    TEMLO   TEMMID    TEMHI TEMNOSO2')
    //F2890  214  FORMAT (1X,'YEAR  MSLUSER    MSLLO   MSLMID    MSLHI')
    //F2891  220  FORMAT (1X,I4,',',e18.10,',',e18.10,7(',',e18.10),',',I6,',',e18.10)
    //F2892  221  FORMAT (1X,I4,',',e18.10,',',e18.10,',',e18.10,',',',,,,,,',I6)
    //F2893  222  FORMAT (1X,I4,',',21(e18.10,','),I6)
    //F2894  223  FORMAT (1X,I4,9F8.1)
    //F2895  224  FORMAT (1X,I4,9F8.1,F10.2)
    //F2896  225  FORMAT (1X,I4,8F9.2)
    //F2897  226  FORMAT (1X,I4,5F9.3)
    //F2898  227  FORMAT (1X,I4,5F9.1)
    //F2899  23   FORMAT (1X,'** INPUT EMISSIONS **')
    //F2900  231  FORMAT (4X,'BALANCED EMISSIONS FOR CH4 & N2O : SO2 EMISSIONS', &
    //F2901       ' RELATIVE TO 1990')
    //F2902  24   FORMAT (1X,'** CARBON CYCLE DETAILS **')
    //F2903  241  FORMAT (1X,'CONCENTRATIONS ARE UNCORRECTED MODEL OUTPUT', &
    //F2904       ' : LEVCO2 =',I2)
    //F2905 ! 25   FORMAT (1X,'FEEDBACKS ** TEMPERATURE * GPP :',F7.4,' * RESP :'
    //F2906 !     +,F7.4,' * LITT OXDN :',F7.4,'  **  FERTIL :',F7.4)
    //F2907  28   FORMAT (1X,F8.1,7F8.2,4f8.1)
    //F2908 !
    //F2909  30   FORMAT (1X,'  ')
    //F2910  31   FORMAT (1X,'****************************************************')
    //F2911  47   FORMAT (1X,'** DECADAL CONTRIBUTIONS TO', &
    //F2912        ' GLOBAL RADIATIVE FORCING **')
    //F2913  48   FORMAT (1X,'   (DELTA-Q in W/m**2 : PERCENTAGES IN BRACKETS : ', &
    //F2914        'BASED ON END-OF-YEAR FORCING VALUES)')
    //F2915 !
    //F2916  50   FORMAT (1X,'  INTERVAL')
    //F2917  53   FORMAT (1X,'STRAT H2O FROM CH4 : DELQH2O/DELQCH4 =',F6.3)
    //F2918  55   FORMAT (1X,'** GAS BY GAS DELTA-Q FROM',I5' : MIDYEAR VALUES **')
    //F2919  56   FORMAT (1X,'(BASED ON MID-YEAR FORCING VALUES : QEXTRA', &
    //F2920       ' NOT INCLUDED)')
    //F2921  561  FORMAT (1X,'CH4tot INCLUDES STRATH2O : TROPO3 INCLUDES CH4', &
    //F2922       ' COMPONENT')
    //F2923 5611  FORMAT (1X,'QAERMN IS THE SUM OF NITRATE AND MINERAL DUST', &
    //F2924       ' AEROSOL FORCING')
    //F2925  562  FORMAT (1X,'HALOtot INCLUDES STRAT O3')
    //F2926  563  FORMAT (1X,'HALOtot (AND QTOTAL) DOES NOT INCLUDE STRAT O3')
    //F2927  57   FORMAT (1X,'YEAR,CO2,CH4tot,N2O, HALOtot,','TROPOZ,SO4DIR,SO4IND,BIOAER,FOC+FBC,QAERMN,QLAND,',&
    //F2928  	   ' TOTAL, YEAR,CH4-O3,',&
    //F2929  		' STRATO3, MONTDIR,QKYOTO')
    //F2930  573   FORMAT (1X,'YEAR,CO2,CH4tot,N2O, HALOtot,','TROPOZ,SO4DIR,SO4IND,BIOAER,FOC+FBC,QAERMN,QLAND,',&
    //F2931  	   ' TOTAL, YEAR,CH4-O3,',&
    //F2932  		' STRATO3, MONTDIR,QKYOTO,BC,OC,QEXTRA')
    //F2933  571  FORMAT (1X,I4,12(',',e18.10),',',I4,10(',',e18.10))
    //F2934  58   FORMAT (1X,'** GAS BY GAS DELTA-Q FROM 1765 : MIDYEAR VALUES **')
    //F2935 !
    //F2936  60   FORMAT (1X,'1990 DIRECT AEROSOL FORCING          =',F6.3,'W/m**2')
    //F2937  61   FORMAT (1X,'1990 INDIRECT AEROSOL FORCING        =',F6.3,'W/m**2')
    //F2938  62   FORMAT (1X,'1990 BIOMASS AEROSOL FORCING         =',F6.3,'W/m**2')
    //F2939  63   FORMAT (1X,'1990 FOSSIL ORG C + BLACK C FORCING  =',F6.3,'W/m**2')
    //F2940 !
    //F2941  756  FORMAT (/1X,'NO EXTRA FORCING ADDED')
    //F2942  757  FORMAT (/1X,'EXTRA GLOBAL MEAN FORCING ADDED FROM', &
    //F2943       ' QEXTRA.IN OVER',I5,' TO',I5,' INCLUSIVE')
    //F2944  758  FORMAT (/1X,'EXTRA HEMISPHERIC FORCINGS ADDED FROM', &
    //F2945       ' QEXTRA.IN OVER',I5,' TO',I5,' INCLUSIVE')
    //F2946  759  FORMAT (/1X,'EXTRA NHO, NHL, SHO, SHL FORCINGS ADDED FROM', &
    //F2947       ' QEXTRA.IN OVER',I5,' TO',I5,' INCLUSIVE')
    //F2948  760  FORMAT (2X,F10.3,' W/m**2 SUBTRACTED FROM ALL VALUES')
    //F2949  761  FORMAT (2X,'FORCING SCALED BY',F7.3,' AFTER OFFSET')
    //F2950  762  FORMAT (/1X,'QEXTRA FORCING USED ALONE')
    //F2951 !
    //F2952  800  FORMAT (1X,'LOW CONC CASE  : NETDEF(80s) = 1.80GtC/yr', &
    //F2953       ' : GIFFORD FERTILIZATION FACTOR =',F6.3)
    //F2954  801  FORMAT (1X,'MID CONC CASE  : NETDEF(80s) = 1.10GtC/yr', &
    //F2955       ' : GIFFORD FERTILIZATION FACTOR =',F6.3)
    //F2956  802  FORMAT (1X,'HIGH CONC CASE : NETDEF(80s) = 0.40GtC/yr', &
    //F2957       ' : GIFFORD FERTILIZATION FACTOR =',F6.3)
    //F2958  803  FORMAT (1X,'USER CONC CASE : NETDEF(80s) =',F5.2, &
    //F2959       'GtC/yr : GIFFORD FERTILIZATION FACTOR =',F6.3)
    //F2960  804  FORMAT (/1X,'ALL CASES USE 1980s MEAN OCEAN FLUX OF 2.0GtC/yr')
    //F2961  805  FORMAT (1X,'DETAILED CARBON CYCLE OUTPUT IS FOR LEVCO2 CASE ONLY')
    //F2962  806  FORMAT (1X,'METHANE OXIDATION TERM INCLUDED IN EMISSIONS')
    //F2963  807  FORMAT (1X,'METHANE OXIDATION TERM NOT INCLUDED IN EMISSIONS')
    //F2964  8071 FORMAT (1X,'NOTE: CORRECTION TO MATCH OBSERVED IN 2000 NOT', &
    //F2965       ' APPLIED IN THIS SECTION')
    //F2966  808  FORMAT (/1X,'*** ERROR : D80SIN SET TOO LOW AT',F7.3, &
    //F2967       ' : RESET AT'F7.3,' ***')
    //F2968  809  FORMAT (/1X,'*** ERROR : D80SIN SET TOO HIGH AT',F7.3, &
    //F2969       ' : RESET AT'F7.3,' ***')
    //F2970  810  FORMAT (/77X,'ENDYEAR')
    //F2971  811  FORMAT (/77X,'MIDYEAR')
    //F2972  812  FORMAT (1X,'YEAR, ETOTAL,  EFOSS, CH4OXN,   NETD, GROSSD,  OFLUX,',' ABFRAC, PLANT C, HLITT,    SOIL,',&
    //F2973  '    CONC,  DEL-M,  YEAR') 
    //F2974  813  FORMAT (1X,I4,',',e18.10,11(',',e18.10),',',I6)
    //F2975 !
    //F2976  871  FORMAT (/1X,'CO2 CONC INPUT OVERWRITES CO2 EMISSIONS :', &
    //F2977       ' POST-1990 CO2 FORCING SCALED UP BY',f5.1,'%')
    //F2978  872  FORMAT (/1X,'CO2 CONC INPUT OVERWRITES CO2 EMISSIONS :', &
    //F2979       ' OTHER GAS EMISSIONS AS SPECIFIED IN GAS.EMK')
    //F2980  873  FORMAT (/1X,'CO2 CONC INPUT OVERWRITES CO2 EMISSIONS :', &
    //F2981       ' SO2 EMISSIONS AS SPECIFIED IN GAS.EMK')
    //F2982 !
    //F2983  900  FORMAT (1X,I5)
    //F2984  901  FORMAT (1X,2I5)
    //F2985  902  FORMAT (1X,I5,F10.0)
    //F2986  903  FORMAT (1X,I5,2F10.0)
    //F2987  904  FORMAT (1X,I5,4F10.0)
    //F2988 !
    //F2989  914  FORMAT (/1X,'DIFF/L SENSITIVITY CASE : RLO =',F6.3,' : XLAML =', &
    //F2990       F10.4,' : XLAMO =',F10.4)
    //F2991  915  FORMAT (/1X,'GLOBAL SENSITIVITY CASE : INITIAL XLAM =',F10.4)
    //F2992  916  FORMAT (1X,'  **  WARNING, XLAML<0.0 : USE SMALLER XKLO **')
    //F2993  930  FORMAT (/1X,'LOW CLIMATE MODEL PARAMETERS')
    //F2994  931  FORMAT (/1X,'MID CLIMATE MODEL PARAMETERS')
    //F2995  932  FORMAT (/1X,'HIGH CLIMATE MODEL PARAMETERS')
    //F2996  933  FORMAT (/1X,'USER CLIMATE MODEL PARAMETERS')
    //F2997  934  FORMAT (/2X,'YEAR       GHG     ESO21     ESO22     ESO23', &
    //F2998       '       ALL')
    //F2999  935  FORMAT (1X,I5,5F10.4)
    //F3000  936  FORMAT (1X,'REF T',40X,F10.4)
    //F3001  937  format (1x,'PROFILE: ',A20)
    //F3002 !
    //F3003  4240 FORMAT (I10)
    //F3004  4241 FORMAT (F10.0)
    //F3005  4242 FORMAT (1X,I5,20F10.0)
    //F3006  4243 FORMAT (I2)
    //F3007  4445 FORMAT(1X,I5,2F10.0)
    //F3008  4446 FORMAT(1X,2I5)
    //F3009  4447 FORMAT(1X,I5,4F10.0,20X,F10.0)
    //F3010  4448 FORMAT(1X,I5,2F10.0,20X,3F10.0)
    //F3011 !
    //F3012       END
    //F3013 !
    //F3014 !********************************************************************
    //F3015 !
    //F3016       BLOCK DATA
    // Assignments corresponding to these DATA statements at top of climat()
    //F3017 !
    //F3018       parameter (iTp=740)
    //F3019 !
    //F3020       COMMON/CLIM/IC,IP,KC,DT,DZ,FK,HM,Q2X,QXX,PI,T,TE,TEND,W0,XK,XKLO, &
    //F3021       XKNS,XLAM,FL(2),FO(2),FLSUM,FOSUM,HEM(2),P(40),TEM(40),TO(2,40), &
    //F3022       AL,BL,CL,DTH,DTZ,DZ1,XLL,WWW,XXX,YYY,RHO,SPECHT,HTCONS,Y(4)
    //F3023 !
    //F3024       COMMON/COBS/COBS(0:236)
    //F3025 !
    //F3026       COMMON/CARB/CCO2(4,224:iTp),EDGROSS(4,226:iTp),EF(226:iTp+1), &
    //F3027       REGROW(4,226:iTp),PL(4,226:iTp),HL(4,226:iTp),SOIL(4,226:iTp), &
    //F3028       TTT(226:iTp),ESUM(226:iTp),ETOT(4,226:iTp),EDNET90(4), &
    //F3029       FOC(4,226:iTp),CO2(0:iTp),CO2SAVE(0:iTp)
    //F3030 !
    //F3031       COMMON/CAR/EL1,EL2,EL3,TINV0(5),TINV(4,5),A(3,5),AA(4,5), &
    //F3032       BCO2(4),BTGPP,BTRESP,BTHUM,GAMP,GPP0,RESP0,QA0,U0,C0,B340(4), &
    //F3033       PHI,RG,TAUP,TAUH,TAUS,THP,THS,THH0,THS0,THPL,G1,G2,G3,FACTOR, &
    //F3034       EL21,EL32,XX1,XX2,XX3,XX4,XX5,XX6,DEE1,DEE2,DEE3,DEE4,DEE5,DEE6, &
    //F3035       FL1,FL2,FL3,XL,GAMH,GAMS,QS0,BTSOIL,FERTTYPE,TOTEM,CONVTERP, &
    //F3036       R(4),CPART(4,5),DELMASS(4,226:iTp),ABFRAC(4,226:iTp)
    //F3037 !
    //F3038       DATA FL(1)/0.420/,FL(2)/0.210
    //F3039 !
    //F3040       DATA RHO/1.026/,SPECHT/0.9333/,HTCONS/4.1856/
    //F3041 !
    //F3042 !  INITIALISE CARBON CYCLE MODEL PARAMETERS.
    //F3043 !  FIRST SPECIFY CUMULATIVE EMISSIONS TRANSITION POINTS.
    //F3044 !
    //F3045       DATA EL1/141./,EL2/565./,EL3/2500./
    //F3046 !
    //F3047       DATA DEE1/0.25/,DEE2/0.5/,DEE3/1.0/,DEE4/2.0/ &
    //F3048       ,DEE5/4.0/,DEE6/8.0/
    //F3049 !
    //F3050 !  THESE ARE THE INVERSE DECAY TIMES AND MULTIPLYING CONSTANTS
    //F3051 !
    //F3052       DATA (TINV0(J),J=1,5)/0.0,0.0030303,0.0125,0.05,0.625/
    //F3053       DATA (A(1,J),J=1,5)/0.131,0.216,0.261,0.294,0.098/
    //F3054       DATA (A(2,J),J=1,5)/0.142,0.230,0.335,0.198,0.095/
    //F3055       DATA (A(3,J),J=1,5)/0.166,0.363,0.304,0.088,0.079/
    //F3056 !
    //F3057         end
    outfile8.close();

    setGlobals( &CARB, &TANDSL, &CONCS, &NEWCONCS, 
              &STOREDVALS, &NEWPARAMS, &BCOC, 
              &METH1, &CAR, &FORCE, &JSTART,
              &QADD, &HALOF );
}
