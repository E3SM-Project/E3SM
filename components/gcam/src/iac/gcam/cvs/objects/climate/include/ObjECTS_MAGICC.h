#ifndef _ObjECTS_MAGICC_h_
#define _ObjECTS_MAGICC_h_

/*
 *  ObjECTS_MAGICC.f90.h
 *  magicc++
 *
 *  Created by d3x290-local on 10/7/09.
 *  Copyright 2009 DOE Pacific Northwest Lab. All rights reserved.
 *
 */

// iTp is used extensively in array declarations, so it's special
#define iTp 740

#include "climate/include/MAGICC_array.h"

//#define DEBUG_MAGICC++
#define __func__ __FUNCTION__


// Shared structures

/*typedef*/ struct AREAS_block {
    //F 200       COMMON /AREAS/FNO,FNL,FSO,FSL
    float   FNO,    /* Area fraction, northern hemisphere ocean */
            FNL,    /* Area fraction, northern hemisphere land  */ 
            FSO,    /* Area fraction, southern hemisphere ocean */
            FSL;    /* Area fraction, southern hemisphere land  */
} ;

typedef struct {
    //F       COMMON/BCOC/FBC1990, FOC1990, FSO2_dir1990,FSO2_ind1990, aBCUnitForcing, aOCUnitForcing, &
    float FBC1990, FSO2_dir1990, FSO2_ind1990, aBCUnitForcing, aOCUnitForcing;
    //F              aBCBaseEmissions, aOCBaseEmissions !sjs
    float aBCBaseEmissions, aOCBaseEmissions;
} BCOC_block;

typedef struct {
    //F       COMMON/CAR/EL1,EL2,EL3,TINV0(5),TINV(4,5),A(3,5),AA(4,5), &
    float EL1, EL2, EL3/* , TINV0[ 5+1 ], TINV[ 4+1 ][ 5+1 ], A[ 3+1 ][ 5+1 ], AA[ 4+1 ][ 5+1 ] */;
    magicc_array TINV0, TINV, A, AA;
    //F       BCO2(4),BTGPP,BTRESP,BTHUM,GAMP,GPP0,RESP0,QA0,U0,C0,B340(4), &
    float BCO2[ 4+1 ], BTGPP, BTRESP, BTHUM, GAMP, GPP0, RESP0, QA0, U0, C0, B340[ 4+1 ];
    //F       PHI,RG,TAUP,TAUH,TAUS,THP,THS,THH0,THS0,THPL,G1,G2,G3,FACTOR, &
    float PHI, RG, TAUP, TAUH, TAUS, THP, THS, THH0, THS0, THPL, G1, G2, G3, FACTOR;
    //F       EL21,EL32,XX1,XX2,XX3,XX4,XX5,XX6,DEE1,DEE2,DEE3,DEE4,DEE5,DEE6, &
    float EL21, EL32, XX1, XX2, XX3, XX4, XX5, XX6, DEE1, DEE2, DEE3, DEE4, DEE5, DEE6;
    //F       FL1,FL2,FL3,XL,GAMH,GAMS,QS0,BTSOIL,FERTTYPE,TOTEM,CONVTERP, &
    float FL1, FL2, FL3, XL, GAMH, GAMS, QS0, BTSOIL;
    int FERTTYPE, TOTEM, CONVTERP;
    //F       R(4),CPART(4,5),DELMASS(4,226:iTp),ABFRAC(4,226:iTp)
    float R[ 4+1 ], CPART[ 4+1 ][ 5+1 ] /*, DELMASS[ 4+1 ][ iTp+1 ], ABFRAC[ 4+1 ][ iTp+1 ] */;
    magicc_array DELMASS, ABFRAC;
} CAR_block;

typedef struct {
    //F 152       COMMON/CARB/CCO2(4,224:iTp),EDGROSS(4,226:iTp),EF(226:iTp+1), &
    /* float CCO2[ 4+1 ][ iTp+1 ], EDGROSS[ 4+1 ][ iTp+1 ], EF[ iTp+1 ]; */
    magicc_array CCO2, EDGROSS, EF;
    //F 153       REGROW(4,226:iTp),PL(4,226:iTp),HL(4,226:iTp),SOIL(4,226:iTp), &
    magicc_array REGROW, PL, HL, SOIL;
    //F 154       TTT(226:iTp),ESUM(226:iTp),ETOT(4,226:iTp),EDNET90(4), &
    float EDNET90[ 4+1 ];
    magicc_array TTT, ESUM, ETOT;
    //F 155       FOC(4,226:iTp),CO2(0:iTp),CO2SAVE(0:iTp)
    float CO2[ iTp+1 ], CO2SAVE[ iTp+1 ];
    magicc_array FOC;
} CARB_block;

typedef struct {
    //F 150       COMMON/COBS/COBS(0:236)
    float COBS[ 236+1 ];    // Observed CO2 -IPCC DATA SET BY ENTING AND WIGLEY
} COBS_block;

typedef struct {
    //F 133       COMMON/CLIM/IC,IP,KC,DT,DZ,FK,HM,Q2X,QXX,PI,T,TE,TEND,W0,XK,XKLO, &
    int IC, IP, KC;
    float DT, DZ, FK, HM, Q2X, QXX, PI, T, TE, TEND, W0, XK, XKLO;
    //F 134       XKNS,XLAM,FL(2),FO(2),FLSUM,FOSUM,HEM(2),P(40),TEM(40),TO(2,40), &
    float XKNS, XLAM, FL[ 2+1 ], FO[ 2+1 ], FLSUM, FOSUM, HEM[ 2+1 ], P[ 40+1 ], TEM[ 40+1 ], 
        TO[ 2+1 ][ 40+1 ]; // OCEAN TEMP CHANGE [ HEMISPHERE ][ LEVEL ]
    //F 135       AL,BL,CL,DTH,DTZ,DZ1,XLL,WWW,XXX,YYY,RHO,SPECHT,HTCONS,Y(4)
    float   AL, BL, CL, DTH, DTZ, DZ1, XLL, WWW, XXX, YYY, RHO, SPECHT, HTCONS, 
            Y[ 4+1 ];   /* 1=northern ocean, 2=southern ocean, 3=northern land, 4=southern land land temps */
} CLIM_block;

typedef struct {
    //F 139       COMMON/CONCS/CH4(0:iTp),CN2O(0:iTp),ECH4(226:iTp+1), &
    float CH4[ iTp+1 ], CN2O[ iTp+1 ] /* , ECH4[ iTp+1+1 ] */;
    magicc_array ECH4;
    //F 140       EN2O(226:iTp+1),ECO(226:iTp+1),COE(iTp+1),EVOC(226:iTp+1), &
    magicc_array EN2O, ECO, COE, EVOC;
    //F 141       ENOX(226:iTp+1),ESO2(0:iTp+1),ESO2SUM(226:iTp+1), &
    magicc_array ENOX, ESO2, ESO2SUM;
    //F 142       ESO21(226:iTp+1),ESO22(226:iTp+1),ESO23(226:iTp+1), &
    magicc_array ESO21, ESO22, ESO23;
    //F 143       EBC(226:iTp+1), EOC(226:iTp+1) ! sjs- add BC-OC
    magicc_array EBC, EOC;
} CONCS_block;

typedef struct {
    //F 212       COMMON /CORREN/CORREN1,CORREN2,CORREN3,CORREN4,CORREN
    float CORREN1, CORREN2, CORREN3, CORREN4, CORREN;
} CORREN_block;

typedef struct {
    //F 187       COMMON /CO2READ/ICO2READ,XC(226:iTp),CO2SCALE,qtot86,LEVCO2
    float ICO2READ, XC[ iTp+1 ], CO2SCALE, qtot86;
    int LEVCO2;
} CO2READ_block;

typedef struct {
    //F 188       COMMON /DSENS/IXLAM,XLAML,XLAMO,ADJUST
    int IXLAM;
    float XLAML, XLAMO, ADJUST;
} DSENS_block;

typedef struct {
    //F 175       COMMON /FORCE/qco2(0:iTp),qm(0:iTp),qn(0:iTp),QCFC(0:iTp), &
    float QCO2[ iTp+1 ], QM[ iTp+1 ], QN[ iTp+1 ], QCFC[ iTp+1 ];
    //F 176       QMONT(0:iTp),QOTHER(0:iTp),QSTRATOZ(0:iTp),QCH4O3(0:iTp), &
    float QMONT[ iTp+1 ], QOTHER[ iTp+1 ], QSTRATOZ[ iTp+1 ], QCH4O3[ iTp+1 ];
    //F 177       CFC12(0:iTp), QCH4H2O(0:iTp),QBC(0:iTp),QOC(0:iTp)
    float CFC12[ iTp+1 ], QCH4H2O[ iTp+1 ], QBC[ iTp+1 ], QOC[ iTp+1 ];
} FORCE_block;

typedef struct {
    //F 215       COMMON /HALOF/QCF4_ar(0:iTp),QC2F6_ar(0:iTp),qSF6_ar(0:iTp), &
    float QCF4_ar[ iTp+1 ], QC2F6_ar[ iTp+1 ], qSF6_ar[ iTp+1 ];
    //F 216        Q125_ar(0:iTp),Q134A_ar(0:iTp), &
    float Q125_ar[ iTp+1 ], Q134A_ar[ iTp+1 ];
    //F 217        Q143A_ar(0:iTp),Q227_ar(0:iTp),Q245_ar(0:iTp)
    float Q143A_ar[ iTp+1 ], Q227_ar[ iTp+1 ], Q245_ar[ iTp+1 ];
} HALOF_block;

typedef struct {
    //F 196       COMMON /ICE/T1990,G1990,SEN,SENG,SENA,ERRG,ERRA, &
    float T1990, G1990, SEN, SENG, SENA, ERRG, ERRA;
    //F 197       DMG,DMA,SENI,SENP,SENS,DSENI,DSENP,DSENS,ICE,MODEL, &
    float DMG, DMA, SENI, SENP, SENS, DSENI, DSENP, DSENS;
    int ICE, MODEL;
    //F 198       NEWGSIC,IXG,VZERO,XG
    float NEWGSIC, IXG, VZERO, XG;
} ICE_block;

typedef struct {
    //F       COMMON /JSTART/JSTART,FOSSHIST(0:236),QKYMAG(0:iTp),IGHG, &
    int JSTART;
    float FOSSHIST[ 236+1 ], QKYMAG[ iTp+1 ];
    int IGHG;
    //F       QCH4OZ,QFOC(0:iTp),ICO2CORR,TROZSENS
    float QCH4OZ, QFOC[ iTp+1 ], ICO2CORR, TROZSENS;
} JSTART_block;

typedef struct {
    //F 129       common /Limits/KEND
    int KEND;
} Limits_block;

typedef struct {
    //F       COMMON /METH1/emeth(226:iTp),imeth,ch4l(225:iTp),ch4b(225:iTp), &
    float IMETH /*, emeth[ iTp+1 ], ch4l[ iTp+1 ], ch4b[ iTp+1 ] */;
    magicc_array emeth, ch4l, ch4b;
    //F       ch4h(225:iTp),ef4(226:iTp),StratH2O,TCH4(iTp),iO3feed, &
    magicc_array ch4h, ef4;
    float STRATH2O, TCH4[ iTp+1 ], IO3FEED;
    //F       ednet(226:iTp+1),DUSER,FUSER,CORRUSER,CORRMHI,CORRMMID,CORRMLO
    float /* ednet[ iTp+1 ], */ DUSER, FUSER, CORRUSER, CORRMHI, CORRMMID, CORRMLO;
    magicc_array ednet;
} METH1_block;    

typedef struct {
    //F       COMMON /METH2/LEVCH4,ch4bar90,QQQN2O
    int LEVCH4;
    float ch4bar90, QQQN2O;
} METH2_block;

typedef struct {
    //F       COMMON /METH3/TCH4CON,TAUINIT,SCH4,DELSS,DELTAU, &
    float TCH4CON, TAUINIT, SCH4, DELSS, DELTAU;
    //F       ANOX,ACO,AVOC,DELANOX,DELACO,DELAVOC,ICH4FEED
    float ANOX, ACO, AVOC, DELANOX, DELACO, DELAVOC, ICH4FEED;
} METH3_block;

typedef struct {
    //F       COMMON /METH4/GAM,TAUOTHER,BBCH4,CM00
    float GAM, TAUOTHER, BBCH4, CM00;
} METH4_block;

typedef struct {
    //F 145       COMMON/NEWCONCS/CF4(iTp),C2F6(iTp),C125(iTp),C134A(iTp), &
    float CF4[ iTp ], C2F6[ iTp ], C125[ iTp ], C134A[ iTp ];
    //F 146       C143A(iTp),C227(iTp),C245(iTp),CSF6(iTp), &
    float C143A[ iTp ], C227[ iTp ], C245[ iTp ], CSF6[ iTp+1 ];
    //F 147       ECF4(226:iTp+1),EC2F6(226:iTp+1),E125(226:iTp+1),E134A(226:iTp+1), &
    magicc_array ECF4, EC2F6, E125, E134A;
    //F 148       E143A(226:iTp+1),E227(226:iTp+1),E245(226:iTp+1),ESF6(226:iTp+1)
    magicc_array E143A, E227, E245, ESF6;
} NEWCONCS_block;    

typedef struct {
    float aNewClimSens, aNewBTsoil, DT2XUSER, aNewBTGPP, aNewBTHumus, 
        aNewDUSER, aNewFUSER, aNewSO2dir1990, aNewSO2ind1990;
} NEWPARAMS_block;

typedef struct {
    //F 206       COMMON /NSIM/NSIM,NCLIM,ISCENGEN,TEMEXP(2,40),IWNHOFF,IWSHOFF, &
    int NSIM, NCLIM, ISCENGEN;
    float TEMEXP[ 2+1 ][ 40+1 ];
    int IWNHOFF, IWSHOFF; 
    //F 207       WTHRESH
    float WTHRESH;
} NSIM_block;

typedef struct {
    //F       COMMON/OZ/OZ00CH4,OZCH4,OZNOX,OZCO,OZVOC
    float OZ00CH4, OZCH4, OZNOX, OZCO, OZVOC;
} OZ_block;

typedef struct {
    //F 202       COMMON /QADD/IQREAD,OrgIQREAD,JQFIRST,JQLAST,QEX(0:iTp),QEXNH(0:iTp), &
    int IQREAD, OrgIQREAD, JQFIRST, JQLAST;
    float QEX[ iTp+1 ], QEXNH[ iTp+1 ];
    //F 203       QEXSH(0:iTp),QEXNHO(0:iTp),QEXNHL(0:iTp),QEXSHO(0:iTp), &
    float QEXSH[ iTp+1 ], QEXNHO[ iTp+1 ], QEXNHL[ iTp+1 ], QEXSHO[ iTp+1 ];
    //F 204       QEXSHL(0:iTp),IOLDTZ
    float QEXSHL[ iTp+1 ];
    int IOLDTZ;
} QADD_block;

typedef struct {
    //F 191       COMMON /QSPLIT/QNHO,QNHL,QSHO,QSHL,QGLOBE(0:iTp), &
    float QNHO, QNHL, QSHO, QSHL, QGLOBE[ iTp+1 ];
    //F 192       QQNHO(0:iTp),QQNHL(0:iTp),QQSHO(0:iTp),QQSHL(0:iTp), &
    float QQNHO[ iTp+1 ], QQNHL[ iTp+1 ], QQSHO[ iTp+1 ], QQSHL[ iTp+1 ];
    //F 193       QQQNHO(0:iTp),QQQNHL(0:iTp),QQQSHO(0:iTp),QQQSHL(0:iTp), &
    float QQQNHO[ iTp+1 ], QQQNHL[ iTp+1 ], QQQSHO[ iTp+1 ], QQQSHL[ iTp+1 ];
    //F 194       EHistBC(iTp),EHistOC(iTp) ! Vars to store read-in BC-OC history.
    float EHistBC[ iTp+1 ], EHistOC[ iTp+1 ];
} QSPLIT_block;

typedef struct {
    //F 229       COMMON/STOREDVALS/ TEMUSER(iTp),QSO2SAVE(0:iTp+1),QDIRSAVE(0:iTp+1), &
    //F 230       KYRREF
    float TEMUSER[ iTp+1 ], QSO2SAVE[ iTp+1+1 ], QDIRSAVE[ iTp+1+1 ]
        , KYRREF
    ;
} STOREDVALS_block;

typedef struct {
    //F 186       common /Sulph/S90DIR,S90IND,S90BIO,enat,ES1990,ECO90,FOC90,IFOC
    float S90DIR, S90IND, S90BIO, ENAT, ES1990, ECO90, FOC90;
    int IFOC;
} Sulph_block;

typedef struct {
    //F 157       COMMON/TANDSL/TEQU(iTp),TGAV(iTp),TNHO(iTp), &
    float   TEQU[ iTp+1 ], 
            TGAV[ iTp+1 ],   /* Global mean temperature */
            TNHO[ iTp+1 ];   /* Mean temp, northern hemisphere ocean */
    //F 158       TSHO(iTp),TNHL(iTp),TSHL(iTp),TDEEP(iTp),TNHAV(iTp),TSHAV(iTp), &
    float   TSHO[ iTp+1 ],   /* Mean temp, southern hemisphere ocean */
            TNHL[ iTp+1 ],   /* Mean temp, northern hemisphere land */
            TSHL[ iTp+1 ],   /* Mean temp, southern hemisphere land */
            TDEEP[ iTp+1 ], 
            TNHAV[ iTp+1 ],  /* Mean temp, northern hemisphere */
            TSHAV[ iTp+1 ];  /* Mean temp, southern hemisphere */
    //F 159       TLAND(iTp),TOCEAN(iTp),TOCN(40),TOCNPREV(40), &
    float   TLAND[ iTp+1 ],    /* Mean temp, land */
            TOCEAN[ iTp+1 ],   /* Mean temp, ocean */
            TOCN[ 40+1 ],    // AREA WEIGHTED MEAN OF HEMIS TEMPS
            TOCNPREV[ 40+1 ];
    //F 160       SIP,SGP,SAP,SLI(iTp),SLG(iTp),SLA(iTp),EX(0:iTp),SLT(iTp), &
    float SLI[ iTp+1 ], SLG[ iTp+1 ], SLA[ iTp+1 ], 
            EX[ iTp+1 ],        /* Thermal expansion contribution to sea level change */
            SLT[ iTp+1 ];
    //F 161       QTOT(0:iTp),QGH(0:iTp),QOZ(0:iTp),QBIO(0:iTp),SLO(iTp), &
    float QTOT[ iTp+1 ], QGH[ iTp+1 ], QOZ[ iTp+1 ], QBIO[ iTp+1 ], SLO[ iTp+1 ];
    //F 162       QSO2(0:iTp+1),QDIR(0:iTp+1),QLAND(0:iTp),QMN(0:iTp+1)
    float QSO2[ iTp+1 ], QDIR[ iTp+1 ], QLAND[ iTp+1 ], QMN[ iTp+1 ];
} TANDSL_block;

typedef struct {
    //F       common /TauNitr/TN2000,BBN2O,SN2O,CN00,NOFFSET
    float TN2000, BBN2O, SN2O, CN00, NOFFSET;
} TauNitr_block;

typedef struct {
    //F 189       COMMON /VARW/Z(40),W(2),DW(2),TO0(2),TP0(2),WNH(iTp),WSH(iTp), &
    float Z[ 40+1 ], W[ 2+1 ], DW[ 2+1 ], TO0[ 2+1 ], TP0[ 2+1 ], WNH[ iTp+1 ], WSH[ iTp+1 ];
    //F 190       TW0NH,TW0SH,IVARW,KEYDW
    float TW0NH, TW0SH, IVARW;
    int KEYDW;
} VARW_block;

// Function prototypes
void CLIMAT();
void tslcalc( int N, Limits_block* Limits, CLIM_block* CLIM, CONCS_block* CONCS, CARB_block* CARB,
             TANDSL_block* TANDSL, VARW_block* VARW, QSPLIT_block* QSPLIT, ICE_block* ICE, 
             NSIM_block* NSIM, std::ofstream* outfile8 );
void init( Limits_block* Limits, CLIM_block* CLIM, CONCS_block* CONCS, TANDSL_block* TANDSL, FORCE_block* FORCE, 
          Sulph_block* Sulph, VARW_block* VARW, ICE_block* ICE, AREAS_block* AREAS, NSIM_block* NSIM,
          OZ_block* OZ, NEWCONCS_block* NEWCONCS, CARB_block* CARB, CAR_block* CAR, METH1_block* METH1,
          METH2_block* METH2, METH3_block* METH3, METH4_block* METH4, CO2READ_block* CO2READ, JSTART_block* JSTART,
          CORREN_block* CORREN, HALOF_block* HALOF, COBS_block* COBS, TauNitr_block* TauNitr );
void interp( int NVAL, int ISTART, int IY[], float X[], magicc_array* Y, int KEND );
void deltaq( Limits_block* Limits, OZ_block* OZ, CLIM_block* CLIM, CONCS_block* CONCS,
            NEWCONCS_block* NEWCONCS, CARB_block* CARB, TANDSL_block* TANDSL, CAR_block* CAR,
            METH1_block* METH1, FORCE_block* FORCE, METH2_block* METH2, METH3_block* METH3,
            METH4_block* METH4, TauNitr_block* TauNitr, Sulph_block* Sulph, NSIM_block* NSIM, 
            CO2READ_block* CO2READ, JSTART_block* JSTART, CORREN_block* CORREN, HALOF_block* HALOF, COBS_block* COBS );
void initcar( const int NN, const float D80, const float F80, COBS_block* COBS, 
             CARB_block* CARB, CAR_block* CAR );
void halocarb( const int N, float C0, float E, float* C1, float* Q, float TAU00, float TAUCH4 );
void history( const int JJJ, float* CO2, float* CH4, float* CN2O, float* eso2, float* eso21,
             float* CF4, float* C2F6, float* C125, float* C134A, float* C143A, float* C227, float* C245, float* CSF6,
             COBS_block* COBS, float ES1990 );
void methane( float ICH4F, float CPREV, float E, float DEN, float DEC, float DEV,
             float* CONC, float TAU00, float* TAUOUT, float S, float AANOX, 
             float AACO, float AAVOC, float TEMP, METH4_block* METH4 );
void nitrous( float C, float CP, float CPP, float E, float* C1, TauNitr_block* TauNitr );
void carbon( const int MM, float TEM, float EFOSS, float ENETDEF, float CPP, float CPREV, float C,
            float PL, float HU, float SO, float REGRO, float ETOT,
            float* PL1, float* HU1, float* SO1, float* REGRO1, float* ETOT1,
            float* SUMEM, float* FLUX, float* DELM, float* EGROSSD, float* C1,
            CAR_block* CAR ); 
void sulphate( const int JY, float ESO2, float ESO21, float ECO, float* QSO2, 
              float* QDIR, float* QFOC, float* QMN, Sulph_block* Sulph );
void lamcalc( float Q, float FNHL, float FSHL, float XK, float XKH, float DT2X, 
             float A, float* LAMOBEST, float* LAMLBEST );
void runmod( Limits_block* Limits, CLIM_block* CLIM, CONCS_block* CONCS, CARB_block* CARB, TANDSL_block* TANDSL,
            CO2READ_block* CO2READ, Sulph_block* Sulph, DSENS_block* DSENS, VARW_block* VARW, QSPLIT_block* QSPLIT,
            AREAS_block* AREAS, QADD_block* QADD, BCOC_block* BCOC, FORCE_block* FORCE, NSIM_block* NSIM,
            OZ_block* OZ, NEWCONCS_block* NEWCONCS, CAR_block* CAR, METH1_block* METH1, METH2_block* METH2, 
            METH3_block* METH3, METH4_block* METH4, TauNitr_block* TauNitr,
            JSTART_block* JSTART, CORREN_block* CORREN, HALOF_block* HALOF, COBS_block* COBS, ICE_block* ICE, std::ofstream* outfile8 );
void split( const float QGLOBE, const float A, const float BN, const float BS, float* QNO, float* QNL, 
           float* QSO, float* QSL, AREAS_block* AREAS );

void setGlobals( CARB_block* CARB, TANDSL_block* TANDSL, CONCS_block* CONCS, NEWCONCS_block* NEWCONCS, 
                STOREDVALS_block* STOREDVALS, NEWPARAMS_block* NEWPARAMS, BCOC_block* BCOC, 
                METH1_block* METH1, CAR_block* CAR, FORCE_block* FORCE, JSTART_block* JSTART,
                QADD_block* QADD, HALOF_block* HALOF );
void setLocals( CARB_block* CARB, TANDSL_block* TANDSL, CONCS_block* CONCS, NEWCONCS_block* NEWCONCS, 
                STOREDVALS_block* STOREDVALS, NEWPARAMS_block* NEWPARAMS, BCOC_block* BCOC, 
                METH1_block* METH1, CAR_block* CAR, FORCE_block* FORCE, JSTART_block* JSTART,
                QADD_block* QADD, HALOF_block* HALOF );

// Externally called methods

float getSLR( const int inYear );
float GETFORCING( const int iGasNumber, const int inYear );
float GETGHGCONC(int, int);
float GETGMTEMP(int);
float GETCARBONRESULTS(int, int);
void SETPARAMETERVALUES(int, float);
void overrideParameters( NEWPARAMS_block* NEWPARAMS, CAR_block* CAR, METH1_block* METH1, BCOC_block* BCOC );

// Internal helper methods

void openfile_read( std::ifstream* infile, const std::string& f, bool echo );
void skipline( std::ifstream* infile, bool echo );
float read_csv_value( std::ifstream* infile, bool echo );
float read_and_discard( std::ifstream* infile, bool echo );
void openfile_write( std::ofstream* outfile, const std::string& f, bool echo );


#endif // _ObjECTS_MAGICC.h_
