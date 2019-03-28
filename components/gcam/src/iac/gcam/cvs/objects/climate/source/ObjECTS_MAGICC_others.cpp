#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <cassert>
#include <limits>


#include "climate/include/ObjECTS_MAGICC.h"


using namespace std;

void f_change( const char* name, int entering )
{
/*     static int printmargin = 0;
     const char* spaces = "     ";
     static char lastname[20];
     static int callcount = 0;
     
     if( entering ) {
     printmargin += 1;
     if( strcmp( name, lastname ) ) {
     if( callcount > 1 ) {
     cout << " (" << callcount << " times)";
     } 
     cout << endl;
     for( int i=0; i<printmargin; i++)
     cout << spaces;
     cout << name;
     
     strcpy( lastname, name );   // moving to new function--reset things
     callcount = 0;
     
     } else {
     callcount++;
     }
     } else {    // exiting, little to do
     printmargin -= 1;
     }
 */
}

void f_enter( const char* name )
{
    f_change( name, 1 );
}

void f_exit( const char* name )
{
    f_change( name, 0 );
}

//F3058 !
//F3059 !********************************************************************
//F3060 !
//F3061       SUBROUTINE INIT
void init( Limits_block* Limits, CLIM_block* CLIM, CONCS_block* CONCS, TANDSL_block* TANDSL, FORCE_block* FORCE, 
          Sulph_block* Sulph, VARW_block* VARW, ICE_block* ICE, AREAS_block* AREAS, NSIM_block* NSIM,
          OZ_block* OZ, NEWCONCS_block* NEWCONCS, CARB_block* CARB, CAR_block* CAR, METH1_block* METH1,
          METH2_block* METH2, METH3_block* METH3, METH4_block* METH4, CO2READ_block* CO2READ, JSTART_block* JSTART,
          CORREN_block* CORREN, HALOF_block* HALOF, COBS_block* COBS, TauNitr_block* TauNitr )
{
    //    std::cout << "SUBROUTINE INIT" << endl;
    f_enter( __func__ );
    
    //F3062       IMPLICIT REAL*4 (a-h,o-z), Integer (I-N)
    //F3063 !
    //F3064       parameter (iTp=740)
    //F3065 !
    //F3066       common /Limits/KEND
    //F3067 !
    //F3068       COMMON/CLIM/IC,IP,KC,DT,DZ,FK,HM,Q2X,QXX,PI,T,TE,TEND,W0,XK,XKLO, &
    //F3069       XKNS,XLAM,FL(2),FO(2),FLSUM,FOSUM,HEM(2),P(40),TEM(40),TO(2,40), &
    //F3070       AL,BL,CL,DTH,DTZ,DZ1,XLL,WWW,XXX,YYY,RHO,SPECHT,HTCONS,Y(4)
    //F3071 !
    //F3072       COMMON/CONCS/CH4(0:iTp),CN2O(0:iTp),ECH4(226:iTp+1), &
    //F3073       EN2O(226:iTp+1),ECO(226:iTp+1),COE(iTp+1),EVOC(226:iTp+1), &
    //F3074       ENOX(226:iTp+1),ESO2(0:iTp+1),ESO2SUM(226:iTp+1), &
    //F3075       ESO21(226:iTp+1),ESO22(226:iTp+1),ESO23(226:iTp+1), &
    //F3076       EBC(226:iTp+1), EOC(226:iTp+1) ! sjs- add BC-OC
    //F3077 !
    //F3078       COMMON/TANDSL/TEQU(iTp),TGAV(iTp),TNHO(iTp), &
    //F3079       TSHO(iTp),TNHL(iTp),TSHL(iTp),TDEEP(iTp),TNHAV(iTp),TSHAV(iTp), &
    //F3080       TLAND(iTp),TOCEAN(iTp),TOCN(40),TOCNPREV(40), &
    //F3081       SIP,SGP,SAP,SLI(iTp),SLG(iTp),SLA(iTp),EX(0:iTp),SLT(iTp), &
    //F3082       QTOT(0:iTp),QGH(0:iTp),QOZ(0:iTp),QBIO(0:iTp),SLO(iTp), &
    //F3083       QSO2(0:iTp+1),QDIR(0:iTp+1),QLAND(0:iTp),QMN(0:iTp+1)
    //F3084 !
    //F3085       COMMON /FORCE/qco2(0:iTp),qm(0:iTp),qn(0:iTp),QCFC(0:iTp), &
    //F3086       QMONT(0:iTp),QOTHER(0:iTp),QSTRATOZ(0:iTp),QCH4O3(0:iTp), &
    //F3087       CFC12(0:iTp), QCH4H2O(0:iTp),QBC(0:iTp),QOC(0:iTp)
    //F3088 !
    //F3089       common /Sulph/S90DIR,S90IND,S90BIO,ENAT,ES1990,ECO90,FOC90,IFOC
    //F3090       COMMON /VARW/Z(40),W(2),DW(2),TO0(2),TP0(2),WNH(iTp),WSH(iTp), &
    //F3091       TW0NH,TW0SH,IVARW,KEYDW
    //F3092 !
    //F3093       COMMON /ICE/T1990,G1990,SEN,SENG,SENA,ERRG,ERRA, &
    //F3094       DMG,DMA,SENI,SENP,SENS,DSENI,DSENP,DSENS,ICE,MODEL, &
    //F3095       NEWGSIC,IXG,VZERO,XG
    //F3096 !
    //F3097       COMMON /AREAS/FNO,FNL,FSO,FSL
    //F3098 !
    //F3099       COMMON /NSIM/NSIM,NCLIM,ISCENGEN,TEMEXP(2,40),IWNHOFF,IWSHOFF, &
    //F3100       WTHRESH
    //F3101 !
    //F3102 !  THE PROGRAM USES VARIOUS COUNTERS TO KEEP TRACK OF TIME.
    //F3103 !   IC BEGINS WITH IC=1 TO IDENTIFY THE YEAR 1765. THUS IC
    //F3104 !   =226 IS THE YEAR 1990, IC=336 IS THE YEAR 2100, ETC.
    //F3105 !   CONC AND FORCING ARRAYS GIVE VALUES AT THE END OF THE
    //F3106 !   YEAR. THUS CONC(1) IS THE VALUE AT THE END OF 1765, ETC.
    //F3107 !   TIME (T) IS COUNTED FROM T=0 AT THE MIDDLE OF 1765, SO
    //F3108 !   THAT THE MIDDLE OF 1990 IS T=225.0. TEMP AND SEALEVEL
    //F3109 !   OUTPUT VALUES ARE AVERAGES OVER CALENDAR YEARS WITH THE
    //F3110 !   VALUES BEING THOSE CALCULATED AT THE MIDPOINTS OF YEARS.
    //F3111 !   TEMP(1) IS THEREFORE THE VALUE FOR THE MIDDLE OF 1765
    //F3112 !   CORRESPONDING TO T=0.0 AND IC=1. EMISSIONS ARE TOTALS
    //F3113 !   OVER CALENDAR YEARS, SO THE E(1) WOULD BE THE TOTAL FOR
    //F3114 !   1765, E(226) THE TOTAL FOR 1990, ETC.
    //F3115 !
    //F3116 !  ****************************************************************
    //F3117 !
    //F3118       FNL=FL(1)/2.0
    AREAS->FNL = CLIM->FL[ 1 ] / 2.0;
    //F3119       FNO=(1.0-FL(1))/2.0
    AREAS->FNO = (1.0 - CLIM->FL[ 1 ] ) / 2.0;
    //F3120       FSL=FL(2)/2.0
    AREAS->FSL = CLIM->FL[ 2 ] / 2.0;
    //F3121       FSO=(1.0-FL(2))/2.0
    AREAS->FSO = (1.0 - CLIM->FL[ 2 ] ) / 2.0;
    //F3122       FLSUM=2.0*(FNL+FSL)
    CLIM->FLSUM = 2.0 * ( AREAS->FNL + AREAS->FSL );
    //F3123       FOSUM=2.0*(FNO+FSO)
    CLIM->FOSUM = 2.0 * ( AREAS->FNO + AREAS->FSO );
    //F3124 !
    //F3125       IWNHOFF=0
    NSIM->IWNHOFF = 0;
    //F3126       IWSHOFF=0
    NSIM->IWSHOFF = 0;
    //F3127       IP=0
    CLIM->IP = 0;
    //F3128       IC=1
    CLIM->IC = 1;
    //F3129       KC=1
    CLIM->KC = 1;
    //F3130       T=0.0
    CLIM->T = 0.0;
    //F3131       DZ=100.
    CLIM->DZ = 100.0;
    //F3132       DZ1=DZ/2.
    CLIM->DZ1 = CLIM->DZ / 2.0;
    //F3133       DTH=DT/HM
    CLIM->DTH = CLIM->DT / CLIM->HM;
    //F3134       DTZ=DT/DZ
    CLIM->DTZ = CLIM->DT / CLIM->DZ;
    //F3135       XXX=XK/DZ1
    CLIM->XXX = CLIM->XK / CLIM->DZ1;
    //F3136       YYY=XK/DZ
    CLIM->YYY = CLIM->XK / CLIM->DZ;
    //F3137       AL=-DTZ*YYY
    CLIM->AL = -CLIM->DTZ * CLIM->YYY;
    //F3138 !
    //F3139 !  INITIALIZE TO(I,L) = OCEAN TEMP CHANGE IN HEMISPHERE "I"
    //F3140 !   AT LEVEL "L", AND TOCN(L) = AREA WEIGHTED MEAN OF HEMIS TEMPS.
    //F3141 !
    //F3142       DO L=1,40
    for( int L=1; L<=40; L++) {
        //F3143       TOCN(L)=0.0
        TANDSL->TOCN[ L ] = 0.0;
        //F3144         DO I=1,2
        for( int I=1; I<=2; I++) {
            //F3145         TO(I,L)=0.0
            CLIM->TO[ I ][ L ] = 0.0;
            //F3146         END DO
        }
        //F3147       END DO
    }
    //F3148 !
    //F3149       Y(1)=0.0
    //F3150       Y(2)=0.0
    //F3151       Y(3)=0.0
    //F3152       Y(4)=0.0
    CLIM->Y[ 1 ] = CLIM->Y[ 2 ] = CLIM->Y[ 3 ] = CLIM->Y[ 4 ] = 0.0;
    //F3153       HEM(1)=0.0
    //F3154       HEM(2)=0.0
    CLIM->HEM[ 1 ] = CLIM->HEM[ 2 ] = 0.0;
    //F3155 !
    //F3156 !  Z(I) DEPTH FROM BOTTOM OF MIXED LAYER FOR VARIABLE W
    //F3157 !
    //F3158       DO I=2,40
    for( int I=2; I<=40; I++) {
        //F3159       Z(I)=(I-2)*DZ+0.5*DZ
        VARW->Z[ I ] = ( I-2 ) * CLIM->DZ + 0.5 * CLIM->DZ;
        //F3146         END DO
    }
    //F3160       ENDDO
    //F3161 !
    //F3162 !  DEFINE INITIAL TEMP PROFILE (TEM(I)) AND PRESSURE (P(I)).
    //F3163 !
    //F3164       ZED=HM/2.
    float ZED = CLIM->HM / 2.0;
    //F3165       P(1)=0.0098*(0.1005*ZED+10.5*(EXP(-ZED/3500.)-1.))
    CLIM->P[ 1 ] = 0.0098 * ( 0.1005 * ZED + 10.5 * ( exp( float( -ZED/3500.0 ) )-1.0 ) );
    //F3166       TEM(1)=20.98-13.12/HM-0.04025*HM
    CLIM->TEM[ 1 ] = 20.98 - 13.12 / CLIM->HM - 0.04025 * CLIM->HM;
    //F3167       DO I=2,40
    for( int I=2; I<=40; I++ ) {
        //F3168       ZED=HM+(I-1)*DZ-DZ/2.
        ZED = CLIM->HM + ( I-1 ) * CLIM->DZ - CLIM->DZ / 2.0;
        //F3169       P(I)=0.0098*(0.1005*ZED+10.5*(EXP(-ZED/3500.)-1.))
        CLIM->P[ I ] = 0.0098*( 0.1005*ZED+10.5*( exp( float( -ZED/3500.0 ) )-1.0 ) );
        //F3170       ZZ=ZED/100.
        const float ZZ = ZED / 100.0;
        //F3171       IF(ZED.LE.130.)THEN
        if( ZED <= 130.0 ) {
            //F3172       TEM(I)=18.98-4.025*ZZ
            CLIM->TEM[ I ] = 18.98-4.025*ZZ;
            //F3173       ELSE IF(ZED.LE.500.0)THEN
        } else if( ZED <= 500 ) {
            //F3174       TEM(I)=17.01-2.829*ZZ+0.228*ZZ*ZZ
            CLIM->TEM[ I ] = 17.01 - 2.829 * ZZ + 0.228 * ZZ * ZZ;
            //F3175       ELSE IF(ZED.LT.2500.)THEN
        } else if( ZED < 2500 ) {
            //F3176       TEM(I)=EXP(EXP(1.007-0.0665*ZZ))
            CLIM->TEM[ I ] = exp( float( exp( float( 1.007-0.0665 * ZZ ) ) ) );
            //F3177       ELSE
        } else
            //F3178       TEM(I)=2.98-0.052*ZZ
            CLIM->TEM[ I ] = 2.98-0.052*ZZ;
        //F3179       ENDIF
        //F3180       END DO
    }
    //F3181 !
    //F3182 !  DEFINE THEORETICAL INITIAL TEMP PROFILE
    //F3183 !
    //F3184       TO0(1)=17.2      ! INITIAL MIXED LAYER TEMPE
    //F3185       TO0(2)=17.2      ! DITTO
    VARW->TO0[ 1 ] = VARW->TO0[ 2 ] = 17.2;
    //F3186       TP0(1)=1.0       ! INITIAL TEMP OF POLAR SINKING WATER
    //F3187       TP0(2)=1.0       ! DITTO
    VARW->TP0[ 1 ] = VARW->TP0[ 2 ] = 1.0;
    //F3188 !
    //F3189       DO I=1,2
    for( int I=1; I<=2; I++)
        //F3190       DO L=2,40
        for( int L=2; L<=40; L++ ) {
            //F3191       TEMEXP(I,L)=TP0(I)+(TO0(I)-TP0(I))*EXP(-W0*Z(L)/XK)
            NSIM->TEMEXP[ I ][ L ] = VARW->TP0[ I ]+( VARW->TO0[ I ]-VARW->TP0[ I ] ) 
            * exp( float( -CLIM->W0*VARW->Z[ L ]/CLIM->XK ) );
        }
    //F3192       END DO
    //F3193       END DO
    //F3194 !
    //F3195 !   SET INITIAL VALUES FOR USE WITH VARIABLE W
    //F3196 !
    //F3197       W(1)=W0
    //F3198       W(2)=W0
    VARW->W[ 1 ] = VARW->W[ 2 ] = CLIM->W0;
    //F3199       DW(1)=0.0
    //F3200       DW(2)=0.0
    VARW->DW[ 1 ] = VARW->DW[ 2 ] = 0.0;
    //F3201 !
    //F3202 !  DEFINE INITIAL TEMP AND SEA LEVEL COMPONENTS.
    //F3203 !
    //F3204       TGAV(1)=0.0
    TANDSL->TGAV[ 1 ] = 0.0;
    //F3205       TDEEP(1)=0.0
    TANDSL->TDEEP[ 1 ] = 0.0;
    //F3206       SLI(1)=0.0
    TANDSL->SLI[ 1 ] = 0.0;
    //F3207       SLG(1)=0.0
    TANDSL->SLG[ 1 ] = 0.0;
    //F3208       SLA(1)=0.0
    TANDSL->SLA[ 1 ] = 0.0;
    //F3209       SLT(1)=0.0
    TANDSL->SLT[ 1 ] = 0.0;
    //F3210       EX(0)=0.0
    TANDSL->EX[ 0 ] = 0.0;
    //F3211 !
    //F3212 !  CALCULATE NEW RADIATIVE FORCINGS EVERY FOURTH LOOP (I.E., WHEN
    //F3213 !   NSIM=1,5,9,13,17) WHEN NEW SO2 EMISSIONS ARE USED.
    //F3214 !
    //F3215       IF(NCLIM.EQ.1.OR.ISCENGEN.EQ.9)THEN
    if( NSIM->NCLIM == 1 || NSIM->ISCENGEN == 9 ) {
        //F3216 !
        //F3217         CALL DELTAQ
        deltaq( Limits, OZ, CLIM, CONCS,
               NEWCONCS, CARB, TANDSL, CAR,
               METH1, FORCE, METH2, METH3,
               METH4, TauNitr, Sulph, NSIM, CO2READ, JSTART,
               CORREN, HALOF, COBS );
        //F3218 !
        //F3219 !  INITIALISE QTOT ETC AT START OF 1765.
        //F3220 !  THIS ENSURES THAT ALL FORCINGS ARE ZERO AT THE MIDPOINT OF 1765.
        //F3221 !
        //F3222         QTOT(0)      = -qtot(1)
        TANDSL->QTOT[ 0 ]       = -TANDSL->QTOT[ 1 ];
        //F3223         qso2(0)      = -qso2(1)
        TANDSL->QSO2[ 0 ]       = -TANDSL->QSO2[ 1 ];
        //F3224         qdir(0)      = -qdir(1)
        TANDSL->QDIR[ 0 ]       = -TANDSL->QDIR[ 1 ];
        //F3225         qco2(0)      = -qco2(1)
        FORCE->QCO2[ 0 ]        = -FORCE->QCO2[ 1 ];
        //F3226         qm(0)        = -qm(1)
        FORCE->QM[ 0 ]          = -FORCE->QM[ 1 ];
        //F3227         qn(0)        = -qn(1)
        FORCE->QN[ 0 ]          = -FORCE->QN[ 1 ];
        //F3228         qcfc(0)      = -qcfc(1)
        FORCE->QCFC[ 0 ]        = -FORCE->QCFC[ 1 ];
        //F3229         QCH4O3(0)    = -QCH4O3(1)
        FORCE->QCH4O3[ 0 ]      = -FORCE->QCH4O3[ 1 ];
        //F3230         qgh(0)       = -qgh(1)
        TANDSL->QGH[ 0 ]        = -TANDSL->QGH[ 1 ];
        //F3231         QBIO(0)      = -QBIO(1)
        TANDSL->QBIO[ 0 ]       = -TANDSL->QBIO[ 1 ];
        //F3232         QLAND(0)     = -QLAND(1)
        TANDSL->QLAND[ 0 ]       = -TANDSL->QLAND[ 1 ];
        //F3233         QMN(0)       = -QMN(1)
        TANDSL->QMN[ 0 ]       = -TANDSL->QMN[ 1 ];
        //F3234       ENDIF
    }
    //F3235 !
    //F3236       RETURN
    //F3237       END
    f_exit( __func__ );
} // init
//F3238 !
//F3239 !  *******************************************************************
//F3240 !
//F3241       SUBROUTINE TSLCALC(N)
void tslcalc( int N, Limits_block* Limits, CLIM_block* CLIM, CONCS_block* CONCS, CARB_block* CARB,
             TANDSL_block* TANDSL, VARW_block* VARW, QSPLIT_block* QSPLIT, ICE_block* ICE, 
             NSIM_block* NSIM, std::ofstream* outfile8 )
{
    //    std::cout << "SUBROUTINE TSLCALC" << endl;
    f_enter( __func__ );
    //F3242       IMPLICIT REAL*4 (a-h,o-z), Integer (I-N)
    //F3243 !
    //F3244       parameter (iTp=740)
    //F3245 !
    //F3246       DIMENSION SLI1(iTp),SLG1(iTp),SLA1(iTp),SLO1(iTp),ZALL1(iTp), &
    //F3247       EX1(iTp),ZALL(iTp), &
    //F3248       SLI2(iTp),SLG2(iTp),SLA2(iTp),SLO2(iTp),ZALL2(iTp),EX2(iTp), &
    //F3249       SLI3(iTp),SLG3(iTp),SLA3(iTp),SLO3(iTp),ZALL3(iTp),EX3(iTp), &
    //F3250       SLI4(iTp),SLG4(iTp),SLA4(iTp),SLO4(iTp),ZALL4(iTp),EX4(iTp)
    float SLI1[ iTp+1 ], SLG1[ iTp+1 ], SLA1[ iTp+1 ], SLO1[ iTp+1 ], ZALL1[ iTp+1 ], EX1[ iTp+1 ], ZALL[ iTp+1 ];
    float SLI2[ iTp+1 ], SLG2[ iTp+1 ], SLA2[ iTp+1 ], SLO2[ iTp+1 ], ZALL2[ iTp+1 ], EX2[ iTp+1 ];
    float SLI3[ iTp+1 ], SLG3[ iTp+1 ], SLA3[ iTp+1 ], SLO3[ iTp+1 ], ZALL3[ iTp+1 ], EX3[ iTp+1 ];
    float SLI4[ iTp+1 ], SLG4[ iTp+1 ], SLA4[ iTp+1 ], SLO4[ iTp+1 ], ZALL4[ iTp+1 ], EX4[ iTp+1 ];
    //F3251 !
    //F3252       common /Limits/KEND
    //F3253 !
    //F3254       COMMON/CLIM/IC,IP,KC,DT,DZ,FK,HM,Q2X,QXX,PI,T,TE,TEND,W0,XK,XKLO, &
    //F3255       XKNS,XLAM,FL(2),FO(2),FLSUM,FOSUM,HEM(2),P(40),TEM(40),TO(2,40), &
    //F3256       AL,BL,CL,DTH,DTZ,DZ1,XLL,WWW,XXX,YYY,RHO,SPECHT,HTCONS,Y(4)
    //F3257 !
    //F3258       COMMON/CONCS/CH4(0:iTp),CN2O(0:iTp),ECH4(226:iTp+1), &
    //F3259       EN2O(226:iTp+1),ECO(226:iTp+1),COE(iTp+1),EVOC(226:iTp+1), &
    //F3260       ENOX(226:iTp+1),ESO2(0:iTp+1),ESO2SUM(226:iTp+1), &
    //F3261       ESO21(226:iTp+1),ESO22(226:iTp+1),ESO23(226:iTp+1), &
    //F3262       EBC(226:iTp+1), EOC(226:iTp+1) ! sjs- add BC-OC
    //F3263 !
    //F3264       COMMON/CARB/CCO2(4,224:iTp),EDGROSS(4,226:iTp),EF(226:iTp+1), &
    //F3265       REGROW(4,226:iTp),PL(4,226:iTp),HL(4,226:iTp),SOIL(4,226:iTp), &
    //F3266       TTT(226:iTp),ESUM(226:iTp),ETOT(4,226:iTp),EDNET90(4), &
    //F3267       FOC(4,226:iTp),CO2(0:iTp),CO2SAVE(0:iTp)
    //F3268 !
    //F3269       COMMON/TANDSL/TEQU(iTp),TGAV(iTp),TNHO(iTp), &
    //F3270       TSHO(iTp),TNHL(iTp),TSHL(iTp),TDEEP(iTp),TNHAV(iTp),TSHAV(iTp), &
    //F3271       TLAND(iTp),TOCEAN(iTp),TOCN(40),TOCNPREV(40), &
    //F3272       SIP,SGP,SAP,SLI(iTp),SLG(iTp),SLA(iTp),EX(0:iTp),SLT(iTp), &
    //F3273       QTOT(0:iTp),QGH(0:iTp),QOZ(0:iTp),QBIO(0:iTp),SLO(iTp), &
    //F3274       QSO2(0:iTp+1),QDIR(0:iTp+1),QLAND(0:iTp),QMN(0:iTp+1)
    //F3275 !
    //F3276       COMMON /VARW/Z(40),W(2),DW(2),TO0(2),TP0(2),WNH(iTp),WSH(iTp), &
    //F3277       TW0NH,TW0SH,IVARW,KEYDW
    //F3278       COMMON /QSPLIT/QNHO,QNHL,QSHO,QSHL,QGLOBE(0:iTp), &
    //F3279       QQNHO(0:iTp),QQNHL(0:iTp),QQSHO(0:iTp),QQSHL(0:iTp), &
    //F3280       QQQNHO(0:iTp),QQQNHL(0:iTp),QQQSHO(0:iTp),QQQSHL(0:iTp), &
    //F3281       EHistBC(iTp),EHistOC(iTp) ! Vars to store read-in BC-OC history.
    //F3282 
    //F3283 !
    //F3284       COMMON /ICE/T1990,G1990,SEN,SENG,SENA,ERRG,ERRA, &
    //F3285       DMG,DMA,SENI,SENP,SENS,DSENI,DSENP,DSENS,ICE,MODEL, &
    //F3286       NEWGSIC,IXG,VZERO,XG
    //F3287 !
    //F3288       COMMON /NSIM/NSIM,NCLIM,ISCENGEN,TEMEXP(2,40),IWNHOFF,IWSHOFF, &
    //F3289       WTHRESH
    //F3290 !
    //F3291       QQNHO(N)  = QNHO
    QSPLIT->QQNHO[ N ] = QSPLIT->QNHO;
    //F3292       QQNHL(N)  = QNHL
    QSPLIT->QQNHL[ N ] = QSPLIT->QNHL;
    //F3293       QQSHO(N)  = QSHO
    QSPLIT->QQSHO[ N ] = QSPLIT->QSHO;
    //F3294       QQSHL(N)  = QSHL
    QSPLIT->QQSHL[ N ] = QSPLIT->QSHL;
    //F3295       TNHO(N)  = Y(1)
    TANDSL->TNHO[ N ] = CLIM->Y[ 1 ];
    //F3296       TNHL(N)  = Y(3)
    TANDSL->TNHL[ N ] = CLIM->Y[ 3 ];
    //F3297       TSHO(N)  = Y(2)
    TANDSL->TSHO[ N ] = CLIM->Y[ 2 ];
    //F3298       TSHL(N)  = Y(4)
    TANDSL->TSHL[ N ] = CLIM->Y[ 4 ];
    //F3299       TNHAV(N) = FO(1)*Y(1)+FL(1)*Y(3)
    TANDSL->TNHAV[ N ] = CLIM->FO[ 1 ] * CLIM->Y[ 1 ] + CLIM->FL[ 1 ] * CLIM->Y[ 3 ];
    //F3300       TSHAV(N) = FO(2)*Y(2)+FL(2)*Y(4)
    TANDSL->TSHAV[ N ] = CLIM->FO[ 2 ] * CLIM->Y[ 2 ] + CLIM->FL[ 2 ] * CLIM->Y[ 4 ];
    //F3301       TLAND(N) = (FL(1)*Y(3)+FL(2)*Y(4))/FLSUM
    TANDSL->TLAND[ N ] = ( CLIM->FL[ 1 ] * CLIM->Y[ 3 ] + CLIM->FL[ 2 ] * CLIM->Y[ 4 ] ) / CLIM->FLSUM;
    //F3302       TOCEAN(N)= (FO(1)*Y(1)+FO(2)*Y(2))/FOSUM
    TANDSL->TOCEAN[ N ] = ( CLIM->FO[ 1 ] * CLIM->Y[ 1 ] + CLIM->FO[ 2 ] * CLIM->Y[ 2 ] ) / CLIM->FOSUM;
    //F3303       TGAV(N)  = (TNHAV(N)+TSHAV(N))/2.
    TANDSL->TGAV[ N ] = ( TANDSL->TNHAV[ N ] + TANDSL->TSHAV[ N ] ) / 2.0;
    //F3304 !
    //F3305 !  CALCULATE MEAN OCEAN TEMP CHANGE PROFILE AND ITS INCREMENTAL
    //F3306 !   CHANGE OVER ONE YEAR
    //F3307 !
    //F3308       DO L=1,40
    for( int L=1; L<=40; L++ ) {
        //F3309         TOCNPREV(L)=TOCN(L)
        TANDSL->TOCNPREV[ L ] = TANDSL->TOCN[ L ];
        //F3310         TOCN(L)=(FO(1)*TO(1,L)+FO(2)*TO(2,L))/FOSUM
        TANDSL->TOCN[ L ] = ( CLIM->FO[ 1 ] * CLIM->TO[ 1 ][ L ] + CLIM->FO[ 2 ] * CLIM->TO[ 2 ][ L ] ) / CLIM->FOSUM;
        //F3311       END DO
    }
    //F3312       TDEEP(N)=TOCN(40)
    TANDSL->TDEEP[ N ] = TANDSL->TOCN[ 40 ];
    //F3313 !
    //F3314 !  INCREMENTAL THERMAL EXPN CONTRIBUTION TO SEA LEVEL CHANGE.
    //F3315 !
    //F3316       EXPAN=EX(N-1)
    float EXPAN = TANDSL->EX[ N-1 ];
    //F3317 !
    //F3318       ZLAYER=HM
    float ZLAYER = CLIM->HM;
    //F3319 !
    //F3320       DO I=1,40
    for( int I=1; I<=40; I++ ) {
        //F3321 !
        //F3322         DELTOC=TOCN(I)-TOCNPREV(I)
        float DELTOC = TANDSL->TOCN[ I ] - TANDSL->TOCNPREV[ I ];
        //F3323         DELTBAR=(TOCN(I)+TOCNPREV(I))/2.0
        float DELTBAR = ( TANDSL->TOCN[ I ] + TANDSL->TOCNPREV[ I ] ) / 2.0;
        //F3324 !
        //F3325         TL1=TEM(I)+DELTBAR
        float TL1 = CLIM->TEM[ I ] + DELTBAR;
        //F3326         TL2=TL1*TL1
        float TL2 = TL1 * TL1;
        //F3327         TL3=TL2*TL1/6000.
        float TL3 = TL2 * TL1 / 6000.0;
        //F3328         PPL=P(I)
        float PPL = CLIM->P[ I ];
        //F3329         COEFE=52.55+28.051*PPL-0.7108*PPL*PPL+TL1*(12.9635-1.0833*PPL)- &
        //F3330          TL2*(0.1713-0.019263*PPL)+TL3*(10.41-1.1338*PPL)
        float COEFE = 52.55 + 28.051 * PPL - 0.7108 * PPL * PPL + TL1 * ( 12.9635 - 1.0833*PPL ) - TL2 * ( 0.1713 - 0.019263*PPL ) + TL3 * ( 10.41 - 1.1338*PPL );
        //F3331         DLR=COEFE*DELTOC*ZLAYER/10000.
        float DLR = COEFE * DELTOC * ZLAYER / 10000.0;
        //F3332         ZLAYER=DZ
        ZLAYER = CLIM->DZ;
        //F3333         EXPAN=EXPAN+DLR
        EXPAN += DLR;
        //F3334 !
        //F3335       END DO
    }
    //F3336 !
    //F3337       EX(N)=EXPAN
    TANDSL->EX[ N ] = EXPAN;
    //F3338 !
    //F3339 !  ICE MELT CONTRIBUTIONS TO SEA LEVEL RISE.
    //F3340 !  FOR TAR MODELS, KEY TEMPERATURE PROBABLY SHOULD BE CHANGE
    //F3341 !   FROM 1880, BUT CHANGE FROM PRE-INDUISTRIAL TIMES WAS USED.
    //F3342 !  FOR MAGICC (MODEL=0) T1990=TBASE
    //F3343 !
    //F3344       IF(N.LE.226)THEN
    float TBAR = 0.0;
    static float TCUM = 0.0;
    if( N <= 226 ) {
        //F3345         TBAR = 0.0
        //F3346         TCUM = 0.0
        TBAR = TCUM = 0.0;
        //F3347         SLT(N)=EX(N)
        TANDSL->SLT[ N ] = TANDSL->EX[ N ];
        //F3348       ENDIF
    }
    //F3349 !
    static float TBASE = 0.0, XX = 0.0, GS1990 = 0.0, B19901 = 0.0, B19902 = 0.0, B19903 = 0.0, B19904 = 0.0;
    static float BZERO1 = 0.0, BZERO2 = 0.0, BZERO3 = 0.0, BZERO4 = 0.0;
    static float GSPREV1 = 0.0, GSPREV2 = 0.0, GSPREV3 = 0.0, GSPREV4 = 0.0;
    static float VZ1 = 0.0, VZ2 = 0.0, VZ3 = 0.0, VZ4 = 0.0;
    float GS, GS1, GS2, GS3, GS4;
    GS = GS1 = GS2 = GS3 = GS4 = 0.0;
    //F3350       IF(N.EQ.226)THEN
    if( N == 226 ) {
        //F3351         TBASE=TGAV(N)
        TBASE = TANDSL->TGAV[ N ];
        //F3352         XX=G1990
        XX = ICE->G1990;
        //F3353         GS1990=0.934*XX-0.01165*XX*XX        ! CM
        GS1990 = 0.934 * XX - 0.01165 * XX * XX;
        //F3354         B19901=(0.934-0.0233*XX)*SEN*0.6                 ! NEW CODE
        B19901 = ( 0.934 - 0.0233 * XX ) * ICE->SEN * 0.6;
        //F3355         B19902=(0.934-0.0233*XX)*SEN                     ! NEW CODE
        B19902 = ( 0.934 - 0.0233 * XX ) * ICE->SEN;
        //F3356         B19903=(0.934-0.0233*XX)*SEN*1.4                 ! NEW CODE
        B19903 = ( 0.934 - 0.0233 * XX ) * ICE->SEN * 1.4;
        //F3357         B19904=(0.934-0.0233*XX)*SEN*(1.0+(ICE-2)*0.4)   ! NEW CODE
        B19904 = ( 0.934 - 0.0233 * XX ) * ICE->SEN * ( 1.0 + ( ICE->ICE - 2 ) * 0.4 );
        //F3358 !
        //F3359 !  YG IS TO ALLOW BZERO SCALING TO BE TURNED OFF
        //F3360 !
        //F3361         YG=XG                                            ! NEW CODE
        float YG = ICE->XG;
        //F3362         IF(IXG.EQ.0)YG=1.0                               ! NEW CODE
        if( ICE->IXG == 0 ) YG = 1.0;
        //F3363 !
        //F3364 !  ERROR BOUNDS ON VZERO CHANGED FROM +/-10 TO FOLLOW AR4 (JUNE 2008)
        //F3365 !
        //F3366         VZ1=VZERO-11.                                    ! NEW CODE
        VZ1 = ICE->VZERO - 11.0;
        //F3367         VZ2=VZERO                                        ! NEW CODE
        VZ2 = ICE->VZERO;
        //F3368         VZ3=VZERO+15.                                    ! NEW CODE
        VZ3 = ICE->VZERO + 15.0;
        //F3369         VZ4=VZERO+(ICE-2)*11.                            ! NEW CODE
        VZ4 = ICE->VZERO + ( ICE->ICE - 2 ) * 11.0;
        //F3370         IF(ICE.EQ.3)VZ4=VZ4+4.
        if( ICE->ICE == 3 ) VZ4 += 4.0;
        //F3371         BZERO1=B19901/((1.0-GS1990/VZ1)**YG)             ! NEW CODE
        BZERO1 = B19901 / pow ( static_cast<float> (1.0) - GS1990 / VZ1, YG );
        //F3372         BZERO2=B19902/((1.0-GS1990/VZ2)**YG)             ! NEW CODE
        BZERO2 = B19902 / pow ( static_cast<float> (1.0) - GS1990 / VZ2, YG );
        //F3373         BZERO3=B19903/((1.0-GS1990/VZ3)**YG)             ! NEW CODE
        BZERO3 = B19903 / pow ( static_cast<float> (1.0) - GS1990 / VZ3, YG );
        //F3374         BZERO4=B19904/((1.0-GS1990/VZ4)**YG)             ! NEW CODE
        BZERO4 = B19904 / pow ( static_cast<float> (1.0) - GS1990 / VZ4, YG );
        //F3375         GSPREV1=GS1990                                   ! NEW CODE
        GSPREV1 = GS1990;
        //F3376         GSPREV2=GS1990                                   ! NEW CODE
        GSPREV2 = GS1990;
        //F3377         GSPREV3=GS1990                                   ! NEW CODE
        GSPREV3 = GS1990;
        //F3378         GSPREV4=GS1990                                   ! NEW CODE
        GSPREV4 = GS1990;
        //F3379 !        WRITE(8,*)SEN,VZ4,B19904,BZERO4,GSPREV4
        //F3380       ENDIF
    }
    //F3381       IF(MODEL.EQ.0)T1990=TBASE
    if( ICE->MODEL == 0 ) ICE->T1990 = TBASE;
    //F3382 !
    //F3383 !  NEED TCUM = INTEGRAL OF TEMP CHANGE FROM MID 1990
    //F3384 !
    //F3385       IF(N.GE.226)THEN
    if( N >= 226 ) {
        //F3386 !        TBAR=(TGAV(N)+TGAV(N-1))/2.0-TGAV(226)
        //F3387         TBAR=(TGAV(N)+TGAV(N-1))/2.0
        TBAR = ( TANDSL->TGAV[ N ] + TANDSL->TGAV[ N-1 ] ) / 2.0;
        //F3388         IF(N.EQ.226)TBAR=0.0
        if( N == 226 ) TBAR = 0.0;
        //F3389         TCUM=TCUM+TBAR
        TCUM += TBAR;
        //F3390         DTB=0.15
        const float DTB = 0.15;
        //F3391         AAA=T1990-TBASE
        const float AAA = ICE->T1990 - TBASE;
        //F3392         DYR=FLOAT(N-226)
        const float DYR = float( N - 226 );
        //F3393         BBB=AAA*DYR+TCUM
        const float BBB = AAA * DYR + TCUM;
        //F3394 !
        //F3395 !  NEW GSIC. NOTE THAT LO, MID, HIGH AND USER CASES MUST ALL BE
        //F3396 !   CARRIED THRU TOGETHER SINCE THEY CANNOT BE CALCULATED BY
        //F3397 !   SIMPLY SCALING THE CENTRAL VALUE AS PREVIOUSLY.
        //F3398 !
        //F3399         IF(NEWGSIC.EQ.1)THEN                                     ! NEW CODE
        if( ICE->NEWGSIC == 1 ) {
            //F3400 !
            //F3401           FF1=BZERO1*(DTB+AAA+TGAV(N))                           ! NEW CODE
            const float FF1 = BZERO1 * ( DTB + AAA + TANDSL->TGAV[ N ] );
            //F3402           FF2=BZERO2*(DTB+AAA+TGAV(N))                           ! NEW CODE
            const float FF2 = BZERO2 * ( DTB + AAA + TANDSL->TGAV[ N ] );
            //F3403           FF3=BZERO3*(DTB+AAA+TGAV(N))                           ! NEW CODE
            const float FF3 = BZERO3 * ( DTB + AAA + TANDSL->TGAV[ N ] );
            //F3404           FF4=BZERO4*(DTB+AAA+TGAV(N))                           ! NEW CODE
            const float FF4 = BZERO4 * ( DTB + AAA + TANDSL->TGAV[ N ] );
            //F3405           X1=1.0-GSPREV1/VZ1                                     ! NEW CODE
            const float X1 = 1.0 - GSPREV1 / VZ1;
            //F3406           X2=1.0-GSPREV2/VZ2                                     ! NEW CODE
            const float X2 = 1.0 - GSPREV2 / VZ2;
            //F3407           X3=1.0-GSPREV3/VZ3                                     ! NEW CODE
            const float X3 = 1.0 - GSPREV3 / VZ3;
            //F3408           X4=1.0-GSPREV4/VZ4                                     ! NEW CODE
            const float X4 = 1.0 - GSPREV4 / VZ4;
            //F3409           DEL1=FF1*(X1**XG)/(1.0+0.5*FF1*XG*(X1**(XG-1.0))/VZ1)  ! NEW CODE
            const float DEL1 = FF1 * pow( X1, ICE->XG ) / ( 1.0 + 0.5 * FF1 * ICE->XG * ( pow( X1, ICE->XG- static_cast<float> (1.0) ) ) / VZ1 ); 
            //F3410           DEL2=FF2*(X2**XG)/(1.0+0.5*FF2*XG*(X2**(XG-1.0))/VZ2)  ! NEW CODE
            const float DEL2 = FF2 * pow( X2, ICE->XG ) / ( 1.0 + 0.5 * FF2 * ICE->XG * ( pow( X2, ICE->XG- static_cast<float> (1.0) ) ) / VZ2 ); 
            //F3411           DEL3=FF3*(X3**XG)/(1.0+0.5*FF3*XG*(X3**(XG-1.0))/VZ3)  ! NEW CODE
            const float DEL3 = FF3 * pow( X3, ICE->XG ) / ( 1.0 + 0.5 * FF3 * ICE->XG * ( pow( X3, ICE->XG- static_cast<float> (1.0) ) ) / VZ3 ); 
            //F3412           DEL4=FF4*(X4**XG)/(1.0+0.5*FF4*XG*(X4**(XG-1.0))/VZ4)  ! NEW CODE
            const float DEL4 = FF4 * pow( X4, ICE->XG ) / ( 1.0 + 0.5 * FF4 * ICE->XG * ( pow( X4, ICE->XG- static_cast<float> (1.0) ) ) / VZ4 ); 
            //F3413           GS1=GSPREV1+DEL1                                       ! NEW CODE
            GS1 = GSPREV1 + DEL1;
            //F3414           GS2=GSPREV2+DEL2                                       ! NEW CODE
            GS2 = GSPREV2 + DEL2;
            //F3415           GS3=GSPREV3+DEL3                                       ! NEW CODE
            GS3 = GSPREV3 + DEL3;
            //F3416           GS4=GSPREV4+DEL4                                       ! NEW CODE
            GS4 = GSPREV4 + DEL4;
            //F3417 !
            //F3418           GSPREV1=GS1                                            ! NEW CODE
            GSPREV1 = GS1;
            //F3419           GSPREV2=GS2                                            ! NEW CODE
            GSPREV2 = GS2;
            //F3420           GSPREV3=GS3                                            ! NEW CODE
            GSPREV3 = GS3;
            //F3421           GSPREV4=GS4                                            ! NEW CODE
            GSPREV4 = GS4;
            //F3422         ELSE                                                     ! NEW CODE
        } else {
            //F3423 !
            //F3424 !  GSICs, GS = MELT CONTRIB FROM 1880, DGS = 1-SIGMA UNCERT.,
            //F3425 !   SEN = GSIC SENSITIVITY.
            //F3426 !
            //F3427           GU=G1990+DTB*SEN*DYR+SEN*BBB         ! THIS IS IN CM
            float GU = ICE->G1990 + DTB * ICE->SEN * DYR + ICE->SEN * BBB;
            //F3428 !
            //F3429 !  TRAP TO AVOID UNREALISTIC BEHAVIOR IN GS FOR LARGE GU. THIS
            //F3430 !   IMPOSES AN UPPER BOUND ON THE CENTRAL GU VALUE OF 40.1 CM, WITH
            //F3431 !   A CORRESP. UPPER BOUND ON GS OF 18.7201 CM.
            //F3432 !
            //F3433           IF(GU.GT.40.0858)GU=40.0858          ! GU IS IN CM
            if( GU > 40.0858 ) GU = 40.0858;
            //F3434           GUM=GU/100.                          ! THIS IS IN M
            const float GUM = GU / 100.0;
            //F3435           GSM=0.934*GUM-1.165*GUM*GUM          ! AREA CORRECTED, IN M
            const float GSM = 0.934 * GUM - 1.165 * GUM * GUM;
            //F3436           GS=GSM*100.                          ! BACK TO CM
            GS = GSM * 100;
            //F3437 !
            //F3438         ENDIF                                                    ! NEW CODE
        }
        //F3439 !
        //F3440         GREF=GS-GS1990
        const float GREF = GS - GS1990;
        //F3441         DGS=0.40*GREF
        float DGS = 0.40 * GREF;
        //F3442         IF(NEWGSIC.EQ.1)DGS=(GS3-GS1)/2.0                        ! NEW CODE
        if( ICE->NEWGSIC == 1 ) DGS = ( GS3 - GS1 ) / 2.0;
        //F3443 !
        //F3444 !  GREENLAND AND ANTARCTICA  (ZGR AND ZAN, FROM 1990)
        //F3445 !
        //F3446         ZGR=SENG*BBB
        const float ZGR = ICE->SENG * BBB;
        //F3447         DZGR1=ERRG*DMG*BBB
        const float DZGR1 = ICE->ERRG * ICE->DMG * BBB;
        //F3448         DZGR2=0.1*ZGR
        const float DZGR2 = 0.1 * ZGR;
        //F3449 !
        //F3450         ZAN=SENA*BBB
        const float ZAN = ICE->SENA * BBB;
        //F3451         DZAN=ERRA*DMA*BBB
        const float DZAN = ICE->ERRA * ICE->DMA * BBB;
        //F3452 !
        //F3453 !  ICE MELT UNCERTAINTY TERM
        //F3454 !
        //F3455         DHV=SQRT(DGS**2+DZGR1**2+DZGR2**2+DZAN**2)
        const float DHV = sqrt( float( pow( DGS, 2) + pow( DZGR1, 2 ) + pow( DZGR2, 2 ) + pow( DZAN, 2 ) ) );
        //F3456 !
        //F3457 !  OTHER CONTRIBUTORS
        //F3458 !
        //F3459         ZI=SENI*DYR
        const float ZI = ICE->SENI * DYR;
        //F3460         DZI=DSENI*DYR
        const float DZI = ICE->DSENI * DYR;
        //F3461         ZP=SENP*DYR
        const float ZP = ICE->SENP * DYR;
        //F3462         DZP=DSENP*DYR
        const float DZP = ICE->DSENP * DYR;
        //F3463         ZS=SENS*DYR
        const float ZS = ICE->SENS * DYR;
        //F3464         DZS=DSENS*DYR
        const float DZS = ICE->DSENS * DYR;
        //F3465 !
        //F3466 !  TOTAL NON-EXPANSION SEA LEVEL RISE, FROM 1880 BUT IGNORING
        //F3467 !   ALL CHANGES OVER 1880-1990 EXCEPT GSICs
        //F3468 !
        //F3469         ZTOT=GS+ZGR+ZAN+ZI+ZP+ZS
        float ZTOT = GS + ZGR + ZAN + ZI + ZP + ZS;
        //F3470         IF(NEWGSIC.EQ.1)ZTOT=GS2+ZGR+ZAN+ZI+ZP+ZS       ! REVISED CODE
        if( ICE->NEWGSIC == 1 ) ZTOT = GS2 + ZGR + ZAN + ZI + ZP + ZS;
        //F3471 !
        //F3472 !  UNCERTAINTY TERM (2-SIGMA)
        //F3473 !
        //F3474         DZTOT=2.0*DHV+DZI+DZP+DZS
        const float DZTOT = 2.0 * DHV + DZI + DZP + DZS;
        //F3475 !
        //F3476 !  ABOVE GIVES MID YEAR VALUES.
        //F3477 !
        //F3478         SLI(N)=GS
        TANDSL->SLI[ N ] = GS;
        //F3479         IF(NEWGSIC.EQ.1)SLI(N)=GS2           ! NEW CODE
        if( ICE->NEWGSIC == 1 ) TANDSL->SLI[ N ] = GS2;
        //F3480         SLG(N)=ZGR
        TANDSL->SLG[ N ] = ZGR;
        //F3481         SLA(N)=ZAN
        TANDSL->SLA[ N ] = ZAN;
        //F3482         SLO(N)=ZI+ZP+ZS
        TANDSL->SLO[ N ] = ZI + ZP + ZS;
        //F3483         ZALL(N)=ZTOT
        ZALL[ N ] = ZTOT;
        //F3484 !
        //F3485 !  CASES:
        //F3486 !   ICE MELT ONLY CALCULATED FOR NSIM=1,2,3 OR 4.
        //F3487 !   NCLIM=1 = LOW = LOW EXPAN + LOW MELT
        //F3488 !
        //F3489        IF(NSIM.EQ.1)THEN
        if( NSIM->NSIM == 1 ) {
            //F3490           SLI(N)=GS-2.0*DGS
            TANDSL->SLI[ N ] = GS - 2.0 * DGS;
            //F3491 !                       
            //F3492 ! THE FACTOR 2.0 IS BECAUSE THIS IS THE 2-SIGMA UNCERTAINTY
            //F3493 !
            //F3494           IF(NEWGSIC.EQ.1)SLI(N)=GS1          ! NEW CODE
            if( ICE->NEWGSIC == 1 ) TANDSL->SLI[ N ] = GS1;
            //F3495           SLG(N)=ZGR-2.0*SQRT(DZGR1**2+DZGR2**2)
            TANDSL->SLG[ N ] = ZGR - 2.0 * sqrt( float( pow( DZGR1, 2) + pow( DZGR2, 2 ) ) );
            //F3496           SLA(N)=ZAN-2.0*DZAN
            TANDSL->SLA[ N ] = ZAN - 2.0 * DZAN;
            //F3497           SLO(N)=0.0
            TANDSL->SLO[ N ] = 0.0;
            //F3498           ZALL(N)=ZTOT-DZTOT
            ZALL[ N ] = ZTOT - DZTOT;
            //F3499 !
            //F3500           SLI1(N) =SLI(N)
            SLI1[ N ] = TANDSL->SLI[ N ];
            //F3501           SLG1(N) =SLG(N)
            SLG1[ N ] = TANDSL->SLG[ N ];
            //F3502           SLA1(N) =SLA(N)
            SLA1[ N ] = TANDSL->SLA[ N ];
            //F3503           SLO1(N) =SLO(N)
            SLO1[ N ] = TANDSL->SLO[ N ];
            //F3504           ZALL1(N)=ZALL(N)
            ZALL1[ N ] = ZALL[ N ];
            //F3505           EX1(N)=EX(N)
            EX1[ N ] = TANDSL->EX[ N ];
            //F3506         ENDIF
        }
        //F3507 !
        //F3508 !     NCLIM=2 = MID = MID EXPAN + MID MELT
        //F3509 !
        //F3510         IF(NSIM.EQ.2)THEN
        if( NSIM->NSIM == 2 ) {
            //F3511           SLI(N)=GS
            TANDSL->SLI[ N ] = GS;
            //F3512           IF(NEWGSIC.EQ.1)SLI(N)=GS2          ! NEW CODE
            if( ICE->NEWGSIC == 1 ) TANDSL->SLI[ N ] = GS2;
            //F3513           SLG(N)=ZGR
            TANDSL->SLG[ N ] = ZGR;
            //F3514           SLA(N)=ZAN
            TANDSL->SLA[ N ] = ZAN;
            //F3515           SLO(N)=ZI+ZP+ZS
            TANDSL->SLO[ N ] = ZI + ZP + ZS;
            //F3516           ZALL(N)=ZTOT
            ZALL[ N ] = ZTOT;
            //F3517 !
            //F3518           SLI2(N) =SLI(N)
            SLI2[ N ] = TANDSL->SLI[ N ];
            //F3519           SLG2(N) =SLG(N)
            SLG2[ N ] = TANDSL->SLG[ N ];
            //F3520           SLA2(N) =SLA(N)
            SLA2[ N ] = TANDSL->SLA[ N ];
            //F3521           SLO2(N) =SLO(N)
            SLO2[ N ] = TANDSL->SLO[ N ];
            //F3522           ZALL2(N)=ZALL(N)
            ZALL2[ N ] = ZALL[ N ];
            //F3523           EX2(N)=EX(N)
            EX2[ N ] = TANDSL->EX[ N ];
            //F3524         ENDIF
        }
        //F3525 !
        //F3526 !     NCLIM=3 = HIGH = HIGH EXPAN + HIGH MELT
        //F3527 !
        //F3528         IF(NSIM.EQ.3)THEN
        if( NSIM->NSIM == 3 ) {
            //F3529           SLI(N)=GS+DGS                       ! PREVIOUSLY +2.0*DGS (??)
            TANDSL->SLI[ N ] = GS + DGS;
            //F3530           IF(NEWGSIC.EQ.1)SLI(N)=GS3          ! NEW CODE
            if( ICE->NEWGSIC == 1 ) TANDSL->SLI[ N ] = GS3;
            //F3531           SLG(N)=ZGR+2.0*SQRT(DZGR1**2+DZGR2**2)
            TANDSL->SLG[ N ] = ZGR + 2.0 * sqrt( float( pow( DZGR1, 2 ) + pow( DZGR2, 2 ) ) );
            //F3532           SLA(N)=ZAN+2.0*DZAN
            TANDSL->SLA[ N ] = ZAN + 2.0 * DZAN;
            //F3533           SLO(N)=2.0*(ZI+ZP+ZS)
            TANDSL->SLO[ N ] = 2.0 * ( ZI + ZP + ZS );
            //F3534           ZALL(N)=ZTOT+DZTOT
            ZALL[ N ] = ZTOT + DZTOT;
            //F3535 !
            //F3536           SLI3(N) =SLI(N)
            SLI3[ N ] = TANDSL->SLI[ N ];
            //F3537           SLG3(N) =SLG(N)
            SLG3[ N ] = TANDSL->SLG[ N ];
            //F3538           SLA3(N) =SLA(N)
            SLA3[ N ] = TANDSL->SLA[ N ];
            //F3539           SLO3(N) =SLO(N)
            SLO3[ N ] = TANDSL->SLO[ N ];
            //F3540           ZALL3(N)=ZALL(N)
            ZALL3[ N ] = ZALL[ N ];
            //F3541           EX3(N)=EX(N)
            EX3[ N ] = TANDSL->EX[ N ];
            //F3542         ENDIF
        }
        //F3543 !
        //F3544 !  NCLIM=4 = USER = USER EXPAN + USER MELT
        //F3545 !   NOTE THAT NSIM=1 RESULTS ABOVE (LOW CASE) ARE OVER-WRITTEN
        //F3546 !    BY USER CASE IF ISCENGEN=9
        //F3547 !   PARAMETER 'ICE' DETERMINES WHETHER USER HAS CHOSEN TO USE
        //F3548 !    LOW, MID OR HIGH ICE MELT
        //F3549 !
        //F3550         IF(NSIM.EQ.4.OR.ISCENGEN.EQ.9)THEN
        if( NSIM->NSIM == 4 || NSIM->ISCENGEN == 9 ) {
            //F3551           SCALE=(ICE-2)*2.0
            const float SCALE = ( ICE->ICE - 2 ) * 2.0;
            //F3552 !
            //F3553           SLI(N)=GS+SCALE*DGS
            TANDSL->SLI[ N ] = GS + SCALE * DGS;
            //F3554           IF(NEWGSIC.EQ.1)SLI(N)=GS4          ! NEW CODE
            if( ICE->NEWGSIC == 1 ) TANDSL->SLI[ N ] = GS4;
            //F3555           SLG(N)=ZGR+SCALE*SQRT(DZGR1**2+DZGR2**2)
            TANDSL->SLG[ N ] = ZGR + SCALE * sqrt( float( pow( DZGR1, 2 ) + pow( DZGR2, 2 ) ) );
            //F3556           SLA(N)=ZAN+SCALE*DZAN
            TANDSL->SLA[ N ] = ZAN + SCALE * DZAN;
            //F3557           SLO(N)=(ICE-1)*(ZI+ZP+ZS)
            TANDSL->SLO[ N ] = ( ICE->ICE - 1 ) * ( ZI + ZP + ZS );
            //F3558           ZALL(N)=ZTOT+SCALE*DZTOT/2.0
            ZALL[ N ] = ZTOT + SCALE * DZTOT / 2.0;
            //F3559 !
            //F3560           SLI4(N) =SLI(N)
            SLI4[ N ] = TANDSL->SLI[ N ];
            //F3561           SLG4(N) =SLG(N)
            SLG4[ N ] = TANDSL->SLG[ N ];
            //F3562           SLA4(N) =SLA(N)
            SLA4[ N ] = TANDSL->SLA[ N ];
            //F3563           SLO4(N) =SLO(N)
            SLO4[ N ] = TANDSL->SLO[ N ];
            //F3564           ZALL4(N)=ZALL(N)
            ZALL4[ N ] = ZALL[ N ];
            //F3565           EX4(N)=EX(N)
            EX4[ N ] = TANDSL->EX[ N ];
            //F3566         ENDIF
        }
        //F3567 !
        //F3568 !  NOW TRANSFER APPROPRIATE MELT CASE RESULTS FOR NSIM VALUES
        //F3569 !   GREATER THAN 4.
        //F3570 !
        //F3571         IF(NSIM.EQ.5.OR.NSIM.EQ.9.OR.NSIM.EQ.13.OR.NSIM.EQ.17)THEN
        if( NSIM->NSIM == 5 || NSIM->NSIM == 9 || NSIM->NSIM == 13 || NSIM->NSIM == 17 ) {
            //F3572           SLI(N) =SLI1(N)
            TANDSL->SLI[ N ] = SLI1[ N ];
            //F3573           SLG(N) =SLG1(N)
            TANDSL->SLG[ N ] = SLG1[ N ];
            //F3574           SLA(N) =SLA1(N)
            TANDSL->SLA[ N ] = SLA1[ N ];
            //F3575           SLO(N) =SLO1(N)
            TANDSL->SLO[ N ] = SLO1[ N ];
            //F3576           ZALL(N)=ZALL1(N)
            ZALL[ N ] = ZALL1[ N ];
            //F3577           EX(N)=EX1(N)
            TANDSL->EX[ N ] = EX1[ N ];
            //F3578         ENDIF
        }
        //F3579 !
        //F3580         IF(NSIM.EQ.6.OR.NSIM.EQ.10.OR.NSIM.EQ.14.OR.NSIM.EQ.18)THEN
        if( NSIM->NSIM == 6 || NSIM->NSIM == 10 || NSIM->NSIM == 14 || NSIM->NSIM == 18 ) {
            //F3581           SLI(N) =SLI2(N)
            TANDSL->SLI[ N ] = SLI2[ N ];
            //F3582           SLG(N) =SLG2(N)
            TANDSL->SLG[ N ] = SLG2[ N ];
            //F3583           SLA(N) =SLA2(N)
            TANDSL->SLA[ N ] = SLA2[ N ];
            //F3584           SLO(N) =SLO2(N)
            TANDSL->SLO[ N ] = SLO2[ N ];
            //F3585           ZALL(N)=ZALL2(N)
            ZALL[ N ] = ZALL2[ N ];
            //F3586           EX(N)=EX2(N)
            TANDSL->EX[ N ] = EX2[ N ];
            //F3587         ENDIF
        }
        //F3588 !
        //F3589         IF(NSIM.EQ.7.OR.NSIM.EQ.11.OR.NSIM.EQ.15.OR.NSIM.EQ.19)THEN
        if( NSIM->NSIM == 7 || NSIM->NSIM == 11 || NSIM->NSIM == 15 || NSIM->NSIM == 19 ) {
            //F3590           SLI(N) =SLI3(N)
            TANDSL->SLI[ N ] = SLI3[ N ];
            //F3591           SLG(N) =SLG3(N)
            TANDSL->SLG[ N ] = SLG3[ N ];
            //F3592           SLA(N) =SLA3(N)
            TANDSL->SLA[ N ] = SLA3[ N ];
            //F3593           SLO(N) =SLO3(N)
            TANDSL->SLO[ N ] = SLO3[ N ];
            //F3594           ZALL(N)=ZALL3(N)
            ZALL[ N ] = ZALL3[ N ];
            //F3595           EX(N)=EX3(N)
            TANDSL->EX[ N ] = EX3[ N ];
            //F3596         ENDIF
        }
        //F3597 !
        //F3598         IF(NSIM.EQ.8.OR.NSIM.EQ.12.OR.NSIM.EQ.16.OR.NSIM.EQ.20)THEN
        if( NSIM->NSIM == 8 || NSIM->NSIM == 12 || NSIM->NSIM == 16 || NSIM->NSIM == 20 ) {
            //F3599           SLI(N) =SLI4(N)
            TANDSL->SLI[ N ] = SLI4[ N ];
            //F3600           SLG(N) =SLG4(N)
            TANDSL->SLG[ N ] = SLG4[ N ];
            //F3601           SLA(N) =SLA4(N)
            TANDSL->SLA[ N ] = SLA4[ N ];
            //F3602           SLO(N) =SLO4(N)
            TANDSL->SLO[ N ] = SLO4[ N ];
            //F3603           ZALL(N)=ZALL4(N)
            ZALL[ N ] = ZALL4[ N ];
            //F3604           EX(N)=EX4(N)
            TANDSL->EX[ N ] = EX4[ N ];
            //F3605         ENDIF
        }
        //F3606 !
        //F3607         SLT(N)=ZALL(N)+EX(N)
        TANDSL->SLT[ N ] = ZALL[ N ] + TANDSL->EX[ N ];
        //F3608       ENDIF
    }
    //F3609 !
    //F3610       RETURN
    //F3611       END
    f_exit( __func__ );
}
//F3612 !
//F3613 !  *******************************************************************
//F3614 !
//F3615       SUBROUTINE SPLIT(QGLOBE,A,BN,BS,QNO,QNL,QSO,QSL)
void split( const float QGLOBE, const float A, const float BN, const float BS, float* QNO, float* QNL, 
           float* QSO, float* QSL, AREAS_block* AREAS )
{
    // std::cout << "SUBROUTINE SPLIT" << endl;    
    f_enter( __func__ );
    //F3616       IMPLICIT REAL*4 (a-h,o-z), Integer (I-N)
    //F3617 !
    //F3618       COMMON /AREAS/FNO,FNL,FSO,FSL
    //F3619 !
    //F3620 !  Q VALUES ARE FORCINGS OVER AREAS IN W/M**2, F VALUES ARE THESE
    //F3621 !   MULTIPLIED BY AREA FRACTIONS. 
    //F3622 !   The resulting fractions times the appropriate area fractions
    //F3623 !   (FNO,FNL,FSO,FSL) sum to the global forcing. sjs
    //F3624 !
    //F3625 !  FIRST SPLIT QGLOBE INTO NH AND SH
    //F3626 !
    //F3627       Q=QGLOBE
    const float Q = QGLOBE;
    //F3628       QS=2.*Q/(A+1.)
    const float QS = 2.0 * Q / ( A + 1.0 );
    //F3629       QN=A*QS
    const float QN = A * QS;
    //F3630 !
    //F3631 !  NOW SPLIT NH AND SH INTO LAND AND OCEAN
    //F3632 !
    //F3633       FFN=2.0*(BN*FNL+FNO)
    const float FFN = 2.0 * ( BN * AREAS->FNL + AREAS->FNO );
    //F3634       QNO=QN/FFN
    *QNO = QN / FFN;
    //F3635       QNL=BN*QNO
    *QNL = BN * *QNO;
    //F3636 !
    //F3637       FFS=2.0*(BS*FSL+FSO)
    const float FFS = 2.0 * ( BS * AREAS->FSL + AREAS->FSO );
    //F3638       QSO=QS/FFS
    *QSO = QS / FFS;
    //F3639       QSL=BS*QSO
    *QSL = BS * *QSO;
    //F3640 !
    //F3641       RETURN
    //F3642       END
    f_exit( __func__ );
} // split
//F3643 !
//F3644 !  *******************************************************************
//F3645 !
//F3646       SUBROUTINE RUNMOD
void runmod( Limits_block* Limits, CLIM_block* CLIM, CONCS_block* CONCS, CARB_block* CARB, TANDSL_block* TANDSL,
            CO2READ_block* CO2READ, Sulph_block* Sulph, DSENS_block* DSENS, VARW_block* VARW, QSPLIT_block* QSPLIT,
            AREAS_block* AREAS, QADD_block* QADD, BCOC_block* BCOC, FORCE_block* FORCE, NSIM_block* NSIM,
            OZ_block* OZ, NEWCONCS_block* NEWCONCS, CAR_block* CAR, METH1_block* METH1, METH2_block* METH2, 
            METH3_block* METH3, METH4_block* METH4, TauNitr_block* TauNitr,
            JSTART_block* JSTART, CORREN_block* CORREN, HALOF_block* HALOF, COBS_block* COBS, ICE_block* ICE, std::ofstream* outfile8 )
{
    //    std::cout << "SUBROUTINE RUNMOD" << endl;    
    f_enter( __func__ );
    //F3647       IMPLICIT REAL*4 (a-h,o-z), Integer (I-N)
    //F3648 !
    //F3649       parameter (iTp=740)
    //F3650 !
    //F3651       DIMENSION AA(40),BB(40),A(40),B(40),C(40),D(40),XLLGLOBE(2), &
    //F3652       XLLDIFF(2)
    float AA[ 40+1 ], BB[ 40+1 ], A[ 40+1 ], B[ 40+1 ], C[ 40+1 ], D[ 40+1 ], XLLGLOBE[ 2+1 ], XLLDIFF[ 2+1 ];
    //F3653 !
    //F3654       common /Limits/KEND
    //F3655 !
    //F3656       COMMON/CLIM/IC,IP,KC,DT,DZ,FK,HM,Q2X,QXX,PI,T,TE,TEND,W0,XK,XKLO, &
    //F3657       XKNS,XLAM,FL(2),FO(2),FLSUM,FOSUM,HEM(2),P(40),TEM(40),TO(2,40), &
    //F3658       AL,BL,CL,DTH,DTZ,DZ1,XLL,WWW,XXX,YYY,RHO,SPECHT,HTCONS,Y(4)
    //F3659 !
    //F3660       COMMON/CONCS/CH4(0:iTp),CN2O(0:iTp),ECH4(226:iTp+1), &
    //F3661       EN2O(226:iTp+1),ECO(226:iTp+1),COE(iTp+1),EVOC(226:iTp+1), &
    //F3662       ENOX(226:iTp+1),ESO2(0:iTp+1),ESO2SUM(226:iTp+1), &
    //F3663       ESO21(226:iTp+1),ESO22(226:iTp+1),ESO23(226:iTp+1), &
    //F3664       EBC(226:iTp+1), EOC(226:iTp+1) ! sjs- add BC-OC
    //F3665 !
    //F3666       COMMON/CARB/CCO2(4,224:iTp),EDGROSS(4,226:iTp),EF(226:iTp+1), &
    //F3667       REGROW(4,226:iTp),PL(4,226:iTp),HL(4,226:iTp),SOIL(4,226:iTp), &
    //F3668       TTT(226:iTp),ESUM(226:iTp),ETOT(4,226:iTp),EDNET90(4), &
    //F3669       FOC(4,226:iTp),CO2(0:iTp),CO2SAVE(0:iTp)
    //F3670 !
    //F3671       COMMON/TANDSL/TEQU(iTp),TGAV(iTp),TNHO(iTp), &
    //F3672       TSHO(iTp),TNHL(iTp),TSHL(iTp),TDEEP(iTp),TNHAV(iTp),TSHAV(iTp), &
    //F3673       TLAND(iTp),TOCEAN(iTp),TOCN(40),TOCNPREV(40), &
    //F3674       SIP,SGP,SAP,SLI(iTp),SLG(iTp),SLA(iTp),EX(0:iTp),SLT(iTp), &
    //F3675       QTOT(0:iTp),QGH(0:iTp),QOZ(0:iTp),QBIO(0:iTp),SLO(iTp), &
    //F3676       QSO2(0:iTp+1),QDIR(0:iTp+1),QLAND(0:iTp),QMN(0:iTp+1)
    //F3677 !
    //F3678       COMMON /CO2READ/ICO2READ,XC(226:iTp),CO2SCALE,qtot86,LEVCO2
    //F3679 !
    //F3680       common /Sulph/S90DIR,S90IND,S90BIO,ENAT,ES1990,ECO90,FOC90,IFOC
    //F3681       COMMON /DSENS/IXLAM,XLAML,XLAMO,ADJUST
    //F3682       COMMON /VARW/Z(40),W(2),DW(2),TO0(2),TP0(2),WNH(iTp),WSH(iTp), &
    //F3683       TW0NH,TW0SH,IVARW,KEYDW
    //F3684       COMMON /QSPLIT/QNHO,QNHL,QSHO,QSHL,QGLOBE(0:iTp), &
    //F3685       QQNHO(0:iTp),QQNHL(0:iTp),QQSHO(0:iTp),QQSHL(0:iTp), &
    //F3686       QQQNHO(0:iTp),QQQNHL(0:iTp),QQQSHO(0:iTp),QQQSHL(0:iTp), &
    //F3687       EHistBC(iTp),EHistOC(iTp) ! Vars to store read-in BC-OC history.
    //F3688 
    //F3689       COMMON /AREAS/FNO,FNL,FSO,FSL
    //F3690 !
    //F3691       COMMON /QADD/IQREAD,OrgIQREAD,JQFIRST,JQLAST,QEX(0:iTp),QEXNH(0:iTp), &
    //F3692       QEXSH(0:iTp),QEXNHO(0:iTp),QEXNHL(0:iTp),QEXSHO(0:iTp), &
    //F3693       QEXSHL(0:iTp),IOLDTZ
    //F3694 ! BCOC params to set read-in forcing
    //F3695       REAL*4 FBC1990, FOC1990, FSO2_dir1990,FSO2_ind1990, aBCUnitForcing, aOCUnitForcing, &
    //F3696              aBCBaseEmissions, aOCBaseEmissions !sjs
    //F3697       COMMON/BCOC/FBC1990, FOC1990, FSO2_dir1990,FSO2_ind1990, aBCUnitForcing, aOCUnitForcing, &
    //F3698              aBCBaseEmissions, aOCBaseEmissions !sjs
    //F3699       COMMON /FORCE/qco2(0:iTp),qm(0:iTp),qn(0:iTp),QCFC(0:iTp), &
    //F3700       QMONT(0:iTp),QOTHER(0:iTp),QSTRATOZ(0:iTp),QCH4O3(0:iTp), &
    //F3701       CFC12(0:iTp), QCH4H2O(0:iTp),QBC(0:iTp),QOC(0:iTp)
    //F3702 !
    //F3703       COMMON /NSIM/NSIM,NCLIM,ISCENGEN,TEMEXP(2,40),IWNHOFF,IWSHOFF, &
    //F3704       WTHRESH
    //F3705 !
    //F3706   11  CONTINUE
    do {    
        //F3707 !
        //F3708 !  INCREMENT COUNTER AND ADD TIME STEP.
        //F3709 !   T CORRESPONDS TO THE TIME AT WHICH NEW VALUES ARE TO
        //F3710 !   BE CALCULATED
        //F3711 !
        //F3712       IP=IC
        CLIM->IP = CLIM->IC;
        //F3713       T=T+DT
        CLIM->T += CLIM->DT;
        //F3714 !
        //F3715       IC=INT(T+1.49)
        CLIM->IC = int( CLIM->T + 1.49 );
        //F3716 !
        //F3717 !  PROTECT AGAINST CUMULATIVE ROUND-OFF ERROR IN T.
        //F3718 !   (IF T IS VERY NEAR AN INTEGER, ROUND UP OR DOWN TO INTEGER.)
        //F3719 !
        //F3720       DIFF=ABS(T-INT(T))
        float XT, DIFF = fabs( CLIM->T - int( CLIM->T ) );
        //F3721       IF(DIFF.LE.0.5)THEN
        if( DIFF <= 0.5 ) {
            //F3722         XT=FLOAT(INT(T))
            XT = float( int( CLIM->T ) );
            //F3723       ELSE
        } else {
            //F3724         DIFF=1.0-DIFF
            DIFF = 1.0 - DIFF;
            //F3725         XT=FLOAT(INT(T)+1)
            XT = float( int( CLIM->T ) +1 );
            //F3726       ENDIF
        }
        //F3727       IF(DIFF.LT.0.01)T=XT
        if( DIFF < 0.01 ) CLIM->T = XT;
        //F3728 !
        //F3729 !  AS SOON AS A NEW YEAR IS ENCOUNTERED (I.E. IC IS
        //F3730 !   INCREMENTED), CALCULATE NEW END OF YEAR VALUES FOR
        //F3731 !   CONCENTRATIONS AND RADIATIVE FORCING.
        //F3732 !  THIS IS ONLY NECESSARY FOR FIRST PASS THROUGH NUMSIMS LOOP
        //F3733 !   WHICH WILL SET ALL VALUES OF FORCING COMPONENT ARRAYS.
        //F3734 !
        //F3735 !      IF(NCLIM.EQ.1.OR.ISCENGEN.EQ.9)THEN
        //F3736 !
        //F3737        IF(IC.GT.IP) CALL DELTAQ
        if( CLIM->IC > CLIM->IP ) deltaq( Limits, OZ, CLIM, CONCS,
                                         NEWCONCS, CARB, TANDSL, CAR,
                                         METH1, FORCE, METH2, METH3,
                                         METH4, TauNitr, Sulph, NSIM, 
                                         CO2READ, JSTART, CORREN, HALOF, COBS );
        //F3738 !
        //F3739 !      ENDIF
        //F3740 !
        //F3741 !  INTERPOLATE FORCING FROM VALUES AT ENDS OF YEARS TO
        //F3742 !   VALUE AT MIDPOINT OF TIME STEP.
        //F3743 !
        //F3744 !  ***********************************************************
        //F3745 !
        //F3746 !  REVISED ALGORITHM FOR INTERPOLATION (15 APR, 1994)
        //F3747 !
        //F3748 !  FIRST CALCULATE FRACTION OF YEAR THAT TIME CORRESPONDS TO AND
        //F3749 !   IDENTIFY THE INTEGER (JC) THAT IS THE END OF THE YEAR IN
        //F3750 !   WHICH THE TIME LIES.
        //F3751 !
        //F3752       T1=T-DT/2.0
        const float T1 = CLIM->T - CLIM->DT / 2.0;
        //F3753       JC=INT(T1+1.49)
        const int JC = int( T1 + 1.49 );
        //F3754       FRAC=T1+0.5-INT(T1+0.5)
        float FRAC = T1 + 0.5 - int( T1 + 0.5 );
        //F3755       IF(FRAC.LT.0.001)FRAC=1.0
        if( FRAC < 0.001 ) FRAC = 1.0;
        //F3756 !
        //F3757 !  CALCULATE GREENHOUSE GAS AND BIOMASS AEROSOL COMPONENTS OF
        //F3758 !   FORCING AT START (0) AND END (1) OF YEAR, AND AT START OF
        //F3759 !   PREVIOUS YEAR (P).
        //F3760 !
        //F3761       JPREV=0
        int JPREV = 0;
        //F3762       IF(JC.GE.2)JPREV=JC-2
        if( JC >= 2 ) JPREV = JC - 2;
        //F3763 !
        //F3764       QGHP  = QGH(JPREV)
        float QGHP = TANDSL->QGH[ JPREV ];
        //F3765       QGH0  = QGH(JC-1)
        float QGH0 = TANDSL->QGH[ JC-1 ];
        //F3766       QGH1  = QGH(JC)
        float QGH1 = TANDSL->QGH[ JC ];
        //F3767 !
        //F3768       QBIOP = QBIO(JPREV)
        float QBIOP = TANDSL->QBIO[ JPREV ];
        //F3769       QBIO0 = QBIO(JC-1)
        float QBIO0 = TANDSL->QBIO[ JC-1 ];
        //F3770       QBIO1 = QBIO(JC)
        float QBIO1 = TANDSL->QBIO[ JC ];
        //F3771 !
        //F3772       QLANDP= QLAND(JPREV)
        float QLANDP = TANDSL->QLAND[ JPREV ];
        //F3773       QLAND0= QLAND(JC-1)
        float QLAND0 = TANDSL->QLAND[ JC-1 ];
        //F3774       QLAND1= QLAND(JC)
        float QLAND1 = TANDSL->QLAND[ JC ];
        //F3775 !
        //F3776       QMNP  = QMN(JPREV)
        float QMNP = TANDSL->QMN[ JPREV ];
        //F3777       QMN0  = QMN(JC-1)
        float QMN0 = TANDSL->QMN[ JC-1 ];
        //F3778       QMN1  = QMN(JC)
        float QMN1 = TANDSL->QMN[ JC ];
        //F3779 !
        //F3780 !  RELABEL DIRECT AEROSOL AND OZONE FORCINGS
        //F3781 !
        //F3782       QDIRP = QDIR(JPREV)
        const float QDIRP = TANDSL->QDIR[ JPREV ];
        //F3783       QDIR0 = QDIR(JC-1)
        const float QDIR0 = TANDSL->QDIR[ JC-1 ];
        //F3784       QDIR1 = QDIR(JC)
        const float QDIR1 = TANDSL->QDIR[ JC ];
        //F3785       QOZP  = QOZ(JPREV)
        const float QOZP = TANDSL->QOZ[ JPREV ];
        //F3786       QOZ0  = QOZ(JC-1)
        const float QOZ0 = TANDSL->QOZ[ JC-1 ];
        //F3787       QOZ1  = QOZ(JC)
        const float QOZ1 = TANDSL->QOZ[ JC ];
        //F3788 !
        //F3789 !  CALCULATE INDIRECT AEROSOL COMPONENT.
        //F3790 !
        //F3791       QINDP  =QSO2(JPREV)-QDIRP
        const float QINDP = TANDSL->QSO2[ JPREV ] - QDIRP;
        //F3792       QIND0  =QSO2(JC-1)-QDIR0
        const float QIND0 = TANDSL->QSO2[ JC-1 ] - QDIR0;
        //F3793       QIND1  =QSO2(JC)-QDIR1
        const float QIND1 = TANDSL->QSO2[ JC ] - QDIR1;
        //F3794 !
        //F3795 !  ***********************************************************
        //F3796 !
        //F3797 !  SPLIT AEROSOL & TROP O3 FORCING INTO LAND AND OCEAN IN NH AND SH
        //F3798 !
        //F3799 !  A IS THE NH/SH FORCING RATIO
        //F3800 !  BN IS THE LAND/OCEAN FORCING RATIO IN THE NH
        //F3801 !  BS IS THE LAND/OCEAN FORCING RATIO IN THE SH
        //F3802 !
        //F3803       ADIR = 4.0
        const float ADIR = 4.0;
        //F3804       AIND = 2.0
        const float AIND = 2.0;
        //F3805       AOZ  =99.0
        const float AOZ = 99.0;
        //F3806       BNDIR= 9.0
        const float BNDIR = 9.0;
        //F3807       BNIND= 9.0
        const float BNIND = 9.0;
        //F3808       BNOZ = 9.0
        const float BNOZ = 9.0;
        //F3809       BSDIR= 9.0
        const float BSDIR = 9.0;
        //F3810       BSIND= 9.0
        const float BSIND = 9.0;
        //F3811       BSOZ = 9.0
        const float BSOZ = 9.0;
        //F3812 !
        //F3813       CALL SPLIT(QINDP,AIND,BNIND,BSIND,QINDNOP,QINDNLP, &
        //F3814       QINDSOP,QINDSLP)
        float QINDNOP = 0.0, QINDNLP = 0.0, QINDSOP = 0.0, QINDSLP = 0.0;
        split( QINDP, AIND, BNIND, BSIND, &QINDNOP, &QINDNLP, &QINDSOP, &QINDSLP, AREAS );
        //F3815       CALL SPLIT(QIND0,AIND,BNIND,BSIND,QINDNO0,QINDNL0, &
        //F3816       QINDSO0,QINDSL0)
        float QINDNO0 = 0.0, QINDNL0 = 0.0, QINDSO0 = 0.0, QINDSL0 = 0.0;
        split( QIND0, AIND, BNIND, BSIND, &QINDNO0, &QINDNL0, &QINDSO0, &QINDSL0, AREAS );
        //F3817       CALL SPLIT(QIND1,AIND,BNIND,BSIND,QINDNO1,QINDNL1, &
        //F3818       QINDSO1,QINDSL1)
        float QINDNO1 = 0.0, QINDNL1 = 0.0, QINDSO1 = 0.0, QINDSL1 = 0.0;
        split( QIND1, AIND, BNIND, BSIND, &QINDNO1, &QINDNL1, &QINDSO1, &QINDSL1, AREAS );
        //F3819 !
        //F3820       CALL SPLIT(QDIRP,ADIR,BNDIR,BSDIR,QDIRNOP,QDIRNLP, &
        //F3821       QDIRSOP,QDIRSLP)
        float QDIRNOP = 0.0, QDIRNLP = 0.0, QDIRSOP = 0.0, QDIRSLP = 0.0;
        split( QDIRP, ADIR, BNDIR, BSDIR, &QDIRNOP, &QDIRNLP, &QDIRSOP, &QDIRSLP, AREAS );
        //F3822       CALL SPLIT(QDIR0,ADIR,BNDIR,BSDIR,QDIRNO0,QDIRNL0, &
        //F3823       QDIRSO0,QDIRSL0)
        float QDIRNO0 = 0.0, QDIRNL0 = 0.0, QDIRSO0 = 0.0, QDIRSL0 = 0.0;
        split( QDIR0, ADIR, BNDIR, BSDIR, &QDIRNO0, &QDIRNL0, &QDIRSO0, &QDIRSL0, AREAS );
        //F3824       CALL SPLIT(QDIR1,ADIR,BNDIR,BSDIR,QDIRNO1,QDIRNL1, &
        //F3825      QDIRSO1,QDIRSL1)
        float QDIRNO1 = 0.0, QDIRNL1 = 0.0, QDIRSO1 = 0.0, QDIRSL1 = 0.0;
        split( QDIR1, ADIR, BNDIR, BSDIR, &QDIRNO1, &QDIRNL1, &QDIRSO1, &QDIRSL1, AREAS );
        //F3826 !
        //F3827       CALL SPLIT(QOZP,AOZ,BNOZ,BSOZ,QOZNOP,QOZNLP, &
        //F3828       QOZSOP,QOZSLP)
        float QOZNOP = 0.0, QOZNLP = 0.0, QOZSOP = 0.0, QOZSLP = 0.0;
        split( QOZP, AOZ, BNOZ, BSOZ, &QOZNOP, &QOZNLP, &QOZSOP, &QOZSLP, AREAS );
        //F3829       CALL SPLIT(QOZ0,AOZ,BNOZ,BSOZ,QOZNO0,QOZNL0, &
        //F3830       QOZSO0,QOZSL0)
        float QOZNO0 = 0.0, QOZNL0 = 0.0, QOZSO0 = 0.0, QOZSL0 = 0.0;
        split( QOZ0, AOZ, BNOZ, BSOZ, &QOZNO0, &QOZNL0, &QOZSO0, &QOZSL0, AREAS );
        //F3831       CALL SPLIT(QOZ1,AOZ,BNOZ,BSOZ,QOZNO1,QOZNL1, &
        //F3832       QOZSO1,QOZSL1)
        float QOZNO1 = 0.0, QOZNL1 = 0.0, QOZSO1 = 0.0, QOZSL1 = 0.0;
        split( QOZ1, AOZ, BNOZ, BSOZ, &QOZNO1, &QOZNL1, &QOZSO1, &QOZSL1, AREAS );
        //F3833 !
        //F3834       QSNHOP =QDIRNOP+QINDNOP+QOZNOP
        float QSNHOP = QDIRNOP + QINDNOP + QOZNOP;
        //F3835       QSNHLP =QDIRNLP+QINDNLP+QOZNLP
        float QSNHLP = QDIRNLP + QINDNLP + QOZNLP;
        //F3836       QSSHOP =QDIRSOP+QINDSOP+QOZSOP
        float QSSHOP = QDIRSOP + QINDSOP + QOZSOP;
        //F3837       QSSHLP =QDIRSLP+QINDSLP+QOZSLP
        float QSSHLP = QDIRSLP + QINDSLP + QOZSLP;
        //F3838 !
        //F3839       QSNHO0 =QDIRNO0+QINDNO0+QOZNO0
        float QSNHO0 = QDIRNO0 + QINDNO0 + QOZNO0;
        //F3840       QSNHL0 =QDIRNL0+QINDNL0+QOZNL0
        float QSNHL0 = QDIRNL0 + QINDNL0 + QOZNL0;
        //F3841       QSSHO0 =QDIRSO0+QINDSO0+QOZSO0
        float QSSHO0 = QDIRSO0 + QINDSO0 + QOZSO0;
        //F3842       QSSHL0 =QDIRSL0+QINDSL0+QOZSL0
        float QSSHL0 = QDIRSL0 + QINDSL0 + QOZSL0;
        //F3843 !
        //F3844       QSNHO1 =QDIRNO1+QINDNO1+QOZNO1
        float QSNHO1 = QDIRNO1 + QINDNO1 + QOZNO1;
        //F3845       QSNHL1 =QDIRNL1+QINDNL1+QOZNL1
        float QSNHL1 = QDIRNL1 + QINDNL1 + QOZNL1;
        //F3846       QSSHO1 =QDIRSO1+QINDSO1+QOZSO1
        float QSSHO1 = QDIRSO1 + QINDSO1 + QOZSO1;
        //F3847       QSSHL1 =QDIRSL1+QINDSL1+QOZSL1
        float QSSHL1 = QDIRSL1 + QINDSL1 + QOZSL1;
        //F3848 !
        //F3849 !  ***********************************************************
        //F3850 !
        //F3851 !  IF EXTRA FORCINGS ADDED THROUGH QEXTRA.IN, CALC CORRESP
        //F3852 !   'P', '0' AND '1' COMPONENTS FOR THESE.  NOTE THAT THESE
        //F3853 !   DATA ARE INPUT AS ANNUAL-MEAN VALUES, SO THEY MUST
        //F3854 !   BE APPLIED EQUALLY AT THE START AND END OF THE YEAR,
        //F3855 !   I.E., Q0=Q1=Q(JC).
        //F3856 !
        //F3857       QEXNHOP=0.0
        //F3858       QEXSHOP=0.0
        //F3859       QEXNHLP=0.0
        //F3860       QEXSHLP=0.0
        //F3861       QEXNHO0=0.0
        //F3862       QEXSHO0=0.0
        //F3863       QEXNHL0=0.0
        //F3864       QEXSHL0=0.0
        //F3865       QEXNHO1=0.0
        //F3866       QEXSHO1=0.0
        //F3867       QEXNHL1=0.0
        //F3868       QEXSHL1=0.0
        float QEXNHOP, QEXSHOP, QEXNHLP, QEXSHLP, QEXNHO0, QEXSHO0, QEXNHL0, QEXSHL0, QEXNHO1, QEXSHO1, QEXNHL1, QEXSHL1;
        QEXNHOP = QEXSHOP = QEXNHLP = QEXSHLP = QEXNHO0 = QEXSHO0 = QEXNHL0 = QEXSHL0 = QEXNHO1 = QEXSHO1 = QEXNHL1 = QEXSHL1 = 0.0;
        //F3869 !
        //F3870       IF(IQREAD.GE.1)THEN
        if( QADD->IQREAD >= 1 ) {
            //F3871         IF((JC.GE.JQFIRST).AND.(JC.LE.JQLAST))THEN
            if( JC >= QADD->JQFIRST && JC <= QADD->JQLAST ) {
                //F3872           QEXNHOP=QEXNHO(JC-1)
                QEXNHOP = QADD->QEXNHO[ JC-1 ];
                //F3873           QEXSHOP=QEXSHO(JC-1)
                QEXSHOP = QADD->QEXSHO[ JC-1 ];
                //F3874           QEXNHLP=QEXNHL(JC-1)
                QEXNHLP = QADD->QEXNHL[ JC-1 ];
                //F3875           QEXSHLP=QEXSHL(JC-1)
                QEXSHLP = QADD->QEXSHL[ JC-1 ];
                //F3876           QEXNHO0=QEXNHO(JC)
                QEXNHO0 = QADD->QEXNHO[ JC ];
                //F3877           QEXSHO0=QEXSHO(JC)
                QEXSHO0 = QADD->QEXSHO[ JC ];
                //F3878           QEXNHL0=QEXNHL(JC)
                QEXNHL0 = QADD->QEXNHL[ JC ];
                //F3879           QEXSHL0=QEXSHL(JC)
                QEXSHL0 = QADD->QEXSHL[ JC ];
                //F3880           QEXNHO1=QEXNHO(JC)
                QEXNHO1 = QADD->QEXNHO[ JC ];
                //F3881           QEXSHO1=QEXSHO(JC)
                QEXSHO1 = QADD->QEXSHO[ JC ];
                //F3882           QEXNHL1=QEXNHL(JC)
                QEXNHL1 = QADD->QEXNHL[ JC ];
                //F3883           QEXSHL1=QEXSHL(JC)
                QEXSHL1 = QADD->QEXSHL[ JC ];
                //F3884         ENDIF
            }
            //F3885       ENDIF
        }
        //F3886 
        //F3887 ! If read-in, calc explicit BCOC forcing, split into hemispheres and add to QEXTRA here
        //F3888      IF( IFOC.GE.3 .AND. IQREAD.NE.2 )THEN
        if( Sulph->IFOC >= 3 && QADD->IQREAD != 2 ) {
            //F3889         ! Save total BC and OC Forcing
            //F3890         QBC(JC) = EHistBC(JC) * aBCUnitForcing
            FORCE->QBC[ JC ] = QSPLIT->EHistBC[ JC ] * BCOC->aBCUnitForcing;
            //F3891         QOC(JC) = EHistOC(JC) * aOCUnitForcing
            FORCE->QOC[ JC ] = QSPLIT->EHistOC[ JC ] * BCOC->aOCUnitForcing;
            //F3892 
            //F3893 !  A IS THE NH/SH FORCING RATIO
            //F3894 !  BN IS THE LAND/OCEAN FORCING RATIO IN THE NH
            //F3895 !  BS IS THE LAND/OCEAN FORCING RATIO IN THE SH
            //F3896 !
            //F3897       A_BCOC = 8.0
            const float A_BCOC = 8.0;
            //F3898       BN_BCOC = 9.0
            const float BN_BCOC = 9.0;
            //F3899       BS_BCOC = 9.0
            const float BS_BCOC = 9.0;
            //F3900  
            //F3901         QBCOCP = ( EHistBC(JPREV) * aBCUnitForcing + EHistOC(JPREV) * aOCUnitForcing )
            //UNUSED const float QBCOCP = ( QSPLIT->EHistBC[ JPREV ] * BCOC->aBCUnitForcing + QSPLIT->EHistOC[ JPREV ] * BCOC->aOCUnitForcing );
            //F3902         QBCOC0 = ( EHistBC(JC-1) * aBCUnitForcing + EHistOC(JC-1) * aOCUnitForcing )
            //UNUSED const float QBCOC0 = ( QSPLIT->EHistBC[ JC-1 ] * BCOC->aBCUnitForcing + QSPLIT->EHistOC[ JC-1 ] * BCOC->aOCUnitForcing );
            //F3903         QBCOC1 = QBC(JC) + QOC(JC)
            float QBCOC1 = FORCE->QBC[ JC ] + FORCE->QOC[ JC ];
            //F3904         QBCOCP = QBC(JPREV) + QOC(JPREV)
            float QBCOCP = FORCE->QBC[ JPREV ] + FORCE->QOC[ JPREV ];
            //F3905         QBCOC0 = QBC(JC-1) + QOC(JC-1)
            float QBCOC0 = FORCE->QBC[ JC-1 ] + FORCE->QOC[ JC-1 ];
            //F3906         
            //F3907         ! If 1990 or later, then use values from input file
            //F3908         ! Input values are in units of Tg
            //F3909         IF ( JC .GT. 225 ) THEN
            if( JC > 225 ) {
                //F3910             QBC(JC) = EBC(JC) * aBCUnitForcing
                FORCE->QBC[ JC ] = CONCS->EBC.getval( JC ) * BCOC->aBCUnitForcing;
                //F3911             QOC(JC) = EOC(JC) * aOCUnitForcing
                FORCE->QOC[ JC ] = CONCS->EOC.getval( JC ) * BCOC->aOCUnitForcing;
                //F3912             
                //F3913             QBCOCP = QBC(JPREV) + QOC(JPREV)
                QBCOCP = FORCE->QBC[ JPREV ] + FORCE->QOC[ JPREV ];
                //F3914             QBCOC0 = QBC(JC-1) + QOC(JC-1)
                QBCOC0 = FORCE->QBC[ JC-1 ] + FORCE->QOC[ JC-1 ];
                //F3915             QBCOC1 = QBC(JC) + QOC(JC)
                QBCOC1 = FORCE->QBC[ JC ] + FORCE->QOC[ JC ];
                //F3916         ENDIF
            }
            //F3917 
            //F3918 		! Split into N and S Hem and ocean same as SO2 Dir
            //F3919         CALL SPLIT(QBCOCP,A_BCOC,BN_BCOC,BS_BCOC,QBCOCNOP,QBCOCNLP,QBCOCSOP,QBCOCSLP)
            float QBCOCNOP = 0.0, QBCOCNLP = 0.0, QBCOCSOP = 0.0,QBCOCSLP = 0.0;
            split( QBCOCP, A_BCOC, BN_BCOC, BS_BCOC, &QBCOCNOP, &QBCOCNLP, &QBCOCSOP, &QBCOCSLP, AREAS );
            //F3920         CALL SPLIT(QBCOC0,A_BCOC,BN_BCOC,BS_BCOC,QBCOCNO0,QBCOCNL0,QBCOCSO0,QBCOCSL0)
            float QBCOCNO0 = 0.0, QBCOCNL0 = 0.0, QBCOCSO0 = 0.0,QBCOCSL0 = 0.0;
            split( QBCOC0, A_BCOC, BN_BCOC, BS_BCOC, &QBCOCNO0, &QBCOCNL0, &QBCOCSO0, &QBCOCSL0, AREAS );
            //F3921         CALL SPLIT(QBCOC1,A_BCOC,BN_BCOC,BS_BCOC,QBCOCNO1,QBCOCNL1,QBCOCSO1,QBCOCSL1)
            float QBCOCNO1 = 0.0, QBCOCNL1 = 0.0, QBCOCSO1 = 0.0,QBCOCSL1 = 0.0;
            split( QBCOC1, A_BCOC, BN_BCOC, BS_BCOC, &QBCOCNO1, &QBCOCNL1, &QBCOCSO1, &QBCOCSL1, AREAS );
            //F3922 
            //F3923         QEXNHOP = QEXNHOP + QBCOCNOP
            QEXNHOP += QBCOCNOP;
            //F3924         QEXNHLP = QEXNHLP + QBCOCNLP
            QEXNHLP += QBCOCNLP;
            //F3925         QEXSHOP = QEXSHOP + QBCOCSOP
            QEXSHOP += QBCOCSOP;
            //F3926         QEXSHLP = QEXSHLP + QBCOCSLP
            QEXSHLP += QBCOCSLP;
            //F3927         QEXNHO0 = QEXNHO0 + QBCOCNO0
            QEXNHO0 += QBCOCNO0;
            //F3928         QEXNHL0 = QEXNHL0 + QBCOCNL0
            QEXNHL0 += QBCOCNL0;
            //F3929         QEXSHO0 = QEXSHO0 + QBCOCSO0
            QEXSHO0 += QBCOCSO0;
            //F3930         QEXSHL0 = QEXSHL0 + QBCOCSL0
            QEXSHL0 += QBCOCSL0;
            //F3931         QEXNHO1 = QEXNHO1 + QBCOCNO1
            QEXNHO1 += QBCOCNO1;
            //F3932         QEXNHL1 = QEXNHL1 + QBCOCNL1
            QEXNHL1 += QBCOCNL1;
            //F3933         QEXSHO1 = QEXSHO1 + QBCOCSO1
            QEXSHO1 += QBCOCSO1;
            //F3934         QEXSHL1 = QEXSHL1 + QBCOCSL1
            QEXSHL1 += QBCOCSL1;
            //F3935     
            //F3936       ENDIF
        }
        //F3937 
        //F3938 !
        //F3939 !  IF EXTRA NH AND SH FORCING INPUT THROUGH QEXTRA.IN, ADD
        //F3940 !   TO AEROSOL FORCING
        //F3941 !  IF IQREAD=2, USE EXTRA FORCING ONLY
        //F3942 !
        //F3943       IQR=1
        int IQR = 1;
        //F3944       IF(IQREAD.EQ.2)IQR=0
        if( QADD->IQREAD == 2 ) IQR = 0;
        //F3945 !
        //F3946       QSNHOP=IQR*QSNHOP+QEXNHOP
        QSNHOP = IQR * QSNHOP + QEXNHOP;
        //F3947       QSSHOP=IQR*QSSHOP+QEXSHOP
        QSSHOP = IQR * QSSHOP + QEXSHOP;
        //F3948       QSNHLP=IQR*QSNHLP+QEXNHLP
        QSNHLP = IQR * QSNHLP + QEXNHLP;
        //F3949       QSSHLP=IQR*QSSHLP+QEXSHLP
        QSSHLP=IQR * QSSHLP + QEXSHLP;
        //F3950       QGHP  =IQR*QGHP
        QGHP = IQR * QGHP;
        //F3951       QBIOP =IQR*QBIOP
        QBIOP = IQR * QBIOP;
        //F3952       QLANDP=IQR*QLANDP
        QLANDP = IQR * QLANDP;
        //F3953       QMNP  =IQR*QMNP
        QMNP = IQR * QMNP;
        //F3954 !
        //F3955       QSNHO0=IQR*QSNHO0+QEXNHO0
        QSNHO0 = IQR * QSNHO0 + QEXNHO0;
        //F3956       QSSHO0=IQR*QSSHO0+QEXSHO0
        QSSHO0 = IQR * QSSHO0 + QEXSHO0;
        //F3957       QSNHL0=IQR*QSNHL0+QEXNHL0
        QSNHL0 = IQR * QSNHL0 + QEXNHL0;
        //F3958       QSSHL0=IQR*QSSHL0+QEXSHL0
        QSSHL0 = IQR * QSSHL0 + QEXSHL0;
        //F3959       QGH0  =IQR*QGH0
        QGH0 = IQR * QGH0;
        //F3960       QBIO0 =IQR*QBIO0
        QBIO0 = IQR * QBIO0;
        //F3961       QLAND0=IQR*QLAND0
        QLAND0=IQR * QLAND0;
        //F3962       QMN0  =IQR*QMN0
        QMN0 =IQR * QMN0;
        //F3963 !
        //F3964       QSNHO1=IQR*QSNHO1+QEXNHO1
        QSNHO1 = IQR * QSNHO1 + QEXNHO1;
        //F3965       QSSHO1=IQR*QSSHO1+QEXSHO1
        QSSHO1 = IQR * QSSHO1 + QEXSHO1;
        //F3966       QSNHL1=IQR*QSNHL1+QEXNHL1
        QSNHL1 = IQR * QSNHL1 + QEXNHL1;
        //F3967       QSSHL1=IQR*QSSHL1+QEXSHL1
        QSSHL1 = IQR * QSSHL1 + QEXSHL1;
        //F3968       QGH1  =IQR*QGH1
        QGH1 = IQR * QGH1;
        //F3969       QBIO1 =IQR*QBIO1
        QBIO1 = IQR * QBIO1;
        //F3970       QLAND1=IQR*QLAND1
        QLAND1 = IQR * QLAND1;
        //F3971       QMN1  =IQR*QMN1
        QMN1 = IQR * QMN1;
        //F3972 !
        //F3973 !  CALCULATE FORCING COMPONENTS FOR MIDPOINT OF INTERVAL.
        //F3974 !
        //F3975       QSNHLM =(QSNHL0+QSNHL1)/2.0
        const float QSNHLM = ( QSNHL0 + QSNHL1 ) / 2.0;
        //F3976       QSSHLM =(QSSHL0+QSSHL1)/2.0
        const float QSSHLM = ( QSSHL0 + QSSHL1 ) / 2.0;
        //F3977       QSNHOM =(QSNHO0+QSNHO1)/2.0
        const float QSNHOM = ( QSNHO0 + QSNHO1 ) / 2.0;
        //F3978       QSSHOM =(QSSHO0+QSSHO1)/2.0
        const float QSSHOM = ( QSSHO0 + QSSHO1 ) / 2.0;
        //F3979       QGHM   =(QGH0  +QGH1  )/2.0
        const float QGHM = ( QGH0 + QGH1 ) / 2.0;
        //if (JC>235) cout << "QGHM " << JC << " " <<  QGHM << " " << QGH0 << " " << QGH1 << " " << IQR << endl; //debug
        
        //F3980       QBIOM  =(QBIO0 +QBIO1 )/2.0
        const float QBIOM = ( QBIO0 + QBIO1 ) / 2.0;
        //F3981       QLANDM =(QLAND0+QLAND1)/2.0
        const float QLANDM = ( QLAND0 + QLAND1 ) / 2.0;
        //F3982       QMNM   =(QMN0+QMN1)/2.0
        const float QMNM = ( QMN0 + QMN1 ) / 2.0;
        //F3983 !
        //F3984 !  ************************************************************
        //F3985 !
        //F3986 !   PUT TOTAL FORCINGS AT MIDPOINT OF YEAR JC INTO ARRAYS.
        //F3987 !
        //F3988       QQQNHL(JC)=QSNHLM+QGHM+QBIOM+QLANDM+QMNM
        QSPLIT->QQQNHL[ JC ] = QSNHLM + QGHM + QBIOM + QLANDM + QMNM;
        //        if (JC>235) cout << "QQQNHL " << JC << " " <<  QSPLIT->QQQNHL[ JC ] << " " << QSNHLM << " " << QGHM << " " << QBIOM << " " << QLANDM << " " << QMNM << endl; //debug
        //F3989       QQQNHO(JC)=QSNHOM+QGHM+QBIOM+QLANDM+QMNM
        QSPLIT->QQQNHO[ JC ] = QSNHOM + QGHM + QBIOM + QLANDM + QMNM;
        //F3990       QQQSHL(JC)=QSSHLM+QGHM+QBIOM+QLANDM+QMNM
        QSPLIT->QQQSHL[ JC ] = QSSHLM + QGHM + QBIOM + QLANDM + QMNM;
        //F3991       QQQSHO(JC)=QSSHOM+QGHM+QBIOM+QLANDM+QMNM
        QSPLIT->QQQSHO[ JC ] = QSSHOM + QGHM + QBIOM + QLANDM + QMNM;
        //F3992       FNHL=QQQNHL(JC)*FNL
        const float FNHL = QSPLIT->QQQNHL[ JC ] * AREAS->FNL;
        //F3993       FNHO=QQQNHO(JC)*FNO
        const float FNHO = QSPLIT->QQQNHO[ JC ] * AREAS->FNO;
        //F3994       FSHL=QQQSHL(JC)*FSL
        const float FSHL = QSPLIT->QQQSHL[ JC ] * AREAS->FSL;
        //F3995       FSHO=QQQSHO(JC)*FSO
        const float FSHO = QSPLIT->QQQSHO[ JC ] * AREAS->FSO;
        //F3996       QGLOBE(JC)=FNHL+FNHO+FSHL+FSHO
        QSPLIT->QGLOBE[ JC ] = FNHL + FNHO + FSHL + FSHO;
        //F3997       TEQU(JC)=TE*QGLOBE(JC)/Q2X
        TANDSL->TEQU[ JC ] = CLIM->TE * QSPLIT->QGLOBE[ JC ] / CLIM->Q2X;
        //F3998 !
        //F3999 !  CALCULATE FORCING INCREMENTS OVER YEAR IN WHICH TIME STEP LIES.
        //F4000 !   DONE FOR 1ST AND 2ND HALVES OF YEAR SEPARATELY.
        //F4001 !
        //F4002       DSNHL0M =(QSNHLM-QSNHL0)*2.0
        const float DSNHL0M = ( QSNHLM - QSNHL0 ) * 2.0;
        //F4003       DSSHL0M =(QSSHLM-QSSHL0)*2.0
        const float DSSHL0M = ( QSSHLM - QSSHL0 ) * 2.0;
        //F4004       DSNHO0M =(QSNHOM-QSNHO0)*2.0
        const float DSNHO0M = ( QSNHOM - QSNHO0 ) * 2.0;
        //F4005       DSSHO0M =(QSSHOM-QSSHO0)*2.0
        const float DSSHO0M = ( QSSHOM - QSSHO0 ) * 2.0;
        //F4006       DGH0M   =(QGHM  -QGH0  )*2.0
        const float DGH0M = ( QGHM - QGH0 ) * 2.0;
        //F4007       DBIO0M  =(QBIOM -QBIO0 )*2.0
        const float DBIO0M = ( QBIOM - QBIO0 ) * 2.0;
        //F4008       DLAND0M =(QLANDM-QLAND0)*2.0
        const float DLAND0M = ( QLANDM - QLAND0 ) * 2.0;
        //F4009       DMN0M   =(QMNM-QMN0)*2.0
        const float DMN0M = ( QMNM - QMN0 ) * 2.0;
        //F4010 !
        //F4011       DSNHLM1 =(QSNHL1-QSNHLM)*2.0
        const float DSNHLM1 = ( QSNHL1 - QSNHLM ) * 2.0;
        //F4012       DSSHLM1 =(QSSHL1-QSSHLM)*2.0
        const float DSSHLM1 = ( QSSHL1 - QSSHLM ) * 2.0;
        //F4013       DSNHOM1 =(QSNHO1-QSNHOM)*2.0
        const float DSNHOM1 = ( QSNHO1 - QSNHOM ) * 2.0;
        //F4014       DSSHOM1 =(QSSHO1-QSSHOM)*2.0
        const float DSSHOM1 = ( QSSHO1 - QSSHOM ) * 2.0;
        //F4015       DGHM1   =(QGH1  -QGHM  )*2.0
        const float DGHM1 = ( QGH1 - QGHM ) * 2.0;
        //F4016       DBIOM1  =(QBIO1 -QBIOM )*2.0
        const float DBIOM1 = ( QBIO1 - QBIOM ) * 2.0;
        //F4017       DLANDM1 =(QLAND1-QLANDM)*2.0
        const float DLANDM1 = ( QLAND1 - QLANDM ) * 2.0;
        //F4018       DMNM1   =(QMN1-QMNM)*2.0
        const float DMNM1 = ( QMN1 - QMNM ) * 2.0;
        //F4019 !
        //F4020 !  NOW CALCULATE THE FORCING VALUES AT THE MIDPOINT OF THE TIME STEP.
        //F4021 !
        //F4022       IF(FRAC.LE.0.5)THEN
        float QSNHL, QSSHL, QSNHO, QSSHO, QGHG, QBIOG, QLANDG, QMNG;
        if( FRAC <= 0.5 ) {
            //F4023         QSNHL =QSNHL0+FRAC*DSNHL0M
            QSNHL = QSNHL0 + FRAC * DSNHL0M;
            //F4024         QSSHL =QSSHL0+FRAC*DSSHL0M
            QSSHL = QSSHL0 + FRAC * DSSHL0M;
            //F4025         QSNHO =QSNHO0+FRAC*DSNHO0M
            QSNHO = QSNHO0 + FRAC * DSNHO0M;
            //F4026         QSSHO =QSSHO0+FRAC*DSSHO0M
            QSSHO = QSSHO0 + FRAC * DSSHO0M;
            //F4027         QGHG  =QGH0  +FRAC*DGH0M
            QGHG = QGH0 + FRAC * DGH0M;
            //F4028         QBIOG =QBIO0 +FRAC*DBIO0M
            QBIOG = QBIO0 + FRAC * DBIO0M;
            //F4029         QLANDG=QLAND0+FRAC*DLAND0M
            QLANDG = QLAND0 + FRAC * DLAND0M;
            //F4030         QMNG  =QMN0+FRAC*DMN0M
            QMNG = QMN0 + FRAC * DMN0M;
            //F4031       ELSE
        } else {
            //F4032         QSNHL =QSNHLM+(FRAC-0.5)*DSNHLM1
            QSNHL = QSNHLM + ( FRAC - 0.5 ) * DSNHLM1;
            //F4033         QSSHL =QSSHLM+(FRAC-0.5)*DSSHLM1
            QSSHL = QSSHLM + ( FRAC - 0.5 ) * DSSHLM1;
            //F4034         QSNHO =QSNHOM+(FRAC-0.5)*DSNHOM1
            QSNHO = QSNHOM + ( FRAC - 0.5 ) * DSNHOM1;
            //F4035         QSSHO =QSSHOM+(FRAC-0.5)*DSSHOM1
            QSSHO = QSSHOM + ( FRAC - 0.5 ) * DSSHOM1;
            //F4036         QGHG  =QGHM  +(FRAC-0.5)*DGHM1
            QGHG = QGHM + ( FRAC - 0.5 ) * DGHM1;
            //F4037         QBIOG =QBIOM +(FRAC-0.5)*DBIOM1
            QBIOG = QBIOM + ( FRAC - 0.5 ) * DBIOM1;
            //F4038         QLANDG=QLANDM+(FRAC-0.5)*DLANDM1
            QLANDG = QLANDM + ( FRAC - 0.5 ) * DLANDM1;
            //F4039         QMNG  =QMNM+(FRAC-0.5)*DMNM1
            QMNG = QMNM + ( FRAC - 0.5 ) * DMNM1;
            //F4040       ENDIF
        }
        //F4041 !
        //F4042 !  COMBINE GREENHOUSE, (SO4 AEROSOL + TROP O3 + QEXTRA) AND BIO
        //F4043 !   AEROSOL FORCINGS
        //F4044 !
        //F4045       QNHL=QGHG+QSNHL+QBIOG+QLANDG+QMNG
        QSPLIT->QNHL = QGHG + QSNHL + QBIOG + QLANDG + QMNG;
        //F4046       QNHO=QGHG+QSNHO+QBIOG+QLANDG+QMNG
        QSPLIT->QNHO = QGHG + QSNHO + QBIOG + QLANDG + QMNG;
        //F4047       QSHL=QGHG+QSSHL+QBIOG+QLANDG+QMNG
        QSPLIT->QSHL = QGHG + QSSHL + QBIOG + QLANDG + QMNG;
        //F4048       QSHO=QGHG+QSSHO+QBIOG+QLANDG+QMNG
        QSPLIT->QSHO = QGHG + QSSHO + QBIOG + QLANDG + QMNG;
        //F4049 !
        //F4050 !  **********************************************************
        //F4051 !
        //F4052       DO 10 II=1,2
        for( int II=1; II<=2; II++ ) {
            //F4053 !
            //F4054 !  *****  START OF DIFFERENTIAL SENSITIVITY TERMS  *****
            //F4055 !
            //F4056         IF(IXLAM.EQ.1)THEN
            if( DSENS->IXLAM == 1 ) {
                //F4057           WWW=FO(II)*(XKLO+FL(II)*XLAML)
                CLIM->WWW = CLIM->FO[ II ] * ( CLIM->XKLO + CLIM->FL[ II ] * DSENS->XLAML );
                //F4058           XLLDIFF(II)=ADJUST*(XLAMO+XLAML*XKLO*FL(II)/WWW)
                XLLDIFF[ II ] = DSENS->ADJUST * ( DSENS->XLAMO + DSENS->XLAML * CLIM->XKLO * CLIM->FL[ II ] / CLIM->WWW );
                //F4059           XLL=XLLDIFF(II)
                CLIM->XLL = XLLDIFF[ II ];
                //F4060         ELSE
            } else {
                //F4061           WWW=FO(II)*(XKLO+FL(II)*XLAM)
                CLIM->WWW = CLIM->FO[ II ] * ( CLIM->XKLO + CLIM->FL[ II ] * CLIM->XLAM );
                //F4062           XLLGLOBE(II)=ADJUST*XLAM*(1.+XKLO*FL(II)/WWW)
                XLLGLOBE[ II ] = DSENS->ADJUST * CLIM->XLAM * ( 1.0 + CLIM->XKLO * CLIM->FL[ II ] / CLIM->WWW );
                //F4063           XLL=XLLGLOBE(II)
                CLIM->XLL = XLLGLOBE[ II ];
                //F4064         ENDIF
            }
            //F4065 !
            //F4066 !  *****  END OF DIFFERENTIAL SENSITIVITY TERMS  *****
            //F4067 !
            //F4068         CL=-DTZ*(YYY+W(II))
            CLIM->CL = -CLIM->DTZ * ( CLIM->YYY + VARW->W[ II ] );
            //F4069         BL=1.-AL-CL
            CLIM->BL = 1.0 - CLIM->AL - CLIM->CL;
            //F4070 !
            //F4071         A(1)=1.+DTH*(XXX+W(II)*PI+XLL/FK)
            A[ 1 ] = 1.0 + CLIM->DTH * ( CLIM->XXX + VARW->W[ II ] * CLIM->PI + CLIM->XLL / CLIM->FK );
            //F4072         B(1)=-DTH*(XXX+W(II))
            B[ 1 ] = -CLIM->DTH * ( CLIM->XXX + VARW->W[ II ] );
            //F4073         A(2)=-DTZ*XXX
            A[ 2 ] = -CLIM->DTZ * CLIM->XXX;
            //F4074         C(2)=CL
            C[ 2 ] = CLIM->CL;
            //F4075         B(2)=1.-A(2)-C(2)
            B[ 2 ] = 1.0 - A[ 2 ] - C[ 2 ];
            //F4076         D(2)=TO(II,2)
            D[ 2 ] = CLIM->TO[ II ][ 2 ];
            //F4077 !
            //F4078 !  TERMS FOR VARIABLE W
            //F4079 !
            //F4080         DTZW=DTZ*DW(II)
            const float DTZW = CLIM->DTZ * VARW->DW[ II ];
            //F4081 !
            //F4082       IF(IVARW.GE.1)THEN
            if( VARW->IVARW >= 1 ) {
                //F4083         IF(IOLDTZ.EQ.1)THEN
                if( QADD->IOLDTZ == 1 )
                    //F4084           D(2)=D(2)+DTZW*(TEMEXP(II,3)-TEMEXP(II,2))
                    D[ 2 ] += DTZW * ( NSIM->TEMEXP[ II ][ 3 ] - NSIM->TEMEXP[ II ][ 2 ] );
                //F4085         ELSE
                else
                    //F4086           D(2)=D(2)+DTZW*(TEM(3)-TEM(2))
                    D[ 2 ] += DTZW * ( CLIM->TEM[ 3 ] - CLIM->TEM[ 2 ] );
                
                //F4087         ENDIF
                //F4088       ENDIF
            }
            //F4089 !
            //F4090         DO L=3,39
            for( int L=3; L<=39; L++ ) {
                //F4091           A(L)=AL
                A[ L ] = CLIM->AL;
                //F4092           B(L)=BL
                B[ L ] = CLIM->BL;
                //F4093           C(L)=CL
                C[ L ] = CLIM->CL;
                //F4094           D(L)=TO(II,L)
                D[ L ] = CLIM->TO[ II ][ L ];
                //F4095 !
                //F4096 !  TERMS FOR VARIABLE W
                //F4097 !
                //F4098 !      TEMEXP(I,L)=TP0(I)+(TO0(I)-TP0(I))*EXP(-W0*Z(L)/XK)
                //F4099         IF(IVARW.GE.1)THEN
                if( VARW->IVARW >= 1 ) {
                    //F4100           IF(IOLDTZ.EQ.1)THEN
                    if( QADD->IOLDTZ == 1 )
                        //F4101             D(L)=D(L)+DTZW*(TEMEXP(II,L+1)-TEMEXP(II,L))
                        D[ L ] += DTZW * ( NSIM->TEMEXP[ II ][ L+1 ] - NSIM->TEMEXP[ II ][ L ] );
                    //F4102           ELSE
                    else
                        //F4103             D(L)=D(L)+DTZW*(TEM(L+1)-TEM(L))
                        D[ L ] += DTZW * ( CLIM->TEM[ L+1 ] - CLIM->TEM[ L ] );
                    //F4104           ENDIF
                    //F4105         ENDIF
                }
                //F4106 !
                //F4107         END DO
            }
            //F4108 !
            //F4109         A(40)=AL
            A[ 40 ] = CLIM->AL;
            //F4110         B(40)=1.-CL
            B[ 40 ] = 1.0 - CLIM->CL;
            //F4111         D(40)=TO(II,40)+TO(II,1)*PI*DTZ*W(II)
            D[ 40 ] = CLIM->TO[ II ][ 40 ] + CLIM->TO[ II ][ 1 ] * CLIM->PI * CLIM->DTZ * VARW->W[ II ];
            //F4112 !
            //F4113 !  TERMS FOR VARIABLE W
            //F4114 !
            //F4115         IF(IVARW.GE.1)THEN
            if( VARW->IVARW >= 1 ) {
                //F4116           IF(IOLDTZ.EQ.1)THEN
                if( QADD->IOLDTZ == 1 ) {
                    //F4117             D(40)=D(40)+DTZW*(TP0(II)-TEMEXP(II,40))
                    D[ 40 ] += DTZW * ( VARW->TP0[ II ] - NSIM->TEMEXP[ II ][ 40 ] );
                    //F4118           ELSE
                } else {
                    //F4119             D(40)=D(40)+DTZW*(TP0(II)-TEM(40))
                    D[ 40 ] += DTZW * ( VARW->TP0[ II ] - CLIM->TEM[ 40 ] );
                }
                //F4120           ENDIF
                //F4121         ENDIF
            }
            //F4122 !
            //F4123 !  FORL is the land forcing term
            //F4124 !
            //F4125         if(ii.eq.1) then
            float FORL, HEAT;
            if( II == 1 ) {
                //F4126           forl = qnhl*xklo*fl(ii)/www
                FORL = QSPLIT->QNHL * CLIM->XKLO * CLIM->FL[ II ] / CLIM->WWW;
                //F4127           heat = (qnho+hem(ii)+forl)*dth/fk
                HEAT = ( QSPLIT->QNHO + CLIM->HEM[ II ] + FORL ) * CLIM->DTH / CLIM->FK;
                //F4128         else
            } else {
                //F4129           forl = qshl*xklo*fl(ii)/www
                FORL = QSPLIT->QSHL * CLIM->XKLO * CLIM->FL[ II ] / CLIM->WWW;
                //F4130           heat = (qsho+hem(ii)+forl)*dth/fk
                HEAT = ( QSPLIT->QSHO + CLIM->HEM[ II ] + FORL ) * CLIM->DTH / CLIM->FK;
                //F4131         endif
            }
            //F4132 !
            //F4133         D(1)=TO(II,1)+HEAT
            D[ 1 ] = CLIM->TO[ II ][ 1 ] + HEAT;
            //F4134 !
            //F4135 !  TERMS FOR VARIABLE W
            //F4136 !
            //F4137         IF(IVARW.GE.1)THEN
            if( VARW->IVARW >= 1 ) {
                //F4138           IF(IOLDTZ.EQ.1)THEN
                if( QADD->IOLDTZ == 1 )
                    //F4139             D(1)=D(1)+DTH*DW(II)*(TEMEXP(II,2)-TP0(II))
                    D[ 1 ] += CLIM->DTH * VARW->DW[ II ] * ( NSIM->TEMEXP[ II ][ 2 ] - VARW->TP0[ II ] );
                //F4140           ELSE
                else
                    //F4141             D(1)=D(1)+DTH*DW(II)*(TEM(2)-TP0(II))
                    D[ 1 ] += CLIM->DTH * VARW->DW[ II ] * ( CLIM->TEM[ 2 ] - VARW->TP0[ II ] );
                //F4142           ENDIF
                //F4143         ENDIF
            }
            //F4144 !
            //F4145 !  THIS IS THE OLD BAND SUBROUTINE
            //F4146 !
            //F4147         AA(1)=-B(1)/A(1)
            AA[ 1 ] = -B[ 1 ] / A[ 1 ];
            //F4148         BB(1)=D(1)/A(1)
            BB[ 1 ] = D[ 1 ] / A[ 1 ];
            //F4149         DO L=2,39
            for( int L=2; L<=39; L++ ) {
                //F4150           VV=A(L)*AA(L-1)+B(L)
                float VV = A[ L ] * AA[ L-1 ] + B[ L ];
                //F4151           AA(L)=-C(L)/VV
                AA[ L ] = -C[ L ] / VV;
                //F4152           BB(L)=(D(L)-A(L)*BB(L-1))/VV
                BB[ L ] = ( D[ L ] - A[ L ] * BB[ L-1 ] ) / VV;
                //F4153         END DO
            }
            //F4154         TO(II,40)=(D(40)-A(40)*BB(39))/(A(40)*AA(39)+B(40))
            CLIM->TO[ II ][ 40 ] = ( D[ 40 ] - A[ 40 ] * BB[ 39 ] ) / ( A[ 40 ] * AA[ 39 ] + B[ 40 ] );
            //F4155         DO I=1,39
            for( int I=1; I<=39; I++ ) {
                //F4156           L=40-I
                int L = 40 - I;
                //F4157           TO(II,L)=AA(L)*TO(II,L+1)+BB(L)
                CLIM->TO[ II ][ L ] = AA[ L ] * CLIM->TO[ II ][ L+1 ] + BB[ L ];
                //F4158         END DO
            }
            //F4159 !
            //F4160   10  CONTINUE
        }
        //F4161 !
        //F4162 !  Y(1,2,3,4) ARE NH OCEAN, SH OCEAN, NH LAND & SH LAND TEMPS.
        //F4163 !
        //F4164       Y(1)=TO(1,1)*ADJUST
        CLIM->Y[ 1 ] = CLIM->TO[ 1 ][ 1 ] * DSENS->ADJUST;
        //F4165       Y(2)=TO(2,1)*ADJUST
        CLIM->Y[ 2 ] = CLIM->TO[ 2 ][ 1 ] * DSENS->ADJUST;
        //F4166 !
        //F4167 !  DIFFERENTIAL SENSITIVITY TERMS
        //F4168 !
        //F4169       IF(IXLAM.EQ.1)THEN
        if( DSENS->IXLAM == 1 ) {
            //F4170         Y(3)=(FL(1)*qnhl+XKLO*Y(1))/(FL(1)*XLAML+XKLO)
            CLIM->Y[ 3 ] = ( CLIM->FL[ 1 ] * QSPLIT->QNHL + CLIM->XKLO * CLIM->Y[ 1 ] ) / ( CLIM->FL[ 1 ] * DSENS->XLAML + CLIM->XKLO );
            //F4171         Y(4)=(FL(2)*qshl+XKLO*Y(2))/(FL(2)*XLAML+XKLO)
            CLIM->Y[ 4 ] = ( CLIM->FL[ 2 ] * QSPLIT->QSHL + CLIM->XKLO * CLIM->Y[ 2 ] ) / ( CLIM->FL[ 2 ] * DSENS->XLAML + CLIM->XKLO );
            //F4172       ELSE
        } else {
            //F4173         Y(3)=(FL(1)*qnhl+XKLO*Y(1))/(FL(1)*XLAM+XKLO)
            CLIM->Y[ 3 ] = ( CLIM->FL[ 1 ] * QSPLIT->QNHL + CLIM->XKLO * CLIM->Y[ 1 ] ) / ( CLIM->FL[ 1 ] * CLIM->XLAM + CLIM->XKLO );
            //F4174         Y(4)=(FL(2)*qshl+XKLO*Y(2))/(FL(2)*XLAM+XKLO)
            CLIM->Y[ 4 ] = ( CLIM->FL[ 2 ] * QSPLIT->QSHL + CLIM->XKLO * CLIM->Y[ 2 ] ) / ( CLIM->FL[ 2 ] * CLIM->XLAM + CLIM->XKLO );
            //F4175       ENDIF
        }
        //F4176 !
        //F4177 !  *****  END OF DIFFERENTIAL SENSITIVITY TERMS  *****
        //F4178 !
        //F4179       HEM(1)=(XKNS/FO(1))*(Y(2)-Y(1))
        CLIM->HEM[ 1 ] = ( CLIM->XKNS / CLIM->FO[ 1 ] ) * ( CLIM->Y[ 2 ] - CLIM->Y[ 1 ] );
        //F4180       HEM(2)=(XKNS/FO(2))*(Y(1)-Y(2))
        CLIM->HEM[ 2 ] = ( CLIM->XKNS / CLIM->FO[ 2 ] ) * ( CLIM->Y[ 1 ] - CLIM->Y[ 2 ] );
        //F4181 !
        //F4182 !  VARIABLE W TERMS
        //F4183 !
        //F4184       IF(IVARW.EQ.0)THEN
        if( VARW->IVARW == 0 ) {
            //F4185         W(1)=W0
            VARW->W[ 1 ] = CLIM->W0;
            //F4186         DW(1)=0.0
            VARW->DW[ 1 ] = 0.0;
            //F4187         W(2)=W0
            VARW->W[ 2 ] = CLIM->W0;
            //F4188         DW(2)=0.0
            VARW->DW[ 2 ] = 0.0;
            //F4189       ENDIF
        }
        //F4190 !
        //F4191       OCEANT=(FO(1)*Y(1)+FO(2)*Y(2))/(FO(1)+FO(2))
        const float OCEANT = ( CLIM->FO[ 1 ] * CLIM->Y[ 1 ] + CLIM->FO[ 2 ] * CLIM->Y[ 2 ] ) / ( CLIM->FO[ 1 ] + CLIM->FO[ 2 ] );
        //F4192       GLOBET=(FO(1)*Y(1)+FL(1)*Y(3)+FO(2)*Y(2)+FL(2)*Y(4))/2.0
        const float GLOBET = ( CLIM->FO[ 1 ] * CLIM->Y[ 1 ] + CLIM->FL[ 1 ] * CLIM->Y[ 3 ] + CLIM->FO[ 2 ] * CLIM->Y[ 2 ]  + CLIM->FL[ 2 ] * CLIM->Y[ 4 ] ) / 2.0;
        //F4193 !
        //F4194 !  ALTERNATIVE WAYS TO DEFINE W(t) (SEE TOP OF CODE FOR DETAILS).
        //F4195 !
        //F4196       AW=1.0
        float AW = 1.0;
        //F4197       IF((IVARW.EQ.2).AND.(KEYDW.GE.4))AW=1.0-WTHRESH/W0
        if( VARW->IVARW == 2 && VARW->KEYDW >= 4 ) AW = 1.0 - NSIM->WTHRESH / CLIM->W0;
        //F4198 !
        //F4199       IF((KEYDW.EQ.1).OR.(KEYDW.EQ.4))THEN
        float TNKEY, TSKEY;
        if( VARW->KEYDW == 1 || VARW->KEYDW == 4 )
            //F4200         TNKEY=GLOBET
            //F4201         TSKEY=GLOBET
            TNKEY = TSKEY = GLOBET;
        //F4202       ENDIF
        //F4203       IF((KEYDW.EQ.2).OR.(KEYDW.EQ.5))THEN
        if( VARW->KEYDW == 2 || VARW->KEYDW == 5 )
            //F4204         TNKEY=OCEANT
            //F4205         TSKEY=OCEANT
            TNKEY = TSKEY = OCEANT;
        //F4206       ENDIF
        //F4207       IF(KEYDW.EQ.3)THEN
        if( VARW->KEYDW == 3 ) {
            //F4208         TNKEY=Y(1)
            TNKEY = CLIM->Y[ 1 ];
            //F4209         TSKEY=Y(2)
            TSKEY = CLIM->Y[ 2 ];
            //F4210       ENDIF
        }
        //F4211 !
        //F4212 !  DROP FULL W TO 0.1 BUT RECOVER IF TEMPERATURE RECOVERS
        //F4213 !
        //F4214       IF(IVARW.EQ.1)THEN
        if( VARW->IVARW == 1 ) {
            //F4215         DW(1)=-AW*W0*TNKEY/TW0NH
            VARW->DW[ 1 ] = -AW * CLIM->W0 * TNKEY / VARW->TW0NH;
            //F4216         W(1)=W0+DW(1)
            VARW->W[ 1 ] = CLIM->W0 + VARW->DW[ 1 ];
            //F4217         IF(W(1).LT.0.1)W(1)=0.1
            if( VARW->W[ 1 ] < 0.1 ) VARW->W[ 1 ] = 0.1;
            //F4218 !
            //F4219         DW(2)=-AW*W0*TSKEY/TW0SH
            VARW->DW[ 2 ] = -AW * CLIM->W0 * TSKEY / VARW->TW0SH;
            //F4220         W(2)=W0+DW(2)
            VARW->W[ 2 ] = CLIM->W0 + VARW->DW[ 2 ];
            //F4221         IF(W(2).LT.0.1)W(1)=0.1
            if( VARW->W[ 2 ] < 0.1 ) VARW->W[ 1 ] = 0.1;
            //F4222       ENDIF
        }
        //F4223 !
        //F4224 !  DROP TO EITHER FULL OR ACTIVE W TO WTHRESH AND STAY THERE
        //F4225 !
        //F4226       IF(IVARW.EQ.2)THEN
        if( VARW->IVARW == 2 ) {
            //F4227         DW(1)=-AW*W0*TNKEY/TW0NH
            VARW->DW[ 1 ] = -AW * CLIM->W0 * TNKEY / VARW->TW0NH;
            //F4228         W(1)=W0+DW(1)
            VARW->W[ 1 ] = CLIM->W0 + VARW->DW[ 1 ];
            //F4229         IF((W(1).LT.WTHRESH).OR.(IWNHOFF.EQ.1))THEN
            if( VARW->W[ 1 ] < NSIM->WTHRESH || NSIM->IWNHOFF == 1 ) {
                //F4230           IWNHOFF=1
                NSIM->IWNHOFF = 1;
                //F4231           W(1)=WTHRESH
                VARW->W[ 1 ] = NSIM->WTHRESH;
                //F4232           DW(1)=WTHRESH-W0
                VARW->DW[ 1 ] = NSIM->WTHRESH - CLIM->W0;
                //F4233         ENDIF
            }
            //F4234 !
            //F4235         DW(2)=-AW*W0*TSKEY/TW0SH
            VARW->DW[ 2 ] = -AW * CLIM->W0 * TSKEY / VARW->TW0SH;
            //F4236         W(2)=W0+DW(2)
            VARW->W[ 2 ] = CLIM->W0 + VARW->DW[ 2 ];
            //F4237         IF((W(2).LT.WTHRESH).OR.(IWSHOFF.EQ.1))THEN
            if( VARW->W[ 2 ] < NSIM->WTHRESH || NSIM->IWNHOFF == 1 ) {
                //F4238           IWSHOFF=1
                NSIM->IWNHOFF = 1;
                //F4239           W(2)=WTHRESH
                VARW->W[ 2 ] = NSIM->WTHRESH;
                //F4240           DW(2)=WTHRESH-W0
                VARW->DW[ 2 ] = NSIM->WTHRESH - CLIM->W0;
                //F4241         ENDIF
            }
            //F4242       ENDIF
        }
        //F4243 !
        //F4244 !  SPACE FOR CODE TO USE INPUT WNH AND WSH TIME SERIES
        //F4245 !
        //F4246 !      IF(IVARW.EQ.3)THEN
        //F4247 !        if(jc.ge.246) w(1)=wthresh
        //F4248 !        if(jc.ge.246) w(2)=wthresh
        //F4249 !        if(jc.ge.246) dw(1)=0.0
        //F4250 !        if(jc.ge.246) dw(2)=0.0
        //F4251         WNH(JC)=W(1)
        VARW->WNH[ JC ] = VARW->W[ 1 ];
        //F4252         WSH(JC)=W(2)
        VARW->WSH[ JC ] = VARW->W[ 2 ];
        //F4253 !
        //F4254 !  HAVING CALCULATED VALUES AT NEW TIME 'T', DECIDE WHETHER OR NOT
        //F4255 !   TEMPS ARE ANNUAL VALUES (I.E., CORRESPOND TO T=MIDPOINT OF
        //F4256 !   YEAR).  IF SO, GO TO TSLCALC, CALCULATE GLOBAL MEAN TEMP (TGAV)
        //F4257 !   AND INSERT INTO ARRAY, AND CALCULATE SEA LEVEL RISE COMPONENTS.
        //F4258 !  NOTE THAT, BECAUSE FORCING VALUES ARE GIVEN AT ENDS OF YEARS
        //F4259 !   AND TIMES ARE INTEGERS AT MIDPOINTS OF YEARS, A ONE YEAR TIME
        //F4260 !   STEP WILL STILL GIVE MID-YEAR VALUES FOR CALCULATED TEMPERATURES.
        //F4261 !
        //F4262       KP=KC
        int KP = CLIM->KC;
        //F4263       KC=INT(T+1.01)
        CLIM->KC = int( CLIM->T + 1.01 );
        //F4264       IF(KC.GT.KP)CALL TSLCALC(KC)
        if( CLIM->KC > KP ) tslcalc( CLIM->KC, Limits, CLIM, CONCS, CARB,
                                    TANDSL, VARW, QSPLIT, ICE, NSIM, outfile8 );
        //F4265 !
        //F4266       IF(T.GE.TEND)RETURN
        //F4267       GO TO  11
    } while ( CLIM->T < CLIM->TEND );
    //F4268       END
    f_exit( __func__ );
} // runmod
//F4269 !
//F4270 !  *******************************************************************
//F4271 !
//F4272       SUBROUTINE DELTAQ
void deltaq( Limits_block* Limits, OZ_block* OZ, CLIM_block* CLIM, CONCS_block* CONCS,
            NEWCONCS_block* NEWCONCS, CARB_block* CARB, TANDSL_block* TANDSL, CAR_block* CAR,
            METH1_block* METH1, FORCE_block* FORCE, METH2_block* METH2, METH3_block* METH3,
            METH4_block* METH4, TauNitr_block* TauNitr, Sulph_block* Sulph, NSIM_block* NSIM, 
            CO2READ_block* CO2READ, JSTART_block* JSTART, CORREN_block* CORREN, HALOF_block* HALOF, COBS_block* COBS )
{
    f_enter( __func__ );
    //F4273       IMPLICIT REAL*4 (a-h,o-z), Integer (I-N)
    //F4274 !
    //F4275       parameter (iTp=740)
    //F4276 !
    //F4277       common /Limits/KEND
    //F4278 !
    //F4279       COMMON/OZ/OZ00CH4,OZCH4,OZNOX,OZCO,OZVOC
    //F4280 !
    //F4281       COMMON/CLIM/IC,IP,KC,DT,DZ,FK,HM,Q2X,QXX,PI,T,TE,TEND,W0,XK,XKLO, &
    //F4282       XKNS,XLAM,FL(2),FO(2),FLSUM,FOSUM,HEM(2),P(40),TEM(40),TO(2,40), &
    //F4283       AL,BL,CL,DTH,DTZ,DZ1,XLL,WWW,XXX,YYY,RHO,SPECHT,HTCONS,Y(4)
    //F4284 !
    //F4285       COMMON/CONCS/CH4(0:iTp),CN2O(0:iTp),ECH4(226:iTp+1), &
    //F4286       EN2O(226:iTp+1),ECO(226:iTp+1),COE(iTp+1),EVOC(226:iTp+1), &
    //F4287       ENOX(226:iTp+1),ESO2(0:iTp+1),ESO2SUM(226:iTp+1), &
    //F4288       ESO21(226:iTp+1),ESO22(226:iTp+1),ESO23(226:iTp+1), &
    //F4289       EBC(226:iTp+1), EOC(226:iTp+1) ! sjs- add BC-OC
    //F4290 !
    //F4291       COMMON/NEWCONCS/CF4(iTp),C2F6(iTp),C125(iTp),C134A(iTp), &
    //F4292       C143A(iTp),C227(iTp),C245(iTp),CSF6(iTp), &
    //F4293       ECF4(226:iTp+1),EC2F6(226:iTp+1),E125(226:iTp+1),E134A(226:iTp+1), &
    //F4294       E143A(226:iTp+1),E227(226:iTp+1),E245(226:iTp+1),ESF6(226:iTp+1)
    //F4295 !
    //F4296       COMMON/CARB/CCO2(4,224:iTp),EDGROSS(4,226:iTp),EF(226:iTp+1), &
    //F4297       REGROW(4,226:iTp),PL(4,226:iTp),HL(4,226:iTp),SOIL(4,226:iTp), &
    //F4298       TTT(226:iTp),ESUM(226:iTp),ETOT(4,226:iTp),EDNET90(4), &
    //F4299       FOC(4,226:iTp),CO2(0:iTp),CO2SAVE(0:iTp)
    //F4300 !
    //F4301       COMMON/TANDSL/TEQU(iTp),TGAV(iTp),TNHO(iTp), &
    //F4302       TSHO(iTp),TNHL(iTp),TSHL(iTp),TDEEP(iTp),TNHAV(iTp),TSHAV(iTp), &
    //F4303       TLAND(iTp),TOCEAN(iTp),TOCN(40),TOCNPREV(40), &
    //F4304       SIP,SGP,SAP,SLI(iTp),SLG(iTp),SLA(iTp),EX(0:iTp),SLT(iTp), &
    //F4305       QTOT(0:iTp),QGH(0:iTp),QOZ(0:iTp),QBIO(0:iTp),SLO(iTp), &
    //F4306       QSO2(0:iTp+1),QDIR(0:iTp+1),QLAND(0:iTp),QMN(0:iTp+1)
    //F4307 !
    //F4308       COMMON/CAR/EL1,EL2,EL3,TINV0(5),TINV(4,5),A(3,5),AA(4,5), &
    //F4309       BCO2(4),BTGPP,BTRESP,BTHUM,GAMP,GPP0,RESP0,QA0,U0,C0,B340(4), &
    //F4310       PHI,RG,TAUP,TAUH,TAUS,THP,THS,THH0,THS0,THPL,G1,G2,G3,FACTOR, &
    //F4311       EL21,EL32,XX1,XX2,XX3,XX4,XX5,XX6,DEE1,DEE2,DEE3,DEE4,DEE5,DEE6, &
    //F4312       FL1,FL2,FL3,XL,GAMH,GAMS,QS0,BTSOIL,FERTTYPE,TOTEM,CONVTERP, &
    //F4313       R(4),CPART(4,5),DELMASS(4,226:iTp),ABFRAC(4,226:iTp)
    //F4314 !
    //F4315       COMMON /METH1/emeth(226:iTp),imeth,ch4l(225:iTp),ch4b(225:iTp), &
    //F4316       ch4h(225:iTp),ef4(226:iTp),StratH2O,TCH4(iTp),iO3feed, &
    //F4317       ednet(226:iTp+1),DUSER,FUSER,CORRUSER,CORRMHI,CORRMMID,CORRMLO
    //F4318 !
    //F4319       COMMON /FORCE/qco2(0:iTp),qm(0:iTp),qn(0:iTp),QCFC(0:iTp), &
    //F4320       QMONT(0:iTp),QOTHER(0:iTp),QSTRATOZ(0:iTp),QCH4O3(0:iTp), &
    //F4321       CFC12(0:iTp), QCH4H2O(0:iTp),QBC(0:iTp),QOC(0:iTp)
    //F4322 !
    //F4323       COMMON /METH2/LEVCH4,ch4bar90,QQQN2O
    //F4324 !
    //F4325       COMMON /METH3/TCH4CON,TAUINIT,SCH4,DELSS,DELTAU, &
    //F4326       ANOX,ACO,AVOC,DELANOX,DELACO,DELAVOC,ICH4FEED
    //F4327 !
    //F4328       COMMON /METH4/GAM,TAUOTHER,BBCH4,CM00
    //F4329       common /TauNitr/TN2000,BBN2O,SN2O,CN00,NOFFSET
    //F4330       common /Sulph/S90DIR,S90IND,S90BIO,ENAT,ES1990,ECO90,FOC90,IFOC
    //F4331 !
    //F4332       COMMON /NSIM/NSIM,NCLIM,ISCENGEN,TEMEXP(2,40),IWNHOFF,IWSHOFF, &
    //F4333       WTHRESH
    //F4334 !
    //F4335       COMMON /CO2READ/ICO2READ,XC(226:iTp),CO2SCALE,qtot86,LEVCO2
    //F4336 !
    //F4337       COMMON /JSTART/JSTART,FOSSHIST(0:236),QKYMAG(0:iTp),IGHG, &
    //F4338       QCH4OZ,QFOC(0:iTp),ICO2CORR,TROZSENS
    //F4339 !
    //F4340       COMMON /CORREN/CORREN1,CORREN2,CORREN3,CORREN4,CORREN
    //F4341 !
    //F4342       SAVE T00LO,T00MID,T00HI,T00USER
    //F4343 ! sjs -- change to make MAGICC  work. need to save these vars
    static float T00LO, T00MID, T00HI, T00USER;
    //F4344 
    //F4345 ! sjs -- add storage for halocarbon variables
    //F4346       COMMON /HALOF/QCF4_ar(0:iTp),QC2F6_ar(0:iTp),qSF6_ar(0:iTp), &
    //F4347        Q125_ar(0:iTp),Q134A_ar(0:iTp), &
    //F4348        Q143A_ar(0:iTp),Q227_ar(0:iTp),Q245_ar(0:iTp)
    //F4349 
    //F4350 ! sjs-- g95 seems to have optomized away these local variables, so put them in common block
    //F4351      COMMON /TEMPSTOR/DQOZPP, DQOZ
    static float /* DQOZPP,*/ DQOZ = 0.0; // DQOZPP unused
    static float QOZ1 = 0.0;
    const float fffrac = 0.18;
    float TAUCH4 = 0.0;
    
    //F4352 
    //F4353 !
    //F4354 !  THIS SUBROUTINE IS ONLY ENTERED WHEN THE IC YEAR COUNT
    //F4355 !   INCREMENTS.  CONCENTRATIONS AND RADIATIVE FORCINGS ARE
    //F4356 !   THEN CALCULATED FOR THE END OF THE IC YEAR.  SINCE THE
    //F4357 !   IC INCREMENT MAY BE GREATER THAN 1, A LOOP IS NEEDED TO
    //F4358 !   FILL IN THE MISSING IC VALUES.
    //F4359 !
    //F4360       QLAND90=-0.2
    const float QLAND90 = -0.2;
    static float TX=0.0, DELT90=0.0, DELT00=0.0;
    //F4361 !
    //F4362       DO 10 J=IP+1,IC
    for( int J=CLIM->IP+1; J<=CLIM->IC; J++ ) {
        //F4363 !
        //F4364 !********************************************************
        //F4365 !
        //F4366 !  LAND ALBEDO CHANGE FORCING
        //F4367 !
        //F4368       IF(J.LE.226)THEN
        if( J <= 226 )
            //F4369         QLAND(J)=QLAND90*FLOAT(J)/226.0
            TANDSL->QLAND[ J ] = QLAND90 * float(J)/226.0;
        //F4370       ELSE
        else
            //F4371         QLAND(J)=QLAND90
            TANDSL->QLAND[ J ] = QLAND90;
        //F4372       ENDIF
        //F4373 !
        //F4374 !  *******************************************************
        //F4375 !
        //F4376 !  START OF LONG SET OF IF STATEMENTS (MODIFIED AUGUST, 2000).
        //F4377 !
        //F4378 !  CARBON CYCLE MODEL STILL RUNS FROM 1990, BUT A CORRECTION
        //F4379 !   IS APPLIED FOR 2000 ONWARDS TO ENSURE CONSISTENCY WITH
        //F4380 !   OBSERVED VALUES IN 2000.
        //F4381 !
        //F4382 !  HISTORY
        //F4383 !  FIRST ACCESS HISTORY DATA
        //F4384 !
        //F4385       IF(J.LE.JSTART)THEN
        if( J <= JSTART->JSTART ) {
            //F4386 !
            //F4387 !  NOTE : FOR PROTOCOL GASES, ONLY 1990 CONCENTRATION VALUE IS GIVEN
            //F4388 !   (FOR ALL J.LE.JSTART) SINCE THIS IS ALL THAT IS NEEDED TO
            //F4389 !   INITIALIZE FUTURE CONCS
            //F4390 !
            //F4391         if(j.ge.226)then
            float ejkeep, ej1keep;
            if( J >= 226 ) {
                //F4392           ejkeep=eso2(j)
                ejkeep = CONCS->ESO2.getval( J );
                //F4393           ej1keep=eso2(j+1)
                ej1keep = CONCS->ESO2.getval( J+1 );
                //F4394         endif
            }
            //F4395 !
            //F4396         CALL HISTORY(J,CO2(J),CH4(J),CN2O(J),eso2(J),eso2(j+1), &
            //F4397         CF4(J),C2F6(J),C125(J),C134A(J),C143A(J),C227(J),C245(J), &
            //F4398         CSF6(J))
            history( J, &CARB->CO2[ J ], &CONCS->CH4[ J ], &CONCS->CN2O[ J ], CONCS->ESO2.getptr( J ), CONCS->ESO2.getptr( J+1 ), 
                    &NEWCONCS->CF4[ J ], &NEWCONCS->C2F6[ J ], &NEWCONCS->C125[ J ],
                    &NEWCONCS->C134A[ J ], &NEWCONCS->C143A[ J ], &NEWCONCS->C227[ J ], &NEWCONCS->C245[ J ], 
                    &NEWCONCS->CSF6[ J ], COBS, Sulph->ES1990 );
            //F4399 !
            //F4400         if(j.ge.226)then
            if( J >= 226 ) {
                //F4401           eso2(j)=ejkeep
                CONCS->ESO2.setval( ejkeep, J );
                //F4402           eso2(j+1)=ej1keep
                CONCS->ESO2.setval( ej1keep, J+1 );
                //F4403         endif
            }
            //F4404 !
            //F4405         IF(J.EQ.JSTART)THEN
            if( J == JSTART->JSTART ) {
                //F4406           CM00=CH4(JSTART-1)
                METH4->CM00 = CONCS->CH4[ JSTART->JSTART-1 ];
                //F4407           CN00=CN2O(JSTART-1)
                TauNitr->CN00 = CONCS->CN2O[ JSTART->JSTART-1 ];
                //F4408         ENDIF
            }
            //F4409 !
            //F4410         if(j.eq.jstart) then
            if( J == JSTART->JSTART ) {
                //F4411           ch4l(jstart-1) = ch4(jstart-1)
                METH1->ch4l.setval( CONCS->CH4[ JSTART->JSTART-1 ], JSTART->JSTART-1 );
                //F4412           ch4b(jstart-1) = ch4(jstart-1)
                METH1->ch4b.setval( CONCS->CH4[ JSTART->JSTART-1 ], JSTART->JSTART-1 );
                //F4413           ch4h(jstart-1) = ch4(jstart-1)
                METH1->ch4h.setval( CONCS->CH4[ JSTART->JSTART-1 ], JSTART->JSTART-1 );
                //F4414           ch4l(jstart) = ch4(jstart)
                METH1->ch4l.setval( CONCS->CH4[ JSTART->JSTART ], JSTART->JSTART );
                //F4415           ch4b(jstart) = ch4(jstart)
                METH1->ch4b.setval( CONCS->CH4[ JSTART->JSTART ], JSTART->JSTART );
                //F4416           ch4h(jstart) = ch4(jstart)
                METH1->ch4h.setval( CONCS->CH4[ JSTART->JSTART ], JSTART->JSTART );
                //F4417         endif
            }
            //F4418 !
            //F4419 !  Calculate additional CO2 emissions due to CH4 oxidation for
            //F4420 !   jstart year.
            //F4421 !  First specify fossil fuel fraction of total CH4 emissions
            //F4422 !   (updated to 0.18 in August 2000 based on TAR).
            //F4423 !
            //F4424         if(j.ge.226)then
            if( J >= 226 ) {
                //F4425           fffrac = 0.18
                // fffrac = 0.18; // this is done @ definition above, actually
                //F4426           emeth(j) = fffrac*0.0020625*(ch4(j)-700.)/TAUINIT
                METH1->emeth.setval( fffrac * 0.0020625 * ( CONCS->CH4[ J ]-700.0 ) / METH3->TAUINIT, J ); 
                //F4427           if(imeth.eq.1) then
                if( METH1->IMETH == 1 )
                    //F4428             ef4(j) = ef(j)+emeth(j)
                    METH1->ef4.setval( CARB->EF.getval( J ) + METH1->emeth.getval( J ), J );
                //F4429           else
                else
                    //F4430             ef4(j) = ef(j)
                    METH1->ef4.setval( CARB->EF.getval( J ), J );
                //F4431           endif
                
                //F4432         endif
            }
            //F4433 !
            //F4434       ENDIF
        }
        //F4435 !
        //F4436 !  *******************************************************
        //F4437 !
        //F4438 !  FOR CONCS BEYOND THE END OF 1990, CALL THE VARIOUS EMISSIONS
        //F4439 !   TO CONCENTRATIONS MODELS.  NOTE THAT THE INPUT EMISSIONS
        //F4440 !   ARRAYS ARE NUMBERED AS FOR OTHER VARIABLES ; I.E. E(226)
        //F4441 !   IS THE 1990 VALUE. EMISSIONS VALUES UP TO AND INCLUDING E(225)
        //F4442 !   ARE NOT USED.
        //F4443 !
        //F4444 !  FOR CALCULATING FEEDBACKS IN CARBON AND METHANE, NEED TO USE
        //F4445 !   EXTRAPOLATED VALUES FOR TEMPERATURE AND CO2 CONCENTRATION.
        //F4446 !   RELEVANT VALUES ARE FROM START YEAR FOR MODEL PROJECTIONS.
        //F4447 !
        //F4448       IF(J.GE.226)THEN
        if( J >= 226 ) {
            //F4449         TX=2.0*TGAV(J-1)-TGAV(J-2)
            TX = 2.0 * TANDSL->TGAV[ J-1 ] - TANDSL->TGAV[ J-2 ];
            //F4450         IF(J.EQ.226)DELT90=TX
            if( J == 226 ) DELT90 = TX;
            //F4451         IF(J.EQ.236)DELT00=TX
            if( J == 236 ) DELT00 = TX;
            //F4452       ENDIF
        }
        //F4453 !
        //F4454       IF(J.GT.JSTART)THEN
        if( J > JSTART->JSTART ) {
            //F4455 !
            //F4456 !  SET INITIAL (YEAR 2000) LIFETIMES
            //F4457 !
            //F4458         IF(j.eq.jstart+1) then
            if( J == JSTART->JSTART+1 ) {
                //F4459           t00lo=TAUINIT-DELTAU
                T00LO = METH3->TAUINIT - METH3->DELTAU;
                //F4460           t00mid=TAUINIT
                T00MID = METH3->TAUINIT;
                //F4461           t00hi=TAUINIT+DELTAU
                T00HI = METH3->TAUINIT + METH3->DELTAU;
                //F4462           if(LEVCH4.eq.1)t00user=t00lo
                if( METH2->LEVCH4 == 1 ) T00USER = T00LO;
                //F4463           if(LEVCH4.eq.2)t00user=t00mid
                if( METH2->LEVCH4 == 2 ) T00USER = T00MID;
                //F4464           if(LEVCH4.eq.3)t00user=t00hi
                if( METH2->LEVCH4 == 3 ) T00USER = T00HI;
                //F4465           if(LEVCH4.eq.4)t00user=TCH4CON
                if( METH2->LEVCH4 == 4 ) T00USER = METH3->TCH4CON;
                //F4466         ENDIF
            }
            //F4467 !
            //F4468 !  *******************************************************
            //F4469 !
            //F4470 !  CH4
            //F4471 !  PRATHER'S TAR METHOD INCORPORATED, AUGUST 2000
            //F4472 !  FOR METHANE CONC PROJECTIONS, NEED TO USE EMISSIONS CONSISTENT WITH
            //F4473 !   THE LIFETIME. THIS IS DONE USING CORRECTION FACTORS CALCULATED
            //F4474 !   IN THE MAIN PROGRAM. NOTE THAT THE INPUT EMISSIONS HAVE ALREADY
            //F4475 !   BEEN OFFSET BY THE AMOUNT APPROPRIATE TO THE USER-SPECIFIED
            //F4476 !   LIFETIME (I.E., BY CORRUSER).
            //F4477 !  THIS WAS CORRECTED ON 97/12/13.
            //F4478 !
            //F4479         DENOX=ENOX(J)-ENOX(236)
            const float DENOX = CONCS->ENOX.getval( J ) - CONCS->ENOX.getval( 236 );
            //F4480         DECO=ECO(J) -ECO(236)
            const float DECO = CONCS->ECO.getval( J ) - CONCS->ECO.getval( 236 );
            //F4481         DEVOC=EVOC(J)-EVOC(236)
            const float DEVOC = CONCS->EVOC.getval( J ) - CONCS->EVOC.getval( 236 );
            //F4482 !
            //F4483 !  ESTIMATED TEMPERATURE CHANGE FROM 2000
            //F4484 !
            //F4485         DELTAT=TX-DELT00
            const float DELTAT = TX - DELT00;
            //F4486 !
            //F4487 !  LOW LIFETIME
            //F4488 !
            //F4489         EECH4 = ECH4(J)-CORRUSER+CORRMLO
            float EECH4 = CONCS->ECH4.getval( J ) - METH1->CORRUSER + METH1->CORRMLO;
            //F4490 !
            //F4491         SSLO=SCH4-DELSS
            const float SSLO = METH3->SCH4 - METH3->DELSS;
            //F4492         ANOXLO=ANOX+DELANOX
            const float ANOXLO = METH3->ANOX + METH3->DELANOX;
            //F4493         ACOLO=ACO-DELACO
            const float ACOLO = METH3->ACO + METH3->DELACO;
            //F4494         AVOCLO=AVOC-DELAVOC
            const float AVOCLO = METH3->AVOC + METH3->DELAVOC;
            //F4495 !
            //F4496         CALL METHANE(ICH4FEED,CH4L(J-1),EECH4,DENOX,DECO,DEVOC,CH4L(J), &
            //F4497         T00LO,TAULO,SSLO,ANOXLO,ACOLO,AVOCLO,DELTAT)
            float TAULO = 0.0;
            methane( METH3->ICH4FEED, METH1->ch4l.getval( J-1 ), EECH4, DENOX, DECO, DEVOC, METH1->ch4l.getptr( J ),
                    T00LO, &TAULO, SSLO, ANOXLO, ACOLO, AVOCLO, DELTAT, METH4 );
            // Note that TAULO has a value return in it, but is never used
            //F4498 !
            //F4499 !  MID (BEST) LIFETIME
            //F4500 !
            //F4501         EECH4 = ECH4(J)-CORRUSER+CORRMMID
            EECH4 = CONCS->ECH4.getval( J ) - METH1->CORRUSER + METH1->CORRMMID;
            //F4502 !
            //F4503         CALL METHANE(ICH4FEED,CH4B(J-1),EECH4,DENOX,DECO,DEVOC,CH4B(J), &
            //F4504         T00MID,TAUBEST,SCH4,ANOX,ACO,AVOC,DELTAT)
            float TAUBEST = 0.0;
            methane( METH3->ICH4FEED, METH1->ch4b.getval( J-1 ), EECH4, DENOX, DECO, DEVOC, METH1->ch4b.getptr( J ),
                    T00MID, &TAUBEST, METH3->SCH4, METH3->ANOX, METH3->ACO, METH3->AVOC, DELTAT, METH4 );
            // Note that TAUBEST has a value return in it, but is never used
            //F4505 !
            //F4506 !  HIGH LIFETIME
            //F4507 !
            //F4508         EECH4 = ECH4(J)-CORRUSER+CORRMHI
            EECH4 = CONCS->ECH4.getval( J ) - METH1->CORRUSER + METH1->CORRMHI;
            //F4509 !
            //F4510         SSHI=SCH4+DELSS
            float SSHI = METH3->SCH4 + METH3->DELSS;
            //F4511         ANOXHI=ANOX-DELANOX
            float ANOXHI = METH3->ANOX + METH3->DELACO;
            //F4512         ACOHI=ACO+DELACO
            float ACOHI = METH3->ACO + METH3->DELACO;
            //F4513         AVOCHI=AVOC+DELAVOC
            float AVOCHI = METH3->AVOC + METH3->DELAVOC;
            //F4514 !
            //F4515         CALL METHANE(ICH4FEED,CH4H(J-1),EECH4,DENOX,DECO,DEVOC,CH4H(J), &
            //F4516         T00HI,TAUHI,SSHI,ANOXHI,ACOHI,AVOCHI,DELTAT)
            float TAUHI = 0.0;
            methane( METH3->ICH4FEED, METH1->ch4h.getval( J-1 ), EECH4, DENOX, DECO, DEVOC, METH1->ch4h.getptr( J ),
                    T00HI, &TAUHI, SSHI, ANOXHI, ACOHI, AVOCHI, DELTAT, METH4 );
            // Note that TAUHI has a value return in it, but is never used
            //F4517 !
            //F4518 !  USER LIFETIME (ONE OF ABOVE, OR CONSTANT AT SPECIFIED 1990 VALUE)
            //F4519 !
            //F4520         EECH4 = ECH4(J)
            EECH4 = CONCS->ECH4.getval( J );
            //F4521 !
            //F4522 !  SET MODEL PARAMETERS FOR USER CASE
            //F4523 !
            //F4524         IF(LEVCH4.EQ.1)THEN
            float SSUSER, ANOXUSER, ACOUSER, AVOCUSER;
            if( METH2->LEVCH4 == 1 ) {
                //F4525           SSUSER=SCH4-DELSS
                SSUSER = METH3->SCH4 - METH3->DELSS;
                //F4526           ANOXUSER=ANOX+DELANOX
                ANOXUSER = METH3->ANOX + METH3->DELANOX;
                //F4527           ACOUSER=ACO-DELACO
                ACOUSER = METH3->ACO - METH3->DELACO;
                //F4528           AVOCUSER=AVOC-DELAVOC
                AVOCUSER = METH3->AVOC - METH3->DELAVOC;
                //F4529         ENDIF
            }
            //F4530 !
            //F4531         IF(LEVCH4.EQ.2)THEN
            if( METH2->LEVCH4 == 2 ) {
                //F4532           SSUSER=SCH4
                SSUSER = METH3->SCH4;
                //F4533           ANOXUSER=ANOX
                ANOXUSER = METH3->ANOX;
                //F4534           ACOUSER=ACO
                ACOUSER = METH3->ACO;
                //F4535           AVOCUSER=AVOC
                AVOCUSER = METH3->AVOC;
                //F4536         ENDIF
            }
            //F4537 !
            //F4538         IF(LEVCH4.EQ.3)THEN
            if( METH2->LEVCH4 == 3 ) {
                //F4539           SSUSER=SCH4+DELSS
                SSUSER = METH3->SCH4 + METH3->DELSS;
                //F4540           ANOXUSER=ANOX-DELANOX
                ANOXUSER = METH3->ANOX - METH3->DELANOX;
                //F4541           ACOUSER=ACO+DELACO
                ACOUSER = METH3->ACO + METH3->DELACO;
                //F4542           AVOCUSER=AVOC+DELAVOC
                AVOCUSER = METH3->AVOC + METH3->DELAVOC;
                //F4543         ENDIF
            }
            //F4544 !
            //F4545         IF(LEVCH4.EQ.4)THEN
            if( METH2->LEVCH4 == 4 ) {
                //F4546           SSUSER=0.0
                //F4547           ANOXUSER=0.0
                //F4548           ACOUSER=0.0
                //F4549           AVOCUSER=0.0
                SSUSER = ANOXUSER = ACOUSER = AVOCUSER = 0.0;
                //F4550         ENDIF
            }
            //F4551 !
            //F4552         CALL METHANE(ICH4FEED,CH4(J-1),EECH4,DENOX,DECO,DEVOC,CH4(J), &
            //F4553         T00USER,TAUCH4,SSUSER,ANOXUSER,ACOUSER,AVOCUSER,DELTAT)
            TAUCH4 = 0.0;
            methane( METH3->ICH4FEED, CONCS->CH4[ J-1 ], EECH4, DENOX, DECO, DEVOC, &CONCS->CH4[ J ],
                    T00USER, &TAUCH4, SSUSER, ANOXUSER, ACOUSER, AVOCUSER, DELTAT, METH4 );
            //F4554 !
            //F4555 !  SAVE USER-MODEL METHANE LIFETIME. TCH4(J) = CHEMICAL (OH)
            //F4556 !   LIFETIME. THIS IS THE SAME AS ......
            //F4557 !   TCH4EFF(J)=CH4BAR/(ECH4(J)/BBCH4-DELCH4-SOILSINK-STRATSINK)
            //F4558 !
            //F4559         TCH4(J)=TAUCH4
            METH1->TCH4[ J ] = TAUCH4;
            //F4560 !
            //F4561 ! Methane oxidation source: based on user methane projection only.
            //F4562 !  User projection determined by choice of LEVCH4 (1,2,3,or 4).
            //F4563 !
            //F4564         CH4BAR=(CH4(J-1)+CH4(J))/2
            float CH4BAR = ( CONCS->CH4[ J-1 ] + CONCS->CH4[ J ] ) / 2;
            //F4565         EMETH(J) = FFFRAC*0.0020625*(CH4BAR-700.)/TAUCH4
            METH1->emeth.setval( fffrac * 0.0020625 * ( CH4BAR-700.0 ) / TAUCH4, J );
            //F4566         IF(IMETH.EQ.1)THEN
            if( METH1->IMETH == 1 )
                //F4567           EF4(J) = EF(J)+EMETH(J)
                METH1->ef4.setval( CARB->EF.getval( J ) + METH1->emeth.getval( J ), J);
            //F4568         ELSE
            //F4569           EF4(J) = EF(J)
            else METH1->ef4.setval( CARB->EF.getval( J ), J );
            //F4570         ENDIF
            //F4571 !
            //F4572       ENDIF
        }
        //F4573 !
        //F4574 !  *******************************************************
        //F4575 !
        //F4576       IF(J.GT.JSTART)THEN
        if( J > JSTART->JSTART ) {
            //F4577 !
            //F4578 !  N2O
            //F4579 !  N2O CONCs
            //F4580 !
            //F4581 !  NOFFSET IS THE TIME IT TAKES FOR N2O TO MIX TO THE STRATOSPHERE
            //F4582 !
            //F4583         J1=J-1
            int J1 = J - 1;
            //F4584         J2=J-NOFFSET
            int J2 = J - TauNitr->NOFFSET;
            //F4585         J3=J-NOFFSET-1
            int J3 = J - TauNitr->NOFFSET - 1;
            //F4586 !
            //F4587         CALL NITROUS(CN2O(J1),CN2O(J2),CN2O(J3),EN2O(J),CN2O(J))
            nitrous( CONCS->CN2O[ J1 ], CONCS->CN2O[ J2 ], CONCS->CN2O[ J3 ], CONCS->EN2O.getval( J ), &CONCS->CN2O[ J ], TauNitr );
            //F4588 !
            //F4589       ENDIF
        }
        //F4590 !
        //F4591 !  *******************************************************
        //F4592 !
        //F4593 !  HALOCARBS, ETC.
        //F4594 !  CONC CALCULATIONS FROM 1991 ONWARDS FOR KYOTO PROTOCOL HFCs,
        //F4595 !   PFCs AND SF6. NOTE : FOR PRESENT VERSION, ALTHO CONCS ARE
        //F4596 !   CALCULATED, THEY ARE NOT OUTPUT. ONLY THE TOTAL FORCING IS USED.
        //F4597 !  HALO FORCINGS FOR 1990 AND EARLIER ARE GIVEN IN QHALOS.IN, BROKEN
        //F4598 !   DOWN INTO MONTREAL GASES, MAGICC KYOTO GASES, OTHER GASES AND 
        //F4599 !   STRAT OZONE. FOR 1991+, QHALOS.IN GIVES FORCINGS FOR MONTREAL
        //F4600 !   GASES, OTHER GASES AND STRAT OZONE.
        //F4601 !  EFFECT OF OH CHANGES ON LIFETIMES ADDED 000909
        //F4602 !
        //F4603        T0=TAUINIT
        float TM, T0 = METH3->TAUINIT;
        //F4604 !
        //F4605       IF(J.GT.226)THEN
        if( J > 226 ) {
            //F4606         IF(J.LE.JSTART)THEN
            if( J <= JSTART->JSTART ) {
                //F4607           TM=T0
                TM = T0;
                //F4608         ELSE
            } else {
                //F4609 !
                //F4610 !  PREVIOUS VERSION USED TAUBEST HERE. IT IS MORE CONSISTENT TO USE
                //F4611 !   TAUCH4, WHICH IS THE USER LIFETIME. THIS WILL GIVE HALOCARBON
                //F4612 !   LIFETIMES CONSISTENT WITH USER METHANE LIFETIME.
                //F4613 !
                //F4614           TM=TAUCH4
                TM = TAUCH4;
                //F4615         ENDIF
            }
            //F4616 !
            float QCF4 = 0.0, QC2F6 = 0.0, Q125 = 0.0, Q134A = 0.0, Q143A = 0.0, 
            Q227 = 0.0, Q245 = 0.0, QSF6 = 0.0;
            //F4617         CALL HALOCARB(1,CF4(J-1),  ECF4(J), CF4(J)  ,QCF4 ,T0,TM)
            halocarb( 1, NEWCONCS->CF4[ J-1 ], NEWCONCS->ECF4.getval( J ), &NEWCONCS->CF4[ J ], &QCF4, T0, TM );
            //F4618         CALL HALOCARB(2,C2F6(J-1), EC2F6(J),C2F6(J) ,QC2F6,T0,TM)
            halocarb( 2, NEWCONCS->C2F6[ J-1 ], NEWCONCS->EC2F6.getval( J ), &NEWCONCS->C2F6[ J ], &QC2F6, T0, TM );
            //F4619         CALL HALOCARB(3,C125(J-1), E125(J), C125(J) ,Q125 ,T0,TM)
            halocarb( 3, NEWCONCS->C125[ J-1 ], NEWCONCS->E125.getval( J ), &NEWCONCS->C125[ J ], &Q125, T0, TM );
            //F4620         CALL HALOCARB(4,C134A(J-1),E134A(J),C134A(J),Q134A,T0,TM)
            halocarb( 4, NEWCONCS->C134A[ J-1 ], NEWCONCS->E134A.getval( J ), &NEWCONCS->C134A[ J ], &Q134A, T0, TM );
            //F4621         CALL HALOCARB(5,C143A(J-1),E143A(J),C143A(J),Q143A,T0,TM)
            halocarb( 5, NEWCONCS->C143A[ J-1 ], NEWCONCS->E143A.getval( J ), &NEWCONCS->C143A[ J ], &Q143A, T0, TM );
            //F4622         CALL HALOCARB(6,C227(J-1), E227(J), C227(J) ,Q227 ,T0,TM)
            halocarb( 6, NEWCONCS->C227[ J-1 ], NEWCONCS->E227.getval( J ), &NEWCONCS->C227[ J ], &Q227, T0, TM );
            //F4623         CALL HALOCARB(7,C245(J-1), E245(J), C245(J) ,Q245 ,T0,TM)
            halocarb( 7, NEWCONCS->C245[ J-1 ], NEWCONCS->E245.getval( J ), &NEWCONCS->C245[ J ], &Q245, T0, TM );
            //F4624         CALL HALOCARB(8,CSF6(J-1), ESF6(J), CSF6(J) ,QSF6 ,T0,TM)
            halocarb( 8, NEWCONCS->CSF6[ J-1 ], NEWCONCS->ESF6.getval( J ), &NEWCONCS->CSF6[ J ], &QSF6, T0, TM );
            //F4625 !
            //F4626         QKYMAG(J)=QCF4+QC2F6+Q125+Q134A+Q143A+Q227+Q245+QSF6
            JSTART->QKYMAG[ J ] = QCF4 + QC2F6 + Q125 + Q134A + Q143A + Q227 + Q245 + QSF6;
            //F4627 !
            //F4628 ! sjs -- Save halocarbon forcing
            //F4629 
            //F4630 	QCF4_ar(J) = QCF4
            HALOF->QCF4_ar[ J ] = QCF4;
            //F4631 	QC2F6_ar(J) = QC2F6
            HALOF->qSF6_ar[ J ] = QC2F6;
            //F4632 	qSF6_ar(J) = QSF6
            HALOF->qSF6_ar[ J ] = QSF6;
            //F4633 	Q125_ar(J) = Q125
            HALOF->Q125_ar[ J ] = Q125;
            //F4634 	Q134A_ar(J) = Q134A
            HALOF->Q134A_ar[ J ] = Q134A;
            //F4635 	Q143A_ar(J) = Q143A
            HALOF->Q143A_ar[ J ] = Q143A;
            //F4636 	Q227_ar(J) = Q227
            HALOF->Q227_ar[ J ] = Q227;
            //F4637 	Q245_ar(J) = Q245
            HALOF->Q245_ar[ J ] = Q245;
            //F4638 
            //F4639       ENDIF
        }
        //F4640 !
        //F4641 !  *******************************************************
        //F4642 !
        //F4643 !  CO2
        //F4644 !  CARBON CYCLE MODEL CALL. NOTE THAT THIS IS CALLED FROM 1991
        //F4645 !   ONWARDS, IRRESPECTIVE OF JSTART VALUE
        //F4646 !
        //F4647       IF(J.GT.226)THEN
        if( J > 226 ) {
            //F4648 !
            //F4649 !  NC=1, BCO2 BASED ON D80S=1.8 TO GIVE LOWER BOUND CONCS
            //F4650 !  NC=2, BCO2 BASED ON D80S=1.1 TO GIVE BEST GUESS CONCS
            //F4651 !  NC=3, BCO2 BASED ON D80S=0.4 TO GIVE UPPER BOUND CONCS
            //F4652 !  NC=4, BCO2 BASED ON USER SELECTED D80S
            //F4653 !  ALL CASES USE F80S =2.0
            //F4654 !
            //F4655 !  BEST GUES CHANGED TO 1.5 (MAY 2008)
            //F4656 !
            //F4657         DO 444 NC = 1,4
            for( int NC = 1; NC <=4; NC++ ) {
                //F4658 !
                //F4659 !  CALL INITCAR TO INITIALIZE CARBON CYCLE MODEL
                //F4660 !
                //F4661         IF(J.EQ.227)THEN
                if( J == 227 ) {
                    //F4662           FIN=2.0
                    float FIN = 2.0;
                    //F4663           DIN=2.9-0.7*NC
                    float DIN = 2.9 - 0.7 * NC;
                    //F4664 !         DIN=2.5-0.7*NC
                    //F4665           IF(NC.EQ.4)THEN
                    if( NC == 4 ) {
                        //F4666             FIN=FUSER
                        FIN = METH1->FUSER;
                        //F4667             DIN=DUSER
                        DIN = METH1->DUSER;
                        //F4668           ENDIF
                    }
                    //F4669           CALL INITCAR(NC,DIN,FIN)
                    initcar( NC, DIN, FIN, COBS, CARB, CAR );
                    //F4670         ENDIF
                }
                //F4671 !
                //F4672 !  IF IDRELAX.NE.O, OVERWRITE EDNET(J) FOR 1990 TO 1990+DRELAX WITH
                //F4673 !   A LINEAR INTERP BETWEEN THE 1990 VALUE BASED ON BALANCING THE
                //F4674 !   1980S-MEAN CARBON BUDGET TO THE APPROPRIATE VALUE OBTAINED FROM
                //F4675 !   THE EMISSIONS INPUT FILE.
                //F4676 !
                //F4677         IDRELAX=10
                const int IDRELAX = 10;
                //F4678         IF(IDRELAX.NE.0)THEN
                if( IDRELAX != 0 ) {
                    //F4679           EDNET(226)=EDNET90(NC)
                    METH1->ednet.setval( CARB->EDNET90[ NC ], 226 );
                    //F4680           JDEND=226+IDRELAX
                    const int JDEND = 226 + IDRELAX;
                    //F4681           DELED=(EDNET(JDEND)-EDNET(226))/FLOAT(IDRELAX)
                    const float DELED = ( METH1->ednet.getval( JDEND ) - METH1->ednet.getval( 226 ) ) / float( IDRELAX );
                    //F4682           DO JD=227,JDEND-1
                    for( int JD=227; JD<=JDEND-1; JD++) {
                        //F4683           EDNET(JD)=EDNET(226)+DELED*(JD-226)
                        METH1->ednet.setval( METH1->ednet.getval( 226 ) + DELED * ( JD-226 ), JD );
                    }
                    //F4684           END DO
                    //F4685         ENDIF
                }
                //F4686 !
                //F4687 !  Note: for temp feedback on CO2, temp from default or user carbon cycle
                //F4688 !        model is used in all four nc cases.
                //F4689 !        Strictly in (eg) the upper bound case one should
                //F4690 !        use the corresponding temp. However the upper bound CO2 is not
                //F4691 !        generally passed to the climate model so this temp is not
                //F4692 !        calculated. The error in making this approx must be small
                //F4693 !        as it is only a second order effect.
                //F4694 !  Note: this also applies to the methane model.
                //F4695 !
                //F4696         TEMP=TX-DELT90
                const float TEMP = TX - DELT90;
                //F4697 !
                //F4698         CALL CARBON(NC,TEMP,EF4(J),EDNET(J),CCO2(NC,J-3),CCO2(NC,J-2), &
                //F4699         CCO2(NC,J-1), &
                //F4700         PL(NC,J-1),HL(NC,J-1),SOIL(NC,J-1),REGROW(NC,J-1),ETOT(NC,J-1), &
                //F4701         PL(NC,J)  ,HL(NC,J)  ,SOIL(NC,J)  ,REGROW(NC,J)  ,ETOT(NC,J)  , &
                //F4702         ESUM(J),FOC(NC,J),DELMASS(NC,J),EDGROSS(NC,J),CCO2(NC,J))
                carbon( NC, TEMP, METH1->ef4.getval( J ), METH1->ednet.getval( J ), CARB->CCO2.getval( NC, J-3 ), CARB->CCO2.getval( NC, J-2 ),
                       CARB->CCO2.getval( NC, J-1 ),
                       CARB->PL.getval( NC, J-1 ), CARB->HL.getval( NC, J-1 ), CARB->SOIL.getval( NC, J-1 ),  CARB->REGROW.getval( NC, J-1 ),  CARB->ETOT.getval( NC, J-1 ),
                       CARB->PL.getptr( NC, J ), CARB->HL.getptr( NC, J ), CARB->SOIL.getptr( NC, J ),  CARB->REGROW.getptr( NC, J ),  CARB->ETOT.getptr( NC, J ),
                       CARB->ESUM.getptr( J ), CARB->FOC.getptr( NC, J ), CAR->DELMASS.getptr( NC, J ), CARB->EDGROSS.getptr( NC, J ), CARB->CCO2.getptr( NC, J ),
                       CAR );
                //F4703 !
                //F4704   444   CONTINUE
            } // for
            //F4705 !
            //F4706 !  SELECT CO2 VALUES (CO2(J)) TO CARRY ON TO FORCING :
            //F4707 !   THE PARTICULAR CARBON CYCLE MODEL OUTPUT THAT IS CARRIED ON IS
            //F4708 !   DETERMINED BY THE SPECIFIED VALUE OF LEVCO2.
            //F4709 !  NOTE THAT, IF ICO2CORR=1 (THE DEFAULT VALUE) ALL CO2 ARRAYS ARE
            //F4710 !   CORRECTED (CHANGE MADE ON 6/10/03). THE ARRAYS ARE CORRECTED TO
            //F4711 !   AGREE WITH OBSERVATIONS THROUGH JSTART. THE ARRAY CO2(J) HAS
            //F4712 !   ALREADY BEEN SPECIFIED AS HISTORICAL OBSERVED DATA THROUGH
            //F4713 !   J=JSTART. WHEN J=JSTART A CORRECTION FACTOR IS CALCULATED AND
            //F4714 !   THIS IS APPLIED TO ALL SUBSEQUENT YEARS.
            //F4715 !
            //F4716         IF(J.EQ.JSTART)THEN
            if( J == JSTART->JSTART ) {
                //F4717           CORREN1=CO2(J)-CCO2(1,J)
                CORREN->CORREN1 = CARB->CO2[ J ] - CARB->CCO2.getval( 1, J );
                //F4718           CORREN2=CO2(J)-CCO2(2,J)
                CORREN->CORREN2 = CARB->CO2[ J ] - CARB->CCO2.getval( 2, J );
                //F4719           CORREN3=CO2(J)-CCO2(3,J)
                CORREN->CORREN3 = CARB->CO2[ J ] - CARB->CCO2.getval( 3, J );
                //F4720           CORREN4=CO2(J)-CCO2(4,J)
                CORREN->CORREN4 = CARB->CO2[ J ] - CARB->CCO2.getval( 4, J );
                //F4721           CORREN=CO2(J)-CCO2(LEVCO2,J)
                CORREN->CORREN = CARB->CO2[ J ] - CARB->CCO2.getval( CO2READ->LEVCO2, J );
                //F4722         ENDIF
            }
            //F4723 !
            //F4724         IF(J.GE.JSTART)THEN
            if( J >= JSTART->JSTART ) {
                //F4725           CO2(J)=CCO2(LEVCO2,J)
                CARB->CO2[ J ] = CARB->CCO2.getval( CO2READ->LEVCO2, J );
                //F4726           IF(ICO2CORR.EQ.1)THEN
                if( JSTART->ICO2CORR == 1 )
                    //F4727             CO2(J)=CCO2(LEVCO2,J)+CORREN
                    CARB->CO2[ J ] = CARB->CCO2.getval( CO2READ->LEVCO2, J ) + CORREN->CORREN;
                //F4728           ENDIF
                //F4729         ENDIF
            }
            //F4730 !
            //F4731 !  FEEDBACK PERCENTILE VALUES
            //F4732 !
            //F4733 !       IF(J.GE.JSTART)THEN
            //F4734 !         CO2(J)=CO2(J)*(1.0-0.0855*(J-236)/100.0)          !! 10%
            //F4735 !         CO2(J)=CO2(J)*(1.0-0.0414*(J-236)/100.0)          !! 30%
            //F4736 !         CO2(J)=CO2(J)*(1.0+0.0650*((J-236)/100.0)**1.5)   !! 70%
            //F4737 !         CO2(J)=CO2(J)*(1.0+0.2340*((J-236)/100.0)**1.5)   !! 90%
            //F4738 !       ENDIF
            //F4739 !
            //F4740 !  OVERWRITE CARBON CYCLE CO2 CONCS FOR YEARS.GE.1990 WITH INPUT
            //F4741 !   DATA FROM CO2INPUT.DAT IF ICO2READ.GE.1.
            //F4742 !
            //F4743         IF(ICO2READ.GE.1.AND.ICO2READ.LE.4) co2(J)=xc(J)
            if( CO2READ->ICO2READ >= 1 && CO2READ->ICO2READ <= 4 ) CARB->CO2[ J ] = CO2READ->XC[ J ];
            //F4744 !
            //F4745       ENDIF
        } // if
        //F4746 !
        //F4747 !      NSAVE=4
        //F4748 !      IF(ISCENGEN.EQ.9)NSAVE=1
        //F4749 !      IF(NSIM.EQ.NSAVE)THEN
        //F4750       IF(NSIM.LE.4.AND.NCLIM.EQ.4)THEN
        if( NSIM->NSIM <= 4 && NSIM->NCLIM == 4 ) {
            //F4751         CO2SAVE(J)=CO2(J)
            CARB->CO2SAVE[ J ] = CARB->CO2[ J ];
            //F4752       ENDIF
        }
        //F4753 !
        //F4754 !  **************************************************
        //F4755 !
        //F4756 !  END OF LONG SEQUENCE OF IF STATEMENTS
        //F4757 !
        //F4758 !  *******************************************************
        //F4759 !
        //F4760 !  HALOCARBON FORCING (IOLDHALO OPTION DELETED, 8 OCT 2000)
        //F4761 !
        //F4762       QCFC(J)=QMONT(J)+QKYMAG(J)+QOTHER(J)+QSTRATOZ(J)
        FORCE->QCFC[ J ] = FORCE->QMONT[ J ] + JSTART->QKYMAG[ J ] + FORCE->QOTHER[ J ] + FORCE->QSTRATOZ[ J ];
        //F4763 !
        //F4764 !  SUBTRACT STRAT OZONE FORCING IF IO3FEED=0 (I.E., NFB CASE)
        //F4765 !
        //F4766       IF(IO3FEED.EQ.0)QCFC(J)=QCFC(J)-QSTRATOZ(J)
        if( METH1->IO3FEED == 0 ) FORCE->QCFC[ J ] = FORCE->QCFC[ J ] - FORCE->QSTRATOZ[ J ];
        //F4767 !
        //F4768 !  *******************************************************
        //F4769 !
        //F4770 !   METHANE FORCING
        //F4771 !
        //F4772       QCH4=0.036*(SQRT(CH4(J))-SQRT(CH4(1)))
        const float QCH4 = 0.036 * ( sqrt( CONCS->CH4[ J ] )-sqrt( CONCS->CH4[ 1 ] ) );
        //F4773 !
        //F4774       XM=CH4(J)/1000.
        float XM = CONCS->CH4[ J ] / 1000.0;
        //F4775       WW=CN2O(1)/1000.
        float WW = CONCS->CN2O[ 1 ] / 1000.0;
        //F4776       AB=0.636*((XM*WW)**0.75) + 0.007*XM*((XM*WW)**1.52)
        float AB = 0.636 * pow( XM*WW, static_cast<float> (0.75) ) + 0.007 * XM * pow( XM*WW, static_cast<float> (1.52) );
        //F4777 !
        //F4778       XM0=CH4(1)/1000.
        const float XM0 = CONCS->CH4[ 1 ] / 1000.0;
        //F4779       AB0=0.636*((XM0*WW)**0.75) + 0.007*XM0*((XM0*WW)**1.52)
        const float AB0 = 0.636 * pow( XM0*WW, static_cast<float> (0.75) ) + 0.007 * XM0 * pow( XM0*WW, static_cast<float> (1.52) );
        //F4780 !
        //F4781       QMeth=QCH4+0.47*ALOG((1.+AB0)/(1.+AB))
        const float QMeth = QCH4 + 0.47 * log( float( ( 1.+AB0 ) / ( 1.+AB ) ) );
        //F4782 !
        //F4783 !  QH2O IS THE ADDITIONAL INDIRECT CH4 FORCING DUE TO PRODUCTION
        //F4784 !   OF STRAT H2O FORM CH4 OXIDATION.  QCH4OZ IS THE ENHANCEMENT
        //F4785 !   OF QCH4 DUE TO CH4-INDUCED OZONE PRODUCTION.
        //F4786 !  THE DENOMINATOR IN QCH4OZ HAS BEEN CHOSEN TO GIVE 0.08W/m**2
        //F4787 !   FORCING IN MID 1990. IT DIFFERS SLIGHTLY FROM THE VALUE
        //F4788 !   ORIGINALLY ADVISED BY PRATHER (WHICH WAS 353). THE 0.0025
        //F4789 !   CORRECTION TERM IS TO MAKE THE CHANGE RELATIVE TO MID 1765.
        //F4790 !
        //F4791       QH2O=STRATH2O*QCH4
        const float QH2O = METH1->STRATH2O * QCH4;
        //F4792       
        //F4793 ! Save this value to Kyoto Forcing can be calculated
        //F4794       QCH4H2O(J) = QH2O
        FORCE->QCH4H2O[ J ] = QH2O;
        //F4795 !
        //F4796 !  QCH4OZ IS THE TROP OZONE FORCING FROM CH4 CHANGES.
        //F4797 !   HISTORY : USE THE TAR LOGARITHMIC RELATIONSHIP, SCALED TO
        //F4798 !    OZ00CH4 (THE CENTRAL TAR CH4-RELATED FORCING IN 2000).
        //F4799 !
        //F4800       IF(J.LE.235)THEN
        if( J <= 235 ) {
            //F4801         AAA=OZ00CH4/ALOG(1760./700.)
            const float AAA = OZ->OZ00CH4 / log( float( 1760.0 / 700.0 ) );
            //F4802         QCH4OZ=AAA*ALOG(CH4(J)/700.0)
            JSTART->QCH4OZ = AAA * log( float( CONCS->CH4[ J ] / 700.0 ) );
            //F4803       ELSE
        } else {
            //F4804         QCH4OZ=OZ00CH4+TROZSENS*OZCH4*ALOG(CH4(J)/CH4(235))
            JSTART->QCH4OZ = OZ->OZ00CH4 + JSTART->TROZSENS * OZ->OZCH4 * log( float( CONCS->CH4[ J ] / CONCS->CH4[235] ) );
            //F4805       ENDIF
        }
        //F4806 !
        //F4807       QM(J) = qMeth+qH2O+QCH4OZ
        FORCE->QM[ J ] = QMeth + QH2O + JSTART->QCH4OZ;
        //F4808       QCH4O3(J)=QCH4OZ
        FORCE->QCH4O3[ J ] = JSTART->QCH4OZ;
        //F4809 !
        //F4810 !  *******************************************************
        //F4811 !
        //F4812 !  TROPOSPHERIC OZONE FORCING NOT ASSOCIATED WITH CH4.
        //F4813 !   HISTORY : SCALE WITH FOSSIL CO2 EMISSIONS AS A PROXY FOR
        //F4814 !    THE INFLUENCE OF THE REACTIVE GASES. EMISSIONS HISTORY
        //F4815 !    IS SMOOTHED VERSION OF MARLAND'S DATA TO 1990, WITH
        //F4816 !    SRES TO 2000. SCALING FACTOR CHOSEN TO MAKE TOTAL TROP
        //F4817 !    OZONE FORCING = 0.35 W/m**2 AT BEGINNING OF 2000.
        //F4818 !   FUTURE : USE TAR RELATIONSHIP AND REACTIVE GAS EMISSIONS.
        //F4819 !
        //F4820       QREF=0.33-OZ00CH4
        const float QREF = 0.33 - OZ->OZ00CH4;
        //F4821       IF(J.LE.235)THEN
        if( J <= 235 ) {
            //F4822         FOSS0=FOSSHIST(1)
            const float FOSS0 = JSTART->FOSSHIST[ 1 ];
            //F4823         QOZ(J)=QREF*(FOSSHIST(J)-FOSS0)/(FOSSHIST(235)-FOSS0)
            TANDSL->QOZ[ J ] = QREF * ( JSTART->FOSSHIST[ J ] - FOSS0 ) / ( JSTART->FOSSHIST[ 235 ] - FOSS0 );
            //F4824         IF(J.EQ.234)QOZ1=QOZ(J)
            if( J == 234 ) QOZ1 = TANDSL->QOZ[ J ];
            //F4825         IF(J.EQ.235)DQOZ=QOZ(J)-QOZ1
            if( J == 235 ) DQOZ = TANDSL->QOZ[ J ] - QOZ1;
            //F4826       ELSE
        } else {
            //F4827         DDEN=ENOX(J)-ENOX(236)
            const float DDEN = CONCS->ENOX.getval( J ) - CONCS->ENOX.getval( 236 );
            //F4828         DDEC=ECO(J) -ECO(236)
            const float DDEC = CONCS->ECO.getval( J ) - CONCS->ECO.getval( 236 );
            //F4829         DDEV=EVOC(J)-EVOC(236)
            const float DDEV = CONCS->EVOC.getval( J ) - CONCS->EVOC.getval( 236 );
            //F4830         QOZ(J)=QREF+DQOZ &
            //F4831         +TROZSENS*(OZNOX*DDEN+OZCO*DDEC+OZVOC*DDEV)
            TANDSL->QOZ[ J ] = QREF + DQOZ + JSTART->TROZSENS * ( OZ->OZNOX * DDEN + OZ->OZCO * DDEC + OZ->OZVOC * DDEV );
            //F4832       ENDIF
        }
        //F4833 !
        //F4834       IF(IGHG.LE.1)QOZ(J)=0.0
        if( JSTART->IGHG <= 1 ) TANDSL->QOZ[ J ] = 0.0;
        //F4835 !
        //F4836 !  *******************************************************
        //F4837 !
        //F4838 !   NITROUS OXIDE FORCING
        //F4839 !
        //F4840       QN2O=QQQN2O*(SQRT(CN2O(J))-SQRT(CN2O(1)))
        const float QN2O = METH2->QQQN2O * ( sqrt( CONCS->CN2O[ J ] ) - sqrt( CONCS->CN2O[ 1 ] ) );
        //F4841 !
        //F4842       XM=XM0
        XM = XM0;
        //F4843       WW=CN2O(J)/1000.
        WW = CONCS->CN2O[ J ] / 1000.0;
        //F4844       AB=0.636*((XM*WW)**0.75) + 0.007*(XM*(XM*WW)**1.52)
        AB = 0.636 * pow( XM*WW, static_cast<float> (0.75) ) + 0.007 * XM * pow( XM*WW, static_cast<float> (1.52) );
        //F4845       QN(J)=QN2O+0.47*ALOG((1.+AB0)/(1.+AB))
        FORCE->QN[ J ] = QN2O + 0.47 * log( float( ( 1.+AB0 ) / ( 1.0+AB ) ) );
        //F4846 !
        //F4847 !  *******************************************************
        //F4848 !
        //F4849 !   TOTAL GREENHOUSE FORCING
        //F4850 !
        //F4851       QCO2(J)=QXX*ALOG(CO2(J)/CO2(1))
        FORCE->QCO2[ J ] = CLIM->QXX * log( float( CARB->CO2[ J ] / CARB->CO2[ 1 ] ) );
        //        if(J>235) cout << "QCO2 " << J << " " << FORCE->QCO2[J] << " " << CLIM->QXX << " " << CARB->CO2[J] << " " << CARB->CO2[1] << endl; //debug
        //F4852       QGH(J)=QCO2(J)+QM(J)+QN(J)+QCFC(J)
        TANDSL->QGH[ J ] = FORCE->QCO2[ J ] + FORCE->QM[ J ] + FORCE->QN[ J ] + FORCE->QCFC[ J ];
        //        if(J>235) cout << "QGH " << J << " " << TANDSL->QGH[J] << " " << FORCE->QCO2[J] << " " << FORCE->QM[J] << " " << FORCE->QN[J] << " " << FORCE->QCFC[J] << endl; //debug
        //F4853 !
        //F4854       IF(IGHG.EQ.0)QGH(J)=0.0
        if( JSTART->IGHG == 0 ) TANDSL->QGH[ J ] = 0.0;
        //F4855       IF(IGHG.EQ.1)QGH(J)=QCO2(J)
        if( JSTART->IGHG == 1 ) TANDSL->QGH[ J ] = FORCE->QCO2[ J ];
        //F4856 !
        //F4857 !  *******************************************************
        //F4858 !
        //F4859 !   CALCULATE BIOMASS BURNING AEROSOL TERM (SUM OF ORGANIC PLUS
        //F4860 !    BLACK CARBON).
        //F4861 !   HISTORY : ASSUME LINEAR RAMP UP OVER 1765-1990. THIS VERY
        //F4862 !    ROUGHLY FOLLOWS THE GROSS DEFORESTATION HISTORY
        //F4863 !   FUTURE : (1991 ONWARDS) SCALED WITH GROSS DEFORESTATION
        //F4864 !    OUTPUT FROM CARBON CYCLE CALCULATIONS, BUT CONSTRAINED
        //F4865 !    NEVER TO GO ABOVE ZERO. (NOTE THAT WITH THE LINEAR
        //F4866 !    RELATIONSHIP THIS COULD HAPPEN IF EDGROSS BECAME LESS
        //F4867 !    THAN ZERO, WHICH CAN HAPPEN WITH THE PRESENT CARBON CYCLE
        //F4868 !    MODEL. NOTE TOO THAT EDGROSS LESS THAN ZERO IS POSSIBLE:
        //F4869 !    THIS CORRRESPONDS TO NET REFORESTATION. CLEARLY, HOWEVER,
        //F4870 !    THIS WOULD NOT LEAD TO POSITIVE QBIO FORCING.)
        //F4871 !
        //F4872       IF(J.LE.226)QBIO(J)=S90BIO*FLOAT(J)/226.0
        if( J <= 226 ) TANDSL->QBIO[ J ] = Sulph->S90BIO * float( J ) / 226.0;
        //F4873       IF(J.GT.226)THEN
        if( J > 226 ) {
            //F4874         EDG=EDGROSS(LEVCO2,J)
            float EDG = CARB->EDGROSS.getval( CO2READ->LEVCO2, J );
            //F4875         IF(EDG.LT.0.0)EDG=0.0
            if( EDG < 0.0 ) EDG = 0.0;
            //F4876         QBIO(J)=S90BIO*EDG/EDGROSS(LEVCO2,226)
            TANDSL->QBIO[ J ] = Sulph->S90BIO * EDG / CARB->EDGROSS.getval( CO2READ->LEVCO2, 226 );
            //F4877 ! sjs - Check if EDGROSS(1990) is < 0, if so, ramp 1990 forcing down to zero by 2050
            //F4878         IF(EDGROSS(LEVCO2,226) .LT. 0) THEN
            if( CARB->EDGROSS.getval( CO2READ->LEVCO2, 226 ) < 0.0 ) {
                //F4879            QBIO(J)=S90BIO*(1.0 - (FLOAT(J)-226)/60.0)
                TANDSL->QBIO[ J ] = Sulph->S90BIO * ( 1.0 - ( float( J ) - 226 ) / 60.0 );
                //F4880            IF (QBIO(J) .GT. 0 ) QBIO(J) = 0
                if( TANDSL->QBIO[ J ] > 0.0 ) TANDSL->QBIO[ J ] = 0.0;
                //F4881         ENDIF
            }
            //F4882       ENDIF
        }
        //F4883 
        //F4884 !
        //F4885 !  *******************************************************
        //F4886 !
        //F4887 !   GLOBAL SULPHATE, FOSSIL ORGANIC AND FOSSIL BLACK CARBON
        //F4888 !    AEROSOL FORCING.
        //F4889 !   NOT : BOTH QSO2 AND QDIR OUTPUTS INCLUDE QFOC
        //F4890 !
        //F4891       CALL SULPHATE(J,ESO2(J),ESO2(J+1),COE(J),QSO2(J),QDIR(J),QFOC(J), &
        //F4892       QMN(J))
        sulphate( J, CONCS->ESO2.getval( J ), CONCS->ESO2.getval( J+1 ), CONCS->COE.getval( J ),
                 &TANDSL->QSO2[ J ], &TANDSL->QDIR[ J ], &JSTART->QFOC[ J ], &TANDSL->QMN[ J ], Sulph );
        //F4893 !
        //F4894 !   TOTAL GLOBAL FORCING INCLUDING AEROSOLS. NOTE THAT QSO2(J)
        //F4895 !    ALREADY INCLUDES QFOC(J)
        //F4896 !
        //F4897       qtot(J) = QGH(J)+qso2(J)+qoz(J)+QBIO(J)+qland(j)+QMN(J)
        TANDSL->QTOT[ J ] = TANDSL->QGH[ J ] + TANDSL->QSO2[ J ] + TANDSL->QOZ[ J ] + TANDSL->QBIO[ J ] + TANDSL->QLAND[ J ] + TANDSL->QMN[ J ];
        //F4898 !
        //F4899       if(j.eq.86)qtot86=qtot(J)
        if( J == 86 ) CO2READ->qtot86 = TANDSL->QTOT[ J ];
        //F4900 !
        //F4901 !  ***************************************************************
        //F4902 !
        //F4903 !  SWITCH TO OVERWRITE CO2 CONCS CALCULATED BY  CARBON CYCLE MODEL.
        //F4904 !   ICO2READ=1, THEN DELTA-QTOT FROM 1990 IS A DIRECT MULTIPLE
        //F4905 !    OF THE CO2 FORCING WITH SCALING FACTOR FROM CFG FILE.
        //F4906 !   ICO2READ=2, QTOT IS THE SUM OF THE NEW CO2 FORCING AND OTHER
        //F4907 !    FORCINGS AS DETERMINED BY THE GAS.EMK EMISSIONS INPUTS.
        //F4908 !   ICO2READ=3, QTOT IS THE SUM OF THE NEW CO2 FORCING AND SULPHATE
        //F4909 !    FORCING AS DETERMINED BY THE EMISSIONS INPUTS.
        //F4910 !
        //F4911       IF(ICO2READ.EQ.1.AND.J.GT.226)THEN
        if( CO2READ->ICO2READ == 1 && J > 226 ) {
            //F4912         QTOT(J) =QTOT(226)+CO2SCALE*(QXX*ALOG(CO2(J)/CO2(226)))
            TANDSL->QTOT[ J ] = TANDSL->QTOT[ 226 ] + CO2READ->CO2SCALE * ( CLIM->QXX * log( float( CARB->CO2[ J ]/CARB->CO2[ 226 ] ) ) );
            //F4913         QSO2(J) =QSO2(226)
            TANDSL->QSO2[ J ] = TANDSL->QSO2[ 226 ];
            //F4914         QDIR(J) =QDIR(226)
            TANDSL->QDIR[ J ] = TANDSL->QDIR[ 226 ];
            //F4915         QOZ(J)  =QOZ(226)
            TANDSL->QOZ[ J ] = TANDSL->QOZ[ 226 ];
            //F4916         QBIO(J) =QBIO(226)
            TANDSL->QBIO[ J ] = TANDSL->QBIO[ 226 ];
            //F4917         QLAND(J)=QLAND(226)
            TANDSL->QLAND[ J ] = TANDSL->QLAND[ 226 ];
            //F4918         QMN(J)  =QMN(226)
            TANDSL->QMN[ J ] = TANDSL->QMN[ 226 ];
            //F4919         QGH(J)  =QGH(226)+CO2SCALE*(QXX*ALOG(CO2(J)/CO2(226)))
            TANDSL->QGH[ J ] = TANDSL->QGH[ 226 ] + CO2READ->CO2SCALE * ( CLIM->QXX * log( float( CARB->CO2[ J ]/CARB->CO2[ 226 ] ) ) );
            //F4920       ENDIF
        }
        //F4921 !
        //F4922 !  IF ICO2READ=2, THE CHANGE REQUIRED IS ONLY FOR CO2 AND THIS HAS
        //F4923 !   ALREADY BEEN MADE.
        //F4924 !
        //F4925       IF((ICO2READ.EQ.3.OR.ICO2READ.EQ.4).AND.J.GT.226)THEN
        if( ( CO2READ->ICO2READ == 3 || CO2READ->ICO2READ == 4 ) && J > 226 ) {
            //F4926         QTOT(J) =QTOT(226)+CO2SCALE*(QXX*ALOG(CO2(J)/CO2(226)))
            TANDSL->QTOT[ J ] = TANDSL->QTOT[ 226 ] + CO2READ->CO2SCALE * ( CLIM->QXX * log( float( CARB->CO2[ J ]/CARB->CO2[ 226 ] ) ) );
            //F4927         QTOT(J) =QTOT(J)+(QSO2(J)-QSO2(226))
            TANDSL->QTOT[ J ] = TANDSL->QTOT[ J ] + ( TANDSL->QSO2[ J ] - TANDSL->QSO2[ 226 ] );
            //F4928         QOZ(J)  =QOZ(226)
            TANDSL->QOZ[ J ] = TANDSL->QOZ[ 226 ];
            //F4929         QBIO(J) =QBIO(226)
            TANDSL->QBIO[ J ] = TANDSL->QBIO[ 226 ];
            //F4930         QLAND(J)=QLAND(226)
            TANDSL->QLAND[ J ] = TANDSL->QLAND[ 226 ];
            //F4931         QMN(J)  =QMN(226)
            TANDSL->QMN[ J ] = TANDSL->QMN[ 226 ];
            //F4932         QGH(J)  =QGH(226)+CO2SCALE*(QXX*ALOG(CO2(J)/CO2(226)))
            TANDSL->QGH[ J ] = TANDSL->QGH[ 226 ] + CO2READ->CO2SCALE * ( CLIM->QXX * log( float( CARB->CO2[ J ]/CARB->CO2[ 226 ] ) ) );
            //F4933       ENDIF
        }
        //F4934 !
        //F4935       IF(ICO2READ.EQ.5.AND.J.GT.236)THEN
        if( CO2READ->ICO2READ == 5 && J > 236 ) {
            //F4936         QTOT(J) =QTOT(236)+QXX*ALOG(CO2(J)/CO2(236))
            TANDSL->QTOT[ J ] = TANDSL->QTOT[ 236 ] + CLIM->QXX * log( float( CARB->CO2[ J ]/CARB->CO2[ 236 ] ) );
            //F4937         QSO2(J) =QSO2(236)
            TANDSL->QSO2[ J ] = TANDSL->QSO2[ 236 ];
            //F4938         QDIR(J) =QDIR(236)
            TANDSL->QDIR[ J ] = TANDSL->QDIR[ 236 ];
            //F4939         QOZ(J)  =QOZ(236)
            TANDSL->QOZ[ J ] = TANDSL->QOZ[ 236 ];
            //F4940         QBIO(J) =QBIO(236)
            TANDSL->QBIO[ J ] = TANDSL->QBIO[ 236 ];
            //F4941         QLAND(J)=QLAND(236)
            TANDSL->QLAND[ J ] = TANDSL->QLAND[ 236 ];
            //F4942         QMN(J)  =QMN(236)
            TANDSL->QMN[ J ] = TANDSL->QMN[ 236 ];
            //F4943         QGH(J)  =QGH(236)+QXX*ALOG(CO2(J)/CO2(236))
            TANDSL->QGH[ J ] = TANDSL->QGH[ 236 ] + CO2READ->CO2SCALE * ( CLIM->QXX * log( float( CARB->CO2[ J ]/CARB->CO2[ 236 ] ) ) );
            //F4944       ENDIF
        }
        //F4945 !
        //F4946       IF(J.EQ.225)THEN
        if( J == 225 ) {
            //F4947         QTOT89 =QTOT(J)
            //UNUSED const float     QTOT89 = TANDSL->QTOT[ J ];
            //F4948         QSO289 =QSO2(J)
            //UNUSED const float     QSO289 = TANDSL->QSO2[ J ];
            //F4949         QDIR89 =QDIR(J)
            //UNUSED const float     QDIR89 = TANDSL->QDIR[ J ];
            //F4950         QOZ89  =QOZ(J)
            //UNUSED const float     QOZ89 = TANDSL->QOZ[ J ];
            //F4951         QBIO89 =QBIO(J)
            //UNUSED const float     QBIO89 = TANDSL->QBIO[ J ];
            //F4952         QLAND89=QLAND(J)
            //UNUSED const float     QLAND89 = TANDSL->QLAND[ J ];
            //F4953         QMN89  =QMN(J)
            //UNUSED const float     QMN89 = TANDSL->QMN[ J ];
            //F4954         QGH89  =QGH(J)
            //UNUSED const float     QGH89 = TANDSL->QGH[ J ];
            //F4955       ENDIF
        }
        //F4956 !
        //F4957       IF(J.EQ.226)THEN
        if( J == 226 ) {
            //F4958         QTOT90 =QTOT(J)
            //UNUSED const float QTOT90 = TANDSL->QTOT[ J ];
            //F4959         QSO290 =QSO2(J)
            //UNUSED const float QSO290 = TANDSL->QSO2[ J ];
            //F4960         QDIR90 =QDIR(J)
            //UNUSED const float QDIR90 = TANDSL->QDIR[ J ];
            //F4961         QOZ90  =QOZ(J)
            //UNUSED const float QOZ90 = TANDSL->QOZ[ J ];
            //F4962         QBIO90 =QBIO(J)
            //UNUSED const float QBIO90 = TANDSL->QBIO[ J ];
            //F4963         QLAND90=QLAND(J)
            //UNUSED const float QLAND90 = TANDSL->QLAND[ J ];
            //F4964         QMN90  =QMN(J)
            //UNUSED const float QMN90 = TANDSL->QMN[ J ];
            //F4965         QGH90  =QGH(J)
            //UNUSED const float QGH90 = TANDSL->QGH[ J ];
            //F4966         QTOTM =(QTOT89+QTOT90)/2.0
            //UNUSED const float QTOTM = ( QTOT89 + QTOT90 ) / 2.0;
            //F4967         QSO2M =(QSO289+QSO290)/2.0
            //UNUSED const float QSO2M = ( QSO289 + QSO290 ) / 2.0;
            //F4968         QDIRM =(QDIR89+QDIR90)/2.0
            //UNUSED const float QDIRM = ( QDIR89 + QDIR90 ) / 2.0;
            //F4969         QOZM  =(QOZ89+QOZ90)/2.0
            //UNUSED const float QOZM = ( QOZ89 + QOZ90 ) / 2.0;
            //F4970         QBIOM =(QBIO89+QBIO90)/2.0
            //UNUSED const float QBIOM = ( QBIO89 + QBIO90 ) / 2.0;
            //F4971         QLANDM=(QLAND89+QLAND90)/2.0
            //UNUSED const float QLANDM = ( QLAND89 + QLAND90 ) / 2.0;
            //F4972         QMNM  =(QMN89+QMN90)/2.0
            //UNUSED const float QMNM = ( QMN89 + QMN90 ) / 2.0;
            //F4973         QGHM  =(QGH89+QGH90)/2.0
            //UNUSED const float QGHM = ( QGH89 + QGH90 ) / 2.0;
            //F4974       ENDIF
        }
        //F4975 !
        //F4976   10  CONTINUE
    } // for
    //F4977 !
    //F4978       RETURN
    //F4979       END
    f_exit( __func__ );
} // deltaq
//F4980 !
//F4981 !  *******************************************
//F4982 !
//F4983       SUBROUTINE HALOCARB(N,C0,E,C1,Q,TAU00,TAUCH4)
void halocarb( const int N, float C0, float E, float* C1, float* Q, float TAU00, float TAUCH4 )
{
    f_enter( __func__ );
    //F4984       IMPLICIT REAL*4 (a-h,o-z), Integer (I-N)
    //F4985 !
    //F4986 !  N   GASNAME AND FORMULA
    //F4987 !  1   CF4           (CF4)
    //F4988 !  2   C2F6         (C2F6)
    //F4989 !  3   HFC125    (CHF2CF3)
    //F4990 !  4   HFC134a   (CH2FCF3)
    //F4991 !  5   HFC143a    (CH3CF3)
    //F4992 !  6   HFC227ea    (C3HF7)
    //F4993 !  7   HFC245ca   (C3H3F5)
    //F4994 !  8   SF6           (SF6)
    //F4995 !
    //F4996 !  THE FOLLOWING CENTRED DIFFERENCE FORMULA IS ONLY GOOD FOR
    //F4997 !   TAU.GT.(dt/2), SO IT FAILS IF TAU.LT.0.5yr.  AN ALTERNATIVE
    //F4998 !   IS TO USE THE EXACT SOLUTION OVER THE ONE YEAR INCREMENT.
    //F4999 !   THIS IS ...
    //F5000 !     XXX=E(J,K)/B(J)
    //F5001 !     EX=EXP(-1.0/TAU(J))
    //F5002 !     C(K) = TAU(J)*XXX*(1.0-EX) + C(K-1)*EX
    //F5003 !
    //F5004       DIMENSION B(10),TAU(10),ANFB(10)
    //F5005 !
    //F5006       B(1)=15.10
    //F5007       B(2)=23.68
    //F5008       B(3)=20.17
    //F5009       B(4)=17.14
    //F5010       B(5)=14.12
    //F5011       B(6)=28.57
    //F5012       B(7)=22.52
    //F5013       B(8)=25.05
    float B[ 8+1 ] = { 0, 15.10, 23.68, 20.17, 17.14, 14.12, 28.57, 22.52, 25.05 };
    //F5014 !
    //F5015       TAU(1)=50000.0
    //F5016       TAU(2)=10000.0
    //F5017       TAU(3)=   29.0
    //F5018       TAU(4)=   13.8
    //F5019       TAU(5)=   52.0
    //F5020       TAU(6)=   33.0
    //F5021       TAU(7)=    6.6
    //F5022       TAU(8)= 3200.0
    float TAU[ 8+1 ] = { 0, 50000.0, 10000.0, 29.0, 13.8, 52.0, 33.0, 6.6, 3200.0 };
    //F5023 !
    //F5024       ANFB(1)=0.08
    //F5025       ANFB(2)=0.26
    //F5026       ANFB(3)=0.23
    //F5027       ANFB(4)=0.15
    //F5028       ANFB(5)=0.13
    //F5029       ANFB(6)=0.30
    //F5030       ANFB(7)=0.23
    //F5031       ANFB(8)=0.52
    float ANFB[ 8+1 ] = { 0, 0.08, 0.26, 0.23, 0.15, 0.13, 0.30, 0.23, 0.52 };
    //F5032 !
    //F5033 !  LIFETIME CHANGE FACTOR (ADDED 000909)
    //F5034 !
    //F5035       IF((N.LE.2).OR.(N.EQ.8))THEN
    //F5036         FACTOR=1.0
    float FACTOR;
    if ( N <= 2 || N == 8 ) FACTOR = 1.0;
    //F5037       ELSE
    //F5038         FACTOR=TAUCH4/TAU00
    else FACTOR = TAUCH4 / TAU00;
    //F5039       ENDIF
    //F5040       TAU(N)=TAU(N)*FACTOR
    TAU[ N ] = TAU[ N ] * FACTOR;
    //F5041 !
    //F5042       XXX=E/B(N)
    float XXX = E / B[ N ];
    //F5043       TT=1.0/(2.0*TAU(N))
    float TT = 1.0 / ( 2.0 * TAU[ N ] );
    //F5044 !
    //F5045       IF(TAU(N).GE.1.0)THEN
    if( TAU[ N ] >= 1.0 ) {
        //F5046         UU=1.-TT
        float UU = 1.0 - TT;
        //F5047         C1=(C0*UU+XXX)/(1.+TT)
        *C1 = ( C0 * UU + XXX ) / ( 1.0 + TT );
        //F5048       ELSE
    } else {
        //F5049         EX=EXP(-1.0/TAU(N))
        float EX = exp( float( -1.0 / TAU[ N ] ) );
        //F5050         C1=TAU(N)*XXX*(1.0-EX)+C0*EX
        *C1 = TAU[ N ] * XXX * ( 1.0 - EX ) + C0 * EX;
        //F5051       ENDIF
    }
    //F5052 !
    //F5053       Q=C1*ANFB(N)/1000.
    *Q = *C1 * ANFB[ N ] / 1000.0;
    //F5054 !
    //F5055       RETURN
    //F5056       END
    f_exit( __func__ );
} // halocarb
//F5057 !
//F5058 !  *******************************************
//F5059 !
//F5060       SUBROUTINE INITCAR(NN,D80,F80)
void initcar( const int NN, const float D80, const float F80, COBS_block* COBS, 
             CARB_block* CARB, CAR_block* CAR )
{
    f_enter( __func__ );
    //F5061       IMPLICIT REAL*4 (a-h,o-z), Integer (I-N)
    //F5062 !
    //F5063       parameter (iTp=740)
    //F5064 !
    //F5065       INTEGER FERTTYPE,TOTEM,CONVTERP
    //F5066 !
    //F5067       COMMON/COBS/COBS(0:236)
    //F5068 !
    //F5069       COMMON/CARB/CCO2(4,224:iTp),EDGROSS(4,226:iTp),EF(226:iTp+1), &
    //F5070       REGROW(4,226:iTp),PL(4,226:iTp),HL(4,226:iTp),SOIL(4,226:iTp), &
    //F5071       TTT(226:iTp),ESUM(226:iTp),ETOT(4,226:iTp),EDNET90(4), &
    //F5072       FOC(4,226:iTp),CO2(0:iTp),CO2SAVE(0:iTp)
    //F5073 !
    //F5074       COMMON/CAR/EL1,EL2,EL3,TINV0(5),TINV(4,5),A(3,5),AA(4,5), &
    //F5075       BCO2(4),BTGPP,BTRESP,BTHUM,GAMP,GPP0,RESP0,QA0,U0,C0,B340(4), &
    //F5076       PHI,RG,TAUP,TAUH,TAUS,THP,THS,THH0,THS0,THPL,G1,G2,G3,FACTOR, &
    //F5077       EL21,EL32,XX1,XX2,XX3,XX4,XX5,XX6,DEE1,DEE2,DEE3,DEE4,DEE5,DEE6, &
    //F5078       FL1,FL2,FL3,XL,GAMH,GAMS,QS0,BTSOIL,FERTTYPE,TOTEM,CONVTERP, &
    //F5079       R(4),CPART(4,5),DELMASS(4,226:iTp),ABFRAC(4,226:iTp)
    //F5080 !
    //F5081 !  FIRST INITIALISE PARAMETERS THAT DEPEND ON MM
    //F5082 !
    //F5083 !  BCO2 ***************************
    //F5084 !
    //F5085 !  COMPLETELY REVISED ON 6/12/95
    //F5086 !  REVISED AGAIN ON 8/8/95. CUBIC RETAINED AS BASIC APPROXIMATION
    //F5087 !   FORMULA FOR BCO2, BUT CORRECTION ADDED FOR D80<0.4 OR >1.8.
    //F5088 !   FORMULA RAPIDLY LOSES ACCURACY OUTSIDE RANGE SHOWN BELOW.
    //F5089 !  FORMULA ADDED TO ACCOUNT FOR CASES WHERE F80S.NE.2.0 USES
    //F5090 !   FACT THAT BCO2 FOR (F80),(D80) IS APPROX THE SAME AS FOR
    //F5091 !   (F80-DEL),(D80-DEL)
    //F5092 !
    //F5093       DD =(D80-1.1)-(F80-2.0)
    const float DD = ( D80 - 1.1 ) - ( F80 - 2.0 );
    //F5094       DD2=DD*DD
    const float DD2 = DD * DD;
    //F5095       DD3=DD*DD2
    const float DD3 = DD * DD2;
    //F5096       BCO2(NN)=0.27628+DD*0.301304+DD2*0.060673+DD3*0.012383
    CAR->BCO2[ NN ] = 0.27628 + DD*0.301304 + DD2*0.060673 + DD3*0.012383;
    //F5097       XD=ABS(DD)-0.7
    const float XD = fabs( DD ) - 0.7;
    //F5098       IF(XD.GT.0.0) BCO2(NN)=BCO2(NN)-(0.0029*XD-.022*XD*XD)
    if( XD > 0.0 ) CAR->BCO2[ NN ] -= ( 0.0029*XD - 0.022 * XD * XD );
    //F5099 !
    //F5100 !  CORRECTION FOR F80.NE.2.0.  ERROR IN BCO2 LESS THAN 0.0014
    //F5101 !   FOR 1.0<F80<3.0
    //F5102 !
    //F5103       BCO2(NN)=BCO2(NN)-0.0379*(F80-2.0)
    CAR->BCO2[ NN ] -= 0.0379 * ( F80-2.0 );
    //F5104 !
    //F5105 !  COMPARISON OF FIT TO TRUE BCO2 FOR F80=2.0
    //F5106 !
    //F5107 !     D80 BCO2(TRUE)  BCO2(FIT)  TRUE-FIT
    //F5108 !     0.0   0.00365    0.00414   -0.00049
    //F5109 !     0.1   0.02438    0.02438    0.00000
    //F5110 !     0.4   0.09085    0.09085    0.00000
    //F5111 !     0.8   0.19092    0.19102   -0.00010
    //F5112 !     1.1   0.27628    0.27628    0.00000
    //F5113 !     1.4   0.37237    0.37247   -0.00010
    //F5114 !     1.8   0.52117    0.52117    0.00000
    //F5115 !     2.2   0.69994    0.69997   -0.00003
    //F5116 !     2.6   0.91868    0.91830    0.00038
    //F5117 !     3.0   1.19350    1.18092    0.01258
    //F5118 !
    //F5119 !  EDNET ETC. **********************
    //F5120 !
    //F5121 !  LATEST (8/8/95) VALUES. NOT YET INSERTED.
    //F5122 !
    //F5123 !  YEAR     EDNET    REGROW     PLANT     HLITT      SOIL      ETOT
    //F5124 !  1989      .724     1.251   707.895    84.195  1416.528   284.230
    //F5125 !  1990      .762     1.254   708.285    84.300  1416.655   289.707
    //F5126 !
    //F5127 !  EDNET90 CORRECTED ON 8/8/95. LEADING TERMS ONLY IN  OTHER ITEMS
    //F5128 !   CORRECTED.
    //F5129 !
    //F5130       EDNET90(NN)    =0.762+1.1185*DD-0.0020*DD2
    CARB->EDNET90[ NN ] = 0.762 + 1.1185 * DD - 0.0020 * DD2;
    //F5131 !
    //F5132       REGROW(NN,226) =1.254+0.584*DD-0.0173*DD2
    CARB->REGROW.setval( 1.254 + 0.584 * DD - 0.0173 * DD2, NN, 226 );
    //F5133       EDGROSS(NN,226)=EDNET90(NN)+REGROW(NN,226)
    CARB->EDGROSS.setval( CARB->EDNET90[ NN ] + CARB->REGROW.getval( NN, 226 ), NN, 226 );
    //F5134 !
    //F5135       PL(NN,226)     = 708.285-6.131*DD+0.137*DD2
    CARB->PL.setval(  708.285 - 6.131 * DD + 0.137 * DD2, NN, 226 );
    //F5136       HL(NN,226)     =  84.300+4.569*DD-0.039*DD2
    CARB->HL.setval( 84.300 + 4.569 * DD - 0.039 * DD2, NN, 226 );
    //F5137       SOIL(NN,226)   =1416.655+1.564*DD-0.099*DD2
    CARB->SOIL.setval( 1416.655 + 1.564 * DD - 0.099 * DD2, NN, 226 );
    //F5138 !
    //F5139 !  NEXT THREE ITEMS NOT UPDATED
    //F5140 !
    //F5141       FOC(NN,226)    =   2.24
    CARB->FOC.setval( 2.24, NN, 226 );
    //F5142       DELMASS(NN,226)=   3.609
    CAR->DELMASS.setval( 3.609, NN, 226 );
    //F5143       ABFRAC(NN,226) =   0.514-0.084*(D80-1.0)+0.013*(D80-1.0)**2
    CAR->ABFRAC.setval( 0.514 - 0.084 * ( D80-1.0 ) + 0.013 * pow( D80-1.0, 2 ), NN, 226 );
    //F5144 !
    //F5145 !  ******************************
    //F5146 !
    //F5147 !  SPECIFY INIT, END-1988, END-1989 AND END-1990 CO2 CONCS.
    //F5148 !
    //F5149       C0          =COBS(0)
    CAR->C0 = COBS->COBS[ 0 ];
    //F5150       CCO2(NN,224)=COBS(224)
    CARB->CCO2.setval( COBS->COBS[ 224 ], NN, 224);
    //F5151       CCO2(NN,225)=COBS(225)
    CARB->CCO2.setval( COBS->COBS[ 225 ], NN, 225 );
    //F5152       CCO2(NN,226)=COBS(226)
    CARB->CCO2.setval( COBS->COBS[ 226 ], NN, 226 );
    //F5153 !
    //F5154       ETOT(NN,226) = 290.797
    CARB->ETOT.setval( 290.797, NN, 226 );
    //F5155 !
    //F5156       PL0  = 750.0
    const float PL0 = 750.0;
    //F5157       HL0  =  80.0
    const float HL0 = 80.0;
    //F5158       SOIL0=1450.0
    const float SOIL0 = 1450.0;
    //F5159 !
    //F5160       FERTTYPE=2
    CAR->FERTTYPE = 2;
    //F5161       TOTEM   =1
    CAR->TOTEM = 1;
    //F5162       CONVTERP=1
    CAR->CONVTERP = 1;
    //F5163       FACTOR  =2.123
    CAR->FACTOR = 2.123;
    //F5164 !
    //F5165 !  ADJUST INVERSE DECAY TIMES ACCORDING TO 1980S-MEAN OCEAN FLUX.
    //F5166 !   PSI CONTROLS THE MAGNITUDE OF THE FLUX INTO THE OCEAN.
    //F5167 !   (IF PSI IS LESS THAN 1, FLUX IS LESS THAN IN THE ORIGINAL
    //F5168 !   MAIER-REIMER AND HASSELMANN MODEL.)
    //F5169 !  A SIMPLE BUT ACCURATE APPROXIMATE EMPIRICAL EXPRESSION IS USED
    //F5170 !   TO ESTIMATE PSI AS A FUNCTION OF F80SIN. THE PSI-F80SIN RELATION
    //F5171 !   DEPENDS ON THE ASSUMED HISTORY OF OBSERVED CONCENTRATION
    //F5172 !   CHANGES. THUS, DIFFERENT PSI-F80SIN RELATIONSHIPS APPLY TO
    //F5173 !   DIFFERENT CO2 HISTORIES.
    //F5174 !
    //F5175 !  REVISED ON 3/16/95, BUT APSI ONLY (HENCE OK FOR F80=2.0 ONLY)
    //F5176 !  REVISED AGAIN ON 4/2/95
    //F5177 !  REVISED AGAIN ON 6/12/95
    //F5178 !  REVISED AGAIN ON 8/8/95
    //F5179 !
    //F5180       FF  =F80-2.0
    const float FF = F80 - 2.0;
    //F5181       FF2 =FF*FF
    const float FF2 = FF * FF;
    //F5182       APSI=1.029606
    const float APSI = 1.029606;
    //F5183       BPSI=0.873692
    const float BPSI = 0.873692;
    //F5184       CPSI=0.165084
    const float CPSI = 0.165084;
    //F5185       PSI =APSI+BPSI*FF+CPSI*FF2
    const float PSI = APSI + BPSI * FF + CPSI * FF2;
    //F5186 !
    //F5187 !  PSI COMPARISON
    //F5188 !
    //F5189 !      F80   PSITRUE    PSIEST
    //F5190 !      1.0   .320730   .320998
    //F5191 !      1.5   .634031   .634031
    //F5192 !      2.0  1.029606  1.029606
    //F5193 !      2.5  1.507723  1.502733
    //F5194 !      3.0  2.074348  2.068382
    //F5195 !
    //F5196       DO J=1,5
    for( int J=1; J<=5; J++ ) {
        //F5197       TINV(NN,J)=TINV0(J)*PSI
        CAR->TINV.setval( CAR->TINV0.getval( J ) * PSI, NN, J );
        //F5198       END DO
    }
    //F5199 !
    //F5200 !  SPECIFY 1990 (J=226) PARTIAL CONCS AND CONVOLUTION CONSTANTS
    //F5201 !   REVISED ON 3/16/95
    //F5202 !   REVISED AGAIN ON 4/2/95
    //F5203 !   REVISED AGAIN ON 6/12/95 (AA NOT CHANGED, CPART VERY MINOR)
    //F5204 !   REVISED AGAIN ON 8/8/95 (AA AND CPART CHANGED)
    //F5205 !
    //F5206       CPART(NN,1)=18.40290
    CAR->CPART[ NN ][ 1 ] = 18.40290;
    //F5207       CPART(NN,2)=25.51023
    CAR->CPART[ NN ][ 2 ] = 25.51023;
    //F5208       CPART(NN,3)=22.51496
    CAR->CPART[ NN ][ 3 ] = 22.51496;
    //F5209       CPART(NN,4)= 9.64078
    CAR->CPART[ NN ][ 4 ] = 9.64078;
    //F5210       CPART(NN,5)= 0.38709
    CAR->CPART[ NN ][ 5 ] = 0.38709;
    //F5211 !
    //F5212       AA(NN,1)= 0.13486
    CAR->AA.setval( 0.13486, NN, 1 );
    //F5213       AA(NN,2)= 0.22091
    CAR->AA.setval( 0.22091, NN, 2 );
    //F5214       AA(NN,3)= 0.28695
    CAR->AA.setval( 0.28695, NN, 3 );
    //F5215       AA(NN,4)= 0.26033
    CAR->AA.setval( 0.26033, NN, 4 );
    //F5216       AA(NN,5)= 0.09695
    CAR->AA.setval( 0.09695, NN, 5 );
    //F5217 !
    //F5218 !  SPECIFY OR CALCULATE OTHER MAIN MODEL PARAMETERS
    //F5219 !
    //F5220 !  GPP IS THE PART OF GPP THAT IS NOT IMMEDIATELY RESPIRED BY
    //F5221 !   LEAVES AND GROUND VEG
    //F5222 !  RESP IS THE PART OF RESPIRATION THAT COMES FROM THE TREES
    //F5223 !  RG=RESP0/GPP0
    //F5224 !
    //F5225       GPP0 =76.0
    CAR->GPP0   = 76.0;
    //F5226       RESP0=14.0
    CAR->RESP0  = 14.0;
    //F5227       PHI  =0.98
    CAR->PHI    = 0.98;
    //F5228       XL   =0.05
    CAR->XL     = 0.05;
    //F5229       G1   =0.35
    CAR->G1     = 0.35;
    //F5230       G2   =0.60
    CAR->G2     = 0.60;
    //F5231       GAMP =0.70
    CAR->GAMP   = 0.70;
    //F5232       GAMH =0.05
    CAR->GAMH   = 0.05;
    //F5233 !
    //F5234       RG  =RESP0/GPP0
    CAR->RG     = CAR->RESP0 / CAR->GPP0;
    //F5235       G3  =1.0-G1-G2
    CAR->G3     = 1.0 - CAR->G1 - CAR->G2;
    //F5236       GAMS=1.0-GAMP-GAMH
    CAR->GAMS   = 1.0 - CAR->GAMP - CAR->GAMH;
    //F5237 !
    //F5238 !  TEMPERATURE FEEDBACK TERMS SPECIFIED IN MAG3GAS.CFG
    //F5239 !
    //F5240       THPL=G1*GPP0-RESP0
    CAR->THPL   = CAR->G1 * CAR->GPP0 - CAR->RESP0;
    //F5241       TAUP=PL0/THPL
    CAR->TAUP   = PL0 / CAR->THPL;
    //F5242       THP =1./(2.*TAUP)
    CAR->THP    = 1.0 / ( 2.0 * CAR->TAUP );
    //F5243       TAUH=HL0/(G2*GPP0+PHI*THPL)
    CAR->TAUH   = HL0 / ( CAR->G2 * CAR->GPP0 + CAR->PHI * CAR->THPL );
    //F5244       THH0=0.5/TAUH
    CAR->THH0   = 0.5 / CAR->TAUH;
    //F5245       TAUS=SOIL0/(GPP0-RESP0-(1.0-XL)*HL0/TAUH)
    CAR->TAUS   = SOIL0 / ( CAR->GPP0 - CAR->RESP0 - ( 1.0-CAR->XL ) * HL0/CAR->TAUH );
    //F5246       THS0=0.5/TAUS
    CAR->THS0   = 0.5 / CAR->TAUS;
    //F5247 !
    //F5248 !  QA0 IS THE INITIAL HLITT DECOMPOSITION FLUX TO ATMOSPHERE.
    //F5249 !
    //F5250       QA0=(1.0-XL)*HL0/TAUH
    CAR->QA0    = ( 1.0 - CAR->XL ) * HL0 / CAR->TAUH;
    //F5251 !
    //F5252 !  QS0 IS THE INITIAL HLITT FLUX TO SOIL : U0 DITTO SOIL TO ATMOS
    //F5253 !
    //F5254       QS0=XL*HL0/TAUH
    CAR->QS0    = CAR->XL * HL0 / CAR->TAUH;
    //F5255       U0 =SOIL0/TAUS
    CAR->U0     = SOIL0 / CAR->TAUS;
    //F5256 !
    //F5257 !      FOC(NN,0) =0.0
    //F5258 !      EF(0)     =0.0
    //F5259 !      ESUM(0)   =0.0
    //F5260 !      ETOT(NN,0)=0.0
    //F5261 !
    //F5262       FL1 =EL1/100.
    CAR->FL1    = CAR->EL1 / 100.0;
    //F5263       FL2 =EL2/100.
    CAR->FL2    = CAR->EL2 / 100.0;
    //F5264       FL3 =EL3/100.
    CAR->FL3    = CAR->EL3 / 100.0;
    //F5265       EL21=FL2-FL1
    CAR->EL21   = CAR->FL2 - CAR->FL1;
    //F5266       EL32=FL3-FL2
    CAR->EL32   = CAR->FL3 - CAR->FL2;
    //F5267       XX1 =FL1-DEE1
    CAR->XX1    = CAR->FL1 - CAR->DEE1;
    //F5268       XX2 =FL1+DEE2
    CAR->XX2    = CAR->FL1 + CAR->DEE2;
    //F5269       XX3 =FL2-DEE3
    CAR->XX3    = CAR->FL2 - CAR->DEE3;
    //F5270       XX4 =FL2+DEE4
    CAR->XX4    = CAR->FL2 + CAR->DEE4;
    //F5271       XX5 =FL3-DEE5
    CAR->XX5    = CAR->FL3 - CAR->DEE5;
    //F5272       XX6 =FL3+DEE6
    CAR->XX6    = CAR->FL3 + CAR->DEE6;
    //F5273 !
    //F5274       RETURN
    //F5275       END
    f_exit( __func__ );
} // initcar
//F5276 !
//F5277 !  *******************************************
//F5278 !
//F5279       SUBROUTINE CARBON(MM,TEM,EFOSS,ENETDEF,CPP,CPREV,C, &
//F5280       PL ,HU ,SO ,REGRO ,ETOT , &
//F5281       PL1,HU1,SO1,REGRO1,ETOT1, &
//F5282       SUMEM1,FLUX,DELM,EGROSSD,C1)
void carbon( const int MM, float TEM, float EFOSS, float ENETDEF, float CPP, float CPREV, float C,
            float PL, float HU, float SO, float REGRO, float ETOT,
            float* PL1, float* HU1, float* SO1, float* REGRO1, float* ETOT1,
            float* SUMEM1, float* FLUX, float* DELM, float* EGROSSD, float* C1,
            CAR_block* CAR )
{
    f_enter( __func__ );
    //F5283       IMPLICIT REAL*4 (a-h,o-z), Integer (I-N)
    //F5284 !
    //F5285       parameter (iTp=740)
    //F5286 !
    //F5287       DIMENSION A1(4,5)
    float A1[ 4+1 ][ 5+1 ];
    //F5288 !
    //F5289       INTEGER FERTTYPE,TOTEM,CONVTERP
    //F5290 !
    //F5291       COMMON/CAR/EL1,EL2,EL3,TINV0(5),TINV(4,5),A(3,5),AA(4,5), &
    //F5292       BCO2(4),BTGPP,BTRESP,BTHUM,GAMP,GPP0,RESP0,QA0,U0,C0,B340(4), &
    //F5293       PHI,RG,TAUP,TAUH,TAUS,THP,THS,THH0,THS0,THPL,G1,G2,G3,FACTOR, &
    //F5294       EL21,EL32,XX1,XX2,XX3,XX4,XX5,XX6,DEE1,DEE2,DEE3,DEE4,DEE5,DEE6, &
    //F5295       FL1,FL2,FL3,XL,GAMH,GAMS,QS0,BTSOIL,FERTTYPE,TOTEM,CONVTERP, &
    //F5296       R(4),CPART(4,5),DELMASS(4,226:iTp),ABFRAC(4,226:iTp)
    //F5297 !
    //F5298 !  TERRESTRIAL CARBON CYCLE SECTION
    //F5299 !
    //F5300 !  FOR CALCULATING FEEDBACKS, NEED TO USE EXTRAPOLATED CO2
    //F5301 !   CONCENTRATION FOR THE MIDPOINT OF YEAR 'J' AFTER 1990.
    //F5302 !   FORMULA BELOW GIVES QUADRATIC EXTRAPOLATED RESULT.
    //F5303 !
    //F5304 !      CBAR=C+0.5*(C-CPREV)
    //F5305       CBAR=(3.0*CPP-10.0*CPREV+15.0*C)/8.0
    float CBAR = (3.0 * CPP - 10.0 * CPREV + 15.0 * C ) / 8.0;
    //F5306 !
    //F5307 !  NOW CALCULATE GROSS DEFOR CORRESP TO INPUT NET DEFOR
    //F5308 !
    //F5309       REGRO1=(REGRO*(TAUP-0.5)+GAMP*ENETDEF)/(TAUP+0.5-GAMP)
    *REGRO1 = ( REGRO * ( CAR->TAUP-0.5 ) + CAR->GAMP * ENETDEF ) / ( CAR->TAUP+0.5-CAR->GAMP );
    //F5310       EGROSSD=ENETDEF+REGRO1
    *EGROSSD = ENETDEF + *REGRO1;
    //F5311 !
    //F5312 !  TEMPERATURE TERMS
    //F5313 !
    //F5314       FG=EXP(BTGPP*TEM)
    const float FG = exp( float( CAR->BTGPP * TEM ) );
    //F5315       FR=EXP(BTRESP*TEM)
    const float FR = exp( float( CAR->BTRESP * TEM ) );
    //F5316       FH=EXP(BTHUM*TEM)
    const float FH = exp( CAR->BTHUM * TEM ); 
    //F5317       FS=EXP(BTSOIL*TEM)
    const float FS = exp( float( CAR->BTSOIL * TEM ) );
    //F5318 !
    //F5319 !  DEFORESTATION TERMS
    //F5320 !
    //F5321       GDP=GAMP*EGROSSD
    const float GDP = CAR->GAMP * *EGROSSD;
    //F5322       GDH=GAMH*EGROSSD
    const float GDH = CAR->GAMH * *EGROSSD;
    //F5323       GDS=GAMS*EGROSSD
    const float GDS = CAR->GAMS * *EGROSSD;
    //F5324 !
    //F5325 !  CO2 FERTILIZATION TERMS. LOG FORM FIRST THEN FORM USED BY ENTING,
    //F5326 !   THEN FORM USED BY GIFFORD. ALSO CALCULATE BCO2 AT C=CCC=340PPMV
    //F5327 !
    //F5328       IF(FERTTYPE.EQ.1)THEN
    float Y1;
    if( CAR->FERTTYPE == 1 ) {
        //F5329         Y1=BCO2(MM)*ALOG(CBAR/C0)+1.0
        Y1 = CAR->BCO2[ MM ] * log( float( CBAR/CAR->C0 ) ) + 1.0;
        //F5330         B340(MM)=BCO2(MM)
        CAR->B340[ MM ] = CAR->BCO2[ MM ];
        //F5331       ELSE
    } else {
        //F5332 !
        //F5333 !        CONCBASE=80.
        //F5334 !        GINF=2.4
        //F5335 !        BEE=C0*(GINF-1.0)-GINF*CONCBASE
        //F5336 !        GEE=GINF*(CBAR-CONCBASE)/(CBAR+BEE)
        //F5337 !        Y1=BCO2(MM)*(GEE-1.0)+1.0
        //F5338 !
        //F5339         CB=31.
        const float CB = 31.0;
        //F5340         R(MM)=(1.+BCO2(MM)*ALOG(680./C0))/(1.+BCO2(MM)*ALOG(340./C0))
        CAR->R[ MM ] = ( 1.0 + CAR->BCO2[ MM ] * log( float( 680.0/CAR->C0 ) ) ) / ( 1.0 + CAR->BCO2[ MM ] * log( float( 340.0/CAR->C0 ) ) );
        //F5341         AR=680.-CB
        const float AR = 680.0 - CB;
        //F5342         BR=340.-CB
        const float BR = 340.0 - CB;
        //F5343         BEE=999.9
        float BEE = 999.9;
        //F5344         IF(R(MM).NE.1.)BEE=(AR/BR-R(MM))/(R(MM)-1.)/AR
        if( CAR->R[ MM ] != 1.0 ) BEE = ( AR/BR - CAR->R[ MM ]) / ( CAR->R[ MM ] - 1.0 ) / AR;
        //F5345         DR=CBAR-CB
        const float DR = CBAR - CB;
        //F5346         CR=C0-CB
        const float CR = CAR->C0 - CB;
        //F5347         Y1=1.0
        Y1 = 1.0;    // weird
        //F5348         CCC=340.
        const float CCC = 340.0;
        //F5349         B340(MM)=0.0
        CAR->B340[ MM ] = 0.0;
        //F5350         IF(R(MM).NE.1.)THEN
        if( CAR->R[ MM ] != 1.0 ) {
            //F5351           Y1=(1./CR+BEE)/(1./DR+BEE)
            Y1 = ( 1.0/CR + BEE ) / ( 1.0/DR + BEE );
            //F5352           B340(MM)=(1./CR+BEE)*CCC/(1.+BEE*(CCC-CB))**2
            CAR->B340[ MM ] = ( 1.0/CR + BEE ) * CCC / pow( 1.0+BEE*( CCC-CB ), 2 );
            //F5353         ENDIF
        }
        //F5354       ENDIF
    }
    //F5355 !
    //F5356       GPP=GPP0*Y1*FG
    const float GPP = CAR->GPP0 * Y1 * FG;
    //F5357 !
    //F5358       PGPP=GPP*G1
    const float PGPP = GPP * CAR->G1;
    //F5359       DPGPP=PGPP-GPP0*G1
    //UNUSED const float DPGPP = PGPP - CAR->GPP0 * CAR->G1; 
    //F5360       HGPP=GPP*G2
    const float HGPP = GPP * CAR->G2;
    //F5361       DHGPP=HGPP-GPP0*G2
    //UNUSED const float DHGPP = HGPP - CAR->GPP0 * CAR->G2; 
    //F5362       SGPP=GPP*G3
    const float SGPP = GPP * CAR->G3;
    //F5363       DSGPP=SGPP-GPP0*G3
    //UNUSED const float DSGPP = SGPP - CAR->GPP0 * CAR->G3; 
    //F5364       RESP=RESP0*Y1*FR
    const float RESP = CAR->RESP0 * Y1 * FR;
    //F5365 !
    //F5366       DRESP=RESP-RESP0
    //UNUSED const float DRESP = RESP - CAR->RESP0; 
    //F5367 !
    //F5368 !  NEW PLANT MASS
    //F5369 !
    //F5370       PTERM=PGPP-RESP-GDP
    const float PTERM = PGPP - RESP - GDP;
    //F5371       PL1=(PL*(1.0-THP)+PTERM)/(1.0+THP)
    *PL1 = ( PL * ( 1.0-CAR->THP ) + PTERM ) / ( 1.0 + CAR->THP );
    //F5372 !
    //F5373 !  NEW HLITT MASS
    //F5374 !
    //F5375       THH=THH0*FH
    const float THH = CAR->THH0 * FH;
    //F5376       Y2=THP*(PL+PL1)
    const float Y2 = CAR->THP * ( PL + *PL1 );
    //F5377       HTERM=HGPP+PHI*Y2-GDH
    const float HTERM = HGPP + CAR->PHI * Y2 - GDH;
    //F5378       HU1=(HU*(1.0-THH)+HTERM)/(1.0+THH)
    *HU1 = ( HU * ( 1.0 - THH ) + HTERM ) / ( 1.0 + THH );
    //F5379 !
    //F5380 !  NEW SOIL MASS
    //F5381 !
    //F5382       THS=THS0*FS
    CAR->THS = CAR->THS0 * FS;
    //F5383       Y3=THH*(HU+HU1)
    const float Y3 = THH * ( HU + *HU1 );
    //F5384       STERM=SGPP+(1.0-PHI)*Y2+XL*Y3-GDS
    const float STERM = SGPP + ( 1.0 - CAR->PHI ) * Y2 + CAR->XL * Y3 - GDS;
    //F5385       SO1=(SO*(1.0-THS)+STERM)/(1.0+THS)
    *SO1 = ( SO * ( 1.0 - CAR->THS ) + STERM ) / ( 1.0 + CAR->THS );
    //F5386 !
    //F5387 !  FERTILIZATION FEEDBACK FLUX
    //F5388 !
    //F5389       BFEED1=GPP0-RESP0-(GPP-RESP)
    //UNUSED const float BFEED1 = CAR->GPP0 - CAR->RESP0 - ( GPP - RESP );
    //F5390 !
    //F5391 !  FEEDBACK FLUX DUE TO ACTIVE TERRESTRIAL BIOMASS
    //F5392 !
    //F5393       HUBAR=(HU+HU1)/2.0
    //UNUSED const float HUBAR = ( HU + *HU1 ) / 2.0;
    //F5394       SOBAR=(SO+SO1)/2.0
    //UNUSED const float SOBAR = ( SO + *SO1 ) / 2.0;
    //F5395       QA1=FH*(1.0-XL)*HUBAR/TAUH
    //UNUSED const float QA1 = FH * ( 1.0 - CAR->XL ) * HUBAR / CAR->TAUH;
    //F5396       U1=FS*SOBAR/TAUS
    //UNUSED const float U1 = FS * SOBAR / CAR->TAUS;
    //F5397       AFEED1=QA1+U1-(QA0+U0)
    //UNUSED const float AFEED1 = QA1 + U1 - ( CAR->QA0 + CAR->U0 );
    //F5398 !
    //F5399 !  TOTAL BIOMASS CHANGE
    //F5400 !
    //F5401       DELSOIL=HU1+SO1-HU-SO
    const float DELSOIL = *HU1 + *SO1 - HU - SO;
    //F5402       DELB=PL1-PL+DELSOIL
    const float DELB = *PL1 - PL + DELSOIL;
    //F5403 !
    //F5404 !  ANN SUM OF EMISSIONS, OCEAN FLUX AND CUM SUM OF EMISSIONS.
    //F5405 !
    //F5406       SUMEM1=EFOSS-DELB
    *SUMEM1 = EFOSS - DELB;
    //F5407       FLUX=SUMEM1-FACTOR*DELC
    static float DELC = 0.0;    // this is exceedingly weird -- a static var -- see note in documentation
    *FLUX = *SUMEM1 - CAR->FACTOR * DELC;
    //F5408       IF(TOTEM.EQ.1)ETOT1=ETOT+SUMEM1
    if( CAR->TOTEM == 1 ) *ETOT1 = ETOT + *SUMEM1;
    //F5409       IF(TOTEM.NE.1)ETOT1=ETOT+EFOSS+EGROSSD
    if( CAR->TOTEM != 1 ) *ETOT1 = ETOT + EFOSS + *EGROSSD;
    //F5410 !      IF(TOTEM.NE.1)ETOT1=ETOT+EFOSS-DELB
    //F5411 !
    //F5412 !  CALCULATE CONVOLUTION CONSTANTS AT THE END OF YEAR I. THIS
    //F5413 !    MAY BE DONE IN TWO WAYS, DETERMINED BY CONVTERP. NORMALLY
    //F5414 !    CONVTERP=1 SHOULD BE USED.
    //F5415 !
    //F5416       EE=ETOT1/100.
    const float EE = *ETOT1 / 100.0;
    //F5417 !
    //F5418       IF(CONVTERP.EQ.1)THEN
    if( CAR->CONVTERP == 1 ) {
        //F5419 !
        //F5420         DO K=1,5
        for( int K=1; K<=5; K++ ) {
            //F5421         D21=(A(2,K)-A(1,K))/EL21
            const float D21 = ( CAR->A.getval( 2, K ) - CAR->A.getval( 1, K ) ) / CAR->EL21;
            //F5422         D32=(A(3,K)-A(2,K))/EL32
            const float D32 = ( CAR->A.getval( 3, K ) - CAR->A.getval( 2, K ) ) / CAR->EL32;
            //F5423         IF(EE.LE.XX1)THEN
            if( EE <= CAR->XX1 ) {
                //F5424           A1(MM,K)=A(1,K)
                A1[ MM ][ K ] = CAR->A.getval( 1, K );
                //F5425 !
                //F5426         ELSE IF(EE.LE.XX2)THEN
            } else if( EE <= CAR->XX2 ) {
                //F5427           DY=D21*DEE2
                const float DY = D21 * CAR->DEE2;
                //F5428           X1=DEE1+DEE2
                const float X1 = CAR->DEE1 + CAR->DEE2;
                //F5429           X12=X1*X1
                const float X12 = X1 * X1;
                //F5430           AX=(D21-2.0*DY/X1)/X12
                const float AX = ( D21 - 2.0 * DY / X1 ) / X12;
                //F5431           BX=(3.*DY/X1-D21)/X1
                const float BX = ( 3.0 * DY/X1 - D21 ) / X1;
                //F5432           U=EE-XX1
                const float U = EE - CAR->XX1;
                //F5433           U2=U*U
                const float U2 = U * U;
                //F5434           A1(MM,K)=A(1,K)+BX*U2+AX*U2*U
                A1[ MM ][ K ] = CAR->A.getval( 1, K ) + BX * U2 + AX * U2 * U;
                //F5435 !
                //F5436         ELSE IF(EE.LE.XX3)THEN
            } else if( EE <= CAR->XX3 ) {
                //F5437           U=EE-FL1
                const float U = EE - CAR->FL1;
                //F5438           A1(MM,K)=A(1,K)+D21*U
                A1[ MM ][ K ] = CAR->A.getval( 1, K ) + D21 * U;
                //F5439 !
                //F5440         ELSE IF(EE.LE.XX4)THEN
            } else if( EE <= CAR->XX4 ) {
                //F5441           Y0=A(2,K)-D21*DEE3
                const float Y0 = CAR->A.getval( 2, K ) - D21 * CAR->DEE3;
                //F5442           DY=D21*DEE3+D32*DEE4
                const float DY = D21 * CAR->DEE3 + D32 * CAR->DEE4;
                //F5443           DD=D21+D32
                const float DD = D21 + D32;
                //F5444           X1=DEE3+DEE4
                const float X1 = CAR->DEE3 + CAR->DEE4;
                //F5445           X12=X1*X1
                const float X12 = X1 * X1;
                //F5446           AX=(DD-2.0*DY/X1)/X12
                const float AX = ( DD-2.0 * DY / X1 ) / X12;
                //F5447           BX=(3.*DY/X1-DD-D21)/X1
                const float BX = ( 3.0 * DY/X1 - DD - D21 ) / X1;
                //F5448           U=EE-XX3
                const float U = EE - CAR->XX3;
                //F5449           U2=U*U
                const float U2 = U * U;
                //F5450           A1(MM,K)=Y0+D21*U+BX*U2+AX*U2*U
                A1[ MM ][ K ] = Y0 + D21 * U + BX * U2 + AX * U2 * U;
                //F5451 !
                //F5452         ELSE IF(EE.LE.XX5)THEN
            } else if ( EE <= CAR->XX5 ) {
                //F5453           U=EE-FL2
                const float U = EE - CAR->FL2;
                //F5454           A1(MM,K)=A(2,K)+D32*U
                A1[ MM ][ K ] = CAR->A.getval( 2, K ) + D32 * U;
                //F5455 !
                //F5456         ELSE IF(EE.LE.XX6)THEN
            } else if( EE <= CAR->XX6 ) {
                //F5457           Y0=A(3,K)-D32*DEE5
                const float Y0 = CAR->A.getval( 3, K ) - D32 * CAR->DEE5;
                //F5458           DY=D32*DEE5
                const float DY = D32 * CAR->DEE5;
                //F5459           X1=DEE5+DEE6
                const float X1 = CAR->DEE5 * CAR->DEE6;
                //F5460           X12=X1*X1
                const float X12 = X1 * X1;
                //F5461           AX=(D32-2.0*DY/X1)/X12
                const float AX = ( D32 - 2.0 * DY / X1 ) / X12;
                //F5462           BX=(3.*DY/X1-2.0*D32)/X1
                const float BX = ( 3.0 * DY / X1 - 2.0 * D32 ) / X1;
                //F5463           U=EE-XX5
                const float U = EE - CAR->XX5;
                //F5464           U2=U*U
                const float U2 = U * U;
                //F5465           A1(MM,K)=Y0+D32*U+BX*U2+AX*U2*U
                A1[ MM ][ K ] = Y0 + D32 * U + BX * U2 + AX * U2 * U;
                //F5466 !
                //F5467         ELSE
                //F5468           A1(MM,K)=A(3,K)
            } else A1[ MM ][ K ] = CAR->A.getval( 3, K );
            //F5469         ENDIF
            //F5470 !
            //F5471         END DO
        }
        //F5472 !
        //F5473       ELSE
    } else {
        //F5474 !
        //F5475         XXX1=0.75*FL1
        const float XXX1 = 0.75 * CAR->FL1;
        //F5476         XXX2=0.75*FL2
        const float XXX2 = 0.75 * CAR->FL2;
        //F5477         XXX3=1.25*FL2
        const float XXX3 = 1.25 * CAR->FL2;
        //F5478         XXX4=1.25*FL3
        const float XXX4 = 1.25 * CAR->FL3;
        //F5479         DX12=XXX2-XXX1
        const float DX12 = XXX2 - XXX1;
        //F5480         DX23=XXX4-XXX3
        const float DX23 = XXX4 - XXX3;
        //F5481 !
        //F5482         DO K=1,5
        for( int K=1; K<=5; K++ ) {
            //F5483         IF(EE.LE.XXX1)THEN
            if( EE <= XXX1 ) {
                //F5484           A1(MM,K)=A(1,K)
                A1[ MM ][ K ] = CAR->A.getval( 1, K );
                //F5485         ELSE IF(EE.LE.XXX2)THEN
            } else if( EE <= XXX2 ) {
                //F5486           Z1=EE-XXX1
                const float Z1 = EE - XXX1;
                //F5487           Z2=EE-XXX2
                const float Z2 = EE - XXX2;
                //F5488           DY12=A(2,K)-A(1,K)
                const float DY12 = CAR->A.getval( 2, K ) - CAR->A.getval( 1, K );
                //F5489           A1(MM,K)=A(1,K)+(DY12/DX12**3)*Z1*Z1*(Z1-3.0*Z2)
                A1[ MM ][ K ] = CAR->A.getval( 1, K ) + ( DY12 / pow( DX12, 3 ) ) * Z1 * Z1 * ( Z1-3.0*Z2 );
                //F5490         ELSE IF(EE.LE.XXX3)THEN
            } else if( EE <= XXX3 ) {
                //F5491           A1(MM,K)=A(2,K)
                A1[ MM ][ K ] = CAR->A.getval( 2, K );
                //F5492         ELSE IF(EE.LE.XXX4)THEN
            } else if( EE <= XXX4 ) {
                //F5493           Z1=EE-XXX3
                const float Z1 = EE - XXX3;
                //F5494           Z2=EE-XXX4
                const float Z2 = EE - XXX4;
                //F5495           DY23=A(3,K)-A(2,K)
                const float DY23 = CAR->A.getval( 3, K ) - CAR->A.getval( 2, K );
                //F5496           A1(MM,K)=A(2,K)+(DY23/DX23**3)*Z1*Z1*(Z1-3.0*Z2)
                A1[ MM ][ K ] = CAR->A.getval( 2, K ) + ( DY23/pow( DX23, 3 ) ) * Z1 * Z1 * ( Z1-3.0*Z2 );
                //F5497         ELSE
                //F5498           A1(MM,K)=A(3,K)
            } else A1[ MM ][ K ] = CAR->A.getval( 3, K );
            //F5499         ENDIF
            //F5500         END DO
        }
        //F5501 !
        //F5502       ENDIF
    }
    //F5503 !
    //F5504 !  ******************************************************************
    //F5505 !
    //F5506 !  CALCULATE NEW PARTIAL CONCS AND NEW CONCENTRATION.
    //F5507 !
    //F5508       DELC=0.0
    DELC = 0.0;
    //F5509       DO J=1,5
    for( int J=1; J<=5; J++ ) {
        //F5510       DELA=A1(MM,J)-AA(MM,J)
        float DELA = A1[ MM ][ J ] - CAR->AA.getval( MM, J );
        //F5511       ABAR=(A1(MM,J)+AA(MM,J))/2.0
        float ABAR = ( A1[ MM ][ J ] + CAR->AA.getval( MM, J ) ) / 2.0;
        //F5512       Z=(TINV(MM,J)-DELA/ABAR)/2.0
        float Z = ( CAR->TINV.getval( MM, J ) - DELA / ABAR ) / 2.0;
        //F5513       DEL=(ABAR*SUMEM1/FACTOR-2.0*Z*CPART(MM,J))/(1.0+Z)
        float DEL = ( ABAR * *SUMEM1/CAR->FACTOR - 2.0 * Z * CAR->CPART[ MM ][ J ] )/( 1.0+Z );
        //F5514       CPART(MM,J)=CPART(MM,J)+DEL
        CAR->CPART[ MM ][ J ] += DEL;
        //F5515       DELC=DELC+DEL
        DELC += DEL;
        //F5516       FLUX=SUMEM1-FACTOR*DELC
        *FLUX = *SUMEM1 - CAR->FACTOR * DELC;
        //F5517       AA(MM,J)=A1(MM,J)
        CAR->AA.setval( A1[ MM ][ J ], MM, J );
        //F5518       END DO
    }
    //F5519       C1=C+DELC
    *C1 = C + DELC;
    //F5520       CBAR=(C1+C)/2.0
    CBAR = ( *C1 + C ) / 2.0;
    //F5521 !
    //F5522 !  CALCULATE ANNUAL CHANGES IN ATMOSPHERIC MASS
    //F5523 !
    //F5524       DELM=FACTOR*(C1-C)
    *DELM = CAR->FACTOR * ( *C1 - C );
    //F5525 !
    //F5526       RETURN
    //F5527       END
    
    f_exit( __func__ );
} // carbon
//F5528 !
//F5529 !  *******************************************
//F5530 !
//F5531       SUBROUTINE HISTORY(JJJ,CO2,CH4,CN2O,eso2,eso21, &
//F5532       CF4,C2F6,C125,C134A,C143A,C227,C245,CSF6)
void history( const int JJJ, float* CO2, float* CH4, float* CN2O, float* eso2, float* eso21,
             float* CF4, float* C2F6, float* C125, float* C134A, float* C143A, float* C227, float* C245, 
             float* CSF6, COBS_block* COBS, float ES1990 )
{
    f_enter( __func__ );
    //F5533       IMPLICIT REAL*4 (a-h,o-z), Integer (I-N)
    //F5534 !
    //F5535 !  THIS SUBROUTINE CALCULATES THE CONCS UP TO AND INCLUDING
    //F5536 !   1990 USING FITS TO OBSERVED DATA.
    //F5537 !  Concs are end of year values. Conc(1) is end of 1765 etc.
    //F5538 !  SO2 emissions values are for whole year, assumed to apply to
    //F5539 !   midpoint of year.
    //F5540 !  FOR THE KYOTO PROTOCOL GASES, ONLY THE 1990 VALUE IS GIVEN,
    //F5541 !   SINCE CURRENT VERSION OF CODE DOES NOT ALLOW CONCS FOR THESE
    //F5542 !   TO BE OUTPUT. THE 1990 VALUE IS NEEDED TO INITIALIZE CONC
    //F5543 !   CALCULATIONS FOR 1991 ONWARDS.
    //F5544 !
    //F5545       COMMON/COBS/COBS(0:236)
    //F5546       common /Sulph/S90DIR,S90IND,S90BIO,ENAT,ES1990,ECO90,FOC90,IFOC
    //F5547 !
    //F5548 !  ****************** HALOCARBONS AND RELATED SPECIES
    //F5549 !
    //F5550 !  SINCE PRE-1990 VALUES OF THESE 'HALO' CONCS ARE NOT USED, 
    //F5551 !   SUBROUTINE HISTORY SETS ALL THESE PRE-1990 VALUES TO THEIR
    //F5552 !   1990 LEVEL
    //F5553 !
    //F5554       CF4    =  69.7
    *CF4     = 69.7;
    //F5555       C2F6   =   3.6
    *C2F6    = 3.6;
    //F5556       C125   =    .0
    //F5557       C134A  =    .0
    //F5558       C143A  =    .0
    //F5559       C227   =    .0
    //F5560       C245   =    .0
    *C125 = *C134A = *C143A = *C227 = *C245 = 0.0;
    //F5561       CSF6   =   3.2
    *CSF6    = 3.2;
    //F5562 !
    //F5563 !  ****************** CO2
    //F5564 !
    //F5565 !  CO2 : END OF YEAR CONCS ARE SPECIFIED IN A DATA STATEMENT IN
    //F5566 !   BLOCK DATA, IN ARRAY COBS. THIS ARRAY IS THE IPCC DATA SET
    //F5567 !   GENERATED BY ENTING AND WIGLEY FOR THE IPCC CONC STABILIZATION
    //F5568 !   EXERCISE. IMPLEMENTED ON DEC 30, 1993, UPDATED MAR 11, 1995
    //F5569 !
    //F5570       CO2=COBS(JJJ)
    *CO2 = COBS->COBS[ JJJ ];
    //F5571 !
    //F5572 !  ****************** CH4
    //F5573 !
    //F5574       Y=JJJ+1764.0
    const float Y = JJJ + 1764.0;
    float Y1, YY, CONC0, A, B, D, YY2;
    //F5575 !
    //F5576 !  Y CORRESPS TO END OF YEAR. E.G. IF JJJ=236, Y=2000.0, SO CONC
    //F5577 !   OUTPUT IS VALUE AT END OF YEAR 2000.
    //F5578 !
    //F5579 !  NEW CH4 (AUG. 2000). BEGINS WITH 700ppbv AT START OF 1750 (I.E.,
    //F5580 !   Y=1749.0) AND GOES TO 1100ppbv AT END OF 1940 (I.E., Y=1940.0).
    //F5581 !
    //F5582       IF(Y.LE.1940.0)THEN
    if( Y <= 1940.0 ) {
        //F5583 !
        //F5584         Y1=Y-1749.
        Y1 = Y - 1749.0;
        //F5585         CH4=700.0+400.0*Y1*Y1/(191.0**2)
        *CH4 = 700.0+400.0*Y1*Y1/(191.0*191.0);
        //F5586 !
        //F5587       ELSE
    } else {
        //F5588 !
        //F5589         IF((Y.LE.1970.0).AND.(Y.GT.1940.0))THEN
        if( Y <= 1970.0 && Y > 1940.0 ) {
            //F5590           YY=Y-1940.
            //F5591           CONC0=1100.0
            //F5592           A=4.1885
            //F5593           B=0.26643
            //F5594           D=-0.0010469
            YY = Y - 1940.0; CONC0 = 1100.0; A = 4.1885; B = 0.26643; D = -0.0010469;
            //F5595         ENDIF
        }
        //F5596 !
        //F5597         IF((Y.LE.1980.0).AND.(Y.GT.1970.0))THEN
        if( Y <= 1980.0 && Y > 1970.0 ) {
            //F5598           YY=Y-1969.
            //F5599           CONC0=1420.0
            //F5600           A=17.0
            //F5601           B=-0.5
            //F5602           D=0.03
            YY = Y - 1969.0; CONC0 = 1420.0; A = 17.0; B = -0.5; D = 0.03;
            //F5603         ENDIF
        }
        //F5604 !
        //F5605         IF((Y.LE.1990.0).AND.(Y.GT.1980.0))THEN
        if( Y <= 1990.0 && Y > 1980.0 ) {
            //F5606           YY=Y-1979.
            //F5607           CONC0=1570.0
            //F5608           A=16.0
            //F5609           B=-0.3
            //F5610           D=0.0
            YY = Y - 1979.0; CONC0 = 1570.0; A = 16.0; B = -0.3; D = 0.0;
            //F5611         ENDIF
        }
        //F5612 !
        //F5613         IF((Y.LE.2001.0).AND.(Y.GT.1990.0))THEN
        if( Y <= 2001.0 && Y > 1990.0 ) {
            //F5614           YY=Y-1989.
            //F5615           CONC0=1700.0
            //F5616           A=10.0
            //F5617           B=-1.0
            //F5618           D=0.06
            YY = Y - 1989.0; CONC0 = 1700.0; A = 10.0; B = -1.0; D = 0.06;
            //F5619         ENDIF
        }
        //F5620 !
        //F5621         YY2=YY*YY
        YY2 = YY * YY;
        //F5622         CH4=CONC0+A*YY+B*YY2+D*YY*YY2
        *CH4=CONC0 + A * YY + B * YY2 + D * YY * YY2;
        //F5623 !
        //F5624       ENDIF
    }
    //F5625 !
    //F5626 !  ****************** N2O
    //F5627 !
    //F5628       Y=JJJ+1764.0
    //UNNECESSARY Y=JJJ+1764.0; // This is already assigned above
    //F5629 !
    //F5630 !  Y CORRESPS TO END OF YEAR. E.G. IF JJJ=236, Y=2000.0, SO CONC
    //F5631 !   OUTPUT IS VALUE AT END OF YEAR 2000.
    //F5632 !
    //F5633 !  NEW N2O (AUG. 2000). BEGINS WITH 270ppbv AT START OF 1750 (I.E.,
    //F5634 !   Y=1749.0) AND GOES TO 290ppbv AT END OF 1950 (I.E., Y=1950.0).
    //F5635 !
    //F5636       IF(Y.LE.1950.0)THEN
    if( Y <= 1950.0 ) {
        //F5637 !
        //F5638         Y1=Y-1749.
        Y1 = Y - 1749.0;
        //F5639         CN2O=270.0+20.0*Y1*Y1/(201.0**2)
        *CN2O = 270.0+20.0*Y1*Y1/(201.0*201.0);
        //F5640 !
        //F5641       ELSE
    } else {
        //F5642 !
        //F5643         IF((Y.LE.1970.0).AND.(Y.GT.1950.0))THEN
        if( Y <= 1970.0 && Y > 1950.0 ) {
            //F5644 !
            //F5645           YY=Y-1950.
            //F5646           CONC0=290.0
            //F5647           A=0.199
            //F5648           B=-0.0083435
            //F5649           D=0.00061685
            YY = Y - 1950.0; CONC0 = 290.0; A = 0.199; B = -0.0083435; D = 0.00061685;
            //F5650         ENDIF
        }
        //F5651 !
        //F5652         IF((Y.LE.1980.0).AND.(Y.GT.1970.0))THEN
        if( Y <= 1980.0 && Y > 1970.0 ) {
            //F5653           YY=Y-1969.
            //F5654           CONC0=295.0
            //F5655           A=0.55
            //F5656           B=0.005
            //F5657           D=0.0
            YY = Y - 1969.0; CONC0 = 295.0; A = 0.55; B = 0.005; D = 0.0;
            //F5658         ENDIF
        }
        //F5659 !
        //F5660         IF((Y.LE.1990.0).AND.(Y.GT.1980.0))THEN
        if( Y <= 1990.0 && Y > 1980.0 ) {
            //F5661           YY=Y-1979.
            //F5662           CONC0=301.0
            //F5663           A=0.65
            //F5664           B=0.005
            //F5665           D=0.0
            YY = Y - 1979.0; CONC0 = 301.0; A = 0.65; B = 0.005; D = 0.0;
            //F5666         ENDIF
        }
        //F5667 !
        //F5668         IF((Y.LE.2001.0).AND.(Y.GT.1990.0))THEN
        if( Y <= 2001.0 && Y > 1990.0 ) {
            //F5669           YY=Y-1989.
            //F5670           CONC0=308.0
            //F5671           A=0.75
            //F5672           B=0.01
            //F5673           D=-0.0005
            YY = Y - 1989.0; CONC0 = 308.0; A = 0.75; B = 0.01; D = -0.0005;
            //F5674         ENDIF
        }
        //F5675 !
        //F5676         YY2=YY*YY
        YY2 = YY*YY;
        //F5677         CN2O=CONC0+A*YY+B*YY2+D*YY*YY2
        *CN2O = CONC0+A*YY+B*YY2+D*YY*YY2;
        //F5678 !
        //F5679       ENDIF
    }
    //F5680 !
    //F5681 !  ****************** SO2 EMISSIONS
    //F5682 !
    //F5683 !  SO2 EMISSIONS : NEED TO CALC JJJ & JJJ+1 VALUES.
    //F5684 !   eso2 IS THE VALUE FOR YEAR J, NOMINALLY A MID-YEAR VALUE.
    //F5685 !   eso21 IS THE VALUE FOR YEAR J+1, NEEDED TO CALCULATE AN
    //F5686 !   EFFECTIVE END OF YEAR VALUE IN DETERMINING FORCING.
    //F5687 !  EXTENDED TO 2000 WITH SRES VALUES (AUG. 2000)
    //F5688 !
    //F5689       DO I=0,1
    for( int I=0; I<=1; I++) {
        //F5690       J=jjj+I
        int J = JJJ+I;
        //F5691       ymid = J+1764.
        int ymid = J+1764;
        float ee;
        //F5692       if(ymid.lt.1860.) then
        //F5693         ee = 0.0
        if( ymid < 1860 ) ee = 0.0;
        //F5694       else if(ymid.lt.1953) then
        //F5695         ee = 35.0*(ymid-1860.)/93.0
        else if ( ymid < 1953 ) ee = 35.0*(ymid-1860)/93.0;
        //F5696       else if(ymid.lt.1973) then
        //F5697         ee = 35.0+33.0*(ymid-1953.)/20.
        else if ( ymid < 1973 )  ee = 35.0+33.0*( ymid-1953 )/20.0;
        //F5698       else IF(YMID.LT.1990)THEN
        //F5699         ee = 68.0+(ES1990-68.0)*(ymid-1973.)/17.
        else if ( ymid < 1990 ) ee = 68.0+( ES1990-68.0 )*( ymid-1973 )/17.0;
        //F5700       else IF(YMID.LT.2000)THEN
        //F5701         ee = ES1990-1.876*(ymid-1990.)/10.
        else if ( ymid < 2000 ) ee = ES1990-1.876*( ymid-1990 )/10.0;
        //F5702       endif
        //F5703       IF(I.EQ.0)eso2=ee
        if( I==0 ) *eso2=ee;
        //F5704       IF(I.EQ.1)eso21=ee
        if( I==1 ) *eso21=ee;
        //F5705       END DO
    }
    //F5706 !
    //F5707       RETURN
    //F5708       END
    f_exit( __func__ );
} // history
//F5709 !
//F5710 !  *******************************************
//F5711 !
//F5712 !     CALL METHANE(CH4B(J-1),EECH4,DENOX,DECO,DEVOC,CH4B(J),
//F5713 !  &  T00MID,TAUBEST,SCH4,ANOX,ACO,AVOC)
//F5714 !
//F5715       SUBROUTINE METHANE(ICH4F,CPREV,E,DEN,DEC,DEV,CONC, &
//F5716       TAU00,TAUOUT,S,AANOX,AACO,AAVOC,TEMP)
void methane( float ICH4F, float CPREV, float E, float DEN, float DEC, float DEV,
             float* CONC, float TAU00, float* TAUOUT, float S, float AANOX, 
             float AACO, float AAVOC, float TEMP, METH4_block* METH4 )
{
    f_enter( __func__ );
    //F5717       IMPLICIT REAL*4 (a-h,o-z), Integer (I-N)
    //F5718 !
    //F5719 !  ********************************************************************
    //F5720 !
    //F5721 !  METHANE MODEL USING PRATHER'S METHOD (FOR TAR)
    //F5722 !  MODIFIED TO BEGIN IN 2000
    //F5723 !  TAR METHOD :
    //F5724 !   dOH/0H = S*dB/B + AANOX*dEN +AACO*dEC +AAVOC*dEV   (Table 4.11)
    //F5725 !    WHERE S = -0.32 (Table 4.11) (Table 4.2 GIVES -0.34 ??)
    //F5726 !   dTAU/TAU = GAM * dOH/OH (GAM=-1.145, FROM DATA IN Table 4.2)
    //F5727 !    GAM SPECIFIED IN METH3 COMMON BLOCK
    //F5728 !
    //F5729 !  NOTE THAT CONC AND TAU VALUES ARE END-OF-YEAR VALUES
    //F5730 !
    //F5731 !  ********************************************************************
    //F5732 !
    //F5733       COMMON /METH4/GAM,TAUOTHER,BBCH4,CM00
    //F5734       Real TauSave(10)
    float TauSave[ 10+1 ];
    //F5735 !  ************************************************************
    //F5736 !
    //F5737 !  METHANE CONC PROJECTION.
    //F5738 !
    //F5739 !  FIRST ITERATION
    //F5740 !
    //F5741       TAU0=TAU00
    //UNUSED float TAU0 = TAU00;
    //F5742       B=CPREV*BBCH4
    float B = CPREV * METH4->BBCH4;
    //F5743       B00=CM00*BBCH4
    float B00 = METH4->CM00 * METH4->BBCH4;
    //F5744       AAA=EXP(GAM*(AANOX*DEN+AACO*DEC+AAVOC*DEV))
    float AAA = exp( float( METH4->GAM*( AANOX * DEN + AACO * DEC + AAVOC * DEV ) ) );
    //F5745       X=GAM*S
    float X = METH4->GAM * S;
    //F5746       U=TAU00*AAA
    float U = TAU00 * AAA;
    //F5747 	TauSave(1) = U
    TauSave[ 1 ] = U;
    //F5748 !
    //F5749 !  FIRST ITERATION
    //F5750 !
    //F5751       BBAR=B
    float BBAR = B;
    //F5752       TAUBAR=U*((BBAR/B00)**X)
    float TAUBAR = U * ( pow( BBAR/B00, X ) );
    //F5753       IF(ICH4F.EQ.1)TAUBAR=TAU00/(TAU00/TAUBAR+0.0316*TEMP)
    if( ICH4F == 1 ) TAUBAR = TAU00 / ( TAU00/TAUBAR + 0.0316*TEMP );
    //F5754       DB1=E-BBAR/TAUBAR-BBAR/TAUOTHER
    float DB1 = E - BBAR / TAUBAR - BBAR / METH4->TAUOTHER;
    //F5755       B1=B+DB1
    float B1 = B + DB1;
    //F5756 !
    //F5757 !  SECOND ITERATION
    //F5758 !
    //F5759       BBAR=(B+B1)/2.0
    BBAR = ( B + B1 ) / 2.0;
    //F5760       TAUBAR=U*((BBAR/B00)**X)
    TAUBAR = U * pow( BBAR/B00, X );
    //F5761       TAUBAR=TAUBAR*(1.0-0.5*X*DB1/B)
    TAUBAR *= ( 1.0 - 0.5 * X * DB1/B );
    //F5762       IF(ICH4F.EQ.1)TAUBAR=TAU00/(TAU00/TAUBAR+0.0316*TEMP)
    if( ICH4F == 1 ) TAUBAR = TAU00 / (TAU00 / TAUBAR + 0.0316 * TEMP );
    //F5763       DB2=E-BBAR/TAUBAR-BBAR/TAUOTHER
    float DB2 = E - BBAR / TAUBAR - BBAR / METH4->TAUOTHER;
    //F5764       B2=B+DB2
    float B2 = B + DB2;
    //F5765 !
    //F5766 !  THIRD ITERATION
    //F5767 !
    //F5768       BBAR=(B+B2)/2.0
    BBAR = ( B + B2 ) / 2.0;
    //F5769       TAUBAR=U*((BBAR/B00)**X)
    TAUBAR = U * pow( BBAR/B00, X );
    //F5770       TAUBAR=TAUBAR*(1.0-0.5*X*DB2/B)
    TAUBAR *= ( 1.0 - 0.5 * X * DB2/B );
    //F5771       IF(ICH4F.EQ.1)TAUBAR=TAU00/(TAU00/TAUBAR+0.0316*TEMP)
    if( ICH4F == 1 ) TAUBAR = TAU00 / (TAU00 / TAUBAR + 0.0316 * TEMP );
    //F5772       DB3=E-BBAR/TAUBAR-BBAR/TAUOTHER
    float DB3 = E - BBAR / TAUBAR - BBAR / METH4->TAUOTHER;
    //F5773       B3=B+DB3
    float B3 = B + DB3;
    //F5774     TauSave(2) = U*((BBAR/B00)**X)
    TauSave[ 2 ] = U * pow( BBAR/B00, X );
    //F5775 	TauSave(3) = TAUBAR
    TauSave[ 3 ] = TAUBAR;
    //F5776 !
    //F5777 !  FOURTH ITERATION
    //F5778 !
    //F5779       BBAR=(B+B3)/2.0
    BBAR = ( B + B3 ) / 2.0;
    //F5780       TAUBAR=U*((BBAR/B00)**X)
    TAUBAR = U * pow( BBAR/B00, X );
    //F5781       TAUBAR=TAUBAR*(1.0-0.5*X*DB3/B)
    TAUBAR *= ( 1.0 - 0.5 * X * DB3/B );
    //F5782       IF(ICH4F.EQ.1)TAUBAR=TAU00/(TAU00/TAUBAR+0.0316*TEMP)
    if( ICH4F == 1 ) TAUBAR = TAU00 / (TAU00 / TAUBAR + 0.0316 * TEMP );
    //F5783       DB4=E-BBAR/TAUBAR-BBAR/TAUOTHER
    float DB4 = E - BBAR / TAUBAR - BBAR / METH4->TAUOTHER;
    //F5784       B4=B+DB4
    float B4 = B + DB4;
    //F5785 	TauSave(4) = U*((BBAR/B00)**X)
    TauSave[ 4 ] = U * pow( BBAR/B00, X );
    //F5786 	TauSave(5) = TAUBAR
    TauSave[ 5 ] = TAUBAR;
    //F5787 !
    //F5788 !  LIFETIME AND CONCENTRATION AT END OF STEP
    //F5789 !
    //F5790       TAUOUT=TAUBAR
    *TAUOUT = TAUBAR;
    //F5791       CONC=B4/BBCH4
    *CONC = B4 / METH4->BBCH4;
    //F5792 !
    //F5793       RETURN
    //F5794       END
    f_exit( __func__ );
} // methane
//F5795 !
//F5796 !  ********************************************************************
//F5797 !
//F5798       SUBROUTINE NITROUS(C,CP,CPP,E,C1)
void nitrous( float C, float CP, float CPP, float E, float* C1, TauNitr_block* TauNitr )
{
    f_enter( __func__ );
    //F5799       IMPLICIT REAL*4 (a-h,o-z), Integer (I-N)
    //F5800       common /TauNitr/TN2000,BBN2O,SN2O,CN00,NOFFSET
    //F5801 !
    //F5802 !  *******************************************************
    //F5803 !
    //F5804 !  N2O CONC PROJECTIONS.
    //F5805 !
    //F5806       B=C*BBN2O
    float B = C * TauNitr->BBN2O;
    //F5807       B00=CN00*BBN2O
    float B00 = TauNitr->CN00 * TauNitr->BBN2O;
    //F5808       BBARPREV=0.5*(CP+CPP)*BBN2O
    float BBARPREV = 0.5 * ( CP + CPP ) * TauNitr->BBN2O;
    //F5809       S=SN2O
    float S = TauNitr->SN2O;
    //F5810 !
    //F5811 !  FIRST ITERATION
    //F5812 !
    //F5813       BBAR=B
    float BBAR = B;
    //F5814       TAUBAR=TN2000*((BBAR/B00)**S)
    float TAUBAR = TauNitr->TN2000 * pow( BBAR/B00, S );
    //F5815       DB1=E-BBARPREV/TAUBAR
    float DB1 = E - BBARPREV / TAUBAR;
    //F5816       B1=B+DB1
    float B1 = B + DB1;
    //F5817 !
    //F5818 !  NOTE : TAUBAR IS MIDYR VALUE; B1 IS ENDYR VALUE
    //F5819 !
    //F5820 !  SECOND ITERATION
    //F5821 !
    //F5822       BBAR=(B+B1)/2.0
    BBAR = ( B + B1 ) / 2.0;
    //F5823       TAUBAR=TN2000*((BBAR/B00)**S)
    TAUBAR = TauNitr->TN2000 * pow( BBAR/B00, S );
    //F5824       DB2=E-BBARPREV/TAUBAR
    float DB2 = E - BBARPREV / TAUBAR;
    //F5825       B2=B+DB2
    float B2 = B + DB2;
    //F5826 !
    //F5827 !  THIRD ITERATION
    //F5828 !
    //F5829       BBAR=(B+B2)/2.0
    BBAR = ( B + B2 ) / 2.0;
    //F5830       TAUBAR=TN2000*((BBAR/B00)**S)
    TAUBAR = TauNitr->TN2000 * pow( BBAR/B00, S );
    //F5831       DB3=E-BBARPREV/TAUBAR
    float DB3 = E - BBARPREV / TAUBAR;
    //F5832       B3=B+DB3
    float B3 = B + DB3;
    //F5833 !
    //F5834 !  FOURTH ITERATION
    //F5835 !
    //F5836       BBAR=(B+B3)/2.0
    BBAR = ( B + B3 ) / 2.0;
    //F5837       TAUBAR=TN2000*((BBAR/B00)**S)
    TAUBAR = TauNitr->TN2000 * pow( BBAR/B00, S );
    //F5838       DB4=E-BBARPREV/TAUBAR
    float DB4 = E - BBARPREV / TAUBAR;
    //F5839       B4=B+DB4
    float B4 = B + DB4;
    //F5840       C1=B4/BBN2O
    *C1 = B4 / TauNitr->BBN2O;
    //F5841 !
    //F5842       RETURN
    //F5843       END
    f_exit( __func__ );
} // nitrous
//F5844 !
//F5845 !  ***************************************
//F5846 !
//F5847       SUBROUTINE INTERP(N,ISTART,IY,X,Y)
void interp( int N, int ISTART, int IY[], float X[], magicc_array* Y, int KEND )
{
    f_enter( __func__ );
    //F5848       IMPLICIT REAL*4 (a-h,o-z), Integer (I-N)
    //F5849 !
    //F5850       parameter (iTp =700)
    //F5851 !
    //F5852       common /Limits/KEND
    //F5853 !
    //F5854       DIMENSION IY(100),X(100),Y(226:iTp+1)
    //F5855 !
    //F5856       IEND=ISTART+IY(N)
    const int IEND = ISTART + IY[ N ];
    //F5857       DO I=0,IY(N)-1
    for( int I=0; I<=IY[ N ]-1; I++ ) {
        //F5858         DO K=1,N
        for( int K=1; K<=N; K++ ) {
            //F5859           IF(I.GE.IY(K).AND.I.LT.IY(K+1))THEN
            if( I >= IY[ K ] && I < IY[ K+1 ] ) {
                //F5860             J=I+ISTART
                int J = I + ISTART;
                //F5861             Y(J)=X(K)+(I-IY(K))*(X(K+1)-X(K))/(IY(K+1)-IY(K))
                (*Y).setval( X[ K ] + (I-IY[ K ] ) * ( X[ K+1 ]-X[ K ] ) / ( IY[ K+1 ]-IY[ K ] ), J );
                //F5862           ENDIF
            }
            //F5863         END DO
        } // for K
        //F5864       END DO
    } // for i
    //F5865       Y(IEND)=X(N)
    (*Y).setval( X[ N ], IEND );
    //F5866 !
    //F5867 ! If last year in profile (relative to 1990) not KEND, then assume
    //F5868 !  constant emissions from last year specified in emissions profile
    //F5869 !  to KEND.
    //F5870 !
    //F5871       if(iy(n).lt.KEND) then
    if( IY[ N ] < KEND )
        //F5872         do i=iend+1,KEND
        for( int i=IEND+1; i<=KEND; i++) {
            //F5873           y(i) = x(n)
            (*Y).setval( X[ N ], i );
        }
    //F5874         end do
    //F5875       end if
    //F5876 !
    //F5877       RETURN
    //F5878       END
    f_exit( __func__ );
} // interp
//F5879 !
//F5880 !*************************************************
//F5881 !
//F5882       SUBROUTINE SULPHATE(JY,ESO2,ESO21,ECO,QSO2,QDIR,QFOC,QMN)
void sulphate( const int JY, float ESO2, float ESO21, float ECO, float* QSO2, 
              float* QDIR, float* QFOC, float* QMN, Sulph_block* Sulph )
{
    f_enter( __func__ );
    //F5883       IMPLICIT REAL*4 (a-h,o-z), Integer (I-N)
    //F5884 !
    //F5885       parameter (iTp=740)
    //F5886 !
    //F5887       common /Sulph/S90DIR,S90IND,S90BIO,ENAT,ES1990,ECO90,FOC90,IFOC
    //F5888 !
    //F5889 !  DIRECT AND INDIRECT SULPHATE FORCING
    //F5890 !
    //F5891 !  Tall stack effect factor
    //F5892 !
    //F5893       if(jy.lt.186) then
    float f;
    if( JY < 186 ) {
        //F5894         f = 0.7
        f = 0.7;
        //F5895       else if(jy.lt.206) then
    } else if( JY < 206 ) {
        //F5896         f = 0.7 + 0.3*(jy-186)/20.
        f = 0.7 + 0.3 * ( JY-186 ) / 20.0;
        //F5897       else
    } else {
        //F5898         f = 1.0
        f = 1.0;
        //F5899       endif
    }
    //F5900 !
    //F5901       ky = jy + 1
    float f1;
    const int ky = JY + 1;
    //F5902       if(ky.lt.186) then
    if( ky < 186 ) {
        //F5903         f1 = 0.7
        f1 = 0.7;
        //F5904       else if(ky.lt.206) then
    } else if( ky < 206 ) {
        //F5905         f1 = 0.7 + 0.3*(ky-186)/20.
        f1 = 0.7 + 0.3 * ( ky-186 ) / 20.0;
        //F5906       else
    } else {
        //F5907         f1 = 1.0
        f1 = 1.0;
        //F5908       endif
    }
    //F5909 !
    //F5910 !  Calculate end of year emissions and ditto corrected for tall
    //F5911 !   stack effect
    //F5912 !
    //F5913       eraw=(eso2+eso21)/2.0
    //UNUSED const float eraw = ( ESO2 + ESO21) /  2.0;
    //F5914       e = (f*eso2+f1*eso21)/2.
    const float e = ( f * ESO2 + f1 * ESO21) / 2.0;
    //F5915 !
    //F5916 !  Calculate global forcing at end of year.
    //F5917 !
    //F5918 !  Initialise SO2 parameters.  ORIGINALLY USED VALUES BASED ON
    //F5919 !   CHARLSON ET AL. ALTERED TO FIT IPCC94 ON OCT 16, 1994.
    //F5920 !   FURTHER MODIFIED ON FEB 24, 1995  FOR IPCC95
    //F5921 !  MOVED FROM SUBROUTINE INIT TO SUBROUTINE SULPHATE ON 3/10/95
    //F5922 !   WITH ASO2 AND BSO2 ELIMINATED
    //F5923 !
    //F5924       qdir   = e*s90dir/ES1990
    *QDIR = e * Sulph->S90DIR / Sulph->ES1990;
    //F5925       qindir = s90ind*(alog(1.0+e/ENAT))/(alog(1.0+ES1990/ENAT))
    const float qindir = Sulph->S90IND * ( log( float( 1.0+e/Sulph->ENAT ) ) ) / ( log( float( 1.0+Sulph->ES1990/Sulph->ENAT ) ) );
    //F5926 !
    //F5927       qso2   =  qdir+qindir
    *QSO2 = *QDIR + qindir;
    //F5928 !
    //F5929 !********************************************************
    //F5930 !
    //F5931 !  FOSSIL ORGANIC PLUS BLACK CARBON
    //F5932 !   HISTORY : FOR SIMPLICITY, SCALE WITH SO2 EMISSIONS
    //F5933 !   FUTURE : IF IFOC=0, QFOC CONSTANT AT 1990 LEVEL
    //F5934 !            IF IFOC=1, SCALE WITH SO2 EMISSIONS
    //F5935 !            IF IFOC=2, SCALE WITH CO EMISSIONS
    //F5936 !
    //F5937       IF(JY.LE.226)QFOC=E*FOC90/ES1990
    if( JY <= 226 ) *QFOC = e * Sulph->FOC90 / Sulph->ES1990;
    //F5938       IF(JY.GT.226)THEN
    if( JY > 226 ) {
        //F5939         IF(IFOC.EQ.0)QFOC=FOC90
        if( Sulph->IFOC == 0 ) *QFOC = Sulph->FOC90;
        //F5940         IF(IFOC.EQ.1)QFOC=E*FOC90/ES1990
        if( Sulph->IFOC == 1 ) *QFOC = e * Sulph->FOC90 / Sulph->ES1990;
        //F5941         IF(IFOC.EQ.2)QFOC=ECO*FOC90/ECO90
        if( Sulph->IFOC == 2 ) *QFOC = ECO * Sulph->FOC90 / Sulph->ECO90;
        //F5942       ENDIF
    }
    //F5943 !
    //F5944 !********************************************************
    //F5945 !
    //F5946 !  ADD NITRATE AND MINERAL DUST TO QDIR AND QSO2
    //F5947 !
    //F5948       QNO390=-0.1
    const float QNO390 = -0.1;
    //F5949       QMIN90=-0.1
    const float QMIN90 = -0.1;
    //F5950       IF(JY.LE.226)THEN
    float QNO3, QMIN;
    if( JY <= 226 ) {
        //F5951         QNO3=QNO390*FLOAT(JY)/226.0
        QNO3 = QNO390 * float( JY ) / 226.0;
        //F5952         QMIN=QMIN90*FLOAT(JY)/226.0
        QMIN = QMIN90 * float( JY ) / 226.0;
        //F5953       ELSE
    } else {
        //F5954         QNO3=QNO390
        QNO3 = QNO390;
        //F5955         QMIN=QMIN90
        QMIN = QMIN90;
        //F5956       ENDIF
    }
    //F5957       QMN=QNO3+QMIN
    *QMN = QNO3 + QMIN;
    //F5958 !
    //F5959 !********************************************************
    //F5960 !
    //F5961 !  NOTE : BECAUSE THESE FORCINGS MUST BE SPLIT INTO NH/SH AND
    //F5962 !   LAND/OCEAN, QDIR AND QSO2 ARE COMBINED HERE WITH QFOC.
    //F5963 !   QFOC IS STILL TRANSFERRED TO MAIN PROGRAM IN CASE IT NEEDS
    //F5964 !   TO BE SPLIT OFF LATER
    //F5965 !
    //F5966       QDIR=QDIR+QFOC
    *QDIR += *QFOC;
    //F5967       QSO2=QSO2+QFOC
    *QSO2 += *QFOC;
    //F5968 !
    //F5969       RETURN
    //F5970       END
    f_exit( __func__ );
} // sulfate
//F5971 !
//F5972 !  ******************************************************************
//F5973 !
//F5974       SUBROUTINE LAMCALC(Q,FNHL,FSHL,XK,XKH,DT2X,A,LAMOBEST,LAMLBEST)
void lamcalc( float Q, float FNHL, float FSHL, float XK, float XKH, float DT2X, 
             float A, float* LAMOBEST, float* LAMLBEST )
{
    f_enter( __func__ );
    //F5975       IMPLICIT REAL*4 (a-h,o-z), Integer (I-N)
    //F5976 !
    //F5977 ! Revision history:
    //F5978 !  950215 : CONVERTED TO SUBROUTINE FOR STAG.FOR
    //F5979 !  950208 : ITERATION ALGORITHM IMPROVED YET AGAIN
    //F5980 !  950206 : ITERATION ALGORITHM IMPROVED AGAIN
    //F5981 !  950205 : MATRIX INVERSION CORRECTED
    //F5982 !  950204 : ITERATION ALGORITHM IMPROVED
    //F5983 !  950203 : FIRST VERSION OF PROGRAM WRITTEN
    //F5984 !
    //F5985 !  THIS SUBROUTINE CALCULATES LAND AND OCEAN FEEDBACK
    //F5986 !   PARAMETER VALUES GIVEN THE GLOBAL DT2X AND THE LAND TO
    //F5987 !   OCEAN EQUILIBRIUM WARMING RATIO (A).
    //F5988 !  BOTH THE INPUT VALUES OF FNHL AND FSHL, AND THE INPUT XK
    //F5989 !   (XKLO IN MAIN) AND XKH (XKNS IN MAIN) ARE DOUBLE WHAT
    //F5990 !   ARE USED HERE.
    //F5991 !
    //F5992       DIMENSION AEST(100),DIFF(100)
    float AEST[ 100+1 ], DIFF[ 100+1 ];
    //F5993 !
    //F5994       REAL LAMO(100),LAML(100),LAMOBEST,LAMLBEST,LAM,KLO,KNS
    float LAMO[ 100+1 ], LAML[ 100+1 ];
    //F5995 !
    //F5996       KLO=XK/2.0
    const float KLO = XK / 2.0;
    //F5997       KNS=XKH/2.0
    const float KNS = XKH / 2.0;
    //F5998       FNL=FNHL/2.0
    const float FNL = FNHL / 2.0;
    //F5999       FSL=FSHL/2.0
    const float FSL = FSHL / 2.0;
    //F6000 !
    //F6001       IMAX=40
    const int IMAX = 40;
    //F6002       DLAMO=1.0
    float DLAMO = 1.0;
    //F6003       DIFFLIM=0.001
    const float DIFFLIM = 0.001;
    //F6004 !
    //F6005       FNO=0.5-FNL
    const float FNO = 0.5 - FNL;
    //F6006       FSO=0.5-FSL
    const float FSO = 0.5 - FSL;
    //F6007       FL=FNL+FSL
    const float FL = FNL + FSL;
    //F6008       FO=FNO+FSO
    const float FO = FNO + FSO;
    //F6009       FRATIO=FO/FL
    const float FRATIO = FO / FL;
    //F6010 !
    //F6011       DT2XO=DT2X/(FO+A*FL)
    //UNUSED const float DT2XO = DT2X / ( FO + A * FL );
    //F6012       DT2XL=A*DT2XO
    //UNUSED const float DT2XL = A * DT2XO;
    //F6013       LAM=Q/DT2X
    const float LAM= Q / DT2X;
    //F6014       LAMO(1)=LAM
    LAMO[ 1 ] = LAM;
    //F6015       LAMO(2)=LAM+DLAMO
    LAMO[ 2 ] = LAM + DLAMO;
    //F6016 !
    //F6017       IFLAG=0
    int IFLAG = 0;
    //F6018       DO 1 I=1,IMAX
    int I;
    for( I=1; I<=IMAX; I++ ) {
        //F6019 !
        //F6020       LAML(I)=LAM+FRATIO*(LAM-LAMO(I))/A
        LAML[ I ] = LAM + FRATIO * ( LAM-LAMO[ I ] ) / A;
        //F6021 !
        //F6022 !  SOLVE FOR NH/SH OCEAN/LAND TEMPS
        //F6023 !   FIRST SPECIFY COEFFICIENT MATRIX : A(I,J)
        //F6024 !
        //F6025       A11= FNO*LAMO(I)+KLO+KNS
        float A11 = FNO * LAMO[ I ] + KLO + KNS;
        //F6026       A12=-KLO
        float A12 = -KLO;
        //F6027       A13=-KNS
        float A13 = -KNS;
        //F6028       A14= 0.0
        //UNUSED float A14 = 0.0;
        //F6029       A22= FNL*LAML(I)+KLO
        float A22 = FNL * LAML[ I ] + KLO;
        //F6030       A23= 0.0
        //UNUSED float A23 = 0.0;
        //F6031       A24= 0.0
        //UNUSED float A24 = 0.0;
        //F6032       A33= FSO*LAMO(I)+KLO+KNS
        float A33 = FSO * LAMO[ I ] + KLO + KNS;
        //F6033       A34= A12
        //UNUSED float A34 = A12;
        //F6034       A44= FSL*LAML(I)+KLO
        float A44 = FSL * LAML[ I ] + KLO;
        //F6035 !
        //F6036 !  CALCULATE INVERSE OF COEFFICIENT MATRIX : B(I,J)
        //F6037 !   FIRST DETERMINE DETERMINANT OF A(I,J) MATRIX
        //F6038 !
        //F6039       C1 = A11*A22-A12*A12
        float C1 = A11*A22-A12*A12;
        //F6040       C2 = A33*A44-A12*A12
        float C2 = A33*A44-A12*A12;
        //F6041       C3 = A22*A13
        float C3 = A22*A13;
        //F6042       C4 = A44*A13
        float C4 = A44*A13;
        //F6043       DET= C1*C2-C3*C4
        float DET= C1*C2-C3*C4;
        //F6044 !
        //F6045       B11= A22*C2/DET
        float B11= A22*C2/DET;
        //F6046       B12=-A12*C2/DET
        float B12=-A12*C2/DET;
        //F6047       B13=-A44*C3/DET
        float B13=-A44*C3/DET;
        //F6048       B14= A12*C3/DET
        float B14= A12*C3/DET;
        //F6049       B22= (A11*C2-A13*C4)/DET
        float B22= (A11*C2-A13*C4)/DET;
        //F6050       B23= A12*C4/DET
        float B23= A12*C4/DET;
        //F6051       B24=-A12*A12*A13/DET
        float B24=-A12*A12*A13/DET;
        //F6052       B33= A44*C1/DET
        float B33= A44*C1/DET;
        //F6053       B34=-A12*C1/DET
        float B34=-A12*C1/DET;
        //F6054       B44= (A33*C1-A13*C3)/DET
        float B44= (A33*C1-A13*C3)/DET;
        //F6055 !
        //F6056 !  CALCULATE ESTIMATED NH/SH OCEAN/LAND EQUILIBRIUM TEMPS
        //F6057 !
        //F6058       TNO=(B11*FNO+B12*FNL+B13*FSO+B14*FSL)*Q
        float  TNO=(B11*FNO+B12*FNL+B13*FSO+B14*FSL)*Q;
        //F6059       TNL=(B12*FNO+B22*FNL+B23*FSO+B24*FSL)*Q
        float TNL=(B12*FNO+B22*FNL+B23*FSO+B24*FSL)*Q;
        //F6060       TSO=(B13*FNO+B23*FNL+B33*FSO+B34*FSL)*Q
        float TSO=(B13*FNO+B23*FNL+B33*FSO+B34*FSL)*Q;
        //F6061       TSL=(B14*FNO+B24*FNL+B34*FSO+B44*FSL)*Q
        float TSL=(B14*FNO+B24*FNL+B34*FSO+B44*FSL)*Q;
        //F6062 !
        //F6063 !  CALCULATE ESTIMATED OCEAN-MEAN AND LAND-MEAN TEMPS
        //F6064 !
        //F6065       DT2XLE=(TNL*FNL+TSL*FSL)/FL
        float DT2XLE=(TNL*FNL+TSL*FSL)/FL;
        //F6066       DT2XOE=(TNO*FNO+TSO*FSO)/FO
        float DT2XOE=(TNO*FNO+TSO*FSO)/FO;
        //F6067 !
        //F6068 !  CALCULATE ESTIMATED N.H. AND S.H. TEMPERATURES
        //F6069 !
        //F6070       TNH=(TNL*FNL+TNO*FNO)/0.5
        //UNUSED float TNH=(TNL*FNL+TNO*FNO)/0.5;
        //F6071       TSH=(TSL*FSL+TSO*FSO)/0.5
        //UNUSED float TSH=(TSL*FSL+TSO*FSO)/0.5;
        //F6072 !
        //F6073 !  CALCULATE ESTIMATED VALUE OF A
        //F6074 !
        //F6075       AEST(I)=DT2XLE/DT2XOE
        AEST[ I ] = DT2XLE/DT2XOE;
        //F6076       AAAA=AEST(I)
        //UNUSED float AAA = AEST[ I ];
        //F6077       DIFF(I)=A-AEST(I)
        DIFF[ I ] = A - AEST[ I ];
        //F6078 !
        //F6079 !  TEST DIFF TO DECIDE WHETHER TO END ITERATION LOOP
        //F6080 !
        //F6081       IF(ABS(DIFF(I)).LT.DIFFLIM)GO TO 2
        if( fabs( DIFF[ I ] ) < DIFFLIM ) break;
        //F6082 !
        //F6083       IF(I.GE.2)THEN
        if( I >= 2 ) {
            //F6084         DD=DIFF(I)*DIFF(I-1)
            float DD = DIFF[ I ] * DIFF[ I-1 ];
            //F6085 !
            //F6086         IF(DD.LT.0.0)THEN
            if( DD < 0.0 ) {
                //F6087           IFLAG=1
                IFLAG = 1;
                //F6088         ELSE
            } else {
                //F6089           IF(ABS(DIFF(I)).GT.ABS(DIFF(I-1)))DLAMO=-DLAMO
                if( fabs( DIFF[ I ] ) > fabs( DIFF[ I-1 ]) ) DLAMO = -DLAMO;
                //F6090           LAMO(I+1)=LAMO(I)+DLAMO
                LAMO[ I+1 ] = LAMO[ I ] + DLAMO;
                //F6091         ENDIF
            }
            //F6092 !
            //F6093         IF(IFLAG.EQ.1)THEN
            if( IFLAG == 1 ) {
                //F6094           IF(DD.LT.0.0)THEN
                if( DD < 0.0 ) {
                    //F6095             RATIO=(LAMO(I)-LAMO(I-1))/(DIFF(I)-DIFF(I-1))
                    float RATIO = ( LAMO[ I ] - LAMO[ I-1 ] ) / ( DIFF[ I ] - DIFF[ I-1 ] );
                    //F6096             LAMO(I+1)=LAMO(I)-RATIO*DIFF(I)
                    LAMO[ I+1 ] = LAMO[ I ] - RATIO * DIFF[ I ];
                    //F6097           ELSE
                } else {
                    //F6098             RATIO=(LAMO(I)-LAMO(I-2))/(DIFF(I)-DIFF(I-2))
                    float RATIO = ( LAMO[ I ] - LAMO[ I-2 ] ) / ( DIFF[ I ] - DIFF[ I-2 ] );
                    //F6099             LAMO(I+1)=LAMO(I)-RATIO*DIFF(I)
                    LAMO[ I+1 ] = LAMO[ I ] - RATIO * DIFF[ I ];
                    //F6100           ENDIF
                }
                //F6101         ENDIF
            }
            //F6102 !
            //F6103       ENDIF
        }
        //F6104 !
        //F6105    1  CONTINUE
    } // for
    //F6106    2  CONTINUE
    //F6107 !
    //F6108       LAMOBEST=LAMO(I)
    *LAMOBEST = LAMO[ I ];
    //F6109       LAMLBEST=LAML(I)
    *LAMLBEST = LAML[ I ];
    //F6110 !
    //F6111       RETURN
    //F6112       END
    f_exit( __func__ );
} //F lamcalc
//F6113 
//F6114 !*****************************************************************************************
//F6115 ! sjs New routines to get values from MAGICC 
//F6116 !*****************************************************************************************
//F6117 

/*  These functions are called by MAGICC and need a way to extract values from data structures.
 For now, we use some globals that are initialized with a call from CLIMAT.
 */

/* Fundamental difference from Fortran: we're going to get call by GCAM w/o
 CLIMAT having initialized things first. So take care of it here. */
CARB_block* G_CARB = new CARB_block;
TANDSL_block* G_TANDSL = new TANDSL_block;
CONCS_block* G_CONCS = new CONCS_block;
NEWCONCS_block* G_NEWCONCS = new NEWCONCS_block;
STOREDVALS_block* G_STOREDVALS = new STOREDVALS_block;
METH1_block* G_METH1 = new METH1_block;
CAR_block* G_CAR = new CAR_block;
FORCE_block* G_FORCE = new FORCE_block;
JSTART_block* G_JSTART = new JSTART_block;
QADD_block* G_QADD = new QADD_block;
HALOF_block* G_HALOF = new HALOF_block;
NEWPARAMS_block* G_NEWPARAMS = new NEWPARAMS_block;
BCOC_block* G_BCOC = new BCOC_block;




void setLocals( CARB_block* CARB, TANDSL_block* TANDSL, CONCS_block* CONCS, NEWCONCS_block* NEWCONCS, 
                STOREDVALS_block* STOREDVALS, NEWPARAMS_block* NEWPARAMS, BCOC_block* BCOC, 
                METH1_block* METH1, CAR_block* CAR, FORCE_block* FORCE, JSTART_block* JSTART,
                QADD_block* QADD, HALOF_block* HALOF )
{
    f_enter( __func__ );
/*    G_CARB = CARB;
    G_TANDSL = TANDSL;
    G_CONCS = CONCS;
    G_NEWCONCS = NEWCONCS;
    G_STOREDVALS = STOREDVALS; */
    *NEWPARAMS = *G_NEWPARAMS;
    *BCOC = *G_BCOC;
/*    G_METH1 = METH1;
    G_CAR = CAR;
    G_FORCE = FORCE;
    G_JSTART = JSTART;
    G_QADD = QADD;
    G_HALOF = HALOF; */
    f_exit( __func__ );
}

void setGlobals( CARB_block* CARB, TANDSL_block* TANDSL, CONCS_block* CONCS, NEWCONCS_block* NEWCONCS, 
               STOREDVALS_block* STOREDVALS, NEWPARAMS_block* NEWPARAMS, BCOC_block* BCOC, 
               METH1_block* METH1, CAR_block* CAR, FORCE_block* FORCE, JSTART_block* JSTART,
               QADD_block* QADD, HALOF_block* HALOF )
{
    f_enter( __func__ );
    *G_CARB = *CARB;
     *G_TANDSL = *TANDSL;
     *G_CONCS = *CONCS;
     *G_NEWCONCS = *NEWCONCS;
     *G_STOREDVALS = *STOREDVALS;
    *G_NEWPARAMS = *NEWPARAMS;
    *G_BCOC = *BCOC;
    *G_METH1 = *METH1;
     *G_CAR = *CAR;
     *G_FORCE = *FORCE;
     *G_JSTART = *JSTART;
     *G_QADD = *QADD;
     *G_HALOF = *HALOF;
    f_exit( __func__ );
}


//F6118       FUNCTION getCO2Conc( inYear )
float getCO2Conc( int inYear )
{
    f_enter( __func__ );
    assert( G_CARB != NULL && G_TANDSL != NULL && G_CONCS != NULL && G_NEWCONCS != NULL &&
           G_STOREDVALS != NULL && G_NEWPARAMS != NULL && G_BCOC != NULL && 
           G_METH1 != NULL && G_CAR != NULL && G_FORCE != NULL && G_JSTART != NULL && 
           G_QADD != NULL && G_HALOF != NULL );
    //F6119       IMPLICIT REAL*4 (a-h,o-z), Integer (I-N)
    //F6120 ! Expose subroutine co2Conc to users of this DLL
    //F6121 !DEC$ATTRIBUTES DLLEXPORT::getCO2Conc
    //F6122 
    //F6123       parameter (iTp=740)
    //F6124 
    //F6125       COMMON/CARB/CCO2(4,224:iTp),EDGROSS(4,226:iTp),EF(226:iTp+1), &
    //F6126       REGROW(4,226:iTp),PL(4,226:iTp),HL(4,226:iTp),SOIL(4,226:iTp), &
    //F6127       TTT(226:iTp),ESUM(226:iTp),ETOT(4,226:iTp),EDNET90(4), &
    //F6128       FOC(4,226:iTp),co2(0:iTp),CO2SAVE(0:iTp)
    //F6129 
    //F6130 	  REAL*4 getCO2Conc
    //F6131 
    //F6132       IYR = inYear-1990+226
    const int IYR = inYear - 1990 + 226;
    //F6133 
    //F6134       getCO2Conc = CO2( IYR )
    return( G_CARB->CO2[ IYR ] );
    //F6135 
    //F6136       RETURN 
    //F6137 	  END
    f_exit( __func__ );
}
//F6138 	    
//F6139       FUNCTION getSLR( inYear )
float getSLR( const int inYear )
{
    f_enter( __func__ );
    assert( G_CARB != NULL && G_TANDSL != NULL && G_CONCS != NULL && G_NEWCONCS != NULL &&
           G_STOREDVALS != NULL && G_NEWPARAMS != NULL && G_BCOC != NULL && 
           G_METH1 != NULL && G_CAR != NULL && G_FORCE != NULL && G_JSTART != NULL && 
           G_QADD != NULL && G_HALOF != NULL );
    //F6140       IMPLICIT REAL*4 (a-h,o-z), Integer (I-N)
    //F6141 ! Expose subroutine co2Conc to users of this DLL
    //F6142 !DEC$ATTRIBUTES DLLEXPORT::getCO2Conc
    //F6143 
    //F6144       parameter (iTp=740)
    //F6145 
    //F6146       COMMON/TANDSL/TEQU(iTp),TGAV(iTp),TNHO(iTp), &
    //F6147      TSHO(iTp),TNHL(iTp),TSHL(iTp),TDEEP(iTp),TNHAV(iTp),TSHAV(iTp), &
    //F6148      TLAND(iTp),TOCEAN(iTp),TOCN(40),TOCNPREV(40), &
    //F6149      SIP,SGP,SAP,SLI(iTp),SLG(iTp),SLA(iTp),EX(0:iTp),SLT(iTp), &
    //F6150      QTOT(0:iTp),QGH(0:iTp),QOZ(0:iTp),QBIO(0:iTp),SLO(iTp), &
    //F6151      QSO2(0:iTp+1),QDIR(0:iTp+1),QLAND(0:iTp),QMN(0:iTp+1)
    //F6152 
    //F6153 	  REAL*4 getSLR
    //F6154 
    //F6155       IYR = inYear-1990+226
    const int IYR = inYear - 1990 + 226;
    //F6156       ST1=SLT(IYR)
    const float ST1 = G_TANDSL->SLT[ IYR ];
    //F6157       SO1=SLO(IYR)
    const float SO1 = G_TANDSL->SLO[ IYR ];
    //F6158       SLRAW1=ST1-SO1
    const float SLRAW1 = ST1 - SO1;
    //F6159 
    //F6160       getSLR = SLRAW1
    return SLRAW1;
    //F6161 
    //F6162       RETURN 
    //F6163 	  END
    f_exit( __func__ );
}
//F6164 
//F6165       FUNCTION getGHGConc( ghgNumber, inYear )
float GETGHGCONC( int ghgNumber, int inYear )
{
    f_enter( __func__ );
    assert( G_CARB != NULL && G_TANDSL != NULL && G_CONCS != NULL && G_NEWCONCS != NULL &&
           G_STOREDVALS != NULL && G_NEWPARAMS != NULL && G_BCOC != NULL && 
           G_METH1 != NULL && G_CAR != NULL && G_FORCE != NULL && G_JSTART != NULL && 
           G_QADD != NULL && G_HALOF != NULL );
    //F6166       IMPLICIT REAL*4 (a-h,o-z), Integer (I-N)
    //F6167 ! Expose subroutine ghgConc to users of this DLL
    //F6168 !DEC$ATTRIBUTES DLLEXPORT::getGHGConc
    //F6169 
    //F6170       parameter (iTp=740)
    //F6171 
    //F6172       COMMON/CARB/CCO2(4,224:iTp),EDGROSS(4,226:iTp),EF(226:iTp+1), &
    //F6173       REGROW(4,226:iTp),PL(4,226:iTp),HL(4,226:iTp),SOIL(4,226:iTp), &
    //F6174       TTT(226:iTp),ESUM(226:iTp),ETOT(4,226:iTp),EDNET90(4), &
    //F6175       FOC(4,226:iTp),co2(0:iTp),CO2SAVE(0:iTp)
    //F6176 
    //F6177       COMMON/CONCS/CH4(0:iTp),CN2O(0:iTp),ECH4(226:iTp+1), &
    //F6178       EN2O(226:iTp+1),ECO(226:iTp+1),COE(iTp+1),EVOC(226:iTp+1), &
    //F6179       ENOX(226:iTp+1),ESO2(0:iTp+1),ESO2SUM(226:iTp+1), &
    //F6180       ESO21(226:iTp+1),ESO22(226:iTp+1),ESO23(226:iTp+1), &
    //F6181       EBC(226:iTp+1), EOC(226:iTp+1) ! sjs- add BC-OC
    //F6182 !
    //F6183       COMMON/NEWCONCS/CF4(iTp),C2F6(iTp),C125(iTp),C134A(iTp), &
    //F6184      C143A(iTp),C227(iTp),C245(iTp),CSF6(iTp), &
    //F6185      ECF4(226:iTp+1),EC2F6(226:iTp+1),E125(226:iTp+1),E134A(226:iTp+1), &
    //F6186      E143A(226:iTp+1),E227(226:iTp+1),E245(226:iTp+1),ESF6(226:iTp+1)
    //F6187 
    //F6188 	  REAL*4 getGHGConc
    //F6189       INTEGER ghgNumber
    //F6190 
    //F6191       IYR = inYear-1990+226
    const int IYR = inYear - 1990 + 226;
    float returnValue;
    //F6192 	  
    //F6193 	  ! For consistency, make sure indices are same here as in getForcing
    //F6194       select case (ghgNumber)
    switch( ghgNumber ) {
            //F6195       case(1); getGHGConc = CO2( IYR )
        case 1: returnValue = G_CARB->CO2[ IYR ]; break;
            //F6196       case(2); getGHGConc = CH4( IYR )
        case 2: returnValue = G_CONCS->CH4[ IYR ]; break;
            //F6197       case(3); getGHGConc = CN2O( IYR )
        case 3: returnValue = G_CONCS->CN2O[ IYR ]; break;
            //F6198       case(4); getGHGConc = C2F6( IYR )
        case 4: returnValue = G_NEWCONCS->C2F6[ IYR ]; break;
            //F6199       case(5); getGHGConc = C125( IYR )
        case 5: returnValue = G_NEWCONCS->C125[ IYR ]; break;
            //F6200       case(6); getGHGConc = C134A( IYR )
        case 6: returnValue = G_NEWCONCS->C134A[ IYR ]; break;
            //F6201       case(7); getGHGConc = C143A( IYR )
        case 7: returnValue = G_NEWCONCS->C143A[ IYR ]; break;
            //F6202       case(8); getGHGConc = C245( IYR )
        case 8: returnValue = G_NEWCONCS->C245[ IYR ]; break;
            //F6203       case(9); getGHGConc = CSF6( IYR )
        case 9: returnValue = G_NEWCONCS->CSF6[ IYR ]; break;
            //F6204       case(10); getGHGConc = CF4( IYR )
        case 10: returnValue = G_NEWCONCS->CF4[ IYR ]; break;
            //F6205       case(11); getGHGConc = C227( IYR )
        case 11: returnValue = G_NEWCONCS->C227[ IYR ]; break;
            //F6206       case default; getGHGConc = -1.0
        default: returnValue = std::numeric_limits<float>::max();
                cerr << __func__ << " undefined gas " << ghgNumber << flush;
            //F6207       end select;
    }
    //F6208       
    //F6209 
    //F6210       RETURN 
    return( returnValue );
    //F6211 	  END
    f_exit( __func__ );
}
//F6212 	  
//F6213 ! Returns mid-year forcing for a given gas
//F6214       FUNCTION getForcing( iGasNumber, inYear )
float GETFORCING( const int iGasNumber, const int inYear )
{
    f_enter( __func__ );
    assert( G_CARB != NULL && G_TANDSL != NULL && G_CONCS != NULL && G_NEWCONCS != NULL &&
           G_STOREDVALS != NULL && G_NEWPARAMS != NULL && G_BCOC != NULL && 
           G_METH1 != NULL && G_CAR != NULL && G_FORCE != NULL && G_JSTART != NULL && 
           G_QADD != NULL && G_HALOF != NULL );
    //F6215       IMPLICIT REAL*4 (a-h,o-z), Integer (I-N)
    //F6216 ! Expose subroutine getForcing to users of this DLL
    //F6217 !DEC$ATTRIBUTES DLLEXPORT::getForcing
    //F6218 
    //F6219       parameter (iTp=740)
    //F6220 
    //F6221       COMMON /FORCE/qco2(0:iTp),qm(0:iTp),qn(0:iTp),QCFC(0:iTp), &
    //F6222       QMONT(0:iTp),QOTHER(0:iTp),QSTRATOZ(0:iTp),QCH4O3(0:iTp), &
    //F6223       CFC12(0:iTp), QCH4H2O(0:iTp),QBC(0:iTp),QOC(0:iTp)
    //F6224 !
    //F6225       COMMON /HALOF/QCF4_ar(0:iTp),QC2F6_ar(0:iTp),qSF6_ar(0:iTp), &
    //F6226        Q125_ar(0:iTp),Q134A_ar(0:iTp), &
    //F6227        Q143A_ar(0:iTp),Q227_ar(0:iTp),Q245_ar(0:iTp)
    //F6228 
    //F6229       COMMON/TANDSL/TEQU(iTp),TGAV(iTp),TNHO(iTp), &
    //F6230       TSHO(iTp),TNHL(iTp),TSHL(iTp),TDEEP(iTp),TNHAV(iTp),TSHAV(iTp), &
    //F6231       TLAND(iTp),TOCEAN(iTp),TOCN(40),TOCNPREV(40), &
    //F6232       SIP,SGP,SAP,SLI(iTp),SLG(iTp),SLA(iTp),EX(0:iTp),SLT(iTp), &
    //F6233       QTOT(0:iTp),QGH(0:iTp),QOZ(0:iTp),QBIO(0:iTp),SLO(iTp), &
    //F6234       QSO2(0:iTp+1),QDIR(0:iTp+1),QLAND(0:iTp),QMN(0:iTp+1)
    //F6235 
    //F6236       COMMON /JSTART/JSTART,FOSSHIST(0:236),QKYMAG(0:iTp),IGHG, &
    //F6237      QCH4OZ,QFOC(0:iTp),ICO2CORR,TROZSENS
    //F6238 
    //F6239       COMMON/STOREDVALS/ TEMUSER(iTp),QSO2SAVE(0:iTp+1),QDIRSAVE(0:iTp+1), &
    //F6240       KYRREF
    //F6241 
    //F6242 
    //F6243       COMMON /QADD/IQREAD,OrgIQREAD,JQFIRST,JQLAST,QEX(0:iTp),QEXNH(0:iTp), &
    //F6244      QEXSH(0:iTp),QEXNHO(0:iTp),QEXNHL(0:iTp),QEXSHO(0:iTp), &
    //F6245      QEXSHL(0:iTp),IOLDTZ
    //F6246 
    //F6247       COMMON /METH1/emeth(226:iTp),imeth,ch4l(225:iTp),ch4b(225:iTp), &
    //F6248      ch4h(225:iTp),ef4(226:iTp),StratH2O,TCH4(iTp),iO3feed, &
    //F6249      ednet(226:iTp+1),DUSER,FUSER,CORRUSER,CORRMHI,CORRMMID,CORRMLO
    //F6250 
    //F6251 	  REAL*4 getForcing
    //F6252         
    //F6253       IYRQALL=1990 ! Duplicated hard coded value from main routine
    const int IYRQALL = 1990;
    //F6254       IYR = inYear-1990+226
    const int IYR = inYear - 1990 + 226;
    //F6255       IYRP = IYR - 1
    const int IYRP = IYR - 1;
    //F6256       M00=IYRQALL-1764
    const int M00 = IYRQALL - 1764;
    //F6257       M01=M00-1
    const int M01 = M00 - 1;
    //F6258       
    //F6259 ! Calculate mid-year forcing components
    //F6260         QQQCO2 = (QCO2(IYR)+QCO2(IYRP))/2.
    const float QQQCO2 = ( G_FORCE->QCO2[ IYR ] + G_FORCE->QCO2[ IYRP ] ) / 2.0;
    //F6261         QQQM   = (QM(IYR)+QM(IYRP))/2.
    /* const */ float QQQM = ( G_FORCE->QM[ IYR ] + G_FORCE->QM[ IYRP ] ) / 2.0;
    //F6262         QQQN   = (QN(IYR)+QN(IYRP))/2.
    const float QQQN = ( G_FORCE->QN[ IYR ] + G_FORCE->QN[ IYRP ] ) / 2.0;
    //F6263         QQQCFC = (QCFC(IYR)+QCFC(IYRP))/2.
    const float QQQCFC = ( G_FORCE->QCFC[ IYR ] + G_FORCE->QCFC[ IYRP ] ) / 2.0;
    //F6264         QQQOZ  = (QOZ(IYR)+QOZ(IYRP))/2.
    /* const */ float QQQOZ = ( G_TANDSL->QOZ[ IYR ] + G_TANDSL->QOZ[ IYRP ] ) / 2.0;
    //F6265         QQQFOCR  = (QFOC(IYR)     +QFOC(IYRP))     /2.
    const float QQQFOCR = ( G_JSTART->QFOC[ IYR ] + G_JSTART->QFOC[ IYRP ] ) / 2.0;
    //F6266 
    //F6267         QQQSO2 = 0.0
    float QQQSO2 = 0.0;
    //F6268         QQQDIR = 0.0
    float QQQDIR = 0.0;
    //F6269         IF(inYear.GT.1860)THEN
    if( inYear > 1860 ) {
        //F6270           QQQSO2 = (QSO2SAVE(IYR)+QSO2SAVE(IYRP))/2.
        QQQSO2 = ( G_STOREDVALS->QSO2SAVE[ IYR ] + G_STOREDVALS->QSO2SAVE[ IYRP ] ) / 2.0;
        //F6271           QQQDIR = (QDIRSAVE(IYR)+QDIRSAVE(IYRP))/2.
        QQQDIR = ( G_STOREDVALS->QDIRSAVE[ IYR ] + G_STOREDVALS->QDIRSAVE[ IYRP ] ) / 2.0;
        //F6272         ENDIF
    }
    //F6273          QQQIND = QQQSO2-QQQDIR
    //UNUSED const float QQQIND = QQQSO2 - QQQDIR;
    //F6274          DELQFOC = (QFOC(IYR)+QFOC(IYRP))/2.-QQQFOCR
    const float DELQFOC = ( G_JSTART->QFOC[ IYR ] + G_JSTART->QFOC[ IYRP ] ) / 2.0;
    //F6275 !
    //F6276          QQQCO2 = (QCO2(IYR)+QCO2(IYRP))/2.
    //UNNECESSARY const float QQQCO2 = ( FORCE->QCO2[ IYR ] + FORCE->QCO2[ IYRP ] ) / 2.0;
    //F6277          QQQM   = (QM(IYR)+QM(IYRP))/2.
    //UNNECESSARY const float QQQM = ( FORCE->QM[ IYR ] + FORCE->QM[ IYRP ] ) / 2.0;
    //F6278          QQQN   = (QN(IYR)+QN(IYRP))/2.
    //UNNECESSARY const float QQQN = ( FORCE->QN[ IYR ] + FORCE->QN[ IYRP ] ) / 2.0;
    //F6279          QQQCFC = (QCFC(IYR)+QCFC(IYRP))/2.
    //UNNECESSARY const float QQQCFC = ( FORCE->QCFC[ IYR ] + FORCE->QCFC[ IYRP ] ) / 2.0;
    //F6280          QQQOZ  = (QOZ(IYR)+QOZ(IYRP))/2.
    //UNNECESSARY const float QQQOZ = ( TANDSL->QOZ[ IYR ] + TANDSL->QOZ[ IYRP ] ) / 2.0;
    //F6281          QQQFOC = (QFOC(IYR)+QFOC(IYRP))/2.
    //UNNECESSARY const float QQQFOC = ( JSTART->QFOC[ M00 ] + JSTART->QFOC[ M01 ] ) / 2.0;
    //F6282          QQQMN  = (QMN(IYR)+QMN(IYRP))/2.
    const float QQQMN = ( G_TANDSL->QMN[ IYR ] + G_TANDSL->QMN[ IYRP ] ) / 2.0;
    //F6283          
    //F6284          QQQEXTRA = ( QEXNH(IYR)+QEXSH(IYR)+QEXNHO(IYR)+QEXNHL(IYR) + &
    //F6285                       QEXNH(IYRP)+QEXSH(IYRP)+QEXNHO(IYRP)+QEXNHL(IYRP) )/2.
    float QQQEXTRA = ( G_QADD->QEXNH[ IYR ] + G_QADD->QEXSH[ IYR ] + G_QADD->QEXNHO[ IYR ] + G_QADD->QEXNHL[ IYR ] + 
                      G_QADD->QEXNH[ IYRP ] + G_QADD->QEXSH[ IYRP ] + G_QADD->QEXNHO[ IYRP ] + G_QADD->QEXNHL[ IYRP ]  ) / 2.0;
    //F6286 !
    //F6287 ! NOTE SPECIAL CASE FOR QOZ BECAUSE OF NONLINEAR CHANGE OVER 1990
    //F6288 !
    //F6289          IF(IYR.EQ.226)QQQOZ=QOZ(IYR)
    if( IYR == 226 ) QQQOZ = G_TANDSL->QOZ[ IYR ];
    //F6290 !
    //F6291          QQQLAND= (QLAND(IYR)+QLAND(IYRP))/2.
    const float QQQLAND = ( G_TANDSL->QLAND[ IYR ] + G_TANDSL->QLAND[ IYRP ] ) / 2.0;
    //F6292          QQQBIO = (QBIO(IYR)+QBIO(IYRP))/2.
    const float QQQBIO = ( G_TANDSL->QBIO[ IYR ] + G_TANDSL->QBIO[ IYRP ] ) / 2.0;
    //F6293          QQQTOT = QQQCO2+QQQM+QQQN+QQQCFC+QQQSO2+QQQBIO+QQQOZ+QQQLAND &
    //F6294          +QQQMN
    float QQQTOT = QQQCO2 + QQQM + QQQN + QQQCFC + QQQSO2 + QQQBIO + QQQOZ + QQQLAND + QQQMN;
    //F6295 !
    //F6296          QQCH4O3= (QCH4O3(IYR)+QCH4O3(IYRP))/2.
    const float QQCH4O3 = ( G_FORCE->QCH4O3[ IYR ] + G_FORCE->QCH4O3[ IYRP ] ) / 2.0;
    //F6297          QQQM   = QQQM-QQCH4O3
    QQQM -= QQCH4O3;
    //F6298          QQQOZ  = QQQOZ+QQCH4O3
    QQQOZ += QQCH4O3;
    //F6299          QQQD   = QQQDIR-QQQFOC
    //UNUSED const float QQQD = QQQDIR - QQQFOCR;    //CHANGE since QQQFOC = QQQFOCR
    //F6300  
    //F6301          QQQSTROZ= (QSTRATOZ(IYR)+QSTRATOZ(IYRP))/2.
    float QQQSTROZ = ( G_FORCE->QSTRATOZ[ IYR ] + G_FORCE->QSTRATOZ[ IYRP ] ) / 2.0;
    //F6302          IF(IO3FEED.EQ.0)QQQSTROZ=0.0 
    if( G_METH1->IO3FEED == 0 ) QQQSTROZ = 0.0;
    //F6303 !
    //F6304          QQQKYMAG = (QKYMAG(IYR)+QKYMAG(IYRP))/2.
    //UNUSED const float QQQKYMAG = ( JSTART->QKYMAG[ IYR ] + JSTART->QKYMAG[ IYRP ] ) / 2.0;
    //F6305          QQQMONT  = (QMONT(IYR) +QMONT(IYRP)) /2.
    const float QQQMONT = ( G_FORCE->QMONT[ IYR ] + G_FORCE->QMONT[ IYRP ] ) / 2.0;
    //F6306          QQQOTHER = (QOTHER(IYR)+QOTHER(IYRP))/2.
    const float QQQOTHER = ( G_FORCE->QOTHER[ IYR ] + G_FORCE->QOTHER[ IYRP ] ) / 2.0;
    //F6307          QQQKYOTO = QQQKYMAG+QQQOTHER
    //UNUSED const float QQQKYOTO = QQQKYMAG + QQQOTHER;
    //F6308 !
    //F6309          QQQStratCH4H2O = (QCH4H2O(IYR)+QCH4H2O(IYRP))/2.	! Strat H2O forcing from CH4
    const float QQQStratCH4H2O = ( G_FORCE->QCH4H2O[ IYR ] + G_FORCE->QCH4H2O[ IYRP ] ) / 2.0;
    //F6310 
    //F6311          QQQBC = ( QBC(IYR) + QBC(IYRP) )/2.
    const float QQQBC = ( G_FORCE->QBC[ IYR ] + G_FORCE->QBC[ IYRP ] ) / 2.0;
    //F6312          QQQOC = ( QOC(IYR) + QOC(IYRP) )/2.
    const float QQQOC = ( G_FORCE->QOC[ IYR ] + G_FORCE->QOC[ IYRP ] ) / 2.0;
    //F6313  
    //F6314  	     QQQTOT = QQQTOT + QQQBC + QQQOC
    QQQTOT += ( QQQBC + QQQOC );
    //F6315  	     QQQEXTRA = QQQEXTRA - (QQQBC + QQQOC)
    QQQEXTRA -= ( QQQBC + QQQOC );
    //F6316  	     
    //F6317 	  ! For consistency, make sure indices are same here as in getGHGConc for gases that overlap
    //F6318       select case (iGasNumber)
    float returnValue;
    switch ( iGasNumber ) {
            //F6319       case(0); getForcing = QQQTOT	! Total anthropogenic forcing
        case 0: returnValue = QQQTOT; break;
            //F6320       case(1); getForcing = (QCO2(IYR)+QCO2(IYRP))/2.
        case 1: returnValue = QQQCO2;  break; //CHANGE  why recalculate this?
            //F6321       case(2); getForcing = (qm(IYR)+qm(IYRP))/2. - QQQStratCH4H2O - QQCH4O3! CH4 forcing, subtract indirect components so are just reporting just CH4 forcing
        case 2: returnValue = ( G_FORCE->QM[ IYR ] + G_FORCE->QM[ IYRP ] ) / 2.0 - QQQStratCH4H2O - QQCH4O3;  break;
            //F6322       case(3); getForcing = (qn(IYR)+qn(IYRP))/2.  ! N2O forcing
        case 3: returnValue = QQQN; break; //CHANGE  why recalculate this?
            //F6323       case(4); getForcing = (QC2F6_ar(IYR)+QC2F6_ar(IYRP))/2.
        case 4: returnValue = ( G_HALOF->QC2F6_ar[ IYR ] + G_HALOF->QC2F6_ar[ IYRP ] ) / 2.0; break;
            //F6324       case(5); getForcing = (Q125_ar(IYR)+Q125_ar(IYRP))/2.
        case 5: returnValue = ( G_HALOF->Q125_ar[ IYR ] + G_HALOF->Q125_ar[ IYRP ] ) / 2.0; break;
            //F6325       case(6); getForcing = (Q134A_ar(IYR)+Q134A_ar(IYRP))/2.
        case 6: returnValue = ( G_HALOF->Q134A_ar[ IYR ] + G_HALOF->Q134A_ar[ IYRP ] ) / 2.0; break;
            //F6326       case(7); getForcing = (Q143A_ar(IYR)+Q143A_ar(IYRP))/2.
        case 7: returnValue = ( G_HALOF->Q143A_ar[ IYR ] + G_HALOF->Q143A_ar[ IYRP ] ) / 2.0; break;
            //F6327       case(8); getForcing = (Q245_ar(IYR)+Q245_ar(IYRP))/2.
        case 8: returnValue = ( G_HALOF->Q245_ar[ IYR ] + G_HALOF->Q245_ar[ IYRP ] ) / 2.0; break;
            //F6328       case(9); getForcing = (qSF6_ar(IYR)+qSF6_ar(IYRP))/2.
        case 9: returnValue = ( G_HALOF->qSF6_ar[ IYR ] + G_HALOF->qSF6_ar[ IYRP ] ) / 2.0; break;
            //F6329       case(10); getForcing = (QCF4_ar(IYR)+QCF4_ar(IYRP))/2.
        case 10: returnValue = ( G_HALOF->QCF4_ar[ IYR ] + G_HALOF->QCF4_ar[ IYRP ] ) / 2.0; break;
            //F6330       case(11); getForcing = (Q227_ar(IYR)+Q227_ar(IYRP))/2.
        case 11: returnValue = ( G_HALOF->Q227_ar[ IYR ] + G_HALOF->Q227_ar[ IYRP ] ) / 2.0; break;
            //F6331       case(12); getForcing = (QOTHER(IYR)+QOTHER(IYRP))/2.	! Other halo forcing (exogenous input)
        case 12: returnValue = QQQOTHER; break; //CHANGE  why recalculate this?
            //F6332       case(13); getForcing = QQQSO2 - DELQFOC ! Total SO2 forcing. Note QSO2 and QDIR includes FOC
        case 13: returnValue = QQQSO2 - DELQFOC; break;
            //F6333       case(14); getForcing = QQQDIR	- DELQFOC ! SO2 direct forcing only. Note QSO2 and QDIR includes FOC
        case 14: returnValue = QQQDIR - DELQFOC; break;
            //F6334       case(15); getForcing = QQQOZ ! Tropospheric Ozone forcing, including CH4 component
        case 15: returnValue = QQQOZ; break;
            //F6335       case(16); getForcing = (QCH4O3(IYR)+QCH4O3(IYRP))/2.	! Trop O3 change due to CH4
        case 16: returnValue = QQCH4O3; break; //CHANGE  why recalculate this?
            //F6336       case(17); getForcing = (QCH4H2O(IYR)+QCH4H2O(IYRP))/2.	! Strat H2O forcing from CH4
        case 17: returnValue = QQQStratCH4H2O; break; //CHANGE  why recalculate this?
            //F6337       case(18); getForcing = (QMONT(IYR)+QMONT(IYRP))/2.	! Montreal Protocol Gases forcing
        case 18: returnValue = QQQMONT; break; //CHANGE  why recalculate this?
            //F6338       case(19); getForcing = (QSTRATOZ(IYR)+QSTRATOZ(IYRP))/2.	! Stratospheric Ozone Forcing due to CFC emissions changes
        case 19: returnValue = QQQSTROZ; break; //CHANGE  why recalculate this?
            //F6339       case(20); getForcing = QQQBIO  ! MAGICC biomass burning aerosol forcing
        case 20: returnValue = QQQBIO; break;
            //F6340       case(21); getForcing = (QFOC(IYR)+QFOC(IYRP))/2. ! MAGICC internal fossil BC+OC
        case 21: returnValue = ( G_JSTART->QFOC[ IYR ] + G_JSTART->QFOC[ IYRP ] ) / 2.0; break;
            //F6341       case(22); getForcing = QQQLAND ! Land Surface Albedo forcing
        case 22: returnValue = QQQLAND; break;
            //F6342       case(23); getForcing = QQQMN	! Mineral and nitrous oxide aerosol forcing
        case 23: returnValue = QQQMN; break;
            //F6343       case(24); getForcing = QQQBC	! Custom BC forcing (total BC; land + combustion)
        case 24: returnValue = QQQBC; break;
            //F6344       case(25); getForcing = QQQOC	! Custom OC forcing (total OC; land + combustion)
        case 25: returnValue = QQQOC; break;
            //F6345       ! Note, QEXTRA is not included in total (since this is not always anthropogenic).
            //F6346       case(26); getForcing = QQQEXTRA	! User (exogenous) input forcing
        case 26: returnValue = QQQEXTRA; break;
            //F6347       ! RCPForcing. This is the older TAR definition of total forcing , exclusive of Albedo, Mineral and nitrous oxide aerosol forcing.
            //F6348       case(27); getForcing = QQQTOT - (QQQLAND + QQQMN)	! User (exogenous) input forcing
        case 27: returnValue = QQQTOT - ( QQQLAND + QQQMN ); break;
            // (Modified from fortran for GCAM reporting) Return Fossil Fuel BC/OC
        case 28: returnValue =  QQQFOCR; break;
            //F6349       case default; getForcing = -1.0
        default: returnValue = std::numeric_limits<float>::max();
                    cerr << __func__ << " undefined gas " << iGasNumber << flush;
    }
    //F6350       end select;
    //F6351       
    //F6352       RETURN 
    return returnValue;
    //F6353 	  END
    f_exit( __func__ );
}
//F6354 	  
//F6355       FUNCTION getGMTemp( inYear )
float GETGMTEMP( int inYear )
{
    f_enter( __func__ );
    assert( G_CARB != NULL && G_TANDSL != NULL && G_CONCS != NULL && G_NEWCONCS != NULL &&
           G_STOREDVALS != NULL && G_NEWPARAMS != NULL && G_BCOC != NULL && 
           G_METH1 != NULL && G_CAR != NULL && G_FORCE != NULL && G_JSTART != NULL && 
           G_QADD != NULL && G_HALOF != NULL );
    //F6356       IMPLICIT REAL*4 (a-h,o-z), Integer (I-N)
    //F6357 ! Expose subroutine gmTemp to users of this DLL
    //F6358 !DEC$ATTRIBUTES DLLEXPORT::gmTemp
    //F6359 
    //F6360       parameter (iTp=740)
    //F6361 
    //F6362       COMMON/STOREDVALS/ TEMUSER(iTp),QSO2SAVE(0:iTp+1),QDIRSAVE(0:iTp+1), &
    //F6363       KYRREF
    //F6364 
    //F6365       COMMON/TANDSL/TEQU(iTp),TGAV(iTp),TNHO(iTp), &
    //F6366       TSHO(iTp),TNHL(iTp),TSHL(iTp),TDEEP(iTp),TNHAV(iTp),TSHAV(iTp), &
    //F6367       TLAND(iTp),TOCEAN(iTp),TOCN(40),TOCNPREV(40), &
    //F6368       SIP,SGP,SAP,SLI(iTp),SLG(iTp),SLA(iTp),EX(0:iTp),SLT(iTp), &
    //F6369       QTOT(0:iTp),QGH(0:iTp),QOZ(0:iTp),QBIO(0:iTp),SLO(iTp), &
    //F6370       QSO2(0:iTp+1),QDIR(0:iTp+1),QLAND(0:iTp),QMN(0:iTp+1)
    //F6371 
    //F6372 	  REAL*4 getGMTemp
    //F6373 
    //F6374       KREF  = KYRREF-1764
    const int KREF = G_STOREDVALS->KYRREF - 1764;
    //F6375       IYR = inYear-1990+226
    const int IYR = inYear - 1990 + 226;
    //F6376       getGMTemp = TEMUSER(IYR)+TGAV(226)
    return( G_STOREDVALS->TEMUSER[ IYR ] + G_TANDSL->TGAV[ 226 ] );
    //F6377 
    //F6378       RETURN 
    //F6379 	  END
    f_exit( __func__ );
}
//F6380 
//F6381 ! Routine to pass in new values of parameters from calling program (e.g. ObjECTS) - sjs	  
//F6382     SUBROUTINE setParameterValues( index, value )
void SETPARAMETERVALUES( int index, float value )
{
    f_enter( __func__ );
    
    assert( G_CARB != NULL && G_TANDSL != NULL && G_CONCS != NULL && G_NEWCONCS != NULL &&
           G_STOREDVALS != NULL && G_NEWPARAMS != NULL && G_BCOC != NULL && 
           G_METH1 != NULL && G_CAR != NULL && G_FORCE != NULL && G_JSTART != NULL && 
           G_QADD != NULL && G_HALOF != NULL );
    
    //F6383       IMPLICIT REAL*4 (a-h,o-z), Integer (I-N)
    //F6384 ! Expose subroutine co2Conc to users of this DLL
    //F6385 !DEC$ATTRIBUTES DLLEXPORT::setParameterValues
    //F6386 
    //F6387 	  REAL*4 value
    //F6388       
    //F6389 	  REAL*4 aNewClimSens, aNewBTsoil, aNewBTGPP,aNewBTHumus,aNewDUSER,aNewFUSER, aNewSO2dir1990, aNewSO2ind1990
    //F6390       COMMON/NEWPARAMS/aNewClimSens, aNewBTsoil, DT2XUSER,aNewBTGPP,aNewBTHumus,aNewDUSER,aNewFUSER, &
    //F6391       					aNewSO2dir1990, aNewSO2ind1990
    //F6392       REAL*4 FBC1990, FOC1990, FSO2_dir1990,FSO2_ind1990, aBCUnitForcing, aOCUnitForcing, &
    //F6393              aBCBaseEmissions, aOCBaseEmissions !sjs
    //F6394       COMMON/BCOC/FBC1990, FOC1990, FSO2_dir1990,FSO2_ind1990, aBCUnitForcing, aOCUnitForcing, &
    //F6395              aBCBaseEmissions, aOCBaseEmissions !sjs
    //F6396 
    //F6397       select case (index)
    switch( index ) {
            //F6398       case(1); aNewClimSens = value
        case 1: G_NEWPARAMS->aNewClimSens = value; break;
            //F6399       case(2); aNewBTsoil = value
        case 2: G_NEWPARAMS->aNewBTsoil = value; break;
            //F6400       case(3); aNewBTHumus = value
        case 3: G_NEWPARAMS->aNewBTHumus = value; break;
            //F6401       case(4); aNewBTGPP = value
        case 4: G_NEWPARAMS->aNewBTGPP = value; break;
            //F6402       case(5); aNewDUSER = value
        case 5: G_NEWPARAMS->aNewDUSER = value; break;
            //F6403       case(6); aNewFUSER = value
        case 6: G_NEWPARAMS->aNewFUSER = value; break;
            //F6404       case(7); aNewSO2dir1990 = value
        case 7: G_NEWPARAMS->aNewSO2dir1990 = value; break;
            //F6405       case(8); aNewSO2ind1990 = value
        case 8: G_NEWPARAMS->aNewSO2ind1990 = value; break;
            //F6406       case(9); aBCUnitForcing = value
        case 9: G_BCOC->aBCUnitForcing = value; break;
            //F6407       case(10); aOCUnitForcing = value
        case 10: G_BCOC->aOCUnitForcing = value; break;
            //F6408       case default; 
            //F6409       end select;
    }
    //F6410 
    //F6411       RETURN 
    //F6412 	  END
    f_exit( __func__ );
}
//F6413 
//F6414 ! Routine to overide MAGICC parameeters with new values if these have been read-in
//F6415     SUBROUTINE overrideParameters( )
void overrideParameters( NEWPARAMS_block* NEWPARAMS, CAR_block* CAR, METH1_block* METH1, BCOC_block* BCOC )
{
    f_enter( __func__ );
    assert( G_CARB != NULL && G_TANDSL != NULL && G_CONCS != NULL && G_NEWCONCS != NULL &&
           G_STOREDVALS != NULL && G_NEWPARAMS != NULL && G_BCOC != NULL && 
           G_METH1 != NULL && G_CAR != NULL && G_FORCE != NULL && G_JSTART != NULL && 
           G_QADD != NULL && G_HALOF != NULL );
    //F6416       IMPLICIT REAL*4 (a-h,o-z), Integer (I-N)
    //F6417 
    //F6418       parameter (iTp=740)
    //F6419 
    //F6420 	  REAL*4 aNewClimSens, aNewBTsoil, aNewBTGPP,aNewBTHumus,aNewDUSER,aNewFUSER, aNewSO2dir1990, aNewSO2ind1990
    //F6421       COMMON/NEWPARAMS/aNewClimSens, aNewBTsoil, DT2XUSER,aNewBTGPP,aNewBTHumus,aNewDUSER,aNewFUSER, &
    //F6422       					aNewSO2dir1990, aNewSO2ind1990
    //F6423 
    //F6424       COMMON/CAR/EL1,EL2,EL3,TINV0(5),TINV(4,5),A(3,5),AA(4,5), &
    //F6425      BCO2(4),BTGPP,BTRESP,BTHUM,GAMP,GPP0,RESP0,QA0,U0,C0,B340(4), &
    //F6426      PHI,RG,TAUP,TAUH,TAUS,THP,THS,THH0,THS0,THPL,G1,G2,G3,FACTOR, &
    //F6427      EL21,EL32,XX1,XX2,XX3,XX4,XX5,XX6,DEE1,DEE2,DEE3,DEE4,DEE5,DEE6, &
    //F6428      FL1,FL2,FL3,XL,GAMH,GAMS,QS0,BTSOIL,FERTTYPE,TOTEM,CONVTERP, &
    //F6429      R(4),CPART(4,5),DELMASS(4,226:iTp),ABFRAC(4,226:iTp)
    //F6430 
    //F6431       COMMON /METH1/emeth(226:iTp),imeth,ch4l(225:iTp),ch4b(225:iTp), &
    //F6432      ch4h(225:iTp),ef4(226:iTp),StratH2O,TCH4(iTp),iO3feed, &
    //F6433      ednet(226:iTp+1),DUSER,FUSER,CORRUSER,CORRMHI,CORRMMID,CORRMLO
    //F6434 
    //F6435       REAL*4 FBC1990, FOC1990, FSO2_dir1990,FSO2_ind1990, aBCUnitForcing, aOCUnitForcing, &
    //F6436              aBCBaseEmissions, aOCBaseEmissions !sjs
    //F6437       COMMON/BCOC/FBC1990, FOC1990, FSO2_dir1990,FSO2_ind1990, aBCUnitForcing, aOCUnitForcing, &
    //F6438              aBCBaseEmissions, aOCBaseEmissions !sjs
    //F6439 
    //F6440       IF(aNewClimSens.GT.0)THEN
    //F6441         DT2XUSER   = aNewClimSens
    //F6442       ENDIF
    if( NEWPARAMS->aNewClimSens > 0) NEWPARAMS->DT2XUSER = NEWPARAMS->aNewClimSens;
    //F6443 
    //F6444       IF(aNewBTsoil.GT.0)THEN
    //F6445         BTSOIL   = aNewBTsoil
    //F6446       ENDIF
    if( NEWPARAMS->aNewBTsoil > 0) CAR->BTSOIL = NEWPARAMS->aNewBTsoil;
    //F6447       
    //F6448       IF(aNewBTHumus.GT.0)THEN
    //F6449         BTHUM   = aNewBTHumus
    //F6450       ENDIF
    if( NEWPARAMS->aNewBTHumus > 0) CAR->BTHUM = NEWPARAMS->aNewBTHumus;
    //F6451       
    //F6452       IF(aNewBTGPP.GT.0)THEN
    //F6453         BTGPP   = aNewBTGPP
    //F6454       ENDIF
    if( NEWPARAMS->aNewBTGPP > 0) CAR->BTGPP = NEWPARAMS->aNewBTGPP;
    //F6455       
    //F6456       IF(aNewDUSER.GT.0)THEN
    //F6457         DUSER   = aNewDUSER
    //F6458       ENDIF
    if( NEWPARAMS->aNewDUSER > 0) METH1->DUSER = NEWPARAMS->aNewDUSER;
    //F6459       
    //F6460       IF(aNewFUSER.GT.0)THEN
    //F6461         FUSER   = aNewFUSER
    //F6462       ENDIF
    if( NEWPARAMS->aNewFUSER > 0) METH1->FUSER = NEWPARAMS->aNewFUSER;
    //F6463       
    //F6464       IF(aNewSO2dir1990.LT.0)THEN
    //F6465         FSO2_dir1990   = aNewSO2dir1990
    //F6466       ENDIF
    if( NEWPARAMS->aNewSO2dir1990 < 0) BCOC->FSO2_dir1990 = NEWPARAMS->aNewSO2dir1990;
    //F6467       
    //F6468       IF(aNewSO2ind1990.LT.0)THEN
    //F6469         FSO2_ind1990   = aNewSO2ind1990
    //F6470       ENDIF
    if( NEWPARAMS->aNewSO2ind1990 < 0) BCOC->FSO2_ind1990 = NEWPARAMS->aNewSO2ind1990;
    //F6471 
    //F6472       RETURN 
    //F6473 	  END
    f_exit( __func__ );
} // overrideParameters
//F6474 	    
//F6475 ! Returns climate results forcing for a given gas
//F6476       FUNCTION getCarbonResults( iResultNumber, inYear )
float GETCARBONRESULTS( int iResultNumber, int inYear )
{
    f_enter( __func__ );
    assert( G_CARB != NULL && G_TANDSL != NULL && G_CONCS != NULL && G_NEWCONCS != NULL &&
           G_STOREDVALS != NULL && G_NEWPARAMS != NULL && G_BCOC != NULL && 
           G_METH1 != NULL && G_CAR != NULL && G_FORCE != NULL && G_JSTART != NULL && 
           G_QADD != NULL && G_HALOF != NULL );
    //F6477       IMPLICIT REAL*4 (a-h,o-z), Integer (I-N)
    //F6478 ! Expose subroutine getCarbonResults to users of this DLL
    //F6479 !DEC$ATTRIBUTES DLLEXPORT::getCarbonResults
    //F6480 
    //F6481       parameter (iTp=740)
    //F6482 
    //F6483       COMMON/CARB/CCO2(4,224:iTp),EDGROSS(4,226:iTp),EF(226:iTp+1), &
    //F6484      REGROW(4,226:iTp),PL(4,226:iTp),HL(4,226:iTp),SOIL(4,226:iTp), &
    //F6485      TTT(226:iTp),ESUM(226:iTp),ETOT(4,226:iTp),EDNET90(4), &
    //F6486      FOC(4,226:iTp),co2(0:iTp),CO2SAVE(0:iTp)
    //F6487 !
    //F6488       COMMON/CAR/EL1,EL2,EL3,TINV0(5),TINV(4,5),A(3,5),AA(4,5), &
    //F6489      BCO2(4),BTGPP,BTRESP,BTHUM,GAMP,GPP0,RESP0,QA0,U0,C0,B340(4), &
    //F6490      PHI,RG,TAUP,TAUH,TAUS,THP,THS,THH0,THS0,THPL,G1,G2,G3,FACTOR, &
    //F6491      EL21,EL32,XX1,XX2,XX3,XX4,XX5,XX6,DEE1,DEE2,DEE3,DEE4,DEE5,DEE6, &
    //F6492      FL1,FL2,FL3,XL,GAMH,GAMS,QS0,BTSOIL,FERTTYPE,TOTEM,CONVTERP, &
    //F6493      R(4),CPART(4,5),DELMASS(4,226:iTp),ABFRAC(4,226:iTp)
    //F6494  
    //F6495       COMMON /METH1/emeth(226:iTp),imeth,ch4l(225:iTp),ch4b(225:iTp), &
    //F6496      ch4h(225:iTp),ef4(226:iTp),StratH2O,TCH4(iTp),iO3feed, &
    //F6497      ednet(226:iTp+1),DUSER,FUSER,CORRUSER,CORRMHI,CORRMMID,CORRMLO
    //F6498 
    //F6499 	  REAL*4 getCarbonResults
    //F6500 	  REAL*4 NetDef, GrossDef
    float NetDef, GrossDef, TOTE;
    //F6501 
    //F6502         
    //F6503       IYR = inYear-1990+226
    const int IYR = inYear - 1990 + 226;
    //F6504 
    //F6505 !	  Branch for years > 1990
    //F6506 
    //F6507       IF ( inYear .ge. 1990 ) THEN
    if( inYear >= 1990 ) {
        //F6508       IF(IMETH.EQ.0)THEN
        if( G_METH1->IMETH == 0.0 )
            //F6509         TOTE=EF(IYR)+EDNET(IYR)
            TOTE = G_CARB->EF.getval( IYR ) + G_METH1->ednet.getval( IYR );
        //F6510       ELSE
        //F6511         TOTE=EF(IYR)+EDNET(IYR)+EMETH(IYR)
        else 
            TOTE = G_CARB->EF.getval( IYR ) + G_METH1->ednet.getval( IYR ) + G_METH1->emeth.getval( IYR );
        //F6512       ENDIF
        //F6513 	    NetDef = EDNET(IYR)
        NetDef = G_METH1->ednet.getval( IYR );
        //F6514 	    GrossDef = EDGROSS(4,IYR)
        GrossDef = G_CARB->EDGROSS.getval( 4, IYR );
    } else {
        //F6515 	  ELSE
        //F6516         TOTE = -1.0
        //F6517 		NetDef = -1.0
        //F6518 		GrossDef = -1.0
        TOTE = NetDef = GrossDef = -1.0;
        //F6519 	  ENDIF
    }
    //F6520 !
    //F6521       ECH4OX=EMETH(IYR)
    float ECH4OX = G_METH1->emeth.getval( IYR );
    //F6522       IF(IMETH.EQ.0)ECH4OX=0.0
    if( G_METH1->IMETH == 0.0 ) ECH4OX = 0.0;
    //F6523       
    //F6524       getCarbonResults = - 1.0
    float returnValue;
    //F6525       
    //F6526       select case (iResultNumber)
    switch( iResultNumber ) {
            //F6527       case(0); getCarbonResults = TOTE    ! Total emissions (fossil + netDef + Oxidation)
        case 0: returnValue = TOTE; break;
            //F6528       case(1); getCarbonResults = EF(IYR) ! Fossil Emissions as used by MAGICC
        case 1: returnValue = G_CARB->EF.getval( IYR ); break;
            //F6529       case(2); getCarbonResults = NetDef  ! Net Deforestation
        case 2: returnValue = NetDef; break;
            //F6530       case(3); getCarbonResults = GrossDef  ! Gross Deforestation
        case 3: returnValue = GrossDef; break;
            //F6531       case(4); getCarbonResults = FOC(4,IYR)  ! Ocean Flux
        case 4: returnValue = G_CARB->FOC.getval( 4, IYR ); break;
            //F6532       case(5); getCarbonResults = PL(4,IYR) ! Plant Carbon
        case 5: returnValue = G_CARB->PL.getval( 4, IYR ); break;
            //F6533       case(6); getCarbonResults = HL(4,IYR) ! Carbon in Litter
        case 6: returnValue = G_CARB->HL.getval( 4, IYR ); break;
            //F6534       case(7); getCarbonResults = SOIL(4,IYR) ! Carbon in Soils
        case 7: returnValue = G_CARB->SOIL.getval( 4, IYR ); break;
            //F6535       case(8); getCarbonResults = DELMASS(4,IYR)  ! Atmospheric Increase
        case 8: returnValue = G_CAR->DELMASS.getval( 4, IYR ); break;
            //F6536       case(9); getCarbonResults = ECH4OX  ! Oxidation Addition to Atmosphere
        case 9: returnValue = ECH4OX; break;
            //F6537       case(10); IF(inYear .ge. 1990 ) getCarbonResults = EF(IYR)+ECH4OX-(FOC(4,IYR)+DELMASS(4,IYR)) ! Net Terrestrial Uptake
        case 10: if( inYear >= 1990 ) returnValue = G_CARB->EF.getval( IYR ) + ECH4OX - (G_CARB->FOC.getval( 4, IYR ) + G_CAR->DELMASS.getval( 4, IYR )); break;
            //F6538       case default; getCarbonResults = -1.0
        default: returnValue = std::numeric_limits<float>::max();
                cerr << __func__ << " undefined result " << iResultNumber << flush;;
            //F6539       end select;
    }
    //F6540 
    //F6541       RETURN 
    return( returnValue );
    //F6542 	  END
    f_exit( __func__ );
}
//F6543 

