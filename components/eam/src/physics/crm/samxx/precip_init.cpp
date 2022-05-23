#include "precip_init.h"

void precip_init() {
  YAKL_SCOPE( rho           , :: rho );
  YAKL_SCOPE( tabs0         , :: tabs0 );
  YAKL_SCOPE( pres          , :: pres );
  YAKL_SCOPE( accrrc        , :: accrrc ); 
  YAKL_SCOPE( accrsi        , :: accrsi );
  YAKL_SCOPE( accrsc        , :: accrsc );
  YAKL_SCOPE( coefice       , :: coefice );
  YAKL_SCOPE( evaps1        , :: evaps1 );
  YAKL_SCOPE( evaps2        , :: evaps2 );
  YAKL_SCOPE( accrgi        , :: accrgi );
  YAKL_SCOPE( accrgc        , :: accrgc );
  YAKL_SCOPE( evapg1        , :: evapg1 );
  YAKL_SCOPE( evapg2        , :: evapg2 );
  YAKL_SCOPE( evapr1        , :: evapr1 );
  YAKL_SCOPE( evapr2        , :: evapr2 );
  YAKL_SCOPE( gams1         , :: gams1 );
  YAKL_SCOPE( gams2         , :: gams2 );
  YAKL_SCOPE( gamr1         , :: gamr1 );
  YAKL_SCOPE( gamr2         , :: gamr2 );
  YAKL_SCOPE( gamg1         , :: gamg1 );
  YAKL_SCOPE( gamg2         , :: gamg2 );
  YAKL_SCOPE( ncrms         , :: ncrms );

  gam3  = gammafff(3.0             );
  gamr1 = gammafff(3.0+b_rain      );
  gamr2 = gammafff((5.0+b_rain)/2.0);
  gamr3 = gammafff(4.0+b_rain      );
  gams1 = gammafff(3.0+b_snow      );
  gams2 = gammafff((5.0+b_snow)/2.0);
  gams3 = gammafff(4.0+b_snow      );
  gamg1 = gammafff(3.0+b_grau      );
  gamg2 = gammafff((5.0+b_grau)/2.0);
  gamg3 = gammafff(4.0+b_grau      );

  if(round(gam3) != 2) {
    std::cout << "cannot compute gamma-function in precip_init." << std::endl;
    std::exit(-1);
  }

  // for (int k=0; k<nzm; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    real pratio = sqrt(1.29 / rho(k,icrm));
    real rrr1=393.0/(tabs0(k,icrm)+120.0)*pow((tabs0(k,icrm)/273.0),1.5);
    real rrr2=pow((tabs0(k,icrm)/273.0),1.94)*(1000.0/pres(k,icrm));
    real estw = 100.0*esatw_crm(tabs0(k,icrm));
    real esti = 100.0*esati_crm(tabs0(k,icrm));
    
    // accretion by snow:
    real coef1 = 0.25 * pi * nzeros * a_snow * gams1 * pratio/pow((pi * rhos * nzeros/rho(k,icrm) ) , ((3.0+b_snow)/4.0));
    real coef2 = exp(0.025*(tabs0(k,icrm) - 273.15));
    accrsi(k,icrm) =  coef1 * coef2 * esicoef;
    accrsc(k,icrm) =  coef1 * esccoef;
    coefice(k,icrm) =  coef2;

    // evaporation of snow:
    coef1  =(lsub/(tabs0(k,icrm)*rv)-1.0)*lsub/(therco*rrr1*tabs0(k,icrm));
    coef2  = rv*tabs0(k,icrm)/(diffelq*rrr2*esti);
    evaps1(k,icrm)  =  0.65*4.0*nzeros/sqrt(pi*rhos*nzeros)/(coef1+coef2)/sqrt(rho(k,icrm));
    evaps2(k,icrm)  =  0.49*4.0*nzeros*gams2*sqrt(a_snow/(muelq*rrr1))/pow((pi*rhos*nzeros) , ((5.0+b_snow)/8.0)) / 
                       (coef1+coef2) * pow(rho(k,icrm) , ((1.0+b_snow)/8.0))*sqrt(pratio);
      
    // accretion by graupel:
    coef1 = 0.25*pi*nzerog*a_grau*gamg1*pratio/pow((pi*rhog*nzerog/rho(k,icrm)) , ((3.0+b_grau)/4.0));
    coef2 = exp(0.025*(tabs0(k,icrm) - 273.15));
    accrgi(k,icrm) =  coef1 * coef2 * egicoef;
    accrgc(k,icrm) =  coef1 * egccoef;

    // evaporation of graupel:
    coef1  =(lsub/(tabs0(k,icrm)*rv)-1.0)*lsub/(therco*rrr1*tabs0(k,icrm));
    coef2  = rv*tabs0(k,icrm)/(diffelq*rrr2*esti);
    evapg1(k,icrm)  = 0.65*4.0*nzerog/sqrt(pi*rhog*nzerog)/(coef1+coef2)/sqrt(rho(k,icrm));
    evapg2(k,icrm)  = 0.49*4.0*nzerog*gamg2*sqrt(a_grau/(muelq*rrr1))/pow((pi * rhog * nzerog) , ((5.0+b_grau)/8.0)) / 
                      (coef1+coef2) * pow(rho(k,icrm) , ((1.0+b_grau)/8.0))*sqrt(pratio);

    // accretion by rain:
    accrrc(k,icrm)=  0.25 * pi * nzeror * a_rain * gamr1 * pratio/pow((pi * rhor * nzeror / rho(k,icrm)) , ((3+b_rain)/4.))* erccoef;

    // evaporation of rain:
    coef1  =(lcond/(tabs0(k,icrm)*rv)-1.0)*lcond/(therco*rrr1*tabs0(k,icrm));
    coef2  = rv*tabs0(k,icrm)/(diffelq * rrr2 * estw);
    evapr1(k,icrm)  =  0.78 * 2.0 * pi * nzeror / 
    sqrt(pi * rhor * nzeror) / (coef1+coef2) / sqrt(rho(k,icrm));
    evapr2(k,icrm)  =  0.31 * 2.0 * pi  * nzeror * gamr2 * 0.89 * sqrt(a_rain/(muelq*rrr1))/
                        pow((pi * rhor * nzeror) , ((5.0+b_rain)/8.0)) / 
                        (coef1+coef2) * pow(rho(k,icrm) , ((1.0+b_rain)/8.0))*sqrt(pratio);
  });
}


