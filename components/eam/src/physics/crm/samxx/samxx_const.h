
#pragma once

#include "YAKL.h"
#include "Array.h"
#include <iomanip>

using yakl::c::Bounds;
using yakl::c::SimpleBounds;
using yakl::min;
using yakl::max;
using yakl::abs;
using yakl::c::parallel_for;
using yakl::SArray;
using yakl::ScalarLiveOut;

template <class T>
void DEBUG(T var) {
  std::cout << var.myname << ": " << std::setprecision(16) << std::scientific << yakl::intrinsics::sum(var) << std::endl;
}
template <class T>
void DEBUG_LOG(T var) {
  std::cout << var.myname << ": " << std::setprecision(16) << std::scientific << yakl::intrinsics::count(var) << std::endl;
}
template <class T>
void DEBUG_SCALAR(T var) {
  std::cout << std::setprecision(16) << std::scientific << var << std::endl;
}

typedef double real;

int  constexpr crm_nx     = CRM_NX;
int  constexpr crm_ny     = CRM_NY;
int  constexpr crm_nz     = CRM_NZ;
int  constexpr crm_nx_rad = CRM_NX_RAD;
int  constexpr crm_ny_rad = CRM_NY_RAD;
real constexpr crm_dx     = CRM_DX;
real constexpr crm_dy     = crm_dx;
real constexpr crm_dt     = CRM_DT;

int  constexpr YES3D = YES3DVAL  ; // Domain dimensionality: 1 - 3D, 0 - 2D
int  constexpr nx_gl = crm_nx    ; // Number of grid points in X
int  constexpr ny_gl = crm_ny    ; // Number of grid points in Y
int  constexpr nz_gl = crm_nz    ; // Number of pressure (scalar) levels
int  constexpr nsubdomains_x  = 1; // No of subdomains in x
int  constexpr nsubdomains_y  = 1; // No of subdomains in y
int  constexpr navgmom_x = -1    ; 
int  constexpr navgmom_y = -1    ; 

int  constexpr nx = nx_gl/nsubdomains_x;
int  constexpr ny = ny_gl/nsubdomains_y;
int  constexpr nz = nz_gl+1            ;                     //  note that nz_gl = crm_nz
int  constexpr nzm = nz-1              ;                     //  note that nzm   = crm_nz
int  constexpr nsubdomains = nsubdomains_x * nsubdomains_y;
int  constexpr RUN3D = ny_gl > 1;
int  constexpr RUN2D = ! RUN3D;
int  constexpr nxp1 = nx + 1;
int  constexpr nyp1 = ny + 1 * YES3D;
int  constexpr nxp2 = nx + 2;
int  constexpr nyp2 = ny + 2 * YES3D;
int  constexpr nxp3 = nx + 3;
int  constexpr nyp3 = ny + 3 * YES3D;
int  constexpr nxp4 = nx + 4;
int  constexpr nyp4 = ny + 4 * YES3D;
int  constexpr dimx1_u = -1          ;
int  constexpr dimx2_u = nxp3        ;
int  constexpr dimy1_u = 1-(2)*YES3D ;
int  constexpr dimy2_u = nyp2        ;
int  constexpr dimx1_v = -1          ;
int  constexpr dimx2_v = nxp2        ;
int  constexpr dimy1_v = 1-2*YES3D   ;
int  constexpr dimy2_v = nyp3        ;
int  constexpr dimx1_w = -1          ;
int  constexpr dimx2_w = nxp2        ;
int  constexpr dimy1_w = 1-(2)*YES3D ;
int  constexpr dimy2_w = nyp2        ;
int  constexpr dimx1_s = -2          ;
int  constexpr dimx2_s = nxp3        ;
int  constexpr dimy1_s = 1-(3)*YES3D ;
int  constexpr dimy2_s = nyp3        ;
int  constexpr dimx1_d = 0           ;
int  constexpr dimx2_d = nxp1        ;
int  constexpr dimy1_d = 1-YES3D     ;
int  constexpr dimy2_d = nyp1        ;
int  constexpr dimx1_sstxy = 0       ;
int  constexpr dimx2_sstxy = nx      ;
int  constexpr dimy1_sstxy = 1-YES3D ;
int  constexpr dimy2_sstxy = ny      ;
int  constexpr ncols = nx*ny;
int  constexpr nadams = 3;

int  constexpr dimx_u     = dimx2_u - dimx1_u   + 1;
int  constexpr dimx_v     = dimx2_v - dimx1_v   + 1;
int  constexpr dimx_w     = dimx2_w - dimx1_w   + 1;
int  constexpr dimx_s     = dimx2_s - dimx1_s   + 1;
int  constexpr dimx_d     = dimx2_d - dimx1_d   + 1;
int  constexpr dimx_sstxy = dimx2_sstxy - dimx1_sstxy + 1;

int  constexpr dimy_u     = dimy2_u - dimy1_u   + 1;
int  constexpr dimy_v     = dimy2_v - dimy1_v   + 1;
int  constexpr dimy_w     = dimy2_w - dimy1_w   + 1;
int  constexpr dimy_s     = dimy2_s - dimy1_s   + 1;
int  constexpr dimy_d     = dimy2_d - dimy1_d   + 1;
int  constexpr dimy_p     = ny      - (1-YES3D) + 1;
int  constexpr dimy_tk2   = nyp1    - (1-YES3D) + 1;
int  constexpr dimy_sstxy = ny      - (1-YES3D) + 1;
int  constexpr dimy_fcory = ny      - 0         + 1;

int  constexpr offx_u     = 1 - dimx1_u  ;
int  constexpr offx_v     = 1 - dimx1_v  ;
int  constexpr offx_w     = 1 - dimx1_w  ;
int  constexpr offx_s     = 1 - dimx1_s  ;
int  constexpr offx_d     = 1 - dimx1_d  ;
int  constexpr offx_sstxy = 1 - dimx1_sstxy;
int  constexpr offx_p     = 1 - 0;

int  constexpr offy_u     = 1 - dimy1_u  ;
int  constexpr offy_v     = 1 - dimy1_v  ;
int  constexpr offy_w     = 1 - dimy1_w  ;
int  constexpr offy_s     = 1 - dimy1_s  ;
int  constexpr offy_d     = 1 - dimy1_d  ;
int  constexpr offy_p     = 1 - (1-YES3D);
int  constexpr offy_tk2   = 1 - (1-YES3D);
int  constexpr offy_sstxy = 1 - (1-YES3D);
int  constexpr offy_fcory = 1 - 0        ;



real constexpr SHR_CONST_PI      = 3.14159265358979323846;  // pi
real constexpr SHR_CONST_CDAY    = 86400.0;      // sec in calendar day ~ sec
real constexpr SHR_CONST_SDAY    = 86164.0;      // sec in siderial day ~ sec
real constexpr SHR_CONST_OMEGA   = 2.0*SHR_CONST_PI/SHR_CONST_SDAY; // earth rot ~ rad/sec
real constexpr SHR_CONST_REARTH  = 6.37122e6  ;  // radius of earth ~ m
real constexpr SHR_CONST_G       = 9.80616    ;  // acceleration of gravity ~ m/s^2
real constexpr SHR_CONST_STEBOL  = 5.67e-8    ;  // Stefan-Boltzmann constant ~ W/m^2/K^4
real constexpr SHR_CONST_BOLTZ   = 1.38065e-23;  // Boltzmann's constant ~ J/K/molecule
real constexpr SHR_CONST_AVOGAD  = 6.02214e26 ;  // Avogadro's number ~ molecules/kmole
real constexpr SHR_CONST_RGAS    = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ;       // Universal gas constant ~ J/K/kmole
real constexpr SHR_CONST_MWDAIR  = 28.966;                       // molecular weight dry air ~ kg/kmole
real constexpr SHR_CONST_MWWV    = 18.016;                       // molecular weight water vapor
real constexpr SHR_CONST_RDAIR   = SHR_CONST_RGAS/SHR_CONST_MWDAIR       ; // Dry air gas constant     ~ J/K/kg
real constexpr SHR_CONST_RWV     = SHR_CONST_RGAS/SHR_CONST_MWWV         ; // Water vapor gas constant ~ J/K/kg
real constexpr SHR_CONST_ZVIR    = (SHR_CONST_RWV/SHR_CONST_RDAIR)-1.0; // RWV/RDAIR - 1.0
real constexpr SHR_CONST_KARMAN  = 0.4                  ;       // Von Karman constant
real constexpr SHR_CONST_PSTD    = 101325.0             ;       // standard pressure ~ pascals
real constexpr SHR_CONST_PDB     = 0.0112372            ;       // ratio of 13C/12C in Pee Dee Belemnite (C isotope standard)
real constexpr SHR_CONST_TKTRIP  = 273.16               ;       // triple point of fresh water        ~ K
real constexpr SHR_CONST_TKFRZ   = 273.15               ;       // freezing T of fresh water          ~ K
real constexpr SHR_CONST_TKFRZSW = SHR_CONST_TKFRZ - 1.8;       // freezing T of salt water  ~ K
real constexpr SHR_CONST_ZSRFLYR = 3.0                  ;       // ocn surf layer depth for diurnal SST cal ~ m
real constexpr SHR_CONST_RHODAIR = SHR_CONST_PSTD/(SHR_CONST_RDAIR*SHR_CONST_TKFRZ);
real constexpr SHR_CONST_RHOFW   = 1.000e3  ;                   // density of fresh water     ~ kg/m^3
real constexpr SHR_CONST_RHOSW   = 1.026e3  ;                   // density of sea water       ~ kg/m^3
real constexpr SHR_CONST_RHOICE  = 0.917e3  ;                   // density of ice             ~ kg/m^3
real constexpr SHR_CONST_CPDAIR  = 1.00464e3;                   // specific heat of dry air   ~ J/kg/K
real constexpr SHR_CONST_CPWV    = 1.810e3  ;                   // specific heat of water vap ~ J/kg/K
real constexpr SHR_CONST_CPVIR   = (SHR_CONST_CPWV/SHR_CONST_CPDAIR)-1.0; // CPWV/CPDAIR - 1.0
real constexpr SHR_CONST_CPFW    = 4.188e3  ;                   // specific heat of fresh h2o ~ J/kg/K
real constexpr SHR_CONST_CPSW    = 3.996e3  ;                   // specific heat of sea h2o   ~ J/kg/K
real constexpr SHR_CONST_CPICE   = 2.11727e3;                   // specific heat of fresh ice ~ J/kg/K
real constexpr SHR_CONST_LATICE  = 3.337e5  ;                   // latent heat of fusion      ~ J/kg
real constexpr SHR_CONST_LATVAP  = 2.501e6  ;                   // latent heat of evaporation ~ J/kg
real constexpr SHR_CONST_LATSUB  = SHR_CONST_LATICE + SHR_CONST_LATVAP;
real constexpr SHR_CONST_CONDICE = 2.1;                         // thermal conductivity of ice ~ W/m/K
real constexpr SHR_CONST_KAPPA_LAND_ICE = SHR_CONST_CONDICE / (SHR_CONST_RHOICE*SHR_CONST_CPICE);
real constexpr SHR_CONST_TF0    = 6.22e-2                     ; // The freezing temperature at zero pressure in
real constexpr SHR_CONST_DTF_DP = -7.43e-8                    ; // The coefficient for the term proportional to the (limited)
real constexpr SHR_CONST_DTF_DS = -5.63e-2                    ; //The coefficient for the term proportional to salinity in
real constexpr SHR_CONST_DTF_DPDS = -1.74e-10                 ; // The coefficient for the term proportional to salinity times
real constexpr SHR_CONST_OCN_REF_SAL = 34.7                   ; // ocn ref salinity (psu)
real constexpr SHR_CONST_ICE_REF_SAL =  4.0                   ; // ice ref salinity (psu)
real constexpr SHR_CONST_SPVAL        = 1.0e30                ; // special missing value
real constexpr SHR_CONST_SPVAL_TOLMIN = 0.99 * SHR_CONST_SPVAL; // min spval tolerance
real constexpr SHR_CONST_SPVAL_TOLMAX = 1.01 * SHR_CONST_SPVAL; // max spval tolerance
real constexpr SHR_CONST_SPVAL_AERODEP= 1.e29                 ; // special aerosol deposition
real constexpr SHR_CONST_VSMOW_18O   = 2005.2e-6              ; // 18O/16O in VMSOW
real constexpr SHR_CONST_VSMOW_17O   = 379.e-6                ; // 18O/16O in VMSOW
real constexpr SHR_CONST_VSMOW_16O   = 0.997628               ; // 16O/Tot in VMSOW
real constexpr SHR_CONST_VSMOW_D   = 155.76e-6                ; // 2H/1H in VMSOW
real constexpr SHR_CONST_VSMOW_T   = 1.85e-6                  ; // 3H/1H in VMSOW
real constexpr SHR_CONST_VSMOW_H   = 0.99984426               ; // 1H/Tot in VMSOW
real constexpr SHR_CONST_RSTD_H2ODEV   = 1.0                  ; // Rstd Dev Use

real constexpr cp    = SHR_CONST_CPDAIR ;
real constexpr ggr   = SHR_CONST_G      ;
real constexpr lcond = SHR_CONST_LATVAP ;
real constexpr lfus  = SHR_CONST_LATICE ;
real constexpr lsub  = lcond + lfus     ;
real constexpr rgas  = SHR_CONST_RDAIR  ;
real constexpr rv    = SHR_CONST_RGAS/SHR_CONST_MWWV ;


real constexpr diffelq = 2.21e-05 ;    // Diffusivity of water vapor, m2/s
real constexpr therco  = 2.40e-02 ;    // Thermal conductivity of air, J/m/s/K
real constexpr muelq   = 1.717e-05;    // Dynamic viscosity of air
real constexpr fac_cond = lcond/cp;
real constexpr fac_fus  = lfus/cp ; 
real constexpr fac_sub  = lsub/cp ; 
real constexpr epsv = 0.61e0;     // = (1-eps)/eps, where eps= Rv/Ra
real constexpr pi = 3.141592653589793 ;  // sine, cosine, cosine, sine, 3.14159 !


int  constexpr nsgs_fields = 1;         // total number of prognostic sgs vars
int  constexpr nsgs_fields_diag = 2;    // total number of diagnostic sgs vars
bool constexpr do_sgsdiag_bound = true; // exchange boundaries for diagnostics fields
int  constexpr nmicro_fields = 2;
int  constexpr index_water_vapor = 0;
int  constexpr index_cloud_ice = 0;

real constexpr rhor = 1000.; // Density of water, kg/m3
real constexpr rhos = 100.;  // Density of snow, kg/m3
real constexpr rhog = 400.;  // Density of graupel, kg/m3
real constexpr tbgmin = 253.16;    // Minimum temperature for cloud water., K
real constexpr tbgmax = 273.16;    // Maximum temperature for cloud ice, K
real constexpr tprmin = 268.16;    // Minimum temperature for rain, K
real constexpr tprmax = 283.16;    // Maximum temperature for snow+graupel, K
real constexpr tgrmin = 223.16;    // Minimum temperature for snow, K
real constexpr tgrmax = 283.16;    // Maximum temperature for graupel, K
real constexpr a_rain = 842.; // Coeff.for rain term vel
real constexpr b_rain = 0.8;  // Fall speed exponent for rain
real constexpr a_snow = 4.84; // Coeff.for snow term vel
real constexpr b_snow = 0.25; // Fall speed exponent for snow
real constexpr a_grau = 94.5; // Lin (1983) (rhog=400)
real constexpr b_grau = 0.5;  // Fall speed exponent for graupel
real constexpr erccoef = 1.0;   // Rain/Cloud water collection efficiency
real constexpr esccoef = 1.0;   // Snow/Cloud water collection efficiency
real constexpr esicoef = 0.1;   // Snow/cloud ice collection efficiency
real constexpr egccoef = 1.0;   // Graupel/Cloud water collection efficiency
real constexpr egicoef = 0.1;   // Graupel/Cloud ice collection efficiency
real constexpr nzeror = 8.e6;   // Intercept coeff. for rain
real constexpr nzeros = 3.e6;   // Intersept coeff. for snow
real constexpr nzerog = 4.e6;   // Intersept coeff. for graupel
real constexpr qp_threshold = 1.e-8; // minimal rain/snow water content

real constexpr qcw0 = 1.e-3;
real constexpr qci0 = 1.e-4;
real constexpr alphaelq = 1.e-3;
real constexpr betaelq = 1.e-3;

real constexpr crm_accel_coef = 1.0/( (real) nx * (real) ny );

int constexpr plev = PLEV;

real constexpr UMAX = 0.5*crm_dx/crm_dt;
real constexpr wmin = 2.0;
real constexpr cwp_threshold = .001;
int  constexpr perturb_seed_scale = 1000;

typedef yakl::Array<real,1,yakl::memDevice,yakl::styleC> umgReal1d;
typedef yakl::Array<real,2,yakl::memDevice,yakl::styleC> umgReal2d;
typedef yakl::Array<real,3,yakl::memDevice,yakl::styleC> umgReal3d;
typedef yakl::Array<real,4,yakl::memDevice,yakl::styleC> umgReal4d;
typedef yakl::Array<real,5,yakl::memDevice,yakl::styleC> umgReal5d;
typedef yakl::Array<real,6,yakl::memDevice,yakl::styleC> umgReal6d;
typedef yakl::Array<real,7,yakl::memDevice,yakl::styleC> umgReal7d;

typedef yakl::Array<int,1,yakl::memDevice,yakl::styleC> umgInt1d;
typedef yakl::Array<int,2,yakl::memDevice,yakl::styleC> umgInt2d;
typedef yakl::Array<int,3,yakl::memDevice,yakl::styleC> umgInt3d;


typedef yakl::Array<real,1,yakl::memDevice,yakl::styleC> real1d;
typedef yakl::Array<real,2,yakl::memDevice,yakl::styleC> real2d;
typedef yakl::Array<real,3,yakl::memDevice,yakl::styleC> real3d;
typedef yakl::Array<real,4,yakl::memDevice,yakl::styleC> real4d;
typedef yakl::Array<real,5,yakl::memDevice,yakl::styleC> real5d;
typedef yakl::Array<real,6,yakl::memDevice,yakl::styleC> real6d;
typedef yakl::Array<real,7,yakl::memDevice,yakl::styleC> real7d;

typedef yakl::Array<int,1,yakl::memDevice,yakl::styleC> int1d;
typedef yakl::Array<int,2,yakl::memDevice,yakl::styleC> int2d;
typedef yakl::Array<int,3,yakl::memDevice,yakl::styleC> int3d;

typedef yakl::Array<bool,1,yakl::memDevice,yakl::styleC> bool1d;
typedef yakl::Array<bool,2,yakl::memDevice,yakl::styleC> bool2d;
typedef yakl::Array<bool,3,yakl::memDevice,yakl::styleC> bool3d;

typedef yakl::Array<bool,1,yakl::memHost,yakl::styleC> boolHost1d;
typedef yakl::Array<bool,2,yakl::memHost,yakl::styleC> boolHost2d;
typedef yakl::Array<bool,3,yakl::memHost,yakl::styleC> boolHost3d;

typedef yakl::Array<real,1,yakl::memHost,yakl::styleC> realHost1d;
typedef yakl::Array<real,2,yakl::memHost,yakl::styleC> realHost2d;
typedef yakl::Array<real,3,yakl::memHost,yakl::styleC> realHost3d;
typedef yakl::Array<real,4,yakl::memHost,yakl::styleC> realHost4d;
typedef yakl::Array<real,5,yakl::memHost,yakl::styleC> realHost5d;
typedef yakl::Array<real,6,yakl::memHost,yakl::styleC> realHost6d;
typedef yakl::Array<real,7,yakl::memHost,yakl::styleC> realHost7d;

typedef yakl::Array<int,1,yakl::memHost,yakl::styleC> intHost1d;
typedef yakl::Array<int,2,yakl::memHost,yakl::styleC> intHost2d;
typedef yakl::Array<int,3,yakl::memHost,yakl::styleC> intHost3d;
typedef yakl::Array<int,4,yakl::memHost,yakl::styleC> intHost4d;
typedef yakl::Array<int,5,yakl::memHost,yakl::styleC> intHost5d;
typedef yakl::Array<int,6,yakl::memHost,yakl::styleC> intHost6d;
typedef yakl::Array<int,7,yakl::memHost,yakl::styleC> intHost7d;

