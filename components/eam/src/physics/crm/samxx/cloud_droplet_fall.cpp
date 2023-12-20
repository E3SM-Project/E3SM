#ifdef MMF_SEDIMENTATION
#include "cloud_droplet_fall.h"
// Sedimentation of cloud water droplets based on parameterization 
// from Rogers and Yau (1989) and GCSS WG1 DYCOMS2_RF2 case
void cloud_droplet_fall() {
  YAKL_SCOPE( dtn           , :: dtn );
  YAKL_SCOPE( rho           , :: rho );
  YAKL_SCOPE( adz           , :: adz );
  YAKL_SCOPE( dz            , :: dz );
  YAKL_SCOPE( t             , :: t );
  YAKL_SCOPE( micro_field   , :: micro_field );
  YAKL_SCOPE( qcl           , :: qcl );
  // YAKL_SCOPE( qifall        , :: qifall );
  // YAKL_SCOPE( precflux      , :: precflux );

  real4d fz("fz",nz,ny,nx,ncrms);

  parallel_for( SimpleBounds<4>(nz,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    fz(k,j,i,icrm) = 0.;
  });

  // These parameters are used in MICRO_M2005 cloud optics routines.
  real constexpr rho_water = rhor;
  // real constexpr rho_snow = rhos;
  // real constexpr rho_cloud_ice = 917.;

  real coef_cl = 1.19e8 * pow( 3./(4.*3.1415*rho_water*Nc0*1.e6) , (2./3.) );

  // Compute cloud ice flux using flux limited advection scheme, as in
  // chapter 6 of Finite Volume Methods for Hyperbolic Problems by R.J.
  // LeVeque, Cambridge University Press, 2002). 
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int kc=min(nzm-1,k+1);
    int kb=max(0,k-1);
  
    // CFL number based on grid spacing interpolated to interface i,j,k-1/2
    real coef = dtn / ( 0.5 * (adz(kb,icrm)+adz(k,icrm)) * dz(icrm) );

    // Compute cloud liquid density in this cell and the ones above/below.
    // Since cloud liquid is falling, the above cell is u (upwind),
    // this cell is c (center) and the one below is d (downwind)
    real qclu = rho(kc,icrm)*qcl(kc,j,i,icrm);
    real qclc = rho(k ,icrm)*qcl(k ,j,i,icrm);
    real qcld = rho(kb,icrm)*qcl(kb,j,i,icrm); 

    // From Rogers and Yau: leading coefficient (assumed to be uniform in space)
    // is computed above and depends on (rho*qcl)^(2/3)
    // a small offset of 1.e-12 prevents issues with raising zero to a fractional power
    real vt_liq = coef_cl * pow(qclc+1.e-12,(2./3.)) * exp( 5.*pow(log(sigmag_fixed),2) );

    // Use monotonized centered difference flux limiter in flux correction
    real tmp_phi;
    if (qclc==qcld) {
      tmp_phi = 0.;
    } else {
      real tmp_theta = (qclu-qclc)/(qclc-qcld);
      tmp_phi = min( 0.5*(1.+tmp_theta), 2.);
      tmp_phi = max(0.,min(tmp_phi,2.*tmp_theta));
    }

    // Compute limited flux - since falling cloud liquid is a 1D advection 
    // problem, this flux-limited advection scheme is monotonic
    fz(k,j,i,icrm) = -vt_liq * ( qclc - 0.5*(1.-coef*vt_liq)*tmp_phi*(qclc-qcld) );
  });
  
  // Apply cloud sedimentation tendnecies to prognostic variables
  // NOTE: This only works for schemes that advect water vapor 
  // and cloud liquid mass together as a single species
  // (for P3 we will need a different approach)
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {    
    real coef = dtn/(dz(icrm)*adz(k,icrm)*rho(k,icrm));
    // The cloud liquid increment is the difference of the fluxes.
    real dqcl = coef*(fz(k,j,i,icrm)-fz(k+1,j,i,icrm));
    // Add this increment to both non-precipitating and total water.
    micro_field(0,k,j,i,icrm)  = micro_field(0,k,j,i,icrm)  + dqcl;

    // Add divergence of latent heat flux to liquid-ice static energy.
    // The latent heat flux induced by the falling cloud liquid enters
    // the LISE budget in the same way as the precipitation.  
    t(k,j,i,icrm)  = t(k,j,i,icrm) - fac_cond * dqcl;

    // // Include this effect in the total moisture budget diagnostics
    // qifall(k,icrm) = qifall(k,icrm) + dqcl; // combine cloud liquid and ice sedimentation in qifall
    // precflux(k,icrm) = precflux(k,icrm) - fz(k,j,i,icrm)*dtn/dz(icrm);

    // // Add divergence to liquid-ice static energy budget diagnostics
    // tlatqi(k,icrm) = tlatqi(k,icrm) - lat_heat; // combined cloud liquid and ice sedimentation in tlatqi
  });

}
#endif
