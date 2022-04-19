#include "advect_all_scalars.h"

void advect_all_scalars() {

  real2d dummy("dummy",nz,ncrms);
  real1d esmt_offset("esmt_offset", ncrms);
#ifdef MMF_ESMT
  YAKL_SCOPE( u_esmt  , :: u_esmt);
  YAKL_SCOPE( v_esmt  , :: v_esmt);
  real1d esmt_min("esmt_min",ncrms);
  yakl::memset(esmt_min,1.0e20);
#endif

  // advection of scalars :
  advect_scalar(t,dummy,dummy);

  // Advection of microphysics prognostics:
  for (int k=0; k<nmicro_fields; k++) {
    if ( k==index_water_vapor || (docloud && flag_precip(k)!=1) || (doprecip && flag_precip(k)==1) ) {
      advect_scalar(micro_field,k,mkadv,k,mkwle,k);  
    }
  }

  // Advection of sgs prognostics:
  if (dosgs && advect_sgs) {
    for (int k=0; k<nsgs_fields; k++) {
      advect_scalar(sgs_field,k,dummy,dummy);
    }
  }

  micro_precip_fall();

#ifdef MMF_ESMT
    // the esmt_offset simply ensures that the scalar momentum
    // tracers are positive definite during the advection calculation
    parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
       esmt_min(icrm) = min(min(u_esmt(k,j+offy_s,i+offx_s,icrm), v_esmt(k,j+offy_s,i+offx_s,icrm)), esmt_min(icrm));
    });

    parallel_for( SimpleBounds<4>(nzm,dimy_s,dimx_s,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
       esmt_offset(icrm)  = abs(esmt_min(icrm)) + 50.;
       u_esmt(k,j,i,icrm) = u_esmt(k,j,i,icrm) + esmt_offset(icrm);
       v_esmt(k,j,i,icrm) = v_esmt(k,j,i,icrm) + esmt_offset(icrm);
    });

    // advection of scalar momentum tracers
    advect_scalar(u_esmt,dummy,dummy);
    advect_scalar(v_esmt,dummy,dummy);

    parallel_for( SimpleBounds<4>(nzm,dimy_s,dimx_s,ncrms), YAKL_LAMBDA (int k, int j, int i, int icrm) {
      u_esmt(k,j,i,icrm) = u_esmt(k,j,i,icrm) - esmt_offset(icrm);
      v_esmt(k,j,i,icrm) = v_esmt(k,j,i,icrm) - esmt_offset(icrm);
   });
#endif

}
