#include "advect_all_scalars.h"

void advect_all_scalars() {

  real2d dummy("dummy",nz,ncrms);

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
}
