#include "uvw.h"

void uvw() {
  YAKL_SCOPE( u             , :: u );
  YAKL_SCOPE( v             , :: v );
  YAKL_SCOPE( w             , :: w );
  YAKL_SCOPE( dudt          , :: dudt );
  YAKL_SCOPE( dvdt          , :: dvdt );
  YAKL_SCOPE( dwdt          , :: dwdt );
  YAKL_SCOPE( nc            , :: nc );
  YAKL_SCOPE( ncrms         , :: ncrms );

  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    u(k,j+offy_u,i+offx_u,icrm) = dudt(nc-1,k,j,i,icrm);
    v(k,j+offy_v,i+offx_v,icrm) = dvdt(nc-1,k,j,i,icrm);
    w(k,j+offy_w,i+offx_w,icrm) = dwdt(nc-1,k,j,i,icrm);
  });

}
