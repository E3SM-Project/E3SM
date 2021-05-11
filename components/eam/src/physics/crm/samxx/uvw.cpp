#include "uvw.h"

void uvw() {
  auto &u             = :: u;
  auto &v             = :: v;
  auto &w             = :: w;
  auto &dudt          = :: dudt;
  auto &dvdt          = :: dvdt;
  auto &dwdt          = :: dwdt;
  auto &nc            = :: nc;
  auto &ncrms         = :: ncrms;

  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    u(k,j+offy_u,i+offx_u,icrm) = dudt(nc-1,k,j,i,icrm);
    v(k,j+offy_v,i+offx_v,icrm) = dvdt(nc-1,k,j,i,icrm);
    w(k,j+offy_w,i+offx_w,icrm) = dwdt(nc-1,k,j,i,icrm);
  });

}
