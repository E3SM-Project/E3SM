
#include "shear_prod2D.h"

void shear_prod2D(real4d &def2) {
  YAKL_SCOPE( dx             , :: dx );
  YAKL_SCOPE( dz             , :: dz );
  YAKL_SCOPE( adz            , :: adz );
  YAKL_SCOPE( adzw           , :: adzw );
  YAKL_SCOPE( u              , :: u );
  YAKL_SCOPE( v              , :: v );
  YAKL_SCOPE( w              , :: w );
  YAKL_SCOPE( u0             , :: u0 );
  YAKL_SCOPE( v0             , :: v0 );
  YAKL_SCOPE( ncrms          , :: ncrms );

  // for (int k=0; k<nzm; k++) {
  //    for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(nzm,nx,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
    real rdx0 = 1.0/dx;
    int j = 0;
    int kb, kc, ib, ic;
    real rdz, rdzw_up, rdzw_dn, rdx, rdx_up, rdx_dn;

    if (k>=1 && k <= nzm-2) {

      kb = k-1;
      kc = k+1;
      rdz = 1.0/(dz(icrm)*adz(k,icrm));
      rdzw_up = 1.0/(dz(icrm)*adzw(kc,icrm));
      rdzw_dn = 1.0/(dz(icrm)*adzw(k,icrm));
      rdx = rdx0*sqrt(dx*rdz);
      rdx_up = rdx0 *sqrt(dx*rdzw_up);
      rdx_dn = rdx0 *sqrt(dx*rdzw_dn);
      ib = i-1;
      ic = i+1;
      real tmp1 = ((u(k,j+offy_u,ic+offx_u,icrm)-u(k,j+offy_u,i+offx_u,icrm))*rdx);
      real tmp2 = ((w(kc,j+offy_w,i+offx_w,icrm)-w(k,j+offy_w,i+offx_w,icrm))*rdz);
      real tmp3 = (v(k,j+offy_v,ic+offx_v,icrm)-v(k,j+offy_v,i+offx_v,icrm))*rdx;
      real tmp4 = (v(k,j+offy_v,i+offx_v,icrm)-v(k,j+offy_v,ib+offx_v,icrm))*rdx;
      real tmp5 = ( (u(kc,j+offy_u,ic+offx_u,icrm)-u0(kc,icrm)-u(k,j+offy_u,ic+offx_u,icrm)+u0(k,icrm))*rdzw_up +
                    (w(kc,j+offy_w,ic+offx_w,icrm)-w(kc,j+offy_w,i+offx_w,icrm))*rdx_up);
      real tmp6 = ( (u(kc,j+offy_u,i+offx_u,icrm)-u0(kc,icrm)-u(k,j+offy_u,i+offx_u,icrm)+u0(k,icrm))*rdzw_up +
                    (w(kc,j+offy_w,i+offx_w,icrm)-w(kc,j+offy_w,ib+offx_w,icrm))*rdx_up);
      real tmp7 = ( (u(k,j+offy_u,ic+offx_u,icrm)-u0(k,icrm)-u(kb,j+offy_u,ic+offx_u,icrm)+u0(kb,icrm))*rdzw_dn +
                    (w(k,j+offy_w,ic+offx_w,icrm)-w(k,j+offy_w,i+offx_w,icrm))*rdx_dn);
      real tmp8 = ( (u(k,j+offy_u,i+offx_u,icrm)-u0(k,icrm)-u(kb,j+offy_u,i+offx_u,icrm)+u0(kb,icrm))*rdzw_dn +
                    (w(k,j+offy_w,i+offx_w,icrm)-w(k,j+offy_w,ib+offx_w,icrm))*rdx_dn);
      real tmp9  = ( (v(kc,j+offy_v,i+offx_v,icrm)-v0(kc,icrm)-v(k,j+offy_v,i+offx_v,icrm)+v0(k,icrm))*rdzw_up);
      real tmp10 = ( (v(k,j+offy_v,i+offx_v,icrm)-v0(k,icrm)-v(kb,j+offy_v,i+offx_v,icrm)+v0(kb,icrm))*rdzw_dn);
      def2(k,j,i,icrm) = 2.0 * ( tmp1*tmp1 + tmp2*tmp2 ) + 0.5 * ( tmp3*tmp3 + tmp4*tmp4 + tmp5*tmp5 + tmp6 *tmp6  +
                                                                   tmp7*tmp7 + tmp8*tmp8 + tmp9*tmp9 + tmp10*tmp10 );

    } else if (k==0) {

      kc = k+1;
      rdz = 1.0/(dz(icrm)*adz(k,icrm));
      rdzw_up = 1.0/(dz(icrm)*adzw(kc,icrm));
      rdx = rdx0*sqrt(dx*rdz);
      rdx_up = rdx0 *sqrt(dx*rdzw_up);
      ib = i-1;
      ic = i+1;
      real tmp1 = ((u(k,j+offy_u,ic+offx_u,icrm)-u(k,j+offy_u,i+offx_u,icrm))*rdx);
      real tmp2 = ((w(kc,j+offy_w,i+offx_w,icrm)-w(k,j+offy_w,i+offx_w,icrm))*rdz);
      real tmp3 = ( (v(k,j+offy_v,ic+offx_v,icrm)-v(k,j+offy_v,i+offx_v,icrm))*rdx);
      real tmp4 = ( (v(k,j+offy_v,i+offx_v,icrm)-v(k,j+offy_v,ib+offx_v,icrm))*rdx);
      real tmp5 = ( (v(kc,j+offy_v,i+offx_v,icrm)-v0(kc,icrm)-v(k,j+offy_v,i+offx_v,icrm)+v0(k,icrm))*rdzw_up);
      real tmp6 = ( (u(kc,j+offy_u,ic+offx_u,icrm)-u0(kc,icrm)-u(k,j+offy_u,ic+offx_u,icrm)+u0(k,icrm))*rdzw_up +
                    (w(kc,j+offy_w,ic+offx_w,icrm)-w(kc,j+offy_w,i+offx_w,icrm))*rdx_up);
      real tmp7 = ( (u(kc,j+offy_u,i+offx_u,icrm)-u0(kc,icrm)-u(k,j+offy_u,i+offx_u,icrm)+u0(k,icrm))*rdzw_up +
                    (w(kc,j+offy_w,i+offx_w,icrm)-w(kc,j+offy_w,ib+offx_w,icrm))*rdx_up);
      def2(k,j,i,icrm) = 2.0* ( tmp1*tmp1 + tmp2*tmp2 ) + 0.5 * ( tmp3*tmp3 + tmp4*tmp4 ) + tmp5*tmp5 +
                                                          0.5 * ( tmp6*tmp6 + tmp7*tmp7 );

    } else if (k==nzm-1) {

      kc = k+1;
      kb = k-1;
      rdz = 1.0/(dz(icrm)*adz(k,icrm));
      rdzw_dn = 1.0/(dz(icrm)*adzw(k,icrm));
      rdx = rdx0*sqrt(dx*rdz);
      rdx_dn = rdx0 *sqrt(dx*rdzw_dn);
      ib = i-1;
      ic = i+1;
      real tmp1 = ( (u(k,j+offy_u,ic+offx_u,icrm)-u(k,j+offy_u,i+offx_u,icrm))*rdx);
      real tmp2 = ( (w(kc,j+offy_w,i+offx_w,icrm)-w(k,j+offy_w,i+offx_w,icrm))*rdz);
      real tmp3 = ( (v(k,j+offy_v,ic+offx_v,icrm)-v(k,j+offy_v,i+offx_v,icrm))*rdx);
      real tmp4 = ( (v(k,j+offy_v,i+offx_v,icrm)-v(k,j+offy_v,ib+offx_v,icrm))*rdx);
      real tmp5 = ( (v(k,j+offy_v,i+offx_v,icrm)-v0(k,icrm)-v(kb,j+offy_v,i+offx_v,icrm)+v0(kb,icrm))*rdzw_dn);
      real tmp6 = ( (u(k,j+offy_u,ic+offx_u,icrm)-u0(k,icrm)-u(kb,j+offy_u,ic+offx_u,icrm)+u0(kb,icrm))*rdzw_dn);
      real tmp7 = ( (w(k,j+offy_w,ic+offx_w,icrm)-w(k,j+offy_w,i+offx_w,icrm))*rdx_dn);
      real tmp8 = ( (u(k,j+offy_u,i+offx_u,icrm)-u0(k,icrm)-u(kb,j+offy_u,i+offx_u,icrm)+u0(kb,icrm))*rdzw_dn+
                    (w(k,j+offy_w,i+offx_w,icrm)-w(k,j+offy_w,ib+offx_w,icrm))*rdx_dn);
      def2(k,j,i,icrm) = 2.0 * ( tmp1*tmp1 + tmp2*tmp2 ) + 0.5 * ( tmp3*tmp3 + tmp4*tmp4 ) + tmp5*tmp5
                                                         + 0.5 * ( tmp6*tmp6 + tmp7*tmp7 + tmp8*tmp8 );

    }
  });
}


