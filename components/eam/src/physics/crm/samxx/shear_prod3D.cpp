
#include "shear_prod3D.h"

void shear_prod3D(real4d &def2) {
  YAKL_SCOPE( dx    , ::dx );
  YAKL_SCOPE( dy    , ::dy );
  YAKL_SCOPE( dz    , ::dz );
  YAKL_SCOPE( adz   , ::adz );
  YAKL_SCOPE( adzw  , ::adzw );
  YAKL_SCOPE( u     , ::u );
  YAKL_SCOPE( v     , ::v );
  YAKL_SCOPE( w     , ::w );
  YAKL_SCOPE( u0    , ::u0 );
  YAKL_SCOPE( v0    , ::v0 );
  YAKL_SCOPE( ncrms , ::ncrms );
  
  // for (int k=0; k<nzm; k++) {
  //    for (int j=0; j<ny; j++) {
  //      for (int i=0; i<nx; i++) {
  //        for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    real rdx0 = 1.0/dx;
    real rdy0 = 1.0/dy;
    real rdz, rdzw_up, rdzw_dn, rdx, rdx_up, rdx_dn, rdy, rdy_up, rdy_dn;
    int kb, kc, ib, ic, jb, jc;

    if(k>=1 && k <= nzm-2) {

      kb = k-1;
      kc = k+1;
      rdz = 1.0/(dz(icrm)*adz(k,icrm));
      rdzw_up = 1.0/(dz(icrm)*adzw(kc,icrm));
      rdzw_dn = 1.0/(dz(icrm)*adzw(k,icrm));
      rdx = rdx0*sqrt(dx*rdz);
      rdy = rdy0*sqrt(dy*rdz);
      rdx_up = rdx0 *sqrt(dx*rdzw_up);
      rdy_up = rdy0 *sqrt(dy*rdzw_up);
      rdx_dn = rdx0 *sqrt(dx*rdzw_dn);
      rdy_dn = rdy0 *sqrt(dy*rdzw_dn);
      jb = j-YES3D;
      jc = j+YES3D;
      ib = i-1;
      ic = i+1;
      real tmp1 = ((u(k,j+offy_u,ic+offx_u,icrm)-u(k,j+offy_u,i+offx_u,icrm))*rdx);
      real tmp2 = ((v(k,jc+offy_v,i+offx_v,icrm)-v(k,j+offy_v,i+offx_v,icrm))*rdy);
      real tmp3 = ((w(kc,j+offy_w,i+offx_w,icrm)-w(k,j+offy_w,i+offx_w,icrm))*rdz);
      real tmp4 = (u(k,jc+offy_u,ic+offx_u,icrm)-u(k,j+offy_v,ic+offx_v,icrm))*rdy +
                  (v(k,jc+offy_v,ic+offx_v,icrm)-v(k,jc+offy_v,i+offx_v,icrm))*rdx;
      real tmp5 = (u(k,jc+offy_u,i+offx_u,icrm)-u(k,j+offy_v,i+offx_v,icrm))*rdy +
                  (v(k,jc+offy_v,i+offx_v,icrm)-v(k,jc+offy_v,ib+offx_v,icrm))*rdx;
      real tmp6 = (u(k,j+offy_u,ic+offx_u,icrm)-u(k,jb+offy_v,ic+offx_v,icrm))*rdy +
                  (v(k,j+offy_v,ic+offx_v,icrm)-v(k,j+offy_v,i+offx_v,icrm))*rdx;
      real tmp7 = (u(k,j+offy_u,i+offx_u,icrm)-u(k,jb+offy_v,i+offx_v,icrm))*rdy +
                  (v(k,j+offy_v,i+offx_v,icrm)-v(k,j+offy_v,ib+offx_v,icrm))*rdx;
      def2(k,j,i,icrm) = 2.0 * ( tmp1*tmp1 + tmp2*tmp2 + tmp3*tmp3 ) + 
                         0.25 * ( tmp4*tmp4 + tmp5*tmp5 + tmp6*tmp6 + tmp7*tmp7);

      tmp1 = (u(kc,j+offy_u,ic+offx_u,icrm)-u0(kc,icrm)-u(k,j+offy_u,ic+offx_u,icrm)+u0(k,icrm))*rdzw_up +
             (w(kc,j+offy_w,ic+offx_w,icrm)-w(kc,j+offy_w,i+offx_w,icrm))*rdx_up;
      tmp2 = (u(kc,j+offy_u,i+offx_u,icrm)-u0(kc,icrm)-u(k,j+offy_u,i+offx_u,icrm)+u0(k,icrm))*rdzw_up +
             (w(kc,j+offy_w,i+offx_w,icrm)-w(kc,j+offy_w,ib+offx_w,icrm))*rdx_up;
      tmp3 = (u(k,j+offy_u,ic+offx_u,icrm)-u0(k,icrm)-u(kb,j+offy_u,ic+offx_u,icrm)+u0(kb,icrm))*rdzw_dn +
             (w(k,j+offy_w,ic+offx_w,icrm)-w(k,j+offy_w,i+offx_w,icrm))*rdx_dn;
      tmp4 = (u(k,j+offy_u,i+offx_u,icrm)-u0(k,icrm)-u(kb,j+offy_u,i+offx_u,icrm)+u0(kb,icrm))*rdzw_dn +
             (w(k,j+offy_w,i+offx_w,icrm)-w(k,j+offy_w,ib+offx_w,icrm))*rdx_dn;
      def2(k,j,i,icrm) = def2(k,j,i,icrm) + .25 * ( tmp1*tmp1 + tmp2*tmp2 + tmp3*tmp3 + tmp4*tmp4 );

      tmp1 = (v(kc,jc+offy_v,i+offx_v,icrm)-v0(kc,icrm)-v(k,jc+offy_v,i+offx_v,icrm)+v0(k,icrm))*rdzw_up +
             (w(kc,jc+offy_w,i+offx_w,icrm)-w(kc,j+offy_w,i+offx_w,icrm))*rdy_up;
      tmp2 = (v(kc,j+offy_v,i+offx_v,icrm)-v0(kc,icrm)-v(k,j+offy_v,i+offx_v,icrm)+v0(k,icrm))*rdzw_up +
             (w(kc,j+offy_w,i+offx_w,icrm)-w(kc,jb+offy_w,i+offx_w,icrm))*rdy_up;
      tmp3 = (v(k,jc+offy_v,i+offx_v,icrm)-v0(k,icrm)-v(kb,jc+offy_v,i+offx_v,icrm)+v0(kb,icrm))*rdzw_dn +
             (w(k,jc+offy_w,i+offx_w,icrm)-w(k,j+offy_w,i+offx_w,icrm))*rdy_dn;
      tmp4 = (v(k,j+offy_v,i+offx_v,icrm)-v0(k,icrm)-v(kb,j+offy_v,i+offx_v,icrm)+v0(kb,icrm))*rdzw_dn +
             (w(k,j+offy_w,i+offx_w,icrm)-w(k,jb+offy_w,i+offx_w,icrm))*rdy_dn;
      def2(k,j,i,icrm) = def2(k,j,i,icrm) + .25 * ( tmp1*tmp1 + tmp2*tmp2 + tmp3*tmp3 + tmp4*tmp4 );

    }  else if (k==0) {

      kc = k+1;
      rdz = 1.0/(dz(icrm)*adz(k,icrm));
      rdzw_up = 1.0/(dz(icrm)*adzw(kc,icrm));
      rdx = rdx0*sqrt(dx*rdz);
      rdy = rdy0*sqrt(dy*rdz);
      rdx_up = rdx0 *sqrt(dx*rdzw_up);
      rdy_up = rdy0 *sqrt(dy*rdzw_up);
      jb = j-YES3D;
      jc = j+YES3D;
      ib = i-1;
      ic = i+1;

      real tmp1  = (u(k,j+offy_u,ic+offx_u,icrm)-u(k,j+offy_u,i+offx_u,icrm))*rdx;
      real tmp2  = (v(k,jc+offy_v,i+offx_v,icrm)-v(k,j+offy_v,i+offx_v,icrm))*rdy;
      real tmp3  = (w(kc,j+offy_w,i+offx_w,icrm)-w(k,j+offy_w,i+offx_w,icrm))*rdz;
      real tmp4  = (u(k,jc+offy_u,ic+offx_u,icrm)-u(k,j+offy_u,ic+offx_u,icrm))*rdy + 
                   (v(k,jc+offy_v,ic+offx_v,icrm)-v(k,jc+offy_v,i+offx_v,icrm))*rdx;
      real tmp5  = (u(k,jc+offy_u,i+offx_u,icrm)-u(k,j+offy_u,i+offx_u,icrm))*rdy + 
                   (v(k,jc+offy_v,i+offx_v,icrm)-v(k,jc+offy_v,ib+offx_v,icrm))*rdx  ;
      real tmp6  = (u(k,j+offy_u,ic+offx_u,icrm)-u(k,jb+offy_u,ic+offx_u,icrm))*rdy + 
                   (v(k,j+offy_v,ic+offx_v,icrm)-v(k,j+offy_v,i+offx_v,icrm))*rdx  ;
      real tmp7  = (u(k,j+offy_u,i+offx_u,icrm)-u(k,jb+offy_u,i+offx_u,icrm))*rdy + 
                   (v(k,j+offy_v,i+offx_v,icrm)-v(k,j+offy_v,ib+offx_v,icrm))*rdx    ;
      real tmp8  = (v(kc,jc+offy_v,i+offx_v,icrm)-v0(kc,icrm)-v(k,jc+offy_v,i+offx_v,icrm)+v0(k,icrm))*rdzw_up +
                   (w(kc,jc+offy_w,i+offx_w,icrm)-w(kc,j+offy_w,i+offx_w,icrm))*rdy_up;
      real tmp9  = (v(kc,j+offy_v,i+offx_v,icrm)-v0(kc,icrm)-v(k,j+offy_v,i+offx_v,icrm)+v0(k,icrm))*rdzw_up +
                   (w(kc,j+offy_w,i+offx_w,icrm)-w(kc,jb+offy_w,i+offx_w,icrm))*rdy_up;
      real tmp10 = (u(kc,j+offy_u,ic+offx_u,icrm)-u0(kc,icrm)-u(k,j+offy_u,ic+offx_u,icrm)+u0(k,icrm))*rdzw_up +
                   (w(kc,j+offy_w,ic+offx_w,icrm)-w(kc,j+offy_w,i+offx_w,icrm))*rdx_up;
      real tmp11 = (u(kc,j+offy_u,i+offx_u,icrm)-u0(kc,icrm)-u(k,j+offy_u,i+offx_u,icrm)+u0(k,icrm))*rdzw_up +
                   (w(kc,j+offy_w,i+offx_w,icrm)-w(kc,j+offy_w,ib+offx_w,icrm))*rdx_up;

      def2(k,j,i,icrm) = 2.0 * ( tmp1*tmp1 + tmp2*tmp2 + tmp3*tmp3 ) + 
                        0.25 * ( tmp4*tmp4 + tmp5*tmp5 + tmp6*tmp6 + tmp7*tmp7 ) +
                         0.5 * ( tmp8*tmp8 + tmp9*tmp9 ) + 0.5 * ( tmp10*tmp10 + tmp11*tmp11 );

    } else if (k==nzm-1) {

      kc = k+1;
      kb = k-1;
      rdz = 1.0/(dz(icrm)*adz(k,icrm));
      rdzw_dn = 1.0/(dz(icrm)*adzw(k,icrm));
      rdx = rdx0*sqrt(dx*rdz);
      rdx = rdy0*sqrt(dy*rdz);
      rdx_dn = rdx0 *sqrt(dx*rdzw_dn);
      rdy_dn = rdy0 *sqrt(dy*rdzw_dn);
      jb = j-1*YES3D;
      jc = j+1*YES3D;
      ib = i-1;
      ic = i+1;

      real tmp1  = (u(k,j+offy_u,ic+offx_u,icrm)-u(k,j+offy_u,i+offx_u,icrm))*rdx;
      real tmp2  = (v(k,jc+offy_v,i+offx_v,icrm)-v(k,j+offy_v,i+offx_v,icrm))*rdy;
      real tmp3  = (w(kc,j+offy_w,i+offx_w,icrm)-w(k,j+offy_w,i+offx_w,icrm))*rdz;
      real tmp4  = (u(k,jc+offy_u,ic+offx_u,icrm)-u(k,j+offy_u,ic+offx_u,icrm))*rdy + 
                   (v(k,jc+offy_v,ic+offx_v,icrm)-v(k,jc+offy_v,i+offx_v,icrm))*rdx;
      real tmp5  = (u(k,jc+offy_u,i+offx_u,icrm)-u(k,j+offy_u,i+offx_u,icrm))*rdy + 
                   (v(k,jc+offy_v,i+offx_v,icrm)-v(k,jc+offy_v,ib+offx_v,icrm))*rdx  ;
      real tmp6  = (u(k,j+offy_u,ic+offx_u,icrm)-u(k,jb+offy_u,ic+offx_u,icrm))*rdy + 
                   (v(k,j+offy_v,ic+offx_v,icrm)-v(k,j+offy_v,i+offx_v,icrm))*rdx  ;
      real tmp7  = (u(k,j+offy_u,i+offx_u,icrm)-u(k,jb+offy_u,i+offx_u,icrm))*rdy + 
                   (v(k,j+offy_v,i+offx_v,icrm)-v(k,j+offy_v,ib+offx_v,icrm))*rdx    ;
      real tmp8  = (v(k,jc+offy_v,i+offx_v,icrm)-v0(k,icrm)-v(kb,jc+offy_v,i+offx_v,icrm)+v0(kb,icrm))*rdzw_dn +
                   (w(k,jc+offy_w,i+offx_w,icrm)-w(k,j+offy_w,i+offx_w,icrm))*rdy_dn;
      real tmp9  = (v(k,j+offy_v,i+offx_v,icrm)-v0(k,icrm)-v(kb,j+offy_v,i+offx_v,icrm)+v0(kb,icrm))*rdzw_dn +
                   (w(k,j+offy_w,i+offx_w,icrm)-w(k,jb+offy_w,i+offx_w,icrm))*rdy_dn  ;
      real tmp10 = (u(k,j+offy_u,ic+offx_u,icrm)-u0(k,icrm)-u(kb,j+offy_u,ic+offx_u,icrm)+u0(kb,icrm))*rdzw_dn +
                   (w(k,j+offy_w,ic+offx_w,icrm)-w(k,j+offy_w,i+offx_w,icrm))*rdx_dn;
      real tmp11 = (u(k,j+offy_u,i+offx_u,icrm)-u0(k,icrm)-u(kb,j+offy_u,i+offx_u,icrm)+u0(kb,icrm))*rdzw_dn +
                   (w(k,j+offy_w,i+offx_w,icrm)-w(k,j+offy_w,ib+offx_w,icrm))*rdx_dn  ;

      def2(k,j,i,icrm) = 2.0 * ( tmp1*tmp1 + tmp2*tmp2 + tmp3*tmp3 ) +
                        0.25 * ( tmp4*tmp4 + tmp5*tmp5 + tmp6*tmp6 + tmp7*tmp7 ) +
                         0.5 * ( tmp8*tmp8 + tmp9*tmp9 ) + 0.5 * ( tmp10*tmp10 + tmp11*tmp11 );

    }
  });
}


