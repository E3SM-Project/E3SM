
#include "kurant.h"
#include "vars.h"
#include "samxx_utils.h"

void kurant () {
  YAKL_SCOPE( w     , ::w );
  YAKL_SCOPE( u     , ::u );
  YAKL_SCOPE( v     , ::v );
  YAKL_SCOPE( dt    , ::dt );
  YAKL_SCOPE( dx    , ::dx );
  YAKL_SCOPE( dy    , ::dy );
  YAKL_SCOPE( dz    , ::dz );
  YAKL_SCOPE( adzw  , ::adzw );
  YAKL_SCOPE( ncrms , ::ncrms );
  YAKL_SCOPE( tabs  , ::tabs );
  YAKL_SCOPE( qv    , ::qv );
  YAKL_SCOPE( qcl   , ::qcl );
  YAKL_SCOPE( qci   , ::qci );
  YAKL_SCOPE( micro_field, :: micro_field );
  YAKL_SCOPE( longitude0 , :: longitude0 );
  YAKL_SCOPE( latitude0  , :: latitude0 );
  YAKL_SCOPE( microphysics_scheme, :: microphysics_scheme );

  int constexpr max_ncycle = 4;
  real cfl;

  real2d wmax("wmax" ,nz ,ncrms);   // max vertical velocity
  real2d umax("umax" ,nzm ,ncrms);  // max horizontal wind magnitude
  real2d tmpMax("tmpMax",nzm,ncrms);

  ncycle = 1;

  parallel_for( SimpleBounds<2>(nz,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    wmax(k,icrm) = 0.0;
  });
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    umax(k,icrm) = 0.0;
  });

  parallel_for( SimpleBounds<4>(nz,ny,nx,ncrms) , YAKL_DEVICE_LAMBDA (int k, int j, int i, int icrm) {
    real tmp = fabs(w(k,j+offy_w,i+offx_w,icrm));
    yakl::atomicMax(wmax(k,icrm),tmp);
  });
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_DEVICE_LAMBDA (int k, int j, int i, int icrm) {
    real utmp = u(k,j+offy_u,i+offx_u,icrm);
    real vtmp = v(k,j+offy_v,i+offx_v,icrm);
    real tmp = sqrt(utmp*utmp +YES3D*vtmp*vtmp);
    yakl::atomicMax(umax(k,icrm),tmp);
  });


  cfl = 0.0;
  // for (int k=0; k<nzm; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    real tmp1 = umax(k,icrm)*dt*sqrt(1.0/(dx*dx) + YES3D*1.0/(dy*dy));
    real dztemp = dz(icrm)*adzw(k,icrm);
    real tmp2 = wmax(k,icrm)*dt/dztemp;
    real tmp3 = wmax(k+1,icrm)*dt/dztemp;
    tmpMax(k,icrm) = max(max(tmp1,tmp2),tmp3);
  });

  yakl::ParallelMax<real,yakl::memDevice> pmax( nzm*ncrms );
  real cfl_loc = pmax(tmpMax.data());
  cfl = max(cfl,cfl_loc);


  if(cfl != cfl) {
    std::cout << "\nkurant() - cfl is NaN." << std::endl;
    finalize();
    exit(-1);
  }

  if (is_same_str(turbulence_scheme, "smag") == 0) { kurant_sgs(cfl); }

  ncycle = max(ncycle,max(1,static_cast<int>(ceil(cfl/0.7))));

#ifdef MMF_FIXED_SUBCYCLE
  ncycle = max_ncycle;
#endif

  //----------------------------------------------------------------------------
  // if ncycle too large, print some debugging info and exit
  //----------------------------------------------------------------------------
  if(ncycle > max_ncycle) {

    real2d tamax("tamax" ,nzm ,ncrms);  // max absolute temperature
    real2d qvmax("qvmax" ,nzm ,ncrms);  // max specific humidity
    real2d qcmax("qcmax" ,nzm ,ncrms);  // max liq cloud water
    real2d qimax("qimax" ,nzm ,ncrms);  // max ice cloud water

    parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
      tamax(k,icrm) = 0.0;
      qvmax(k,icrm) = 0.0;
      qcmax(k,icrm) = 0.0;
      qimax(k,icrm) = 0.0;
    });

    parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_DEVICE_LAMBDA (int k, int j, int i, int icrm) {
      real tmp;
      tmp = tabs(k,j,i,icrm); yakl::atomicMax(tamax(k,icrm),tmp);
      tmp = qv(k,j,i,icrm);   yakl::atomicMax(qvmax(k,icrm),tmp);
      tmp = qcl(k,j,i,icrm);  yakl::atomicMax(qcmax(k,icrm),tmp);
      tmp = qci(k,j,i,icrm);  yakl::atomicMax(qimax(k,icrm),tmp);
    });

    std::cout << "\nkurant() - the number of cycles exceeded max_ncycle = "<< max_ncycle << std::endl;

    for (int icrm=0; icrm<ncrms; icrm++) {
      for (int k=0; k<nzm; k++) {
        std::cout<<"  "
        <<"  icrm:"<<icrm
        <<"  k:"<<k
        <<"  lat: "<<latitude0(icrm)
        <<"  lon: "<<longitude0(icrm)
        <<"  wmax: "<<wmax(k,icrm)
        <<"  umax: "<<umax(k,icrm)
        <<"  tamax: "<<tamax(k,icrm)
        <<"  qvmax: "<<qvmax(k,icrm)
        <<"  qcmax: "<<qcmax(k,icrm)
        <<"  qimax: "<<qimax(k,icrm)
        << std::endl;
      }
    }

    finalize();
    exit(-1);
  }
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
}


