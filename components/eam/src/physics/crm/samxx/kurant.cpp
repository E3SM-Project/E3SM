
#include "kurant.h"
#include "vars.h"

void kurant () {
  auto &w     = ::w;
  auto &u     = ::u;
  auto &v     = ::v;
  auto &dt    = ::dt;
  auto &dx    = ::dx;
  auto &dy    = ::dy;
  auto &dz    = ::dz;
  auto &adzw  = ::adzw;
  auto &ncrms = ::ncrms;

  int constexpr max_ncycle = 4;
  real cfl;

  real2d wm    ("wm"   ,nz ,ncrms);
  real2d uhm   ("uhm"  ,nz ,ncrms);
  real2d tmpMax("uhMax",nzm,ncrms);

  ncycle = 1;
  parallel_for( SimpleBounds<2>(nz,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    wm(k,icrm) = 0.0;
    uhm(k,icrm) = 0.0;
  });

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    real tmp;
    tmp = fabs(w(k,j+offy_w,i+offx_w,icrm));
    yakl::atomicMax(wm(k,icrm),tmp);

    real utmp = u(k,j+offy_u,i+offx_u,icrm);
    real vtmp = v(k,j+offy_v,i+offx_v,icrm);
    tmp = sqrt(utmp*utmp +YES3D*vtmp*vtmp);
    yakl::atomicMax(uhm(k,icrm),tmp);
  });


  cfl = 0.0;
  // for (int k=0; k<nzm; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    real tmp1 = uhm(k,icrm)*dt*sqrt(1.0/(dx*dx) + YES3D*1.0/(dy*dy));
    real dztemp = dz(icrm)*adzw(k,icrm);
    real tmp2 = wm(k,icrm)*dt/dztemp;
    real tmp3 = wm(k+1,icrm)*dt/dztemp;
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

  kurant_sgs(cfl);

  ncycle = max(ncycle,max(1,static_cast<int>(ceil(cfl/0.7))));

  if(ncycle > max_ncycle) {
    std::cout << "\nkurant() - the number of cycles exceeded max_ncycle = "<< max_ncycle << std::endl;
    exit(-1);
  }
}


