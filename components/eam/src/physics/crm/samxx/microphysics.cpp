#include "microphysics.h"

void precip_fall(int hydro_type, real4d &omega) {
  auto &rho           = :: rho;
  auto &adz           = :: adz;
  auto &dtn           = :: dtn;
  auto &dz            = :: dz;
  auto &micro_field   = :: micro_field;
  auto &rhow          = :: rhow;
  auto &qpfall        = :: qpfall;
  auto &tlat          = :: tlat;
  auto &precflux      = :: precflux;
  auto &precsfc       = :: precsfc;
  auto &precssfc      = :: precssfc;
  auto &prec_xy       = :: prec_xy;
  auto &t             = :: t;
  auto &vrain         = :: vrain;
  auto &vsnow         = :: vsnow;
  auto &vgrau         = :: vgrau;
  auto &crain         = :: crain;
  auto &csnow         = :: csnow;
  auto &cgrau         = :: cgrau;
  auto &tabs          = :: tabs;
  auto &a_pr          = :: a_pr;
  auto &a_gr          = :: a_gr;
  auto &ncrms         = :: ncrms;

  real constexpr eps = 1.e-10;
  bool constexpr nonos = true;

  real4d mx("mx",nzm,ny,nx,ncrms);
  real4d mn("mn",nzm,ny,nx,ncrms);
  real4d lfac("lfac",nz,ny,nx,ncrms);
  real4d www("www",nz,ny,nx,ncrms);
  real4d fz("fz",nz,ny,nx,ncrms);
  real4d wp("wp",nzm,ny,nx,ncrms);
  real4d tmp_qp("tmp_qp",nzm,ny,nx,ncrms);
  real2d irhoadz("irhoadz",nzm,ncrms);
  real2d iwmax("iwmax",nzm,ncrms);
  real2d rhofac("rhofac",nzm,ncrms);

  // for (int k=0; k<nzm; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    rhofac(k,icrm) = sqrt(1.29/rho(k,icrm));
    irhoadz(k,icrm) = 1.0/(rho(k,icrm)*adz(k,icrm));
    int kb = max(0,k-1);
    real wmax       = dz(icrm)*adz(kb,icrm)/dtn;   // Velocity equivalent to a cfl of 1.0.
    iwmax(k,icrm)   = 1.0/wmax;
  });

  //  Add sedimentation of precipitation field to the vert. vel.
  real prec_cfl = 0.0;
  real4d prec_cfl_arr("prec_cfl_arr",nzm,ny,nx,ncrms);

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    if (hydro_type == 0) {
      lfac(k,j,i,icrm) = fac_cond;
    }
    else if (hydro_type == 1) {
      lfac(k,j,i,icrm) = fac_sub;
    }
    else if (hydro_type == 2) {
      lfac(k,j,i,icrm) = fac_cond + (1.0-omega(k,j,i,icrm))*fac_fus;
    }
    else if (hydro_type == 3) {
      lfac(k,j,i,icrm) = 0.0;
    }
    real tmp;
    tmp = term_vel_qp(icrm,i,j,k,micro_field(1,k,j+offy_s,i+offx_s,icrm), 
                      vrain, vsnow, vgrau, crain, csnow, cgrau, rho(k,icrm),
                      tabs(k,j,i,icrm), a_pr, a_gr);
    wp(k,j,i,icrm)=rhofac(k,icrm)*tmp;
    tmp = wp(k,j,i,icrm)*iwmax(k,icrm);
    prec_cfl_arr(k,j,i,icrm) = tmp;
    wp(k,j,i,icrm) = -wp(k,j,i,icrm)*rhow(k,icrm)*dtn/dz(icrm);
    if (k == 0) {
      fz(nz-1,j,i,icrm)=0.0;
      www(nz-1,j,i,icrm)=0.0;
      lfac(nz-1,j,i,icrm)=0.0;
    }
  });

  yakl::ParallelMax<real,yakl::memDevice> pmax( nzm*ny*nx*ncrms );
  real prec_cfl_loc  = pmax( prec_cfl_arr.data() );
  prec_cfl = max(prec_cfl , prec_cfl_loc);

  // If maximum CFL due to precipitation velocity is greater than 0.9,
  // take more than one advection step to maintain stability.
  int nprec;
  if (prec_cfl > 0.9) {
    nprec = ceil(prec_cfl/0.9);
    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      // wp already includes factor of dt, so reduce it by a
      // factor equal to the number of precipitation steps.
      wp(k,j,i,icrm) = wp(k,j,i,icrm)/nprec;
    });
  } else {
    nprec = 1;
  }

  //if(nprec > 1){std::cout << nprec << std::endl;}

#ifdef MMF_FIXED_SUBCYCLE
    nprec = 4;
#endif
  

  for(int iprec = 1; iprec<=nprec; iprec++) {
    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      tmp_qp(k,j,i,icrm) = micro_field(1,k,j+offy_s,i+offx_s,icrm); // Temporary array for qp in this column
    });

    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      if (nonos) {
        int kc=min(nzm-1,k+1);
        int kb=max(0,k-1);
        mx(k,j,i,icrm)=max(tmp_qp(kb,j,i,icrm),max(tmp_qp(kc,j,i,icrm),tmp_qp(k,j,i,icrm)));
        mn(k,j,i,icrm)=min(tmp_qp(kb,j,i,icrm),min(tmp_qp(kc,j,i,icrm),tmp_qp(k,j,i,icrm)));
      }
      // Define upwind precipitation flux
      fz(k,j,i,icrm)=tmp_qp(k,j,i,icrm)*wp(k,j,i,icrm);
    });

    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int kc=k+1;
      tmp_qp(k,j,i,icrm)=tmp_qp(k,j,i,icrm)-(fz(kc,j,i,icrm)-fz(k,j,i,icrm))*irhoadz(k,icrm); //Update temporary qp
    });

    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      // Also, compute anti-diffusive correction to previous
      // (upwind) approximation to the flux
      int kb=max(0,k-1);
      // The precipitation velocity is a cell-centered quantity,
      // since it is computed from the cell-centered
      // precipitation mass fraction.  Therefore, a reformulated
      // anti-diffusive flux is used here which accounts for
      // this and results in reduced numerical diffusion.
      www(k,j,i,icrm) = 0.5*(1.0+wp(k,j,i,icrm)*irhoadz(k,icrm))*(tmp_qp(kb,j,i,icrm)*wp(kb,j,i,icrm) - 
                        tmp_qp(k,j,i,icrm)*wp(k,j,i,icrm)); // works for wp(k)<0
    });

    if (nonos) {
      // for (int k=0; k<nzm; k++) {
      //   for (int j=0; j<ny; j++) {
      //     for (int i=0; i<nx; i++) {
      //       for (int icrm=0; icrm<ncrms; icrm++) {
      parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
        int kc=min(nzm-1,k+1);
        int kb=max(0,k-1);
        mx(k,j,i,icrm)=max(tmp_qp(kb,j,i,icrm),max(tmp_qp(kc,j,i,icrm),max(tmp_qp(k,j,i,icrm),mx(k,j,i,icrm))));
        mn(k,j,i,icrm)=min(tmp_qp(kb,j,i,icrm),min(tmp_qp(kc,j,i,icrm),min(tmp_qp(k,j,i,icrm),mn(k,j,i,icrm))));
        kc=min(nzm-1,k+1);
        mx(k,j,i,icrm)=rho(k,icrm)*adz(k,icrm)*(mx(k,j,i,icrm)-tmp_qp(k,j,i,icrm))/(pn(www(kc,j,i,icrm)) + 
                                                                                    pp(www(k,j,i,icrm))+eps);
        mn(k,j,i,icrm)=rho(k,icrm)*adz(k,icrm)*(tmp_qp(k,j,i,icrm)-mn(k,j,i,icrm))/(pp(www(kc,j,i,icrm)) + 
                                                                                    pn(www(k,j,i,icrm))+eps);
      });

      // for (int k=0; k<nzm; k++) {
      //   for (int j=0; j<ny; j++) {
      //     for (int i=0; i<nx; i++) {
      //       for (int icrm=0; icrm<ncrms; icrm++) {
      parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
        int kb=max(0,k-1);
        // Add limited flux correction to fz(k).
        fz(k,j,i,icrm) = fz(k,j,i,icrm) + pp(www(k,j,i,icrm))*min(1.0,min(mx(k,j,i,icrm), mn(kb,j,i,icrm))) -
                                        pn(www(k,j,i,icrm))*min(1.0,min(mx(kb,j,i,icrm),mn(k,j,i,icrm))); // Anti-diffusive flux
      });
    }

    // Update precipitation mass fraction and liquid-ice static
    // energy using precipitation fluxes computed in this column.

    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int kc=k+1;
      // Update precipitation mass fraction.
      // Note that fz is the total flux, including both the
      // upwind flux and the anti-diffusive correction.
      real flagstat = 1.0;
      micro_field(1,k,j+offy_s,i+offx_s,icrm)=micro_field(1,k,j+offy_s,i+offx_s,icrm)-(fz(kc,j,i,icrm)-
                                                                                       fz(k,j,i,icrm))*irhoadz(k,icrm);
      real tmp = -(fz(kc,j,i,icrm)-fz(k,j,i,icrm))*irhoadz(k,icrm)*flagstat;  // For qp budget
      yakl::atomicAdd(qpfall(k,icrm),tmp);
      real lat_heat = -(lfac(kc,j,i,icrm)*fz(kc,j,i,icrm)-lfac(k,j,i,icrm)*fz(k,j,i,icrm))*irhoadz(k,icrm);
      t(k,j+offy_s,i+offx_s,icrm)=t(k,j+offy_s,i+offx_s,icrm)-lat_heat;
      yakl::atomicAdd(tlat(k,icrm),-lat_heat);
      tmp = fz(k,j,i,icrm)*flagstat;
      yakl::atomicAdd(precflux(k,icrm),-tmp);
      if (k == 0) {
        precsfc(j,i,icrm) = precsfc(j,i,icrm) - fz(0,j,i,icrm)*flagstat; 
        precssfc(j,i,icrm) = precssfc(j,i,icrm) - fz(0,j,i,icrm)*(1.0-omega(0,j,i,icrm))*flagstat;
        prec_xy(j,i,icrm) = prec_xy(j,i,icrm) - fz(0,j,i,icrm)*flagstat;
      }
    });

    if (iprec < nprec) {
      // Re-compute precipitation velocity using new value of qp.
      // for (int j=0; j<ny; j++) {
      //  for (int i=0; i<nx; i++) {
      //    for (int k=0; k<nzm; k++) {
      //       for (int icrm=0; icrm<ncrms; icrm++) {
      parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
        real tmp = term_vel_qp(icrm,i,j,k,micro_field(1,k,j+offy_s,i+offx_s,icrm), 
                               vrain, vsnow, vgrau, crain, csnow, cgrau, rho(k,icrm),
                               tabs(k,j,i,icrm), a_pr, a_gr);
        wp(k,j,i,icrm) = rhofac(k,icrm)*tmp;
        // Decrease precipitation velocity by factor of nprec
        wp(k,j,i,icrm) = -wp(k,j,i,icrm)*rhow(k,icrm)*dtn/dz(icrm)/nprec;
        // Note: Don't bother checking CFL condition at each
        // substep since it's unlikely that the CFL will
        // increase very much between substeps when using
        // monotonic advection schemes.
        if (k == 0) {
          fz(nz-1,j,i,icrm)=0.0;
          www(nz-1,j,i,icrm)=0.0;
          lfac(nz-1,j,i,icrm)=0.0;
        }
      });
    }
  } // iprec loop
}


void micro_precip_fall() {
  auto &tabs  = ::tabs;
  auto &a_pr  = ::a_pr;
  auto &ncrms = ::ncrms;

  real4d omega("omega", nzm, ny, nx, ncrms);

  crain = b_rain / 4.0;
  csnow = b_snow / 4.0;
  cgrau = b_grau / 4.0;
  vrain = a_rain * gamr3 / 6.0 / pow((pi * rhor * nzeror), crain);
  vsnow = a_snow * gams3 / 6.0 / pow((pi * rhos * nzeros), csnow);
  vgrau = a_grau * gamg3 / 6.0 / pow((pi * rhog * nzerog), cgrau);

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    omega(k,j,i,icrm) = max(0.0,min(1.0,(tabs(k,j,i,icrm)-tprmin)*a_pr));
  });

  precip_fall(2,omega);
}


void micro_flux() {
  auto &fluxbmk       = :: fluxbmk;
  auto &fluxbq        = :: fluxbq;
  auto &fluxtmk       = :: fluxtmk;
  auto &fluxtq        = :: fluxtq;
  auto &ncrms         = :: ncrms;

  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(ny,nx,ncrms) , YAKL_LAMBDA (int j, int i, int icrm) {
    fluxbmk(index_water_vapor,j,i,icrm) = fluxbq(j,i,icrm);
    fluxtmk(index_water_vapor,j,i,icrm) = fluxtq(j,i,icrm);
  });
}

void micro_diagnose() {
  auto &qv            = :: qv;
  auto &micro_field   = :: micro_field;
  auto &qn            = :: qn;
  auto &tabs          = :: tabs;
  auto &a_bg          = :: a_bg;
  auto &a_pr          = :: a_pr;
  auto &qcl           = :: qcl;
  auto &qci           = :: qci;
  auto &qpl           = :: qpl;
  auto &qpi           = :: qpi;
  auto &ncrms         = :: ncrms;

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    qv(k,j,i,icrm) = micro_field(0,k,j+offy_s,i+offx_s,icrm) - qn(k,j,i,icrm);
    real omn = max(0.0,min(1.0,(tabs(k,j,i,icrm)-tbgmin)*a_bg));
    qcl(k,j,i,icrm) = qn(k,j,i,icrm)*omn;
    qci(k,j,i,icrm) = qn(k,j,i,icrm)*(1.0-omn);
    real omp = max(0.0,min(1.0,(tabs(k,j,i,icrm)-tprmin)*a_pr));
    qpl(k,j,i,icrm) = micro_field(1,k,j+offy_s,i+offx_s,icrm)*omp;
    qpi(k,j,i,icrm) = micro_field(1,k,j+offy_s,i+offx_s,icrm)*(1.0-omp);
  });
}

void micro_proc() {
  if (doprecip && icycle == 1) {
    precip_init();
  }

  if (docloud) {
    cloud(micro_field, 0, micro_field, 1);
    if (doprecip) {
      precip_proc(micro_field,0,micro_field,1);
    }
    micro_diagnose();
  }

  if (dosmoke) {
    micro_diagnose();
  }
}

void micro_init() {
  auto &fluxbmk = ::fluxbmk;
  auto &fluxtmk = ::fluxtmk;
  auto &mkwle   = ::mkwle;
  auto &mkwsb   = ::mkwsb;
  auto &mkadv   = ::mkadv;
  auto &mkdiff  = ::mkdiff;
  auto &qpsrc   = ::qpsrc;
  auto &qpevp   = ::qpevp;
  auto &ncrms   = ::ncrms;

  a_bg = 1.0/(tbgmax-tbgmin);
  a_pr = 1.0/(tprmax-tprmin);
  a_gr = 1.0/(tgrmax-tgrmin);

  if (nrestart == 0) {
    // for (int l=0; l<nmicro_fields; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nmicro_fields,ny,nx,ncrms) , YAKL_LAMBDA (int l, int j, int i, int icrm) {
      fluxbmk(l,j,i,icrm) = 0.0;
      fluxtmk(l,j,i,icrm) = 0.0;
    });
    if (docloud || dosmoke) {
      micro_diagnose();
    }
  }

  // for (int l=0; l<nmicro_fields; k++) {
  //  for (int k=0; k<nz; k++) {
  //    for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(nmicro_fields,nz,ncrms) , YAKL_LAMBDA (int l, int k, int icrm) {
    mkwle (l,k,icrm) = 0.0;
    mkwsb (l,k,icrm) = 0.0;
    mkadv (l,k,icrm) = 0.0;
    mkdiff(l,k,icrm) = 0.0;
  });
  
  // for (int k=0; k<nz; k++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nz,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    qpsrc(k,icrm) = 0.0;
    qpevp(k,icrm) = 0.0;
  });
}

