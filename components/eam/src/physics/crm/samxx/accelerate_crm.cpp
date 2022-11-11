
#include "accelerate_crm.h"

void accelerate_crm(int nstep, int nstop, bool &ceaseflag) {
  YAKL_SCOPE( t                  , ::t);
  YAKL_SCOPE( qcl                , ::qcl);
  YAKL_SCOPE( qci                , ::qci);
  YAKL_SCOPE( qv                 , ::qv);
  YAKL_SCOPE( u                  , ::u);
  YAKL_SCOPE( v                  , ::v);
  YAKL_SCOPE( t0                 , ::t0);
  YAKL_SCOPE( u0                 , ::u0);
  YAKL_SCOPE( v0                 , ::v0);
  YAKL_SCOPE( q0                 , ::q0);
  YAKL_SCOPE( crm_accel_uv       , ::crm_accel_uv);
  YAKL_SCOPE( use_crm_accel      , ::use_crm_accel);
  YAKL_SCOPE( micro_field        , ::micro_field);
  YAKL_SCOPE( ncrms              , ::ncrms);
  YAKL_SCOPE( crm_accel_factor   , ::crm_accel_factor);

  real ttend_threshold = 5.0;  // 5K, following UP-CAM implementation
  real tmin = 50.0;  // should never get below 50K in crm, following UP-CAM implementation
  int idx_qt = index_water_vapor;

  real2d ubaccel("ubaccel", nzm, ncrms);
  real2d vbaccel("vbaccel", nzm, ncrms);
  real2d tbaccel("tbaccel", nzm, ncrms);
  real2d qtbaccel("qtbaccel", nzm, ncrms);
  real2d ttend_acc("ttend_acc", nzm, ncrms);
  real2d qtend_acc("qtend_acc", nzm, ncrms);
  real2d utend_acc("utend_acc", nzm, ncrms);
  real2d vtend_acc("vtend_acc", nzm, ncrms);
  real2d qpoz("qpoz", nzm, ncrms);
  real2d qneg("qneg", nzm, ncrms);

  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // Compute the average among horizontal columns for each variable
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  // for (int k=0; k<nzm; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    tbaccel(k,icrm) = 0.0;
    qtbaccel(k,icrm) = 0.0;
    if (crm_accel_uv) {
      ubaccel(k,icrm) = 0.0;
      vbaccel(k,icrm) = 0.0;
    }
  });

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    // calculate tendency * dtn
    yakl::atomicAdd( tbaccel(k,icrm) , t(k,j+offy_s,i+offx_s,icrm) * crm_accel_coef );
    yakl::atomicAdd( qtbaccel(k,icrm) , (qcl(k,j,i,icrm) + qci(k,j,i,icrm) + qv(k,j,i,icrm)) * crm_accel_coef );
    if (crm_accel_uv) {
      yakl::atomicAdd( ubaccel(k,icrm) , u(k,j+offy_u,i+offx_u,icrm) * crm_accel_coef );
      yakl::atomicAdd( vbaccel(k,icrm) , v(k,j+offy_v,i+offx_v,icrm) * crm_accel_coef );
    }
  });

  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //!! Compute the accelerated tendencies
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ScalarLiveOut<bool> ceaseflag_liveout(false);

  // for (int k=0; k<nzm; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    ttend_acc(k,icrm) = tbaccel(k,icrm) - t0(k,icrm);
    qtend_acc(k,icrm) = qtbaccel(k,icrm) - q0(k,icrm);
    if (crm_accel_uv) {
      utend_acc(k,icrm) = ubaccel(k,icrm) - u0(k,icrm);
      vtend_acc(k,icrm) = vbaccel(k,icrm) - v0(k,icrm);
    }
    if (abs(ttend_acc(k,icrm)) > ttend_threshold) {
      ceaseflag_liveout = true;
    }
  });
  ceaseflag = ceaseflag_liveout.hostRead();


  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //!! Make sure it isn't insane
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (ceaseflag) { // special case for dT/dt too large
    // MSA will not be applied here or for the remainder of the CRM integration.
    // nstop must be updated to ensure the CRM integration duration is unchanged.
    // 
    // The effective MSA timestep is dt_a = crm_dt * (1 + crm_accel_factor). When
    // ceaseflag is triggered at nstep, we've taken (nstep - 1) previous steps of
    // size crm_dt * (1 + crm_accel_factor). The current step, and all future
    // steps, will revert to size crm_dt. Therefore, the total crm integration
    // time remaining after this step is
    //     time_remaining = crm_run_time - (nstep - 1)* dt_a + crm_dt
    //     nsteps_remaining = time_remaining / crm_dt
    //     updated nstop = nstep + nsteps_remaining
    // Because we set nstop = crm_run_time / dt_a in crm_accel_nstop, subbing
    // crm_run_time = nstop * dt_a and working through algebra yields 
    //     updated nstop = nstop + (nstop - nstep + 1) * crm_accel_factor.
    std::cout << "accelerate_crm: mean-state acceleration not applied this step";
    std::cout << "crm: nstop increased from " << nstop << " to " << round(nstop+(nstop-nstep+1)*crm_accel_factor);
    nstop = nstop + (nstop - nstep + 1)*crm_accel_factor; // only can happen once
    return;
  }

  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //!! Apply the accelerated tendencies
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    // don't let T go negative!
    t(k,j+offy_s,i+offx_s,icrm) = max(tmin, t(k,j+offy_s,i+offx_s,icrm) + crm_accel_factor * ttend_acc(k,icrm));
    if (crm_accel_uv) {
      u(k,j+offy_u,i+offx_u,icrm) = u(k,j+offy_u,i+offx_u,icrm) + crm_accel_factor * utend_acc(k,icrm); 
      v(k,j+offy_v,i+offx_v,icrm) = v(k,j+offy_v,i+offx_v,icrm) + crm_accel_factor * vtend_acc(k,icrm); 
    }
    micro_field(idx_qt,k,j+offy_s,i+offx_s,icrm) = 
        micro_field(idx_qt,k,j+offy_s,i+offx_s,icrm) + crm_accel_factor * qtend_acc(k,icrm);
  });

  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //!! Fix negative micro and readjust among separate water species
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  // for (int k=0; k<nzm; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    qpoz(k,icrm) = 0.0;
    qneg(k,icrm) = 0.0;
  });

  // separately accumulate positive and negative qt values in each layer k
  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    if (micro_field(idx_qt,k,j+offy_s,i+offx_s,icrm) < 0.0) {
      yakl::atomicAdd( qneg(k,icrm) , micro_field(idx_qt,k,j+offy_s,i+offx_s,icrm) ); 
    }
    else {
      yakl::atomicAdd( qpoz(k,icrm) , micro_field(idx_qt,k,j+offy_s,i+offx_s,icrm) );
    }
  });

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    real factor;
    if (qpoz(k,icrm) + qneg(k,icrm) <= 0.0) {
      // all moisture depleted in layer
      micro_field(idx_qt,k,j+offy_s,i+offx_s,icrm) = 0.0;
      qv(k,j,i,icrm) = 0.0;
      qcl(k,j,i,icrm) = 0.0;
      qci(k,j,i,icrm) = 0.0;
    } else {
      // Clip qt values at 0 and remove the negative excess in each layer
      // proportionally from the positive qt fields in the layer
      factor = 1.0 + qneg(k,icrm) / qpoz(k,icrm);
      micro_field(idx_qt,k,j+offy_s,i+offx_s,icrm) = max(0.0, micro_field(idx_qt,k,j+offy_s,i+offx_s,icrm) * factor);
      // Partition micro_field == qv + qcl + qci following these rules:
      //    (1) attempt to satisfy purely by adjusting qv
      //    (2) adjust qcl and qci only if needed to ensure positivity
      if (micro_field(idx_qt,k,j+offy_s,i+offx_s,icrm) < 0.0) {
        qv (k,j,i,icrm) = 0.0;
        qcl(k,j,i,icrm) = 0.0;
        qci(k,j,i,icrm) = 0.0;
      } else {
        // deduce qv as residual between qt - qcl - qci
        real qt_res = micro_field(idx_qt,k,j+offy_s,i+offx_s,icrm) - qcl(k,j,i,icrm) - qci(k,j,i,icrm);
        qv(k,j,i,icrm) = max(0.0, qt_res);
        if (qt_res < 0.0) {
          // qv was clipped; need to reduce qcl and qci accordingly
          factor = 1.0 + qt_res / (qcl(k,j,i,icrm) + qci(k,j,i,icrm));
          qcl(k,j,i,icrm) = qcl(k,j,i,icrm) * factor;
          qci(k,j,i,icrm) = qci(k,j,i,icrm) * factor;
        }
      }
    } 
  });
}



void crm_accel_nstop(int &nstop) {
  if(nstop%static_cast<int>((1+crm_accel_factor)) != 0) {
    std::cout << "CRM acceleration unexpected exception:\n";
    std::cout << "(1+crm_accel_factor) does not divide equally into nstop\n";
    std::cout << "nstop = " <<  nstop << std::endl;
    std::cout << "crm_accel_factor = " << crm_accel_factor << std::endl;
    exit(-1);
  } else {
    nstop = nstop / (1 + crm_accel_factor);
  }
}


