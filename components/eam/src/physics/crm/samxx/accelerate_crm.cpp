
#include "accelerate_crm.h"
#include "samxx_utils.h"

void accelerate_crm(int nstep, int nstop, bool &ceaseflag) {
  YAKL_SCOPE( t                , :: t);
  YAKL_SCOPE( tabs             , :: tabs);
  YAKL_SCOPE( qcl              , :: qcl);
  YAKL_SCOPE( qci              , :: qci);
  YAKL_SCOPE( qv               , :: qv);
  YAKL_SCOPE( u                , :: u);
  YAKL_SCOPE( v                , :: v);
  YAKL_SCOPE( t0               , :: t0);
  YAKL_SCOPE( tabs0            , :: tabs0);
  YAKL_SCOPE( u0               , :: u0);
  YAKL_SCOPE( v0               , :: v0);
  YAKL_SCOPE( q0               , :: q0);
  YAKL_SCOPE( crm_accel_uv     , :: crm_accel_uv);
  YAKL_SCOPE( micro_field      , :: micro_field);
  YAKL_SCOPE( ncrms            , :: ncrms);
  YAKL_SCOPE( crm_accel_factor , :: crm_accel_factor);
  YAKL_SCOPE( microphysics_scheme, :: microphysics_scheme );

  real ttend_threshold = 2.0; // reduced from 5 to 2 after P3 implementation
  real tmin = 50.0;  // should never get below 50K in crm, following UP-CAM implementation
  int idx_qv = index_water_vapor;

  real2d ubaccel("ubaccel", nzm, ncrms);
  real2d vbaccel("vbaccel", nzm, ncrms);
  real2d tbaccel("tbaccel", nzm, ncrms);
  real2d qbaccel("qbaccel", nzm, ncrms);
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
    qbaccel(k,icrm) = 0.0;
    if (crm_accel_uv) {
      ubaccel(k,icrm) = 0.0;
      vbaccel(k,icrm) = 0.0;
    }
  });

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_DEVICE_LAMBDA (int k, int j, int i, int icrm) {
    // calculate tendency * dtn
    real q_tmp;
    real t_tmp;
    if (is_same_str(microphysics_scheme, "sam1mom")==0){
      t_tmp = t(k,j+offy_s,i+offx_s,icrm);
      q_tmp = qv(k,j,i,icrm)+qcl(k,j,i,icrm)+qci(k,j,i,icrm);
    }
    if (is_same_str(microphysics_scheme, "p3"     )==0){
      t_tmp = tabs(k,j,i,icrm);
      q_tmp = qv(k,j,i,icrm)+qcl(k,j,i,icrm)+qci(k,j,i,icrm);
    }
    yakl::atomicAdd( tbaccel(k,icrm) , t_tmp * crm_accel_coef );
    yakl::atomicAdd( qbaccel(k,icrm) , q_tmp * crm_accel_coef );
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
    real t0_tmp;
    real q0_tmp;
    if (is_same_str(microphysics_scheme, "sam1mom")==0){
      t0_tmp = t0(k,icrm);
      q0_tmp = q0(k,icrm);
    }
    if (is_same_str(microphysics_scheme, "p3"     )==0){
      t0_tmp = tabs0(k,icrm);
      q0_tmp = q0(k,icrm);
    }
    ttend_acc(k,icrm) = tbaccel(k,icrm) - t0_tmp;
    qtend_acc(k,icrm) = qbaccel(k,icrm) - q0_tmp;
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
    int new_nstop = nstop + (nstop - nstep + 1)*crm_accel_factor;
    std::cout << "accelerate_crm: mean-state acceleration not applied this step"<<std::endl;;
    std::cout << "accelerate_crm: nstop increased from "<<nstop<<" to "<<round(new_nstop)<<std::endl;;
    nstop = new_nstop; // only can happen once
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
    micro_field(0,k,j+offy_s,i+offx_s,icrm) =
        micro_field(0,k,j+offy_s,i+offx_s,icrm) + crm_accel_factor * qtend_acc(k,icrm);
  });

  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //!! Fix negative micro and re-adjust among separate water species
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
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_DEVICE_LAMBDA (int k, int j, int i, int icrm) {
    real tmp_q;
    if (is_same_str(microphysics_scheme, "sam1mom") == 0) {
      tmp_q = micro_field(0,k,j+offy_s,i+offx_s,icrm);
    }
    if (is_same_str(microphysics_scheme, "p3") == 0) {
      tmp_q = micro_field(idx_qv,k,j+offy_s,i+offx_s,icrm)
             +micro_field(idx_qc,k,j+offy_s,i+offx_s,icrm);
             +micro_field(idx_qi,k,j+offy_s,i+offx_s,icrm);
    }
    // real tmp_q = micro_field(0,k,j+offy_s,i+offx_s,icrm);
    if (tmp_q < 0.0) {
      yakl::atomicAdd( qneg(k,icrm) , tmp_q );
    }
    else {
      yakl::atomicAdd( qpoz(k,icrm) , tmp_q );
    }
  });

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    real factor;
    if (qpoz(k,icrm) + qneg(k,icrm) <= 0.0) {

      // all vapor and cloud water is depleted
      if (is_same_str(microphysics_scheme, "sam1mom") == 0) {
        micro_field(0,k,j+offy_s,i+offx_s,icrm) = 0.0;
        qv (k,j,i,icrm) = 0.0;
        qcl(k,j,i,icrm) = 0.0;
        qci(k,j,i,icrm) = 0.0;
      }
      if (is_same_str(microphysics_scheme, "p3") == 0) {
        micro_field(idx_qv,k,j+offy_s,i+offx_s,icrm) = 0.0;
        micro_field(idx_qc,k,j+offy_s,i+offx_s,icrm) = 0.0;
        micro_field(idx_qi,k,j+offy_s,i+offx_s,icrm) = 0.0;
        qv(k,j,i,icrm)  = micro_field(idx_qv,k,j+offy_s,i+offx_s,icrm);
        qcl(k,j,i,icrm) = micro_field(idx_qc,k,j+offy_s,i+offx_s,icrm);
        qci(k,j,i,icrm) = micro_field(idx_qi,k,j+offy_s,i+offx_s,icrm);
      }

    } else {

      // Clip qt values at 0 and remove the negative excess in each layer
      // proportionally from the positive qt fields in the layer
      factor = 1.0 + qneg(k,icrm) / qpoz(k,icrm);
      micro_field(0,k,j+offy_s,i+offx_s,icrm) = max(0.0, micro_field(0,k,j+offy_s,i+offx_s,icrm) * factor);

      if (is_same_str(microphysics_scheme, "sam1mom") == 0) {
        // Partition micro_field == qv + qcl + qci following these rules:
        //    (1) attempt to satisfy purely by adjusting qv
        //    (2) adjust qcl and qci only if needed to ensure positivity
        if (micro_field(0,k,j+offy_s,i+offx_s,icrm) < 0.0) {
          qv (k,j,i,icrm) = 0.0;
          qcl(k,j,i,icrm) = 0.0;
          qci(k,j,i,icrm) = 0.0;
        } else {
          // deduce qv as residual between qt - qcl - qci
          real qt_res = micro_field(0,k,j+offy_s,i+offx_s,icrm) - qcl(k,j,i,icrm) - qci(k,j,i,icrm);
          qv(k,j,i,icrm) = max(0.0, qt_res);
          if (qt_res < 0.0) {
            // qv was clipped; need to reduce qcl and qci accordingly
            factor = 1.0 + qt_res / (qcl(k,j,i,icrm) + qci(k,j,i,icrm));
            qcl(k,j,i,icrm) = qcl(k,j,i,icrm) * factor;
            qci(k,j,i,icrm) = qci(k,j,i,icrm) * factor;
          }
        }
      }

      if (is_same_str(microphysics_scheme, "p3") == 0) {
        // check for locally negative water vapor and re-partition
        // total water (qv+qcl+qci) to avoid negative values
        if (micro_field(idx_qv,k,j+offy_s,i+offx_s,icrm) < 0.0) {
          real qv_res = micro_field(idx_qv,k,j+offy_s,i+offx_s,icrm);
          real qc_tmp = micro_field(idx_qc,k,j+offy_s,i+offx_s,icrm); 
          real qi_tmp = micro_field(idx_qi,k,j+offy_s,i+offx_s,icrm);
          factor = 1.0 + qv_res / (qc_tmp + qi_tmp);
          micro_field(idx_qv,k,j+offy_s,i+offx_s,icrm) = 0.0; 
          micro_field(idx_qc,k,j+offy_s,i+offx_s,icrm) = qc_tmp * factor; 
          micro_field(idx_qi,k,j+offy_s,i+offx_s,icrm) = qi_tmp * factor;
          qv(k,j,i,icrm)  = micro_field(idx_qv,k,j+offy_s,i+offx_s,icrm);
          qcl(k,j,i,icrm) = micro_field(idx_qc,k,j+offy_s,i+offx_s,icrm);
          qci(k,j,i,icrm) = micro_field(idx_qi,k,j+offy_s,i+offx_s,icrm);
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


