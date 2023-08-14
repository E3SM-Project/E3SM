
#include "timeloop.h"

#ifdef MMF_RAD_SORT
void update_sort_idx() {
  YAKL_SCOPE( t                        , :: t );
  YAKL_SCOPE( qv                       , :: qv  );
  YAKL_SCOPE( qcl                      , :: qcl );
  YAKL_SCOPE( qci                      , :: qci );
  YAKL_SCOPE( sort_q                   , :: sort_q );
  YAKL_SCOPE( sort_i                   , :: sort_i );
  YAKL_SCOPE( sort_j                   , :: sort_j );

  int num_idx = ny*nx;

  // initialize arrays used for sorting
  parallel_for( SimpleBounds<3>(ny,nx,ncrms) , YAKL_LAMBDA (int j, int i, int icrm) {
    int idx = i*ny + j;
    sort_q(idx,icrm)=0.; sort_i(idx,icrm)=i; sort_j(idx,icrm)=j;
  });
  // set water condensate+vapor quantity to use for sorting
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int idx = i*ny + j;
    // scale the water species so the sorting favors qcl, qci, and qv in that order
    real tmp_q = qcl(k,j,i,icrm) + 1e-3*qci(k,j,i,icrm) + 1e-6*qv(k,j,i,icrm);
    yakl::atomicAdd( sort_q(idx,icrm) , tmp_q );
  });

  // // print some debugging stuff
  // parallel_for( SimpleBounds<3>(ny,nx,ncrms) , YAKL_LAMBDA ( int j, int i, int icrm) {
  //   int idx = i*ny + j;
  //   int ii = sort_i(idx,icrm);
  //   std::cout << "WHDEBUG - icrm: "     << icrm
  //                      << "   pre-sort "
  //                      << "   i: "      << i
  //                      << "   ii: "     << ii
  //                      << "   sort_q: " << sort_q(idx,icrm)
  //                      << std::endl;
  // });

  // us simple bubble sort to order columns and save the original column indices
  parallel_for( SimpleBounds<1>(ncrms) , YAKL_LAMBDA (int icrm) {
    for(int idx1 = 0; idx1<num_idx; idx1++) {
      for(int idx2 = idx1+1; idx2<num_idx; idx2++) {
        if( sort_q(idx2,icrm) < sort_q(idx1,icrm) ) {
          real tq = sort_q(idx1,icrm);
          real ti = sort_i(idx1,icrm);
          real tj = sort_j(idx1,icrm);
          sort_q(idx1,icrm) = sort_q(idx2,icrm);
          sort_i(idx1,icrm) = sort_i(idx2,icrm);
          sort_j(idx1,icrm) = sort_j(idx2,icrm);
          sort_q(idx2,icrm) = tq;
          sort_i(idx2,icrm) = ti;
          sort_j(idx2,icrm) = tj;
        }
      }
    }
  });

  // // print some debugging stuff
  // parallel_for( SimpleBounds<3>(ny,nx,ncrms) , YAKL_LAMBDA ( int j, int i, int icrm) {
  //   int idx = i*ny + j;
  //   int ii = sort_i(idx,icrm);
  //   std::cout << "WHDEBUG - icrm: "     << icrm
  //                      << "   post-sort "
  //                      << "   i: "      << i
  //                      << "   ii: "     << ii
  //                      << "   sort_q: " << sort_q(idx,icrm)
  //                      << std::endl;
  // });
}
#endif

void timeloop() {
  YAKL_SCOPE( crm_output_subcycle_factor , :: crm_output_subcycle_factor );
  YAKL_SCOPE( t                        , :: t );
  YAKL_SCOPE( qv                       , :: qv  );
  YAKL_SCOPE( qcl                      , :: qcl );
  YAKL_SCOPE( qci                      , :: qci );
  YAKL_SCOPE( crm_rad_qrad             , :: crm_rad_qrad );
  YAKL_SCOPE( dtn                      , :: dtn );
  YAKL_SCOPE( ncrms                    , :: ncrms );
  YAKL_SCOPE( na                       , :: na );
  YAKL_SCOPE( dt3                      , :: dt3 );
  YAKL_SCOPE( use_VT                   , :: use_VT );
  YAKL_SCOPE( use_ESMT                 , :: use_ESMT );
#ifdef MMF_RAD_SORT
  YAKL_SCOPE( CF3D                     , :: CF3D );
  YAKL_SCOPE( crm_rad_temperature      , :: crm_rad_temperature );
  YAKL_SCOPE( tabs                     , :: tabs );
  YAKL_SCOPE( pres                     , :: pres );
  YAKL_SCOPE( crm_rad_qv               , :: crm_rad_qv );
  YAKL_SCOPE( crm_rad_qc               , :: crm_rad_qc );
  YAKL_SCOPE( crm_rad_qi               , :: crm_rad_qi );
  YAKL_SCOPE( crm_rad_cld              , :: crm_rad_cld );
  YAKL_SCOPE( crm_clear_rh             , :: crm_clear_rh );
  YAKL_SCOPE( crm_clear_rh_cnt         , :: crm_clear_rh_cnt );
  // YAKL_SCOPE( sort_q                   , :: sort_q );
  YAKL_SCOPE( sort_i                   , :: sort_i );
  YAKL_SCOPE( sort_j                   , :: sort_j );

  // start by populating rad sort indices
  update_sort_idx();
#endif


  nstep = 0;

  do {
    nstep = nstep + 1;

    //------------------------------------------------------------------
    //  Check if the dynamical time step should be decreased
    //  to handle the cases when the flow being locally linearly unstable
    //------------------------------------------------------------------
    kurant();

    for(int icyc=1; icyc<=ncycle; icyc++) {
      icycle = icyc;
      dtn = dt/ncycle;
      parallel_for( 1 , YAKL_LAMBDA ( int i ) {
        dt3(na-1) = dtn;
      });
      dtfactor = dtn/dt;

      parallel_for( ncrms , YAKL_LAMBDA (int icrm) {
        crm_output_subcycle_factor(icrm) = crm_output_subcycle_factor(icrm)+1;
      });

      //---------------------------------------------
      //    the Adams-Bashforth scheme in time
      abcoefs();

      //---------------------------------------------
      //    initialize stuff:
      zero();

      //-----------------------------------------------------------
      //       Buoyancy term:
      buoyancy();

      //-----------------------------------------------------------
      // variance transport forcing
      if (use_VT) {
        VT_diagnose();
        VT_forcing();
      }

      //------------------------------------------------------------
      //       Large-scale and surface forcing:
      forcing();

      // Apply radiative tendency
      // for (int k=0; k<nzm; k++) {
      //   for (int j=0; j<ny; j++) {
      //     for (int i=0; i<nx; i++) {
      //       for (int icrm=0; icrm<ncrms; icrm++) {
#ifdef MMF_RAD_SORT
      parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
        int i_rad = i / (nx/crm_nx_rad);
        int j_rad = j / (ny/crm_ny_rad);
        int idx = i*ny + j;
        int ii = sort_i(idx,icrm) + offx_s;
        int jj = sort_j(idx,icrm) + offy_s;
        t(k,jj,ii,icrm) = t(k,jj,ii,icrm) + crm_rad_qrad(k,j_rad,i_rad,icrm)*dtn;
      });
#else
      parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
        int i_rad = i / (nx/crm_nx_rad);
        int j_rad = j / (ny/crm_ny_rad);
        int ii = i+offx_s;
        int jj = j+offy_s;
        t(k,jj,ii,icrm) = t(k,jj,ii,icrm) + crm_rad_qrad(k,j_rad,i_rad,icrm)*dtn;
      });
#endif
      //----------------------------------------------------------
      //    suppress turbulence near the upper boundary (spange):
      if (dodamping) { 
        damping();
      }

      //---------------------------------------------------------
      //   Ice fall-out
      if (docloud) { 
        ice_fall();
      }

      //----------------------------------------------------------
      //     Update scalar boundaries after large-scale processes:
      boundaries(3);

      //---------------------------------------------------------
      //     Update boundaries for velocities:
      boundaries(0);

      //-----------------------------------------------
      //     surface fluxes:
      if (dosurface) {
        crmsurface(bflx);
      }

      //-----------------------------------------------------------
      //  SGS physics:
      if (dosgs) {
        sgs_proc();
      }

      //----------------------------------------------------------
      //     Fill boundaries for SGS diagnostic fields:
      boundaries(4);

      //-----------------------------------------------
      //       advection of momentum:
      advect_mom();

      //----------------------------------------------------------
      //  SGS effects on momentum:
      if (dosgs) { 
        sgs_mom();
      }

      //----------------------------------------------------------
      //  Explicit scalar momentum transport scheme (ESMT)
      if (use_ESMT) {
        scalar_momentum_tend();
      }

      //-----------------------------------------------------------
      //       Coriolis force:
      if (docoriolis) {
        coriolis();
      }

      //---------------------------------------------------------
      //       compute rhs of the Poisson equation and solve it for pressure.
      pressure();

      //---------------------------------------------------------
      //       find velocity field at n+1/2 timestep needed for advection of scalars:
      //  Note that at the end of the call, the velocities are in nondimensional form.
      adams();

      //----------------------------------------------------------
      //     Update boundaries for all prognostic scalar fields for advection:
      boundaries(2);

      //---------------------------------------------------------
      //      advection of scalars :
      advect_all_scalars();

      //-----------------------------------------------------------
      //    Convert velocity back from nondimensional form:
      uvw();

      //----------------------------------------------------------
      //     Update boundaries for scalars to prepare for SGS effects:
      boundaries(3);

      //---------------------------------------------------------
      //      SGS effects on scalars :
      if (dosgs) { 
        sgs_scalars();
      }

      //-----------------------------------------------------------
      //       Calculate PGF for scalar momentum tendency

      //-----------------------------------------------------------
      //       Cloud condensation/evaporation and precipitation processes:
      if (docloud || dosmoke) {
        micro_proc();
      }

      //-----------------------------------------------------------
      //       Apply mean-state acceleration
      if (use_crm_accel && !crm_accel_ceaseflag) {
        // Use Jones-Bretherton-Pritchard methodology to accelerate
        // CRM horizontal mean evolution artificially.
        accelerate_crm(nstep, nstop, crm_accel_ceaseflag);
      }

      //-----------------------------------------------------------
      //    Compute diagnostics fields:
      diagnose();

      //----------------------------------------------------------
      // Rotate the dynamic tendency arrays for Adams-bashforth scheme:

      int nn=na;
      na=nc;
      nc=nb;
      nb=nn;
    } // icycle

    post_icycle();

#ifdef MMF_RAD_SORT
  update_sort_idx();

  // average output rad data
  // int num_idx = ny*nx;
  // parallel_for( SimpleBounds<2>(num_idx,ncrms) , YAKL_LAMBDA (int idx, int icrm) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int idx = i*ny + j;
    int ii = sort_i(idx,icrm);
    int jj = sort_j(idx,icrm);
    int i_rad = i / (nx/crm_nx_rad);
    int j_rad = j / (ny/crm_ny_rad);
    // for (int k=0; k<nzm; k++) {
      yakl::atomicAdd(crm_rad_temperature(k,j_rad,i_rad,icrm) , tabs(k,jj,ii,icrm));
      real tmp = max(0.0,qv(k,jj,ii,icrm));
      yakl::atomicAdd(crm_rad_qv(k,j_rad,i_rad,icrm) , tmp);
      yakl::atomicAdd(crm_rad_qc(k,j_rad,i_rad,icrm) , qcl(k,jj,ii,icrm));
      yakl::atomicAdd(crm_rad_qi(k,j_rad,i_rad,icrm) , qci(k,jj,ii,icrm));
      if (qcl(k,jj,ii,icrm) + qci(k,jj,ii,icrm) > 0) {
        yakl::atomicAdd(crm_rad_cld(k,j_rad,i_rad,icrm) , CF3D(k,jj,ii,icrm));
      } else {
        real qsat_tmp;
        qsatw_crm(tabs(k,jj,ii,icrm),pres(k,icrm),qsat_tmp);
        real rh_tmp = qv(k,jj,ii,icrm)/qsat_tmp;
        yakl::atomicAdd(crm_clear_rh(k,icrm) , rh_tmp);
        yakl::atomicAdd(crm_clear_rh_cnt(k,icrm),1);
      }
    // }
  });
// #else
//   parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
//     // Reduced radiation method allows for fewer radiation calculations
//     // by collecting statistics and doing radiation over column groups
//     int i_rad = i / (nx/crm_nx_rad);
//     int j_rad = j / (ny/crm_ny_rad);
//     real qsat_tmp;
//     real rh_tmp;

//     yakl::atomicAdd(crm_rad_temperature(k,j_rad,i_rad,icrm) , tabs(k,j,i,icrm));
//     real tmp = max(0.0,qv(k,j,i,icrm));
//     yakl::atomicAdd(crm_rad_qv(k,j_rad,i_rad,icrm) , tmp);
//     yakl::atomicAdd(crm_rad_qc(k,j_rad,i_rad,icrm) , qcl(k,j,i,icrm));
//     yakl::atomicAdd(crm_rad_qi(k,j_rad,i_rad,icrm) , qci(k,j,i,icrm));
//     if (qcl(k,j,i,icrm) + qci(k,j,i,icrm) > 0) {
//       yakl::atomicAdd(crm_rad_cld(k,j_rad,i_rad,icrm) , CF3D(k,j,i,icrm));
//     } else {
//       qsatw_crm(tabs(k,j,i,icrm),pres(k,icrm),qsat_tmp);
//       rh_tmp = qv(k,j,i,icrm)/qsat_tmp;
//       yakl::atomicAdd(crm_clear_rh(k,icrm) , rh_tmp);
//       yakl::atomicAdd(crm_clear_rh_cnt(k,icrm),1);
//     }
//   });
#endif

  } while (nstep < nstop);

}
