
#include "timeloop.h"

void timeloop() {
  YAKL_SCOPE( crm_output_subcycle_factor , :: crm_output_subcycle_factor );
  YAKL_SCOPE( t                        , :: t );
  YAKL_SCOPE( crm_rad_qrad             , :: crm_rad_qrad );
  YAKL_SCOPE( dtn                      , :: dtn );
  YAKL_SCOPE( ncrms                    , :: ncrms );
  YAKL_SCOPE( na                       , :: na );
  YAKL_SCOPE( dt3                      , :: dt3 );
  YAKL_SCOPE( use_VT                   , :: use_VT );
  YAKL_SCOPE( use_ESMT                 , :: use_ESMT );

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
      parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
        int i_rad = i / (nx/crm_nx_rad);
        int j_rad = j / (ny/crm_ny_rad);
        t(k,j+offy_s,i+offx_s,icrm) = t(k,j+offy_s,i+offx_s,icrm) + crm_rad_qrad(k,j_rad,i_rad,icrm)*dtn;
      });

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

  } while (nstep < nstop);

}
