#include <mam4xx/mam4.hpp>
#include <mam_coupling.hpp>

namespace scream::impl {

// performs gas phase chemistry calculations on a single level of a single
// atmospheric column
KOKKOS_INLINE_FUNCTION
void gas_phase_chemistry(Real zm, Real zi, Real phis, Real q[mam_coupling::gas_pcnst()]) {
  constexpr Real rga = 1.0/haero::gravity;
  constexpr Real m2km = 0.01; // converts m -> km
                              //
  // FIXME: The following things are chemical mechanism dependent! See mam4xx/src/mam4xx/gas_chem_mechanism.hpp)
  constexpr int gas_pcnst = mam_coupling::gas_pcnst();
  constexpr int rxntot = 7; // number of chemical reactions
  constexpr int extcnt = 9; // number of species with external forcing

  // fetch the zenith angle (not its cosine!) in degrees for this column.
  // FIXME: RRTMGP does this in Fortran--shall we also for now?
  Real zenith = 0.0; // [deg]

  // xform geopotential height from m to km and pressure from Pa to mb
  Real zsurf = rga * phis;
  Real zintr = m2km * zi;
  Real zmid = m2km * (zm + zsurf);
  Real zint = m2km * (zi + zsurf);

  // ... set atmosphere mean mass
  set_mean_mass(mbar);

  // ... xform water vapor from mmr to vmr and set upper bndy values
  Real qh2o = q[0];
  Real h2ovmr;
  h2o_to_vmr(qh2o, h2ovmr, mbar);

  // ... set the "invariants"
  setinv(invariants, tfld, h2ovmr, vmr, pmid, ncol, lchnk, pbuf);

  // ... set the column densities at the upper boundary
  // FIXME: This is level-independent, but we can probably get away with
  // FIXME: calling it at all levels
  set_ub_col(col_delta, vmr, invariants, pdel, ncol, lchnk);

  // ... set rates for "tabular" and user specified reactions
  Real rxt_rates[3];
  setrxt(reaction_rates, tfld);

  // compute the relative humidity
  // FIXME: We have a better way of doing this in EAMxx
  Real satq;
  qsat(tfld, pmid, satv, satq);
  Real relhum = min_max_bound(0, 1, 0.622 * h2ovmr/satq);

  usrrxt(reaction_rates, tfld, invariants, invariants[indexm]);
  adjrxt(reaction_rates, invariants, invariants[indexm]);

  //===================================
  // Photolysis rates at time = t(n+1)
  //===================================

  // ... set the column density
  // FIXME: Again, level-independent
  setcol(col_delta, col_dens);

  // ... calculate the photodissociation rates
  // FIXME: This looks like an EAM Fortran call.
  Real esfact = 1.0;
  shr_orb_decl(calday, eccen, mvelpp, lambm0, obliqr, delta, esfact);

  // ... compute the extraneous frcing at time = t(n+1)
  setext(extfrc, lchnk, ncol, zintr);

  for (int mm = 0; mm < extcnt; ++mm) {
    if (mm != synoz_ndx) {
      extfrc[mm] /= invariants[indexm];
    }
  }

  // ... Form the washout rates
  sethet(het_rates, pmid, zmid, phis, tfld, cmfdqr, prain, nevapr, delt,
    invariants[indexm], vmr, ncol, lchnk);

  ltrop_sol = 0; // apply solver to all levels

  // save h2so4 before gas phase chem (for later new particle nucleation)
  del_h2so4_gasprod = q[ndx_h2so4];

  //===========================
  // Class solution algorithms
  //===========================

  // ... solve for "Implicit" species
  // FIXME: need to figure out arguments to ported C++ code
  imp_sol(q, reaction_rates, het_rates, extfrc, delt,
    invariants[indexm], ncol, lchnk, ltrop_sol);

  /* I don't think we need to worry about this, do we?
  if (convproc_do_aer) {
    vmr2mmr(vmr, mmr_new, mbar, ncol);
    mmr_new(:ncol,:,:) = 0.5_r8*( mmr(:ncol,:,:)+mmr_new(:ncol,:,:) );
    //RCE - mmr_new = average of mmr values before and after imp_sol
    het_diags(het_rates(:ncol,:,:), mmr_new(:ncol,:,:), pdel(:ncol,:), lchnk, ncol );
  }
  */

  // save h2so4 change by gas phase chem (for later new particle nucleation)
  if (ndx_h2so4 > 0) {
    del_h2so4_gasprod = q[ndx_h2so4] - del_h2so4_gasprod;
  }
}

} // namespace scream::impl
