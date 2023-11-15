#include <mam4xx/gas_chem_mechanism.hpp>
#include <mam4xx/mam4.hpp>

namespace scream::impl {

using mam4::utils::min_max_bound;

// performs gas phase chemistry calculations on a single level of a single
// atmospheric column
KOKKOS_INLINE_FUNCTION
void gas_phase_chemistry(Real zm, Real zi, Real phis, Real temp, Real pmid, Real q[mam4::gas_chemistry::gas_pcnst]) {
  constexpr Real rga = 1.0/haero::Constants::gravity;
  constexpr Real mwdry = 1.0/haero::Constants::molec_weight_dry_air;
  constexpr Real m2km = 0.01; // converts m -> km

  // FIXME: The following things are chemical mechanism dependent! See mam4xx/src/mam4xx/gas_chem_mechanism.hpp)
  constexpr int gas_pcnst = mam4::gas_chemistry::gas_pcnst;
  constexpr int rxntot = 7; // number of chemical reactions
  constexpr int extcnt = 9; // number of species with external forcing
  constexpr int nfs = 8;    // number of "fixed species"

  // fetch the zenith angle (not its cosine!) in degrees for this column.
  // FIXME: For now, we fix the zenith angle. At length, we need to compute it
  // FIXME: from EAMxx's current set of orbital parameters, which requires some
  // FIXME: conversation with the EAMxx team.
  Real zenith = 0.0; // [deg]

  // xform geopotential height from m to km and pressure from Pa to mb
  Real zsurf = rga * phis;
  Real zintr = m2km * zi;
  Real zmid = m2km * (zm + zsurf);
  Real zint = m2km * (zi + zsurf);

  // ... map incoming mass mixing ratios to working array
  Real mmr[gas_pcnst] = {};
  // FIXME: come back and fix this

  // ... set atmosphere mean mass to the molecular weight of dry air
  Real mbar = mwdry;

  // ... Xform from mmr to vmr
  Real vmr[gas_pcnst];
  for (int i = 0; i < gas_pcnst; ++i) {
    vmr[i] = mam4::conversions::vmr_from_mmr(mmr[i], mbar);
  }

  // ... xform water vapor from mmr to vmr and set upper bndy values
  Real qh2o = q[0];
  Real h2ovmr = mam4::conversions::vmr_from_mmr(qh2o, mbar);

  // ... set the "invariants"
  Real invariants[nfs];
  setinv(invariants, temp, h2ovmr, vmr, pmid, ncol, lchnk, pbuf);

  // ... set the column densities at the upper boundary
  // FIXME: This is level-independent, but we can probably get away with
  // FIXME: calling it at all levels
  set_ub_col(col_delta, vmr, invariants, pdel, ncol, lchnk);

  // ... set rates for "tabular" and user specified reactions
  Real rxt_rates[3];
  mam4::gas_chemistry::setrxt(reaction_rates, temp);

  // compute the relative humidity
  // FIXME: We have a better way of doing this in EAMxx
  Real satq;
  qsat(temp, pmid, satv, satq);
  Real relhum = min_max_bound(0, 1, 0.622 * h2ovmr/satq);

  mam4::gas_chemistry::usrrxt(reaction_rates, temp, invariants, invariants[indexm]);
  mam4::gas_chemistry::adjrxt(reaction_rates, invariants, invariants[indexm]);

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
  mam4::mo_setext::setext(extfrc, lchnk, ncol, zintr);

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
  mam4::gas_chemistry::imp_sol(q, reaction_rates, het_rates, extfrc, delt,
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
