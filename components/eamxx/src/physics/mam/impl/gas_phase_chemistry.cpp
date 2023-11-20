#include <mam4xx/mam4.hpp>

namespace scream::impl {

using mam4::utils::min_max_bound;

// performs gas phase chemistry calculations on a single level of a single
// atmospheric column
KOKKOS_INLINE_FUNCTION
void gas_phase_chemistry(Real zm, Real zi, Real phis, Real temp, Real pmid, Real pdel, Real dt,
                         const view_1d& photo_rates, // per-column photolysis rates
                         Real q[mam4::gas_chemistry::gas_pcnst]) {
  constexpr Real rga = 1.0/haero::Constants::gravity;
  constexpr Real mwdry = 1.0/haero::Constants::molec_weight_dry_air;
  constexpr Real m2km = 0.01; // converts m -> km

  // FIXME: The following things are chemical mechanism dependent! See mam4xx/src/mam4xx/gas_chem_mechanism.hpp)
  constexpr int gas_pcnst = mam4::gas_chemistry::gas_pcnst;
  constexpr int rxntot = 7;  // number of chemical reactions
  constexpr int extcnt = 9;  // number of species with external forcing
  constexpr int nfs = 8;     // number of "fixed species"
  constexpr int nabscol = 2; // number of "absorbing column densities"
  constexpr int indexm = 0;  // index of total atm density in invariants array

  constexpr int ndx_h2so4 = 0; // FIXME: get_spc_ndx('H2SO4')
  constexpr int o3_ndx = 0;    // FIXME: get_spc_ndx('O3')
  constexpr int synoz_ndx = 0; // FIXME: get_extfrc_ndx('SYNOZ')

  // fetch the zenith angle (not its cosine!) in degrees for this column.
  // FIXME: For now, we fix the zenith angle. At length, we need to compute it
  // FIXME: from EAMxx's current set of orbital parameters, which requires some
  // FIXME: conversation with the EAMxx team.
  Real zen_angle = 0.0; // [deg]

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
  // setinv(invariants, temp, h2ovmr, vmr, pmid); FIXME

  // ... set the column densities at the upper boundary
  // FIXME: This is level-independent, but we can probably get away with
  // FIXME: calling it at all levels
  Real col_delta[nabscol];
  // set_ub_col(col_delta, vmr, invariants, pdel); FIXME

  // ... set rates for "tabular" and user specified reactions
  Real reaction_rates[rxntot];
  mam4::gas_chemistry::setrxt(reaction_rates, temp);

  // compute the relative humidity
  Real relhum = mam4::conversions::relative_humidity_from_vapor_mixing_ratio(qh2o, temp, pmid);

  // FIXME: we need to figure out the arguments for these functions
  //mam4::gas_chemistry::usrrxt(reaction_rates, temp, invariants, invariants[indexm]);
  //mam4::gas_chemistry::adjrxt(reaction_rates, invariants, invariants[indexm]);

  //===================================
  // Photolysis rates at time = t(n+1)
  //===================================

  // ... set the column density
  Real col_dens[nabscol];
  //setcol(col_delta, col_dens); // FIXME

  // ... compute the extraneous frcing at time = t(n+1)
  Real extfrc[extcnt];
  // FIXME: Same thing: can we do this for each single level?
  // mam4::mo_setext::setext(extfrc, zintr);

  for (int mm = 0; mm < extcnt; ++mm) {
    if (mm != synoz_ndx) {
      extfrc[mm] /= invariants[indexm];
    }
  }

  // ... Form the washout rates
  Real het_rates[gas_pcnst];
  //sethet(het_rates, pmid, zmid, phis, temp, cmfdqr, prain, nevapr, delt,
  //       invariants[indexm], vmr);

  int ltrop_sol = 0; // apply solver to all levels

  // save h2so4 before gas phase chem (for later new particle nucleation)
  Real del_h2so4_gasprod = q[ndx_h2so4];

  //===========================
  // Class solution algorithms
  //===========================

  // copy photolysis rates to reaction_rates

  // ... solve for "Implicit" species
  bool factor[mam4::gas_chemistry::itermax];
  for (int i = 0; i < mam4::gas_chemistry::itermax; ++i) {
    factor[i] = true;
  }
  Real epsilon[mam4::gas_chemistry::clscnt4];
  imp_slv_inti(epsilon);
  Real prod_out[mam4::gas_chemistry::clscnt4], loss_out[mam4::gas_chemistry::clscnt4];
  mam4::gas_chemistry::imp_sol(vmr, reaction_rates, het_rates, extfrc, dt,
    mam4::gas_chemistry::permute_4, mam4::gas_chemistry::clsmap_4, factor,
    epsilon, prod_out, loss_out);

  /* I don't think we need to worry about this, do we?
  if (convproc_do_aer) {
    vmr2mmr(vmr, mmr_new, mbar, ncol);
    mmr_new(:ncol,:,:) = 0.5_r8*( mmr(:ncol,:,:)+mmr_new(:ncol,:,:) );
    //RCE - mmr_new = average of mmr values before and after imp_sol
    het_diags(het_rates(:ncol,:,:), mmr_new(:ncol,:,:), pdel(:ncol,:));
  }
  */

  // save h2so4 change by gas phase chem (for later new particle nucleation)
  if (ndx_h2so4 > 0) {
    del_h2so4_gasprod = q[ndx_h2so4] - del_h2so4_gasprod;
  }
}

} // namespace scream::impl
