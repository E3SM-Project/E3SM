#include <mam4xx/mam4.hpp>

namespace scream::impl {

using mam4::utils::min_max_bound;

// performs gas phase chemistry calculations on a single level of a single
// atmospheric column
KOKKOS_INLINE_FUNCTION
void gas_phase_chemistry(Real zm, Real zi, Real phis, Real temp, Real pmid, Real pdel, Real dt,
                         Real col_dens[mam4::gas_chemistry::nabscol],
                         Real photo_rates[mam4::mo_photo::phtcnt],
                         Real q[mam4::gas_chemistry::gas_pcnst]) { // VMRs
  constexpr Real rga = 1.0/haero::Constants::gravity;
  constexpr Real mwdry = 1.0/haero::Constants::molec_weight_dry_air;
  constexpr Real m2km = 0.01; // converts m -> km

  // FIXME: The following things are chemical mechanism dependent! See mam4xx/src/mam4xx/gas_chem_mechanism.hpp)
  constexpr int gas_pcnst = mam4::gas_chemistry::gas_pcnst; // number of gas phase species
  constexpr int rxntot = mam4::gas_chemistry::rxntot;       // number of chemical reactions
  constexpr int extcnt = mam4::gas_chemistry::extcnt;       // number of species with external forcing
  constexpr int nfs = mam4::gas_chemistry::nfs;             // number of "fixed species"
  constexpr int nabscol = mam4::gas_chemistry::nabscol;     // number of "absorbing column densities"
  constexpr int indexm = 0;  // index of total atm density in invariants array

  constexpr int phtcnt = mam4::mo_photo::phtcnt; // number of photolysis reactions

  constexpr int itermax = mam4::gas_chemistry::itermax;
  constexpr int clscnt4 = mam4::gas_chemistry::clscnt4;
  int permute_4[gas_pcnst]; // = mam4::gas_chemistry::permute_4; FIXME
  int clsmap_4[gas_pcnst];  // = mam4::gas_chemistry::clsmap_4; FIXME

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

  // ... compute the column's invariants
  Real invariants[nfs];
  Real h2ovmr = q[0];
  // setinv(invariants, temp, h2ovmr, q, pmid); FIXME: nœt ported yet

  // ... set rates for "tabular" and user specified reactions
  Real reaction_rates[rxntot];
  mam4::gas_chemistry::setrxt(reaction_rates, temp);

  // set reaction rates based on chemical invariants
  // FIXME: figure out these indices
  int usr_HO2_HO2_ndx = -1, usr_DMS_OH_ndx = -1,
      usr_SO2_OH_ndx = -1, inv_h2o_ndx = -1;
  mam4::gas_chemistry::usrrxt(reaction_rates, temp, invariants, invariants[indexm],
                              usr_HO2_HO2_ndx, usr_DMS_OH_ndx,
                              usr_SO2_OH_ndx, inv_h2o_ndx);
  mam4::gas_chemistry::adjrxt(reaction_rates, invariants, invariants[indexm]);

  //===================================
  // Photolysis rates at time = t(n+1)
  //===================================

  // ... compute the extraneous frcing at time = t(n+1)
  Real extfrc[extcnt];
  // mam4::mo_setext::setext(extfrc, zintr); // FIXME: not ported yet

  for (int mm = 0; mm < extcnt; ++mm) {
    if (mm != synoz_ndx) {
      extfrc[mm] /= invariants[indexm];
    }
  }

  // ... Form the washout rates
  Real het_rates[gas_pcnst];
  // FIXME: not ported yet
  //sethet(het_rates, pmid, zmid, phis, temp, cmfdqr, prain, nevapr, delt,
  //       invariants[indexm], q);

  int ltrop_sol = 0; // apply solver to all levels

  // save h2so4 before gas phase chem (for later new particle nucleation)
  Real del_h2so4_gasprod = q[ndx_h2so4];

  //===========================
  // Class solution algorithms
  //===========================

  // copy photolysis rates into reaction_rates (assumes photolysis rates come first)
  for (int i = 0; i < phtcnt; ++i) {
    reaction_rates[i] = photo_rates[i];
  }

  // ... solve for "Implicit" species
  bool factor[itermax];
  for (int i = 0; i < itermax; ++i) {
    factor[i] = true;
  }
  Real epsilon[clscnt4];
  mam4::gas_chemistry::imp_slv_inti(epsilon);
  Real prod_out[clscnt4], loss_out[clscnt4];
  mam4::gas_chemistry::imp_sol(q, reaction_rates, het_rates, extfrc, dt,
    permute_4, clsmap_4, factor, epsilon, prod_out, loss_out);

  /* I don't think we need to worry about this, do we?
  if (convproc_do_aer) {
    vmr2mmr(q, mmr_new, mbar, ncol);
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
