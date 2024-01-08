#include <mam4xx/mam4.hpp>

namespace scream::impl {

using mam4::utils::min_max_bound;

// performs gas phase chemistry calculations on a single level of a single
// atmospheric column
KOKKOS_INLINE_FUNCTION
void gas_phase_chemistry(Real zm, Real zi, Real phis, Real temp, Real pmid, Real pdel, Real dt,
                         const Real photo_rates[mam4::mo_photo::phtcnt], // in
                         const Real extfrc[mam4::gas_chemistry::extcnt], // in
                         Real q[mam4::gas_chemistry::gas_pcnst], // VMRs, inout
                         Real invariants[mam4::gas_chemistry::nfs]) { // out
  constexpr Real rga = 1.0/haero::Constants::gravity;
  constexpr Real m2km = 0.01; // converts m -> km

  // The following things are chemical mechanism dependent! See mam4xx/src/mam4xx/gas_chem_mechanism.hpp)
  constexpr int gas_pcnst = mam4::gas_chemistry::gas_pcnst; // number of gas phase species
  constexpr int rxntot = mam4::gas_chemistry::rxntot;       // number of chemical reactions
  constexpr int extcnt = mam4::gas_chemistry::extcnt;       // number of species with external forcing
  constexpr int nabscol = mam4::gas_chemistry::nabscol;     // number of "absorbing column densities"
  constexpr int indexm = 0;  // FIXME: index of total atm density in invariants array

  constexpr int phtcnt = mam4::mo_photo::phtcnt; // number of photolysis reactions

  constexpr int itermax = mam4::gas_chemistry::itermax;
  constexpr int clscnt4 = mam4::gas_chemistry::clscnt4;

  // NOTE: vvv these arrays were copied from mam4xx/gas_chem_mechanism.hpp vvv
  constexpr int permute_4[gas_pcnst] = {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                        10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                                        20, 21, 22, 23, 24, 25, 26, 27, 28, 29};
  constexpr int clsmap_4[gas_pcnst] = {1,  2,  3,  4,  5,  6,  7,  8,  9,  10,
                                       11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                                       21, 22, 23, 24, 25, 26, 27, 28, 29, 30};

  // These indices for species are fixed by the chemical mechanism
  // std::string solsym[] = {"O3", "H2O2", "H2SO4", "SO2", "DMS", "SOAG",
  //                         "so4_a1", "pom_a1", "soa_a1", "bc_a1", "dst_a1",
  //                         "ncl_a1", "mom_a1", "num_a1", "so4_a2", "soa_a2",
  //                         "ncl_a2", "mom_a2", "num_a2", "dst_a3", "ncl_a3",
  //                         "so4_a3", "bc_a3", "pom_a3", "soa_a3", "mom_a3",
  //                         "num_a3", "pom_a4", "bc_a4", "mom_a4", "num_a4"};
  constexpr int ndx_h2so4 = 2;
  constexpr int o3_ndx = 0;
  // std::string extfrc_list[] = {"SO2", "so4_a1", "so4_a2", "pom_a4", "bc_a4",
  //                              "num_a1", "num_a2", "num_a3", "num_a4", "SOAG"};
  constexpr int synoz_ndx = -1;

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
  Real h2ovmr = q[0];
  // setinv(invariants, temp, h2ovmr, q, pmid); FIXME: not ported yet

  // ... set rates for "tabular" and user specified reactions
  Real reaction_rates[rxntot];
  mam4::gas_chemistry::setrxt(reaction_rates, temp);

  // set reaction rates based on chemical invariants
  // (indices (ndxes?) are taken from mam4 validation data and translated from
  // 1-based indices to 0-based indices)
  int usr_HO2_HO2_ndx = 1, usr_DMS_OH_ndx = 5,
      usr_SO2_OH_ndx = 3, inv_h2o_ndx = 3;
  mam4::gas_chemistry::usrrxt(reaction_rates, temp, invariants, invariants[indexm],
                              usr_HO2_HO2_ndx, usr_DMS_OH_ndx,
                              usr_SO2_OH_ndx, inv_h2o_ndx);
  mam4::gas_chemistry::adjrxt(reaction_rates, invariants, invariants[indexm]);

  //===================================
  // Photolysis rates at time = t(n+1)
  //===================================

  // compute the rate of change from forcing
  Real extfrc_rates[extcnt]; // [1/cm^3/s]
  for (int mm = 0; mm < extcnt; ++mm) {
    if (mm != synoz_ndx) {
      extfrc_rates[mm] = extfrc[mm] / invariants[indexm];
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

  // initialize error tolerances
  Real epsilon[clscnt4];
  mam4::gas_chemistry::imp_slv_inti(epsilon);

  // solve chemical system implicitly
  Real prod_out[clscnt4], loss_out[clscnt4];
  mam4::gas_chemistry::imp_sol(q, reaction_rates, het_rates, extfrc_rates, dt,
    permute_4, clsmap_4, factor, epsilon, prod_out, loss_out);

  // save h2so4 change by gas phase chem (for later new particle nucleation)
  if (ndx_h2so4 > 0) {
    del_h2so4_gasprod = q[ndx_h2so4] - del_h2so4_gasprod;
  }
}

} // namespace scream::impl
