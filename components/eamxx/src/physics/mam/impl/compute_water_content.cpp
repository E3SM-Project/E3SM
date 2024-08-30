#include <mam4xx/mam4.hpp>

namespace scream::impl {

KOKKOS_INLINE_FUNCTION
void compute_water_content( Real *state_q,   // in
                           const Real *qqcw, // in,
                           Real qv,// in
                           Real temp,// in
                           Real pmid, // in
                           const Real num_k_a[mam4::AeroConfig::num_modes()],// in
                           Real dgncur_a[mam4::AeroConfig::num_modes()], // out
                           Real dgncur_awet[mam4::AeroConfig::num_modes()],// out
                           Real wetdens[mam4::AeroConfig::num_modes()],// out
                           Real qaerwat[mam4::AeroConfig::num_modes()]// out
                           ) {
  constexpr int num_modes    = mam4::AeroConfig::num_modes();
  constexpr int num_aero_ids = mam4::AeroConfig::num_aerosol_ids();
  constexpr int pcnst = mam4::pcnst;


  // get some information about aerosol species
  // FIXME: this isn't great!
  int numptr_amode[num_modes];
  int mam_idx[num_modes][ndrop::nspec_max];
  int mam_cnst_idx[num_modes][ndrop::nspec_max];
  constexpr int maxd_aspectype = mam4::water_uptake::maxd_aspectype;
  int nspec_amode[num_modes], lspectype_amode[maxd_aspectype][num_modes];
  Real specdens_amode[maxd_aspectype], spechygro[maxd_aspectype];
  int lmassptr_amode[ndrop::maxd_aspectype][num_modes];
  mam4::ndrop::get_e3sm_parameters(nspec_amode, lspectype_amode, lmassptr_amode,
                             numptr_amode, specdens_amode, spechygro, mam_idx,
                             mam_cnst_idx);

  Real inv_density[num_modes][AeroConfig::num_aerosol_ids()] = {};
  Real num2vol_ratio_min[num_modes] = {};
  Real num2vol_ratio_max[num_modes] = {};
  Real num2vol_ratio_max_nmodes[num_modes] = {};
  Real num2vol_ratio_min_nmodes[num_modes] = {};
  Real num2vol_ratio_nom_nmodes[num_modes] = {};
  Real dgnmin_nmodes[num_modes] = {};
  Real dgnmax_nmodes[num_modes] = {};
  Real dgnnom_nmodes[num_modes] = {};
  // outputs
  bool noxf_acc2ait[AeroConfig::num_aerosol_ids()] = {};
  int n_common_species_ait_accum = {};
  int ait_spec_in_acc[AeroConfig::num_aerosol_ids()] = {};
  int acc_spec_in_ait[AeroConfig::num_aerosol_ids()] = {};
  Real mean_std_dev_nmodes[num_modes];
  mam4::modal_aero_calcsize::init_calcsize(
      inv_density, num2vol_ratio_min, num2vol_ratio_max,
      num2vol_ratio_max_nmodes, num2vol_ratio_min_nmodes,
      num2vol_ratio_nom_nmodes, dgnmin_nmodes, dgnmax_nmodes, dgnnom_nmodes,
      mean_std_dev_nmodes,
      // outputs
      noxf_acc2ait, n_common_species_ait_accum, ait_spec_in_acc,
      acc_spec_in_ait);

  // extract aerosol tracers for this level into state_q, which is needed
  // for computing dry aerosol properties below
  // FIXME: we should eliminate this index translation stuff

  // compute the dry volume for each mode, and from it the current dry
  // geometric nominal particle diameter.
  // FIXME: We have to do some gymnastics here to set up the calls to
  // FIXME: calcsize. This could be improved.

  for(int imode = 0; imode < num_modes; ++imode) {
    const auto v2nmin = num2vol_ratio_min[imode];
    const auto v2nmax = num2vol_ratio_max[imode];
    const auto dgnmin = dgnmin_nmodes[imode];
    const auto dgnmax = dgnmax_nmodes[imode];
    const auto mean_std_dev = mean_std_dev_nmodes[imode];

    Real dryvol_i, dryvol_c = 0.0;  // interstitial and cloudborne dry volumes
    mam4::modal_aero_calcsize::compute_dry_volume(imode,       // in
                       state_q,     // in
                       qqcw,        // in
                       inv_density, // in
                       lmassptr_amode,
                       dryvol_i, // out
                       dryvol_c);

    // NOTE: there's some disagreement over whether vol2num should be called
    // NOTE: num2vol here, so I'm just adopting the nomenclature used by
    // NOTE: the following call to calcsize)
    Real num2vol_ratio_cur_i=0;
    // Make it non-negative
    auto num_i_k = num_k_a[imode] < 0 ? 0 : num_k_a[imode];
    calcsize::update_diameter_and_vol2num(dryvol_i, num_i_k, v2nmin, v2nmax,
                                        dgnmin, dgnmax, mean_std_dev,
                                        dgncur_a[imode], num2vol_ratio_cur_i);
  }

  // calculate dry aerosol properties
  Real hygro[num_modes], naer[num_modes], dryrad[num_modes], dryvol[num_modes],
      drymass[num_modes], rhcrystal[num_modes], rhdeliques[num_modes],
      specdens_1[num_modes];
  mam4::water_uptake::modal_aero_water_uptake_dryaer(
      nspec_amode, specdens_amode, spechygro, lspectype_amode, state_q,
      dgncur_a, hygro, naer, dryrad, dryvol, drymass, rhcrystal, rhdeliques,
      specdens_1);

  // calculate wet aerosol properties
  Real rh = mam4::conversions::relative_humidity_from_vapor_mixing_ratio(
      qv, temp, pmid);
  Real wetrad[num_modes], wetvol[num_modes], wtrvol[num_modes];
  mam4::water_uptake::modal_aero_water_uptake_wetaer(
      rhcrystal, rhdeliques, dgncur_a, dryrad, hygro, rh, naer, dryvol, wetrad,
      wetvol, wtrvol, dgncur_awet, qaerwat);
  mam4::water_uptake::modal_aero_water_uptake_wetdens(wetvol, wtrvol, drymass,
                                                      specdens_1, wetdens);
}

}  // namespace scream::impl
