#include <mam4xx/mam4.hpp>

namespace scream::impl {

KOKKOS_INLINE_FUNCTION
void compute_water_content(const mam4::Prognostics &progs, int k,
                           Real qv, Real temp, Real pmid,
                           Real dgncur_a[mam4::AeroConfig::num_modes()],
                           Real dgncur_awet[mam4::AeroConfig::num_modes()],
                           Real wetdens[mam4::AeroConfig::num_modes()],
                           Real qaerwat[mam4::AeroConfig::num_modes()]) {
  constexpr int num_modes = mam4::AeroConfig::num_modes();
  constexpr int num_aero_ids = mam4::AeroConfig::num_aerosol_ids();

  // get some information about aerosol species
  // FIXME: this isn't great!
  constexpr int maxd_aspectype = mam4::water_uptake::maxd_aspectype;
  int nspec_amode[num_modes], lspectype_amode[maxd_aspectype][num_modes];
  Real specdens_amode[maxd_aspectype], spechygro[maxd_aspectype];
  mam4::water_uptake::get_e3sm_parameters(nspec_amode, lspectype_amode,
      specdens_amode, spechygro);

  // extract aerosol tracers for this level into state_q, which is needed
  // for computing dry aerosol properties below
  // FIXME: we should eliminate this index translation stuff
  constexpr int nvars = mam4::water_uptake::nvars;
  Real state_q[nvars]; // aerosol tracers for level k
  for (int imode = 0; imode < num_modes; ++imode) {
    int la, lc; // interstitial and cloudborne indices within state_q

    // number mixing ratios
    mam4::convproc::assign_la_lc(imode, -1, la, lc);
    state_q[la] = progs.n_mode_i[imode](k);
    state_q[lc] = progs.n_mode_c[imode](k);
    // aerosol mass mixing ratios
    for (int iaero = 0; iaero < num_aero_ids; ++iaero) {
      mam4::convproc::assign_la_lc(imode, iaero, la, lc);
      auto mode = static_cast<mam4::ModeIndex>(imode);
      auto aero = static_cast<mam4::AeroId>(iaero);
      int ispec = mam4::aerosol_index_for_mode(mode, aero);
      if (ispec != -1) {
        state_q[la] = progs.q_aero_i[imode][ispec](k);
        state_q[lc] = progs.q_aero_c[imode][ispec](k);
      }
    }
  }

  // compute the dry volume for each mode, and from it the current dry
  // geometric nominal particle diameter.
  // FIXME: We have to do some gymnastics here to set up the calls to
  // FIXME: calcsize. This could be improved.
  Real inv_densities[num_modes][num_aero_ids] = {};
  for (int imode = 0; imode < num_modes; ++imode) {
    const int n_spec = mam4::num_species_mode(imode);
    for (int ispec = 0; ispec < n_spec; ++ispec) {
      const int iaer = static_cast<int>(mam4::mode_aero_species(imode, ispec));
      const Real density = mam4::aero_species(iaer).density;
      inv_densities[imode][ispec] = 1.0 / density;
    }
  }
  for (int imode = 0; imode < num_modes; ++imode) {
    Real dryvol_i, dryvol_c; // interstitial and cloudborne dry volumes
    mam4::calcsize::compute_dry_volume_k(k, imode, inv_densities, progs,
        dryvol_i, dryvol_c);

    // NOTE: there's some disagreement over whether vol2num should be called
    // NOTE: num2vol here, so I'm just adopting the nomenclature used by
    // NOTE: the following call to calcsize)
    const mam4::Mode& mode = mam4::modes(imode);
    Real vol2num_min = 1.0/mam4::conversions::mean_particle_volume_from_diameter(
        mode.max_diameter, mode.mean_std_dev);
    Real vol2num_max = 1.0/mam4::conversions::mean_particle_volume_from_diameter(
        mode.min_diameter, mode.mean_std_dev);
    Real vol2num;
    mam4::calcsize::update_diameter_and_vol2num(dryvol_i,
        progs.n_mode_i[imode](k), vol2num_min, vol2num_max,
        mode.min_diameter, mode.max_diameter, mode.mean_std_dev,
        dgncur_a[imode], vol2num);
  }

  // calculate dry aerosol properties
  Real hygro[num_modes], naer[num_modes], dryrad[num_modes],
  dryvol[num_modes], drymass[num_modes],
  rhcrystal[num_modes], rhdeliques[num_modes], specdens_1[num_modes];
  mam4::water_uptake::modal_aero_water_uptake_dryaer(nspec_amode, specdens_amode,
      spechygro, lspectype_amode, state_q, dgncur_a, hygro,
      naer, dryrad, dryvol, drymass, rhcrystal, rhdeliques, specdens_1);

  // calculate wet aerosol properties
  Real rh = mam4::conversions::relative_humidity_from_vapor_mixing_ratio(qv, temp, pmid);
  Real wetrad[num_modes], wetvol[num_modes], wtrvol[num_modes];
  mam4::water_uptake::modal_aero_water_uptake_wetaer(rhcrystal, rhdeliques, dgncur_a,
      dryrad, hygro, rh, naer, dryvol, wetrad, wetvol, wtrvol, dgncur_awet,
      qaerwat);
  mam4::water_uptake::modal_aero_water_uptake_wetdens(wetvol, wtrvol,
      drymass, specdens_1, wetdens);
}

} // namespace scream::impl
