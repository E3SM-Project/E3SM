#ifndef ONLINE_EMISSION_IMPL_HPP
#define ONLINE_EMISSION_IMPL_HPP

namespace scream::mam_coupling {
template <typename S, typename D>
void onlineEmissions<S, D>::init_from_input_file(const ekat::ParameterList &params) {
  // FIXME: move this out to the onlineEmissions struct to avoid extra work
  // by doing it again below?
  const int nspec = online_emis_data.nspec;
  const int ncols = online_emis_data.ncols;
  using ExeSpace = typename KT::ExeSpace;
  using ESU = ekat::ExeSpaceUtils<ExeSpace>;
  const auto policy = ESU::get_default_team_policy(ncols, nspec);
  // Read from input file
  // FIXME: currently reading a single placeholder scalar--should be
  //        ncols-sized array when we know what the input data looks like
  for (int ispec = 0; ispec < nspec; ++ispec) {
    Real init_cond_val =
        params.get<Real>(online_emis_data.root_IC_str + online_emis_data.spec_names[ispec]);
    // TODO: is this overkill?--i.e., would a mirror/deep_copy make more sense?
    Kokkos::parallel_for(
        policy, KOKKOS_LAMBDA(const MemberType &team) {
          const int jcol = team.league_rank(); // column index
          online_emis_data.flux_data(ispec, jcol) = init_cond_val;
        });
  }
} // end init_from_input_file()

template <typename S, typename D>
void onlineEmissions<S, D>::transfer_to_cflux(
    const onlineEmissData &data, const std::map<std::string, int> idx_map,
    view_2d &fluxes) {
  const int nspec = data.nspec;
  const int ncols = data.ncols;
  using ExeSpace = typename KT::ExeSpace;
  using ESU = ekat::ExeSpaceUtils<ExeSpace>;
  const auto policy = ESU::get_default_team_policy(ncols, nspec);

  Kokkos::parallel_for(
      policy, KOKKOS_LAMBDA(const MemberType &team) {
        const int jcol = team.league_rank(); // column index
        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team, nspec), [&](const int ispec) {
              auto s_idx = idx_map.at(data.spec_names[ispec]);
              data.cfluxes(jcol, s_idx) = data.flux_data(ispec, jcol);
            });
      });
  Kokkos::deep_copy(fluxes, data.cfluxes);
}
} // namespace scream::mam_coupling

#endif // ONLINE_EMISSION_IMPL_HPP
