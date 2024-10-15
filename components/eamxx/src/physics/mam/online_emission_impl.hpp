#ifndef ONLINE_EMISSION_IMPL_HPP
#define ONLINE_EMISSION_IMPL_HPP

namespace scream::mam_coupling {
template <typename S, typename D>
void onlineEmissions<S, D>::init_from_input_file(const ekat::ParameterList &params) {
  const int nspec = online_emis_data.nspec;
  const int ncols = online_emis_data.ncols;
  using ExeSpace = typename KT::ExeSpace;
  using ESU = ekat::ExeSpaceUtils<ExeSpace>;
  const auto policy = ESU::get_default_team_policy(ncols, nspec);

  start_timer("EAMxx::onlineEmiss::init_onlineEmiss_data_from_input_file");
  // 1. Read from input file
  // FIXME: currently reading a single placeholder scalar--should be
  //        ncols-sized array
  for (int ispec = 0; ispec < nspec; ++ispec) {
    Real init_cond_val =
        params.get<Real>(online_emis_data.root_IC_str + online_emis_data.spec_names[ispec]);
    // FIXME: is this overkill--i.e., would a mirror/deep_copy make more sense?
    Kokkos::parallel_for(
        policy, KOKKOS_LAMBDA(const MemberType &team) {
          const int jcol = team.league_rank(); // column index
          online_emis_data.flux_data(ispec, jcol) = init_cond_val;
        });
  }
  stop_timer("EAMxx::onlineEmiss::init_onlineEmiss_data_from_input_file");

} // end init_from_input_file()

template <typename S, typename D>
void onlineEmissions<S, D>::transfer_to_cflux(
    const onlineEmissData &data, const std::map<std::string, int> idx_map,
    view_2d &fluxes) {
  // FIXME: move out to onlineEmissions struct?
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
              std::cout << "=========================================" << "\n";
              std::cout << "data.cfluxes(jcol, s_idx) = " << data.cfluxes(jcol, s_idx) << "\n";
              std::cout << "=========================================" << "\n";
            });
      });

  Kokkos::deep_copy(fluxes, data.cfluxes);
}
} // namespace scream::mam_coupling

#endif // ONLINE_EMISSION_IMPL_HPP
