#ifndef ONLINE_EMISSION_HPP
#define ONLINE_EMISSION_HPP

#include "share/util/scream_timing.hpp"

namespace scream::mam_coupling {
template <typename ScalarType, typename DeviceType> struct onlineEmissions {
  using Device = DeviceType;
  using KT = KokkosTypes<Device>;
  using MemberType = typename KT::MemberType;
  static constexpr int pcnst = mam4::aero_model::pcnst;

  struct onlineEmissData {
    // Basic spatial dimensions of the data
    int ncols;
    view_2d flux_data;
    // local copy of main fluxes array
    view_2d cfluxes;
    // FIXME: read this from elsewhere? input?
    const std::vector<std::string> spec_names = {
        "ncl_a1", "ncl_a2", "ncl_a3", "mom_a1", "mom_a2", "mom_a4",
        "num_a1", "num_a2", "num_a3", "num_a4", "dst_a1", "dst_a3"};
    // FIXME: change this when/if the above is dynamically-determined
    int nspec = spec_names.size();
    const std::string root_IC_str = "online_emis_IC_";

    onlineEmissData() = default;
    onlineEmissData(const int ncol_, const int nspec_)
        : ncols(ncol_), nspec(nspec_) {
      init(ncols, nspec, true);
    }

    onlineEmissData(const int ncol_) : ncols(ncol_) {
      init(ncols, nspec, true);
    }

    // overloads of init() in case npsec is not hard-coded, or if you want to
    // control allocation via bool flag
    void init(const int ncol_, const int nspec_, const bool allocate_) {
      ncols = ncol_;
      nspec = nspec_;
      if (allocate_) {
        flux_data = view_2d("onlineEmissData", nspec, ncols);
        cfluxes = view_2d("onlineEmisLocalCflux", ncols, pcnst);
      }
    } // onlineEmissData init
    void init(const int ncol_, const bool allocate_) {
      ncols = ncol_;
      if (allocate_) {
        flux_data = view_2d("onlineEmissData", nspec, ncols);
        cfluxes = view_2d("onlineEmisLocalCflux", ncols, pcnst);
      }
    } // onlineEmissData init
    void init(const int ncol_) {
      ncols = ncol_;
      flux_data = view_2d("onlineEmissData", nspec, ncols);
      cfluxes = view_2d("onlineEmisLocalCflux", ncols, pcnst);
    } // onlineEmissData init
  };  // onlineEmissData

  onlineEmissData online_emis_data;

  // ---------------------------------------------------------------------------
  // Online emissions routines
  // ---------------------------------------------------------------------------
  void init_from_input_file(const ekat::ParameterList &params) {
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
}; // struct onlineEmissions
} // namespace scream::mam_coupling
#endif // ONLINE_EMISSION_HPP
