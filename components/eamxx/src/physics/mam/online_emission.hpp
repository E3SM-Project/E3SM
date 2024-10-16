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
  void init_from_input_file(const ekat::ParameterList &m_params);
  void transfer_to_cflux(const onlineEmissData &data,
                                const std::map<std::string, int> idx_map,
                                view_2d &fluxes);
}; // struct onlineEmissions
} // namespace scream::mam_coupling
#endif // ONLINE_EMISSION_HPP

#include "online_emission_impl.hpp"
