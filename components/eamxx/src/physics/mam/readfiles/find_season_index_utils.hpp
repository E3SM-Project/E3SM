#ifndef EAMXX_MAM_FIND_SEASON_INDEX_UTILS
#define EAMXX_MAM_FIND_SEASON_INDEX_UTILS

#include <ekat/kokkos/ekat_kokkos_utils.hpp>
#include <mam4xx/mam4.hpp>

#include "share/io/scorpio_input.hpp"
#include "share/io/scream_scorpio_interface.hpp"

namespace scream::mam_coupling {

// views for single- and multi-column data

using const_view_1d = typename KT::template view_1d<const Real>;
using view_int_2d   = typename KT::template view_2d<int>;

using view_1d_host     = typename KT::view_1d<Real>::HostMirror;
using view_int_3d_host = typename KT::view_3d<int>::HostMirror;
using view_int_2d_host = typename KT::view_2d<int>::HostMirror;

/**
 * @brief Reads the season index from the given file and computes the season
 * indices based on latitudes.
 *
 * @param[in] season_wes_file The path to the season_wes.nc file.
 * @param[in] clat A 1D view of latitude values in degrees.
 * @param[out] index_season_lai A 2D view to store the computed season indices.
 * Note that indices are in C++ (starting from zero).
 */

inline void find_season_index_reader(const std::string &season_wes_file,
                                     const const_view_1d &clat,
                                     view_int_2d &index_season_lai) {
  const int plon = clat.extent(0);
  scorpio::register_file(season_wes_file, scorpio::Read);

  const int nlat_lai = scorpio::get_dimlen(season_wes_file, "lat");
  const int npft_lai = scorpio::get_dimlen(season_wes_file, "pft");

  view_1d_host lat_lai("lat_lai", nlat_lai);
  view_int_2d_host wk_lai_temp("wk_lai", npft_lai, nlat_lai);
  view_int_3d_host wk_lai("wk_lai", nlat_lai, npft_lai, 12);

  scorpio::read_var(season_wes_file, "lat", lat_lai.data());

  Kokkos::MDRangePolicy<Kokkos::HostSpace::execution_space, Kokkos::Rank<2>>
      policy_wk_lai({0, 0}, {nlat_lai, npft_lai});

  // loop over time to get all 12 instantence of season_wes
  for(int itime = 0; itime < 12; ++itime) {
    scorpio::read_var(season_wes_file, "season_wes", wk_lai_temp.data(), itime);
    // copy data from wk_lai_temp to wk_lai.
    // NOTE: season_wes has different layout that wk_lai
    Kokkos::parallel_for("copy_to_wk_lai", policy_wk_lai,
                         [&](const int j, const int k) {
                           wk_lai(j, k, itime) = wk_lai_temp(k, j);
                         });
    Kokkos::fence();
  }
  scorpio::release_file(season_wes_file);

  index_season_lai = view_int_2d("index_season_lai", plon, 12);
  const view_int_2d_host index_season_lai_host =
      Kokkos::create_mirror_view(index_season_lai);

  const view_1d_host clat_host = Kokkos::create_mirror_view(clat);
  Kokkos::deep_copy(clat_host, clat);

  // Computation is performed on the host
  mam4::mo_drydep::find_season_index(clat_host, lat_lai, nlat_lai, wk_lai,
                                     index_season_lai_host);
  Kokkos::deep_copy(index_season_lai, index_season_lai_host);
}
}  // namespace scream::mam_coupling
#endif  //
