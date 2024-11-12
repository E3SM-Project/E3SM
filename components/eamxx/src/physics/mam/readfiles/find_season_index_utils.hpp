#ifndef EAMXX_MAM_FIND_SEASON_INDEX_UTILS
#define EAMXX_MAM_FIND_SEASON_INDEX_UTILS

#include <ekat/kokkos/ekat_kokkos_utils.hpp>
#include "share/io/scorpio_input.hpp"
#include "share/io/scream_scorpio_interface.hpp"
#include <mam4xx/mam4.hpp>

namespace scream::mam_coupling {

using ExeSpace = typename KT::ExeSpace;
using ESU      = ekat::ExeSpaceUtils<ExeSpace>;

// views for single- and multi-column data
using view_1d       = typename KT::template view_1d<Real>;
using const_view_1d       = typename KT::template view_1d<const Real>;

using view_int_1d       = typename KT::template view_1d<int>;
using view_int_2d       = typename KT::template view_2d<int>;
using view_int_3d       = typename KT::template view_3d<int>;

using view_1d_host = typename KT::view_1d<Real>::HostMirror;
using view_3d_host = typename KT::view_3d<Real>::HostMirror;
using view_int_3d_host = typename KT::view_3d<int>::HostMirror;


inline void find_season_index_reader(const std::string& season_wes_file,
                                     const const_view_1d& clat,
                                     view_int_2d &index_season_lai)
{
  const int  plon= clat.extent(0);
  scorpio::register_file(season_wes_file, scorpio::Read);

  const int nlat_lai = scorpio::get_dimlen(season_wes_file, "lat");
  const int npft_lai = scorpio::get_dimlen(season_wes_file, "pft");
  view_1d_host lat_lai_host("lat_lai_host", nlat_lai);
  view_int_3d_host wk_lai_host("wk_lai_host", nlat_lai, npft_lai, 12);

  scorpio::read_var(season_wes_file, "lat", lat_lai_host.data());
  scorpio::read_var(season_wes_file, "season_wes", wk_lai_host.data());
  scorpio::release_file(season_wes_file);

  view_int_3d wk_lai("wk_lai", nlat_lai, npft_lai, 12);

  view_1d lat_lai("lat_lai", nlat_lai);
  Kokkos::deep_copy(lat_lai, lat_lai_host);
  Kokkos::deep_copy(wk_lai, wk_lai_host);

  // output
  index_season_lai=view_int_2d("index_season_lai", plon,12);
  const auto policy = ESU::get_default_team_policy(plon, 1);
  Kokkos::parallel_for(
        policy, KOKKOS_LAMBDA(const Team &team) {
          const int j = team.league_rank();
          const auto index_season_lai_at_j = ekat::subview(index_season_lai, j);
          mam4::mo_drydep::find_season_index(clat(j), lat_lai, nlat_lai, wk_lai,
                                       index_season_lai_at_j);
  });
}

}  // namespace scream::mam_coupling
#endif  //
