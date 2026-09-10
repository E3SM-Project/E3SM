#ifndef GW_GW_CONVECT_INIT_IMPL_HPP
#define GW_GW_CONVECT_INIT_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU
#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"

#include <ekat_assert.hpp>

namespace scream {
namespace gw {

/*
 * Implementation of gw gw_convect_init. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
typename Functions<S,D>::template view_3d<Real>
Functions<S,D>::gw_convect_read_mfcc_table(
  // Inputs
  const std::string& gw_drag_file,
  const int& npgw)
{
  scorpio::register_file(gw_drag_file, scorpio::FileMode::Read);
  const int PS_dim_size = scorpio::get_dimlen(gw_drag_file, "PS"); // Phase Speed [m/s]
  const int MW_dim_size = scorpio::get_dimlen(gw_drag_file, "MW"); // Mean Wind in Heating [m/s]
  const int HD_dim_size = scorpio::get_dimlen(gw_drag_file, "HD"); // Heating Depth [km]

  // mfcc is stored on disk with dims ordered (PS,MW,HD) -- PS slowest, HD
  // fastest -- not (HD,MW,PS). Guard against that assumption silently
  // breaking if a different table is ever supplied.
  EKAT_REQUIRE_MSG(scorpio::get_var(gw_drag_file,"mfcc").dim_names() ==
                    std::vector<std::string>({"PS","MW","HD"}),
                    "Error! Unexpected dimension order for 'mfcc' in " << gw_drag_file <<
                    ". Expected (PS,MW,HD).\n");

  // Read the raw, full-spectrum table into a plain buffer matching that
  // on-disk order (avoids relying on Kokkos's device-dependent default
  // view layout for the raw fill).
  std::vector<Real> mfcc_raw(static_cast<size_t>(PS_dim_size)*MW_dim_size*HD_dim_size);
  scorpio::read_var(gw_drag_file, "mfcc", mfcc_raw.data());
  scorpio::release_file(gw_drag_file);

  // The file's native phase-speed spectrum spans [-ngwv_file,ngwv_file].
  // Extract only the subset matching the runtime pgwv, centered on zero --
  // mirrors gw_drag.F90's NF90_GET_VAR(..., start=[1,1,ngwv_file-pgwv+1]).
  const int ngwv_file = (PS_dim_size - 1) / 2;
  const int pgwv      = (npgw - 1) / 2;
  EKAT_REQUIRE_MSG(pgwv <= ngwv_file,
    "Error! Requested pgwv (" << pgwv << ") exceeds the max wavenumber "
    "available in gw_drag_file (" << ngwv_file << ").\n");
  const int ps_offset = ngwv_file - pgwv;

  // Transpose (PS,MW,HD) -> (HD,MW,PS) and window PS down to npgw, staging
  // into the device view mfcc(hdepth_index, uh_index, wave_index).
  view_3d<Real> mfcc("mfcc", HD_dim_size, MW_dim_size, npgw);
  auto mfcc_h = Kokkos::create_mirror_view(mfcc);
  for (int ps = 0; ps < npgw; ++ps) {
    const int ps_file = ps + ps_offset;
    for (int mw = 0; mw < MW_dim_size; ++mw) {
      for (int hd = 0; hd < HD_dim_size; ++hd) {
        mfcc_h(hd,mw,ps) = mfcc_raw[(static_cast<size_t>(ps_file)*MW_dim_size + mw)*HD_dim_size + hd];
      }
    }
  }
  Kokkos::deep_copy(mfcc, mfcc_h);

  return mfcc;
}

template<typename S, typename D>
void Functions<S,D>::gw_convect_init(
  // Inputs
  const Real& plev_src_wind,
  const uview_3d<const Real>& mfcc_in)
{
  static bool s_convect_init_constructed = false;
  if (!s_convect_init_constructed) {
    new (&s_convect_init) GwConvectInit();
    s_convect_init_constructed = true;
  }

  // Just set k_src_wind to pver. We don't have access to pref_edge
  s_convect_init.k_src_wind = s_common_init.pver - 1;

  // First dimension is maxh.
  s_convect_init.maxh   = mfcc_in.extent(0);

  // Second dimension is -maxuh:maxuh (size 2*maxuh+1).
  s_convect_init.maxuh = (mfcc_in.extent(1) - 1) / 2;

  s_convect_init.mfcc = view_3d<Real>("mfcc", mfcc_in.extent(0), mfcc_in.extent(1), mfcc_in.extent(2));
  Kokkos::deep_copy(s_convect_init.mfcc, mfcc_in);
}

/*
The version above was written for unit tests. The version below
is better suited for running the GW schemes in the full model
*/

template<typename S, typename D>
void Functions<S,D>::gw_convect_init(
  // Inputs
  ekat::ParameterList& params,
  const uview_1d<const Real>& pref_int,
  const uview_3d<const Real>& mfcc_in)
{
  static bool s_convect_init_constructed = false;
  if (!s_convect_init_constructed) {
    new (&s_convect_init) GwConvectInit();
    s_convect_init_constructed = true;
  }

  s_convect_init.gw_convect_eff             = params.get<Real>("gw_convect_eff", s_convect_init.gw_convect_eff);
  s_convect_init.gw_convect_hcf             = params.get<Real>("gw_convect_hcf", s_convect_init.gw_convect_hcf);
  s_convect_init.gw_convect_hdepth_scale    = params.get<Real>("gw_convect_hdepth_scale", s_convect_init.gw_convect_hdepth_scale);
  s_convect_init.gw_convect_hdepth_min      = params.get<Real>("gw_convect_hdepth_min", s_convect_init.gw_convect_hdepth_min);
  s_convect_init.gw_convect_storm_speed_min = params.get<Real>("gw_convect_storm_speed_min", s_convect_init.gw_convect_storm_speed_min);
  s_convect_init.gw_convect_plev_src_wind   = params.get<Real>("gw_convect_plev_src_wind", s_convect_init.gw_convect_plev_src_wind);
  s_convect_init.use_gw_convect_old         = params.get<bool>("use_gw_convect_old", s_convect_init.use_gw_convect_old);
  s_convect_init.mfcc = view_3d<Real>("mfcc", mfcc_in.extent(0), mfcc_in.extent(1), mfcc_in.extent(2));
  Kokkos::deep_copy(s_convect_init.mfcc, mfcc_in);

  // First dimension is maxh
  s_convect_init.maxh   = mfcc_in.extent(0);

  // Second dimension is -maxuh:maxuh (size 2*maxuh+1)
  s_convect_init.maxuh = (mfcc_in.extent(1) - 1) / 2;

  // set index for source wind (small serial scan; do on host)
  auto pref_int_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pref_int);
  s_convect_init.k_src_wind = -1;
  for (int i = 0; i < static_cast<int>(pref_int_h.size()); ++i) {
    if (pref_int_h[i] < s_convect_init.gw_convect_plev_src_wind) {
      s_convect_init.k_src_wind = i;
    }
  }
  EKAT_REQUIRE_MSG(s_convect_init.k_src_wind >= 1,
    "Error! No reference interface pressure level found below plev_src_wind.\n");

}

} // namespace gw
} // namespace scream

#endif
