#ifndef P3_WARM_RAIN_EMULATOR_IMPL_HPP
#define P3_WARM_RAIN_EMULATOR_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

/*
 * Implementation of the merge of emulated warm-rain rates into P3's.
 * Clients should NOT #include this file, #include p3_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::warm_rain_emulator_merge(
  const Pack& nc, const Pack& nr, const Scalar& inv_dt,
  const Scalar& autoconv_fallback_scale, const bool& cloud_self_collection,
  const Pack& emu_qc2qr_autoconv_tend, const Pack& emu_nc2nr_autoconv_tend,
  const Pack& emu_ncautr, const Pack& emu_nc_selfcollect_tend,
  const Pack& emu_qc2qr_accret_tend, const Pack& emu_nc_accret_tend,
  const Pack& emu_nr_selfcollect_tend, const Pack& emu_use_cloud,
  const Pack& emu_use_rain, Pack& qc2qr_autoconv_tend, Pack& nc2nr_autoconv_tend,
  Pack& ncautr, Pack& nc_selfcollect_tend, Pack& qc2qr_accret_tend,
  Pack& nc_accret_tend, Pack& nr_selfcollect_tend, const Mask& context)
{
  const auto use_cloud = emu_use_cloud > sp(0.5) && context;
  const auto use_rain  = emu_use_rain  > sp(0.5) && context;
  const auto fallback  = !use_cloud && context;

  // Where the emulator does not apply: P3, with scaled autoconversion
  qc2qr_autoconv_tend.set(fallback, qc2qr_autoconv_tend * autoconv_fallback_scale);
  nc2nr_autoconv_tend.set(fallback, nc2nr_autoconv_tend * autoconv_fallback_scale);
  ncautr.set(fallback, ncautr * autoconv_fallback_scale);

  if (use_cloud.any()) {
    // These are sinks of qc and nc and sources of qr and nr: never negative
    qc2qr_autoconv_tend.set(use_cloud, max(emu_qc2qr_autoconv_tend, 0));
    nc2nr_autoconv_tend.set(use_cloud, max(emu_nc2nr_autoconv_tend, 0));
    ncautr.set(use_cloud, max(emu_ncautr, 0));
    qc2qr_accret_tend.set(use_cloud, max(emu_qc2qr_accret_tend, 0));
    nc_accret_tend.set(use_cloud, max(emu_nc_accret_tend, 0));

    // Droplet self-collection is a sink of nc (<= 0, added to nc).
    // nc_conservation takes nc + nc_selfcollect_tend*dt as the source of the
    // other nc sinks, so it must not remove more than there is.
    const Pack nc_sc = cloud_self_collection ? emu_nc_selfcollect_tend : Pack(0);
    nc_selfcollect_tend.set(use_cloud, max(min(nc_sc, 0), -max(nc, 0) * inv_dt));
  }

  if (use_rain.any()) {
    // Rain self-collection is a sink of nr (>= 0, subtracted from nr), which
    // must not remove more than there is, with what autoconversion adds.
    // nr_conservation counts nc2nr_autoconv_tend, not ncautr, as that source,
    // so its own check is not enough.
    nr_selfcollect_tend.set(use_rain, min(max(emu_nr_selfcollect_tend, 0),
                                          max(nr, 0) * inv_dt + ncautr));
  }
}

} // namespace p3
} // namespace scream

#endif // P3_WARM_RAIN_EMULATOR_IMPL_HPP
