#ifndef SRF_EMISSION_HPP
#define SRF_EMISSION_HPP

namespace scream::mam_coupling {
namespace {

template <typename ScalarType, typename DeviceType>
struct srfEmissFunctions {
  /* -------------------------------------------------------------------------------------------
   */
  // Surface emissions routines

  static std::shared_ptr<AbstractRemapper> create_horiz_remapper(
      const std::shared_ptr<const AbstractGrid> &model_grid,
      const std::string &spa_data_file, const std::string &map_file,
      const bool use_iop = false);
};  // struct srfEmissFunctions
}  // namespace
}  // namespace scream::mam_coupling
#endif  // SRF_EMISSION_HPP

#include "srf_emission_impl.hpp"