/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "VerticalRemapManager.hpp"
#include "SimulationParams.hpp"
#include "Context.hpp"
#include "Elements.hpp"
#include "Tracers.hpp"
#include "HybridVCoord.hpp"
#include "HommexxEnums.hpp"
#include "RemapFunctor.hpp"
#include "PpmRemap.hpp"

namespace Homme {

struct VerticalRemapManager::Impl {
  Impl(const SimulationParams &params, const Elements &e, const Tracers &t,
       const HybridVCoord &h) {
    using namespace Remap;
    using namespace Remap::Ppm;
    if (params.remap_alg == RemapAlg::PPM_FIXED_PARABOLA) {
      if (params.rsplit != 0) {
        remapper = std::make_shared<RemapFunctor<
            true, PpmVertRemap<PpmFixedParabola>> >(
            params.qsize, e, t, h);
      } else {
        remapper = std::make_shared<RemapFunctor<
            false, PpmVertRemap<PpmFixedParabola>> >(
            params.qsize, e, t, h);
      }
    } else if (params.remap_alg == RemapAlg::PPM_FIXED_MEANS) {
      if (params.rsplit != 0) {
        remapper = std::make_shared<RemapFunctor<
            true, PpmVertRemap<PpmFixedMeans>> >(
            params.qsize, e, t, h);
      } else {
        remapper = std::make_shared<RemapFunctor<
            false, PpmVertRemap<PpmFixedMeans>> >(
            params.qsize, e, t, h);
      }
    } else if (params.remap_alg == RemapAlg::PPM_MIRRORED) {
      if (params.rsplit != 0) {
        remapper = std::make_shared<RemapFunctor<
            true, PpmVertRemap<PpmMirrored>> >(
            params.qsize, e, t, h);
      } else {
        remapper = std::make_shared<RemapFunctor<
            false, PpmVertRemap<PpmMirrored>> >(
            params.qsize, e, t, h);
      }
    } else {
      Errors::runtime_abort(
          "Error in VerticalRemapManager: unknown remap algorithm.\n",
          Errors::err_unknown_option);
    }
  }

  std::shared_ptr<Remap::Remapper> remapper;
};

VerticalRemapManager::VerticalRemapManager() {
  const auto &h = Context::singleton().get<HybridVCoord>();
  const auto &p = Context::singleton().get<SimulationParams>();
  const auto &e = Context::singleton().get<Elements>();
  const auto &t = Context::singleton().get<Tracers>();
  assert(p.params_set);
  p_.reset(new Impl(p, e, t, h));
}

void VerticalRemapManager::run_remap(int np1, int np1_qdp, double dt) const {
  assert(p_);
  assert(p_->remapper);
  p_->remapper->run_remap(np1, np1_qdp, dt);
}

int VerticalRemapManager::requested_buffer_size () const {
  assert (p_);
  assert (p_->remapper);

  return p_->remapper->requested_buffer_size();
}

void VerticalRemapManager::init_buffers(const FunctorsBuffersManager& fbm) {
  assert (p_);
  assert (p_->remapper);

  p_->remapper->init_buffers(fbm);
}

} // namespace Homme
