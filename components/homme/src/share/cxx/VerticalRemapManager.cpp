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
       const HybridVCoord &h)
    : m_hvcoord(h)
    , m_params(params)
    , m_elements(e)
    , m_tracers(t)
  {
    setup_remmapper();
  }

  Impl(const SimulationParams &params, const HybridVCoord &h)
    : m_hvcoord(h)
    , m_params(params)
  {}

  void setup_remmapper ()
  {
    using namespace Remap;
    using namespace Remap::Ppm;
    if (m_params.remap_alg == RemapAlg::PPM_FIXED_PARABOLA) {
      if (m_params.rsplit != 0) {
        remapper = std::make_shared<RemapFunctor<
            true, PpmVertRemap<PpmFixedParabola>> >(
            m_params.qsize, m_elements, m_tracers, m_hvcoord);
      } else {
        remapper = std::make_shared<RemapFunctor<
            false, PpmVertRemap<PpmFixedParabola>> >(
            m_params.qsize, m_elements, m_tracers, m_hvcoord);
      }
    } else if (m_params.remap_alg == RemapAlg::PPM_FIXED_MEANS) {
      if (m_params.rsplit != 0) {
        remapper = std::make_shared<RemapFunctor<
            true, PpmVertRemap<PpmFixedMeans>> >(
            m_params.qsize, m_elements, m_tracers, m_hvcoord);
      } else {
        remapper = std::make_shared<RemapFunctor<
            false, PpmVertRemap<PpmFixedMeans>> >(
            m_params.qsize, m_elements, m_tracers, m_hvcoord);
      }
    } else if (m_params.remap_alg == RemapAlg::PPM_MIRRORED) {
      if (m_params.rsplit != 0) {
        remapper = std::make_shared<RemapFunctor<
            true, PpmVertRemap<PpmMirrored>> >(
            m_params.qsize, m_elements, m_tracers, m_hvcoord);
      } else {
        remapper = std::make_shared<RemapFunctor<
            false, PpmVertRemap<PpmMirrored>> >(
            m_params.qsize, m_elements, m_tracers, m_hvcoord);
      }
    } else {
      Errors::runtime_abort(
          "Error in VerticalRemapManager: unknown remap algorithm.\n",
          Errors::err_unknown_option);
    }
  }

  void setup (const Elements &e, const Tracers &t)
  {
    m_elements = e;
    m_tracers  = t;

    setup_remmapper();
  }

  std::shared_ptr<Remap::Remapper> remapper;

  HybridVCoord     m_hvcoord;
  SimulationParams m_params;
  Elements         m_elements;
  Tracers          m_tracers;

};

VerticalRemapManager::VerticalRemapManager()
  : is_setup(true)
{
  const auto &h = Context::singleton().get<HybridVCoord>();
  const auto &p = Context::singleton().get<SimulationParams>();
  const auto &e = Context::singleton().get<Elements>();
  m_num_elems = e.num_elems();
  const auto &t = Context::singleton().get<Tracers>();
  assert(p.params_set);
  p_.reset(new Impl(p, e, t, h));
}

VerticalRemapManager::VerticalRemapManager(const int num_elems)
  : m_num_elems(num_elems)
  , is_setup(false)
{
  const auto &h = Context::singleton().get<HybridVCoord>();
  const auto &p = Context::singleton().get<SimulationParams>();
  assert(p.params_set);
  p_.reset(new Impl(p, h));
}

void VerticalRemapManager::setup ()
{
  assert(!is_setup);

  const auto &e = Context::singleton().get<Elements>();
  assert(m_num_elems == e.num_elems()); // Sanity check

  const auto &t = Context::singleton().get<Tracers>();
  p_->setup(e,t);

  is_setup = true;
}

void VerticalRemapManager::run_remap(int np1, int np1_qdp, double dt) const {
  assert(is_setup);

  assert(p_);
  assert(p_->remapper);
  p_->remapper->run_remap(np1, np1_qdp, dt);
}

struct TempTagStruct  {};

int VerticalRemapManager::requested_buffer_size () const {
  assert (p_);
  assert (p_->remapper);

  if (is_setup) {
    return p_->remapper->requested_buffer_size();
  }
  // If the struct is not fully setup, we must manually compute
  // buffer request, since internally the Element functor is required.
  else {
#ifdef HOMMEXX_BFB_TESTING
    const bool process_nh_vars = true;
#else
    const bool process_nh_vars = !(p_->m_params.theta_hydrostatic_mode);
#endif
    if (p_->m_params.rsplit == 0 || !process_nh_vars) {
      return 0;
    } else {
      const int iters = 2;
      const int team_size = (iters*m_num_elems > 0 ? iters*m_num_elems : 1);
      TeamUtils<ExecSpace> dummy_tu(Homme::get_default_team_policy<ExecSpace,TempTagStruct>(team_size));

      using temp_type = ExecViewUnmanaged<Scalar*  [NP][NP][NUM_LEV  ]>;
      using phi_type  = ExecViewUnmanaged<Scalar*  [NP][NP][NUM_LEV_P]>;

      const int temp_size = temp_type::shmem_size(dummy_tu.get_num_ws_slots())/sizeof(Real);
      const int phi_size  = phi_type::shmem_size(dummy_tu.get_num_ws_slots())/sizeof(Real);
      return temp_size + phi_size;
    }
  }
}

void VerticalRemapManager::init_buffers(const FunctorsBuffersManager& fbm) {
  assert (p_);
  assert (p_->remapper);

  p_->remapper->init_buffers(fbm);
}

} // namespace Homme
