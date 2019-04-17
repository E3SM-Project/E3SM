#include "atmosphere_dynamics.hpp"
#include "scream_homme_interface.hpp"
#include "Context.hpp"
#include "mpi/MpiContext.hpp"
#include "Types.hpp"
#include "Elements.hpp"
#include "Tracers.hpp"
#include "Dimensions.hpp"
#include "CaarFunctor.hpp"
#include "EulerStepFunctor.hpp"
#include "HyperviscosityFunctor.hpp"

namespace scream
{

HommeDynamics::HommeDynamics (const Comm& comm,const ParameterList& /* params */)
 : m_dynamics_comm (comm)
{
  init_homme1(comm);

  constexpr int NUM_PHYSICAL_LEV = Homme::NUM_PHYSICAL_LEV;
  constexpr int NUM_TIME_LEVELS  = Homme::NUM_TIME_LEVELS;
  const int num_elems = get_homme_param_value<int>("nelemd");

  // Create the identifiers of input and output fields

  // Create layout of scalar/vector/tensor 2d/3d fields
  std::vector<FieldTag> tags_scalar_2d = {FieldTag::Element, FieldTag::GaussPoint, FieldTag::GaussPoint};
  std::vector<FieldTag> tags_vector_2d = {FieldTag::Element, FieldTag::Component, FieldTag::GaussPoint, FieldTag::GaussPoint};

  std::vector<FieldTag> tags_scalar_3d = {FieldTag::Element, FieldTag::GaussPoint, FieldTag::GaussPoint, FieldTag::VerticalLevel};
  std::vector<FieldTag> tags_vector_3d = {FieldTag::Element, FieldTag::Component, FieldTag::GaussPoint, FieldTag::GaussPoint, FieldTag::VerticalLevel};
  std::vector<FieldTag> tags_scalar_state_3d = {FieldTag::Element, FieldTag::TimeLevel, FieldTag::GaussPoint, FieldTag::GaussPoint, FieldTag::VerticalLevel};
  std::vector<FieldTag> tags_vector_state_3d = {FieldTag::Element, FieldTag::TimeLevel, FieldTag::Component, FieldTag::GaussPoint, FieldTag::GaussPoint, FieldTag::VerticalLevel};

  std::vector<int> dims_scalar_2d = {num_elems, NP, NP};
  std::vector<int> dims_scalar_3d_mid = {num_elems, NP, NP, NUM_PHYSICAL_LEV};
  std::vector<int> dims_scalar_3d_int = {num_elems, NP, NP, NUM_PHYSICAL_LEV+1};
  std::vector<int> dims_scalar_state_3d = {num_elems, NUM_TIME_LEVELS, NP, NP, NUM_PHYSICAL_LEV};
  std::vector<int> dims_vector_state_3d = {num_elems, NUM_TIME_LEVELS, 2, NP, NP, NUM_PHYSICAL_LEV};

  FieldLayout v_layout {tags_vector_state_3d,dims_vector_state_3d};
  FieldLayout t_dp_layout {tags_scalar_state_3d,dims_scalar_state_3d};

  m_computed_fields.emplace("horizontal velocity",v_layout,"Dynamics");
  m_computed_fields.emplace("temperature",t_dp_layout,"Dynamics");
  m_computed_fields.emplace("dp",t_dp_layout,"Dynamics");
}

void HommeDynamics::initialize (const std::shared_ptr<const GridsManager> /* grids_manager */)
{
  // We need to set the pointers in the C++ views to the ones contained in the scream
  // Fields *before* they ever get copied/filled. In particular, we need to make sure
  // that the Elements and Tracers structures contain scream Field's views before:
  //  - the Caar, Esf, Hvf, and Remap functors are created (cause they create copies)
  //  - any BoundaryExchange exchanging Elements/Tracers's views is created (cause they subview their views)

  auto& elements = Homme::Context::singleton().get_elements();
  auto& tracers  = Homme::Context::singleton().get_tracers();

  const int num_elems = elements.num_elems();
  constexpr int NUM_PHYSICAL_LEV  = Homme::NUM_PHYSICAL_LEV;
  constexpr int NUM_TIME_LEVELS   = Homme::NUM_TIME_LEVELS;
  constexpr int Q_NUM_TIME_LEVELS = Homme::Q_NUM_TIME_LEVELS;
  // constexpr int QSIZE_D           = Homme::QSIZE_D;
  using Scalar = Homme::Scalar;

  for (auto& it : m_dyn_fields_out) {
    const auto& name = it.first;
    const auto& f    = it.second;

    if (name=="horizontal velocity") {
      // Velocity
      auto& v = elements.m_v;
      auto v_in = f.get_reshaped_view<Real*[NUM_TIME_LEVELS][2][NP][NP][NUM_PHYSICAL_LEV]>();
      using v_type = std::remove_reference<decltype(v)>::type;
      v = v_type (reinterpret_cast<Scalar*>(v_in.data()),num_elems);
    } else if (name=="temperature") {
      // Temperature
      auto& t = elements.m_t;
      auto t_in = f.get_reshaped_view<Real*[NUM_TIME_LEVELS][NP][NP][NUM_PHYSICAL_LEV]>();
      using t_type = std::remove_reference<decltype(t)>::type;
      t = t_type(reinterpret_cast<Scalar*>(t_in.data()),num_elems);
    } else if (name=="dp") {
      // Levels thickness
      auto& dp = elements.m_dp3d;
      auto dp_in = f.template get_reshaped_view<Real*[NUM_TIME_LEVELS][NP][NP][NUM_PHYSICAL_LEV]>();
      using dp_type = std::remove_reference<decltype(dp)>::type;
      dp = dp_type(reinterpret_cast<Scalar*>(dp_in.data()),num_elems);
    } else if (name=="qdp") {
      // Tracers mass
      auto& qdp = tracers.qdp;
      auto qdp_in = f.template get_reshaped_view<Real*[Q_NUM_TIME_LEVELS][QSIZE_D][NP][NP][NUM_PHYSICAL_LEV]>();
      using qdp_type = std::remove_reference<decltype(qdp)>::type;
      qdp = qdp_type(reinterpret_cast<Scalar*>(qdp_in.data()),num_elems);
    } else {
      error::runtime_abort("Error! Unexpected field name. This is an internal error. Please, contact developers.\n");
    }
  }

  // Now that we set the correct pointers inside the Kokkos views, we can finish homme's initialization
  init_homme2_f90 ();
}

void HommeDynamics::register_fields (FieldRepository<Real, device_type>& field_repo) const
{
  for (const auto& fid : m_computed_fields) {
    field_repo.register_field(fid);
  }
}

void HommeDynamics::set_computed_field_impl (const Field<Real, device_type>& f)
{
  m_dyn_fields_out.emplace(f.get_header().get_identifier().name(),f);
}

void HommeDynamics::run (/* what inputs? */)
{
  run_homme_f90 ();
}

void HommeDynamics::finalize (/* what inputs? */)
{
  Homme::Context::singleton().finalize_singleton();
  Homme::MpiContext::singleton().finalize_singleton();
  finalize_homme_f90();
}

} // namespace scream
