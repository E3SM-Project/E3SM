#include "atmosphere_dynamics.hpp"
#include <memory>

// HOMMEXX Includes
#include "CaarFunctor.hpp"
#include "Context.hpp"
#include "Dimensions.hpp"
#include "Elements.hpp"
#include "EulerStepFunctor.hpp"
#include "ExecSpaceDefs.hpp"
#include "Hommexx_Session.hpp"
#include "HyperviscosityFunctor.hpp"
#include "utilities/MathUtils.hpp"
#include "utilities/SubviewUtils.hpp"
#include "TimeLevel.hpp"
#include "Tracers.hpp"
#include "Types.hpp"

// Scream includes
#include "dynamics/homme/homme_dimensions.hpp"
#include "dynamics/homme/homme_dynamics_helpers.hpp"

#include "ekat/ekat_assert.hpp"
#include "dynamics/homme/interface/scream_homme_interface.hpp"

namespace scream
{

HommeDynamics::HommeDynamics (const ekat::Comm& comm, const ekat::ParameterList& params)
 : m_params        (params)
 , m_dynamics_comm (comm)
{
  // Nothing to do here
}

void HommeDynamics::set_grids (const std::shared_ptr<const GridsManager> grids_manager)
{
  // Init prim structures
  // TODO: they should not be inited yet; should we error out if they are?
  //       I'm gonna say 'no', for now, cause it might be a pb with unit tests.
  if (!is_data_structures_inited_f90()) {
    prim_init_data_structures_f90 ();
  }

  // Note: time levels are just an expedient used by Homme to
  //  store temporaries in the RK timestepping schemes.
  //  It is best to have this extra array dimension (rather than,
  //  say, having NTL separate arrays) because of memory locality.
  //  At the end of Homme's timestep, only one of those slices
  //  will be meaningful. The phys-dyn remapper will use Homme's
  //  TimeLevel structure to know exactly where to copy data from/to
  //  during the remap.

  using namespace ShortFieldTagsNames;
  using namespace ekat::units;
  using FID = FieldIdentifier;

  constexpr int NGP  = HOMMEXX_NP;
  constexpr int NVL  = HOMMEXX_NUM_PHYSICAL_LEV;
  constexpr int NTL  = HOMMEXX_NUM_TIME_LEVELS;

  auto Q = kg/kg;
  Q.set_string("kg/kg");

  const auto dgn = "Dynamics";
  m_dyn_grid = grids_manager->get_grid(dgn);

  const int ne = m_dyn_grid->get_num_local_dofs()/(NGP*NGP);

  // Sanity check for the grid. This should *always* pass, since Homme builds the grids
  EKAT_REQUIRE_MSG(get_num_local_elems_f90()==ne,
      "Error! The number of elements computed from the Dynamics grid num_dof()\n"
      "       does not match the number of elements internal in Homme.\n");

  // Create the dyn grid identifiers of input and output fields
  FieldLayout scalar_state_3d_mid { {EL,TL,GP,GP,LEV},      {ne, NTL,    NGP,NGP,NVL} };
  FieldLayout scalar_state_3d_int { {EL,TL,GP,GP,ILEV},     {ne, NTL,    NGP,NGP,NVL+1} };
  FieldLayout vector_state_3d_mid { {EL,TL,CMP,GP,GP,LEV},  {ne, NTL,  2,NGP,NGP,NVL} };

  m_dyn_fids.emplace("v"        , FID("v",          vector_state_3d_mid, m/s,            dgn));
  m_dyn_fids.emplace("vtheta_dp", FID("vtheta_dp",  scalar_state_3d_mid, K,              dgn));
  m_dyn_fids.emplace("phi_i"    , FID("phi_i",      scalar_state_3d_int, Pa*pow(m,3)/kg, dgn));
  m_dyn_fids.emplace("w_i"      , FID("w_i",        scalar_state_3d_int, m/s,            dgn));
  m_dyn_fids.emplace("dp"       , FID("dp",         scalar_state_3d_mid, Pa,             dgn));

  // Create the PD remapper, and register fields
  m_p2d_remapper = grids_manager->create_remapper_from_ref_grid(m_dyn_grid);
  m_p2d_remapper->registration_begins();
  m_p2d_remapper->register_field_from_tgt(m_dyn_fids.at("v"));
  m_p2d_remapper->register_field_from_tgt(m_dyn_fids.at("vtheta_dp"));
  m_p2d_remapper->register_field_from_tgt(m_dyn_fids.at("phi_i"));
  m_p2d_remapper->register_field_from_tgt(m_dyn_fids.at("w_i"));
  m_p2d_remapper->register_field_from_tgt(m_dyn_fids.at("dp"));

  // Create the std::set of required/computed fids
  for (int i=0; i<m_p2d_remapper->get_num_registered_fields(); ++i) {
    const auto& ref_fid = m_p2d_remapper->get_src_field_id(i);
    m_required_fields.insert(ref_fid);
    m_computed_fields.insert(ref_fid);
  }

  // qv is needed to make sure Q is not empty (dyn needs qv to transform T<->Theta),
  // while ps is needed for initial conditions only.
  FieldLayout scalar_3d_mid { {EL,    GP,GP,LEV}, {ne,    NGP,NGP,NVL} };
  FieldIdentifier qv("qv", scalar_3d_mid, Q, dgn);
  m_required_fields.insert(m_p2d_remapper->create_src_fid(qv));

  const int ftype = get_homme_param<int>("ftype");
  EKAT_REQUIRE_MSG(ftype==0 || ftype==2 || ftype==4,
                     "Error! The scream interface to homme *assumes* ftype to be 2 or 4.\n"
                     "       Found " + std::to_string(ftype) + " instead.\n");

  // Input-output groups
  m_inout_groups_req.emplace("TRACERS",grids_manager->get_reference_grid()->name());
}

void HommeDynamics::
set_updated_group (const FieldGroup<Real>& group)
{
  const auto& name = group.m_info->m_group_name;
  EKAT_REQUIRE_MSG(name=="TRACERS",
    "Error! We were not expecting a field group called '" << name << "\n");

  EKAT_REQUIRE_MSG(not group.m_info->empty(),
    "Error! There should be at least one tracer (qv) in the 'TRACERS' group.\n");

  EKAT_REQUIRE_MSG(group.m_info->m_bundled,
      "Error! Homme expects a bundled field for tracers.\n");

  // Store the big Q array
  auto& Q_phys = m_ref_grid_fields["Q"] = *group.m_bundle;
  auto& Q_dyn = m_dyn_grid_fields["Q"];
  Q_dyn = Field<Real>(m_p2d_remapper->create_tgt_fid(group.m_bundle->get_header().get_identifier()));
  Q_dyn.get_header().get_alloc_properties().request_allocation<Homme::Scalar>();
  Q_dyn.allocate_view();

  // Set Q in the remapper
  m_p2d_remapper->register_field(Q_phys,Q_dyn);
  m_p2d_remapper->registration_ends();

  // Now that we have Q, we have the exact count for tracers,
  // and we can use that info to setup tracers stuff in Homme
  const int qsize = group.m_info->size();
  auto& params = Homme::Context::singleton().get<Homme::SimulationParams>();
  auto& tracers = Homme::Context::singleton().get<Homme::Tracers>();
  params.qsize = qsize;
  set_homme_param("qsize",qsize);
  tracers.init(tracers.num_elems(),qsize);

  // Set the dyn grid field's view in Hommexx data structures
  constexpr int NVL  = HOMMEXX_NUM_LEV;
  auto Q_view = Q_dyn.template get_reshaped_view<Homme::Scalar**[NP][NP][NVL]>();
  tracers.Q = decltype(tracers.Q)(Q_view.data(),tracers.num_elems(),qsize);
}

void HommeDynamics::initialize_impl (const util::TimeStamp& /* t0 */)
{
  // Print homme's parameters, so user can see whether something wasn't set right.
  ::Homme::Context::singleton().get<::Homme::SimulationParams>().print();

  // Now that all fields on the ref grid are set, create the fields on the

  // Now that we have all fields set in homme, let's remap the input fields,
  // so that Homme gets the correct Initial Conditions
  m_p2d_remapper->remap(true);

  // Finish homme initialization
  // Homme::initialize_dp3d_from_ps_c();

  prim_init_model_f90 ();
}

void HommeDynamics::register_fields (FieldRepository<Real>& field_repo) const
{
  // Inputs
  for (int i=0; i<m_p2d_remapper->get_num_registered_fields(); ++i) {
    const auto& src = m_p2d_remapper->get_src_field_id(i);
    field_repo.register_field(src);
  }

  // Make sure qv is registered (we need at least that one tracer, for T<->Theta conversions)
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  constexpr int NGP  = HOMMEXX_NP;
  constexpr int NVL  = HOMMEXX_NUM_PHYSICAL_LEV;

  auto Q = kg/kg;
  Q.set_string("kg/kg");

  const int ne = m_dyn_grid->get_num_local_dofs()/(NGP*NGP);
  FieldLayout qv_layout { {EL,    GP,GP,LEV}, {ne,    NGP,NGP,NVL} };
  FieldIdentifier qv("qv", qv_layout, Q,  m_dyn_grid->name());

  field_repo.register_field(m_p2d_remapper->create_src_fid(qv),"TRACERS");
}

void HommeDynamics::run_impl (const Real dt)
{
  try {
    // Remap inputs
    m_p2d_remapper->remap(true);
    // Run hommexx
    prim_run_f90 (dt);
    // Remap outputs
    m_p2d_remapper->remap(false);

    // Get a copy of the current timestamp (at the beginning of the step) and
    // advance it, updating the p3 fields.
    auto ts = timestamp();
    ts += dt;
    for (auto& it : m_ref_grid_fields) {
      it.second.get_header().get_tracking().update_time_stamp(ts);
    }
  } catch (std::exception& e) {
    ekat::error::runtime_abort(e.what());
  } catch (...) {
    ekat::error::runtime_abort("Something went wrong, but we don't know what.\n");
  }
}

void HommeDynamics::finalize_impl (/* what inputs? */)
{
  Homme::Context::singleton().finalize_singleton();
  prim_finalize_f90();
}

void HommeDynamics::set_required_field_impl (const Field<const Real>& f) {
  // Since all inputs are also outputs, we don't store a copy of f here.
  // Instead, we'll store a copy of it when it's given to us during
  // set_computed_field, since we'll have access to a non-const version
  // then, which we need for remapping purposes.

  // Add myself as customer to the field.
  this->add_me_as_customer(f);
}

void HommeDynamics::set_computed_field_impl (const Field<Real>& f) {
  auto& state = Homme::Context::singleton().get<Homme::ElementsState>();

  const int num_elems = state.num_elems();
  using Scalar = Homme::Scalar;

  constexpr int NTL  = HOMMEXX_NUM_TIME_LEVELS;
  constexpr int NVL  = HOMMEXX_NUM_LEV;
  constexpr int NVLI = HOMMEXX_NUM_LEV_P;

  const auto& name = f.get_header().get_identifier().name();

  // Create the dyn grid field
  auto& f_dyn = m_dyn_grid_fields[name];
  f_dyn = Field<Real>(m_dyn_fids.at(name));
  f_dyn.get_header().get_alloc_properties().request_allocation<Homme::Scalar>();
  f_dyn.allocate_view();

  if (name=="v") {
    // Velocity
    auto& v = state.m_v;
    auto v_in = f_dyn.get_reshaped_view<Scalar*[NTL][2][NP][NP][NVL]>();
    using v_type = std::remove_reference<decltype(v)>::type;
    v = v_type (v_in.data(),num_elems);
  } else if (name=="vtheta_dp") {
    // Virtual potential temperature
    auto& vtheta = state.m_vtheta_dp;
    auto vtheta_in = f_dyn.get_reshaped_view<Scalar*[NTL][NP][NP][NVL]>();
    using vtheta_type = std::remove_reference<decltype(vtheta)>::type;
    vtheta = vtheta_type(vtheta_in.data(),num_elems);
  } else if (name=="phi_i") {
    // Geopotential
    auto& phi = state.m_phinh_i;
    auto phi_in = f_dyn.get_reshaped_view<Scalar*[NTL][NP][NP][NVLI]>();
    using phi_type = std::remove_reference<decltype(phi)>::type;
    phi = phi_type(phi_in.data(),num_elems);
  } else if (name=="w_i") {
    // Geopotential
    auto& w = state.m_w_i;
    auto w_in = f_dyn.get_reshaped_view<Scalar*[NTL][NP][NP][NVLI]>();
    using w_type = std::remove_reference<decltype(w)>::type;
    w = w_type(w_in.data(),num_elems);
  } else if (name=="dp") {
    // Levels thickness
    auto& dp = state.m_dp3d;
    auto dp_in = f_dyn.template get_reshaped_view<Scalar*[NTL][NP][NP][NVL]>();
    using dp_type = std::remove_reference<decltype(dp)>::type;
    dp = dp_type(dp_in.data(),num_elems);
  } else {
    EKAT_ERROR_MSG("Error! Unexpected field name. This is an internal error. Please, contact developers.\n");
  }

  // Set the field in the remapper
  m_p2d_remapper->bind_field(f,f_dyn);

  // Store a copy of the field. We need this in order to do some tracking updates
  // at the end of the run call. Other than that, there would be really
  // no need to store a scream field here; we could simply set the view ptr
  // in the Homme's view, and be done with it.
  m_ref_grid_fields.emplace(name,f);

  // Add myself as provider for the field
  this->add_me_as_provider(f);
}

} // namespace scream
