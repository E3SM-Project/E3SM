#include "atmosphere_dynamics.hpp"

// HOMMEXX Includes
#include "CaarFunctor.hpp"
#include "Context.hpp"
#include "Dimensions.hpp"
#include "Elements.hpp"
#include "EulerStepFunctor.hpp"
#include "Hommexx_Session.hpp"
#include "HyperviscosityFunctor.hpp"
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

  constexpr int NGP  = HOMMEXX_NP;
  constexpr int NVL  = HOMMEXX_NUM_PHYSICAL_LEV;
  constexpr int NTL  = HOMMEXX_NUM_TIME_LEVELS;

  auto Q = kg/kg;
  Q.set_string("kg/kg");

  const auto dyn_grid_name = "Dynamics";
  m_dyn_grid = grids_manager->get_grid(dyn_grid_name);
  const int ne = m_dyn_grid->get_num_local_dofs()/(NGP*NGP);

  // Sanity check for the grid. This should *always* pass, since Homme builds the grids
  EKAT_REQUIRE_MSG(get_num_local_elems_f90()==ne,
      "Error! The number of elements computed from the Dynamics grid num_dof()\n"
      "       does not match the number of elements internal in Homme.\n");
  const int nmf = get_homme_param<int>("num momentum forcings");

  // Create the identifiers of input and output fields

  // Inputs
  FieldLayout vector_layout  { {EL,CMP,GP,GP,LEV}, {ne,nmf,NGP,NGP,NVL} };
  FieldLayout scalar_layout  { {EL,    GP,GP,LEV}, {ne,    NGP,NGP,NVL} };

  m_required_fields.emplace("FM", vector_layout, m/pow(s,2), dyn_grid_name);
  m_required_fields.emplace("FT", scalar_layout, K/s,        dyn_grid_name);
  m_required_fields.emplace("qv", scalar_layout, Q,          dyn_grid_name);

  // Inputs-Outputs
  FieldLayout dyn_scalar_3d_mid_layout { {EL,TL,GP,GP,LEV},      {ne, NTL,    NGP,NGP,NVL} };
  FieldLayout dyn_scalar_3d_int_layout { {EL,TL,GP,GP,ILEV},     {ne, NTL,    NGP,NGP,NVL+1} };
  FieldLayout dyn_vector_3d_mid_layout { {EL,TL,CMP,GP,GP,LEV},  {ne, NTL,  2,NGP,NGP,NVL} };

  m_required_fields.emplace("v",          dyn_vector_3d_mid_layout, m/s,            dyn_grid_name);
  m_required_fields.emplace("vtheta_dp",  dyn_scalar_3d_mid_layout, K,              dyn_grid_name);
  m_required_fields.emplace("phi_i",      dyn_scalar_3d_int_layout, Pa*pow(m,3)/kg, dyn_grid_name);
  m_required_fields.emplace("w_i",        dyn_scalar_3d_int_layout, m/s,            dyn_grid_name);
  m_required_fields.emplace("dp",         dyn_scalar_3d_mid_layout, Pa,             dyn_grid_name);
  m_computed_fields.emplace("v",          dyn_vector_3d_mid_layout, m/s,            dyn_grid_name);
  m_computed_fields.emplace("vtheta_dp",  dyn_scalar_3d_mid_layout, K,              dyn_grid_name);
  m_computed_fields.emplace("phi_i",      dyn_scalar_3d_int_layout, Pa*pow(m,3)/kg, dyn_grid_name);
  m_computed_fields.emplace("w_i",        dyn_scalar_3d_int_layout, m/s,            dyn_grid_name);
  m_computed_fields.emplace("dp",         dyn_scalar_3d_mid_layout, Pa,             dyn_grid_name);

  const int ftype = get_homme_param<int>("ftype");
  EKAT_REQUIRE_MSG(ftype==0 || ftype==2 || ftype==4,
                     "Error! The scream interface to homme *assumes* ftype to be 2 or 4.\n"
                     "       Found " + std::to_string(ftype) + " instead.\n");

  // Input-output groups
  m_inout_groups_req.emplace("TRACERS",m_dyn_grid->name());
  m_in_groups_req.emplace("TRACERS TENDENCY",m_dyn_grid->name());
}

// =========================================================================================
void HommeDynamics::
set_required_group (const FieldGroup<const Real>& group)
{
  if (not group.m_info->empty()) {
    const auto& name = group.m_info->m_group_name;

    EKAT_REQUIRE_MSG(name=="TRACERS TENDENCY",
      "Error! We were not expecting a field group called '" << name << "\n");

    EKAT_REQUIRE_MSG(group.m_info->m_bundled,
        "Error! Shoc expects bundled fields for tracers.\n");

    m_dyn_fields_in["FQ"] = *group.m_bundle;
  } else {
    // This must be a Homme-standalone run. Allocate a field manually
    // TODO: add check to make sure this is indeed a standalone run
    constexpr int qsize_d = HOMMEXX_QSIZE_D;
    const auto VAR = ShortFieldTagsNames::VAR;
    auto layout = m_dyn_grid->get_3d_vector_layout(true,VAR,qsize_d);
    auto nondim = ekat::units::Units::nondimensional();
    FieldIdentifier FQ_fid("FQ",layout,nondim,m_dyn_grid->name());
    Field<Real> FQ(FQ_fid);
    FQ.get_header().get_alloc_properties().request_allocation<Homme::Scalar>();
    FQ.allocate_view();
    m_dyn_fields_in.emplace("FQ",FQ);
  }
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
      "Error! Homme expects bundled fields for tracers.\n");

  // Store the big Q array
  m_dyn_fields_out["Q"] = *group.m_bundle;

  // Now that we have Q, we have the exact count for tracers,
  // and we can use that info to setup tracers stuff in Homme
  const int qsize = group.m_info->size();
  auto& params = Homme::Context::singleton().get<Homme::SimulationParams>();
  auto& tracers = Homme::Context::singleton().get<Homme::Tracers>();
  params.qsize = qsize;
  set_homme_param("qsize",qsize);
  tracers.init(tracers.num_elems(),qsize);
}

void HommeDynamics::initialize_impl (const util::TimeStamp& /* t0 */)
{
  // We need to set the pointers in the C++ views to the ones contained in the scream
  // Fields *before* they ever get copied/filled. In particular, we need to make sure
  // that the Elements and Tracers structures contain scream Field's views before:
  //  - the Caar, Esf, Hvf, and Remap functors are created (cause they create copies)
  //  - any BoundaryExchange exchanging Elements/Tracers's views is created (cause they subview their views)
  // If we don't do this, the Elements/Tracers structures *will* have views
  // containing pointers allocated by scream, but all the functors *will not*,
  // since they will contain a copy of the views *before* we had replaced the ptrs.

  auto& state = Homme::Context::singleton().get<Homme::ElementsState>();
  auto& tracers  = Homme::Context::singleton().get<Homme::Tracers>();

  const int num_elems = state.num_elems();
  const int num_tracers = tracers.num_tracers();
  using Scalar = Homme::Scalar;

  constexpr int NTL  = HOMMEXX_NUM_TIME_LEVELS;
  constexpr int NVL  = HOMMEXX_NUM_LEV;
  constexpr int NVLI = HOMMEXX_NUM_LEV_P;

  // Computed fields
  for (ci_string name : {"v", "vtheta_dp", "phi_i", "w_i", "dp", "Q"} ) {
    const auto& f = m_dyn_fields_out.at(name);

    if (name=="v") {
      // Velocity
      auto& v = state.m_v;
      auto v_in = f.get_reshaped_view<Scalar*[NTL][2][NP][NP][NVL]>();
      using v_type = std::remove_reference<decltype(v)>::type;
      v = v_type (v_in.data(),num_elems);
    } else if (name=="vtheta_dp") {
      // Virtual potential temperature
      auto& vtheta = state.m_vtheta_dp;
      auto vtheta_in = f.get_reshaped_view<Scalar*[NTL][NP][NP][NVL]>();
      using vtheta_type = std::remove_reference<decltype(vtheta)>::type;
      vtheta = vtheta_type(vtheta_in.data(),num_elems);
    } else if (name=="phi_i") {
      // Geopotential
      auto& phi = state.m_phinh_i;
      auto phi_in = f.get_reshaped_view<Scalar*[NTL][NP][NP][NVLI]>();
      using phi_type = std::remove_reference<decltype(phi)>::type;
      phi = phi_type(phi_in.data(),num_elems);
    } else if (name=="w_i") {
      // Geopotential
      auto& w = state.m_w_i;
      auto w_in = f.get_reshaped_view<Scalar*[NTL][NP][NP][NVLI]>();
      using w_type = std::remove_reference<decltype(w)>::type;
      w = w_type(w_in.data(),num_elems);
    } else if (name=="dp") {
      // Levels thickness
      auto& dp = state.m_dp3d;
      auto dp_in = f.template get_reshaped_view<Scalar*[NTL][NP][NP][NVL]>();
      using dp_type = std::remove_reference<decltype(dp)>::type;
      dp = dp_type(dp_in.data(),num_elems);
    } else if (name=="Q") {
      // Tracers mass
      auto& Q = tracers.Q;
      auto Q_in = f.template get_reshaped_view<Scalar**[NP][NP][NVL]>();
      using Q_type = std::remove_reference<decltype(Q)>::type;
      Q = Q_type(Q_in.data(),num_elems,num_tracers);
    } else {
      ekat::error::runtime_abort("Error! Unexpected field name. This is an internal error. Please, contact developers.\n");
    }
  }

  // Required fields.
  // NOTE: Homme's Elements store all field as views to non-const data. Why? I'm not
  //       sure why we did that, but I assume it was to allow initialization.
  //       Therefore, we need to cast away the const from the scream input fields.
  // TODO: make Hommexx Elements structure store const views for stuff that is indeed (logically) const.
  auto& forcing = Homme::Context::singleton().get<Homme::ElementsForcing>();
  for (ci_string name : {"FQ", "FM", "FT"}) {
    const auto& f = m_dyn_fields_in.at(name);

    if (name=="FQ") {
      // Tracers forcing
      auto& fq = tracers.fq;
      auto fq_in = f.template get_reshaped_view<const Scalar**[NP][NP][NVL]>();
      using fq_type = std::remove_reference<decltype(fq)>::type;
      auto non_const_ptr = const_cast<Scalar*>(fq_in.data());
      fq = fq_type(non_const_ptr,num_elems,num_tracers);
    } else if (name=="FM") {
      // Momemntum forcing
      auto& fm = forcing.m_fm;
      // Use dynamic extent for second dimension, since preqx and theta have different extents
      auto fm_in = f.template get_reshaped_view<const Scalar**[NP][NP][NVL]>();
      using fm_type = std::remove_reference<decltype(fm)>::type;
      auto non_const_ptr = const_cast<Scalar*>(fm_in.data());
      fm = fm_type(non_const_ptr,num_elems);
    } else if (name=="FT") {
      // Temperature forcing
      auto& ft = forcing.m_ft;
      auto ft_in = f.template get_reshaped_view<const Scalar*[NP][NP][NVL]>();
      using ft_type = std::remove_reference<decltype(ft)>::type;
      auto non_const_ptr = const_cast<Scalar*>(ft_in.data());
      ft = ft_type(non_const_ptr,num_elems);
    } else {
      ekat::error::runtime_abort("Error! Unexpected field name. This is an internal error. Please, contact developers.\n");
    }
  }

  prim_init_model_f90 ();

  ::Homme::Context::singleton().get<::Homme::SimulationParams>().print();
}

void HommeDynamics::register_fields (FieldRepository<Real>& field_repo) const
{
  using Scalar = Homme::Scalar;
  for (const auto& fid : m_computed_fields) {
    field_repo.register_field<Scalar>(fid);
  }
  for (const auto& fid : m_required_fields) {
    if (fid.name()=="qv") {
      field_repo.register_field<Scalar>(fid,"TRACERS");
    } else {
      field_repo.register_field<Scalar>(fid);
    }
  }
}

void HommeDynamics::run_impl (const Real dt)
{
  try {
    prim_run_f90 (dt);

    // Get a copy of the current timestamp (at the beginning of the step) and
    // advance it, updating the p3 fields.
    auto ts = timestamp();
    ts += dt;
    for (auto& it : m_dyn_fields_out) {
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
  // Store a copy of the field. We need this in order to do some tracking checks
  // at the beginning of the run call. Other than that, there would be really
  // no need to store a scream field here; we could simply set the view ptr
  // in the Homme's view, and be done with it.
  m_dyn_fields_in.emplace(f.get_header().get_identifier().name(),f);

  // Add myself as customer to the field
  this->add_me_as_customer(f);
}

void HommeDynamics::set_computed_field_impl (const Field<Real>& f) {
  // Store a copy of the field. We need this in order to do some tracking updates
  // at the end of the run call. Other than that, there would be really
  // no need to store a scream field here; we could simply set the view ptr
  // in the Homme's view, and be done with it.
  const auto& name = f.get_header().get_identifier().name();
  m_dyn_fields_out.emplace(name,f);

  // Add myself as provider for the field
  this->add_me_as_provider(f);
}

} // namespace scream
