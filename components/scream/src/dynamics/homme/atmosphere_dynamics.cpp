#include "atmosphere_dynamics.hpp"
#include "dynamics/homme/scream_homme_interface.hpp"
#include "dynamics/homme/hommexx_dimensions.hpp"
#include "dynamics/homme/homme_dynamics_helpers.hpp"

#include "ekat/ekat_assert.hpp"

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
#include "mpi/MpiContext.hpp"

namespace scream
{

HommeDynamics::HommeDynamics (const ekat::Comm& comm, const ekat::ParameterList& /* params */)
 : m_dynamics_comm (comm)
{
  init_homme1(comm);

  // Make Homme throw rather than abort. In Homme, abort causes finalization of Kokkos,
  // which is bad, since scream still has outstanding views.
  ::Homme::Session::m_throw_instead_of_abort = true;
}

void HommeDynamics::set_grids (const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;

  // The units of mixing ratio Q are technically non-dimensional.
  // Nevertheless, for output reasons, we like to see 'kg/kg'.
  auto Q = kg/kg;
  auto Qdp = Q * Pa;
  Q.set_string("kg/kg");
  Qdp.set_string("kg/kg Pa");

  constexpr int NGP  = HOMMEXX_NP;
  constexpr int QSZ  = HOMMEXX_QSIZE_D;
  constexpr int NVL  = HOMMEXX_NUM_PHYSICAL_LEV;
  constexpr int NTL  = HOMMEXX_NUM_TIME_LEVELS;
  constexpr int QNTL = HOMMEXX_Q_NUM_TIME_LEVELS;

  const auto dyn_grid_name = "SE Dynamics";
  const auto dyn_grid = grids_manager->get_grid(dyn_grid_name);
  const int ne = dyn_grid->get_num_local_dofs()/(NGP*NGP);

  // Sanity check for the grid. This should *always* pass, since Homme builds the grids
  EKAT_REQUIRE_MSG(get_homme_param_value<int>("nelemd")==ne,
                     "Error! The number of elements computed from the Dynamis grid num_dof()\n"
                     "       does not match the number of elements internal in Homme.\n");
  const int nmf = get_homme_param_value<int>("num momentum forcings");

  // Create the identifiers of input and output fields
  using namespace ShortFieldTagsNames;
  using Tags = std::vector<FieldTag>;
  Tags dyn_3d_scalar_tags        {EL,       GP,GP,VL};
  Tags dyn_3d_vector_tags        {EL,   CMP,GP,GP,VL};
  Tags dyn_3d_scalar_state_tags  {EL,TL,    GP,GP,VL};
  Tags dyn_3d_vector_state_tags  {EL,TL,CMP,GP,GP,VL};
  Tags dyn_3d_tracer_state_tags  {EL,TL,VAR,GP,GP,VL};

  // Inputs
  FieldLayout FQ_layout  { {EL,VAR,GP,GP,VL}, {ne,QSZ,NGP,NGP,NVL} };
  FieldLayout FM_layout  { {EL,CMP,GP,GP,VL}, {ne,nmf,NGP,NGP,NVL} };
  FieldLayout FT_layout  { {EL,    GP,GP,VL}, {ne,    NGP,NGP,NVL} };

  m_required_fields.emplace("FQ", FQ_layout, Q,          dyn_grid_name);
  m_required_fields.emplace("FM", FM_layout, m/pow(s,2), dyn_grid_name);
  m_required_fields.emplace("FT", FT_layout, K/s,        dyn_grid_name);

  // Outputs
  FieldLayout dyn_scalar_3d_mid_layout { dyn_3d_scalar_state_tags,  {ne, NTL,    NGP,NGP,NVL} };
  FieldLayout dyn_vector_3d_mid_layout { dyn_3d_vector_state_tags,  {ne, NTL,  2,NGP,NGP,NVL} };
  FieldLayout dyn_tracers_layout       { dyn_3d_tracer_state_tags,  {ne,QNTL,QSZ,NGP,NGP,NVL} };
  m_computed_fields.emplace("v",  dyn_vector_3d_mid_layout, m/s, dyn_grid_name);
  m_computed_fields.emplace("t",  dyn_scalar_3d_mid_layout, K,   dyn_grid_name);
  m_computed_fields.emplace("dp", dyn_scalar_3d_mid_layout, Pa,  dyn_grid_name);
  m_computed_fields.emplace("qdp",dyn_tracers_layout,       Qdp, dyn_grid_name);

  const int ftype = get_homme_param_value<int>("ftype");
  EKAT_REQUIRE_MSG(ftype==0 || ftype==2 || ftype==4,
                     "Error! The scream interface to homme *assumes* ftype to be 2 or 4.\n"
                     "       Found " + std::to_string(ftype) + " instead.\n");
}

void HommeDynamics::initialize (const util::TimeStamp& t0)
{
  m_current_ts = t0;

  // We need to set the pointers in the C++ views to the ones contained in the scream
  // Fields *before* they ever get copied/filled. In particular, we need to make sure
  // that the Elements and Tracers structures contain scream Field's views before:
  //  - the Caar, Esf, Hvf, and Remap functors are created (cause they create copies)
  //  - any BoundaryExchange exchanging Elements/Tracers's views is created (cause they subview their views)
  // If we don't do this, the Elements/Tracers structures *will* have views
  // containing pointers allocated by scream, but all the functors *will not*,
  // since they will contain a copy of the views *before* we had replaced the ptrs.

  auto& elements = Homme::Context::singleton().get_elements();
  auto& tracers  = Homme::Context::singleton().get_tracers();

  const int num_elems = elements.num_elems();
  using Scalar = Homme::Scalar;

  constexpr int NTL  = HOMMEXX_NUM_TIME_LEVELS;
  constexpr int NVL  = HOMMEXX_NUM_LEV;
  constexpr int QNTL = HOMMEXX_Q_NUM_TIME_LEVELS;

  // Computed fields
  for (auto& it : m_dyn_fields_out) {
    const auto& name = it.first;
    const auto& f    = it.second;

    if (name=="v") {
      // Velocity
      auto& v = elements.m_v;
      auto v_in = f.get_reshaped_view<Scalar*[NTL][2][NP][NP][NVL]>();
      using v_type = std::remove_reference<decltype(v)>::type;
      v = v_type (v_in.data(),num_elems);
    } else if (name=="t") {
      // Temperature
      auto& t = elements.m_t;
      auto t_in = f.get_reshaped_view<Scalar*[NTL][NP][NP][NVL]>();
      using t_type = std::remove_reference<decltype(t)>::type;
      t = t_type(t_in.data(),num_elems);
    } else if (name=="dp") {
      // Levels thickness
      auto& dp = elements.m_dp3d;
      auto dp_in = f.template get_reshaped_view<Scalar*[NTL][NP][NP][NVL]>();
      using dp_type = std::remove_reference<decltype(dp)>::type;
      dp = dp_type(dp_in.data(),num_elems);
    } else if (name=="qdp") {
      // Tracers mass
      auto& qdp = tracers.qdp;
      auto qdp_in = f.template get_reshaped_view<Scalar*[QNTL][QSIZE_D][NP][NP][NVL]>();
      using qdp_type = std::remove_reference<decltype(qdp)>::type;
      qdp = qdp_type(qdp_in.data(),num_elems);
    } else {
      ekat::error::runtime_abort("Error! Unexpected field name. This is an internal error. Please, contact developers.\n");
    }
  }

  // Required fields.
  // NOTE: Homme's Elements store all field as views to non-const data. Why? I'm not
  //       sure why we did that, but I assume it was to allow initialization.
  //       Therefore, we need to cast away the const from the scream input fields.
  // TODO: make Hommexx Elements structure store const views for stuff that is indeed (logically) const.
  for (auto& it : m_dyn_fields_in) {
    const auto& name = it.first;
    const auto& f    = it.second;

    if (name=="FQ") {
      // Tracers forcing
      auto& fq = tracers.fq;
      auto fq_in = f.template get_reshaped_view<const Scalar*[QSIZE_D][NP][NP][NVL]>();
      using fq_type = std::remove_reference<decltype(fq)>::type;
      auto non_const_ptr = const_cast<Scalar*>(fq_in.data());
      fq = fq_type(non_const_ptr,num_elems);
    } else if (name=="FM") {
      // Momemntum forcing
      auto& fm = elements.m_fm;
      // Use dynamic extent for second dimension, since preqx and theta have different extents
      auto fm_in = f.template get_reshaped_view<const Scalar**[NP][NP][NVL]>();
      using fm_type = std::remove_reference<decltype(fm)>::type;
      auto non_const_ptr = const_cast<Scalar*>(fm_in.data());
      fm = fm_type(non_const_ptr,num_elems);
    } else if (name=="FT") {
      // Temperature forcing
      auto& ft = elements.m_ft;
      auto ft_in = f.template get_reshaped_view<const Scalar*[NP][NP][NVL]>();
      using ft_type = std::remove_reference<decltype(ft)>::type;
      auto non_const_ptr = const_cast<Scalar*>(ft_in.data());
      ft = ft_type(non_const_ptr,num_elems);
    } else {
      ekat::error::runtime_abort("Error! Unexpected field name. This is an internal error. Please, contact developers.\n");
    }
  }

  // Now that we set the correct pointers inside the Kokkos views, we can finish homme's initialization
  init_homme2_f90 ();
}

void HommeDynamics::register_fields (FieldRepository<Real, device_type>& field_repo) const
{
  using Scalar = Homme::Scalar;
  for (const auto& fid : m_computed_fields) {
    field_repo.register_field<Scalar>(fid);
  }
  for (const auto& fid : m_required_fields) {
    field_repo.register_field<Scalar>(fid);
  }
}

void HommeDynamics::run (const Real dt)
{
  try {
    run_homme_f90 (dt);

    m_current_ts += dt;

    // Update all fields time stamp
    for (auto& it : m_dyn_fields_out) {
      it.second.get_header().get_tracking().update_time_stamp(m_current_ts);
    }
  } catch (std::exception& e) {
    ekat::error::runtime_abort(e.what());
  } catch (...) {
    ekat::error::runtime_abort("Something went wrong, but we don't know what.\n");
  }
}

void HommeDynamics::finalize (/* what inputs? */)
{
  Homme::Context::singleton().finalize_singleton();
  Homme::MpiContext::singleton().finalize_singleton();
  finalize_homme_f90();
}

void HommeDynamics::set_required_field_impl (const Field<const Real, device_type>& f) {
  // Store a copy of the field. We need this in order to do some tracking checks
  // at the beginning of the run call. Other than that, there would be really
  // no need to store a scream field here; we could simply set the view ptr
  // in the Homme's view, and be done with it.
  m_dyn_fields_in.emplace(f.get_header().get_identifier().name(),f);

  // Add myself as customer to the field
  this->add_me_as_customer(f);
}

void HommeDynamics::set_computed_field_impl (const Field<      Real, device_type>& f) {
  // Store a copy of the field. We need this in order to do some tracking updates
  // at the end of the run call. Other than that, there would be really
  // no need to store a scream field here; we could simply set the view ptr
  // in the Homme's view, and be done with it.
  const auto& name = f.get_header().get_identifier().name();
  m_dyn_fields_out.emplace(name,f);

  // Add myself as provider for the field
  this->add_me_as_provider(f);

  // Set extra data specifying whether this state corresponds to tracers,
  // as well as the current time level index (to be used, e.g., during
  // dyn-phys remap, or during I/O)
  int* tl_ptr;
  bool is_tracer;
  if (name=="qdp") {
    tl_ptr = &Homme::Context::singleton().get_time_level().np1_qdp;
    is_tracer = true;
  } else {
    tl_ptr = &Homme::Context::singleton().get_time_level().np1;
    is_tracer = false;
  }
  ekat::any tl, tracer;
  tl.reset<int*>(tl_ptr);
  tracer.reset<bool>(is_tracer);

  // Throw if this data was already set (who dared?)
  f.get_header_ptr()->set_extra_data("Current Time Level",tl,true);
  f.get_header_ptr()->set_extra_data("Is Tracer State",tracer,true);
}

} // namespace scream
