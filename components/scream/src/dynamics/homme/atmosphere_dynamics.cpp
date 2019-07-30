#include "atmosphere_dynamics.hpp"
#include "share/scream_assert.hpp"
#include "dynamics/homme/scream_homme_interface.hpp"
#include "dynamics/homme/hommexx_dimensions.hpp"
#include "Context.hpp"
#include "mpi/MpiContext.hpp"
#include "Types.hpp"
#include "Elements.hpp"
#include "Tracers.hpp"
#include "Dimensions.hpp"
#include "CaarFunctor.hpp"
#include "Hommexx_Session.hpp"
#include "EulerStepFunctor.hpp"
#include "HyperviscosityFunctor.hpp"

namespace scream
{

namespace util
{

// Specialize ScalarProperties on Homme's Scalar type
template<>
struct ScalarProperties<Homme::Scalar> {
  using scalar_type = Homme::Real;
  static constexpr bool is_pack = true;
};

} // namespace util

HommeDynamics::HommeDynamics (const Comm& comm,const ParameterList& /* params */)
 : m_dynamics_comm (comm)
{
  init_homme1(comm);

  // Make Homme throw rather than abort. In Homme, abort causes finalization of Kokkos,
  // which is bad, since scream still has outstanding views.
  ::Homme::Session::m_throw_instead_of_abort = true;
}

void HommeDynamics::set_grids (const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace units;

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

  auto grid = grids_manager->get_grid("Dynamics");
  const int num_dofs = grid->num_dofs();

  const int ne  = num_dofs/(NGP*NGP);

  scream_require_msg(get_homme_param_value<int>("nelemd")==ne,
                     "Error! The number of elements computed from the Dynamis grid num_dof()\n"
                     "       does not match the number of elements internal in Homme.\n");
  const int nmf = get_homme_param_value<int>("num momentum forcings");

  // Create the identifiers of input and output fields
  auto EL = FieldTag::Element;
  auto GP = FieldTag::GaussPoint;
  auto TL = FieldTag::TimeLevel;
  auto CM = FieldTag::Component;
  auto VL = FieldTag::VerticalLevel;
  auto VR = FieldTag::Variable;

  // Create layout of scalar/vector/tensor 2d/3d fields
  FieldLayout scalar2d_layout  { {EL,GP,GP}, {ne,NGP,NGP} };

  FieldLayout scalar_state_3d_mid_layout { {EL,TL,   GP,GP,VL} , {ne,NTL,  NGP,NGP,NVL}};
  FieldLayout vector_state_3d_mid_layout { {EL,TL,CM,GP,GP,VL} , {ne,NTL,2,NGP,NGP,NVL}};

  FieldLayout tracers_state_layout { {EL,TL,VR,GP,GP,VL}, {ne,QNTL,QSZ,NGP,NGP,NVL} };

  FieldLayout q_forcing_layout  { {EL,VR,GP,GP,VL}, {ne,QSZ,NGP,NGP,NVL} };
  FieldLayout m_forcing_layout  { {EL,VR,GP,GP,VL}, {ne,nmf,NGP,NGP,NVL} };
  FieldLayout t_forcing_layout  { {EL,   GP,GP,VL}, {ne,    NGP,NGP,NVL} };

  // Set requirements
  const int ftype = get_homme_param_value<int>("ftype");
  scream_require_msg(ftype==0 || ftype==2 || ftype==4,
                     "Error! The scream interface to homme *assumes* ftype to be 2 or 4.\n"
                     "       Found " + std::to_string(ftype) + " instead.\n");
  m_required_fields.emplace("phis", scalar2d_layout,  pow(m,2)/pow(s,2),"Dynamics");
  m_required_fields.emplace("FQ",   q_forcing_layout, Q,                "Dynamics");
  m_required_fields.emplace("FM",   m_forcing_layout, m/pow(s,2),       "Dynamics");
  m_required_fields.emplace("FT",   t_forcing_layout, K/s,              "Dynamics");

  // Set computed fields
  m_computed_fields.emplace("v",  vector_state_3d_mid_layout,m/s,"Dynamics");
  m_computed_fields.emplace("t",  scalar_state_3d_mid_layout,K,  "Dynamics");
  m_computed_fields.emplace("dp", scalar_state_3d_mid_layout,Pa, "Dynamics");
  m_computed_fields.emplace("qdp",tracers_state_layout,      Qdp,"Dynamics");
}

void HommeDynamics::initialize (const util::TimeStamp& t0)
{
  m_current_ts = t0;

  // We need to set the pointers in the C++ views to the ones contained in the scream
  // Fields *before* they ever get copied/filled. In particular, we need to make sure
  // that the Elements and Tracers structures contain scream Field's views before:
  //  - the Caar, Esf, Hvf, and Remap functors are created (cause they create copies)
  //  - any BoundaryExchange exchanging Elements/Tracers's views is created (cause they subview their views)

  auto& elements = Homme::Context::singleton().get_elements();
  auto& tracers  = Homme::Context::singleton().get_tracers();

  const int num_elems = elements.num_elems();
  using Scalar = Homme::Scalar;

  // Computed fields
  for (auto& it : m_dyn_fields_out) {
    const auto& name = it.first;
    const auto& f    = it.second;

    if (name=="v") {
      // Velocity
      auto& v = elements.m_v;
      auto v_in = f.get_reshaped_view<Scalar*[HOMMEXX_NUM_TIME_LEVELS][2][HOMMEXX_NP][HOMMEXX_NP][HOMMEXX_NUM_LEV]>();
      using v_type = std::remove_reference<decltype(v)>::type;
      v = v_type (v_in.data(),num_elems);
    } else if (name=="t") {
      // Temperature
      auto& t = elements.m_t;
      auto t_in = f.get_reshaped_view<Real*[HOMMEXX_NUM_TIME_LEVELS][HOMMEXX_NP][HOMMEXX_NP][HOMMEXX_NUM_PHYSICAL_LEV]>();
      using t_type = std::remove_reference<decltype(t)>::type;
      t = t_type(reinterpret_cast<Scalar*>(t_in.data()),num_elems);
    } else if (name=="dp") {
      // Levels thickness
      auto& dp = elements.m_dp3d;
      auto dp_in = f.template get_reshaped_view<Real*[HOMMEXX_NUM_TIME_LEVELS][HOMMEXX_NP][HOMMEXX_NP][HOMMEXX_NUM_PHYSICAL_LEV]>();
      using dp_type = std::remove_reference<decltype(dp)>::type;
      dp = dp_type(reinterpret_cast<Scalar*>(dp_in.data()),num_elems);
    } else if (name=="qdp") {
      // Tracers mass
      auto& qdp = tracers.qdp;
      auto qdp_in = f.template get_reshaped_view<Real*[HOMMEXX_NUM_TIME_LEVELS][QSIZE_D][HOMMEXX_NP][HOMMEXX_NP][HOMMEXX_NUM_PHYSICAL_LEV]>();
      using qdp_type = std::remove_reference<decltype(qdp)>::type;
      qdp = qdp_type(reinterpret_cast<Scalar*>(qdp_in.data()),num_elems);
    } else {
      error::runtime_abort("Error! Unexpected field name. This is an internal error. Please, contact developers.\n");
    }
  }

  // Required fields.
  // NOTE: Homme's Elements store all field as views to non-const data (due to initialization),
  //       so we need to cast away the const from the scream input fields.
  // TODO: make Hommexx Elements structure store const views for stuff that is indeed const.
  for (auto& it : m_dyn_fields_in) {
    const auto& name = it.first;
    const auto& f    = it.second;

    if (name=="phis") {
      // Surface geo-potential
      auto& phis = elements.m_phis;
      auto phis_in = f.template get_reshaped_view<const Real*[HOMMEXX_NP][HOMMEXX_NP]>();
      using phis_type = std::remove_reference<decltype(phis)>::type;
      auto non_const_ptr = const_cast<Real*>(phis_in.data());
      phis = phis_type(non_const_ptr,num_elems);
    } else if (name=="FQ") {
      // Tracers forcing
      auto& fq = tracers.fq;
      auto fq_in = f.template get_reshaped_view<const Real*[QSIZE_D][HOMMEXX_NP][HOMMEXX_NP][HOMMEXX_NUM_PHYSICAL_LEV]>();
      using fq_type = std::remove_reference<decltype(fq)>::type;
      auto non_const_ptr = reinterpret_cast<Scalar*>(const_cast<Real*>(fq_in.data()));
      fq = fq_type(non_const_ptr,num_elems);
    } else if (name=="FM") {
      // Momemntum forcing
      auto& fm = elements.m_fm;
      // Use dynamic extent for second dimension, since preqx and theta have different extents
      auto fm_in = f.template get_reshaped_view<const Real**[HOMMEXX_NP][HOMMEXX_NP][HOMMEXX_NUM_PHYSICAL_LEV]>();
      using fm_type = std::remove_reference<decltype(fm)>::type;
      auto non_const_ptr = reinterpret_cast<Scalar*>(const_cast<Real*>(fm_in.data()));
      fm = fm_type(non_const_ptr,num_elems);
    } else if (name=="FT") {
      // Temperature forcing
      auto& ft = elements.m_ft;
      auto ft_in = f.template get_reshaped_view<const Real*[HOMMEXX_NP][HOMMEXX_NP][HOMMEXX_NUM_PHYSICAL_LEV]>();
      using ft_type = std::remove_reference<decltype(ft)>::type;
      auto non_const_ptr = reinterpret_cast<Scalar*>(const_cast<Real*>(ft_in.data()));
      ft = ft_type(non_const_ptr,num_elems);
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
  for (const auto& fid : m_required_fields) {
    field_repo.register_field(fid);
  }
}

void HommeDynamics::run (const double dt)
{
  try {
    run_homme_f90 (dt);

    m_current_ts += dt;

    // Update all fields time stamp
    for (auto& it : m_dyn_fields_out) {
      it.second.get_header().get_tracking().update_time_stamp(m_current_ts);
    }
  } catch (std::exception& e) {
    error::runtime_abort(e.what());
  } catch (...) {
    error::runtime_abort("Something went wrong, but we don't know what.\n");
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
  f.get_header_ptr()->get_tracking().add_customer(weak_from_this());
}

void HommeDynamics::set_computed_field_impl (const Field<      Real, device_type>& f) {
  // Store a copy of the field. We need this in order to do some tracking updates
  // at the end of the run call. Other than that, there would be really
  // no need to store a scream field here; we could simply set the view ptr
  // in the Homme's view, and be done with it.
  m_dyn_fields_out.emplace(f.get_header().get_identifier().name(),f);

  // Add myself as provider for the field
  f.get_header_ptr()->get_tracking().add_provider(weak_from_this());
}

} // namespace scream
