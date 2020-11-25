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
#include "dynamics/homme/homme_inputs_initializer.hpp"

#include "ekat/ekat_assert.hpp"
#include "dynamics/homme/interface/scream_homme_interface.hpp"

namespace scream
{

HommeDynamics::HommeDynamics (const ekat::Comm& comm, const ekat::ParameterList& params)
 : m_params        (params)
 , m_dynamics_comm (comm)
{
  // Init homme context par struct
  if (!is_parallel_inited_f90()) {
    auto comm_f = MPI_Comm_c2f(comm.mpi_comm());
    init_parallel_f90(comm_f);
  }

  // Load homme params from namelist
  if (!is_params_inited_f90()) {
    auto nl_fname = m_params.get<std::string>("Homme namelist filename","namelist.nl");
    const char* nl_fname_c = nl_fname.c_str();
    init_params_f90(nl_fname_c);
  }

  // Init homme geometry structures
  if (!is_geometry_inited_f90()) {
    init_geometry_f90();
  }

  // Init prim structures
  if (!is_data_structures_inited_f90()) {
    prim_init_data_structures_f90 ();
  }
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

  const auto dyn_grid_name = "Dynamics";
  const auto dyn_grid = grids_manager->get_grid(dyn_grid_name);
  const int ne = dyn_grid->get_num_local_dofs()/(NGP*NGP);

  // Sanity check for the grid. This should *always* pass, since Homme builds the grids
  EKAT_REQUIRE_MSG(get_num_owned_elems_f90()==ne,
                     "Error! The number of elements computed from the Dynamis grid num_dof()\n"
                     "       does not match the number of elements internal in Homme.\n");
  const int nmf = get_homme_param<int>("num momentum forcings");

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
  FieldLayout dyn_scalar_3d_int_layout { dyn_3d_scalar_state_tags,  {ne, NTL,    NGP,NGP,NVL+1} };
  FieldLayout dyn_vector_3d_mid_layout { dyn_3d_vector_state_tags,  {ne, NTL,  2,NGP,NGP,NVL} };
  FieldLayout dyn_tracers_layout       { dyn_3d_tracer_state_tags,  {ne,QNTL,QSZ,NGP,NGP,NVL} };
  m_computed_fields.emplace("u",       dyn_vector_3d_mid_layout, m/s, dyn_grid_name);
  m_computed_fields.emplace("vtheta",  dyn_scalar_3d_mid_layout, K,   dyn_grid_name);
  m_computed_fields.emplace("phi",     dyn_scalar_3d_int_layout, Pa*pow(m,3)/kg, dyn_grid_name);
  m_computed_fields.emplace("w",       dyn_scalar_3d_int_layout, m/s, dyn_grid_name);
  m_computed_fields.emplace("dp",      dyn_scalar_3d_mid_layout, Pa,  dyn_grid_name);
  m_computed_fields.emplace("qdp",     dyn_tracers_layout,       Qdp, dyn_grid_name);

  const int ftype = get_homme_param<int>("ftype");
  EKAT_REQUIRE_MSG(ftype==0 || ftype==2 || ftype==4,
                     "Error! The scream interface to homme *assumes* ftype to be 2 or 4.\n"
                     "       Found " + std::to_string(ftype) + " instead.\n");

  auto standalone = m_params.get<bool>("Initialize Inputs", false);
  if (standalone) {
    m_initializer = create_field_initializer<HommeInputsInitializer>();
  }
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
  using Scalar = Homme::Scalar;

  constexpr int NTL  = HOMMEXX_NUM_TIME_LEVELS;
  constexpr int NVL  = HOMMEXX_NUM_LEV;
  constexpr int QNTL = HOMMEXX_Q_NUM_TIME_LEVELS;

  // Computed fields
  for (auto& it : m_dyn_fields_out) {
    const auto& name = it.first;
    const auto& f    = it.second;

    if (name=="u") {
      // Velocity
      auto& v = state.m_v;
      auto v_in = f.get_reshaped_view<Scalar*[NTL][2][NP][NP][NVL]>();
      using v_type = std::remove_reference<decltype(v)>::type;
      v = v_type (v_in.data(),num_elems);
    } else if (name=="vtheta") {
      // Virtual potential temperature
      auto& vtheta = state.m_vtheta_dp;
      auto vtheta_in = f.get_reshaped_view<Scalar*[NTL][NP][NP][NVL]>();
      using vtheta_type = std::remove_reference<decltype(vtheta)>::type;
      vtheta = vtheta_type(vtheta_in.data(),num_elems);
    } else if (name=="phi") {
      // Geopotential
      auto& phi = state.m_phinh_i;
      auto phi_in = f.get_reshaped_view<Scalar*[NTL][NP][NP][NVL]>();
      using phi_type = std::remove_reference<decltype(phi)>::type;
      phi = phi_type(phi_in.data(),num_elems);
    } else if (name=="w") {
      // Geopotential
      auto& w = state.m_w_i;
      auto w_in = f.get_reshaped_view<Scalar*[NTL][NP][NP][NVL]>();
      using w_type = std::remove_reference<decltype(w)>::type;
      w = w_type(w_in.data(),num_elems);
    } else if (name=="dp") {
      // Levels thickness
      auto& dp = state.m_dp3d;
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
  auto& forcing = Homme::Context::singleton().get<Homme::ElementsForcing>();
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

  // Now that we set the correct pointers inside the Kokkos views, we can finish homme's initialization

  auto standalone = m_params.get<bool>("Initialize Inputs", false);
  if (standalone) {
    for (auto& f : m_dyn_fields_in) {
      m_initializer->add_me_as_initializer(f.second);
      std::cout << "Added " << m_initializer->name() << " as initializer for " << f.second.get_header().get_identifier().name() << "\n";
    }
    for (auto f : m_dyn_fields_out) {
      m_initializer->add_me_as_initializer(f.second);
    }
    m_initializer->initialize_fields();
  } else {
    // Some I/O routine must have loaded initial states. Homme needs those values
    // to complete the model initialization. So, leverage Homme functions that
    // copy stuff from cxx views to f90 arrays, to fwd the data I/O loaded to
    // Homme's f90 model initialization routines
    // NOTE: this is not needed in the if branch above, since homme inits fields
    //       in fortran, so the initial states values are already available to model_init
    prim_copy_cxx_to_f90 ();
  }

  prim_init_model_f90 (standalone);

  ::Homme::Context::singleton().get<::Homme::SimulationParams>().print();
}

void HommeDynamics::register_fields (FieldRepository<Real>& field_repo) const
{
  using Scalar = Homme::Scalar;
  for (const auto& fid : m_computed_fields) {
    field_repo.register_field<Scalar>(fid);
  }
  for (const auto& fid : m_required_fields) {
    field_repo.register_field<Scalar>(fid);
  }
}

void HommeDynamics::run_impl (const Real dt)
{
  try {
    prim_run_f90 (dt);

    // Update all fields time stamp
    for (auto& it : m_dyn_fields_out) {
      it.second.get_header().get_tracking().update_time_stamp(this->timestamp());
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

  // Set extra data specifying whether this state corresponds to tracers
  // This is needed by the phys-dyn remapper to figure out which timeleve
  // to use in the homme's TimeLevel structure
  if (name=="qdp") {
    ekat::any tracer;
    tracer.reset<bool>(true);
    // Throw if this data was already set (who dared?)
    f.get_header_ptr()->set_extra_data("Is Tracer State",tracer,true);
  }
}

} // namespace scream
