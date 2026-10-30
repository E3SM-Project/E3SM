#include "eamxx_homme_process_interface.hpp"

// HOMMEXX Includes
#include "Context.hpp"
#include "Elements.hpp"
#include "ElementsState.hpp"
#include "HommexxEnums.hpp"
#include "SimulationParams.hpp"
#include "ElementsForcing.hpp"
#include "EulerStepFunctor.hpp"
#include "ComposeTransport.hpp"
#include "Diagnostics.hpp"
#include "DirkFunctor.hpp"
#include "ForcingFunctor.hpp"
#include "CaarFunctor.hpp"
#include "VerticalRemapManager.hpp"
#include "HyperviscosityFunctor.hpp"
#include "TimeLevel.hpp"
#include "Tracers.hpp"
#include "mpi/ConnectivityHelpers.hpp"
#include "utilities/MathUtils.hpp"
#include "utilities/SubviewUtils.hpp"
#include "ExecSpaceDefs.hpp"
#include "Types.hpp"

// Scream includes
#include "dynamics/homme/physics_dynamics_remapper.hpp"
#include "dynamics/homme/homme_dimensions.hpp"
#include "dynamics/homme/homme_dynamics_helpers.hpp"
#include "dynamics/homme/interface/eamxx_homme_interface.hpp"
#include "share/physics/physics_constants.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"
#include "share/util/eamxx_column_ops.hpp"
#include "share/property_checks/field_lower_bound_check.hpp"

// Ekat includes
#include <ekat_assert.hpp>
#include <ekat_team_policy_utils.hpp>
#include <ekat_subview_utils.hpp>
#include <ekat_pack.hpp>
#include <ekat_math_utils.hpp>
#include <ekat_units.hpp>
#include <ekat_pack_utils.hpp>
#include <ekat_workspace.hpp>

namespace scream
{

HommeDynamics::HommeDynamics (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params), m_phys_grid_pgN(-1)
{
  // This class needs Homme's context, so register as a user
  HommeContextUser::singleton().add_user();

  m_restart_extra_data["homme_nsteps"] = std::make_shared<std::any>(std::make_any<int>(-1));

  if (!is_parallel_inited_f90()) {
    // While we're here, we can init homme's parallel session
    auto fcomm = MPI_Comm_c2f(comm.mpi_comm());
    init_parallel_f90 (fcomm);
  }

  // Set the log filename in the F90 interface
  const char* logname = m_atm_logger->get_logfile_name().c_str();
  set_homme_log_file_name_f90 (&logname);

  m_bfb_hash_nstep = 0;
  if (params.isParameter("bfb_hash"))
    m_bfb_hash_nstep = std::max(0, params.get<int>("bfb_hash"));
}

HommeDynamics::~HommeDynamics ()
{
  // This class is done with Homme. Remove from its users list
  HommeContextUser::singleton().remove_user();
}

void HommeDynamics::create_requests ()
{
  // Grab dynamics, physics, and physicsGLL grids
  const auto dgn = "dynamics";
  m_dyn_grid = m_grids_manager->get_grid(dgn);
  m_phys_grid = m_grids_manager->get_grid("physics");
  m_cgll_grid = m_grids_manager->get_grid("physics_gll");

  m_instance = HommeInstance::get();

  fv_phys_set_grids();

  // Init prim structures
  // TODO: they should not be inited yet; should we error out if they are?
  //       I'm gonna say 'no', for now, cause it might be a pb with unit tests.
  if (!is_data_structures_inited_f90()) {
    // Set moisture in homme base on input file:
    const auto& moisture = m_params.get<std::string>("moisture");
    set_homme_param("moisture",ekat::upper_case(moisture)!="DRY");

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

  constexpr int NGP = HOMMEXX_NP;
  constexpr int NTL = HOMMEXX_NUM_TIME_LEVELS;
  constexpr int QTL = HOMMEXX_Q_NUM_TIME_LEVELS;
  constexpr int N   = HOMMEXX_PACK_SIZE;

  const int nelem = m_dyn_grid->get_num_local_dofs()/(NGP*NGP);
  const int nlev_mid = m_dyn_grid->get_num_vertical_levels();
  const int nlev_int = nlev_mid+1;

  // Sanity check for the grid. This should *always* pass, since Homme builds the grids
  EKAT_REQUIRE_MSG(get_num_local_elems_f90()==nelem,
      "Error! The number of elements computed from the dynamics grid num_dof()\n"
      "       does not match the number of elements internal in Homme.\n");

  const auto& params = Homme::Context::singleton().get<Homme::SimulationParams>();

  /*
     Explanation of what Homme needs, including what it needs for BFB restarts

     From the pt of view of the atmosphere, Homme has v,w,dp,T,phi,Q on physics grid
     as both inputs and outputs.
     Homme stores *prognostic states* on the dyn grid (v, w, dp, vtheta_dp, Qdp, phi),
     which have to be saved to file for BFB restart.
     Additionally, Homme computes forcing as:
       - FT (temperature forcing): FT_dyn=PD( (T_phys_new - T_phys_old)/dt )
       - FM (momentum forcing): FM_dyn=PD( (v_phys_new-v_phys_old)/dt) )
       - FQ (tracers forcing): FQ_dyn =  PD(Q_phys_new) (used as a hard adjustment or in increment calc)
     Here, PD(x) is quantity x remapped from physics to dynamics grid (viceversa for DP(x)).
     So for BFB restart, we need so store T_phys_old, [u,v,w]_phys_old.
     However, Homme computes Q=Qdp/dp at the end of the time step, so for BFB restarts we need to
     save Qdp only, and recompute Q_dyn_old=Qdp/dp.
   */
  // Notes:
  //  - physics will update T_mid, Q, horiz_winds, but not w_i, phi_i, and pseudo_density.
  //    Dyn will have to back out tendencies, as well as convert T->VTheta_dp.
  //    The simplest way is to store T, uv, and Q at the end of the previous dyn step.
  //    To make SCREAM more similar to what EAM does, we store Q_prev on dyn grid,
  //    while we keep uv_prev and T_prev on the phys grid (better for storage as well).
  //  - physics *does* update phi_i after each parametrization, but dyn discards it,
  //    in favor of reconstructing it from the FT, FQ, and the other thermodyn vars.

  // Note: the field T_mid will contain Temperature (at midpoints) at entry/exit.
  //       However, Homme uses VTheta*dp as state variable. Instead of keeping both
  //       T and VTheta*dp on dyn grid, we remap T_mid from phys directly
  //       into the n0 time-slice of Homme's vtheta_dp, and then do the conversion
  //       T_mid->VTheta_dp in place.

  const auto m2 = pow(m,2);
  const auto s2 = pow(s,2);

  // Note: qv is needed to transform T<->Theta

  const auto& pgn = m_phys_grid->name();
  auto pg_scalar2d     = m_phys_grid->get_2d_scalar_layout();
  auto pg_scalar3d_mid = m_phys_grid->get_3d_scalar_layout(LEV);
  auto pg_scalar3d_int = m_phys_grid->get_3d_scalar_layout(ILEV);
  auto pg_vector3d_mid = m_phys_grid->get_3d_vector_layout(LEV,2);
  add_field<Updated> ("horiz_winds",        pg_vector3d_mid, m/s,   pgn,N);
  add_field<Updated> ("T_mid",              pg_scalar3d_mid, K,     pgn,N);
  add_field<Computed>("pseudo_density",     pg_scalar3d_mid, Pa,    pgn,N);
  add_field<Computed>("pseudo_density_dry", pg_scalar3d_mid, Pa,    pgn,N);
  add_field<Updated> ("ps",                 pg_scalar2d    , Pa,    pgn);
  add_field<Updated >("phis",               pg_scalar2d    , m2/s2, pgn);
  add_field<Computed>("p_int",              pg_scalar3d_int, Pa,    pgn,N);
  add_field<Computed>("p_mid",              pg_scalar3d_mid, Pa,    pgn,N);
  add_field<Computed>("p_dry_int",          pg_scalar3d_int, Pa,    pgn,N);
  add_field<Computed>("p_dry_mid",          pg_scalar3d_mid, Pa,    pgn,N);
  add_field<Computed>("omega",              pg_scalar3d_mid, Pa/s,  pgn,N);
  if (params.do_3d_turbulence) {
    add_field<Required>("eddy_diff_heat_horiz", pg_scalar3d_mid, m2/s,  pgn,N);
    add_field<Required>("eddy_diff_mom_horiz",  pg_scalar3d_mid, m2/s,  pgn,N);
    auto pg_shear_components_mid = m_phys_grid->get_3d_vector_layout(LEV,6);
    add_field<Computed>("tke_shear_strain3d_components", pg_shear_components_mid, 1/s, pgn,N);
  }

  add_tracer<Updated >("qv", m_phys_grid, kg/kg, N);
  add_group<Updated>("tracers",pgn,N, MonolithicAlloc::Required);

  // Dynamics grid states
  create_helper_field("v_dyn",        {EL,TL,CMP,GP,GP,LEV}, {nelem,NTL,2,NP,NP,nlev_mid}, dgn);
  create_helper_field("vtheta_dp_dyn",{EL,TL,    GP,GP,LEV}, {nelem,NTL,  NP,NP,nlev_mid}, dgn);
  create_helper_field("dp3d_dyn",     {EL,TL,    GP,GP,LEV}, {nelem,NTL,  NP,NP,nlev_mid}, dgn);
  create_helper_field("w_int_dyn",    {EL,TL,    GP,GP,ILEV},{nelem,NTL,  NP,NP,nlev_int}, dgn);
  create_helper_field("phi_int_dyn",  {EL,TL,    GP,GP,ILEV},{nelem,NTL,  NP,NP,nlev_int}, dgn);
  create_helper_field("ps_dyn",       {EL,TL,    GP,GP},     {nelem,NTL,  NP,NP         }, dgn);
  create_helper_field("phis_dyn",     {EL,       GP,GP},     {nelem,      NP,NP         }, dgn);
  create_helper_field("omega_dyn",    {EL,       GP,GP,LEV}, {nelem,      NP,NP,nlev_mid}, dgn);
  create_helper_field("Qdp_dyn",      {EL,TL,CMP,GP,GP,LEV}, {nelem,QTL,HOMMEXX_QSIZE_D,NP,NP,nlev_mid},dgn);
  // Forcing (the tracers one, FQ_dyn, needs the number of tracers: see set_computed_group_impl)
  create_helper_field("FT_dyn",       {EL,       GP,GP,LEV}, {nelem,      NP,NP,nlev_mid}, dgn);
  create_helper_field("FM_dyn",       {EL,CMP,   GP,GP,LEV}, {nelem,    3,NP,NP,nlev_mid}, dgn);

  // The T and uv at the end of the last Homme step, on the physics grid. They are used to
  // compute the T and uv tendencies, so they need to be inited by HommeModelInit too.
  const int ncols_phys = m_phys_grid->get_num_local_dofs();
  create_helper_field("FT_phys",{COL,LEV},    {ncols_phys,  nlev_mid},pgn);
  create_helper_field("FM_phys",{COL,CMP,LEV},{ncols_phys,2,nlev_mid},pgn);
  m_instance->add_phys_field("FT_phys",m_helper_fields.at("FT_phys"));
  m_instance->add_phys_field("FM_phys",m_helper_fields.at("FM_phys"));
  if (params.do_3d_turbulence) {
    create_helper_field("Km_dyn",       {EL,       GP,GP,LEV}, {nelem,      NP,NP,nlev_mid}, dgn);
    create_helper_field("Kh_dyn",       {EL,       GP,GP,LEV}, {nelem,      NP,NP,nlev_mid}, dgn);
    create_helper_field("grad_Ux_dyn",  {EL,CMP,   GP,GP,LEV}, {nelem,2,    NP,NP,nlev_mid}, dgn);
    create_helper_field("grad_Uy_dyn",  {EL,CMP,   GP,GP,LEV}, {nelem,2,    NP,NP,nlev_mid}, dgn);
    create_helper_field("grad_Uz_dyn",  {EL,CMP,   GP,GP,LEV}, {nelem,2,    NP,NP,nlev_mid}, dgn);
    create_helper_field("shear_strain3d_components_dyn", {EL,CMP,GP,GP,LEV}, {nelem,6,NP,NP,nlev_mid}, dgn);
  }

  share_dyn_fields_with_instance({"v_dyn","vtheta_dp_dyn","dp3d_dyn","w_int_dyn","phi_int_dyn",
                                  "ps_dyn","phis_dyn","omega_dyn","Qdp_dyn","FT_dyn","FM_dyn"});
  if (params.do_3d_turbulence) {
    share_dyn_fields_with_instance({"Km_dyn","Kh_dyn","shear_strain3d_components_dyn"});
  }

  // For BFB restart, we need to read in the state on the dyn grid. The state above has NTL time slices,
  // but only one is really needed for restart. Therefore, we create "dynamic" subfields for
  // the state fields. This allows to save/read only the single slice needed
  // NOTE: the fcn init_time_level_c in Homme should really init also the qdp time levels,
  //       but it doesn't. So let's update them here. Notice that the qdp time levels update
  //       is "idempotent", meaning that it's safe to call it many times. In fact, all it
  //       does it to compute np1 and n0 from nsteps, so if nsteps is unchanged, calling
  //       the fcn twice won't "update" n0 and np1.
  auto& tl = Homme::Context::singleton().get<Homme::TimeLevel>();
  tl.update_tracers_levels(params.qsplit);

  // NOTE: 'true' is for 'dynamic' subfield; the idx of the subfield slice will move
  //       during execution (this class will take care of moving it, by calling
  //       reset_subview_idx on each field).
  add_internal_field (m_helper_fields.at("v_dyn").subfield(1,tl.n0,true),{"RESTART"});
  add_internal_field (m_helper_fields.at("vtheta_dp_dyn").subfield(1,tl.n0,true),{"RESTART"});
  add_internal_field (m_helper_fields.at("dp3d_dyn").subfield(1,tl.n0,true),{"RESTART"});
  add_internal_field (m_helper_fields.at("phi_int_dyn").subfield(1,tl.n0,true),{"RESTART"});
  add_internal_field (m_helper_fields.at("ps_dyn").subfield(1,tl.n0,true),{"RESTART"});
  add_internal_field (m_helper_fields.at("Qdp_dyn").subfield(1,tl.n0_qdp,true),{"RESTART"});
  if (not params.theta_hydrostatic_mode) {
    add_internal_field (m_helper_fields.at("w_int_dyn").subfield(1,tl.n0,true),{"RESTART"});
  }

  // The output manager pulls from the atm process fields. Add
  // helper fields for the case that a user request output.
  add_internal_field (m_helper_fields.at("omega_dyn"));
  add_internal_field (m_helper_fields.at("phis_dyn"),{"RESTART"});

  if (not fv_phys_active()) {
    // Dynamics backs out tendencies from the states, and passes those to Homme.
    // After Homme completes, we remap the updated state to the ref grid.  Thus,
    // is more convenient to use two different remappers: the pd remapper will
    // remap into Homme's forcing views, while the dp remapper will remap from
    // Homme's states.
    m_p2d_remapper = m_grids_manager->create_remapper(m_phys_grid,m_dyn_grid);
    m_d2p_remapper = m_grids_manager->create_remapper(m_dyn_grid,m_phys_grid);

    // Create separate remapper for restart
    m_ic_remapper = m_grids_manager->create_remapper(m_cgll_grid,m_dyn_grid);
  }

  // Layout for 2D (1d horiz X 1d vertical) variable
  const auto& grid_name = m_phys_grid->name();
  // Boundary flux fields for energy and mass conservation checks
  if (has_energy_fixer()) {
    add_field<Computed>("vapor_flux", pg_scalar2d, kg/(m2*s), grid_name);
    add_field<Computed>("water_flux", pg_scalar2d, m/s,     grid_name);
    add_field<Computed>("ice_flux",   pg_scalar2d, m/s,     grid_name);
    add_field<Computed>("heat_flux",  pg_scalar2d, W/m2,    grid_name);
  }

}//set_grids

size_t HommeDynamics::requested_buffer_size_in_bytes() const
{
  using namespace Homme;
  constexpr int num_turb3d_scratch_buffers = 8;
  constexpr int np2 = HOMMEXX_NP*HOMMEXX_NP;

  auto& c       = Context::singleton();
  auto& params  = c.get<SimulationParams>();

  const auto num_elems = c.get<Elements>().num_elems();
  const auto num_tracers = c.get<Tracers>().num_tracers();

  auto& caar = c.create_if_not_there<CaarFunctor>(num_elems,params);
  auto& hvf  = c.create_if_not_there<HyperviscosityFunctor>(num_elems, params);
  auto& ff   = c.create_if_not_there<ForcingFunctor>(num_elems, num_elems, params.qsize);
  auto& diag = c.create_if_not_there<Diagnostics> (num_elems, num_tracers, params.theta_hydrostatic_mode);
  auto& vrm  = c.create_if_not_there<VerticalRemapManager>(num_elems);

  const bool need_dirk = (params.time_step_type==TimeStepType::ttype7_imex ||
                          params.time_step_type==TimeStepType::ttype9_imex ||
                          params.time_step_type==TimeStepType::ttype10_imex  );

  // Request buffer sizes in FunctorsBuffersManager and then
  // return the total bytes using the calculated buffer size.
  auto& fbm  = c.create_if_not_there<FunctorsBuffersManager>();
  fbm.request_size(caar.requested_buffer_size());
  fbm.request_size(hvf.requested_buffer_size());
  fbm.request_size(diag.requested_buffer_size());
  fbm.request_size(ff.requested_buffer_size());
  fbm.request_size(vrm.requested_buffer_size());
  // Functors whose creation depends on the Homme namelist.
  if (params.transport_alg == 0) {
    auto& esf = c.create_if_not_there<EulerStepFunctor>(num_elems);
    fbm.request_size(esf.requested_buffer_size());
  } else {
#ifdef HOMME_ENABLE_COMPOSE
    auto& ct = c.create_if_not_there<ComposeTransport>(num_elems);
    fbm.request_size(ct.requested_buffer_size());
#endif
  }
  if (need_dirk) {
    // Create dirk functor only if needed
    auto& dirk = c.create_if_not_there<DirkFunctor>(num_elems);
    fbm.request_size(dirk.requested_buffer_size());
  }
  fv_phys_requested_buffer_size_in_bytes();

  size_t requested_bytes = fbm.allocated_size()*sizeof(Real);
  if (params.do_3d_turbulence) {
    const size_t ncols = num_elems*np2;
    requested_bytes += num_turb3d_scratch_buffers*sizeof(Real)*ncols*NUM_PHYSICAL_LEV;
  }

  return requested_bytes;
}

void HommeDynamics::init_buffers(const ATMBufferManager &buffer_manager)
{
  EKAT_REQUIRE_MSG(buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(),
                   "Error! Buffers size not sufficient.\n");

  using namespace Homme;
  auto& c = Context::singleton();
  auto& fbm  = c.get<FunctorsBuffersManager>();
  const auto& params = c.get<SimulationParams>();

  // Reset Homme buffer to use AD buffer memory.
  // Internally, homme will actually initialize its own buffers.
  EKAT_REQUIRE(buffer_manager.allocated_bytes()%sizeof(Real)==0); // Sanity check

  const int fbm_size = fbm.allocated_size();
  Real* mem = reinterpret_cast<Real*>(buffer_manager.get_memory());
  fbm.allocate(mem, fbm_size);
  mem += fbm_size;

  if (params.do_3d_turbulence) {
    constexpr int np2 = HOMMEXX_NP*HOMMEXX_NP;
    const int ncols = c.get<Elements>().num_elems()*np2;
    const int scratch_col_size = ncols*NUM_PHYSICAL_LEV;

    auto assign_scratch = [&](HommeDynamics::fixed_view_2d_phys& view) {
      view = HommeDynamics::fixed_view_2d_phys(mem, ncols);
      mem += scratch_col_size;
    };

    assign_scratch(m_w_mid_row_all);
    assign_scratch(m_w_mid_col_all);
    assign_scratch(m_dsdx_Ux_all);
    assign_scratch(m_dsdy_Ux_all);
    assign_scratch(m_dsdx_Uy_all);
    assign_scratch(m_dsdy_Uy_all);
    assign_scratch(m_dsdx_Uz_all);
    assign_scratch(m_dsdy_Uz_all);
  }

  size_t used_mem = (mem - buffer_manager.get_memory())*sizeof(Real);
  EKAT_REQUIRE_MSG(used_mem==requested_buffer_size_in_bytes(),
                   "Error! Used memory != requested memory for HommeDynamics.");
}

void HommeDynamics::initialize_impl (const RunType run_type)
{
  const auto& pgn = m_phys_grid->name();

  // Use common/shorter names for tracers.
  alias_group_in  ("tracers",pgn,"Q");
  alias_group_out ("tracers",pgn,"Q");

  // Grab handles of some Homme data structure
  const auto& c       = Homme::Context::singleton();
  const auto& params  = c.get<Homme::SimulationParams>();

  // Homme overwrites these diagnostic components after each dynamics step when 3D turbulence
  // is enabled, so start from a benign value (to avoid NaNs in an initial output).
  // NOTE: with the pg grid, HommeModelInit has already done this, before the first D->P remap.
  if (not fv_phys_active() and params.do_3d_turbulence) {
    m_helper_fields.at("shear_strain3d_components_dyn").deep_copy(0);
    get_field_out("tke_shear_strain3d_components").deep_copy(0);
  }

  // ------ Sanity checks ------- //

  // Nobody should claim to be a provider for dp.
  // WARNING! If the assumption on 'pseudo_density' ceases to be true, you have to revisit
  //          how you restart homme. In particular, p_mid is restarted from pseudo_density,
  //          as it is read from restart file. If other procs update it, the restarted value
  //          might no longer match the end-of-homme-step value, which is what you need
  //          to compute p_mid. Hence, if this assumption goes away, you need to restart
  //          p_mid by first remapping the restarted dp3d_dyn back to ref grid, and using
  //          that value to compute p_mid. Or, perhaps easier, write p_mid to restart file.
  EKAT_REQUIRE_MSG (
      get_field_out("pseudo_density",pgn).get_header().get_tracking().get_providers().size()==1,
      "Error! Someone other than dynamics is trying to update the pseudo_density.\n");

  if (not fv_phys_active()) {
    // Setup the p2d and d2p remappers

    // ftype==FORCING_0:
    //  1) remap Q_pgn->FQ_dyn
    //  2) compute FQ_dyn=dp3d*(FQ_dyn-Q_dyn_old)/dt
    // ftype!=FORCING_0:
    //  1) remap Q_pgn->FQ_dyn
    // Remap Q directly into FQ, tendency computed in pre_process step
    m_p2d_remapper->register_field(get_group_out("Q",pgn).monolithic_field(),m_helper_fields.at("FQ_dyn"));
    m_p2d_remapper->register_field(m_helper_fields.at("FT_phys"),m_helper_fields.at("FT_dyn"));

    // FM has 3 components on dyn grid, but only 2 on phys grid
    auto& FM_phys = m_helper_fields.at("FM_phys");
    auto& FM_dyn  = m_helper_fields.at("FM_dyn");
    m_p2d_remapper->register_field(FM_phys.get_component(0),FM_dyn.get_component(0));
    m_p2d_remapper->register_field(FM_phys.get_component(1),FM_dyn.get_component(1));

    // NOTE: for states, if/when we can remap subfields, we can remap the corresponding internal fields,
    //       which are subviews of the corresponding helper field at time slice np1
    m_d2p_remapper->register_field(get_internal_field("vtheta_dp_dyn"),get_field_out("T_mid"));
    m_d2p_remapper->register_field(get_internal_field("v_dyn"),get_field_out("horiz_winds"));
    m_d2p_remapper->register_field(get_internal_field("dp3d_dyn"), get_field_out("pseudo_density"));
    m_d2p_remapper->register_field(get_internal_field("ps_dyn"), get_field_out("ps"));
    m_d2p_remapper->register_field(m_helper_fields.at("Q_dyn"),get_group_out("Q",pgn).monolithic_field());
    m_d2p_remapper->register_field(m_helper_fields.at("omega_dyn"), get_field_out("omega"));

    if (params.do_3d_turbulence) {
      // Remap SHOC eddy diffusivities from physics grid to dynamics grid.
      m_p2d_remapper->register_field(get_field_in("eddy_diff_mom_horiz",pgn),m_helper_fields.at("Km_dyn"));
      m_p2d_remapper->register_field(get_field_in("eddy_diff_heat_horiz",pgn),m_helper_fields.at("Kh_dyn"));

      // Remap horizontal/local strain tensor components from dynamics to physics grid.
      m_d2p_remapper->register_field(m_helper_fields.at("shear_strain3d_components_dyn"), get_field_out("tke_shear_strain3d_components"));
    }

    m_p2d_remapper->registration_ends();
    m_d2p_remapper->registration_ends();
  }

  // Print some info about Homme setup
  print_homme_info ();

  // Finalize what HommeModelInit could not do (it has no access to the atm process interface)
  if (run_type==RunType::Initial) {
    initialize_homme_state ();
  } else {
    restart_homme_state ();
  }

  // With the pg grid, HommeModelInit remapped the dyn state to the physics grid, so we can now
  // compute the pressure. For the GLL grid, this is done in initialize/restart_homme_state.
  if (fv_phys_active()) {
    update_pressure (m_phys_grid);
  }

  // Since we just inited them, ensure p_mid/p_int (dry and wet) timestamps are valid
  get_field_out("p_int")    .get_header().get_tracking().update_time_stamp(start_of_step_ts());
  get_field_out("p_mid")    .get_header().get_tracking().update_time_stamp(start_of_step_ts());
  get_field_out("p_dry_int").get_header().get_tracking().update_time_stamp(start_of_step_ts());
  get_field_out("p_dry_mid").get_header().get_tracking().update_time_stamp(start_of_step_ts());

  // Set up field property checks
  // Note: We are seeing near epsilon negative values in a handful of places.
  //       We figured out the source of those values in Homme (see #1842),
  //       but they should be "benign" neg values (due to roundoffs), which we
  //       can safely clip to 0..
  using Interval = FieldWithinIntervalCheck;
  using LowerBound = FieldLowerBoundCheck;

  add_postcondition_check<LowerBound>(get_group_out("Q",pgn).monolithic_field(),m_phys_grid,0,true);
  add_postcondition_check<Interval>(get_field_out("T_mid",pgn),m_phys_grid,100.0, 500.0,false);
  add_postcondition_check<Interval>(get_field_out("horiz_winds",pgn),m_phys_grid,-400.0, 400.0,false);
  add_postcondition_check<Interval>(get_field_out("ps"),m_phys_grid,30000.0, 120000.0,false);

  // Initialize Rayleigh friction variables
  rayleigh_friction_init();

}//initialize_impl

void HommeDynamics::run_impl (const double dt)
{
  try {

    // Note: Homme's step lasts homme_dt*max(dt_remap_factor,dt_tracers_factor), and it must divide dt.
    // We neeed to compute dt/homme_dt, and subcycle homme that many times

    // NOTE: we did not have atm_dt when we inited homme, so we set nsplit=1.
    //       Now we can compute the actual nsplit, and need to update its value
    //       in Hommexx's data structures.
    //       Also, nsplit calculation requires an integer atm timestep, so we
    //       check to ensure that's the case
    EKAT_REQUIRE_MSG (std::abs(dt-std::round(dt))<std::numeric_limits<double>::epsilon()*10,
        "[HommeDynamics] Error! Input timestep departure from integer above tolerance.\n"
        "  - input dt : " << dt << "\n"
        "  - tolerance: " << std::numeric_limits<double>::epsilon()*10 << "\n");

    if (m_bfb_hash_nstep > 0 && start_of_step_ts().get_num_steps() % m_bfb_hash_nstep == 0)
      print_fast_global_state_hash("Hommexx",start_of_step_ts());

    const int dt_int = static_cast<int>(std::round(dt));

    const int nsplit = get_homme_nsplit_f90(dt_int);
    const auto& c = Homme::Context::singleton();
    auto& params = c.get<Homme::SimulationParams>();
    params.nsplit = nsplit;

    Kokkos::fence();
    homme_pre_process (dt);

    for (int subiter=0; subiter<nsplit; ++subiter) {
      Kokkos::fence();
      prim_run_f90(/* nsplit_iteration = */ subiter+1);
    }

    // This is where we will compute the strain term needed for Shear Production of TKE
    if (params.do_3d_turbulence){
      compute_horizontal_derivs_of_car_velocity();
      compute_local_strain_components3d();
    }

    // Update nstep in the restart extra data, so it can be written to restart if needed.
    const auto& tl = c.get<Homme::TimeLevel>();
    std::any_cast<int&>(*m_restart_extra_data["homme_nsteps"]) = tl.nstep;

    // Post process Homme's output, to produce what the rest of Atm expects
    Kokkos::fence();
    homme_post_process (dt);
  } catch (std::exception& e) {
    EKAT_ERROR_MSG(e.what());
  } catch (...) {
    EKAT_ERROR_MSG("Something went wrong, but we don't know what.\n");
  }
}

void HommeDynamics::finalize_impl (/* what inputs? */)
{
  prim_finalize_f90();

  // This class is done needing Homme's context, so remove myself as customer
  Homme::Context::singleton().finalize_singleton();
}


void HommeDynamics::set_computed_group_impl (const FieldGroup& group)
{
  const auto& c = Homme::Context::singleton();
        auto& tracers = c.get<Homme::Tracers>();

  if (group.name()=="tracers") {
    // The groups 'tracers' and 'tracers_mass_dyn' should contain the same fields
    const int qsize = group.size();
    EKAT_REQUIRE_MSG(qsize,
      "Error! There should be at least one tracer (qv) in the tracers group.\n");

    // Set runtime number of tracers in Homme
    auto& params = c.get<Homme::SimulationParams>();
    params.qsize = qsize;           // Set in the CXX data structure
    set_homme_param("qsize",qsize); // Set in the F90 module
    tracers.init(tracers.num_elems(),qsize);

    // Now that we know the number of tracers, create the dyn fields that depend on it
    // NOTE: unlike Qdp_dyn, these are padded neither to HOMMEXX_QSIZE_D, since they are
    //       remapped with fields that have exactly qsize components
    constexpr int NGP = HOMMEXX_NP;
    const int nelem = m_dyn_grid->get_num_local_dofs()/(NGP*NGP);
    const int nlevs = m_dyn_grid->get_num_vertical_levels();
    using namespace ShortFieldTagsNames;
    create_helper_field("FQ_dyn",{EL,CMP,GP,GP,LEV},{nelem,qsize,NGP,NGP,nlevs},m_dyn_grid->name());
    create_helper_field("Q_dyn" ,{EL,CMP,GP,GP,LEV},{nelem,qsize,NGP,NGP,nlevs},m_dyn_grid->name());
    share_dyn_fields_with_instance({"FQ_dyn","Q_dyn"});
  }
}

void HommeDynamics::homme_pre_process (const double dt) {
  // T and uv tendencies are backed out on the ref grid.
  // Homme takes care of turning the FT tendency into a tendency for VTheta_dp.

  using namespace Homme;
  const auto& c = Context::singleton();
  const auto& params = c.get<SimulationParams>();

  const int ncols = m_phys_grid->get_num_local_dofs();
  const int nlevs = m_phys_grid->get_num_vertical_levels();
  const int npacks = ekat::npack<Pack>(nlevs);

  const auto& pgn = m_phys_grid->name();

  // At the beginning of the step, FT and FM store T_prev and V_prev,
  // the temperature and (3d) velocity at the end of the previous
  // Homme step
  auto T  = get_field_in("T_mid",pgn).get_view<const Pack**>();
  auto v  = get_field_in("horiz_winds",pgn).get_view<const Pack***>();
  auto FT = m_helper_fields.at("FT_phys").get_view<Pack**>();
  auto FM = m_helper_fields.at("FM_phys").get_view<Pack***>();

  // If there are other atm procs updating the vertical velocity,
  // then we need to compute forcing for w as well
  Kokkos::parallel_for(KT::RangePolicy(0,ncols*npacks),
                       KOKKOS_LAMBDA(const int& idx) {
    const int icol = idx / npacks;
    const int ilev = idx % npacks;

    // Temperature forcing
    // Note: Homme takes care of converting ft into a forcing for vtheta
    const auto& t_new =  T(icol,ilev);
    const auto& t_old = FT(icol,ilev);
          auto& ft    = FT(icol,ilev);
    ft = (t_new - t_old) / dt;

    // Horizontal velocity forcing
    const auto& u_new =  v(icol,0,ilev);
    const auto& u_old = FM(icol,0,ilev);
          auto& fu    = FM(icol,0,ilev);
    fu = (u_new - u_old) / dt;

    const auto& v_new =  v(icol,1,ilev);
    const auto& v_old = FM(icol,1,ilev);
          auto& fv    = FM(icol,1,ilev);
    fv = (v_new - v_old) / dt;
  });

  const auto ftype = params.ftype;

  if (fv_phys_active()) {
    fv_phys_pre_process();
  } else {
    // Remap FT, FM, and Q->FQ
    m_p2d_remapper->remap_fwd();
  }

  auto& tl = c.get<TimeLevel>();

  // Note: np1_qdp and n0_qdp are 'deduced' from tl.nstep, so the
  //       following call may not even change them (i.e., they are
  //       not updated regardless). So if they were already up-to-date,
  //       the following call will do nothing.
  tl.update_tracers_levels(params.qsplit);

  // At this point, FQ contains Qnew (coming from physics).
  // Depending on ftype, we are going to do different things:
  //  ftype=0: FQ = dp*(Qnew-Qold) / dt
  //  ftype=2: nothing
  if (ftype == ForcingAlg::FORCING_0) {
    // Back out tracers tendency for Qdp
    const auto& tracers = c.get<Tracers>();
    const auto& state = c.get<ElementsState>();
    auto Q    = tracers.Q;
    auto FQ   = tracers.fq;
    auto dp3d = state.m_dp3d;
    const auto n0 = tl.n0;  // The time level where pd coupling remapped into
    constexpr int NVL = HOMMEXX_NUM_LEV;
    const int qsize = params.qsize;
    const auto Q_dyn = m_helper_fields.at("Q_dyn").get_view<Homme::Scalar*****>();
    Kokkos::parallel_for(Kokkos::RangePolicy<>(0,Q.size()),KOKKOS_LAMBDA(const int idx) {
      const int ie = idx / (qsize*NP*NP*NVL);
      const int iq = (idx / (NP*NP*NVL)) % qsize;
      const int ip = (idx / (NP*NVL)) % NP;
      const int jp = (idx / NVL) % NP;
      const int k  =  idx % NVL;

      // fq is currently storing q_new
      const auto& q_prev = Q(ie,iq,ip,jp,k);
      const auto& dp     = dp3d(ie,n0,ip,jp,k);
      auto& fq     = FQ(ie,iq,ip,jp,k);

      fq -= q_prev;
      fq /= dt;
      fq *= dp;
    });
    Kokkos::fence();
  }
}

void HommeDynamics::homme_post_process (const double dt) {
  const auto& pgn = m_phys_grid->name();
  const auto& c = Homme::Context::singleton();
  const auto& params = c.get<Homme::SimulationParams>();
  const auto& tl = c.get<Homme::TimeLevel>();

  // The internal fields are dynamic subfields of the homme states.
  // We need to update the slice they are subviewing
  get_internal_field("v_dyn").get_header().get_alloc_properties().reset_subview_idx(tl.n0);
  get_internal_field("vtheta_dp_dyn").get_header().get_alloc_properties().reset_subview_idx(tl.n0);
  get_internal_field("dp3d_dyn").get_header().get_alloc_properties().reset_subview_idx(tl.n0);
  get_internal_field("phi_int_dyn").get_header().get_alloc_properties().reset_subview_idx(tl.n0);
  get_internal_field("ps_dyn").get_header().get_alloc_properties().reset_subview_idx(tl.n0);
  get_internal_field("Qdp_dyn").get_header().get_alloc_properties().reset_subview_idx(tl.np1_qdp);
  if (not params.theta_hydrostatic_mode) {
    get_internal_field("w_int_dyn").get_header().get_alloc_properties().reset_subview_idx(tl.n0);
  }

  if (fv_phys_active()) {
    fv_phys_post_process();
    // Apply Rayleigh friction to update temperature and horiz_winds
    rayleigh_friction_apply(dt);

    return;
  }

  // Remap outputs to ref grid
  m_d2p_remapper->remap_fwd();

  using ColOps = ColumnOps<DefaultDevice,Real>;
  using PF = PhysicsFunctions<DefaultDevice>;
  using TPF = ekat::TeamPolicyFactory<KT::ExeSpace>;

  // Convert VTheta_dp->T, store T,uv, and possibly w in FT, FM,
  // compute p_int on ref grid.
  const auto dp_view = get_field_out("pseudo_density").get_view<Pack**>();
  const auto p_mid_view = get_field_out("p_mid").get_view<Pack**>();
  const auto p_int_view = get_field_out("p_int").get_view<Pack**>();
  const auto dp_dry_view    = get_field_out("pseudo_density_dry").get_view<Pack**>();
  const auto p_dry_int_view = get_field_out("p_dry_int").get_view<Pack**>();
  const auto p_dry_mid_view = get_field_out("p_dry_mid").get_view<Pack**>();
  const auto Q_view   = get_group_out("Q",pgn).monolithic_field().get_view<Pack***>();

  const auto T_view  = get_field_out("T_mid").get_view<Pack**>();
  const auto T_prev_view = m_helper_fields.at("FT_phys").get_view<Pack**>();

  m_helper_fields.at("FM_phys").deep_copy(get_field_out("horiz_winds"));

  const auto ncols = m_phys_grid->get_num_local_dofs();
  const auto nlevs = m_phys_grid->get_num_vertical_levels();
  const auto npacks= ekat::npack<Pack>(nlevs);

  const auto policy = TPF::get_thread_range_parallel_scan_team_policy(ncols,npacks);

  // Establish the boundary condition for the TOA
  const auto& hvcoord = c.get<Homme::HybridVCoord>();
  const auto ps0 = hvcoord.ps0 * hvcoord.hybrid_ai0;

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const KT::MemberType& team) {
    const int& icol = team.league_rank();

    auto qv  = ekat::subview(Q_view,icol,0);
    // TODO: Here we update the wet and dry pressure coordinates which is the same set of code used
    //       in the update_pressure() subroutine.  A low-priority todo item would clean-up the
    //       interface to call update_pressure here or change that routine so that it can be
    //       called within the top-level Kokkos parallel_for loop.
    // Compute p_int and p_mid
    auto dp = ekat::subview(dp_view,icol);
    auto p_mid = ekat::subview(p_mid_view,icol);
    auto p_int = ekat::subview(p_int_view,icol);

    ColOps::column_scan<true>(team,nlevs,dp,p_int,ps0);
    team.team_barrier();
    ColOps::compute_midpoint_values(team,nlevs,p_int,p_mid);
    team.team_barrier();
    // Compute p_dry_int and p_dry_mid
    auto dp_dry = ekat::subview(dp_dry_view,icol);
    auto p_dry_mid = ekat::subview(p_dry_mid_view,icol);
    auto p_dry_int = ekat::subview(p_dry_int_view,icol);

    Kokkos::parallel_for(Kokkos::TeamVectorRange(team,npacks), [&](const int& jpack) {
      dp_dry(jpack) = dp(jpack) * (1.0 - qv(jpack));
    });
    ColOps::column_scan<true>(team,nlevs,dp_dry,p_dry_int,ps0);
    team.team_barrier();
    ColOps::compute_midpoint_values(team,nlevs,p_dry_int,p_dry_mid);
    team.team_barrier();

    // Convert VTheta_dp->VTheta->Theta->T
    auto T      = ekat::subview(T_view,icol);
    auto T_prev = ekat::subview(T_prev_view,icol);

    Kokkos::parallel_for(Kokkos::TeamVectorRange(team,npacks),
                         [&](const int ilev) {
      // VTheta_dp->VTheta->Theta->T
      auto& T_val = T(ilev);
      T_val /= dp(ilev);
      T_val = PF::calculate_temperature_from_virtual_temperature(T_val,qv(ilev));
      T_val = PF::calculate_T_from_theta(T_val,p_mid(ilev));

      // Store T at end of the dyn timestep (to back out tendencies later)
      T_prev(ilev) = T_val;
    });
  }); //op()

  // Apply Rayleigh friction to update temperature and horiz_winds
  rayleigh_friction_apply(dt);

  if (has_energy_fixer()) {

      get_field_out("vapor_flux").deep_copy(0);
      get_field_out("ice_flux").deep_copy(0);
      get_field_out("water_flux").deep_copy(0);
      get_field_out("heat_flux").deep_copy(0);

  }; //if fixer

}//homme_post_proc

void HommeDynamics::
create_helper_field (const std::string& name,
                     const std::vector<FieldTag>& tags,
                     const std::vector<int>& dims,
                     const std::string& grid)
{
  FieldIdentifier id(name,FieldLayout{tags,dims},ekat::units::none,grid);

  const auto lt = id.get_layout().type();

  // Only request packed field for 3d quantities
  int pack_size = 1;
  if (lt==LayoutType::Scalar3D || lt==LayoutType::Vector3D || lt==LayoutType::Tensor3D) {
    pack_size = HOMMEXX_PACK_SIZE;
  }

  // Create the field. Init with NaN's, so we spot instances of uninited memory usage
  Field f(id);
  f.get_header().get_alloc_properties().request_allocation(pack_size);
  f.allocate_view();
  f.deep_copy(ekat::invalid<Real>());

  m_helper_fields[name] = f;
}

void HommeDynamics::
share_dyn_fields_with_instance (const std::vector<std::string>& names)
{
  for (const auto& n : names) {
    m_instance->add_dyn_field(n,m_helper_fields.at(n));
  }
}

void HommeDynamics::print_homme_info () {
  using TPF = ekat::TeamPolicyFactory<KT::ExeSpace>;

  const auto& c = Homme::Context::singleton();
  auto& params  = c.get<Homme::SimulationParams>();
  auto& tracers = c.get<Homme::Tracers>();

  const int qsize = tracers.num_tracers();

  const auto ncols = m_phys_grid->get_num_local_dofs();
  const auto nlevs = m_phys_grid->get_num_vertical_levels();
  const auto npacks= ekat::npack<Pack>(nlevs);

  const auto default_policy = TPF::get_default_team_policy(ncols,npacks);

  // Print homme's parameters, so user can see whether something wasn't set right.
  // TODO: make Homme::SimulationParams::print accept an ostream.
  std::stringstream msg;
  msg << "\n************** HOMMEXX SimulationParams **********************\n\n";
  msg << "   time_step_type: " << Homme::etoi(params.time_step_type) << "\n";
  msg << "   moisture: " << (params.use_moisture ? "moist" : "dry") << "\n";
  msg << "   remap_alg: " << Homme::etoi(params.remap_alg) << "\n";
  msg << "   test case: " << Homme::etoi(params.test_case) << "\n";
  msg << "   ftype: " << Homme::etoi(params.ftype) << "\n";
  msg << "   theta_adv_form: " << Homme::etoi(params.theta_adv_form) << "\n";
  msg << "   rsplit: " << params.rsplit << "\n";
  msg << "   qsplit: " << params.qsplit << "\n";
  msg << "   qsize: " << qsize << "\n";
  msg << "   limiter_option: " << params.limiter_option << "\n";
  msg << "   state_frequency: " << params.state_frequency << "\n";
  msg << "   dcmip16_mu: " << params.dcmip16_mu << "\n";
  msg << "   nu: " << params.nu << "\n";
  msg << "   nu_p: " << params.nu_p << "\n";
  msg << "   nu_q: " << params.nu_q << "\n";
  msg << "   nu_s: " << params.nu_s << "\n";
  msg << "   nu_top: " << params.nu_top << "\n";
  msg << "   nu_div: " << params.nu_div << "\n";
  msg << "   hypervis_order: " << params.hypervis_order << "\n";
  msg << "   hypervis_subcycle: " << params.hypervis_subcycle << "\n";
  msg << "   horiz_turb_subcycle: " << params.horiz_turb_subcycle << "\n";
  msg << "   hypervis_subcycle_tom: " << params.hypervis_subcycle_tom << "\n";
  msg << "   hypervis_scaling: " << params.hypervis_scaling << "\n";
  msg << "   nu_ratio1: " << params.nu_ratio1 << "\n";
  msg << "   nu_ratio2: " << params.nu_ratio2 << "\n";
  msg << "   use_cpstar: " << (params.use_cpstar ? "yes" : "no") << "\n";
  msg << "   transport_alg: " << params.transport_alg << "\n";
  msg << "   disable_diagnostics: " << (params.disable_diagnostics ? "yes" : "no") << "\n";
  msg << "   theta_hydrostatic_mode: " << (params.theta_hydrostatic_mode ? "yes" : "no") << "\n";
  msg << "   prescribed_wind: " << (params.prescribed_wind ? "yes" : "no") << "\n";

  msg << "\n************** General run info **********************\n\n";
  msg << "   ncols: " << ncols << "\n";
  msg << "   nlevs: " << nlevs << "\n";
  msg << "   npacks: " << npacks << "\n";
  msg << "   league_size: " << default_policy.league_size() << "\n";
  msg << "   team_size: " << default_policy.team_size() << "\n";
  msg << "   concurrent teams: " << KT::ExeSpace().concurrency() / default_policy.team_size() << "\n";

  // TODO: Replace with scale_factor and laplacian_rigid_factor when available.
  //msg << "   rearth: " << params.rearth << "\n";
  msg << "\n**********************************************************\n" << "\n";
  this->log(LogLevel::info,msg.str());
}

void HommeDynamics::restart_homme_state () {
  // Safety checks: RESTART fields *should* have been restarted (and therefore have a valid timestamp)
  for (auto& f : get_internal_fields()) {
    const auto& track = f.get_header().get_tracking();
    if (ekat::contains(track.get_groups_names(),"RESTART")) {
      auto ts = track.get_time_stamp();
      EKAT_REQUIRE_MSG(track.get_time_stamp().is_valid(),
          "Error! Found HommeDynamics internal field not restarted.\n"
          "  - field name: " + f.get_header().get_identifier().name() + "\n");
    }
  }

  using TPF = ekat::TeamPolicyFactory<KT::ExeSpace>;
  using PF = PhysicsFunctions<DefaultDevice>;

  const auto& pgn = m_phys_grid->name();

  // All internal fields should have been read from restart file.
  // We need to remap Q_dyn, v_dyn, w_dyn, T_dyn back to phys grid,
  // to handle the backing out of the tendencies
  // TODO: p2d remapper does not support subfields, so we need to create temps
  //       to remap v_dyn and w_dyn separately, then copy into v3d_prev.
  const auto& c = Homme::Context::singleton();
        auto& params = c.get<Homme::SimulationParams>();

  // NOTE: HommeModelInit has already restarted Homme's state. With the pg grid, it also
  //       already computed the T and uv at the end of the last Homme step.
  if (fv_phys_active()) {
    return;
  }

  const int nlevs = m_phys_grid->get_num_vertical_levels();
  const int ncols = m_phys_grid->get_num_local_dofs();
  const int npacks = ekat::npack<Pack>(nlevs);

  // NOTE: when restarting stuff like T_prev, and other "previous steps" quantities that HommeDynamics
  //       uses for tendencies calculation, we need to compute them in the *exact same way* as they
  //       were computed during the original simulation (in homme_post_process).
  //       E.g., we read vtheta_dp(dyn) from restart file, and need to recompute T_prev. Inside
  //       homme_post_process, we use qv(ref), but that's the qv obtained by remapping qv(dyn)
  //       to ref grid *right after homme ran*. Here, we cannot use qv(ref) as read from restart
  //       file, since that's qv(ref) *at the end of the timestep* in the original simulation.
  //       Therefore, we need to remap the end-of-homme-step qv from dyn to ref grid, and use that one.
  //       Another field we need is dp3d(ref), but Homme *CHECKS* that no other atm proc updates
  //       dp3d(ref), so the p2d-remapped value read from restart file *coincides* with the value at the end
  //       of the last Homme run. So we can safely recompute pressure using p2d(dp_dyn), with dp_dynread from restart.


  m_ic_remapper->register_field(m_helper_fields.at("FT_phys"),get_internal_field("vtheta_dp_dyn"));
  m_ic_remapper->register_field(m_helper_fields.at("FM_phys"),get_internal_field("v_dyn"));
  m_ic_remapper->register_field(get_field_out("pseudo_density",pgn),get_internal_field("dp3d_dyn"));
  auto qv_prev_ref = std::make_shared<Field>();
  auto Q_dyn = m_helper_fields.at("Q_dyn");
  if (params.ftype==Homme::ForcingAlg::FORCING_2) {
    auto Q_old = get_group_in("Q",pgn).monolithic_field();
    m_ic_remapper->register_field(Q_old,Q_dyn);

    // Grab qv_ref_old from Q_old
    *qv_prev_ref = Q_old.get_component(0);
  } else {
    using namespace ShortFieldTagsNames;

    // NOTE: we need the end-of-homme-step qv on the ref grid, to do the Theta->T conversion
    //       to compute T_prev *in the same way as homme_post_process* did in the original run.
    create_helper_field("qv_prev_phys",{COL,LEV},{ncols,nlevs},pgn);
    m_ic_remapper->register_field(m_helper_fields.at("qv_prev_phys"),Q_dyn.get_component(0));

    *qv_prev_ref = m_helper_fields.at("qv_prev_phys");
  }
  m_ic_remapper->registration_ends();
  m_ic_remapper->remap_bwd();
  m_ic_remapper = nullptr; // Can clean up the IC remapper now.

  // Now that we have dp_ref, we can recompute pressure
  update_pressure(m_phys_grid);

  // FT_ref contains vtheta_dp, so convert it to actual temperature
  auto T_prev_view = m_helper_fields.at("FT_phys").get_view<Pack**>();
  auto dp_view     = get_field_out("pseudo_density",pgn).get_view<const Pack**>();
  auto p_mid_view  = get_field_out("p_mid").get_view<Pack**>();
  auto qv_view     = qv_prev_ref->get_view<Pack**>();

  const auto policy = TPF::get_default_team_policy(ncols,npacks);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const KT::MemberType& team){
    const int icol = team.league_rank();

    auto p_mid = ekat::subview(p_mid_view,icol);
    auto qv    = ekat::subview(qv_view,icol);

    Kokkos::parallel_for(Kokkos::TeamVectorRange(team,npacks),
                         [&](const int& ilev) {
      // T_prev as of now contains vtheta_dp. Convert to temperature
      auto& T_prev = T_prev_view(icol,ilev);
      T_prev /= dp_view(icol,ilev);
      T_prev = PF::calculate_temperature_from_virtual_temperature(T_prev,qv(ilev));
      T_prev = PF::calculate_T_from_theta(T_prev,p_mid(ilev));
    });
  });
  Kokkos::fence();

  // Erase also qv_prev_phys (if we created it).
  m_helper_fields.erase("qv_prev_phys");

  // Update the time stamp of the fields we inited in here (to avoid triggering invalid output in IO)
  get_field_out("pseudo_density",pgn).get_header().get_tracking().update_time_stamp(start_of_step_ts());
}

void HommeDynamics::initialize_homme_state () {
  const auto& pgn = m_phys_grid->name();

  // Some Homme structures
  const auto& c = Homme::Context::singleton();
  auto& params = c.get<Homme::SimulationParams>();

  std::any_cast<int&>(*m_restart_extra_data["homme_nsteps"]) = c.get<Homme::TimeLevel>().nstep;

  // Update internal fields time stamp
  for (const auto& it : get_internal_fields()) {
    // Unfortunately, get_internal_fields() returns a list of const fields,
    // so grab the name and grid name, then call get_internal_field(name,grid)
    // It's a bit clunky, but not that bad
    const auto& name = it.get_header().get_identifier().name();
    const auto& grid = it.get_header().get_identifier().get_grid_name();
    auto& f = get_internal_field(name,grid);
    f.get_header().get_tracking().update_time_stamp(start_of_step_ts());
  }

  if (not fv_phys_active()) {
    // Forcings are computed as some version of "value coming in from AD
    // minus value at the end of last HommeDynamics run".
    // At the first time step, we don't have a value at the end of last
    // HommeDynamics run, so init with the initial conditions.
    // NOTE: for FM, we can't deep copy w_int, since w_int and FM_w
    //       have different number of levels. For u,v, we could, but
    //       we cannot (11/2021) subview 2 slices of FM together, so
    //       we'd need to also subview horiz_winds. Since we
    m_helper_fields.at("FT_phys").deep_copy(get_field_in("T_mid",pgn));
    m_helper_fields.at("FM_phys").deep_copy(get_field_out("horiz_winds",pgn));

    // Initialize p_mid/p_int
    update_pressure (m_phys_grid);
  }

  // Can clean up the IC remapper now.
  m_ic_remapper = nullptr;

  // Update the time stamp of the fields we inited in here (to avoid triggering invalid output in IO)
  get_field_out("pseudo_density",pgn).get_header().get_tracking().update_time_stamp(start_of_step_ts());
  get_internal_field("v_dyn").get_header().get_tracking().update_time_stamp(start_of_step_ts());
  get_internal_field("dp3d_dyn").get_header().get_tracking().update_time_stamp(start_of_step_ts());
  get_internal_field("ps_dyn").get_header().get_tracking().update_time_stamp(start_of_step_ts());
  get_internal_field("phis_dyn").get_header().get_tracking().update_time_stamp(start_of_step_ts());
  get_internal_field("vtheta_dp_dyn").get_header().get_tracking().update_time_stamp(start_of_step_ts());
  if (not params.theta_hydrostatic_mode) {
    get_internal_field("w_int_dyn").get_header().get_tracking().update_time_stamp(start_of_step_ts());
  }
}
// =========================================================================================
// Note: Any update to this routine will also need to be made within the homme_post_process
//       routine, which is responsible for updating the pressure every timestep.  There is a
//       TODO item to consolidate how we update the pressure during initialization and run, but
//       for now we have two locations where we do this.
void HommeDynamics::update_pressure(const std::shared_ptr<const AbstractGrid>& grid) {
  using ColOps = ColumnOps<DefaultDevice,Real>;
  using TPF = ekat::TeamPolicyFactory<KT::ExeSpace>;

  const auto ncols = grid->get_num_local_dofs();
  const auto nlevs = grid->get_num_vertical_levels();
  const auto npacks= ekat::npack<Pack>(nlevs);

  const auto& c = Homme::Context::singleton();
  const auto& hvcoord = c.get<Homme::HybridVCoord>();
  const auto ps0 = hvcoord.ps0 * hvcoord.hybrid_ai0;

  const auto& gn = grid->name();
  const auto dp_view    = get_field_out("pseudo_density",gn).get_view<Pack**>();
  const auto p_int_view = get_field_out("p_int",gn).get_view<Pack**>();
  const auto p_mid_view = get_field_out("p_mid",gn).get_view<Pack**>();

  const auto qv_view        = get_field_in("qv",gn).get_view<const Pack**>();
  const auto dp_dry_view    = get_field_out("pseudo_density_dry").get_view<Pack**>();
  const auto p_dry_int_view = get_field_out("p_dry_int").get_view<Pack**>();
  const auto p_dry_mid_view = get_field_out("p_dry_mid").get_view<Pack**>();

  const auto policy = TPF::get_thread_range_parallel_scan_team_policy(ncols,npacks);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const KT::MemberType& team) {
    const int& icol = team.league_rank();

    auto dp = ekat::subview(dp_view,icol);
    auto p_mid = ekat::subview(p_mid_view,icol);
    auto p_int = ekat::subview(p_int_view,icol);

    auto qv        = ekat::subview(qv_view,icol);
    auto dp_dry    = ekat::subview(dp_dry_view,icol);
    auto p_dry_mid = ekat::subview(p_dry_mid_view,icol);
    auto p_dry_int = ekat::subview(p_dry_int_view,icol);

    Kokkos::parallel_for(Kokkos::TeamVectorRange(team,npacks), [&](const int& jpack) {
      dp_dry(jpack) = dp(jpack) * (1.0 - qv(jpack));
    });

    ColOps::column_scan<true>(team,nlevs,dp,p_int,ps0);
    team.team_barrier();
    ColOps::compute_midpoint_values(team,nlevs,p_int,p_mid);

    ColOps::column_scan<true>(team,nlevs,dp_dry,p_dry_int,ps0);
    team.team_barrier();
    ColOps::compute_midpoint_values(team,nlevs,p_dry_int,p_dry_mid);
  });
}
// =========================================================================================

} // namespace scream
