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
#include "dynamics/homme/interface/scream_homme_interface.hpp"
#include "physics/share/physics_constants.hpp"
#include "share/util/scream_common_physics_functions.hpp"
#include "share/util/scream_column_ops.hpp"
#include "share/property_checks/field_lower_bound_check.hpp"

// Ekat includes
#include "ekat/ekat_assert.hpp"
#include "ekat/kokkos//ekat_subview_utils.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_scalar_traits.hpp"
#include "ekat/util/ekat_units.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "ekat/ekat_workspace.hpp"

namespace scream
{

HommeDynamics::HommeDynamics (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params), m_phys_grid_pgN(-1)
{
  // This class needs Homme's context, so register as a user
  HommeContextUser::singleton().add_user();

  ekat::any homme_nsteps;
  homme_nsteps.reset<int>(-1);
  m_restart_extra_data["homme_nsteps"] = homme_nsteps;

  if (!is_parallel_inited_f90()) {
    // While we're here, we can init homme's parallel session
    auto fcomm = MPI_Comm_c2f(comm.mpi_comm());
    init_parallel_f90 (fcomm);
  }

  // Set the log filename in the F90 interface
  const char* logname = m_atm_logger->get_logfile_name().c_str();
  set_homme_log_file_name_f90 (&logname);

  m_bfb_hash_nstep = 0;
  if (params.isParameter("BfbHash"))
    m_bfb_hash_nstep = std::max(0, params.get<int>("BfbHash"));
}

HommeDynamics::~HommeDynamics ()
{
  // This class is done with Homme. Remove from its users list
  HommeContextUser::singleton().remove_user();
}

void HommeDynamics::set_grids (const std::shared_ptr<const GridsManager> grids_manager)
{
  // Grab dynamics, physics, and physicsGLL grids
  const auto dgn = "Dynamics";
  m_dyn_grid = grids_manager->get_grid(dgn);
  m_phys_grid = grids_manager->get_grid("Physics");
  m_cgll_grid = grids_manager->get_grid("Physics GLL");

  fv_phys_set_grids();

  // Init prim structures
  // TODO: they should not be inited yet; should we error out if they are?
  //       I'm gonna say 'no', for now, cause it might be a pb with unit tests.
  if (!is_data_structures_inited_f90()) {
    // Set moisture in homme base on input file:
    const auto& moisture = m_params.get<std::string>("Moisture");
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
      "Error! The number of elements computed from the Dynamics grid num_dof()\n"
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
  auto pg_scalar3d_mid = m_phys_grid->get_3d_scalar_layout(true);
  auto pg_scalar3d_int = m_phys_grid->get_3d_scalar_layout(false);
  auto pg_vector3d_mid = m_phys_grid->get_3d_vector_layout(true,2);
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

  add_tracer<Updated >("qv", m_phys_grid, kg/kg, N);
  add_group<Updated>("tracers",pgn,N, Bundling::Required);

  if (fv_phys_active()) {
    // [CGLL ICs in pg2] Read CGLL IC data even though our in/out format is
    // pgN. I don't want to read these directly from the netcdf file because
    // doing so may conflict with additional IC mechanisms in the AD, e.g.,
    // init'ing a field to a constant.
    const auto& rgn = m_cgll_grid->name();
    auto rg_scalar2d     = m_cgll_grid->get_2d_scalar_layout();
    auto rg_scalar3d_mid = m_cgll_grid->get_3d_scalar_layout(true);
    auto rg_vector3d_mid = m_cgll_grid->get_3d_vector_layout(true,2);

    add_field<Required>("horiz_winds",   rg_vector3d_mid,m/s,   rgn,N);
    add_field<Required>("T_mid",         rg_scalar3d_mid,K,     rgn,N);
    add_field<Required>("ps",            rg_scalar2d    ,Pa,    rgn);
    add_field<Required>("phis",          rg_scalar2d    ,m2/s2, rgn);
    add_group<Required>("tracers",rgn,N, Bundling::Required, DerivationType::Import, "tracers", pgn);
    fv_phys_rrtmgp_active_gases_init(grids_manager);
    // This is needed for the dp_ref init in initialize_homme_state.
    add_field<Computed>("pseudo_density",rg_scalar3d_mid,Pa,    rgn,N);
  }

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
  add_internal_field (m_helper_fields.at("v_dyn").subfield(1,tl.n0,true));
  add_internal_field (m_helper_fields.at("vtheta_dp_dyn").subfield(1,tl.n0,true));
  add_internal_field (m_helper_fields.at("dp3d_dyn").subfield(1,tl.n0,true));
  add_internal_field (m_helper_fields.at("phi_int_dyn").subfield(1,tl.n0,true));
  add_internal_field (m_helper_fields.at("ps_dyn").subfield(1,tl.n0,true));
  add_internal_field (m_helper_fields.at("Qdp_dyn").subfield(1,tl.n0_qdp,true));
  if (not params.theta_hydrostatic_mode) {
    add_internal_field (m_helper_fields.at("w_int_dyn").subfield(1,tl.n0,true));
  }

  // The output manager pulls from the atm process fields. Add
  // helper fields for the case that a user request output.
  add_internal_field (m_helper_fields.at("omega_dyn"));
  add_internal_field (m_helper_fields.at("phis_dyn"));

  if (not fv_phys_active()) {
    // Dynamics backs out tendencies from the states, and passes those to Homme.
    // After Homme completes, we remap the updated state to the ref grid.  Thus,
    // is more convenient to use two different remappers: the pd remapper will
    // remap into Homme's forcing views, while the dp remapper will remap from
    // Homme's states.
    m_p2d_remapper = grids_manager->create_remapper(m_phys_grid,m_dyn_grid);
    m_d2p_remapper = grids_manager->create_remapper(m_dyn_grid,m_phys_grid);
  }

  // Create separate remapper for initial_conditions
  m_ic_remapper = grids_manager->create_remapper(m_cgll_grid,m_dyn_grid);
}

size_t HommeDynamics::requested_buffer_size_in_bytes() const
{
  using namespace Homme;

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

  return fbm.allocated_size()*sizeof(Real);
}

void HommeDynamics::init_buffers(const ATMBufferManager &buffer_manager)
{
  EKAT_REQUIRE_MSG(buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(),
                   "Error! Buffers size not sufficient.\n");

  using namespace Homme;
  auto& c = Context::singleton();
  auto& fbm  = c.get<FunctorsBuffersManager>();

  // Reset Homme buffer to use AD buffer memory.
  // Internally, homme will actually initialize its own buffers.
  EKAT_REQUIRE(buffer_manager.allocated_bytes()%sizeof(Real)==0); // Sanity check

  Real* mem = reinterpret_cast<Real*>(buffer_manager.get_memory());
  fbm.allocate(mem, buffer_manager.allocated_bytes()/sizeof(Real));
  mem += fbm.allocated_size();

  size_t used_mem = (mem - buffer_manager.get_memory())*sizeof(Real);
  EKAT_REQUIRE_MSG(used_mem==requested_buffer_size_in_bytes(),
                   "Error! Used memory != requested memory for HommeDynamics.");
}

void HommeDynamics::initialize_impl (const RunType run_type)
{
  const auto& dgn = m_dyn_grid->name();
  const auto& pgn = m_phys_grid->name();

  // Use common/shorter names for tracers.
  alias_group_in  ("tracers",pgn,"Q");
  alias_group_out ("tracers",pgn,"Q");

  // Grab handles of some Homme data structure
  const auto& c       = Homme::Context::singleton();
  const auto& params  = c.get<Homme::SimulationParams>();

  // Complete Homme prim_init1_xyz sequence
  prim_complete_init1_phase_f90 ();

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
      "Error! Someone other than Dynamics is trying to update the pseudo_density.\n");

  // The groups 'tracers' and 'tracers_mass_dyn' should contain the same fields
  EKAT_REQUIRE_MSG(not get_group_out("Q",pgn).m_info->empty(),
    "Error! There should be at least one tracer (qv) in the tracers group.\n");

  // Create remaining internal fields
  constexpr int NGP  = HOMMEXX_NP;
  const int nelem = m_dyn_grid->get_num_local_dofs()/(NGP*NGP);
  const int ncols = m_phys_grid->get_num_local_dofs();
  const int nlevs = m_dyn_grid->get_num_vertical_levels();
  assert(nlevs == m_dyn_grid->get_num_vertical_levels());
  assert(nlevs == m_cgll_grid->get_num_vertical_levels());
  assert(nlevs == m_phys_grid->get_num_vertical_levels());
  const int qsize = params.qsize;

  using namespace ShortFieldTagsNames;
  create_helper_field("FQ_dyn",{EL,CMP,GP,GP,LEV},{nelem,qsize,NGP,NGP,nlevs},dgn);
  create_helper_field("FT_dyn",{EL,    GP,GP,LEV},{nelem,      NGP,NGP,nlevs},dgn);
  create_helper_field("FM_dyn",{EL,CMP,GP,GP,LEV},{nelem,    3,NGP,NGP,nlevs},dgn);
  create_helper_field("Q_dyn" ,{EL,CMP,GP,GP,LEV},{nelem,qsize,NGP,NGP,nlevs},dgn);

  // Tendencies for temperature and momentum computed on physics grid
  create_helper_field("FT_phys",{COL,LEV},    {ncols,  nlevs},pgn);
  create_helper_field("FM_phys",{COL,CMP,LEV},{ncols,2,nlevs},pgn);

  // Unfortunately, Homme *does* use FM_z even in hydrostatic mode (missing ifdef).
  // Later, it computes diags on w, so if FM_w contains NaN's, repro sum in Homme
  // will crap out. To prevent that, set FM_z=0

  if (fv_phys_active()) {
    fv_phys_initialize_impl();
  } else {
    // Setup the p2d and d2p remappers
    m_p2d_remapper->registration_begins();
    m_d2p_remapper->registration_begins();

    // ftype==FORCING_0:
    //  1) remap Q_pgn->FQ_dyn
    //  2) compute FQ_dyn=dp3d*(FQ_dyn-Q_dyn_old)/dt
    // ftype!=FORCING_0:
    //  1) remap Q_pgn->FQ_dyn
    // Remap Q directly into FQ, tendency computed in pre_process step
    m_p2d_remapper->register_field(*get_group_out("Q",pgn).m_bundle,m_helper_fields.at("FQ_dyn"));
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
    m_d2p_remapper->register_field(m_helper_fields.at("Q_dyn"),*get_group_out("Q",pgn).m_bundle);
    m_d2p_remapper->register_field(m_helper_fields.at("omega_dyn"), get_field_out("omega"));

    m_p2d_remapper->registration_ends();
    m_d2p_remapper->registration_ends();
  }

  // Sets the scream views into the hommexx internal data structures
  init_homme_views ();

  // Import I.C. from the ref grid to the dyn grid.
  if (run_type==RunType::Initial) {
    initialize_homme_state ();
  } else {
    if (m_iop) {
      // We need to reload IOP data after restarting
      m_iop->read_iop_file_data(timestamp());
    }

    restart_homme_state ();
  }

  // Since we just inited them, ensure p_mid/p_int (dry and wet) timestamps are valid
  get_field_out("p_int")    .get_header().get_tracking().update_time_stamp(timestamp());
  get_field_out("p_mid")    .get_header().get_tracking().update_time_stamp(timestamp());
  get_field_out("p_dry_int").get_header().get_tracking().update_time_stamp(timestamp());
  get_field_out("p_dry_mid").get_header().get_tracking().update_time_stamp(timestamp());

  // Complete homme model initialization
  prim_init_model_f90 ();

  if (fv_phys_active()) {
    fv_phys_dyn_to_fv_phys(run_type == RunType::Restart);
    // [CGLL ICs in pg2] Remove the CGLL fields from the process. The AD has a
    // separate fvphyshack-based line to remove the whole CGLL FM. The intention
    // is to clear the view memory on the device, but I don't know if these two
    // steps are enough.
    //   We might end up needing these fields, anyway, for dynamics
    // output. Alternatively, it might be more efficient to make some of the
    // dynamics helper fields Computed and output on the dynamics grid. Worse
    // for I/O but better for device memory.
    const auto& rgn = m_cgll_grid->name();
    for (const auto& f : {"horiz_winds", "T_mid", "ps", "phis", "pseudo_density"})
      remove_field(f, rgn);
    remove_group("tracers", rgn);
    fv_phys_rrtmgp_active_gases_remap(run_type);
  }

  // Set up field property checks
  // Note: We are seeing near epsilon negative values in a handful of places.
  //       We figured out the source of those values in Homme (see #1842),
  //       but they should be "benign" neg values (due to roundoffs), which we
  //       can safely clip to 0..
  using Interval = FieldWithinIntervalCheck;
  using LowerBound = FieldLowerBoundCheck;

  add_postcondition_check<LowerBound>(*get_group_out("Q",pgn).m_bundle,m_phys_grid,0,true);
  add_postcondition_check<Interval>(get_field_out("T_mid",pgn),m_phys_grid,100.0, 500.0,false);
  add_postcondition_check<Interval>(get_field_out("horiz_winds",pgn),m_phys_grid,-400.0, 400.0,false);
  add_postcondition_check<Interval>(get_field_out("ps"),m_phys_grid,40000.0, 120000.0,false);

  // Initialize Rayleigh friction variables
  rayleigh_friction_init();
}

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

    if (m_bfb_hash_nstep > 0 && timestamp().get_num_steps() % m_bfb_hash_nstep == 0)
      print_fast_global_state_hash("Hommexx");

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

    // Update nstep in the restart extra data, so it can be written to restart if needed.
    const auto& tl = c.get<Homme::TimeLevel>();
    auto& nstep = ekat::any_cast<int>(m_restart_extra_data["homme_nsteps"]);
    nstep = tl.nstep;

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

  if (group.m_info->m_group_name=="tracers") {
    // Set runtime number of tracers in Homme
    auto& params = c.get<Homme::SimulationParams>();
    const int qsize = group.m_info->size();
    params.qsize = qsize;           // Set in the CXX data structure
    set_homme_param("qsize",qsize); // Set in the F90 module
    tracers.init(tracers.num_elems(),qsize);
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
    m_p2d_remapper->remap(true);
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

  if (m_iop) {
    apply_iop_forcing(dt);
  }

  if (fv_phys_active()) {
    fv_phys_post_process();
    // Apply Rayleigh friction to update temperature and horiz_winds
    rayleigh_friction_apply(dt);

    return;
  }

  // Remap outputs to ref grid
  m_d2p_remapper->remap(true);

  using ColOps = ColumnOps<DefaultDevice,Real>;
  using PF = PhysicsFunctions<DefaultDevice>;

  // Convert VTheta_dp->T, store T,uv, and possibly w in FT, FM,
  // compute p_int on ref grid.
  const auto dp_view = get_field_out("pseudo_density").get_view<Pack**>();
  const auto p_mid_view = get_field_out("p_mid").get_view<Pack**>();
  const auto p_int_view = get_field_out("p_int").get_view<Pack**>();
  const auto dp_dry_view    = get_field_out("pseudo_density_dry").get_view<Pack**>();
  const auto p_dry_int_view = get_field_out("p_dry_int").get_view<Pack**>();
  const auto p_dry_mid_view = get_field_out("p_dry_mid").get_view<Pack**>();
  const auto Q_view   = get_group_out("Q",pgn).m_bundle->get_view<Pack***>();

  const auto T_view  = get_field_out("T_mid").get_view<Pack**>();
  const auto T_prev_view = m_helper_fields.at("FT_phys").get_view<Pack**>();

  m_helper_fields.at("FM_phys").deep_copy(get_field_out("horiz_winds"));

  const auto ncols = m_phys_grid->get_num_local_dofs();
  const auto nlevs = m_phys_grid->get_num_vertical_levels();
  const auto npacks= ekat::npack<Pack>(nlevs);

  using ESU = ekat::ExeSpaceUtils<KT::ExeSpace>;
  const auto policy = ESU::get_thread_range_parallel_scan_team_policy(ncols,npacks);

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
  });

  // Apply Rayleigh friction to update temperature and horiz_winds
  rayleigh_friction_apply(dt);
}

void HommeDynamics::
create_helper_field (const std::string& name,
                     const std::vector<FieldTag>& tags,
                     const std::vector<int>& dims,
                     const std::string& grid) {
  using namespace ekat::units;
  FieldIdentifier id(name,FieldLayout{tags,dims},Units::nondimensional(),grid);

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
  f.deep_copy(ekat::ScalarTraits<Real>::invalid());

  m_helper_fields[name] = f;
}

void HommeDynamics::init_homme_views () {

  const auto& c = Homme::Context::singleton();
  auto& params  = c.get<Homme::SimulationParams>();
  auto& state   = c.get<Homme::ElementsState>();
  auto& derived = c.get<Homme::ElementsDerivedState>();
  auto& tracers = c.get<Homme::Tracers>();
  auto& forcing = c.get<Homme::ElementsForcing>();

  constexpr int NGP  = HOMMEXX_NP;
  constexpr int NTL  = HOMMEXX_NUM_TIME_LEVELS;
  constexpr int QTL  = HOMMEXX_Q_NUM_TIME_LEVELS;
  constexpr int QSZ  = HOMMEXX_QSIZE_D;
  constexpr int NVL  = HOMMEXX_NUM_LEV;
  constexpr int NVLI = HOMMEXX_NUM_LEV_P;

  const int nelem = m_dyn_grid->get_num_local_dofs()/(NGP*NGP);
  const int qsize = tracers.num_tracers();

  const auto ncols = m_phys_grid->get_num_local_dofs();
  const auto nlevs = m_phys_grid->get_num_vertical_levels();
  const auto npacks= ekat::npack<Pack>(nlevs);

  using ESU = ekat::ExeSpaceUtils<KT::ExeSpace>;
  const auto default_policy = ESU::get_default_team_policy(ncols,npacks);

  // Print homme's parameters, so user can see whether something wasn't set right.
  // TODO: make Homme::SimulationParams::print accept an ostream.
  std::stringstream msg;
  msg << "\n************** HOMMEXX SimulationParams **********************\n\n";
  msg << "   time_step_type: " << Homme::etoi(params.time_step_type) << "\n";
  msg << "   moisture: " << (params.moisture==Homme::MoistDry::DRY ? "dry" : "moist") << "\n";
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

  // ------------ Set views in Homme ------------- //
  // Velocity
  auto v_in = m_helper_fields.at("v_dyn").get_view<Homme::Scalar*[NTL][2][NP][NP][NVL]>();
  using v_type = std::remove_reference<decltype(state.m_v)>::type;
  state.m_v = v_type (v_in.data(),nelem);

  // Virtual potential temperature
  auto vtheta_in = m_helper_fields.at("vtheta_dp_dyn").get_view<Homme::Scalar*[NTL][NP][NP][NVL]>();
  using vtheta_type = std::remove_reference<decltype(state.m_vtheta_dp)>::type;
  state.m_vtheta_dp = vtheta_type(vtheta_in.data(),nelem);

  // Geopotential
  auto phi_in = m_helper_fields.at("phi_int_dyn").get_view<Homme::Scalar*[NTL][NP][NP][NVLI]>();
  using phi_type = std::remove_reference<decltype(state.m_phinh_i)>::type;
  state.m_phinh_i = phi_type(phi_in.data(),nelem);

  // Vertical velocity
  auto w_in = m_helper_fields.at("w_int_dyn").get_view<Homme::Scalar*[NTL][NP][NP][NVLI]>();
  using w_type = std::remove_reference<decltype(state.m_w_i)>::type;
  state.m_w_i = w_type(w_in.data(),nelem);

  // Pseudo-density
  auto dp3d_in = m_helper_fields.at("dp3d_dyn").template get_view<Homme::Scalar*[NTL][NP][NP][NVL]>();
  using dp3d_type = std::remove_reference<decltype(state.m_dp3d)>::type;
  state.m_dp3d = dp3d_type(dp3d_in.data(),nelem);

  // Surface pressure
  auto ps_in = m_helper_fields.at("ps_dyn").template get_view<Real*[NTL][NP][NP]>();
  using ps_type = std::remove_reference<decltype(state.m_ps_v)>::type;
  state.m_ps_v = ps_type(ps_in.data(),nelem);

  // Vertical pressure velocity
  auto omega_in = m_helper_fields.at("omega_dyn").template get_view<Homme::Scalar*[NP][NP][NVL]>();
  using omega_type = std::remove_reference<decltype(derived.m_omega_p)>::type;
  derived.m_omega_p = omega_type(omega_in.data(),nelem);

  // Tracers mixing ratio
  auto q_in = m_helper_fields.at("Q_dyn").template get_view<Homme::Scalar**[NP][NP][NVL]>();
  using q_type = std::remove_reference<decltype(tracers.Q)>::type;
  tracers.Q = q_type(q_in.data(),nelem,qsize);

  // Tracers mass
  auto qdp_in = m_helper_fields.at("Qdp_dyn").template get_view<Homme::Scalar*[QTL][QSZ][NP][NP][NVL]>();
  using qdp_type = std::remove_reference<decltype(tracers.qdp)>::type;
  tracers.qdp = qdp_type(qdp_in.data(),nelem);

  // Tracers forcing
  auto fq_in = m_helper_fields.at("FQ_dyn").template get_view<Homme::Scalar**[NP][NP][NVL]>();
  using fq_type = std::remove_reference<decltype(tracers.fq)>::type;
  tracers.fq = fq_type(fq_in.data(),nelem,qsize);

  // Temperature forcing
  auto ft_in = m_helper_fields.at("FT_dyn").template get_view<Homme::Scalar*[NP][NP][NVL]>();
  using ft_type = std::remove_reference<decltype(forcing.m_ft)>::type;
  forcing.m_ft = ft_type(ft_in.data(),nelem);

  // Momentum forcing
  auto fm_in = m_helper_fields.at("FM_dyn").template get_view<Homme::Scalar**[NP][NP][NVL]>();
  using fm_type = std::remove_reference<decltype(forcing.m_fm)>::type;
  forcing.m_fm = fm_type(fm_in.data(),nelem);

  // Homme has 3 components for FM, but the 3d (the omega forcing) is not computed
  // by EAMxx, so we set FM(3)=0 right away
  m_helper_fields.at("FM_dyn").get_component(2).deep_copy(0);
}

void HommeDynamics::restart_homme_state () {
  // Safety checks: internal fields *should* have been restarted (and therefore have a valid timestamp)
  for (auto& f : get_internal_fields()) {
    auto ts = f.get_header().get_tracking().get_time_stamp();
    EKAT_REQUIRE_MSG(ts.is_valid(),
        "Error! Found HommeDynamics internal field not restarted.\n"
        "  - field name: " + f.get_header().get_identifier().name() + "\n");
  }

  using ESU = ekat::ExeSpaceUtils<KT::ExeSpace>;
  using PF = PhysicsFunctions<DefaultDevice>;

  const auto& dgn = m_dyn_grid->name();
  const auto& pgn = m_phys_grid->name();

  // All internal fields should have been read from restart file.
  // We need to remap Q_dyn, v_dyn, w_dyn, T_dyn back to phys grid,
  // to handle the backing out of the tendencies
  // TODO: p2d remapper does not support subfields, so we need to create temps
  //       to remap v_dyn and w_dyn separately, then copy into v3d_prev.
  const auto& c = Homme::Context::singleton();
        auto& params = c.get<Homme::SimulationParams>();
        auto& tl = c.get<Homme::TimeLevel>();

  // For BFB restarts, set nstep counter in Homme's TimeLevel to match the restarted value.
  const auto& nstep = ekat::any_ptr_cast<int>(m_restart_extra_data["homme_nsteps"]);
  tl.nstep = *nstep;
  set_homme_param("num_steps",*nstep);

  constexpr int NGP = HOMMEXX_NP;
  const int nlevs = m_phys_grid->get_num_vertical_levels();
  const int ncols = m_phys_grid->get_num_local_dofs();
  const int nelem = m_dyn_grid->get_num_local_dofs() / (NGP*NGP);
  const int npacks = ekat::npack<Pack>(nlevs);
  const int qsize = params.qsize;

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

  // Copy all restarted dyn states on all timelevels.
  copy_dyn_states_to_all_timelevels ();

  if (params.theta_hydrostatic_mode) {
    // Nothing read from restart file for w_int, but Homme still does some global reduction on w_int when
    // printing the state, so we need to make sure it doesn't contain NaNs
    m_helper_fields.at("w_int_dyn").deep_copy(0);
  }

  // Restart end-of-homme-step Q as Qdp/dp. That's what Homme does at the end of the timestep,
  // and by writing/loading only Qdp in the restart file, we save space.
  auto Qdp_dyn_view = get_internal_field("Qdp_dyn",dgn).get_view<Pack*****>();
  auto Q_dyn_view = m_helper_fields.at("Q_dyn").get_view<Pack*****>();
  auto dp_dyn_view = get_internal_field("dp3d_dyn",dgn).get_view<Pack****>();
  Kokkos::parallel_for(Kokkos::RangePolicy<>(0,nelem*qsize*NGP*NGP*npacks),
                       KOKKOS_LAMBDA (const int idx) {
    const int ie =  idx / (qsize*NGP*NGP*npacks);
    const int iq = (idx / (NGP*NGP*npacks)) % qsize;
    const int ip = (idx / (NGP*npacks)) % NGP;
    const int jp = (idx / npacks) % NGP;
    const int k  =  idx % npacks;
    Q_dyn_view(ie,iq,ip,jp,k) = Qdp_dyn_view(ie,iq,ip,jp,k) / dp_dyn_view(ie,ip,jp,k);
  });

  if (fv_phys_active()) {
    m_ic_remapper = nullptr;
    return;
  }

  m_ic_remapper->registration_begins();
  m_ic_remapper->register_field(m_helper_fields.at("FT_phys"),get_internal_field("vtheta_dp_dyn"));
  m_ic_remapper->register_field(m_helper_fields.at("FM_phys"),get_internal_field("v_dyn"));
  m_ic_remapper->register_field(get_field_out("pseudo_density",pgn),get_internal_field("dp3d_dyn"));
  auto qv_prev_ref = std::make_shared<Field>();
  auto Q_dyn = m_helper_fields.at("Q_dyn");
  if (params.ftype==Homme::ForcingAlg::FORCING_2) {
    auto Q_old = *get_group_in("Q",pgn).m_bundle;
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
  m_ic_remapper->remap(/*forward = */false);
  m_ic_remapper = nullptr; // Can clean up the IC remapper now.

  // Now that we have dp_ref, we can recompute pressure
  update_pressure(m_phys_grid);

  // FT_ref contains vtheta_dp, so convert it to actual temperature
  auto T_prev_view = m_helper_fields.at("FT_phys").get_view<Pack**>();
  auto dp_view     = get_field_out("pseudo_density",pgn).get_view<const Pack**>();
  auto p_mid_view  = get_field_out("p_mid").get_view<Pack**>();
  auto qv_view     = qv_prev_ref->get_view<Pack**>();

  const auto policy = ESU::get_default_team_policy(ncols,npacks);
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
}

void HommeDynamics::initialize_homme_state () {
  // Some types
  using ColOps = ColumnOps<DefaultDevice,Real>;
  using PF = PhysicsFunctions<DefaultDevice>;
  using ESU = ekat::ExeSpaceUtils<KT::ExeSpace>;
  using EOS = Homme::EquationOfState;

  const auto& rgn = m_cgll_grid->name();

  // Some Homme structures
  const auto& c = Homme::Context::singleton();
  auto& params = c.get<Homme::SimulationParams>();
  const auto& hvcoord = c.get<Homme::HybridVCoord>();

  // Some extents
  const auto ncols = m_cgll_grid->get_num_local_dofs();
  const auto nlevs = m_cgll_grid->get_num_vertical_levels();
  constexpr int NGP = HOMMEXX_NP;
  const int nelem = m_dyn_grid->get_num_local_dofs()/(NGP*NGP);
  const int qsize = params.qsize;
  const int npacks_mid = ekat::npack<Pack>(nlevs);
  const int npacks_int = ekat::npack<Pack>(nlevs+1);

  // Bootstrap dp on phys grid, and let the ic remapper transfer dp on dyn grid
  // NOTE: HybridVCoord already stores hyai and hybi deltas as packed views,
  //       but it used KokkosKernels packs, which are incompatible with ekat::Pack.
  //       If Homme switched to ekat::Pack, you could do the loop below with packs.
  const auto ps0 = hvcoord.ps0;
  const auto dp_ref = get_field_out("pseudo_density",rgn).get_view<Real**>();
  const auto ps_ref = get_field_in("ps",rgn).get_view<const Real*>();
  const auto hyai = hvcoord.hybrid_ai;
  const auto hybi = hvcoord.hybrid_bi;
  const auto policy_dp = ESU::get_default_team_policy(ncols, nlevs);
  Kokkos::parallel_for(policy_dp, KOKKOS_LAMBDA (const KT::MemberType& team) {
    const int icol = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team,nlevs),
                        [&](const int ilev) {
       dp_ref(icol,ilev) = (hyai(ilev+1)-hyai(ilev))*ps0
                         + (hybi(ilev+1)-hybi(ilev))*ps_ref(icol);
    });
    team.team_barrier();
  });

  // Import IC from ref grid to dyn grid
  // NOTE: if/when PD remapper supports remapping directly to/from subfields,
  //       you can use get_internal_field (which have a single time slice) rather than
  //       the helper fields (which have NTL time slices).
  m_ic_remapper->registration_begins();
  m_ic_remapper->register_field(get_field_in("horiz_winds",rgn),get_internal_field("v_dyn"));
  m_ic_remapper->register_field(get_field_out("pseudo_density",rgn),get_internal_field("dp3d_dyn"));
  m_ic_remapper->register_field(get_field_in("ps",rgn),get_internal_field("ps_dyn"));
  m_ic_remapper->register_field(get_field_in("phis",rgn),m_helper_fields.at("phis_dyn"));
  m_ic_remapper->register_field(get_field_in("T_mid",rgn),get_internal_field("vtheta_dp_dyn"));
  m_ic_remapper->register_field(*get_group_in("tracers",rgn).m_bundle,m_helper_fields.at("Q_dyn"));
  m_ic_remapper->registration_ends();
  m_ic_remapper->remap(true);

  // Wheter w_int is computed or not, Homme still does some global reduction on w_int when
  // printing the state, so we need to make sure it doesn't contain NaNs
  m_helper_fields.at("w_int_dyn").deep_copy(0.0);

  // Homme states
  const auto dp_view  = m_helper_fields.at("dp3d_dyn").get_view<Pack*****>();
  const auto vth_view = m_helper_fields.at("vtheta_dp_dyn").get_view<Pack*****>();
  const auto Q_view   = m_helper_fields.at("Q_dyn").get_view<Pack*****>();

  // State time slices
  const auto& tl = c.get<Homme::TimeLevel>();
  const int n0  = tl.n0;
  const int n0_qdp  = tl.n0_qdp;

  ekat::any_cast<int>(m_restart_extra_data["homme_nsteps"]) = tl.nstep;

  const auto phis_dyn_view = m_helper_fields.at("phis_dyn").get_view<const Real***>();
  const auto phi_int_view = m_helper_fields.at("phi_int_dyn").get_view<Pack*****>();
  const auto hyai0 = hvcoord.hybrid_ai0;
  // Need two temporaries, for pi_mid and pi_int
  const auto policy = ESU::get_thread_range_parallel_scan_team_policy(nelem*NGP*NGP,npacks_mid);
  WorkspaceMgr wsm(npacks_int,2,policy);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const KT::MemberType& team) {
    const int ie  =  team.league_rank() / (NGP*NGP);
    const int igp = (team.league_rank() / NGP) % NGP;
    const int jgp =  team.league_rank() % NGP;

    // Compute dp assuming hydrostatic initial state
    auto dp = ekat::subview(dp_view,ie,n0,igp,jgp);

    // Compute p_mid
    auto ws = wsm.get_workspace(team);
    ekat::Unmanaged<WorkspaceMgr::view_1d<Pack> > p_int, p_mid;
    ws.take_many_and_reset<2>({"p_int", "p_mid"}, {&p_int, &p_mid});

    ColOps::column_scan<true>(team,nlevs,dp,p_int,ps0*hyai0);
    team.team_barrier();
    ColOps::compute_midpoint_values(team,nlevs,p_int,p_mid);
    team.team_barrier();

    // Convert T->Theta->VTheta->VTheta*dp in place
    auto T      = ekat::subview(vth_view,ie,n0,igp,jgp);
    auto vTh_dp = ekat::subview(vth_view,ie,n0,igp,jgp);
    auto qv     = ekat::subview(Q_view,ie,0,igp,jgp);
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team,npacks_mid),
                         [&](const int ilev) {
      const auto th = PF::calculate_theta_from_T(T(ilev),p_mid(ilev));
      vTh_dp(ilev) = PF::calculate_virtual_temperature(th,qv(ilev))*dp(ilev);
    });
    team.team_barrier();

    // Init geopotential
    auto dphi   = [&](const int ilev)->Pack {
      return EOS::compute_dphi(vTh_dp(ilev), p_mid(ilev));
    };
    auto phi_int = ekat::subview(phi_int_view,ie,n0,igp,jgp);
    ColOps::column_scan<false>(team,nlevs,dphi,phi_int,phis_dyn_view(ie,igp,jgp));
  });

  // Update internal fields time stamp
  for (const auto& it : get_internal_fields()) {
    // Unfortunately, get_internal_fields() returns a list of const fields,
    // so grab the name and grid name, then call get_internal_field(name,grid)
    // It's a bit clunky, but not that bad
    const auto& name = it.get_header().get_identifier().name();
    const auto& grid = it.get_header().get_identifier().get_grid_name();
    auto& f = get_internal_field(name,grid);
    f.get_header().get_tracking().update_time_stamp(timestamp());
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
    const auto& pgn = m_phys_grid->name();
    m_helper_fields.at("FT_phys").deep_copy(get_field_in("T_mid",pgn));
    m_helper_fields.at("FM_phys").deep_copy(get_field_out("horiz_winds",pgn));
  }

  // For initial runs, it's easier to prescribe IC for Q, and compute Qdp = Q*dp
  auto& tracers = c.get<Homme::Tracers>();
  const auto qdp = tracers.qdp;
  const auto q   = tracers.Q;
  const auto dp  = c.get<Homme::ElementsState>().m_dp3d;
  Kokkos::parallel_for(Kokkos::RangePolicy<>(0,nelem*qsize*NGP*NGP*npacks_mid),
                       KOKKOS_LAMBDA (const int idx) {
    const int ie =  idx / (qsize*NGP*NGP*npacks_mid);
    const int iq = (idx / (NGP*NGP*npacks_mid)) % qsize;
    const int ip = (idx / (NGP*npacks_mid)) % NGP;
    const int jp = (idx / npacks_mid) % NGP;
    const int k  = idx % npacks_mid;

    qdp(ie,n0_qdp,iq,ip,jp,k) = q(ie,iq,ip,jp,k) * dp(ie,n0,ip,jp,k);
  });

  if (not fv_phys_active()) {
    // Initialize p_mid/p_int
    update_pressure (m_phys_grid);
  }

  // If "Instant" averaging type is used for output,
  // an initial output is performed before AD processes
  // are run. If omega_dyn output is requested, it will
  // not have valid computed values for this initial
  // output. Set to zero avoid potential FPE.
  get_internal_field("omega_dyn").deep_copy(0);

  // Copy IC states on all timelevel slices
  copy_dyn_states_to_all_timelevels ();

  // Can clean up the IC remapper now.
  m_ic_remapper = nullptr;
}
// =========================================================================================
void HommeDynamics::
copy_dyn_states_to_all_timelevels () {
  const auto& c = Homme::Context::singleton();

  // State time slices
  const int nm1     = c.get<Homme::TimeLevel>().nm1;
  const int n0      = c.get<Homme::TimeLevel>().n0;
  const int np1     = c.get<Homme::TimeLevel>().np1;
  const int n0_qdp  = c.get<Homme::TimeLevel>().n0_qdp;
  const int np1_qdp = c.get<Homme::TimeLevel>().np1_qdp;

  // States
  const auto dp3d      = m_helper_fields.at("dp3d_dyn");
  const auto ps        = m_helper_fields.at("ps_dyn");
  const auto v         = m_helper_fields.at("v_dyn");
  const auto w_i       = m_helper_fields.at("w_int_dyn");
  const auto phinh_i   = m_helper_fields.at("phi_int_dyn");
  const auto vtheta_dp = m_helper_fields.at("vtheta_dp_dyn");
  const auto qdp       = m_helper_fields.at("Qdp_dyn");

  // Note: it might be somewhat faster to do a single parallel region,
  //       rather than 13 deep copies. However, this is much clearer
  //       to read, and since this method is called only during
  //       initialization, we prefer a more readable method
  dp3d.subfield(1,nm1).deep_copy(dp3d.subfield(1,n0));
  dp3d.subfield(1,np1).deep_copy(dp3d.subfield(1,n0));
  ps.subfield(1,nm1).deep_copy(ps.subfield(1,n0));
  ps.subfield(1,np1).deep_copy(ps.subfield(1,n0));
  v.subfield(1,nm1).deep_copy(v.subfield(1,n0));
  v.subfield(1,np1).deep_copy(v.subfield(1,n0));
  w_i.subfield(1,nm1).deep_copy(w_i.subfield(1,n0));
  w_i.subfield(1,np1).deep_copy(w_i.subfield(1,n0));
  phinh_i.subfield(1,nm1).deep_copy(phinh_i.subfield(1,n0));
  phinh_i.subfield(1,np1).deep_copy(phinh_i.subfield(1,n0));
  vtheta_dp.subfield(1,nm1).deep_copy(vtheta_dp.subfield(1,n0));
  vtheta_dp.subfield(1,np1).deep_copy(vtheta_dp.subfield(1,n0));
  qdp.subfield(1,np1_qdp).deep_copy(qdp.subfield(1,n0_qdp));
}

// =========================================================================================
// Note: Any update to this routine will also need to be made within the homme_post_process
//       routine, which is responsible for updating the pressure every timestep.  There is a
//       TODO item to consolidate how we update the pressure during initialization and run, but
//       for now we have two locations where we do this.
void HommeDynamics::update_pressure(const std::shared_ptr<const AbstractGrid>& grid) {
  using ColOps = ColumnOps<DefaultDevice,Real>;

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

  const auto qv_view        = get_field_in("qv").get_view<const Pack**>();
  const auto dp_dry_view    = get_field_out("pseudo_density_dry").get_view<Pack**>();
  const auto p_dry_int_view = get_field_out("p_dry_int").get_view<Pack**>();
  const auto p_dry_mid_view = get_field_out("p_dry_mid").get_view<Pack**>();

  using ESU = ekat::ExeSpaceUtils<KT::ExeSpace>;
  const auto policy = ESU::get_thread_range_parallel_scan_team_policy(ncols,npacks);
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
