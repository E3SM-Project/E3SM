#include "atmosphere_dynamics.hpp"
#include <string>

// HOMMEXX Includes
#include "Context.hpp"
#include "Elements.hpp"
#include "SimulationParams.hpp"
#include "ElementsForcing.hpp"
#include "EulerStepFunctor.hpp"
#include "Diagnostics.hpp"
#include "DirkFunctor.hpp"
#include "ForcingFunctor.hpp"
#include "CaarFunctor.hpp"
#include "VerticalRemapManager.hpp"
#include "HyperviscosityFunctor.hpp"
#include "ExecSpaceDefs.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "ekat/ekat_workspace.hpp"
#include "TimeLevel.hpp"
#include "Tracers.hpp"
#include "Types.hpp"

// Scream includes
#include "dynamics/homme/homme_dimensions.hpp"
#include "dynamics/homme/homme_dynamics_helpers.hpp"
#include "dynamics/homme/interface/scream_homme_interface.hpp"
#include "physics/share/physics_constants.hpp"
#include "share/util/scream_common_physics_functions.hpp"
#include "share/util//scream_column_ops.hpp"
#include "share/field/field_property_checks/field_lower_bound_check.hpp"

// Ekat includes
#include "ekat/ekat_assert.hpp"
#include "ekat/kokkos//ekat_subview_utils.hpp"

namespace scream
{

HommeDynamics::HommeDynamics (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
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
  using FL  = FieldLayout;

  constexpr int NGP  = HOMMEXX_NP;
  constexpr int NVL  = HOMMEXX_NUM_PHYSICAL_LEV;
  constexpr int NTL  = HOMMEXX_NUM_TIME_LEVELS;

  // Some units
  const auto nondim = Units::nondimensional();
  const auto rho = kg/(m*m*m);
  auto Q = kg/kg;
  Q.set_string("kg/kg");

  const auto dgn = "Dynamics";
  m_dyn_grid = grids_manager->get_grid(dgn);
  m_ref_grid = grids_manager->get_reference_grid();
  const auto& rgn = m_ref_grid->name();

  const int ne = m_dyn_grid->get_num_local_dofs()/(NGP*NGP);
  const int ncols = m_ref_grid->get_num_local_dofs();

  // Sanity check for the grid. This should *always* pass, since Homme builds the grids
  EKAT_REQUIRE_MSG(get_num_local_elems_f90()==ne,
      "Error! The number of elements computed from the Dynamics grid num_dof()\n"
      "       does not match the number of elements internal in Homme.\n");

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

  // Note: qv is needed to make sure Q is not empty (dyn needs qv to transform T<->Theta),
  constexpr int N = sizeof(Homme::Scalar) / sizeof(Real);
  add_field<Updated> ("horiz_winds",   FL({COL,CMP, LEV},{ncols,2,NVL}),  m/s,   rgn,N);
  add_field<Updated> ("T_mid",         FL({COL,     LEV},{ncols,  NVL}),  K,     rgn,N);
  add_field<Updated> ("w_int",         FL({COL,    ILEV},{ncols,  NVL+1}),m/s,   rgn,N);
  add_field<Updated> ("phi_int",       FL({COL,    ILEV},{ncols,  NVL+1}),Pa/rho,rgn,N);
  add_field<Updated> ("pseudo_density",FL({COL,     LEV},{ncols,  NVL}),  Pa,    rgn,N);
  add_field<Updated> ("ps",            FL({COL         },{ncols      }),  Pa,    rgn);
  add_field<Required>("qv",            FL({COL,     LEV},{ncols,  NVL}),  Q,     rgn,"tracers",N);
  add_field<Computed>("p_int",         FL({COL,    ILEV},{ncols,  NVL+1}),Pa,    rgn,N);
  add_field<Computed>("p_mid",         FL({COL,     LEV},{ncols,  NVL}),  Pa,    rgn,N);

  // These are used to back out tendencies
  add_field<Updated>("T_mid_prev",      FL({COL,LEV},    {ncols,  NVL}),  K/s,    rgn,N);
  add_field<Updated>("horiz_winds_prev",FL({COL,CMP,LEV},{ncols,2,NVL}),  m/(s*s),rgn,N);
  add_field<Updated>("w_int_prev",      FL({COL,LEV},    {ncols,  NVL+1}),m/(s*s),rgn,N);

  // Input-output groups.
  // Note: HommeDynamics needs Q (tracers) and FQ (tracers forcing). Depending on ftype
  // (the forcing type), the latter may be computed as FQ=(Qnew-Qold)/dt, where
  // Qnew is the new value of Q (coming from physics), and Qold is the value of Q
  // at the end of the last dyn step.
  // So we register 2 groups:
  //   - Q on ref grid: every atm proc interacts with other atm procs on the ref grid,
  //     so we *must* use this as the main input/output.
  //   - Q on dyn grid: what Homme will use for its computations. Also, this will store
  //     the end-of-step values for the tracers. Homme will make sure that *no other atm proc*
  //     updates this field (or any subfield) since it's needed for exact restart.
  // The forcing FQ is instead stored as an internal field, without storing it in the field manager.
  GroupRequest tracers_ref("tracers",m_ref_grid->name(),N, Bundling::Required);
  GroupRequest tracers_dyn("tracers",m_dyn_grid->name(),N, Bundling::Required,
                           DerivationType::Import, "tracers", m_ref_grid->name());
  add_group<Updated>(tracers_ref);
  add_group<Updated>(tracers_dyn);

  // NOTE: since these fields are internal, and they don't end up in the FieldManager,
  //       we don't care about units.
  // Note: the dyn grid names do not follow the FieldManager conventions,
  //       but rather the var names in Homme. That's ok, since these fields are
  //       not exposed to the ourside world. We use fields (rather than raw views)
  //       to leverage the scream remapper structures.
  const int nelem = m_dyn_grid->get_num_local_dofs()/(NGP*NGP);
  const int nlevs = m_dyn_grid->get_num_vertical_levels();

  create_internal_field("FM",       {EL,   CMP,GP,GP,LEV}, {nelem,    3,NP,NP,nlevs},  m_dyn_grid->name());
  create_internal_field("FT",       {EL,       GP,GP,LEV}, {nelem,      NP,NP,nlevs},  m_dyn_grid->name());
  create_internal_field("v",        {EL,TL,CMP,GP,GP,LEV}, {nelem,NTL,2,NP,NP,nlevs},  m_dyn_grid->name());
  create_internal_field("vtheta_dp",{EL,TL,    GP,GP,LEV}, {nelem,NTL,  NP,NP,nlevs},  m_dyn_grid->name());
  create_internal_field("dp3d",     {EL,TL,    GP,GP,LEV}, {nelem,NTL,  NP,NP,nlevs},  m_dyn_grid->name());
  create_internal_field("w_i",      {EL,TL,    GP,GP,ILEV},{nelem,NTL,  NP,NP,nlevs+1},m_dyn_grid->name());
  create_internal_field("phinh_i",  {EL,TL,    GP,GP,ILEV},{nelem,NTL,  NP,NP,nlevs+1},m_dyn_grid->name());
  create_internal_field("ps",       {EL,TL,    GP,GP},     {nelem,NTL,  NP,NP},        m_dyn_grid->name());

  // For momentum and temperature, we back out tendencies on the ref grid,
  // and then remap them. So we need extra fields for FM and FT on the ref grid,
  // but we don't expose them, since they are internal. And we don't care about units.
  // We need a temp also for tracers, but at the time this fcn is called, we still don't
  // know the size of the tracers group. We can create that temp later.
  create_internal_field("FM",{COL,CMP,LEV},{ncols,3,NVL},m_ref_grid->name());
  create_internal_field("FT",{COL,    LEV},{ncols,  NVL},m_ref_grid->name());

  // Dynamics backs out tendencies from the states, and passes those to Homme.
  // After Homme completes, we remap the updates state to the ref grid.
  // Thus, is more convenient to use two different remappers: the pd remapper
  // will remap into Homme's forcing views, while the dp remapper will remap
  // from Homme's states.
  m_p2d_remapper = grids_manager->create_remapper_from_ref_grid(m_dyn_grid);
  m_d2p_remapper = grids_manager->create_remapper_to_ref_grid(m_dyn_grid);

  // Create separate remapper for Initial Conditions, since for those we
  // remap not into forcings, but directly into states.
  m_ic_remapper_fwd = grids_manager->create_remapper_from_ref_grid(m_dyn_grid);
  m_ic_remapper_bwd = grids_manager->create_remapper_from_ref_grid(m_dyn_grid);
}

void HommeDynamics::
set_computed_group_impl (const FieldGroup<Real>& group)
{
  const auto& name = group.m_info->m_group_name;
  EKAT_REQUIRE_MSG(name=="tracers",
    "Error! We were not expecting a field group called '" << name << "\n");

  EKAT_REQUIRE_MSG(not group.m_info->empty(),
    "Error! There should be at least one tracer (qv) in the '" + name + "' group.\n");

  EKAT_REQUIRE_MSG(group.m_info->m_bundled,
      "Error! Homme expects a bundled field for the group '" + name + "'.\n");

  // Now that we have Q, we have the exact count for tracers,
  // and we can use that info to setup tracers stuff in Homme
  const auto& c = Homme::Context::singleton();
  auto& params = c.get<Homme::SimulationParams>();
  if (params.qsize==0) {
    auto& tracers = c.get<Homme::Tracers>();
    const int qsize = group.m_info->size();
    params.qsize = qsize; // Set in the CXX data structure
    set_homme_param("qsize",qsize); // Set in the F90 module
    tracers.init(tracers.num_elems(),qsize);
  } else {
    EKAT_REQUIRE_MSG (params.qsize==group.m_info->size(),
        "Error! Input tracers group has a size different from the already initialized qsize param in Homme.\n"
        "    - homme qsize: " + std::to_string(params.qsize) + "\n"
        "    - tracers group size: " + std::to_string(group.m_info->size()) + "\n");
  }
}

int HommeDynamics::requested_buffer_size_in_bytes() const
{
  using namespace Homme;
  auto& c = Context::singleton();

  auto& elems   = c.get<Elements>();
  auto& tracers = c.get<Tracers>();
  auto& params  = c.get<SimulationParams>();

  auto& caar = c.create_if_not_there<CaarFunctor>(elems.num_elems(),params);
  auto& esf  = c.create_if_not_there<EulerStepFunctor>(elems.num_elems());
  auto& hvf  = c.create_if_not_there<HyperviscosityFunctor>(elems.num_elems(), params);
  auto& ff   = c.create_if_not_there<ForcingFunctor>(elems.num_elems(), tracers.num_elems(), tracers.num_tracers());
  auto& diag = c.create_if_not_there<Diagnostics> (elems.num_elems(),params.theta_hydrostatic_mode);
  auto& vrm  = c.create_if_not_there<VerticalRemapManager>(elems.num_elems());

  const bool need_dirk = (params.time_step_type==TimeStepType::IMEX_KG243 ||
                          params.time_step_type==TimeStepType::IMEX_KG254 ||
                          params.time_step_type==TimeStepType::IMEX_KG255 ||
                          params.time_step_type==TimeStepType::IMEX_KG355);
  if (need_dirk) {
    // Create dirk functor only if needed
    c.create_if_not_there<DirkFunctor>(elems.num_elems());
  }

  // Request buffer sizes in FunctorsBuffersManager and then
  // return the total bytes using the calculated buffer size.
  auto& fbm  = c.create_if_not_there<FunctorsBuffersManager>();
  fbm.request_size(caar.requested_buffer_size());
  fbm.request_size(esf.requested_buffer_size());
  fbm.request_size(hvf.requested_buffer_size());
  fbm.request_size(diag.requested_buffer_size());
  fbm.request_size(ff.requested_buffer_size());
  fbm.request_size(vrm.requested_buffer_size());
  if (need_dirk) {
    const auto& dirk = c.get<DirkFunctor>();
    fbm.request_size(dirk.requested_buffer_size());
  }

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
  fbm.allocate(buffer_manager.get_memory(), buffer_manager.allocated_bytes()/sizeof(Real));
}

void HommeDynamics::initialize_impl ()
{
  const auto& dgn = m_dyn_grid->name();
  const auto& rgn = m_ref_grid->name();

  // Use common names for some fields.
  alias_group_in  ("tracers",rgn,"Q");
  alias_group_out ("tracers",rgn,"Q");
  alias_group_in  ("tracers",dgn,"Q");
  alias_group_out ("tracers",dgn,"Q");

  // ------ Sanity checks ------- //

  // Nobody should claim to be a provider for dp, w_i, or tracers (on dyn grid).
  // Homme needs Q at the end of the timestep to not be altered by any other atm proc,
  // since it HAS to be saved correctly to allow exact restarts.
  const auto& rho_track = get_field_out("pseudo_density").get_header().get_tracking();
  const auto& w_i_track = get_field_out("w_int").get_header().get_tracking();
  const auto& Q_dyn_track = get_group_out("Q",dgn).m_bundle->get_header().get_tracking();
  EKAT_REQUIRE_MSG (
      rho_track.get_providers().size()==1,
      "Error! Someone other than Dynamics is trying to update the pseudo_density.\n");
  EKAT_REQUIRE_MSG (
      w_i_track.get_providers().size()==1,
      "Error! Someone other than Dynamics is trying to update the vertical velocity.\n");
  EKAT_REQUIRE_MSG (Q_dyn_track.get_providers().size()==1,
      "Error! Only HommeDynamics is allowed to update the tracers field on the dynamics grid.\n"
      "       Homme needs this field to be saved exactly as it is at the end of the dyn call,\n"
      "       in order to allow an exact restart.\n");

  // Create remaining internal fields
  // TODO: tracers_prev on physics grid is only needed if you need to remap dQ.
  const auto& Q_ref_lt = get_group_out("Q",rgn).m_bundle->get_header().get_identifier().get_layout();
  const auto& Q_dyn_lt = get_group_out("Q",dgn).m_bundle->get_header().get_identifier().get_layout();
  create_internal_field("FQ",Q_dyn_lt.tags(),Q_dyn_lt.dims(),m_dyn_grid->name());
  create_internal_field("Q_prev",Q_ref_lt.tags(),Q_ref_lt.dims(),m_ref_grid->name());

  // Sets the scream views into the hommexx internal data structures
  init_homme_views ();

  // Import I.C. from the ref grid to the dyn grid.
  import_initial_conditions ();

  // Update p_int and p_mid. Other models might import these values from
  // SurfaceCoupling during initialization. They will be overwritten
  // during postprocessing.
  update_pressure();

  // Complete homme model initialization
  prim_init_model_f90 ();
}

void HommeDynamics::run_impl (const int dt)
{
  try {
    // Prepare inputs for homme
    Kokkos::fence();
    homme_pre_process (dt);

    // Note: Homme's step lasts homme_dt*max(dt_remap_factor,dt_tracers_factor), and it must divide dt.
    // We neeed to compute dt/homme_dt, and subcycle homme that many times
    const int nsplit = get_homme_nsplit_f90(dt);
    for (int subiter=0; subiter<nsplit; ++subiter) {
      Kokkos::fence();
      prim_run_f90 ();
    }

    // Post process Homme's output, to produce what the rest of Atm expects
    Kokkos::fence();
    homme_post_process ();
  } catch (std::exception& e) {
    EKAT_ERROR_MSG(e.what());
  } catch (...) {
    EKAT_ERROR_MSG("Something went wrong, but we don't know what.\n");
  }
}

void HommeDynamics::finalize_impl (/* what inputs? */)
{
  Homme::Context::singleton().finalize_singleton();
  prim_finalize_f90();
}

void HommeDynamics::homme_pre_process (const int dt) {
  // T and uv tendencies are backed out on the ref grid.
  // Homme takes care of turning the FT tendency into a tendency for VTheta_dp.
  using KT = KokkosTypes<DefaultDevice>;

  constexpr int N = sizeof(Homme::Scalar) / sizeof(Real);
  using Pack = ekat::Pack<Real,N>;

  const int ncols = m_ref_grid->get_num_local_dofs();
  const int nlevs = m_ref_grid->get_num_vertical_levels();
  const int npacks = ekat::PackInfo<N>::num_packs(nlevs);

  auto T  = get_field_out("T_mid").get_view<Pack**>();
  auto FT = get_field_out("T_mid_prev").get_view<Pack**>();
  auto v  = get_field_out("horiz_winds").get_view<Pack***>();
  auto w  = get_field_out("w_int").get_view<Pack**>();
  auto FM = get_field_out("horiz_winds_prev").get_view<Pack***>();

  // If there are other atm procs updating the vertical velocity,
  // then we need to compute forcing for w as well
  const bool has_w_forcing = get_field_out("w_int").get_header().get_tracking().get_providers().size()>1;
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

    if (has_w_forcing) {
      // Vertical velocity forcing
      // Recall: fm(2) stores forcing for w_i at [0,num_int_levels-1],
      //         since w_i at surf is determined with no penetration bc.
      const auto& w_new =  w(icol,ilev);
      const auto& w_old = FM(icol,2,ilev);
            auto& fw    = FM(icol,2,ilev);
      fw = (w_new - w_old) / dt;
    }
  });

  // Remap FT, FM, and Q
  m_p2d_remapper->remap(true);

  using namespace Homme;
  const auto& c = Context::singleton();
  const auto& params = c.get<SimulationParams>();
  auto& tl = c.get<TimeLevel>();
  auto& ff = c.get<ForcingFunctor>();

  const auto ftype = params.ftype;
  const auto& tracers = c.get<Tracers>();
  const auto& state = c.get<ElementsState>();
  auto Q    = tracers.Q;
  auto FQ   = tracers.fq;
  auto Qdp  = tracers.qdp;
  auto dp3d = state.m_dp3d;

  // Note: np1_qdp and n0_qdp are 'deduced' from tl.nstep, so the
  //       following call may not even change them (i.e., they are
  //       not updated regardless). So if they were already up-to-date,
  //       the following call will do nothing.
  tl.update_tracers_levels(params.qsplit);

  const auto n0 = tl.n0;  // The time level where pd coupling remapped into
  constexpr int NVL = HOMMEXX_NUM_LEV;
  const int qsize = params.qsize;
  const auto n0_qdp = tl.n0_qdp;

  // At this point, FQ contains Qnew (coming from physics).
  // Depending on ftype, we are going to do different things:
  //  ftype=0: FQ = dp*(Qnew-Qold) / dt
  //  ftype=2: qdp = dp*Qnew
  switch(ftype) {
    case ForcingAlg::FORCING_DEBUG:
      // Back out tracers tendency for Qdp
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
      break;
    case ForcingAlg::FORCING_2:
      // Hard adjustment of qdp
      ff.tracers_forcing(dt,n0,n0_qdp,true,params.moisture);
      break;
    default:
      EKAT_ERROR_MSG ("Error! Unexpected/unsupported forcing algorithm.\n"
                      "  ftype: " + std::to_string(Homme::etoi(ftype)) + "\n");
  }
}

void HommeDynamics::homme_post_process () {
  const auto& rgn = m_ref_grid->name();

  // Remap outputs to ref grid
  m_d2p_remapper->remap(true);

  using KT = KokkosTypes<DefaultDevice>;
  constexpr int N = sizeof(Homme::Scalar) / sizeof(Real);
  using Pack = ekat::Pack<Real,N>;
  using ColOps = ColumnOps<DefaultDevice,Real>;
  using PF = PhysicsFunctions<DefaultDevice>;

  // Convert VTheta_dp->T, store T,uv, and possibly w in FT, FM,
  // compute p_int on ref grid.
  const auto ps_view = get_field_out("ps").get_view<Real*>();
  const auto dp_view = get_field_out("pseudo_density").get_view<Pack**>();
  const auto p_mid_view = get_field_out("p_mid").get_view<Pack**>();
  const auto p_int_view = get_field_out("p_int").get_view<Pack**>();
  const auto Q_view  = get_group_out("Q",rgn).m_bundle->get_view<Pack***>();

  const auto T_view  = get_field_out("T_mid").get_view<Pack**>();
  const auto v_view  = get_field_out("horiz_winds").get_view<Pack***>();
  const auto w_view  = get_field_out("w_int").get_view<Pack**>();
  const auto T_prev_view = get_field_out("T_mid_prev").get_view<Pack**>();
  const auto v_prev_view = get_field_out("horiz_winds_prev").get_view<Pack***>();
  const auto w_prev_view = get_field_out("w_int_prev").get_view<Pack**>();

  const auto ncols = m_ref_grid->get_num_local_dofs();
  const auto nlevs = m_ref_grid->get_num_vertical_levels();
  const auto npacks= ekat::PackInfo<N>::num_packs(nlevs);

  using ESU = ekat::ExeSpaceUtils<KT::ExeSpace>;
  const auto policy = ESU::get_thread_range_parallel_scan_team_policy(ncols,npacks);

  // If there are other atm procs updating the vertical velocity,
  // then we need to store w_int_old (to compute forcing for w next iteration)
  const bool has_w_forcing = get_field_out("w_int").get_header().get_tracking().get_providers().size()>1;

  // Establish the boundary condition for the TOA
  const auto& c = Homme::Context::singleton();
  const auto& hvcoord = c.get<Homme::HybridVCoord>();
  const auto ps0 = hvcoord.ps0 * hvcoord.hybrid_ai0;

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const KT::MemberType& team) {
    const int& icol = team.league_rank();

    // Compute p_int and p_mid
    auto dp = ekat::subview(dp_view,icol);
    auto p_mid = ekat::subview(p_mid_view,icol);
    auto p_int = ekat::subview(p_int_view,icol);

    ColOps::column_scan<true>(team,nlevs,dp,p_int,ps0);
    team.team_barrier();
    ColOps::compute_midpoint_values(team,nlevs,p_int,p_mid);
    team.team_barrier();

    // Convert VTheta_dp->VTheta->Theta->T
    auto T   = ekat::subview(T_view,icol);
    auto v   = ekat::subview(v_view,icol);
    auto w   = ekat::subview(w_view,icol);
    auto qv  = ekat::subview(Q_view,icol,0);

    auto T_prev = ekat::subview(T_prev_view,icol);
    auto v_prev = ekat::subview(v_prev_view,icol);
    auto w_prev = ekat::subview(w_prev_view,icol);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,npacks),
                         [&](const int ilev) {
      // VTheta_dp->VTheta->Theta->T
      auto& T_val = T(ilev);
      T_val /= dp(ilev);
      T_val = PF::calculate_temperature_from_virtual_temperature(T_val,qv(ilev));
      T_val = PF::calculate_T_from_theta(T_val,p_mid(ilev));

      // Store T, v (and possibly w) at end of the dyn timestep (to back out tendencies later)
      T_prev(ilev) = T_val;
      v_prev(0,ilev) = v(0,ilev);
      v_prev(1,ilev) = v(1,ilev);
      if (has_w_forcing) {
        w_prev(ilev) = w(ilev);
      }
    });
  });
}

void HommeDynamics::
create_internal_field (const std::string& name,
                       const std::vector<FieldTag>& tags,
                       const std::vector<int>& dims,
                       const std::string& grid) {
  using namespace ekat::units;
  FieldIdentifier id(name,FieldLayout{tags,dims},Units::nondimensional(),grid);

  const auto lt = get_layout_type(id.get_layout().tags());

  // Only request packed field for 3d quantities
  int pack_size = 1;
  if (lt==LayoutType::Scalar3D || lt==LayoutType::Vector3D || lt==LayoutType::Tensor3D) {
    pack_size = sizeof(Homme::Scalar) / sizeof(Real);
  }

  field_type f(id);
  f.get_header().get_alloc_properties().request_allocation<Real>(pack_size);
  f.allocate_view();

  m_internal_fields[name+grid] = f;
}

const Field<Real>& HommeDynamics::
get_internal_field (const std::string& name, const std::string& grid) const {
  auto it = m_internal_fields.find(name+grid);
  EKAT_REQUIRE_MSG (it!=m_internal_fields.end(),
      "Error! Internal field '" + name + "' on grid '" + grid + "' not found.\n");
  return it->second;
}

void HommeDynamics::init_homme_views () {
  const auto& dgn = m_dyn_grid->name();

  const auto& c = Homme::Context::singleton();
  auto& state = c.get<Homme::ElementsState>();
  auto& tracers = c.get<Homme::Tracers>();
  auto& forcing = c.get<Homme::ElementsForcing>();

  constexpr int NGP  = HOMMEXX_NP;
  constexpr int NTL  = HOMMEXX_NUM_TIME_LEVELS;
  constexpr int NVL  = HOMMEXX_NUM_LEV;
  constexpr int NVLI = HOMMEXX_NUM_LEV_P;

  const int nelem = m_dyn_grid->get_num_local_dofs()/(NGP*NGP);
  const int qsize = get_group_out("Q",dgn).m_info->size();

  // Print homme's parameters, so user can see whether something wasn't set right.
  if (get_comm().am_i_root()) c.get<Homme::SimulationParams>().print();

  // ------------ Set views in Homme ------------- //
  // Velocity
  auto v_in = get_internal_field("v",dgn).get_view<Homme::Scalar*[NTL][2][NP][NP][NVL]>();
  using v_type = std::remove_reference<decltype(state.m_v)>::type;
  state.m_v = v_type (v_in.data(),nelem);

  // Virtual potential temperature
  auto vtheta_in = get_internal_field("vtheta_dp",dgn).get_view<Homme::Scalar*[NTL][NP][NP][NVL]>();
  using vtheta_type = std::remove_reference<decltype(state.m_vtheta_dp)>::type;
  state.m_vtheta_dp = vtheta_type(vtheta_in.data(),nelem);

  // Geopotential
  auto phi_in = get_internal_field("phinh_i",dgn).get_view<Homme::Scalar*[NTL][NP][NP][NVLI]>();
  using phi_type = std::remove_reference<decltype(state.m_phinh_i)>::type;
  state.m_phinh_i = phi_type(phi_in.data(),nelem);

  // Vertical velocity
  auto w_in = get_internal_field("w_i",dgn).get_view<Homme::Scalar*[NTL][NP][NP][NVLI]>();
  using w_type = std::remove_reference<decltype(state.m_w_i)>::type;
  state.m_w_i = w_type(w_in.data(),nelem);

  // Pseudo-density
  auto dp3d_in = get_internal_field("dp3d",dgn).template get_view<Homme::Scalar*[NTL][NP][NP][NVL]>();
  using dp3d_type = std::remove_reference<decltype(state.m_dp3d)>::type;
  state.m_dp3d = dp3d_type(dp3d_in.data(),nelem);

  // Surface pressure
  auto ps_in = get_internal_field("ps",dgn).template get_view<Real*[NTL][NP][NP]>();
  using ps_type = std::remove_reference<decltype(state.m_ps_v)>::type;
  state.m_ps_v = ps_type(ps_in.data(),nelem);

  // Tracers
  auto q_in = get_group_out("Q",dgn).m_bundle->template get_view<Homme::Scalar**[NP][NP][NVL]>();
  using q_type = std::remove_reference<decltype(tracers.Q)>::type;
  tracers.Q = q_type(q_in.data(),nelem,qsize);

  // Tracers forcing
  auto fq_in = get_internal_field("FQ",dgn).template get_view<Homme::Scalar**[NP][NP][NVL]>();
  using fq_type = std::remove_reference<decltype(tracers.fq)>::type;
  tracers.fq = fq_type(fq_in.data(),nelem,qsize);

  // Temperature forcing
  auto ft_in = get_internal_field("FT",dgn).template get_view<Homme::Scalar*[NP][NP][NVL]>();
  using ft_type = std::remove_reference<decltype(forcing.m_ft)>::type;
  forcing.m_ft = ft_type(ft_in.data(),nelem);

  // Momentum forcing
  auto fm_in = get_internal_field("FM",dgn).template get_view<Homme::Scalar**[NP][NP][NVL]>();
  using fm_type = std::remove_reference<decltype(forcing.m_fm)>::type;
  forcing.m_fm = fm_type(fm_in.data(),nelem);
}

void HommeDynamics::import_initial_conditions () {
  const auto& dgn = m_dyn_grid->name();
  const auto& rgn = m_ref_grid->name();

  // Setup the p2d and d2p remappers
  m_p2d_remapper->registration_begins();
  m_d2p_remapper->registration_begins();

  m_p2d_remapper->register_field(get_field_out("T_mid_prev"),get_internal_field("FT",dgn));
  m_p2d_remapper->register_field(get_field_out("horiz_winds_prev"),get_internal_field("FM",dgn));
  m_p2d_remapper->register_field(*get_group_out("Q",rgn).m_bundle, get_internal_field("FQ",dgn));

  m_d2p_remapper->register_field(get_internal_field("vtheta_dp",dgn),get_field_out("T_mid"));
  m_d2p_remapper->register_field(get_internal_field("v",dgn),get_field_out("horiz_winds"));
  m_d2p_remapper->register_field(get_internal_field("dp3d",dgn), get_field_out("pseudo_density"));
  m_d2p_remapper->register_field(get_internal_field("phinh_i",dgn), get_field_out("phi_int"));
  m_d2p_remapper->register_field(get_internal_field("w_i",dgn), get_field_out("w_int"));
  m_d2p_remapper->register_field(get_internal_field("ps",dgn), get_field_out("ps"));
  m_d2p_remapper->register_field(*get_group_out("Q",dgn).m_bundle,*get_group_out("Q",rgn).m_bundle);

  m_p2d_remapper->registration_ends();
  m_d2p_remapper->registration_ends();

  // Setup, run, and destroy the IC remapper, to remap IC directly into Homme's states
  m_ic_remapper_fwd->registration_begins();
  m_ic_remapper_fwd->register_field(get_field_out("w_int",rgn),get_internal_field("w_i",dgn));
  m_ic_remapper_fwd->register_field(get_field_out("phi_int",rgn),get_internal_field("phinh_i",dgn));
  m_ic_remapper_fwd->register_field(get_field_out("horiz_winds",rgn),get_internal_field("v",dgn));
  m_ic_remapper_fwd->register_field(get_field_out("pseudo_density",rgn),get_internal_field("dp3d",dgn));
  m_ic_remapper_fwd->register_field(get_field_out("ps",rgn),get_internal_field("ps",dgn));
  m_ic_remapper_fwd->register_field(get_field_out("T_mid",rgn),get_internal_field("vtheta_dp",dgn));
  // m_ic_remapper->register_field(get_field_out("Q",rgn),get_internal_field("Q"));
  m_ic_remapper_fwd->registration_ends();
  
  m_ic_remapper_fwd->remap(true);
  m_ic_remapper_fwd = nullptr;

  // We keep a copy of Q_prev on the ref grid, since we remap dQ=Qnew-Qprev from
  // ref grid to dyn grid, and dynamics will then do Qnew_dyn = Qold_dyn+dQ_dyn
  // It would be tempting to init Q_prev on ref grid with the value of Q on the
  // ref grid, but that's not necessarily the same (it depends on the atm procs
  // order in the atm DAG). So we just create a remapper dyn->ref on the fly
  m_ic_remapper_bwd->registration_begins();
  m_ic_remapper_bwd->register_field(get_internal_field("Q_prev",rgn),*get_group_out("Q",dgn).m_bundle);
  m_ic_remapper_bwd->registration_ends();
  m_ic_remapper_bwd->remap(false);
  m_ic_remapper_bwd = nullptr;

  // Convert T->vtheta_dp (in place).
  constexpr int N = sizeof(Homme::Scalar) / sizeof(Real);
  using KT = KokkosTypes<DefaultDevice>;
  using Pack = ekat::Pack<Real,N>;
  using ColOps = ColumnOps<DefaultDevice,Real>;
  using PF = PhysicsFunctions<DefaultDevice>;

  const auto dp3d      = get_internal_field("dp3d",dgn).get_view<Pack*****>();
  const auto ps        = get_internal_field("ps",dgn).get_view<Real****>();
  const auto v         = get_internal_field("v",dgn).get_view<Pack******>();
  const auto w_i       = get_internal_field("w_i",dgn).get_view<Pack*****>();
  const auto phinh_i   = get_internal_field("phinh_i",dgn).get_view<Pack*****>();
  const auto vtheta_dp = get_internal_field("vtheta_dp",dgn).get_view<Pack*****>();
  const auto Q         = get_group_out("tracers",dgn).m_bundle->get_view<Pack*****>();

  const auto& c = Homme::Context::singleton();

  // Time slices indices
  const int n0 = c.get<Homme::TimeLevel>().n0;
  const int nm1 = c.get<Homme::TimeLevel>().nm1;
  const int np1 = c.get<Homme::TimeLevel>().np1;

  // Spatial extents of views
  constexpr int NGP  = HOMMEXX_NP;
  constexpr int NVL  = HOMMEXX_NUM_LEV;
  constexpr int NVLI = HOMMEXX_NUM_LEV_P;
  const int nelem = m_dyn_grid->get_num_local_dofs()/(NGP*NGP);
  const int nlevs = m_dyn_grid->get_num_vertical_levels();

  const auto& hvcoord = c.get<Homme::HybridVCoord>();
  const auto ps0 = hvcoord.ps0 * hvcoord.hybrid_ai0;

  using ESU = ekat::ExeSpaceUtils<KT::ExeSpace>;
  const auto policy = ESU::get_thread_range_parallel_scan_team_policy(nelem*NP*NP,NVL);

  // Need two temporaries, for pi_mid and pi_int
  ekat::WorkspaceManager<Pack,DefaultDevice> wsm(NVLI,2,policy);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const KT::MemberType& team) {
    const int ie  =  team.league_rank() / (NP*NP);
    const int igp = (team.league_rank() / NP) % NP;
    const int jgp =  team.league_rank() % NP;

    auto ws = wsm.get_workspace(team);

    // Compute p_mid
    auto dp = ekat::subview(dp3d,ie,n0,igp,jgp);

    auto p_int = ws.take("p_int");
    auto p_mid = ws.take("p_mid");
    ColOps::column_scan<true>(team,nlevs,dp,p_int,ps0);
    team.team_barrier();
    ColOps::compute_midpoint_values(team,nlevs,p_int,p_mid);
    team.team_barrier();
    
    // Convert T->Theta->Theta*dp->VTheta*dp
    auto T      = ekat::subview(vtheta_dp,ie,n0,igp,jgp);
    auto vTh_dp = ekat::subview(vtheta_dp,ie,n0,igp,jgp);
    auto qv     = ekat::subview(Q,ie,0,igp,jgp);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,NVL),
                         [&](const int ilev) {
      const auto vthdp = dp(ilev)*PF::calculate_theta_from_T(T(ilev),p_mid(ilev));
      vTh_dp(ilev) = PF::calculate_virtual_temperature(vthdp,qv(ilev));

      // Copy states from the n0 timelevel to all the other ones
      dp3d(ie,nm1,igp,jgp,ilev) = dp3d(ie,n0,igp,jgp,ilev);
      dp3d(ie,np1,igp,jgp,ilev) = dp3d(ie,n0,igp,jgp,ilev);
      
      v(ie,nm1,0,igp,jgp,ilev) = v(ie,n0,0,igp,jgp,ilev);
      v(ie,nm1,1,igp,jgp,ilev) = v(ie,n0,1,igp,jgp,ilev);
      v(ie,np1,0,igp,jgp,ilev) = v(ie,n0,0,igp,jgp,ilev);
      v(ie,np1,1,igp,jgp,ilev) = v(ie,n0,1,igp,jgp,ilev);

      vtheta_dp(ie,nm1,igp,jgp,ilev) = vtheta_dp(ie,n0,igp,jgp,ilev);
      vtheta_dp(ie,np1,igp,jgp,ilev) = vtheta_dp(ie,n0,igp,jgp,ilev);

      w_i(ie,nm1,igp,jgp,ilev) = w_i(ie,n0,igp,jgp,ilev);
      w_i(ie,np1,igp,jgp,ilev) = w_i(ie,n0,igp,jgp,ilev);

      phinh_i(ie,nm1,igp,jgp,ilev) = phinh_i(ie,n0,igp,jgp,ilev);
      phinh_i(ie,np1,igp,jgp,ilev) = phinh_i(ie,n0,igp,jgp,ilev);
    });

    // Copy ps (and last interface of w_i and phinh_i)  from the n0 timelevel to all the other ones
    Kokkos::single(Kokkos::PerTeam(team),[&](){
      ps(ie,nm1,igp,jgp) = ps(ie,n0,igp,jgp);
      ps(ie,np1,igp,jgp) = ps(ie,n0,igp,jgp);
      w_i(ie,nm1,igp,jgp,NVLI-1) = w_i(ie,n0,igp,jgp,NVLI-1);
      w_i(ie,np1,igp,jgp,NVLI-1) = w_i(ie,n0,igp,jgp,NVLI-1);

      phinh_i(ie,nm1,igp,jgp,NVLI-1) = phinh_i(ie,n0,igp,jgp,NVLI-1);
      phinh_i(ie,np1,igp,jgp,NVLI-1) = phinh_i(ie,n0,igp,jgp,NVLI-1);
    });

    // Release the scratch mem
    ws.release(p_int);
    ws.release(p_mid);
  });
  Kokkos::fence();

  // Init qdp = q*dp (q comes from initial conditions)
  auto& tracers = c.get<Homme::Tracers>();
  const auto qdp = tracers.qdp;
  const auto q   = tracers.Q;
  const auto dp  = c.get<Homme::ElementsState>().m_dp3d;
  const int n0_qdp = c.get<Homme::TimeLevel>().n0_qdp;
  const int qsize = get_group_out("tracers",dgn).m_info->size();

  Kokkos::parallel_for(Kokkos::RangePolicy<>(0,nelem*qsize*NP*NP*NVL),
                       KOKKOS_LAMBDA (const int idx) {
    const int ie =  idx / (qsize*NP*NP*NVL);
    const int iq = (idx / (NP*NP*NVL)) % qsize;
    const int ip = (idx / (NP*NVL)) % NP;
    const int jp = (idx / NVL) % NP;
    const int k  = idx % NVL;

    qdp(ie,n0_qdp,iq,ip,jp,k) = q(ie,iq,ip,jp,k) * dp(ie,n0,ip,jp,k);
  });
}
// =========================================================================================
void HommeDynamics::update_pressure() {
  using KT = KokkosTypes<DefaultDevice>;
  constexpr int N = sizeof(Homme::Scalar) / sizeof(Real);
  using Pack = ekat::Pack<Real,N>;
  using ColOps = ColumnOps<DefaultDevice,Real>;

  const auto ncols = m_ref_grid->get_num_local_dofs();
  const auto nlevs = m_ref_grid->get_num_vertical_levels();
  const auto npacks= ekat::PackInfo<N>::num_packs(nlevs);

  const auto& c = Homme::Context::singleton();
  const auto& hvcoord = c.get<Homme::HybridVCoord>();
  const auto ps0 = hvcoord.ps0 * hvcoord.hybrid_ai0;

  const auto dp_view    = get_field_out("pseudo_density").get_view<Pack**>();
  const auto p_int_view = get_field_out("p_int").get_view<Pack**>();
  const auto p_mid_view = get_field_out("p_mid").get_view<Pack**>();

  using ESU = ekat::ExeSpaceUtils<KT::ExeSpace>;
  const auto policy = ESU::get_thread_range_parallel_scan_team_policy(ncols,npacks);

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const KT::MemberType& team) {
    const int& icol = team.league_rank();

    auto dp = ekat::subview(dp_view,icol);
    auto p_mid = ekat::subview(p_mid_view,icol);
    auto p_int = ekat::subview(p_int_view,icol);

    ColOps::column_scan<true>(team,nlevs,dp,p_int,ps0);
    team.team_barrier();
    ColOps::compute_midpoint_values(team,nlevs,p_int,p_mid);
  });
}
// =========================================================================================
void HommeDynamics::
check_computed_fields_impl () {
    // Note: We are seeing near epsilon negative values in a handful of places,
    // The strategy is to 
    // 1. First check that no values are sufficiently negative as to require an error.
    //    i.e. below some tolerance.
    // 2. Clip all negative values to zero.
    // TODO: Construct a more robust check that compares the value of Q against an
    // average value or maximum value over each column.  That way we can use a relative
    // error as our threshold, rather than an arbitrary tolerance.
   
    // Grab the pointer to the tracer group. 
    const auto& rgn = m_ref_grid->name();
          auto& Q   = *get_group_out("Q",rgn).m_bundle;
    // Create a local copy of a lower bound check to ensure we are not encountering truly
    // bad negative tracer values.
    Real tol = -1e-20;
    auto lower_bound_check = std::make_shared<FieldLowerBoundCheck<Real>>(tol);
    lower_bound_check->check(Q);
    // Now repair negative tracers using a lower bounds check at 0.0
    auto lower_bound_repair = std::make_shared<FieldLowerBoundCheck<Real>>(0.0);
    lower_bound_repair->repair(Q);
  
}
// =========================================================================================

} // namespace scream
