#include "atmosphere_dynamics.hpp"

// HOMMEXX Includes
#include "Context.hpp"
#include "Elements.hpp"
#include "ElementsForcing.hpp"
#include "ForcingFunctor.hpp"
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

// Ekat includes
#include "ekat/ekat_assert.hpp"
#include "ekat/kokkos//ekat_subview_utils.hpp"

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
  //       T and VTheta on dyn grid, vtheta_dp, we remap T_mid from phys directly
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

  // Dynamics backs out tendencies from the states, and pass those to Homme.
  // After Homme completes, we remap the updates state to the ref grid.
  // Thus, is more convenient to use two different remappers: the pd remapper
  // will remap into Homme's forcing views, while the dp remapper will remap
  // from Homme's states.

  // Helper lambda to make the following lines shorter.
  // NOTE: since these fields are internal, and they don't end up in the FieldManager,
  //       we don't care about units.
  // Note: the dyn grid names do not follow the FieldManager conventions,
  //       but rather the var names in Homme, so they can be easily found there. 
  create_dyn_field("FM",{EL,CMP,GP,GP,LEV},{ne,3,NP,NP,NVL});
  create_dyn_field("FT",{EL,    GP,GP,LEV},{ne,  NP,NP,NVL});

  create_dyn_field("v",        {EL,TL,CMP,GP,GP,LEV}, {ne,NTL,2,NP,NP,NVL});
  create_dyn_field("vtheta_dp",{EL,TL,    GP,GP,LEV}, {ne,NTL,  NP,NP,NVL});
  create_dyn_field("dp3d",     {EL,TL,    GP,GP,LEV}, {ne,NTL,  NP,NP,NVL});
  create_dyn_field("w_i",      {EL,TL,    GP,GP,ILEV},{ne,NTL,  NP,NP,NVL+1});
  create_dyn_field("phinh_i",  {EL,TL,    GP,GP,ILEV},{ne,NTL,  NP,NP,NVL+1});
  create_dyn_field("ps",       {EL,TL,    GP,GP},     {ne,NTL,  NP,NP});

  // For momentum and temperature, we back out tendencies on the ref grid,
  // and then remap them. So we need extra fields for FM and FT on the ref grid,
  // but we don't expose them, since they are internal. And we don't care about units.
  field_type FM(FID("FM",FL({COL,CMP,LEV},{ncols,3,NVL}),nondim,rgn));
  field_type FT(FID("FT",FL({COL,LEV},{ncols,NVL}),nondim,rgn));
  FM.get_header().get_alloc_properties().request_allocation<Real>(N);
  FT.get_header().get_alloc_properties().request_allocation<Real>(N);
  FM.allocate_view();
  FT.allocate_view();
  m_ref_grid_fields["FM"] = FM;
  m_ref_grid_fields["FT"] = FT;
  
  // Input-output groups
  // Note: tracers_prev is tracers after the last dyn step, and is used to back out tracer tendencies.
  //       We want it directy on the dyn grid, and must contain all fields in the tracers group.
  //       So we create it as an 'alias' of "tracers", but on a different grid.
  GroupRequest tracers("tracers",m_ref_grid->name(),N, Bundling::Required);
  add_group<Updated>(tracers);

  // Q_prev is only needed for ftype=0.
  const int ftype = get_homme_param<int>("ftype");
  EKAT_REQUIRE_MSG(ftype==0 || ftype==2,
                     "Error! The scream interface to homme *assumes* ftype to be 0 or 2.\n"
                     "       Found " + std::to_string(ftype) + " instead.\n");
  if (ftype==0) {
    GroupRequest tracers_prev("tracers_prev",m_dyn_grid->name(),N, Bundling::Required,
                              &tracers,Relationship::Alias);
    add_group<Updated>(tracers_prev);
  }

  // Create the PD and DP remappers
  m_p2d_remapper = grids_manager->create_remapper_from_ref_grid(m_dyn_grid);
  m_d2p_remapper = grids_manager->create_remapper_to_ref_grid(m_dyn_grid);

  // Create separate remapper for Initial Conditions, since for those we
  // remap not into forcings, but directly into states.
  m_ic_remapper = grids_manager->create_remapper_from_ref_grid(m_dyn_grid);
}

void HommeDynamics::
set_updated_group (const FieldGroup<Real>& group)
{
  const auto& name = group.m_info->m_group_name;
  const int ftype = get_homme_param<int>("ftype");
  EKAT_REQUIRE_MSG(name=="tracers" || (ftype==0 && name=="tracers_prev"),
    "Error! We were not expecting a field group called '" << name << "\n");

  EKAT_REQUIRE_MSG(not group.m_info->empty(),
    "Error! There should be at least one tracer (qv) in the '" + name + "' group.\n");

  EKAT_REQUIRE_MSG(group.m_info->m_bundled,
      "Error! Homme expects a bundled field for the group '" + name + "'.\n");

  if (m_dyn_grid_fields.find("Q")==m_dyn_grid_fields.end()) {
    // We can use the exact count for the tracers group, to set homme's qsize,
    // initialize the Tracers struct accordingly, and create Q/FQ dyn fields.
    const int qsize = group.m_info->size();
    auto& tracers = Homme::Context::singleton().get<Homme::Tracers>();

    // Now that we have Q, we have the exact count for tracers,
    // and we can use that info to setup tracers stuff in Homme
    auto& params = Homme::Context::singleton().get<Homme::SimulationParams>();
    params.qsize = qsize; // Set in the CXX data structure
    set_homme_param("qsize",qsize); // Set in the F90 module
    tracers.init(tracers.num_elems(),qsize);

    // Create Q and FQ on dyn grid
    using namespace ShortFieldTagsNames;

    const int ngp  = HOMMEXX_NP;
    const int ne = m_dyn_grid->get_num_local_dofs()/(ngp*ngp);
    const int nlevs = m_dyn_grid->get_num_vertical_levels();
    create_dyn_field("Q", {EL,CMP,GP,GP,LEV},{ne,qsize,ngp,ngp,nlevs});
    create_dyn_field("FQ",{EL,CMP,GP,GP,LEV},{ne,qsize,ngp,ngp,nlevs});
  }

  if (name=="tracers") {
    m_ref_grid_fields["Q"] = *group.m_bundle;
  } else {
    m_dyn_grid_fields["Q"] = *group.m_bundle;
  }
}

void HommeDynamics::initialize_impl (const util::TimeStamp& /* t0 */)
{
  const auto& c = Homme::Context::singleton();

  // Print homme's parameters, so user can see whether something wasn't set right.
  c.get<Homme::SimulationParams>().print();

  constexpr int NTL  = HOMMEXX_NUM_TIME_LEVELS;
  constexpr int NVL  = HOMMEXX_NUM_LEV;
  constexpr int NVLI = HOMMEXX_NUM_LEV_P;
  {
    // Now that all dyn fields are created (including tracers), replace Hommexx views
    auto& state = c.get<Homme::ElementsState>();
    auto& tracers = c.get<Homme::Tracers>();
    auto& forcing = c.get<Homme::ElementsForcing>();

    const int num_elems = state.num_elems();
    const int num_tracers = tracers.num_tracers();

    // Velocity
    auto& v = state.m_v;
    auto v_in = m_dyn_grid_fields.at("v").get_reshaped_view<Homme::Scalar*[NTL][2][NP][NP][NVL]>();
    using v_type = std::remove_reference<decltype(v)>::type;
    v = v_type (v_in.data(),num_elems);

    // Virtual potential temperature
    auto& vtheta = state.m_vtheta_dp;
    auto vtheta_in = m_dyn_grid_fields.at("vtheta_dp").get_reshaped_view<Homme::Scalar*[NTL][NP][NP][NVL]>();
    using vtheta_type = std::remove_reference<decltype(vtheta)>::type;
    vtheta = vtheta_type(vtheta_in.data(),num_elems);

    // Geopotential
    auto& phi = state.m_phinh_i;
    auto phi_in = m_dyn_grid_fields.at("phinh_i").get_reshaped_view<Homme::Scalar*[NTL][NP][NP][NVLI]>();
    using phi_type = std::remove_reference<decltype(phi)>::type;
    phi = phi_type(phi_in.data(),num_elems);

    // Vertical velocity
    auto& w = state.m_w_i;
    auto w_in = m_dyn_grid_fields.at("w_i").get_reshaped_view<Homme::Scalar*[NTL][NP][NP][NVLI]>();
    using w_type = std::remove_reference<decltype(w)>::type;
    w = w_type(w_in.data(),num_elems);

    // Pseudo-density
    auto& dp3d = state.m_dp3d;
    auto dp3d_in = m_dyn_grid_fields.at("dp3d").template get_reshaped_view<Homme::Scalar*[NTL][NP][NP][NVL]>();
    using dp3d_type = std::remove_reference<decltype(dp3d)>::type;
    dp3d = dp3d_type(dp3d_in.data(),num_elems);

    // Surface pressure
    auto& ps = state.m_ps_v;
    auto ps_in = m_dyn_grid_fields.at("ps").template get_reshaped_view<Real*[NTL][NP][NP]>();
    using ps_type = std::remove_reference<decltype(ps)>::type;
    ps = ps_type(ps_in.data(),num_elems);

    // Tracers
    auto& q = tracers.Q;
    auto q_in = m_dyn_grid_fields.at("Q").template get_reshaped_view<Homme::Scalar**[NP][NP][NVL]>();
    using q_type = std::remove_reference<decltype(q)>::type;
    q = q_type(q_in.data(),num_elems,num_tracers);

    // Tracers forcing
    auto& fq = tracers.fq;
    auto fq_in = m_dyn_grid_fields.at("FQ").template get_reshaped_view<Homme::Scalar**[NP][NP][NVL]>();
    using fq_type = std::remove_reference<decltype(fq)>::type;
    fq = fq_type(fq_in.data(),num_elems,num_tracers);

    // Temperature forcing
    auto& ft = forcing.m_ft;
    auto ft_in = m_dyn_grid_fields.at("FT").template get_reshaped_view<Homme::Scalar*[NP][NP][NVL]>();
    using ft_type = std::remove_reference<decltype(ft)>::type;
    ft = ft_type(ft_in.data(),num_elems);

    // Momentum forcing
    auto& fm = forcing.m_fm;
    auto fm_in = m_dyn_grid_fields.at("FM").template get_reshaped_view<Homme::Scalar**[NP][NP][NVL]>();
    using fm_type = std::remove_reference<decltype(fm)>::type;
    fm = fm_type(fm_in.data(),num_elems);
  }

  // Setup the p2d and d2p remappers
  m_p2d_remapper->registration_begins();
  m_d2p_remapper->registration_begins();

  m_p2d_remapper->register_field(m_ref_grid_fields.at("T_mid_prev"),m_dyn_grid_fields.at("FT"));
  m_p2d_remapper->register_field(m_ref_grid_fields.at("horiz_winds_prev"),m_dyn_grid_fields.at("FM"));
  m_p2d_remapper->register_field(m_ref_grid_fields.at("Q"), m_dyn_grid_fields.at("FQ"));

  m_d2p_remapper->register_field(m_dyn_grid_fields.at("vtheta_dp"),m_ref_grid_fields.at("T_mid"));
  m_d2p_remapper->register_field(m_dyn_grid_fields.at("v"),m_ref_grid_fields.at("horiz_winds"));
  m_d2p_remapper->register_field(m_dyn_grid_fields.at("dp3d"), m_ref_grid_fields.at("pseudo_density"));
  m_d2p_remapper->register_field(m_dyn_grid_fields.at("phinh_i"), m_ref_grid_fields.at("phi_int"));
  m_d2p_remapper->register_field(m_dyn_grid_fields.at("w_i"), m_ref_grid_fields.at("w_int"));
  m_d2p_remapper->register_field(m_dyn_grid_fields.at("ps"), m_ref_grid_fields.at("ps"));
  m_d2p_remapper->register_field(m_dyn_grid_fields.at("Q"), m_ref_grid_fields.at("Q"));

  m_p2d_remapper->registration_ends();
  m_d2p_remapper->registration_ends();

  // Setup, run, and destroy the IC remapper, to remap IC directly into Homme's states
  m_ic_remapper->registration_begins();
  m_ic_remapper->register_field(m_ref_grid_fields.at("w_int"),m_dyn_grid_fields.at("w_i"));
  m_ic_remapper->register_field(m_ref_grid_fields.at("phi_int"),m_dyn_grid_fields.at("phinh_i"));
  m_ic_remapper->register_field(m_ref_grid_fields.at("horiz_winds"),m_dyn_grid_fields.at("v"));
  m_ic_remapper->register_field(m_ref_grid_fields.at("pseudo_density"),m_dyn_grid_fields.at("dp3d"));
  m_ic_remapper->register_field(m_ref_grid_fields.at("ps"),m_dyn_grid_fields.at("ps"));
  m_ic_remapper->register_field(m_ref_grid_fields.at("T_mid"),m_dyn_grid_fields.at("vtheta_dp"));
  m_ic_remapper->register_field(m_ref_grid_fields.at("Q"),m_dyn_grid_fields.at("Q"));
  m_ic_remapper->registration_ends();
  
  m_ic_remapper->remap(true);
  m_ic_remapper = nullptr;

  // Convert T->vtheta_dp (in place).
  constexpr int N = sizeof(Homme::Scalar) / sizeof(Real);
  using KT = KokkosTypes<DefaultDevice>;
  using Pack = ekat::Pack<Real,N>;
  using ColOps = ColumnOps<DefaultDevice,Real>;
  using PF = PhysicsFunctions<DefaultDevice>;

  const auto ps   = m_dyn_grid_fields.at("ps").get_reshaped_view<Real****>();
  const auto dp3d = m_dyn_grid_fields.at("dp3d").get_reshaped_view<Pack*****>();
  const auto v    = m_dyn_grid_fields.at("v").get_reshaped_view<Pack******>();
  const auto w_i    = m_dyn_grid_fields.at("w_i").get_reshaped_view<Pack*****>();
  const auto phinh_i    = m_dyn_grid_fields.at("phinh_i").get_reshaped_view<Pack*****>();
  const auto vtheta_dp = m_dyn_grid_fields.at("vtheta_dp").get_reshaped_view<Pack*****>();
  const auto Q = m_dyn_grid_fields.at("Q").get_reshaped_view<Pack*****>();

  const auto num_elems = dp3d.extent_int(0);
  const auto nlevs  = m_dyn_grid->get_num_vertical_levels();
  const auto policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(num_elems*NP*NP,NVL);
  const int n0 = c.get<Homme::TimeLevel>().n0;
  const int nm1 = c.get<Homme::TimeLevel>().nm1;
  const int np1 = c.get<Homme::TimeLevel>().np1;

  const auto& hvcoord = c.get<Homme::HybridVCoord>();
  const auto ps0 = hvcoord.ps0 * hvcoord.hybrid_ai0;

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
    ColOps::compute_midpoint_values(team,nlevs,p_int,p_mid);
    
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
  const auto qdp = c.get<Homme::Tracers>().qdp;
  const auto q   = c.get<Homme::Tracers>().Q;
  const auto dp  = c.get<Homme::ElementsState>().m_dp3d;
  const int n0_qdp = c.get<Homme::TimeLevel>().n0_qdp;
  const int num_tracers = qdp.extent_int(2);
  Kokkos::parallel_for(Kokkos::RangePolicy<>(0,num_elems*num_tracers*NP*NP*NVL),
                       KOKKOS_LAMBDA (const int idx) {
    const int ie =  idx / (num_tracers*NP*NP*NVL);
    const int iq = (idx / (NP*NP*NVL)) % num_tracers;
    const int ip = (idx / (NP*NVL)) % NP;
    const int jp = (idx / NVL) % NP;
    const int k  = idx % NVL;

    qdp(ie,n0_qdp,iq,ip,jp,k) = q(ie,iq,ip,jp,k) * dp(ie,n0,ip,jp,k);
  });

  prim_init_model_f90 ();

  // Some checks: nobody should claim to be a provider for dp or w_i.
  EKAT_REQUIRE_MSG (
      m_ref_grid_fields.at("pseudo_density").get_header().get_tracking().get_providers().size()==1,
      "Error! Someone other than Dynamics is trying to update the pseudo_density.\n");
  EKAT_REQUIRE_MSG (
      m_ref_grid_fields.at("w_int").get_header().get_tracking().get_providers().size()==1,
      "Error! Someone other than Dynamics is trying to update the vertical velocity.\n");
}

void HommeDynamics::run_impl (const Real dt)
{
  try {
    // Prepare inputs for homme
    Kokkos::fence();
    homme_pre_process (dt);

    // Run hommexx
    Kokkos::fence();
    prim_run_f90 (dt);

    // Post process Homme's output, to produce what the rest of Atm expects
    Kokkos::fence();
    homme_post_process ();

    // Get a copy of the current timestamp (at the beginning of the step) and
    // advance it, updating the p3 fields.
    auto ts = timestamp();
    ts += dt;
    for (auto& it : m_ref_grid_fields) {
      it.second.get_header().get_tracking().update_time_stamp(ts);
    }
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

void HommeDynamics::set_required_field_impl (const Field<const Real>& f) {
  // Since all inputs are also outputs, we don't store a copy of f here.
  // Instead, we'll store a copy of it when it's given to us during
  // set_computed_field, since we'll have access to a non-const version
  // then, which we need for remapping purposes.
  // The only exception is qv, but qv is available via Q anyays. We only
  // "require" it to guarantee the tracers group is non-empty.

  // Add myself as customer to the field.
  this->add_me_as_customer(f);
}

void HommeDynamics::set_computed_field_impl (const Field<Real>& f) {
  const auto& name = f.get_header().get_identifier().name();

  m_ref_grid_fields.emplace(name,f);

  // Add myself as provider for the field
  this->add_me_as_provider(f);
}

void HommeDynamics::homme_pre_process (const Real dt) {
  // T and uv tendencies are backed out on the ref grid.
  // Homme takes care of turning the FT tendency into a tendency for VTheta_dp.
  using KT = KokkosTypes<DefaultDevice>;

  constexpr int N = sizeof(Homme::Scalar) / sizeof(Real);
  using Pack = ekat::Pack<Real,N>;

  const int ncols = m_ref_grid->get_num_local_dofs();
  const int nlevs = m_ref_grid->get_num_vertical_levels();
  const int npacks = ekat::PackInfo<N>::num_packs(nlevs);

  auto T  = m_ref_grid_fields.at("T_mid").get_reshaped_view<Pack**>();
  auto FT = m_ref_grid_fields.at("T_mid_prev").get_reshaped_view<Pack**>();
  auto v  = m_ref_grid_fields.at("horiz_winds").get_reshaped_view<Pack***>();
  auto w  = m_ref_grid_fields.at("w_int").get_reshaped_view<Pack**>();
  auto FM = m_ref_grid_fields.at("horiz_winds_prev").get_reshaped_view<Pack***>();

  // If there are other atm procs updating the vertical velocity,
  // then we need to compute forcing for w as well
  const bool has_w_forcing = m_ref_grid_fields.at("w_int").get_header().get_tracking().get_providers().size()>1;
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
  const int qsize = Q.extent_int(1);
  const auto n0_qdp = tl.n0_qdp;

  // At this point, FQ contains Qnew (coming from physics).
  // Depending on ftype, we are going to modify it.
  //  ftype=0: FQ = dp*(Qnew-Qold) / dt
  //  ftype=2: FQ = dp*Qnew
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
      break;
    case ForcingAlg::FORCING_2:
      // Hard adjustment of qdp
      Kokkos::parallel_for(Kokkos::RangePolicy<>(0,Q.size()),KOKKOS_LAMBDA(const int idx) {
        const int ie = idx / (qsize*NP*NP*NVL);
        const int iq = (idx / (NP*NP*NVL)) % qsize;
        const int ip = (idx / (NP*NVL)) % NP;
        const int jp = (idx / NVL) % NP;
        const int k  =  idx % NVL;

        // fq is currently storing q_new
              auto& qdp = Qdp(ie,n0_qdp,iq,ip,jp,k);
        const auto& dp  = dp3d(ie,n0,ip,jp,k);
        const auto& fq  = FQ(ie,iq,ip,jp,k);

        qdp = fq*dp;
      });
      Kokkos::fence();

      ff.tracers_forcing(dt,n0,n0_qdp,true,params.moisture);

    default:
      EKAT_ERROR_MSG ("Error! Unexpected/unsupported forcing algorithm.\n"
                      "  ftype: " + std::to_string(Homme::etoi(ftype)) + "\n");
  }
}

void HommeDynamics::homme_post_process () {
  // Remap outputs to ref grid
  m_d2p_remapper->remap(true);

  using KT = KokkosTypes<DefaultDevice>;
  constexpr int N = sizeof(Homme::Scalar) / sizeof(Real);
  using Pack = ekat::Pack<Real,N>;
  using ColOps = ColumnOps<DefaultDevice,Real>;
  using PF = PhysicsFunctions<DefaultDevice>;

  // Convert VTheta_dp->T, store T,uv, and possibly w in FT, FM,
  // compute p_int on ref grid.
  const auto ps_view = m_ref_grid_fields.at("ps").get_reshaped_view<Real*>();
  const auto dp_view = m_ref_grid_fields.at("pseudo_density").get_reshaped_view<Pack**>();
  const auto p_mid_view = m_ref_grid_fields.at("p_mid").get_reshaped_view<Pack**>();
  const auto p_int_view = m_ref_grid_fields.at("p_int").get_reshaped_view<Pack**>();
  const auto Q_view  = m_ref_grid_fields.at("Q").get_reshaped_view<Pack***>();

  const auto T_view  = m_ref_grid_fields.at("T_mid").get_reshaped_view<Pack**>();
  const auto v_view  = m_ref_grid_fields.at("horiz_winds").get_reshaped_view<Pack***>();
  const auto w_view  = m_ref_grid_fields.at("w_int").get_reshaped_view<Pack**>();
  const auto T_prev_view = m_ref_grid_fields.at("T_mid_prev").get_reshaped_view<Pack**>();
  const auto v_prev_view = m_ref_grid_fields.at("horiz_winds_prev").get_reshaped_view<Pack***>();
  const auto w_prev_view = m_ref_grid_fields.at("w_int_prev").get_reshaped_view<Pack**>();

  const auto ncols = m_ref_grid->get_num_local_dofs();
  const auto nlevs = m_ref_grid->get_num_vertical_levels();
  const auto npacks= ekat::PackInfo<N>::num_packs(nlevs);
  const auto policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(ncols,npacks);

  // If there are other atm procs updating the vertical velocity,
  // then we need to store w_int_old (to compute forcing for w next iteration)
  const bool has_w_forcing = m_ref_grid_fields.at("w_int").get_header().get_tracking().get_providers().size()>1;

  constexpr auto P0 = physics::Constants<Real>::P0;
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const KT::MemberType& team) {
    const int& icol = team.league_rank();

    // Compute p_int and p_mid
    auto dp = ekat::subview(dp_view,icol);
    auto p_mid = ekat::subview(p_mid_view,icol);
    auto p_int = ekat::subview(p_int_view,icol);

    ColOps::column_scan<true>(team,nlevs,dp,p_int,P0);
    ColOps::compute_midpoint_values(team,nlevs,p_int,p_mid);
    
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
create_dyn_field (const std::string& name,
                  const std::vector<FieldTag>& tags,
                  const std::vector<int>& dims) {
  using namespace ekat::units;
  FieldIdentifier id(name,FieldLayout{tags,dims},Units::nondimensional(),m_dyn_grid->name());

  const auto lt = get_layout_type(id.get_layout().tags());

  // Only request packed field for 3d quantities
  int pack_size = 1;
  if (lt==LayoutType::Scalar3D || lt==LayoutType::Vector3D || lt==LayoutType::Tensor3D) {
    pack_size = sizeof(Homme::Scalar) / sizeof(Real);
  }

  field_type f(id);
  f.get_header().get_alloc_properties().request_allocation<Real>(pack_size);
  f.allocate_view();

  m_dyn_grid_fields[name] = f;
}

} // namespace scream
