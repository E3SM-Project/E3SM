#include "eamxx_water_tracers_process_interface.hpp"
#include "eamxx_water_tracers_functions.hpp"
#include "../water_isotopes/eamxx_water_isotopes_ocean_flux.hpp"

#include "share/physics/eamxx_common_physics_functions.hpp"

#include <ekat_assert.hpp>
#include <ekat_pack.hpp>
#include <ekat_units.hpp>

namespace scream
{

namespace {

namespace wt = water_tracers;

// Kernels evaluate the functions in eamxx_water_tracers_functions.hpp one
// column at a time, on the bottom level only, so a pack of one is enough.
using Pack1 = ekat::Pack<Real,1>;

} // anonymous namespace

// =========================================================================================
WaterTracers::WaterTracers (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  water_tracers::RegistryOptions opts;
#ifdef EAMXX_HAS_WATER_ISOTOPES
  water_isotopes::WaterIsotopes::enable(opts);
#endif

  const std::string key = "tracers_file";
  EKAT_REQUIRE_MSG(m_params.isParameter(key) and m_params.isType<std::string>(key) and
                   m_params.get<std::string>(key)!="",
      "Error! Process '" + m_params.name() + "' is in atm_procs_list, but its parameter\n"
      "  'tracers_file' is missing or empty. Set it to the water tracer YAML file.\n");

  m_specs = water_tracers::parse_tracers_file(m_params.get<std::string>(key), opts);

  // Reading the isotope formulation options is safe even with no isotope
  // tracers: they are simply never read by compute_isotope_fluxes below.
  m_isotopes.configure(m_params);
  m_isotopes.resolve(m_specs);
}

// =========================================================================================
void WaterTracers::create_requests ()
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  m_grid = m_grids_manager->get_grid("physics");
  const auto& grid_name = m_grid->name();
  m_num_cols = m_grid->get_num_local_dofs();
  m_num_levs = m_grid->get_num_vertical_levels();

  constexpr int ps = SCREAM_PACK_SIZE;
  constexpr auto m2 = pow(m, 2);
  const auto scalar2d     = m_grid->get_2d_scalar_layout();
  const auto scalar3d_mid = m_grid->get_3d_scalar_layout(LEV);

  // Bottom-level atmosphere state (for qv and the air density)
  add_tracer<Required>("qv", m_grid, kg/kg, ps);
  add_field<Required>("T_mid",          scalar3d_mid, K,  grid_name, ps);
  add_field<Required>("p_mid",          scalar3d_mid, Pa, grid_name, ps);
  add_field<Required>("pseudo_density", scalar3d_mid, Pa, grid_name, ps);

  // Ocean surface state from the coupler (atmosphere_surface_coupling_importer)
  add_field<Required>("ocnfrac",   scalar2d, none,  grid_name);
  add_field<Required>("sst",       scalar2d, K,     grid_name);
  add_field<Required>("ocn_ssq",   scalar2d, kg/kg, grid_name);
  add_field<Required>("ocn_re",    scalar2d, none,  grid_name);
  add_field<Required>("ocn_ustar", scalar2d, m/s,   grid_name);

  // The tracers and their diagnostic ocean fluxes. The flux group is not used
  // yet; it is where a later process will find the fluxes to apply.
  bool any_ratio_init = false;
  for (const auto& spec : m_specs) {
    add_tracer<Updated>(spec.name, m_grid, kg/kg, ps, TracerAdvection::DynamicsAndTurbulence);
    add_field<Computed>(spec.name + "_surf_flux", scalar2d, kg/m2/s, grid_name, "tracer_surface_fluxes");
    any_ratio_init |= spec.init==water_tracers::InitType::Ratio;
  }

  // [pgN initial conditions]
  // With a pgN physics grid (e.g. physics_pg2 in CIME ne*pg2 runs), the driver
  // reads initial conditions only for fields on physics_gll and leaves the pgN
  // copies unset (AtmosphereDriver::set_initial_conditions, "fvphyshack").
  // HommeDynamics then remaps the physics_gll "tracers" group to the dynamics
  // grid, and the dynamics state to pgN, at the end of its own initialize_impl.
  // This process runs before homme, so on pgN its qv is still unset during
  // initialize_impl, and anything it wrote to its pgN tracers would be
  // overwritten by homme's remap.
  // Hence, in a pgN run we also request qv and the tracers on physics_gll, set
  // the init: ratio tracers there (see initialize_impl), and homme carries
  // them to the dynamics grid and to pgN, like every other tracer IC. The
  // field manager accepts a tracer registered on both physics grids, and the
  // driver gives the physics_gll copies the IC values before any process is
  // initialized. Like homme, we drop the physics_gll handles at the end of
  // initialize_impl, so the driver can free them.
  // A run is pgN iff a physics_gll grid exists and is not the physics grid
  // itself; in particular, a standalone grid named physics_gll and aliased to
  // physics takes the direct path.
  if (any_ratio_init and m_grids_manager->has_grid("physics_gll") and grid_name!="physics_gll") {
    m_ic_grid = m_grids_manager->get_grid("physics_gll");
    add_tracer<Required>("qv", m_ic_grid, kg/kg, ps);
    for (const auto& spec : m_specs) {
      add_tracer<Updated>(spec.name, m_ic_grid, kg/kg, ps, TracerAdvection::DynamicsAndTurbulence);
    }
  }

  m_sfc_f     = KT::view_1d<Real>("water_tracers_sfc_f",     m_num_cols);
  m_sfc_rho   = KT::view_1d<Real>("water_tracers_sfc_rho",   m_num_cols);
  m_sfc_ustar = KT::view_1d<Real>("water_tracers_sfc_ustar", m_num_cols);
  m_sfc_re    = KT::view_1d<Real>("water_tracers_sfc_re",    m_num_cols);
  m_sfc_ssq   = KT::view_1d<Real>("water_tracers_sfc_ssq",   m_num_cols);
  m_sfc_qv    = KT::view_1d<Real>("water_tracers_sfc_qv",    m_num_cols);
  m_sfc_evap  = KT::view_1d<Real>("water_tracers_sfc_evap",  m_num_cols);
  m_sfc_sst   = KT::view_1d<Real>("water_tracers_sfc_sst",   m_num_cols);
}

// =========================================================================================
void WaterTracers::initialize_impl (const RunType run_type)
{
  // On restarts the tracers come from the restart file. See [pgN initial
  // conditions] in create_requests for why the ratio ICs may go on physics_gll.
  if (run_type==RunType::Initial) {
    apply_ratio_initial_conditions(m_ic_grid ? m_ic_grid->name() : m_grid->name());
  }

  if (m_ic_grid) {
    const auto& ic_grid_name = m_ic_grid->name();
    remove_field("qv", ic_grid_name);
    for (const auto& spec : m_specs) {
      remove_field(spec.name, ic_grid_name);
    }
  }
}

// =========================================================================================
void WaterTracers::apply_ratio_initial_conditions (const std::string& grid_name)
{
  using MDPolicy = Kokkos::MDRangePolicy<KT::ExeSpace, Kokkos::Rank<2>>;

  const auto grid = m_grids_manager->get_grid(grid_name);
  const MDPolicy policy({0, 0}, {grid->get_num_local_dofs(), grid->get_num_vertical_levels()});

  const auto q_v = get_field_in("qv", grid_name).get_view<const Real**>();
  for (const auto& spec : m_specs) {
    if (spec.init!=water_tracers::InitType::Ratio) continue;
    const auto q = get_field_out(spec.name, grid_name).get_view<Real**>();
    const Real r_init = spec.r_init;
    Kokkos::parallel_for("water_tracers_ratio_ic_" + spec.name, policy,
                         KOKKOS_LAMBDA(const int icol, const int ilev) {
      q(icol, ilev) = r_init*q_v(icol, ilev);
    });
  }
  Kokkos::fence();
}

// =========================================================================================
void WaterTracers::run_impl (const double /* dt */)
{
  // Each step below is a column kernel over the bottom level only, so a
  // RangePolicy over columns is used rather than a TeamPolicy. There is one
  // kernel per tracer: tracers are few, and each reads different fields.
  compute_bulk_state();
  compute_generic_fluxes();
  compute_isotope_fluxes();
}

// =========================================================================================
void WaterTracers::compute_bulk_state ()
{
  using PF = PhysicsFunctions<DefaultDevice>;

  const auto& gn  = m_grid->name();
  const int  kbot = m_num_levs - 1;

  const auto qv      = get_field_in("qv",             gn).get_view<const Real**>();
  const auto T_mid   = get_field_in("T_mid",          gn).get_view<const Real**>();
  const auto p_mid   = get_field_in("p_mid",          gn).get_view<const Real**>();
  const auto dp      = get_field_in("pseudo_density", gn).get_view<const Real**>();
  const auto ocnfrac = get_field_in("ocnfrac",        gn).get_view<const Real*>();
  const auto sst     = get_field_in("sst",            gn).get_view<const Real*>();
  const auto ssq     = get_field_in("ocn_ssq",        gn).get_view<const Real*>();
  const auto re      = get_field_in("ocn_re",         gn).get_view<const Real*>();
  const auto ustar   = get_field_in("ocn_ustar",      gn).get_view<const Real*>();

  const auto sfc_f     = m_sfc_f;
  const auto sfc_rho   = m_sfc_rho;
  const auto sfc_ustar = m_sfc_ustar;
  const auto sfc_re    = m_sfc_re;
  const auto sfc_ssq   = m_sfc_ssq;
  const auto sfc_qv    = m_sfc_qv;
  const auto sfc_evap  = m_sfc_evap;
  const auto sfc_sst   = m_sfc_sst;

  Kokkos::parallel_for("water_tracers_bulk_state", KT::RangePolicy(0, m_num_cols),
                       KOKKOS_LAMBDA(const int icol) {
    // Same density the surface coupling exporter sends to the coupler (Sa_dens)
    const Real q   = qv(icol, kbot);
    const Real dz  = PF::calculate_dz(dp(icol, kbot), p_mid(icol, kbot), T_mid(icol, kbot), q);
    const Real rho = PF::calculate_density(dp(icol, kbot), dz);

    Pack1 f(ocnfrac(icol)), u(ustar(icol)), r(re(icol)), s(ssq(icol));
    const Pack1 e = wt::ocean_evaporation(f, u, r, s, Pack1(rho), Pack1(q));

    sfc_f(icol)     = f[0];
    sfc_rho(icol)   = rho;
    sfc_ustar(icol) = u[0];
    sfc_re(icol)    = r[0];
    sfc_ssq(icol)   = s[0];
    sfc_qv(icol)    = q;
    sfc_evap(icol)  = e[0];
    sfc_sst(icol)   = sst(icol);
  });
}

// =========================================================================================
void WaterTracers::compute_generic_fluxes ()
{
  const auto& gn  = m_grid->name();
  const int  kbot = m_num_levs - 1;

  const auto sfc_f     = m_sfc_f;
  const auto sfc_rho   = m_sfc_rho;
  const auto sfc_ustar = m_sfc_ustar;
  const auto sfc_re    = m_sfc_re;
  const auto sfc_ssq   = m_sfc_ssq;
  const auto sfc_qv    = m_sfc_qv;
  const auto sfc_evap  = m_sfc_evap;

  for (const auto& spec : m_specs) {
    if (spec.type!=water_tracers::TracerType::Generic) continue;

    const auto q_i  = get_field_in(spec.name, gn).get_view<const Real**>();
    const auto flux = get_field_out(spec.name + "_surf_flux", gn).get_view<Real*>();
    const auto kind  = spec.ocean_flux;
    const Real r_ocn = spec.r_ocn;

    Kokkos::parallel_for("water_tracers_ocean_flux_" + spec.name, KT::RangePolicy(0, m_num_cols),
                         KOKKOS_LAMBDA(const int icol) {
      const Pack1 f_i = wt::generic_ocean_flux(kind, r_ocn,
          Pack1(sfc_f(icol)), Pack1(sfc_rho(icol)), Pack1(sfc_ustar(icol)), Pack1(sfc_re(icol)),
          Pack1(sfc_ssq(icol)), Pack1(sfc_qv(icol)), Pack1(q_i(icol, kbot)), Pack1(sfc_evap(icol)));
      flux(icol) = (f_i + wt::land_ice_tracer_flux<Pack1>())[0];
    });
  }
}

// =========================================================================================
void WaterTracers::compute_isotope_fluxes ()
{
  // Local copies, so the lambda captures views and values rather than `this`.
  const auto sfc_f     = m_sfc_f;
  const auto sfc_rho   = m_sfc_rho;
  const auto sfc_ustar = m_sfc_ustar;
  const auto sfc_re    = m_sfc_re;
  const auto sfc_ssq   = m_sfc_ssq;
  const auto sfc_qv    = m_sfc_qv;
  const auto sfc_evap  = m_sfc_evap;
  const auto sfc_sst   = m_sfc_sst;
  const auto iso_params = m_isotopes.parameters();
  const int  ilev_bot  = m_num_levs - 1;

  // Tracers may also live on another grid (e.g. physics_gll for initial
  // conditions), so look them up on this process' physics grid explicitly.
  const auto& grid_name = m_grid->name();

  // One kernel per isotope: there are at most a handful, and each reads a
  // different tracer field. No-op when there are no isotope tracers.
  for (const auto& iso : m_isotopes.tracers()) {
    const auto q_i  = get_field_in(iso.name, grid_name).get_view<const Real**>();
    const auto flux = get_field_out(iso.name + "_surf_flux", grid_name).get_view<Real*>();
    const auto species = iso.species;
    const Real r_ocn   = iso.r_ocn;

    Kokkos::parallel_for("water_isotopes_ocean_flux_" + iso.name,
                         KT::RangePolicy(0, m_num_cols),
                         KOKKOS_LAMBDA(const int icol) {
      flux(icol) = water_isotopes::isotope_ocean_flux(
          Pack1(sfc_f(icol)), Pack1(sfc_rho(icol)), Pack1(sfc_ustar(icol)),
          Pack1(sfc_re(icol)), Pack1(sfc_ssq(icol)), Pack1(sfc_qv(icol)),
          Pack1(q_i(icol, ilev_bot)), Pack1(sfc_evap(icol)), Pack1(sfc_sst(icol)),
          r_ocn, species, iso_params)[0];
    });
  }
}

} // namespace scream
