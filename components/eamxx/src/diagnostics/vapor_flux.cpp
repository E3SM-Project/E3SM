#include "diagnostics/vapor_flux.hpp"
#include "physics/share/physics_constants.hpp"

#include <ekat/kokkos/ekat_kokkos_utils.hpp>

namespace scream
{

VaporFluxDiagnostic::
VaporFluxDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereDiagnostic(comm,params)
{
  EKAT_REQUIRE_MSG (params.isParameter("Wind Component"),
      "Error! VaporFluxDiagnostic requires 'Wind Component' in its input parameters.\n");

  const auto& comp = m_params.get<std::string>("Wind Component");
  if (comp=="Zonal") {
    m_component = 0;
  } else if (comp=="Meridional") {
    m_component = 1;
  } else {
    EKAT_ERROR_MSG (
        "Error! Invalid choice for 'Wind Component' in VaporFluxDiagnostic.\n"
        "  - input value: " + comp + "\n"
        "  - valid values: Zonal, Meridional\n");
  }
  m_name = comp + "VapFlux";
}

void VaporFluxDiagnostic::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;

  auto grid  = grids_manager->get_grid("Physics");
  const auto& grid_name = grid->name();
  m_num_cols = grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = grid->get_num_vertical_levels();  // Number of levels per column

  auto scalar2d = grid->get_2d_scalar_layout();
  auto scalar3d = grid->get_3d_scalar_layout(true);
  auto vector3d = grid->get_3d_vector_layout(true,2);

  // The fields required for this diagnostic to be computed
  add_field<Required>("pseudo_density", scalar3d, Pa,    grid_name);
  add_field<Required>("qv",             scalar3d, kg/kg, grid_name);
  add_field<Required>("horiz_winds",    vector3d, m/s,   grid_name);

  // Construct and allocate the diagnostic field
  FieldIdentifier fid (m_name, scalar2d, kg/m/s, grid_name);
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.allocate_view();
}

void VaporFluxDiagnostic::compute_diagnostic_impl()
{
  using PC  = scream::physics::Constants<Real>;
  using KT  = KokkosTypes<DefaultDevice>;
  using MT  = typename KT::MemberType;
  using ESU = ekat::ExeSpaceUtils<typename KT::ExeSpace>;

  constexpr Real g = PC::gravit;

  const auto diag  = m_diagnostic_output.get_view<Real*>();
  const auto qv    = get_field_in("qv").get_view<const Real**>();
  const auto rho   = get_field_in("pseudo_density").get_view<const Real**>();
  const auto wind  = get_field_in("horiz_winds").get_component(m_component).get_view<const Real**>();

  const auto num_levs = m_num_levs;
  const auto policy = ESU::get_default_team_policy(m_num_cols, m_num_levs);
  Kokkos::parallel_for("Compute " + m_name, policy,
                       KOKKOS_LAMBDA(const MT& team) {
    const int icol = team.league_rank();

    auto qv_icol   = ekat::subview(qv,icol);
    auto rho_icol  = ekat::subview(rho,icol);
    auto wind_icol = ekat::subview(wind,icol);

    Kokkos::parallel_reduce(Kokkos::TeamVectorRange(team, num_levs),
                            [&] (const int& ilev, Real& lsum) {
      lsum += wind_icol(ilev) * qv_icol(ilev) * rho_icol(ilev) / g;
    },diag(icol));
    team.team_barrier();
  });
}

} //namespace scream
