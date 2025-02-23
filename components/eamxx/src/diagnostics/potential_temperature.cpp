#include "diagnostics/potential_temperature.hpp"
#include "share/util/eamxx_common_physics_functions.hpp"

namespace scream
{

// =========================================================================================
PotentialTemperatureDiagnostic::PotentialTemperatureDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereDiagnostic(comm,params)
{
  EKAT_REQUIRE_MSG(params.isParameter("Temperature Kind"),
      "Error! PotentialTemperatureDiagnostic requires 'Temperature Kind' in its input parameters.\n");
  
  auto pt_type = params.get<std::string>("Temperature Kind");

  if (pt_type=="Tot"){
    m_ptype = "PotentialTemperature";
  } else if (pt_type=="Liq") {
    m_ptype = "LiqPotentialTemperature";
  } else {
    EKAT_ERROR_MSG (
        "Error! Invalid choice for 'TemperatureKind' in PotentialTemperatureDiagnostic.\n"
        "  - input value: " + pt_type + "\n"
        "  - valid values: Tot, Liq\n");
  }
}

// =========================================================================================
void PotentialTemperatureDiagnostic::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  auto grid  = grids_manager->get_grid("Physics");
  const auto& grid_name = grid->name();
  m_num_cols = grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = grid->get_num_vertical_levels();  // Number of levels per column

  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols,m_num_levs} };

  // The fields required for this diagnostic to be computed
  add_field<Required>("T_mid", scalar3d_layout_mid, K,  grid_name);
  add_field<Required>("p_mid", scalar3d_layout_mid, Pa, grid_name);
  if (m_ptype=="LiqPotentialTemperature") {
    add_field<Required>("qc",    scalar3d_layout_mid, kg/kg,  grid_name);
  }

  // Construct and allocate the diagnostic field
  FieldIdentifier fid (m_ptype, scalar3d_layout_mid, K, grid_name);
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.allocate_view();
}
// =========================================================================================
void PotentialTemperatureDiagnostic::compute_diagnostic_impl()
{
  using PF = scream::PhysicsFunctions<DefaultDevice>;

  bool is_liq = (m_ptype=="LiqPotentialTemperature");

  auto theta = m_diagnostic_output.get_view<Real**>();
  auto T_mid = get_field_in("T_mid").get_view<const Real**>();
  auto p_mid = get_field_in("p_mid").get_view<const Real**>();
  decltype(p_mid) q_mid;
  if (is_liq)
   q_mid = get_field_in("qc").get_view<const Real**>();

  int nlevs = m_num_levs;
  Kokkos::parallel_for("PotentialTemperatureDiagnostic",
                       Kokkos::RangePolicy<>(0,m_num_cols*m_num_levs),
                       KOKKOS_LAMBDA (const int& idx) {
      const int icol = idx / nlevs;
      const int ilev = idx % nlevs;
      auto temp = PF::calculate_theta_from_T(T_mid(icol,ilev),p_mid(icol,ilev));
      if (is_liq) {
        // Liquid potential temperature (consistent with how it is calculated in SHOC)
        theta(icol,ilev) = PF::calculate_thetal_from_theta(temp,T_mid(icol,ilev),q_mid(icol,ilev));
      } else {
        // The total potential temperature
        theta(icol,ilev) = temp;
      }
  });
  Kokkos::fence();
}
// =========================================================================================

} //namespace scream
