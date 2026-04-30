#include "potential_temperature.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"

namespace scream
{

PotentialTemperature::
PotentialTemperature (const ekat::Comm& comm, const ekat::ParameterList& params,
                      const std::shared_ptr<const AbstractGrid>& grid)
 : AbstractDiagnostic(comm,params,grid)
{
  EKAT_REQUIRE_MSG(params.isParameter("temperature_kind"),
      "Error! PotentialTemperature requires 'temperature_kind' in its input parameters.\n");
  
  auto pt_type = params.get<std::string>("temperature_kind");

  if (pt_type=="Tot"){
    m_ptype = "PotentialTemperature";
  } else if (pt_type=="Liq") {
    m_ptype = "LiqPotentialTemperature";
  } else {
    EKAT_ERROR_MSG (
        "Error! Invalid choice for 'TemperatureKind' in PotentialTemperature.\n"
        "  - input value: " + pt_type + "\n"
        "  - valid values: Tot, Liq\n");
  }

  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  // The fields required for this diagnostic to be computed
  m_field_in_names.push_back("T_mid");
  m_field_in_names.push_back("p_mid");
  if (m_ptype=="LiqPotentialTemperature") {
    m_field_in_names.push_back("qc");
  }

  // Construct and allocate the diagnostic field
  auto diag_layout = m_grid->get_3d_scalar_layout(LEV);
  FieldIdentifier fid (m_ptype, diag_layout, K, m_grid->name());
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.allocate_view();
}

void PotentialTemperature::compute_impl()
{
  using KT      = KokkosTypes<DefaultDevice>;
  using PF      = PhysicsFunctions<DefaultDevice>;
  using MDRange = Kokkos::MDRangePolicy<typename KT::ExeSpace,Kokkos::Rank<2>>;

  bool is_liq = (m_ptype=="LiqPotentialTemperature");

  auto theta = m_diagnostic_output.get_view<Real**>();
  auto T_mid = m_fields_in.at("T_mid").get_view<const Real**>();
  auto p_mid = m_fields_in.at("p_mid").get_view<const Real**>();
  auto q_mid = is_liq ? m_fields_in.at("qc").get_view<const Real**>() : decltype(p_mid){};

  int ncols = m_grid->get_num_local_dofs();
  int nlevs = m_grid->get_num_vertical_levels();
  MDRange policy({0,0},{ncols,nlevs});
  auto lambda = KOKKOS_LAMBDA (const int icol, const int ilev) {
    auto temp = PF::calculate_theta_from_T(T_mid(icol,ilev),p_mid(icol,ilev));
    if (is_liq) {
      // Liquid potential temperature (consistent with how it is calculated in SHOC)
      theta(icol,ilev) = PF::calculate_thetal_from_theta(temp,T_mid(icol,ilev),q_mid(icol,ilev));
    } else {
      // The total potential temperature
      theta(icol,ilev) = temp;
    }
  };

  Kokkos::parallel_for("PotentialTemperature", policy, lambda);
}

} //namespace scream
