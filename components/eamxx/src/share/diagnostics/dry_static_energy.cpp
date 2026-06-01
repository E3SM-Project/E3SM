#include "dry_static_energy.hpp"

#include "share/physics/eamxx_common_physics_functions.hpp"
#include "share/core/eamxx_types.hpp"

namespace scream
{

DryStaticEnergy::
DryStaticEnergy (const ekat::Comm& comm, const ekat::ParameterList& params,
                 const std::shared_ptr<const AbstractGrid>& grid)
 : AbstractDiagnostic(comm,params,grid)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  auto m2  = pow(m,2);
  auto s2  = pow(s,2);

  const auto& grid_name = m_grid->name();

  auto scalar2d = m_grid->get_2d_scalar_layout();
  auto scalar3d = m_grid->get_3d_scalar_layout(LEV);

  // The fields required for this diagnostic to be computed
  m_field_in_names.push_back("T_mid");
  m_field_in_names.push_back("height_mid");
  m_field_in_names.push_back("phis");

  // Construct and allocate the diagnostic field
  FieldIdentifier fid (name(), scalar3d, m2/s2, grid_name);
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.allocate_view();
}

void DryStaticEnergy::compute_impl()
{
  using KT      = KokkosTypes<DefaultDevice>;
  using PF      = PhysicsFunctions<DefaultDevice>;
  using MDRange = Kokkos::MDRangePolicy<typename KT::ExeSpace,Kokkos::Rank<2>>;

  const int ncols = m_grid->get_num_local_dofs();
  const int nlevs = m_grid->get_num_vertical_levels();

  const auto& dse   = m_diagnostic_output.get_view<Real**>();
  const auto& T_mid = m_fields_in.at("T_mid").get_view<const Real**>();
  const auto& h_mid = m_fields_in.at("height_mid").get_view<const Real**>();
  const auto& phis  = m_fields_in.at("phis").get_view<const Real*>();

  MDRange policy ({0,0},{ncols,nlevs});
  auto lambda = KOKKOS_LAMBDA(const int icol, const int ilev) {
    dse(icol,ilev) = PF::calculate_dse(T_mid(icol,ilev),h_mid(icol,ilev),phis(icol));
  };

  Kokkos::parallel_for("DryStaticEnergy", policy, lambda);
}

} //namespace scream
