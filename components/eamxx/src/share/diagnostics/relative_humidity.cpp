#include "relative_humidity.hpp"
#include "share/physics/physics_functions.hpp" // also for ETI not on GPUs
#include "share/physics/physics_saturation_impl.hpp"

#include <ekat_pack.hpp>

namespace scream
{

RelativeHumidity::
RelativeHumidity (const ekat::Comm& comm, const ekat::ParameterList& params,
                             const std::shared_ptr<const AbstractGrid>& grid)
 : AbstractDiagnostic(comm,params,grid)
{
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  auto scalar3d = m_grid->get_3d_scalar_layout(LEV);

  m_field_in_names.push_back("T_mid");
  m_field_in_names.push_back("p_dry_mid");
  m_field_in_names.push_back("qv");
  m_field_in_names.push_back("pseudo_density");
  m_field_in_names.push_back("pseudo_density_dry");

  auto diag_layout = m_grid->get_3d_scalar_layout(LEV);
  FieldIdentifier fid (name(), diag_layout, none, m_grid->name());
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.get_header().get_alloc_properties().request_allocation(SCREAM_PACK_SIZE);
  m_diagnostic_output.allocate_view();
}

void RelativeHumidity::compute_impl()
{
  using KT      = KokkosTypes<DefaultDevice>;
  using MDRange = Kokkos::MDRangePolicy<typename KT::ExeSpace,Kokkos::Rank<2>>;
  using physics = scream::physics::Functions<Real, DefaultDevice>;
  using Pack    = ekat::Pack<Real,SCREAM_PACK_SIZE>;

  auto theta     = m_diagnostic_output.get_view<Pack**>();
  auto T_mid     = m_fields_in.at("T_mid").get_view<const Pack**>();
  auto p_dry_mid = m_fields_in.at("p_dry_mid").get_view<const Pack**>();
  auto dp_wet    = m_fields_in.at("pseudo_density").get_view<const Pack**>();
  auto dp_dry    = m_fields_in.at("pseudo_density_dry").get_view<const Pack**>();
  auto qv_mid    = m_fields_in.at("qv").get_view<const Pack**>();
  auto RH        = m_diagnostic_output.get_view<Pack**>();

  const int nlevs = m_grid->get_num_vertical_levels();
  const int ncols = m_grid->get_num_local_dofs();

  MDRange policy({0,0},{ncols,ekat::PackInfo<SCREAM_PACK_SIZE>::num_packs(nlevs)});

  auto lambda = KOKKOS_LAMBDA (const int& icol, const int& ipack) {
    const auto range_pack = ekat::range<Pack>(ipack*Pack::n);
    const auto range_mask = range_pack < nlevs;
    auto qv_sat_l = physics::qv_sat_wet(T_mid(icol,ipack),  p_dry_mid(icol,ipack), true, range_mask, dp_wet(icol,ipack), dp_dry(icol,ipack),
                                         physics::MurphyKoop, "RelativeHumidity::compute_impl");
    RH(icol,ipack) = qv_mid(icol,ipack)/qv_sat_l;
  };
  Kokkos::parallel_for("RelativeHumidity", policy, lambda);
}

} //namespace scream
