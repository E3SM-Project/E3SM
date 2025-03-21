#include "diagnostics/zonal_avg.hpp"

#include "share/field/field_utils.hpp"

namespace scream {

template <typename WeightType>
void ZonalAvgDiag::compute_zonal_sum(const Field &field,
  const WeightType &weight, const Field &lat, const Field &result,
  const ekat::Comm &comm)
{
  auto result_layout = result.get_header().get_identifier().get_layout();
  const int lat_num = result_layout.dim(0);
  const int ncols = field.get_header().get_identifier().get_layout().dim(0);

  const Real lat_delta = sp(180.0) / lat_num;

  auto weight_view = get_view(weight);
  auto field_view = field.get_view<const Real *>();
  auto lat_view = lat.get_view<const Real *>();
  auto result_view = result.get_view<Real *>();
  using KT = ekat::KokkosTypes<DefaultDevice>;
  using TeamPolicy = Kokkos::TeamPolicy<Field::device_t::execution_space>;
  using TeamMember = typename TeamPolicy::member_type;
  using ESU = ekat::ExeSpaceUtils<typename KT::ExeSpace>;
  TeamPolicy team_policy = ESU::get_default_team_policy(lat_num, ncols);
  Kokkos::parallel_for("compute_zonal_sum_" + field.name(), team_policy,
    KOKKOS_LAMBDA(const TeamMember& tm) {
      const int lat_i = tm.league_rank();
      const Real lat_lower = sp(-90.0) + lat_i * lat_delta;
      const Real lat_upper = lat_lower + lat_delta;
      Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, ncols),
        KOKKOS_LAMBDA(int i, Real& val){
          // TODO: check if some conditional is ok here instead of flag
          int flag = (lat_lower <= lat_view(i)) && (lat_view(i) < lat_upper);
          val += flag * weight_view(i) * field_view(i);
        }, result_view(lat_i));
    });

//if(m_comm) {
  // TODO: use device-side MPI calls
  // TODO: the dev ptr causes problems; revisit this later
  // TODO: doing cuda-aware MPI allreduce would be ~10% faster
  Kokkos::fence();
  result.sync_to_host();
  comm.all_reduce(result.template get_internal_view_data<Real, Host>(),
    result_layout.size(), MPI_SUM);
  result.sync_to_dev();
//}
}

ZonalAvgDiag::ZonalAvgDiag(const ekat::Comm &comm,
                           const ekat::ParameterList &params)
    : AtmosphereDiagnostic(comm, params) {
  const auto &field_name = m_params.get<std::string>("field_name");
  m_diag_name = field_name + "_zonal_avg";

  // TODO: consider removing this and replacing with static 180
  m_lat_num = params.isParameter("num_lat_vals") ?
    params.get<int>("num_lat_vals") : 180;
}

void ZonalAvgDiag::set_grids(
    const std::shared_ptr<const GridsManager> grids_manager) {
  const auto &field_name = m_params.get<std::string>("field_name");
  const auto &grid_name = m_params.get<std::string>("grid_name");
  const GridsManager::grid_ptr_type grid = grids_manager->get_grid("Physics");

  add_field<Required>(field_name, grid_name);

  m_area = grid->get_geometry_data("area");
  m_lat = grid->get_geometry_data("lat");
}

void ZonalAvgDiag::initialize_impl(const RunType /*run_type*/) {
  // TODO: auto everything
  using FieldIdentifier = FieldHeader::identifier_type;
  using FieldLayout = FieldIdentifier::layout_type;
  using ShortFieldTagsNames::COL;
  const Field &field = get_fields_in().front();
  const FieldIdentifier &field_id = field.get_header().get_identifier();
  const FieldLayout &field_layout = field_id.get_layout();

  // TODO: support 3D and >3D fields
  EKAT_REQUIRE_MSG(field_layout.rank() >= 1 && field_layout.rank() <= 1,
                   "Error! Field rank not supported by ZonalAvgDiag.\n"
                   " - field name: " +
                       field_id.name() +
                       "\n"
                       " - field layout: " +
                       field_layout.to_string() + "\n");
  EKAT_REQUIRE_MSG(field_layout.tags()[0] == COL,
                   "Error! ZonalAvgDiag diagnostic expects a layout starting "
                   "with the 'COL' tag.\n"
                   " - field name  : " +
                       field_id.name() +
                       "\n"
                       " - field layout: " +
                       field_layout.to_string() + "\n");

  // TODO: replace COL with something else
  FieldLayout diagnostic_layout({COL}, {m_lat_num});
  FieldIdentifier diagnostic_id(m_diag_name, diagnostic_layout,
    field_id.get_units(), field_id.get_grid_name());
  m_diagnostic_output = Field(diagnostic_id);
  m_diagnostic_output.allocate_view();

  // compute zonal area
  FieldIdentifier zonal_area_id("zonal area", diagnostic_layout,
    m_area.get_header().get_identifier().get_units(), field_id.get_grid_name());
  m_zonal_area = Field(zonal_area_id);
  m_zonal_area.allocate_view();
  IdentityField one;
  compute_zonal_sum(m_area, one, m_lat, m_zonal_area, m_comm);
}

void ZonalAvgDiag::compute_diagnostic_impl() {
  const auto &field = get_fields_in().front();

  compute_zonal_sum(field, m_area, m_lat, m_diagnostic_output, m_comm);
  auto zonal_area_view = m_zonal_area.get_view<const Real *>();
  auto output_view = m_diagnostic_output.get_view<Real *>();
  for (int i=0; i < m_lat_num; i++) {
    output_view(i) /= zonal_area_view(i);
  }
}

}  // namespace scream
