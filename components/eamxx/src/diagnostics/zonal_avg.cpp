#include "diagnostics/zonal_avg.hpp"

#include "share/field/field_utils.hpp"

namespace scream {

void ZonalAvgDiag::compute_zonal_sum(const Field &result, const Field &field, const Field &weight,
                                     const Field &lat, const ekat::Comm *comm) {
  auto result_layout       = result.get_header().get_identifier().get_layout();
  const int num_zonal_bins = result_layout.dim(0);
  const int ncols          = field.get_header().get_identifier().get_layout().dim(0);
  const Real lat_delta     = sp(180.0) / num_zonal_bins;

  auto weight_view = weight.get_view<const Real *>();
  auto lat_view    = lat.get_view<const Real *>();
  using KT         = ekat::KokkosTypes<DefaultDevice>;
  using TeamPolicy = Kokkos::TeamPolicy<Field::device_t::execution_space>;
  using TeamMember = typename TeamPolicy::member_type;
  using ESU        = ekat::ExeSpaceUtils<typename KT::ExeSpace>;
  switch (result_layout.rank()) {
  case 1: {
    auto field_view        = field.get_view<const Real *>();
    auto result_view       = result.get_view<Real *>();
    TeamPolicy team_policy = ESU::get_default_team_policy(num_zonal_bins, ncols);
    Kokkos::parallel_for(
        "compute_zonal_sum_" + field.name(), team_policy, KOKKOS_LAMBDA(const TeamMember &tm) {
          const int lat_i      = tm.league_rank();
          const Real lat_lower = sp(-90.0) + lat_i * lat_delta;
          const Real lat_upper = lat_lower + lat_delta;
          Kokkos::parallel_reduce(
              Kokkos::TeamVectorRange(tm, ncols),
              [&](int i, Real &val) {
                // TODO: check if tenary is ok here (if not, multiply by flag)
                int flag = (lat_lower <= lat_view(i)) && (lat_view(i) < lat_upper);
                val += flag ? weight_view(i) * field_view(i) : sp(0.0);
              },
              result_view(lat_i));
        });
  } break;
  case 2: {
    const int d1           = result_layout.dim(1);
    auto field_view        = field.get_view<const Real **>();
    auto result_view       = result.get_view<Real **>();
    TeamPolicy team_policy = ESU::get_default_team_policy(num_zonal_bins * d1, ncols);
    Kokkos::parallel_for(
        "compute_zonal_sum_" + field.name(), team_policy, KOKKOS_LAMBDA(const TeamMember &tm) {
          const int idx        = tm.league_rank();
          const int d1_i       = idx / num_zonal_bins;
          const int lat_i      = idx % num_zonal_bins;
          const Real lat_lower = sp(-90.0) + lat_i * lat_delta;
          const Real lat_upper = lat_lower + lat_delta;
          Kokkos::parallel_reduce(
              Kokkos::TeamVectorRange(tm, ncols),
              [&](int i, Real &val) {
                int flag = (lat_lower <= lat_view(i)) && (lat_view(i) < lat_upper);
                // TODO: check if tenary is ok here (if not, multiply by flag)
                val += flag ? weight_view(i) * field_view(i, d1_i) : sp(0.0);
              },
              result_view(lat_i, d1_i));
        });
  } break;
  case 3: {
    const int d1           = result_layout.dim(1);
    const int d2           = result_layout.dim(2);
    auto field_view        = field.get_view<const Real ***>();
    auto result_view       = result.get_view<Real ***>();
    TeamPolicy team_policy = ESU::get_default_team_policy(num_zonal_bins * d1 * d2, ncols);
    Kokkos::parallel_for(
        "compute_zonal_sum_" + field.name(), team_policy, KOKKOS_LAMBDA(const TeamMember &tm) {
          const int idx        = tm.league_rank();
          const int d1_i       = idx / (num_zonal_bins * d2);
          const int idx2       = idx % (num_zonal_bins * d2);
          const int d2_i       = idx2 / num_zonal_bins;
          const int lat_i      = idx2 % num_zonal_bins;
          const Real lat_lower = sp(-90.0) + lat_i * lat_delta;
          const Real lat_upper = lat_lower + lat_delta;
          Kokkos::parallel_reduce(
              Kokkos::TeamVectorRange(tm, ncols),
              [&](int i, Real &val) {
                int flag = (lat_lower <= lat_view(i)) && (lat_view(i) < lat_upper);
                // TODO: check if tenary is ok here (if not, multiply by flag)
                val += flag ? weight_view(i) * field_view(i, d1_i, d2_i) : sp(0.0);
              },
              result_view(lat_i, d1_i, d2_i));
        });
  } break;
  default:
    EKAT_ERROR_MSG("Error! Unsupported field rank for zonal averages.\n");
  }

  if (comm) {
    // TODO: use device-side MPI calls
    // TODO: the dev ptr causes problems; revisit this later
    // TODO: doing cuda-aware MPI allreduce would be ~10% faster
    Kokkos::fence();
    result.sync_to_host();
    comm->all_reduce(result.template get_internal_view_data<Real, Host>(), result_layout.size(),
                     MPI_SUM);
    result.sync_to_dev();
  }
}

ZonalAvgDiag::ZonalAvgDiag(const ekat::Comm &comm, const ekat::ParameterList &params)
    : AtmosphereDiagnostic(comm, params) {
  const auto &field_name     = m_params.get<std::string>("field_name");
  const auto &num_bins_value = params.get<std::string>("number_of_zonal_bins");
  m_diag_name                = field_name + "_zonal_avg_" + num_bins_value + "_bins";
  m_num_zonal_bins           = std::stoi(num_bins_value);
}

void ZonalAvgDiag::set_grids(const std::shared_ptr<const GridsManager> grids_manager) {
  const auto &field_name = m_params.get<std::string>("field_name");
  const auto &grid_name  = m_params.get<std::string>("grid_name");

  add_field<Required>(field_name, grid_name);
  const GridsManager::grid_ptr_type grid = grids_manager->get_grid(grid_name);
  m_lat                                  = grid->get_geometry_data("lat");
  // area will be scaled in initialize_impl
  m_scaled_area = grid->get_geometry_data("area").clone();
}

void ZonalAvgDiag::initialize_impl(const RunType /*run_type*/) {
  using FieldIdentifier = FieldHeader::identifier_type;
  using FieldLayout     = FieldIdentifier::layout_type;
  using ShortFieldTagsNames::COL, ShortFieldTagsNames::CMP, ShortFieldTagsNames::LEV;
  const Field &field              = get_fields_in().front();
  const FieldIdentifier &field_id = field.get_header().get_identifier();
  const FieldLayout &field_layout = field_id.get_layout();

  EKAT_REQUIRE_MSG(field_layout.rank() >= 1 && field_layout.rank() <= 3,
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

  FieldLayout diagnostic_layout =
      field_layout.clone().strip_dim(COL).prepend_dim({CMP}, {m_num_zonal_bins}, {"bin"});
  FieldIdentifier diagnostic_id(m_diag_name, diagnostic_layout, field_id.get_units(),
                                field_id.get_grid_name());
  m_diagnostic_output = Field(diagnostic_id);
  m_diagnostic_output.allocate_view();

  // allocate zonal area
  const FieldIdentifier &area_id = m_scaled_area.get_header().get_identifier();
  FieldLayout zonal_area_layout({CMP}, {m_num_zonal_bins}, {"bin"});
  FieldIdentifier zonal_area_id("zonal area", zonal_area_layout, area_id.get_units(),
                                area_id.get_grid_name());
  Field zonal_area(zonal_area_id);
  zonal_area.allocate_view();

  // compute zonal area
  FieldLayout ones_layout = area_id.get_layout().clone();
  FieldIdentifier ones_id("ones", ones_layout, area_id.get_units(), area_id.get_grid_name());
  Field ones(ones_id);
  ones.allocate_view();
  ones.deep_copy(1.0);
  compute_zonal_sum(zonal_area, m_scaled_area, ones, m_lat, &m_comm);

  // scale area by 1 / zonal area
  using RangePolicy     = Kokkos::RangePolicy<Field::device_t::execution_space>;
  const Real lat_delta  = sp(180.0) / m_num_zonal_bins;
  const int ncols       = field_layout.dim(0);
  auto lat_view         = m_lat.get_view<const Real *>();
  auto zonal_area_view  = zonal_area.get_view<const Real *>();
  auto scaled_area_view = m_scaled_area.get_view<Real *>();
  Kokkos::parallel_for(
      "scale_area_by_zonal_area_" + field.name(), RangePolicy(0, ncols),
      KOKKOS_LAMBDA(const int &i) {
        const int lat_i = (lat_view(i) + sp(90.0)) / lat_delta;
        scaled_area_view(i) /= zonal_area_view(lat_i);
      });
}

void ZonalAvgDiag::compute_diagnostic_impl() {
  const auto &field = get_fields_in().front();
  compute_zonal_sum(m_diagnostic_output, field, m_scaled_area, m_lat, &m_comm);
}

} // namespace scream
