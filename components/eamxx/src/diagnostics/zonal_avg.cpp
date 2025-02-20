#include "diagnostics/zonal_avg.hpp"

#include "share/field/field_utils.hpp"

namespace scream {

ZonalAvgDiag::ZonalAvgDiag(const ekat::Comm &comm,
                           const ekat::ParameterList &params)
    : AtmosphereDiagnostic(comm, params) {
  const auto &field_name = m_params.get<std::string>("field_name");
  m_diag_name            = field_name + "_zonal_avg";
}

void ZonalAvgDiag::set_grids(
    const std::shared_ptr<const GridsManager> grids_manager) {
  const auto &field_name              = m_params.get<std::string>("field_name");
  const auto &grid_name               = m_params.get<std::string>("grid_name");
  const GridsManager::grid_ptr_type grid = grids_manager->get_grid("Physics");

  add_field<Required>(field_name, grid_name);

  // first clone the area unscaled, we will scale it later in initialize_impl
  m_scaled_area = grid->get_geometry_data("area").clone();

  Field m_lat = grid->get_geometry_data("lat"); // TODO: maybe this is a field in?
  // m_lat.get_view<const Real *>()
}

void ZonalAvgDiag::initialize_impl(const RunType /*run_type*/) {
  // TODO: auto everything
  using FieldIdentifier = FieldHeader::identifier_type;
  using FieldLayout = FieldIdentifier::layout_type;
  using ShortFieldTagsNames::COL;
  const Field &field              = get_fields_in().front();
  const FieldIdentifier &field_id = field.get_header().get_identifier();
  const FieldLayout &layout       = field_id.get_layout();

  EKAT_REQUIRE_MSG(layout.rank() >= 1 && layout.rank() <= 3,
                   "Error! Field rank not supported by HorizAvgDiag.\n"
                   " - field name: " +
                       field_id.name() +
                       "\n"
                       " - field layout: " +
                       layout.to_string() + "\n");
  EKAT_REQUIRE_MSG(layout.tags()[0] == COL,
                   "Error! ZonalAvgDiag diagnostic expects a layout starting "
                   "with the 'COL' tag.\n"
                   " - field name  : " +
                       field_id.name() +
                       "\n"
                       " - field layout: " +
                       layout.to_string() + "\n");

  FieldIdentifier diagnostic_field_id(m_diag_name,
                                      layout.clone().strip_dim(COL),
                              field_id.get_units(), field_id.get_grid_name());
  m_diagnostic_output = Field(diagnostic_field_id);
  m_diagnostic_output.allocate_view();

  // scale the area field
  auto total_area = field_sum<Real>(m_scaled_area, &m_comm);
  m_scaled_area.scale(sp(1.0) / total_area);
}

void ZonalAvgDiag::compute_diagnostic_impl() {
  // TODO: auto everything
  using FieldLayout = FieldIdentifier::layout_type;
  const Field &field  = get_fields_in().front();
  const Field &output = m_diagnostic_output;
  const Field &weight = m_scaled_area;
  // Call the horiz_contraction impl that will take care of everything
  //horiz_contraction<Real>(output, field, m_scaled_area, &m_comm);

  using KT          = ekat::KokkosTypes<DefaultDevice>;
  using RangePolicy = Kokkos::RangePolicy<Field::device_t::execution_space>;
  using TeamPolicy  = Kokkos::TeamPolicy<Field::device_t::execution_space>;
  using TeamMember  = typename TeamPolicy::member_type;
  using ESU         = ekat::ExeSpaceUtils<typename KT::ExeSpace>;

  FieldLayout output_layout = output.get_header().get_identifier().get_layout();
  FieldLayout field_layout = field.get_header().get_identifier().get_layout();

  auto weight_view = weight.get_view<const Real *>();

  const int ncols = field_layout.dim(0);

  auto field_view    = field.get_view<const Real ***>();
  auto output_view   = output.get_view<Real **>();
  const int d1 = field_layout.dim(1);
  const int d2 = field_layout.dim(2);
  TeamPolicy policy       = ESU::get_default_team_policy(d1 * d2, ncols);
  Kokkos::parallel_for(
      output.name(), policy, KOKKOS_LAMBDA(const TeamMember &tm) {
        const int idx = tm.league_rank();
        const int j   = idx / d2;
        const int k   = idx % d2;
        Kokkos::parallel_reduce(
            Kokkos::TeamVectorRange(tm, ncols),
            [&](int i, Real &ac) { ac += weight_view(i) * field_view(i, j, k); },
            output_view(j, k));
      });

//  if(m_comm) {
    // TODO: use device-side MPI calls
    // TODO: the dev ptr causes problems; revisit this later
    // TODO: doing cuda-aware MPI allreduce would be ~10% faster
    Kokkos::fence();
    output.sync_to_host();
    m_comm.all_reduce(output.template get_internal_view_data<Real, Host>(),
                     output_layout.size(), MPI_SUM);
    output.sync_to_dev();
//  }
}

}  // namespace scream
