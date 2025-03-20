#include "diagnostics/zonal_avg.hpp"

#include "share/field/field_utils.hpp"

namespace scream {

ZonalAvgDiag::ZonalAvgDiag(const ekat::Comm &comm,
                           const ekat::ParameterList &params)
    : AtmosphereDiagnostic(comm, params) {
  const auto &field_name = m_params.get<std::string>("field_name");
  m_diag_name = field_name + "_zonal_avg";
}

void ZonalAvgDiag::set_grids(
    const std::shared_ptr<const GridsManager> grids_manager) {
  const auto &field_name = m_params.get<std::string>("field_name");
  const auto &grid_name = m_params.get<std::string>("grid_name");
  const GridsManager::grid_ptr_type grid = grids_manager->get_grid("Physics");

  add_field<Required>(field_name, grid_name);

  // first clone the area unscaled, we will scale it later in initialize_impl
  m_scaled_area = grid->get_geometry_data("area").clone();

  m_lat = grid->get_geometry_data("lat");
}

void ZonalAvgDiag::initialize_impl(const RunType /*run_type*/) {
  // TODO: auto everything
  using FieldIdentifier = FieldHeader::identifier_type;
  using FieldLayout = FieldIdentifier::layout_type;
  using ShortFieldTagsNames::COL;
  const Field &field = get_fields_in().front();
  const FieldIdentifier &field_id = field.get_header().get_identifier();
  const FieldLayout &layout = field_id.get_layout();

  // TODO: support 3D and >3D fields
  EKAT_REQUIRE_MSG(layout.rank() >= 1 && layout.rank() <= 1,
                   "Error! Field rank not supported by ZonalAvgDiag.\n"
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

  // TODO: make this dynamic
  FieldLayout diagnostic_layout({COL}, {2}, {"lat"});
  FieldIdentifier diagnostic_id(m_diag_name, diagnostic_layout,
    field_id.get_units(), field_id.get_grid_name());
  m_diagnostic_output = Field(diagnostic_id);
  m_diagnostic_output.allocate_view();

  // scale the area field
  auto total_area = field_sum<Real>(m_scaled_area, &m_comm);
  m_scaled_area.scale(sp(1.0) / total_area);
}

void ZonalAvgDiag::compute_diagnostic_impl() {
  // TODO: auto everything
  using FieldLayout = FieldIdentifier::layout_type;
  const Field &field = get_fields_in().front();
  const Field &output = m_diagnostic_output;
  const Field &lat = m_lat;

  //using KT = ekat::KokkosTypes<DefaultDevice>;
  using RangePolicy = Kokkos::RangePolicy<Field::device_t::execution_space>;
  //using TeamPolicy = Kokkos::TeamPolicy<Field::device_t::execution_space>;
  //using TeamMember = typename TeamPolicy::member_type;
  //using ESU = ekat::ExeSpaceUtils<typename KT::ExeSpace>;

  FieldLayout output_layout = output.get_header().get_identifier().get_layout();
  FieldLayout field_layout = field.get_header().get_identifier().get_layout();

  const int ncols = field_layout.dim(0);
  auto scaled_area_view = m_scaled_area.get_view<const Real *>();
  auto lat_view = m_lat.get_view<const Real *>();

  switch (field_layout.rank()) {
    case 1: {
      auto field_view    = field.get_view<const Real *>();
      auto output_view   = output.get_view<Real *>();
      Kokkos::parallel_reduce(output.name(), RangePolicy(0, ncols),
        KOKKOS_LAMBDA(const int i, Real& val)
        {
          if (lat_view(i) < 0.0)
            val += scaled_area_view(i) * field_view(i);
        }, output_view(0));
      Kokkos::parallel_reduce(output.name(), RangePolicy(0, ncols),
        KOKKOS_LAMBDA(const int i, Real& val)
        {
          if (lat_view(i) > 0.0)
            val += scaled_area_view(i) * field_view(i);
        }, output_view(1));
    } break;
  }

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
