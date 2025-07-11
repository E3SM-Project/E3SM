#include "diagnostics/histogram.hpp"

namespace scream {


HistogramDiag::HistogramDiag(const ekat::Comm &comm, const ekat::ParameterList &params)
    : AtmosphereDiagnostic(comm, params) {
  const auto &field_name = m_params.get<std::string>("field_name");
  m_bin_configuration    = params.get<std::string>("bin_configuration");
  m_diag_name            = field_name + "_histogram_" + m_bin_configuration;
}

void HistogramDiag::set_grids(const std::shared_ptr<const GridsManager> grids_manager) {
  const auto &field_name = m_params.get<std::string>("field_name");
  const auto &grid_name  = m_params.get<std::string>("grid_name");
  add_field<Required>(field_name, grid_name);
}

void HistogramDiag::initialize_impl(const RunType /*run_type*/) {
  using FieldIdentifier = FieldHeader::identifier_type;
  using FieldLayout     = FieldIdentifier::layout_type;
  using ShortFieldTagsNames::CMP;
  const Field &field              = get_fields_in().front();
  const FieldIdentifier &field_id = field.get_header().get_identifier();

  // extract bin values from configuration
  std::istringstream stream(m_bin_configuration);
  std::string token;
  std::vector<Real> bin_values;
  while (std::getline(stream, token, '_')) {
    bin_values.push_back(std::stod(token));
  }
  int num_bins = bin_values.size()-1;

  // allocate histogram field
  FieldLayout diagnostic_layout({CMP}, {num_bins}, {"bin"});
  FieldIdentifier diagnostic_id(m_diag_name, diagnostic_layout, field_id.get_units(),
                                field_id.get_grid_name(), DataType::IntType);
  m_diagnostic_output = Field(diagnostic_id);
  m_diagnostic_output.allocate_view();

  // allocate field for bin values
  FieldLayout bin_values_layout({CMP}, {num_bins+1}, {"bin"});
  FieldIdentifier bin_values_id(m_diag_name + "_bin_values", bin_values_layout,
                                field_id.get_units(), field_id.get_grid_name());
  m_bin_values = Field(bin_values_id);
  m_bin_values.allocate_view();
  auto bin_values_view = m_bin_values.get_view<Real *>();
  using RangePolicy    = Kokkos::RangePolicy<Field::device_t::execution_space>;
  Kokkos::parallel_for("store_histogram_bin_values_" + field.name(),
      RangePolicy(0, bin_values_layout.dim(0)), [&] (int i) {
        bin_values_view(i) = bin_values[i];
      });
}

void HistogramDiag::compute_diagnostic_impl() {
  const auto &field = get_fields_in().front();
  auto field_layout = field.get_header().get_identifier().get_layout();
  auto histogram_layout = m_diagnostic_output.get_header().get_identifier().get_layout();
  const int num_bins = histogram_layout.dim(0);
  const int field_size = field_layout.size();

  auto bin_values_view = m_bin_values.get_view<Real *>();
  auto histogram_view = m_diagnostic_output.get_view<Int *>();
  auto field_view = field.get_view<const Real *>();
  using KT         = ekat::KokkosTypes<DefaultDevice>;
  using TeamPolicy = Kokkos::TeamPolicy<Field::device_t::execution_space>;
  using TeamMember = typename TeamPolicy::member_type;
  using ESU        = ekat::ExeSpaceUtils<typename KT::ExeSpace>;
  TeamPolicy team_policy = ESU::get_default_team_policy(num_bins, field_size);
  Kokkos::parallel_for("compute_histogram_" + field.name(), team_policy,
      KOKKOS_LAMBDA(const TeamMember &tm) {
        const int bin_i = tm.league_rank();
        const Real bin_lower = bin_values_view(bin_i);
        const Real bin_upper = bin_values_view(bin_i+1);
        Kokkos::parallel_reduce(Kokkos::TeamVectorRange(tm, field_size),
            [&](int i, Int &val) {
              if ((bin_lower <= field_view(i)) && (field_view(i) < bin_upper))
                val++;
            },
            histogram_view(bin_i));
      });

  // TODO: use device-side MPI calls
  // TODO: the dev ptr causes problems; revisit this later
  // TODO: doing cuda-aware MPI allreduce would be ~10% faster
  Kokkos::fence();
  m_diagnostic_output.sync_to_host();
  m_comm.all_reduce(m_diagnostic_output.template get_internal_view_data<Int, Host>(),
                      histogram_layout.size(), MPI_SUM);
  m_diagnostic_output.sync_to_dev();
}

} // namespace scream
