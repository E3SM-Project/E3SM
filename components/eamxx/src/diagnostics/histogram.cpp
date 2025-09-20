#include "diagnostics/histogram.hpp"

#include <ekat_team_policy_utils.hpp>
#include <ekat_string_utils.hpp>


namespace scream {

HistogramDiag::HistogramDiag(const ekat::Comm &comm, const ekat::ParameterList &params)
    : AtmosphereDiagnostic(comm, params) {
  const auto &field_name       = m_params.get<std::string>("field_name");
  const std::string bin_config = m_params.get<std::string>("bin_configuration");
  m_diag_name                  = field_name + "_histogram_" + bin_config;

  // extract bin values from configuration, append end values, and check
  const std::vector<std::string> bin_strings = ekat::split(bin_config, "_");
  m_bin_reals.resize(bin_strings.size()+2);
  m_bin_reals[0] = std::numeric_limits<Real>::lowest();
  for (long unsigned int i=1; i < m_bin_reals.size()-1; i++)
  {
    m_bin_reals[i] = std::stod(bin_strings[i-1]);
    EKAT_REQUIRE_MSG(m_bin_reals[i] > m_bin_reals[i-1],
                     "Error! HistogramDiag bin values must be monotonically "
                     "increasing.\n"
                     " - bin configuration: " + bin_config + "\n");
  }
  m_bin_reals.back() = std::numeric_limits<Real>::max();
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
  const FieldLayout &field_layout = field_id.get_layout();
  EKAT_REQUIRE_MSG(field_layout.rank() >= 1 && field_layout.rank() <= 3,
                   "Error! Field rank not supported by HistogramDiag.\n"
                   " - field name: " +
                       field_id.name() +
                       "\n"
                       " - field layout: " +
                       field_layout.to_string() + "\n");

  // allocate histogram field
  const int num_bins = m_bin_reals.size()-1;
  FieldLayout diagnostic_layout({CMP}, {num_bins}, {"bin"});
  FieldIdentifier diagnostic_id(m_diag_name, diagnostic_layout,
                                FieldIdentifier::Units::nondimensional(),
                                field_id.get_grid_name());
  m_diagnostic_output = Field(diagnostic_id);
  m_diagnostic_output.allocate_view();

  // allocate field for bin values
  FieldLayout bin_values_layout({CMP}, {num_bins+1}, {"bin"});
  FieldIdentifier bin_values_id(m_diag_name + "_bin_values", bin_values_layout,
                                field_id.get_units(), field_id.get_grid_name());
  m_bin_values = Field(bin_values_id);
  m_bin_values.allocate_view();

  // copy bin values into field
  auto bin_values_view = m_bin_values.get_view<Real *>();
  auto bin_values_view_host = Kokkos::create_mirror_view(bin_values_view);
  for (long unsigned int i=0; i < bin_values_layout.dim(0); i++)
    bin_values_view_host(i) = m_bin_reals[i];
  Kokkos::deep_copy(bin_values_view,bin_values_view_host);
}

void HistogramDiag::compute_diagnostic_impl() {
  const auto &field = get_fields_in().front();
  auto field_layout = field.get_header().get_identifier().get_layout();
  auto histogram_layout = m_diagnostic_output.get_header().get_identifier().get_layout();
  const int num_bins = histogram_layout.dim(0);

  auto bin_values_view = m_bin_values.get_view<Real *>();
  auto histogram_view = m_diagnostic_output.get_view<Real *>();
  using KT         = ekat::KokkosTypes<DefaultDevice>;
  using TeamPolicy = Kokkos::TeamPolicy<Field::device_t::execution_space>;
  using TeamMember = typename TeamPolicy::member_type;
  using TPF        = ekat::TeamPolicyFactory<typename KT::ExeSpace>;
  switch (field_layout.rank())
  {
    case 1: {
      const int d1 = field_layout.dim(0);
      auto field_view = field.get_view<const Real *>();
      TeamPolicy team_policy = TPF::get_default_team_policy(num_bins, d1);
      Kokkos::parallel_for("compute_histogram_" + field.name(), team_policy,
          KOKKOS_LAMBDA(const TeamMember &tm) {
            const int bin_i = tm.league_rank();
            const Real bin_lower = bin_values_view(bin_i);
            const Real bin_upper = bin_values_view(bin_i+1);
            Kokkos::parallel_reduce(Kokkos::TeamVectorRange(tm, d1),
                [&](int i, Real &val) {
                  if ((bin_lower <= field_view(i)) && (field_view(i) < bin_upper))
                    val += sp(1.0);
                },
                histogram_view(bin_i));
          });
    } break;
    case 2: {
      const int d1 = field_layout.dim(0);
      const int d2 = field_layout.dim(1);
      auto field_view = field.get_view<const Real **>();
      TeamPolicy team_policy = TPF::get_default_team_policy(num_bins, d1*d2);
      Kokkos::parallel_for("compute_histogram_" + field.name(), team_policy,
          KOKKOS_LAMBDA(const TeamMember &tm) {
            const int bin_i = tm.league_rank();
            const Real bin_lower = bin_values_view(bin_i);
            const Real bin_upper = bin_values_view(bin_i+1);
            Kokkos::parallel_reduce(Kokkos::TeamVectorRange(tm, d1*d2),
                [&](int ind, Real &val) {
                  const int i1 = ind / d2;
                  const int i2 = ind % d2;
                  if ((bin_lower <= field_view(i1,i2)) && (field_view(i1,i2) < bin_upper))
                    val += sp(1.0);
                },
                histogram_view(bin_i));
          });
    } break;
    case 3: {
      const int d1 = field_layout.dim(0);
      const int d2 = field_layout.dim(1);
      const int d3 = field_layout.dim(2);
      auto field_view = field.get_view<const Real ***>();
      TeamPolicy team_policy = TPF::get_default_team_policy(num_bins, d1*d2*d3);
      Kokkos::parallel_for("compute_histogram_" + field.name(), team_policy,
          KOKKOS_LAMBDA(const TeamMember &tm) {
            const int bin_i = tm.league_rank();
            const Real bin_lower = bin_values_view(bin_i);
            const Real bin_upper = bin_values_view(bin_i+1);
            Kokkos::parallel_reduce(Kokkos::TeamVectorRange(tm, d1*d2*d3),
                [&](int ind, Real &val) {
                  const int i1 = ind / (d2*d3);
                  const int ind2 = ind % (d2*d3);
                  const int i2 = ind2 / d3;
                  const int i3 = ind2 % d3;
                  if ((bin_lower <= field_view(i1,i2,i3)) && (field_view(i1,i2,i3) < bin_upper))
                    val += sp(1.0);
                },
                histogram_view(bin_i));
          });
  } break;

  default:
    EKAT_ERROR_MSG("Error! Unsupported field rank for histogram.\n");
  }

  // TODO: use device-side MPI calls
  // TODO: the dev ptr causes problems; revisit this later
  // TODO: doing cuda-aware MPI allreduce would be ~10% faster
  Kokkos::fence();
  m_diagnostic_output.sync_to_host();
  m_comm.all_reduce(m_diagnostic_output.template get_internal_view_data<Real, Host>(),
                      histogram_layout.size(), MPI_SUM);
  m_diagnostic_output.sync_to_dev();
}

} // namespace scream
