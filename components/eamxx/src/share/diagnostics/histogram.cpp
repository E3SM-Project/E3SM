#include "histogram.hpp"

#include "share/field/field_visit.hpp"

#include <ekat_team_policy_utils.hpp>
#include <ekat_string_utils.hpp>


namespace scream {

Histogram::Histogram(const ekat::Comm &comm, const ekat::ParameterList &params,
                     const std::shared_ptr<const AbstractGrid>& grid)
 : AbstractDiagnostic(comm, params, grid)
{
  m_field_name = m_params.get<std::string>("field_name");
  const std::string bin_config = m_params.get<std::string>("bin_configuration");

  // Extract the bin edges from the configuration and check them. As in numpy and
  // xhistogram, n edges make n-1 bins, and values outside the edges are not counted.
  const std::vector<std::string> bin_strings = ekat::split(bin_config, "_");
  EKAT_REQUIRE_MSG(bin_strings.size()>=2,
      "Error! Histogram needs at least two bin edges.\n"
      " - bin configuration: " + bin_config + "\n");
  m_bin_reals.reserve(bin_strings.size());
  for (const auto& s : bin_strings) {
    auto val = std::stod(s);
    EKAT_REQUIRE_MSG(m_bin_reals.empty() or val > m_bin_reals.back(),
        "Error! Histogram bin values must be monotonically increasing.\n"
        " - bin configuration: " + bin_config + "\n");
    m_bin_reals.push_back(val);
  }

  m_field_in_names.push_back(m_field_name);
}

void Histogram::initialize_impl()
{
  using ShortFieldTagsNames::CMP;

  const auto& field    = m_fields_in.at(m_field_in_names.front());
  const auto& field_id = field.get_header().get_identifier();
  const auto& field_layout = field_id.get_layout();
  EKAT_REQUIRE_MSG(field_layout.rank() >= 1 && field_layout.rank() <= 3,
      "Error! Field rank not supported by Histogram.\n"
      " - field name: " + field_id.name() + "\n"
      " - field layout: " + field_layout.to_string() + "\n");

  // allocate histogram field
  const int num_bins = m_bin_reals.size()-1;
  const auto& bin_config = m_params.get<std::string>("bin_configuration");
  auto diag_name = m_field_name + "_histogram_" + bin_config;
  FieldLayout diagnostic_layout({CMP}, {num_bins}, {"bin"});
  FieldIdentifier diagnostic_id(diag_name, diagnostic_layout,
                                ekat::units::none, field_id.get_grid_name());
  add_output(diagnostic_id);

  // allocate field for bin values
  FieldLayout bin_values_layout({CMP}, {num_bins+1}, {"bin"});
  auto bin_values_id = field_id.clone(diag_name + "_bin_values").reset_layout(bin_values_layout);
  m_bin_values = Field(bin_values_id);
  m_bin_values.allocate_view();

  // copy bin values into field
  auto bin_values_view = m_bin_values.get_view<Real *>();
  auto bin_values_view_host = Kokkos::create_mirror_view(bin_values_view);
  for (auto i=0; i < bin_values_layout.dim(0); i++)
    bin_values_view_host(i) = m_bin_reals[i];
  Kokkos::deep_copy(bin_values_view,bin_values_view_host);
}

namespace {

// The kernel, for a field of rank N (N>=1). The mask is read with a LayoutRight view
// if it allows it (LRMask=true), and with a strided one otherwise (e.g., if it is a broadcast of a smaller mask)
template<int N, bool LRMask>
void compute_histogram_impl (const Field& field, const Field& bin_values, const Field& histogram)
{
  using KT         = ekat::KokkosTypes<DefaultDevice>;
  using TeamPolicy = Kokkos::TeamPolicy<Field::device_t::execution_space>;
  using TeamMember = typename TeamPolicy::member_type;
  using TPF        = ekat::TeamPolicyFactory<typename KT::ExeSpace>;

  const auto& field_layout = field.get_header().get_identifier().get_layout();
  const int num_bins = histogram.get_header().get_identifier().get_layout().dim(0);

  auto bin_values_view = bin_values.get_view<Real *>();
  auto histogram_view  = histogram.get_view<Real *>();
  auto field_view      = field.get_view<Field::data_nd_t<const Real,N>>();

  const bool masked = field.has_valid_mask();
  auto mask_view = masked ? details::get_valid_mask_view<N,LRMask>(field)
                          : details::valid_mask_view_t<N,LRMask>{};

  // Each team loops over all the entries of the field, for one bin
  Kokkos::Array<int,N> dims;
  int size = 1;
  for (int k=0; k<N; ++k) {
    dims[k] = field_layout.dim(k);
    size *= dims[k];
  }

  TeamPolicy team_policy = TPF::get_default_team_policy(num_bins, size);
  Kokkos::parallel_for("compute_histogram_" + field.name(), team_policy,
      KOKKOS_LAMBDA(const TeamMember &tm) {
        const int bin_i = tm.league_rank();
        const Real bin_lower = bin_values_view(bin_i);
        const Real bin_upper = bin_values_view(bin_i+1);
        // Each bin is [lower,upper), except the last, which is [lower,upper]
        const bool last_bin = bin_i==num_bins-1;
        Kokkos::parallel_reduce(Kokkos::TeamVectorRange(tm, size),
            [&](int ind, Real &val) {
              const auto idx = details::unravel<N>(ind, dims);
              const Real f = details::at(field_view,idx);
              if ((not masked or details::at(mask_view,idx)!=0) and
                  bin_lower <= f and
                  (f < bin_upper or (last_bin and f == bin_upper)))
                val += sp(1.0);
            },
            histogram_view(bin_i));
      });
}

} // anonymous namespace

void Histogram::compute_impl()
{
  const auto& field = m_fields_in.at(m_field_name);
  auto field_layout = field.get_header().get_identifier().get_layout();
  auto histogram_layout = get().get_header().get_identifier().get_layout();

  const bool lr_mask = details::valid_mask_is_contiguous(field);
  details::dispatch_on_rank(field_layout.rank(), [&](auto rank_c) {
    constexpr int N = decltype(rank_c)::value;
    if constexpr (N>=1 and N<=3) {
      lr_mask ? compute_histogram_impl<N,true> (field,m_bin_values,get())
              : compute_histogram_impl<N,false>(field,m_bin_values,get());
    } else {
      EKAT_ERROR_MSG("Error! Unsupported field rank for histogram.\n");
    }
  });

  // TODO: use device-side MPI calls
  // TODO: the dev ptr causes problems; revisit this later
  // TODO: doing cuda-aware MPI allreduce would be ~10% faster
  Kokkos::fence();
  get().sync_to_host();
  m_comm.all_reduce(get().template get_internal_view_data<Real, Host>(),
                      histogram_layout.size(), MPI_SUM);
  get().sync_to_dev();
}

} // namespace scream
