#include "zonal_avg.hpp"

#include "share/field/field_utils.hpp"

#include <ekat_math_utils.hpp>
#include <ekat_team_policy_utils.hpp>

namespace scream {

// Utility to compute the contraction of a field along its column dimension.
// This is equivalent to f_out = einsum('i,i...k->...k', weight, f_in).
// The implementation is such that:
// - all Field objects must be allocated
// - the first dimension for field and weight and second dimension for
//     bin_to_cols is for the columns (COL)
// - the first dimension for result and bin_to_cols is for the zonal bins (CMP,"bin")
// - field and result must be the same dimension, up to 3
void compute_zonal_sum(const Field &result, const Field &field, const Field &weight,
                       const Field &bin_to_cols, const ekat::Comm *comm)
{
  auto result_layout = result.get_header().get_identifier().get_layout();
  auto bin_to_cols_layout = bin_to_cols.get_header().get_identifier().get_layout();
  const int num_zonal_bins = bin_to_cols_layout.dim(0);
  const int max_ncols_per_bin = bin_to_cols_layout.dim(1);

  auto weight_view = weight.get_view<const Real *>();
  auto bin_to_cols_view = bin_to_cols.get_view<const Int **>();
  using KT         = ekat::KokkosTypes<DefaultDevice>;
  using TeamPolicy = Kokkos::TeamPolicy<Field::device_t::execution_space>;
  using TeamMember = typename TeamPolicy::member_type;
  using TPF        = ekat::TeamPolicyFactory<typename KT::ExeSpace>;
  using cmask1d_t = Field::view_dev_t<const int*>;
  using cmask2d_t = Field::view_dev_t<const int**>;
  using cmask3d_t = Field::view_dev_t<const int***>;

  bool masked = field.has_valid_mask();
  switch (result_layout.rank()) {
  case 1: {
    auto field_view        = field.get_view<const Real *>();
    auto result_view       = result.get_view<Real *>();
    TeamPolicy team_policy = TPF::get_default_team_policy(num_zonal_bins, max_ncols_per_bin);
    auto mask_view = masked ? field.get_valid_mask().get_view<const int*>() : cmask1d_t{};
    Kokkos::parallel_for(
        "compute_zonal_sum_" + field.name(), team_policy,
        KOKKOS_LAMBDA(const TeamMember &tm) {
          const int bin_i = tm.league_rank();
          Kokkos::parallel_reduce(
              Kokkos::TeamVectorRange(tm, 1, 1+bin_to_cols_view(bin_i,0)),
              [&](int lcol_j, Real &val) {
                const int col_i = bin_to_cols_view(bin_i, lcol_j);
                if (not masked or mask_view(col_i)!=0)
                  val += weight_view(col_i) * field_view(col_i);
              },
              result_view(bin_i));
        });
  } break;
  case 2: {
    const int d1           = result_layout.dim(1);
    auto field_view        = field.get_view<const Real **>();
    auto result_view       = result.get_view<Real **>();
    TeamPolicy team_policy = TPF::get_default_team_policy(num_zonal_bins * d1, max_ncols_per_bin);
    auto mask_view = masked ? field.get_valid_mask().get_view<const int**>() : cmask2d_t{};
    Kokkos::parallel_for(
        "compute_zonal_sum_" + field.name(), team_policy,
        KOKKOS_LAMBDA(const TeamMember &tm) {
          const int idx        = tm.league_rank();
          const int d1_i       = idx / num_zonal_bins;
          const int bin_i      = idx % num_zonal_bins;
          Kokkos::parallel_reduce(
              Kokkos::TeamVectorRange(tm, 1, 1+bin_to_cols_view(bin_i,0)),
              [&](int lcol_j, Real &val) {
                const int col_i = bin_to_cols_view(bin_i, lcol_j);
                if (not masked or mask_view(col_i,d1_i)!=0)
                  val += weight_view(col_i) * field_view(col_i, d1_i);
              },
              result_view(bin_i, d1_i));
        });
  } break;
  case 3: {
    const int d1           = result_layout.dim(1);
    const int d2           = result_layout.dim(2);
    auto field_view        = field.get_view<const Real ***>();
    auto result_view       = result.get_view<Real ***>();
    TeamPolicy team_policy = TPF::get_default_team_policy(num_zonal_bins * d1 * d2, max_ncols_per_bin);
    auto mask_view = masked ? field.get_valid_mask().get_view<const int***>() : cmask3d_t{};
    Kokkos::parallel_for(
        "compute_zonal_sum_" + field.name(), team_policy, KOKKOS_LAMBDA(const TeamMember &tm) {
          const int idx        = tm.league_rank();
          const int d1_i       = idx / (num_zonal_bins * d2);
          const int idx2       = idx % (num_zonal_bins * d2);
          const int d2_i       = idx2 / num_zonal_bins;
          const int bin_i      = idx2 % num_zonal_bins;
          Kokkos::parallel_reduce(
              Kokkos::TeamVectorRange(tm, 1, 1+bin_to_cols_view(bin_i,0)),
              [&](int lcol_j, Real &val) {
                const int col_i = bin_to_cols_view(bin_i, lcol_j);
                if (not masked or mask_view(col_i,d1_i, d2_i)!=0)
                  val += weight_view(col_i) * field_view(col_i, d1_i, d2_i);
              },
              result_view(bin_i, d1_i, d2_i));
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

ZonalAvg::ZonalAvg(const ekat::Comm &comm, const ekat::ParameterList &params,
                   const std::shared_ptr<const AbstractGrid> &grid)
  : AbstractDiagnostic(comm, params, grid)
{
  m_field_name = m_params.get<std::string>("field_name");
  m_num_zonal_bins = std::stoi(params.get<std::string>("number_of_zonal_bins"));

  m_field_in_names.push_back(m_field_name);
  m_lat       = m_grid->get_geometry_data("lat");

  // Area will be used as weight in the zonal sum
  m_area = grid->get_geometry_data("area");
}

void ZonalAvg::initialize_impl()
{
  using namespace ShortFieldTagsNames;
  const auto& field    = m_fields_in.at(m_field_name);
  const auto& field_id = field.get_header().get_identifier();
  const auto& field_layout = field_id.get_layout();

  EKAT_REQUIRE_MSG(field_layout.rank() >= 1 && field_layout.rank() <= 3,
      "Error! Field rank not supported by ZonalAvg.\n"
      " - field name: " + field_id.name() + "\n"
      " - field layout: " + field_layout.to_string() + "\n");
  EKAT_REQUIRE_MSG(field_layout.tags()[0] == COL,
      "Error! ZonalAvg diagnostic expects a layout starting with the 'COL' tag.\n"
      " - field name  : " + field_id.name() + "\n"
      " - field layout: " + field_layout.to_string() + "\n");

  FieldLayout bins_layout({CMP}, {m_num_zonal_bins}, {"bin"});

  // allocate column counter
  FieldIdentifier ncols_per_bin_id("number of columns per bin",
      bins_layout, ekat::units::none, m_grid->name(), DataType::IntType);
  Field ncols_per_bin(ncols_per_bin_id);
  ncols_per_bin.allocate_view();
  ncols_per_bin.deep_copy(0);

  // count how many columns are in each zonal bin
  using KT         = ekat::KokkosTypes<DefaultDevice>;
  using TeamPolicy = Kokkos::TeamPolicy<Field::device_t::execution_space>;
  using TeamMember = typename TeamPolicy::member_type;
  using TPF        = ekat::TeamPolicyFactory<typename KT::ExeSpace>;
  const int ncols      = field_layout.dim(COL);
  const Real lat_delta = sp(180.0) / m_num_zonal_bins;
  auto lat_view    = m_lat.get_view<const Real *>();
  auto ncols_per_bin_view = ncols_per_bin.get_view<Int *>();
  const int num_zonal_bins = m_num_zonal_bins; // for use inside lambdas
  TeamPolicy team_policy = TPF::get_default_team_policy(m_num_zonal_bins, ncols);
  Kokkos::parallel_for("count_columns_per_zonal_bin_" + field.name(),
      team_policy, KOKKOS_LAMBDA(const TeamMember &tm) {
        const int bin_i      = tm.league_rank();
        const Real lat_lower = sp(-90.0) + bin_i * lat_delta;
        const Real lat_upper = (bin_i < num_zonal_bins-1)
          ? lat_lower + lat_delta : sp(90.0 + 0.5*lat_delta);
        Kokkos::parallel_reduce(Kokkos::TeamVectorRange(tm, ncols),
            [&](int col_i, Int &val) {
              if ((lat_lower <= lat_view(col_i)) && (lat_view(col_i) < lat_upper))
                val++;
            },
            ncols_per_bin_view(bin_i));
      });

  // determine maximum number of columns per bin & allocate bin to column map
  using RangePolicy     = Kokkos::RangePolicy<Field::device_t::execution_space>;
  Int max_ncols_per_bin = 0;
  Kokkos::parallel_reduce(RangePolicy(0, m_num_zonal_bins),
      KOKKOS_LAMBDA(int bin_i, Int &val) {
        val = ncols_per_bin_view(bin_i) > val ? ncols_per_bin_view(bin_i) : val;
      },
      Kokkos::Max<Int>(max_ncols_per_bin));
  FieldLayout bin_to_cols_layout = bins_layout.clone().append_dim(COL, 1+max_ncols_per_bin);
  FieldIdentifier bin_to_cols_id("columns in each zonal bin",
    bin_to_cols_layout, ekat::units::none, m_grid->name() , DataType::IntType);
  m_bin_to_cols = Field(bin_to_cols_id);
  m_bin_to_cols.allocate_view();

  // compute bin to column map, where the (i,j)-th entry is such that
  //   - for j=0, the entry is the number of columns in the i-th zonal bin
  //   - for j>0, the entry is a column index in the i-th zonal bin
  auto bin_to_cols_view = m_bin_to_cols.get_view<Int **>();
  Kokkos::parallel_for("assign_columns_to_zonal_bins_" + field.name(),
      RangePolicy(0, m_num_zonal_bins), KOKKOS_LAMBDA(int bin_i) {
        const Real lat_lower = sp(-90.0) + bin_i * lat_delta;
        const Real lat_upper = (bin_i < num_zonal_bins-1)
          ? lat_lower + lat_delta : sp(90.0 + 0.5*lat_delta);
        bin_to_cols_view(bin_i, 0) = 0;
        for (int col_i=0; col_i < ncols; col_i++)
        {
          if ((lat_lower <= lat_view(col_i)) && (lat_view(col_i) < lat_upper))
          {
            bin_to_cols_view(bin_i, 0) += 1;
            bin_to_cols_view(bin_i, bin_to_cols_view(bin_i,0)) = col_i;
          }
        }
      });

  // Create the diagnostic
  auto diag_name = m_field_name + "_zonal_avg_" + std::to_string(m_num_zonal_bins) + "_bins";
  auto diagnostic_layout = field_layout.clone().strip_dim(COL).prepend_dim(CMP, m_num_zonal_bins, "bin");
  auto diagnostic_id = field_id.clone(diag_name).reset_layout(diagnostic_layout);
  m_diagnostic_output = Field(diagnostic_id,true);

  if (not field.has_valid_mask()) {
    m_ones = m_area.clone("ones");
    m_ones.deep_copy(1);
    // It's not nondimensional, but we're just using it as a scaling factor field
    FieldIdentifier zonal_area_fid ("zonal area",bins_layout,ekat::units::none,m_grid->name());
    m_zonal_area = Field(zonal_area_fid,true);
    compute_zonal_sum(m_zonal_area,m_ones,m_area,m_bin_to_cols,&m_comm);

    // Create a copy of area which is scaled by 1 / zonal area
    m_scaled_area = m_area.clone("scaled_area");
    auto zonal_area_view  = m_zonal_area.get_view<const Real *>();
    auto scaled_area_view = m_scaled_area.get_view<Real *>();
    Kokkos::parallel_for("scale_area_by_zonal_area_" + field.name(),
        team_policy, KOKKOS_LAMBDA(const TeamMember &tm) {
          const int bin_i = tm.league_rank();
          Kokkos::parallel_for(
            Kokkos::TeamVectorRange(tm, 1, 1+bin_to_cols_view(bin_i,0)),
              [&](int lcol_j) {
                const int col_i = bin_to_cols_view(bin_i, lcol_j);
                scaled_area_view(col_i) /= zonal_area_view(bin_i);
              });
        });
  } else {
    m_ones = field.clone("ones");
    m_ones.deep_copy(1);
    m_ones.set_valid_mask(field.get_valid_mask());
    m_zonal_area = m_diagnostic_output.clone("zonal_area");
    m_diagnostic_output.create_valid_mask();
    m_diagnostic_output.get_header().set_may_be_filled(true);
  }
}

void ZonalAvg::compute_impl()
{
  const auto &field = m_fields_in.begin()->second;

  if (field.has_valid_mask()) {
    // Compute masked field zonal sum, masked area zonal sum, and divide
    compute_zonal_sum(m_diagnostic_output, field, m_area, m_bin_to_cols, &m_comm);
    compute_zonal_sum(m_zonal_area,m_ones,m_area,m_bin_to_cols,&m_comm);

    auto& dmask = m_diagnostic_output.get_valid_mask();
    compute_mask(m_zonal_area,0,Comparison::NE,dmask);
    m_diagnostic_output.scale_inv(m_zonal_area,dmask);

    // TODO: remove when IO stops relying on mask=0 entries being already set to FillValue
    m_diagnostic_output.deep_copy(constants::fill_value<Real>,dmask,true);
  } else {
    // Compute field zonal sum using scaled area as weight
    compute_zonal_sum(m_diagnostic_output, field, m_scaled_area, m_bin_to_cols, &m_comm);
  }
}

} // namespace scream
