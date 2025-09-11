#include "diagnostics/conditional_sampling.hpp"
#include "share/util/eamxx_universal_constants.hpp"
#include <ekat_team_policy_utils.hpp>
#include <string>

namespace scream {

// Utility function to apply conditional sampling logic
KOKKOS_INLINE_FUNCTION
bool evaluate_condition(const Real &condition_val, const int &op_code, const Real &comparison_val) {
  // op_code: 0=eq, 1=ne, 2=gt, 3=ge, 4=lt, 5=le
  switch (op_code) {
    case 0: return condition_val == comparison_val;  // eq or ==
    case 1: return condition_val != comparison_val;  // ne or !=
    case 2: return condition_val > comparison_val;   // gt or >
    case 3: return condition_val >= comparison_val;  // ge or >=
    case 4: return condition_val < comparison_val;   // lt or <
    case 5: return condition_val <= comparison_val;  // le or <=
    default: return false;
  }
}

// Utility function to convert operator string to code
int get_operator_code(const std::string& op) {
  if (op == "eq" || op == "==") return 0;
  if (op == "ne" || op == "!=") return 1;
  if (op == "gt" || op == ">")  return 2;
  if (op == "ge" || op == ">=") return 3;
  if (op == "lt" || op == "<")  return 4;
  if (op == "le" || op == "<=") return 5;
  return -1; // Invalid operator
}

// Utility function to apply conditional sampling for 1D fields (either ncols or nlevs)
void apply_conditional_sampling_1d(
    const Field &output_field, const Field &input_field, const Field &condition_field,
    const std::string &condition_op, const Real &condition_val,
    const Real &fill_value = constants::fill_value<Real>) {

  // if fill_value is 0, we are counting
  const auto is_counting = (fill_value == 0);
  const auto output_v    = output_field.get_view<Real *>();
  const auto mask_v = !is_counting ? output_field.get_header().get_extra_data<Field>("mask_data").get_view<Real *>() : output_v;
  const auto input_v     = input_field.get_view<const Real *>();
  const auto condition_v = condition_field.get_view<const Real *>();

  // Try to get input and condition masks, if present
  bool has_input_mask = input_field.get_header().has_extra_data("mask_data");
  bool has_condition_mask = condition_field.get_header().has_extra_data("mask_data");
  auto input_mask_v = has_input_mask ? input_field.get_header().get_extra_data<Field>("mask_data").get_view<const Real *>() : input_v;
  auto condition_mask_v = has_condition_mask ? condition_field.get_header().get_extra_data<Field>("mask_data").get_view<const Real *>() : condition_v;

  const int n_elements = output_field.get_header().get_identifier().get_layout().dims()[0];

  // Convert operator string to integer code for device use
  const int op_code = get_operator_code(condition_op);

  Kokkos::parallel_for(
      "ConditionalSampling1D", Kokkos::RangePolicy<>(0, n_elements), KOKKOS_LAMBDA(const int &idx) {
        bool input_masked = has_input_mask && (input_mask_v(idx) == 0);
        bool condition_masked = has_condition_mask && (condition_mask_v(idx) == 0);
        if (input_masked || condition_masked) {
          output_v(idx) = fill_value;
          if (!is_counting) mask_v(idx) = 0;
        } else if (evaluate_condition(condition_v(idx), op_code, condition_val)) {
          output_v(idx) = input_v(idx);
          if (!is_counting) mask_v(idx) = 1;
        } else {
          output_v(idx) = fill_value;
          if (!is_counting) mask_v(idx) = 0;
        }
      });
}

// Utility function to apply conditional sampling for 2D fields (ncols x nlevs)
void apply_conditional_sampling_2d(
    const Field &output_field, const Field &input_field, const Field &condition_field,
    const std::string &condition_op, const Real &condition_val,
    const Real &fill_value = constants::fill_value<Real>) {

  // if fill_value is 0, we are counting
  const auto is_counting = (fill_value == 0);

  const auto output_v    = output_field.get_view<Real **>();
  const auto mask_v = !is_counting ? output_field.get_header().get_extra_data<Field>("mask_data").get_view<Real **>() : output_v;
  const auto input_v     = input_field.get_view<const Real **>();
  const auto condition_v = condition_field.get_view<const Real **>();

  // Try to get input and condition masks, if present
  bool has_input_mask = input_field.get_header().has_extra_data("mask_data");
  bool has_condition_mask = condition_field.get_header().has_extra_data("mask_data");
  auto input_mask_v = has_input_mask ? input_field.get_header().get_extra_data<Field>("mask_data").get_view<const Real **>() : input_v;
  auto condition_mask_v = has_condition_mask ? condition_field.get_header().get_extra_data<Field>("mask_data").get_view<const Real **>() : condition_v;

  const int ncols = output_field.get_header().get_identifier().get_layout().dims()[0];
  const int nlevs = output_field.get_header().get_identifier().get_layout().dims()[1];

  // Convert operator string to integer code for device use
  const int op_code = get_operator_code(condition_op);

  Kokkos::parallel_for(
      "ConditionalSampling2D", Kokkos::RangePolicy<>(0, ncols * nlevs),
      KOKKOS_LAMBDA(const int &idx) {
        const int icol = idx / nlevs;
        const int ilev = idx % nlevs;
        bool input_masked = has_input_mask && (input_mask_v(icol, ilev) == 0);
        bool condition_masked = has_condition_mask && (condition_mask_v(icol, ilev) == 0);
        if (input_masked || condition_masked) {
          output_v(icol, ilev) = fill_value;
          if (!is_counting) mask_v(icol, ilev) = 0;
        } else if (evaluate_condition(condition_v(icol, ilev), op_code, condition_val)) {
          output_v(icol, ilev) = input_v(icol, ilev);
          if (!is_counting) mask_v(icol, ilev) = 1;
        } else {
          output_v(icol, ilev) = fill_value;
          if (!is_counting) mask_v(icol, ilev) = 0;
        }
      });
}

// Utility function to apply conditional sampling for 1D fields against level indices
void apply_conditional_sampling_1d_lev(
    const Field &output_field, const Field &input_field,
    const std::string &condition_op, const Real &condition_val,
    const Real &fill_value = constants::fill_value<Real>) {

  // if fill_value is 0, we are counting
  const auto is_counting = (fill_value == 0);

  const auto output_v = output_field.get_view<Real *>();
  const auto mask_v = !is_counting ? output_field.get_header().get_extra_data<Field>("mask_data").get_view<Real *>() : output_v;
  const auto input_v  = input_field.get_view<const Real *>();

  // Try to get input mask, if present
  bool has_input_mask = input_field.get_header().has_extra_data("mask_data");
  auto input_mask_v = has_input_mask ? input_field.get_header().get_extra_data<Field>("mask_data").get_view<const Real *>() : input_v;

  const int n_elements = output_field.get_header().get_identifier().get_layout().dims()[0];

  // Convert operator string to integer code for device use
  const int op_code = get_operator_code(condition_op);

  Kokkos::parallel_for(
      "ConditionalSampling1D_Lev", Kokkos::RangePolicy<>(0, n_elements), 
      KOKKOS_LAMBDA(const int &idx) {
        // For 1D case, the level index is just the element index
        const Real level_idx = static_cast<Real>(idx);
        bool input_masked = has_input_mask && (input_mask_v(idx) == 0);
        if (input_masked) {
          output_v(idx) = fill_value;
          if (!is_counting) mask_v(idx) = 0;
        } else if (evaluate_condition(level_idx, op_code, condition_val)) {
          output_v(idx) = input_v(idx);
          if (!is_counting) mask_v(idx) = 1;
        } else {
          output_v(idx) = fill_value;
          if (!is_counting) mask_v(idx) = 0;
        }
      });
}

// Utility function to apply conditional sampling for 2D fields against level indices
void apply_conditional_sampling_2d_lev(
    const Field &output_field, const Field &input_field,
    const std::string &condition_op, const Real &condition_val,
    const Real &fill_value = constants::fill_value<Real>) {

  // if fill_value is 0, we are counting
  const auto is_counting = (fill_value == 0);

  const auto output_v = output_field.get_view<Real **>();
  const auto mask_v = !is_counting ? output_field.get_header().get_extra_data<Field>("mask_data").get_view<Real **>() : output_v;
  const auto input_v  = input_field.get_view<const Real **>();

  // Try to get input mask, if present
  bool has_input_mask = input_field.get_header().has_extra_data("mask_data");
  auto input_mask_v = has_input_mask ? input_field.get_header().get_extra_data<Field>("mask_data").get_view<const Real **>() : input_v;

  const int ncols = output_field.get_header().get_identifier().get_layout().dims()[0];
  const int nlevs = output_field.get_header().get_identifier().get_layout().dims()[1];

  // Convert operator string to integer code for device use
  const int op_code = get_operator_code(condition_op);

  Kokkos::parallel_for(
      "ConditionalSampling2D_Lev", Kokkos::RangePolicy<>(0, ncols * nlevs),
      KOKKOS_LAMBDA(const int &idx) {
        const int icol = idx / nlevs;
        const int ilev = idx % nlevs;
        // The level index is the second dimension index (ilev)
        const Real level_idx = static_cast<Real>(ilev);
        bool input_masked = has_input_mask && (input_mask_v(icol, ilev) == 0);
        if (input_masked) {
          output_v(icol, ilev) = fill_value;
          if (!is_counting) mask_v(icol, ilev) = 0;
        } else if (evaluate_condition(level_idx, op_code, condition_val)) {
          output_v(icol, ilev) = input_v(icol, ilev);
          if (!is_counting) mask_v(icol, ilev) = 1;
        } else {
          output_v(icol, ilev) = fill_value;
          if (!is_counting) mask_v(icol, ilev) = 0;
        }
      });
}

ConditionalSampling::ConditionalSampling(const ekat::Comm &comm, const ekat::ParameterList &params)
    : AtmosphereDiagnostic(comm, params) {

  m_input_f      = m_params.get<std::string>("input_field");
  m_condition_f  = m_params.get<std::string>("condition_field");
  m_condition_op = m_params.get<std::string>("condition_operator");

  const auto str_condition_v = m_params.get<std::string>("condition_value");
  // TODO: relying on std::stod to throw if bad val is given
  m_condition_v = static_cast<Real>(std::stod(str_condition_v));

  m_diag_name =
      m_input_f + "_where_" + m_condition_f + "_" + m_condition_op + "_" + str_condition_v;
}

void ConditionalSampling::set_grids(const std::shared_ptr<const GridsManager> grids_manager) {
  const auto &gn = m_params.get<std::string>("grid_name");
  const auto g   = grids_manager->get_grid("physics");

  // Special case: if the input field is "count", we don't need to add it
  if (m_input_f != "count") {
    add_field<Required>(m_input_f, gn);
  }

  // Special case: if condition field is "lev", we don't add it as a required field
  // since it's geometric information from the grid, not an actual field
  if (m_condition_f != "lev") {
    add_field<Required>(m_condition_f, gn);
  } else {
    // Store grid info for potential use in count operations
    m_nlevs = g->get_num_vertical_levels();
    m_gn = gn;
  }
}

void ConditionalSampling::initialize_impl(const RunType /*run_type*/) {

  if (m_input_f != "count") {
    auto ifid = get_field_in(m_input_f).get_header().get_identifier();
    FieldIdentifier d_fid(m_diag_name, ifid.get_layout().clone(), ifid.get_units(),
                          ifid.get_grid_name());
    m_diagnostic_output = Field(d_fid);
    m_diagnostic_output.allocate_view();
  } else {
    if (m_condition_f != "lev") {
      auto ifid = get_field_in(m_condition_f).get_header().get_identifier();
      FieldIdentifier d_fid(m_diag_name, ifid.get_layout().clone(), ifid.get_units(),
                          ifid.get_grid_name());
      m_diagnostic_output = Field(d_fid);
      m_diagnostic_output.allocate_view();
    } else {
      using namespace ShortFieldTagsNames;
      using namespace ekat::units;
      const auto nondim = Units::nondimensional();
      FieldIdentifier d_fid(m_diag_name, {{LEV},{m_nlevs}}, nondim, m_gn);
      m_diagnostic_output = Field(d_fid);
      m_diagnostic_output.allocate_view();
    }
  }

  auto ifid = m_diagnostic_output.get_header().get_identifier();
  FieldIdentifier mask_fid(m_diag_name + "_mask", ifid.get_layout().clone(), ifid.get_units(), ifid.get_grid_name());
  Field diag_mask(mask_fid);
  diag_mask.allocate_view();

  const auto var_fill_value = constants::fill_value<Real>;
  m_mask_val = m_params.get<double>("mask_value", var_fill_value);
  if (m_input_f != "count") {
    m_diagnostic_output.get_header().set_extra_data("mask_data", diag_mask);
    m_diagnostic_output.get_header().set_extra_data("mask_value", m_mask_val);
  }
  // Special case: if the input field is "count", let's create a field of 1s
  if (m_input_f == "count") {
    ones = m_diagnostic_output.clone("count_ones");
    ones.deep_copy(1.0);
  }
  
  // Special case: if condition field is "lev", we don't need to check layout compatibility
  // since "lev" is geometric information, not an actual field
  if (m_condition_f == "lev") {
    using namespace ShortFieldTagsNames;
    EKAT_REQUIRE_MSG(ifid.get_layout().tags().back() == LEV,
                   "Error! ConditionalSampling with \"lev\" condition field must have level in layout.\n"
                   " - field name: " + ifid.name() + "\n"
                   " - field layout: " + ifid.get_layout().to_string() + "\n");
  } else {
    const auto cfid = get_field_in(m_condition_f).get_header().get_identifier();
    
    // check that m_input_f and m_condition_f have the same layout
    EKAT_REQUIRE_MSG(ifid.get_layout() == cfid.get_layout(),
                     "Error! ConditionalSampling only supports comparing fields of the same layout.\n"
                     " - input field has layout of " + ifid.get_layout().to_string() + "\n" +
                     " - condition field has layout of " + cfid.get_layout().to_string() + "\n");
  }
}

void ConditionalSampling::compute_diagnostic_impl() {
  Field f;
  if (m_input_f == "count") {
    // Special case: if the input field is "count", we use the diagnostic output as the input
    f = ones;
  } else {
    f = get_field_in(m_input_f);
  }
  const auto &d = m_diagnostic_output;

  // Validate operator
  const int op_code = get_operator_code(m_condition_op);
  EKAT_REQUIRE_MSG(op_code >= 0,
                   "Error! Invalid condition operator: '" + m_condition_op + "'\n"
                   "Valid operators are: eq, ==, ne, !=, gt, >, ge, >=, lt, <, le, <=\n");

  // Get the fill value from constants
  const Real fill_value = (m_input_f == "count") ? 0.0 : m_mask_val;

  // Determine field layout and apply appropriate conditional sampling
  const auto &layout = f.get_header().get_identifier().get_layout();
  const int rank     = layout.rank();

  if (rank > 2) {
      // no support for now, contact devs
      EKAT_ERROR_MSG("Error! ConditionalSampling only supports 1D or 2D field layouts.\n"
                     " - field layout: " + layout.to_string() + "\n"
                     " - field rank: " + std::to_string(rank) + "\n");
  }
  if (m_condition_f == "lev") {
    // Level-based conditional sampling
    if (rank == 1) {
      // 1D field: level index is just the element index
      apply_conditional_sampling_1d_lev(d, f, m_condition_op, m_condition_v, fill_value);
    } else if (rank == 2) {
      // 2D field: level index is the second dimension (ilev)
      apply_conditional_sampling_2d_lev(d, f, m_condition_op, m_condition_v, fill_value);
    }
  } else {
    // Field-based conditional sampling
    const auto &c = get_field_in(m_condition_f);
    if (rank == 1) {
      // 1D field: (ncols) or (nlevs)
      apply_conditional_sampling_1d(d, f, c, m_condition_op, m_condition_v, fill_value);
    } else if (rank == 2) {
      // 2D field: (ncols, nlevs)
      apply_conditional_sampling_2d(d, f, c, m_condition_op, m_condition_v, fill_value);
    }
  }
}

} // namespace scream
