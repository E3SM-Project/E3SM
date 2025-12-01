#ifndef SCREAM_FIELD_UTILS_HPP
#define SCREAM_FIELD_UTILS_HPP

#include "share/util/eamxx_scalar_wrapper.hpp"
#include "share/field/field.hpp"

namespace scream {

// Check that two fields store the same entries.
// NOTE: if the field is padded, padding entries are NOT checked.
bool views_are_equal(const Field& f1, const Field& f2, const ekat::Comm* comm = nullptr);

// Uniform or normal randomization utils
void randomize_uniform (const Field& f, int seed,
                        const ScalarWrapper lb = ScalarWrapper::zero(),
                        const ScalarWrapper ub = ScalarWrapper::one());
void randomize_normal (const Field& f, int seed,
                       const ScalarWrapper mean = ScalarWrapper::zero(),
                       const ScalarWrapper sigma = ScalarWrapper::one());
void randomize_discrete (const Field& f, int seed,
                         const std::vector<ScalarWrapper>& values);
template<typename ST>
void randomize_discrete (const Field& f, int seed,
                         const std::vector<ST>& values)
{
  std::vector<ScalarWrapper> op_values;
  for (auto v : values)
    op_values.push_back(ScalarWrapper(v));
  return randomize_discrete(f,seed,op_values);
}

// Compute a random perturbation of a field for all view entries whose
// level index satisfies the mask.
// Input:
//   - f:             Field to perturbed. Required to have level midpoint
//                    tag as last dimension.
//   - perturb_level: magnitude of RELATIVE perturbation. That is,
//                    field *= (1 + U(1-perturb_level,1+perturb_level))
//   - base_seed:     Seed used for creating the engine input.
//   - level_mask:    A 1d int field (size of the level dimension of f) where f(i0,...,k) is
//                    perturbed if level_mask(k)=1
//   - dof_gids:      Field containing global DoF IDs for columns of f (if applicable)
void perturb (Field& f,
              const ScalarWrapper perturb_level,
              const int base_seed,
              const Field& level_mask,
              const Field& dof_gids = Field());

// Vertical/horizontal contractions of field (possibly averaging)
void vert_contraction (const Field& f_out, const Field& f_in, const Field& weight, bool AVG = false);
void horiz_contraction(const Field& f_out, const Field& f_in, const Field& weight, bool AVG = true,
                       const ekat::Comm* comm = nullptr, const Field& f_tmp = Field());

// Reduce field to a single scalar, and return an opaque type, to allow hiding impl in cpp file.
// NOTE: all calculations are done serially HOST
ScalarWrapper frobenius_norm(const Field& f, const ekat::Comm* comm = nullptr);
ScalarWrapper field_sum(const Field& f, const ekat::Comm* comm = nullptr);
ScalarWrapper field_max(const Field& f, const ekat::Comm* comm = nullptr);
ScalarWrapper field_min(const Field& f, const ekat::Comm* comm = nullptr);

// Prints the value of a field at a certain location, specified by tags and indices.
// If the field layout contains all the location tags, we will slice the field along
// those tags, and print it. E.g., f might be a <COL,LEV> field, and the tags/indices
// refer to a single column, in which case we'll print a whole column worth of data.
void print_field_hyperslab (const Field& f,
                            const std::vector<FieldTag>& tags = {},
                            const std::vector<int>& indices = {},
                            std::ostream& out = std::cout);

// Compute where input field comparese correctly to given value
void compute_mask (const Field& x, const ScalarWrapper value, Comparison CMP, Field& mask);
Field compute_mask (const Field& x, const ScalarWrapper value, Comparison CMP);

// Transpose a field layout
void transpose (const Field& src, Field& tgt);
Field transpose (const Field& src);

} // namespace scream

#endif // SCREAM_FIELD_UTILS_HPP
