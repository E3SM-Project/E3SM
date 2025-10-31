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
template<typename ST>
void randomize_uniform (const Field& f, int seed, const ST& lb, const ST& ub)
{
  return randomize_uniform(f,seed,ScalarWrapper(lb),ScalarWrapper(ub));
}
void randomize_normal (const Field& f, int seed,
                       const ScalarWrapper mean = ScalarWrapper::zero(),
                       const ScalarWrapper sigma = ScalarWrapper::one());
template<typename ST>
void randomize_normal (const Field& f, int seed, const ST& mean, const ST& sigma = 1)
{
  return randomize_normal(f,seed,ScalarWrapper(mean),ScalarWrapper(sigma));
}
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
//   - lb/ub:         Lower/upper bound for the uniform distribution used
//                    to generate the RELATIVE field perturbation. That is,
//                    field_view(i0,...,iN) *= (1 + U(lb,ub))
//   - base_seed:     Seed used for creating the engine input.
//   - level_mask:    Mask (size of the level dimension of f) where f(i0,...,k) is
//                    perturbed if level_mask(k)=true
//   - dof_gids:      Field containing global DoF IDs for columns of f (if applicable)
template<typename MaskType>
void perturb (Field& f,
              const ScalarWrapper lb, const ScalarWrapper ub,
              const int base_seed,
              const MaskType& level_mask,
              const Field& dof_gids = Field())
{
  EKAT_REQUIRE_MSG(f.is_allocated(),
       	           "Error! Cannot perturb the values of a field not yet allocated.\n");

  using namespace ShortFieldTagsNames;
  const auto& fl = f.get_header().get_identifier().get_layout();

  // Field we are perturbing should have a level dimension,
  // and it is required to be the last dimension
  EKAT_REQUIRE_MSG(fl.rank()>0 &&
                   (fl.tags().back() == LEV || fl.tags().back() == ILEV),
                   "Error! Trying to perturb field \""+f.name()+"\", but field "
	                 "does not have LEV or ILEV as last dimension.\n"
                   "  - field name: " + f.name() + "\n"
                   "  - field layout: " + fl.to_string() + "\n");

  if (fl.has_tag(COL)) {
    // If field has a column dimension, it should be the first dimension
    EKAT_REQUIRE_MSG(fl.tag(0) == COL,
                     "Error! Trying to perturb field \""+f.name()+"\", but field "
	                   "does not have COL as first dimension.\n"
                     "  - field name: " + f.name() + "\n"
                     "  - field layout: " + fl.to_string() + "\n");

    const auto& dof_gids_fl = dof_gids.get_header().get_identifier().get_layout();
    EKAT_REQUIRE_MSG(dof_gids_fl.dim(0) == fl.dim(COL),
                     "Error! Field of DoF GIDs should have the same size as "
                     "perturbed field's column dimension.\n"
                     "  - dof_gids dim: " + std::to_string(dof_gids_fl.dim(0)) + "\n"
                     "  - field name: " + f.name() + "\n"
                     "  - field layout: " + fl.to_string() + "\n");
    EKAT_REQUIRE_MSG(dof_gids.data_type() == DataType::IntType,
                     "Error! DoF GIDs field must have \"int\" as data type.\n");
  }

  // Check to see if field has a column dimension
  using namespace ShortFieldTagsNames;
  const bool has_column_dim = fl.has_tag(COL);
  const bool has_lev_dim = fl.has_tag(LEV);

  if (has_column_dim) {
    // Because Column is the partitioned dimension, we must reset the
    // RNG seed to be the same on every column so that a column will
    // have the same value no matter where it exists in an MPI rank's
    // set of local columns.
    const auto gids = dof_gids.get_strided_view<const int*, Host>();

    // Create a field to store perturbation values with layout
    // the same as f, but stripped of column and level dimension.
    auto perturb_fl = fl.clone().strip_dims({COL,LEV});
    FieldIdentifier perturb_fid("perturb_field", perturb_fl, ekat::units::Units::nondimensional(), "");
    Field perturb_f(perturb_fid);
    perturb_f.allocate_view();

    // Loop through columns as reset RNG seed based on GID of column
    for (auto icol=0; icol<fl.dims().front(); ++icol) {
      const auto new_seed = base_seed+gids(icol);

      if (has_lev_dim) {
        // Loop through levels. For each that satisfy the level_mask,
        // apply a random perturbation to f.
        for (auto ilev=0; ilev<fl.dims().back(); ++ilev) {
          if (level_mask(ilev)) {
            randomize_uniform(perturb_f, new_seed+ilev, lb, ub);
            f.subfield(COL, icol).subfield(LEV, ilev).scale(perturb_f);
          }
        }
      } else {
        randomize_uniform(perturb_f, new_seed, lb, ub);
        f.subfield(COL, icol).scale(perturb_f);
      }
    }
  } else {
    // If no Column tag exists, this field is not partitioned.

    // Create a field to store perturbation values with layout
    // the same as f, but stripped of level dimension.
    auto perturb_fl = fl.clone().strip_dim(LEV);
    FieldIdentifier perturb_fid("perturb_field", perturb_fl, ekat::units::Units::nondimensional(), "");
    Field perturb_f(perturb_fid);
    perturb_f.allocate_view();

    if (has_lev_dim) {
      // Loop through levels. For each that satisfy the level_mask,
      // apply a random perturbation to f.
      for (auto ilev=0; ilev<fl.dims().back(); ++ilev) {
        if (level_mask(ilev)) {
          randomize_uniform(perturb_f, base_seed+ilev, lb, ub);
          f.subfield(LEV, ilev).scale(perturb_f);
        }
      }
    } else {
      randomize_uniform(perturb_f, base_seed, lb, ub);
      f.scale(perturb_f);
    }
  }
}

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

void compute_mask (const Field& x, const ScalarWrapper value, Comparison CMP, Field& mask);

template<typename ST>
void compute_mask (const Field& x, const ST value, Comparison CMP, Field& mask) {
  compute_mask(x,ScalarWrapper(value),CMP,mask);
}

template<typename ST>
Field compute_mask (const Field& x, const ST value, Comparison CMP) {
  const auto& fid_x = x.get_header().get_identifier();
  const auto nondim = ekat::units::Units::nondimensional();
  FieldIdentifier fid(x.name()+"_mask",fid_x.get_layout(),nondim,fid_x.get_grid_name(),DataType::IntType);
  Field mask(fid,true);

  compute_mask(x,value,CMP,mask);
  return mask;
}

void transpose (const Field& src, Field& tgt);
Field transpose (const Field& src);

} // namespace scream

#endif // SCREAM_FIELD_UTILS_HPP
