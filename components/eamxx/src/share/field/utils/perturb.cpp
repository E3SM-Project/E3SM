#include "share/field/field_utils.hpp"

namespace scream {

namespace impl {

template<typename ST>
void perturb (Field& f, const ST perturb_level,
              const int base_seed,
              const Field& level_mask,
              const Field& dof_gids)
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

  ST ub = 1 + perturb_level;
  ST lb = 1 - perturb_level;

  Field::view_host_t<const int*> lev_mask_h;
  if (has_lev_dim)
    lev_mask_h = level_mask.get_view<const int*,Host>();

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
      const auto seed = base_seed+gids(icol);

      if (has_lev_dim) {
        // Loop through levels. For each that satisfy the level_mask,
        // apply a random perturbation to f.
        for (auto ilev=0; ilev<fl.dims().back(); ++ilev) {
          if (lev_mask_h(ilev)) {
            randomize_uniform(perturb_f, seed+ilev, lb, ub);
            f.subfield(COL, icol).subfield(LEV, ilev).scale(perturb_f);
          }
        }
      } else {
        randomize_uniform(perturb_f, seed, lb, ub);
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
        if (lev_mask_h(ilev)) {
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

} // namespace impl

void perturb (Field& f, const ScalarWrapper perturb_level,
              const int base_seed,
              const Field& level_mask,
              const Field& dof_gids)
{
  using namespace ShortFieldTagsNames;

  // Perform some sanity checks, before handing off to the impl method
  EKAT_REQUIRE_MSG(f.is_allocated(),
    "[perturb] Error! Trying to perturb the values of a field not yet allocated.\n"
    " - field name: " + f.name() + "\n");

  const auto& fl = f.get_header().get_identifier().get_layout();

  // Field we are perturbing should have a level dimension,
  // and it is required to be the last dimension
  EKAT_REQUIRE_MSG(fl.rank()>0 && (fl.tags().back() == LEV || fl.tags().back() == ILEV),
    "[perturb] Error! Field must have LEV or ILEV as its last dimension.\n"
    " - field name: " + f.name() + "\n"
    " - field layout: " + fl.to_string() + "\n");

  if (fl.has_tag(COL)) {
    // If field has a column dimension, it should be the first dimension
    EKAT_REQUIRE_MSG(fl.tag(0) == COL,
      "[perturb] Error! Fields with COL dimension should have it as their first dimension.\n"
      " - field name: " + f.name() + "\n"
      " - field layout: " + fl.to_string() + "\n");

    const auto& dof_gids_fl = dof_gids.get_header().get_identifier().get_layout();
    EKAT_REQUIRE_MSG(dof_gids_fl.dim(0) == fl.dim(COL),
      "Error! Field of DoF GIDs should have the same size as perturbed field's column dimension.\n"
      " - field name: " + f.name() + "\n"
      " - field layout: " + fl.to_string() + "\n"
      " - dof_gids layout: " + dof_gids_fl.to_string() + "\n");
    EKAT_REQUIRE_MSG(dof_gids.data_type() == DataType::IntType,
      "Error! DoF GIDs field must have \"int\" as data type.\n"
      " - dof gids data type: " + e2str(dof_gids.data_type()));
  }

  switch (f.data_type()) {
    case DataType::IntType:
      impl::perturb(f,perturb_level.as<int>(),base_seed,level_mask,dof_gids);
      break;
    case DataType::FloatType:
      impl::perturb(f,perturb_level.as<float>(),base_seed,level_mask,dof_gids);
      break;
    case DataType::DoubleType:
      impl::perturb(f,perturb_level.as<double>(),base_seed,level_mask,dof_gids);
      break;
    default:
      EKAT_ERROR_MSG ("[perturb] Error! Unsupported field data type.\n");
  }
}

} // namespace scream
