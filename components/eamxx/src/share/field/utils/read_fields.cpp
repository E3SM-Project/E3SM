#include "share/field/field_utils.hpp"
#include "share/field/field_alloc_prop.hpp"
#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"

#include <ekat_string_utils.hpp>
#include <ekat_zip.hpp>

#include <algorithm>
#include <limits>

namespace scream {

namespace {

// Creates an I/O-friendly copy of a field: no padding and no subfield parent.
// If the field already satisfies these conditions it is returned as-is (alias).
// Otherwise a freshly allocated copy with the same identifier is returned.
std::vector<Field> make_io_fields (const std::vector<Field>& fields)
{
  std::vector<Field> io_fields;
  for (const auto& f : fields) {
    const auto& fh  = f.get_header();
    const auto& fap = fh.get_alloc_properties();
    const bool can_alias = fh.get_parent() == nullptr && fap.get_padding() == 0;
    if (can_alias)
      io_fields.emplace_back(f);
    else
      io_fields.emplace_back(fh.get_identifier(),true);
  }
  return io_fields;
}

// Handle possibly non-unlimited "time" dim, validate field extents against
// the file metadata, and changes each var dtype to "real" so that read_var
// can use the native precision.
void setup_file_and_vars (const std::string& filename,
                          const std::vector<Field>& io_fields)
{
  // Mark a non-unlimited "time" dimension so sliced reads work correctly.
  if (!scorpio::has_time_dim(filename) && scorpio::has_dim(filename, "time")) {
    scorpio::mark_dim_as_time(filename, "time");
  }

  for (const auto& f : io_fields) {
    const auto& fdims = f.get_header().get_identifier().get_layout().dims();
    const auto& fname = f.name();

    EKAT_REQUIRE_MSG (scorpio::has_var(filename, fname),
        "Error! Variable not found in input file.\n"
        " - filename: " + filename + "\n"
        " - varname : " + fname + "\n");

    const auto& var = scorpio::get_var(filename, fname);
    std::vector<int> vardims;
    for (auto d : var.dims) {
      switch (d->decomp_rank) {
        case 0: vardims.push_back(d->length); break;
        case 1: vardims.push_back(var.decomp->dim_decomp->offsets.size()); break;
        default:
          EKAT_ERROR_MSG ("[read_fields] Error! Unsupported (and unexpected) decomposition rank\n"
                          " - var name: " + var.name + "\n"
                          " - decomp  : " + var.decomp->name + "\n"
                          " - rank    : " + std::to_string(d->decomp_rank) + "\n");
      }
    }
    EKAT_REQUIRE_MSG (vardims == fdims,
        "Error! Layout mismatch for variable in input file.\n"
        " - filename     : " + filename + "\n"
        " - varname      : " + fname + "\n"
        " - expected dims: " + ekat::join(fdims, ",") + "\n"
        " - file dims    : " + ekat::join(vardims, ",") + "\n");

    // Allow scorpio to cast to the field's native type during read.
    scorpio::change_var_dtype(filename, fname, "real");
  }
}

// Read from file into io_fields, and copy into the fields (no-op for aliasing fields)
void do_read_and_copy (const std::string& filename,
                       const std::vector<Field>& fields,
                       const std::vector<Field>& io_fields,
                       const int time_index)
{
  for (auto [f, io_f] : ekat::zip(fields,io_fields)) {
    const auto& fname = io_f.name();

    switch (io_f.data_type()) {
      case DataType::DoubleType:
        scorpio::read_var(filename, fname, io_f.get_internal_view_data<double, Host>(), time_index);
        break;
      case DataType::FloatType:
        scorpio::read_var(filename, fname, io_f.get_internal_view_data<float, Host>(), time_index);
        break;
      case DataType::IntType:
        scorpio::read_var(filename, fname, io_f.get_internal_view_data<int, Host>(), time_index);
        break;
      default:
        EKAT_ERROR_MSG (
            "Error! Unsupported data type while reading field from file.\n"
            " - filename  : " + filename + "\n"
            " - field name: " + fname + "\n");
    }

    io_f.sync_to_dev();
    f.deep_copy(io_f); // If io_f is aliasing f, this is a no-op
  }
}

} // anonymous namespace

void read_fields (const std::string& filename,
                  const std::vector<Field>& fields,
                  const int time_index)
{
  if (fields.size()==0)
    return;

  scorpio::register_file(filename, scorpio::Read);
  auto io_fields = make_io_fields(fields);
  setup_file_and_vars(filename, io_fields);
  do_read_and_copy(filename, fields, io_fields, time_index);
  scorpio::release_file(filename);
}

void read_fields (const std::string& filename,
                  const std::vector<Field>& fields,
                  const Field& decomp_dim_gids,
                  const ekat::Comm& comm,
                  const int time_index)
{
  if (fields.size()==0)
    return;

  EKAT_REQUIRE_MSG (decomp_dim_gids.rank()==1,
      "Error! Decomposed dim gids field must have rank 1.\n"
      " - field name: " + decomp_dim_gids.name() + "\n"
      " - field rank: " + std::to_string(decomp_dim_gids.rank()) + "\n");
  EKAT_REQUIRE_MSG (decomp_dim_gids.data_type()==DataType::IntType,
      "Error! Decomposed dim gids field must have integer data type.\n"
      " - field name : " + decomp_dim_gids.name() + "\n"
      " - field dtype: " + e2str(decomp_dim_gids.data_type()) + "\n");

  // Compute the global minimum GID so offsets can be 0-based.
  auto global_min_gid = field_min(decomp_dim_gids,&comm).as<int>();
  auto gids = decomp_dim_gids.get_view<const int*,Host>();
  std::vector<scorpio::offset_t> offsets(gids.size());
  for (std::size_t i = 0; i < gids.size(); ++i) {
    offsets[i] = static_cast<scorpio::offset_t>(gids[i] - global_min_gid);
  }

  scorpio::register_file(filename, scorpio::Read);

  // Register the distributed decomposition along the requested dimension.
  const auto& decomp_dim = decomp_dim_gids.get_header().get_identifier().get_layout().name(0);
  scorpio::set_dim_decomp(filename, decomp_dim, offsets);

  auto io_fields = make_io_fields(fields);
  setup_file_and_vars(filename, io_fields);

  do_read_and_copy(filename, fields, io_fields, time_index);
  scorpio::release_file(filename);
}

} // namespace scream
