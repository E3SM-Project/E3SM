#include <catch2/catch.hpp>

#include "share/field/field_reader.hpp"
#include "share/field/field.hpp"
#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"

#include <ekat_comm.hpp>

#include <numeric>

// ============================================================ //
//  Helpers
// ============================================================ //

namespace scream {

namespace {

// Build and allocate a field with a given layout and data type.
Field make_field (const std::string& name,
                  const FieldLayout& layout,
                  const DataType dtype = DataType::RealType)
{
  return Field(FieldIdentifier(name, layout, ekat::units::none, "grid", dtype),true);
}

// Fill a field with sequentially increasing values starting from `start`.
void iota (Field& f, const ScalarWrapper& init)
{
  const int n = f.get_header().get_alloc_properties().get_alloc_size()
              / f.get_header().get_alloc_properties().get_scalar_size();
  switch (f.data_type()) {
    case DataType::IntType:
    {
      auto* ptr = f.get_internal_view_data<int, Host>();
      std::iota(ptr, ptr + n, init.as<int>());
      break;
    }
    case DataType::FloatType:
    {
      auto* ptr = f.get_internal_view_data<float, Host>();
      std::iota(ptr, ptr + n, init.as<float>());
      break;
    }
    case DataType::DoubleType:
    {
      auto* ptr = f.get_internal_view_data<double, Host>();
      std::iota(ptr, ptr + n, init.as<double>());
      break;
    }
    default:
      EKAT_ERROR_MSG ("Unrecognized/unsupported data type.\n"
                      " - field name: " + f.name() + "\n"
                      " - data type : " + e2str(f.data_type()) + "\n");
  }
  f.sync_to_dev();
}

// Write a single variable to an already-open (and enddef'd) scorpio file.
template <typename ST>
void write_field_to_file (const std::string& filename, const Field& f)
{
  f.sync_to_host();
  scorpio::write_var(filename, f.name(),
                     f.get_internal_view_data<const ST, Host>());
}

} // anonymous namespace

// ============================================================ //
//  Test: non-decomposed read
// ============================================================ //

TEST_CASE ("read_fields_no_decomp")
{
  using namespace ShortFieldTagsNames;

  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);

  const std::string filename =
      "field_io_test_no_decomp_np" + std::to_string(comm.size()) + ".nc";

  // Layouts
  const int ncols = 4;
  const int nlevs = 6;
  FieldLayout lay2d({COL, LEV}, {ncols, nlevs});
  FieldLayout lay1d({LEV},      {nlevs});

  // ---- Write phase ---- //
  {
    scorpio::register_file(filename, scorpio::Write);
    scorpio::define_dim(filename, "ncol", ncols);
    scorpio::define_dim(filename, "lev",  nlevs);
    scorpio::define_var(filename, "f2d", {"ncol", "lev"}, "real", false);
    scorpio::define_var(filename, "f1d", {"lev"},          "real", false);
    scorpio::enddef(filename);

    Field f2d = make_field("f2d", lay2d);
    Field f1d = make_field("f1d", lay1d, DataType::IntType);
    iota(f2d,0);
    iota(f1d,100);

    write_field_to_file<Real>(filename, f2d);
    write_field_to_file<int>(filename, f1d);

    scorpio::release_file(filename);
  }

  // ---- Read phase ---- //
  {
    Field f2d = make_field("f2d", lay2d);
    Field f1d = make_field("f1d", lay1d, DataType::IntType);
    f2d.deep_copy(-1);
    f1d.deep_copy(-1);

    read_fields(filename, {f2d, f1d});

    // Verify
    f2d.sync_to_host();
    f1d.sync_to_host();

    auto v2d = f2d.get_view<const Real**, Host>();
    auto v1d = f1d.get_view<const int*,  Host>();

    int idx = 0;
    for (int i = 0; i < ncols; ++i) {
      for (int j = 0; j < nlevs; ++j, ++idx) {
        REQUIRE (v2d(i, j) == Real(idx));
      }
    }
    for (int j = 0; j < nlevs; ++j) {
      REQUIRE (v1d(j) == Real(100 + j));
    }
  }

  scorpio::finalize_subsystem();
}

// ============================================================ //
//  Test: decomposed read (distributed column dimension)
// ============================================================ //

TEST_CASE ("read_fields_with_decomp", "[field_utils][io]")
{
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);

  const std::string filename =
      "field_io_test_decomp_np" + std::to_string(comm.size()) + ".nc";

  // Global sizes: use a multiple of comm.size() so decomp divides evenly
  const int lcols_per_rank = 3;
  const int ncols_global   = lcols_per_rank * comm.size();
  const int nlevs          = 4;

  FieldLayout lay2d({COL, LEV}, {lcols_per_rank, nlevs});
  FieldLayout lay1d({COL},      {lcols_per_rank});

  // GIDs owned by this rank (1-based)
  auto my_gids = make_field("gids",lay1d,DataType::IntType);
  auto my_gids_beg = my_gids.get_internal_view_data<int,Host>();
  auto my_gids_end = my_gids_beg + lcols_per_rank;
  std::iota(my_gids_beg,my_gids_end, 1 + comm.rank() * lcols_per_rank);
  my_gids.sync_to_dev();

  // ---- Write phase (rank 0 writes full data) ---- //
  {
    scorpio::register_file(filename, scorpio::Write);
    scorpio::define_dim(filename, "ncol", ncols_global);
    scorpio::define_dim(filename, "lev",  nlevs);

    // Decomposition offsets
    std::vector<scorpio::offset_t> offsets(lcols_per_rank);
    for (int i = 0; i < lcols_per_rank; ++i) {
      offsets[i] = my_gids_beg[i] - 1;
    }
    scorpio::set_dim_decomp(filename, "ncol", offsets);

    scorpio::define_var(filename, "f2d", {"ncol", "lev"}, "real", false);
    scorpio::define_var(filename, "f1d", {"ncol"},        "real", false);
    scorpio::enddef(filename);

    // Fill: global value at (gcol, lev) = gcol*nlevs + lev
    //       global value at (gcol)      = gcol
    Field f2d = make_field("f2d", lay2d);
    Field f1d = make_field("f1d", lay1d, DataType::IntType);
    f2d.sync_to_host();
    f1d.sync_to_host();

    auto v2d = f2d.get_view<Real**, Host>();
    auto v1d = f1d.get_view<int*,  Host>();
    for (int li = 0; li < lcols_per_rank; ++li) {
      const int gcol = my_gids_beg[li];
      v1d(li) = gcol;
      for (int j = 0; j < nlevs; ++j) {
        v2d(li, j) = Real(gcol * nlevs + j);
      }
    }
    f2d.sync_to_dev();
    f1d.sync_to_dev();

    write_field_to_file<Real>(filename, f2d);
    write_field_to_file<int>(filename, f1d);

    scorpio::release_file(filename);
  }

  // ---- Read phase ---- //
  {
    Field f2d = make_field("f2d", lay2d);
    Field f1d = make_field("f1d", lay1d, DataType::IntType);
    f2d.deep_copy(-1);
    f1d.deep_copy(-1);

    read_fields(filename, {f2d, f1d}, my_gids, comm);

    // Verify local data
    f2d.sync_to_host();
    f1d.sync_to_host();

    auto v2d = f2d.get_view<const Real**, Host>();
    auto v1d = f1d.get_view<const int*,  Host>();

    for (int li = 0; li < lcols_per_rank; ++li) {
      const int gcol = my_gids_beg[li];
      REQUIRE (v1d(li) == gcol);
      for (int j = 0; j < nlevs; ++j) {
        REQUIRE (v2d(li, j) == Real(gcol * nlevs + j));
      }
    }
  }

  scorpio::finalize_subsystem();
}

// ============================================================ //
//  Test: reading a specific time slice
// ============================================================ //

TEST_CASE ("read_fields_time_slice")
{
  using namespace ShortFieldTagsNames;

  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);

  const std::string filename =
      "field_io_test_time_np" + std::to_string(comm.size()) + ".nc";

  const int nlevs  = 5;
  const int ntimes = 3;
  FieldLayout lay1d({LEV}, {nlevs});

  // ---- Write phase ---- //
  {
    scorpio::register_file(filename, scorpio::Write);
    scorpio::define_time(filename, "s");
    scorpio::define_dim(filename, "lev", nlevs);
    scorpio::define_var(filename, "f1d", {"lev"}, "real", true);
    scorpio::enddef(filename);

    Field f1d = make_field("f1d", lay1d);

    for (int t = 0; t < ntimes; ++t) {
      scorpio::update_time(filename, double(t));
      f1d.deep_copy(Real(t * 10));
      write_field_to_file<Real>(filename, f1d);
    }

    scorpio::release_file(filename);
  }

  // ---- Read phase: request middle time slice ---- //
  {
    const int target_t = 1;

    Field f1d = make_field("f1d", lay1d);
    f1d.deep_copy(Real(-1));

    read_fields(filename, {f1d}, target_t);

    f1d.sync_to_host();
    auto v = f1d.get_view<const Real*, Host>();
    for (int j = 0; j < nlevs; ++j) {
      REQUIRE (v(j) == Real(target_t * 10));
    }
  }

  scorpio::finalize_subsystem();
}

} // namespace scream
