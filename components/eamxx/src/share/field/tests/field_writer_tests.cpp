#include <catch2/catch.hpp>

#include "share/field/field_writer.hpp"
#include "share/field/field_reader.hpp"
#include "share/field/field.hpp"
#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"

#include <ekat_comm.hpp>

#include <numeric>

namespace scream {

namespace {

Field make_field (const std::string& name,
                  const FieldLayout& layout,
                  const DataType dtype = DataType::RealType)
{
  return Field(FieldIdentifier(name, layout, ekat::units::none, "grid", dtype),true);
}

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

} // anonymous namespace

TEST_CASE ("write_fields_no_decomp")
{
  using namespace ShortFieldTagsNames;

  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);

  const std::string filename =
      "field_writer_test_no_decomp_np" + std::to_string(comm.size()) + ".nc";

  const int ncols = 4;
  const int nlevs = 6;
  FieldLayout lay2d({COL, LEV}, {ncols, nlevs});
  FieldLayout lay1d({LEV},      {nlevs});

  Field f2d = make_field("f2d", lay2d);
  Field f1d = make_field("f1d", lay1d, DataType::IntType);
  iota(f2d,0);
  iota(f1d,100);

  FieldWriter writer;
  writer.set_fields({f2d,f1d});
  writer.set_file_specs(filename);
  writer.init_scorpio_structures();
  writer.write();

  Field f2d_in = make_field("f2d", lay2d);
  Field f1d_in = make_field("f1d", lay1d, DataType::IntType);
  f2d_in.deep_copy(-1);
  f1d_in.deep_copy(-1);

  read_fields(filename,{f2d_in,f1d_in});

  f2d_in.sync_to_host();
  f1d_in.sync_to_host();
  auto v2d = f2d_in.get_view<const Real**, Host>();
  auto v1d = f1d_in.get_view<const int*, Host>();

  int idx = 0;
  for (int i = 0; i < ncols; ++i) {
    for (int j = 0; j < nlevs; ++j, ++idx) {
      REQUIRE (v2d(i,j)==Real(idx));
    }
  }
  for (int j = 0; j < nlevs; ++j) {
    REQUIRE (v1d(j)==100+j);
  }

  scorpio::finalize_subsystem();
}

TEST_CASE ("write_fields_with_decomp")
{
  using namespace ShortFieldTagsNames;

  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);

  const std::string filename =
      "field_writer_test_decomp_np" + std::to_string(comm.size()) + ".nc";

  const int lcols_per_rank = 3;
  const int ncols_global   = lcols_per_rank * comm.size();
  const int nlevs          = 4;

  FieldLayout lay2d({COL, LEV}, {lcols_per_rank, nlevs});
  FieldLayout lay1d({COL},      {lcols_per_rank});

  auto my_gids = make_field("gids",lay1d,DataType::IntType);
  auto my_gids_beg = my_gids.get_internal_view_data<int,Host>();
  auto my_gids_end = my_gids_beg + lcols_per_rank;
  std::iota(my_gids_beg,my_gids_end, 1 + comm.rank() * lcols_per_rank);
  my_gids.sync_to_dev();

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
      v2d(li,j) = Real(gcol * nlevs + j);
    }
  }
  f2d.sync_to_dev();
  f1d.sync_to_dev();

  FieldWriter writer;
  writer.set_file_specs(filename);
  writer.set_dim_decomp(my_gids,comm);
  writer.set_fields({f2d,f1d});
  writer.init_scorpio_structures();
  writer.write();

  Field f2d_in = make_field("f2d", lay2d);
  Field f1d_in = make_field("f1d", lay1d, DataType::IntType);
  f2d_in.deep_copy(-1);
  f1d_in.deep_copy(-1);

  read_fields(filename,{f2d_in,f1d_in},my_gids,comm);

  f2d_in.sync_to_host();
  f1d_in.sync_to_host();
  auto v2d_in = f2d_in.get_view<const Real**, Host>();
  auto v1d_in = f1d_in.get_view<const int*, Host>();

  for (int li = 0; li < lcols_per_rank; ++li) {
    const int gcol = my_gids_beg[li];
    REQUIRE (v1d_in(li)==gcol);
    for (int j = 0; j < nlevs; ++j) {
      REQUIRE (v2d_in(li,j)==Real(gcol * nlevs + j));
    }
  }

  scorpio::finalize_subsystem();
}

TEST_CASE ("write_fields_time_dep")
{
  using namespace ShortFieldTagsNames;

  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);

  const std::string filename =
      "field_writer_test_time_np" + std::to_string(comm.size()) + ".nc";

  const int nlevs  = 5;
  const int ntimes = 3;
  FieldLayout lay1d({LEV}, {nlevs});

  Field f1d = make_field("f1d", lay1d);

  FieldWriter writer;
  writer.set_fields({f1d});
  writer.set_time_dependent(true,"s","time");
  writer.set_file_specs(filename);
  writer.init_scorpio_structures();

  for (int t = 0; t < ntimes; ++t) {
    f1d.deep_copy(Real(t*10));
    writer.write(double(t));
  }

  const int target_t = 1;
  Field f1d_in = make_field("f1d", lay1d);
  f1d_in.deep_copy(-1);

  read_fields(filename,{f1d_in},target_t);

  f1d_in.sync_to_host();
  auto v = f1d_in.get_view<const Real*, Host>();
  for (int j = 0; j < nlevs; ++j) {
    REQUIRE (v(j)==Real(target_t*10));
  }

  scorpio::finalize_subsystem();
}

TEST_CASE ("write_fields_append_and_time_behavior_mismatch")
{
  using namespace ShortFieldTagsNames;

  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);

  const std::string filename =
      "field_writer_test_append_np" + std::to_string(comm.size()) + ".nc";

  const int nlevs = 4;
  FieldLayout lay({LEV}, {nlevs});

  Field f1 = make_field("f1", lay);
  f1.deep_copy(Real(1));

  {
    FieldWriter writer;
    writer.set_fields({f1});
    writer.set_file_specs(filename);
    writer.init_scorpio_structures();
    writer.write();
  }

  Field f2 = make_field("f2", lay, DataType::IntType);
  iota(f2,7);

  {
    FieldWriter writer;
    writer.set_fields({f2});
    writer.set_file_specs(filename);
    writer.init_scorpio_structures();
    writer.write();
  }

  Field f2_in = make_field("f2", lay, DataType::IntType);
  f2_in.deep_copy(-1);
  read_fields(filename,{f2_in});
  f2_in.sync_to_host();
  auto v2 = f2_in.get_view<const int*, Host>();
  for (int i=0; i<nlevs; ++i) {
    REQUIRE (v2(i)==7+i);
  }

  {
    FieldWriter writer;
    writer.set_fields({f1});
    writer.set_time_dependent(true,"s","time");
    writer.set_file_specs(filename);
    REQUIRE_THROWS(writer.init_scorpio_structures());
  }

  scorpio::finalize_subsystem();
}

TEST_CASE ("write_fields_append_new_time_dep_var")
{
  using namespace ShortFieldTagsNames;

  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);

  const std::string filename =
      "field_writer_test_append_time_np" + std::to_string(comm.size()) + ".nc";

  const int nlevs  = 3;
  const int ntimes = 2;
  FieldLayout lay({LEV}, {nlevs});

  Field f1 = make_field("f1", lay);
  FieldWriter w1;
  w1.set_fields({f1});
  w1.set_time_dependent(true,"s","time");
  w1.set_file_specs(filename);
  w1.init_scorpio_structures();
  for (int t=0; t<ntimes; ++t) {
    f1.deep_copy(Real(10+t));
    w1.write(double(t));
  }

  Field f2 = make_field("f2", lay, DataType::IntType);
  FieldWriter w2;
  w2.set_fields({f2});
  w2.set_time_dependent(true,"s","time");
  w2.set_file_specs(filename);
  w2.init_scorpio_structures();
  iota(f2,100);
  w2.write(double(ntimes));

  Field f2_in = make_field("f2", lay, DataType::IntType);
  f2_in.deep_copy(-1);
  read_fields(filename,{f2_in},ntimes);
  f2_in.sync_to_host();
  auto v = f2_in.get_view<const int*, Host>();
  for (int i=0; i<nlevs; ++i) {
    REQUIRE(v(i)==100+i);
  }

  scorpio::finalize_subsystem();
}

} // namespace scream
