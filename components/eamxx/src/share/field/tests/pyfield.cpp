#include <catch2/catch.hpp>
#include <numeric>

#include "share/field/field.hpp"
#include "share/field/field_pyutils.hpp"
#include "share/field/field_utils.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/embed.h>

#include <filesystem>

namespace py = pybind11;

namespace scream {

TEST_CASE("pyfield", "") {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  auto& ps = PySession::get();
  ps.initialize();

  std::vector<FieldTag> tags = {COL,CMP,LEV};
  const int ncol = 2;
  const int ndim = 3;
  const int nlev = 4;
  std::vector<int> dims = {ncol, ndim, nlev};

  // Get the current file path and add it to sys.path in python
  std::filesystem::path curr_file(__FILE__);
  ps.add_path(curr_file.parent_path().string());

  FieldIdentifier fid ("field", {tags,dims}, m/s,"some_grid");
  Field f1 (fid);
  f1.allocate_view();

  auto f2 = f1.clone();
  auto f2_h = f2.get_view<Real***,Host>();
  for (int icol=0,k=0; icol<ncol; ++icol) {
    for (int idim=0; idim<ndim; ++idim) {
      for (int ilev=0; ilev<nlev; ++ilev, ++k) {
        f2_h(icol,idim,ilev) = k;
      }
    }
  }
  f2.sync_to_dev();

  REQUIRE (not views_are_equal(f1,f2));

  // Use scope, so all py structures are destroyed BEFORE py.finalize()
  {
    auto py_mod = ps.safe_import("pyfield");

    auto f_py = create_py_field<Host>(f1);
    py_mod.attr("set_iota")(f_py);
  }
  REQUIRE (views_are_equal(f1,f2));

  ps.finalize();
}

} // namespace scream
