#include <catch2/catch.hpp>
#include <numeric>

#include "share/field/field.hpp"
#include "share/field/field_pyutils.hpp"
#include "share/field/field_utils.hpp"

#include <pybind11/pybind11.h>

#include <filesystem>

namespace py = pybind11;

namespace scream {

TEST_CASE("pysession", "") {
  auto& ps = PySession::get();

  REQUIRE_THROWS(ps.finalize());
  REQUIRE (not ps.is_initialized());

  ps.initialize();
  REQUIRE (ps.is_initialized());
  ps.finalize();
  REQUIRE (not ps.is_initialized());
}

TEST_CASE("pyfield", "") {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  std::vector<FieldTag> tags = {COL,CMP,LEV};
  const int ncol = 2;
  const int ndim = 3;
  const int nlev = 4;
  std::vector<int> dims = {ncol, ndim, nlev};

  // Get the current file path
  std::filesystem::path curr_file(__FILE__);
  std::filesystem::path curr_src_path = curr_file.parent_path();

  // auto& ps = PySession::get();

  // ps.initialize();
  // ps.add_path(curr_src_path);

  FieldIdentifier fid ("field", {tags,dims}, m/s,"some_grid");
  Field f1 (fid);
  f1.allocate_view();

  auto f2 = f1.clone();
  auto f2_h = f2.get_view<double***,Host>();
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
    py::scoped_interpreter guard;
    auto sys = py::module::import("sys");
    auto sysPath = sys.attr("path");
    // sysPath.append(__FILE__);
    std::cout << "sys.path entries:" << std::endl;
    for (const auto& pathEntry : sysPath) {
        std::cout << py::cast<std::string>(pathEntry) << std::endl;
    }
    // std::cout << "py inited: " << ps.is_initialized() << "\n";
    auto py_mod = py::module::import("pyfield");
    // std::cout << "py inited: " << ps.is_initialized() << "\n";

    // auto f_py = create_py_field<Host>(f1);
    // py_mod.attr("set_iota")(f_py);
    // py_mod.attr("hello")();
  }
  // REQUIRE (views_are_equal(f1,f2));

  // ps.finalize();
}

} // namespace scream
