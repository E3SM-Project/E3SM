#include <eamxx_pysession.hpp>

#include <ekat/ekat_assert.hpp>

#include <filesystem>
#include <pybind11/embed.h>

namespace py = pybind11;

namespace scream {

void PySession::initialize () {
  if (num_customers==0) {
    // Note: if Py interpreter is already inited, we ASSUME someone else
    // is handling the interpreter initialization/finalization
    if (not Py_IsInitialized()) {
      py_guard = std::make_shared<py::scoped_interpreter>();
    }
  }
  ++num_customers;
}
void PySession::finalize () {
  EKAT_REQUIRE_MSG (num_customers>0,
      "Error! Invalid number of customers.\n"
      "  Did you call PySession::finalize() without calling PySession::initialize()?\n");

  --num_customers;
  if (num_customers==0) {
    py_guard.reset();
  }
}

void PySession::add_path (const std::string& path)
{
  EKAT_REQUIRE_MSG (is_initialized(),
      "Error! Cannot modify python's sys.path, since PySession was not initialized yet.\n");

  try {
    // Import the sys module
    py::module sysModule = py::module::import("sys");

    // Get the sys.path list
    py::list sysPath = sysModule.attr("path");

    // Append the new path to sys.path
    sysPath.append(path);
  } catch (const py::error_already_set& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    throw std::runtime_error("Could not modify sys.path. Aborting.");
  }
}

void PySession::add_curr_path ()
{
  auto curr_path = std::filesystem::current_path();
  add_path(curr_path.string());
}

} // namespace scream
