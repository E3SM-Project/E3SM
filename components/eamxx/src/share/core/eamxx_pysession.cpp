#include "share/core/eamxx_pysession.hpp"

#include <ekat_assert.hpp>
#include <ekat_fpe.hpp>

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

pybind11::module PySession::safe_import (const std::string& module_name) const
{
  pybind11::module m;

  // Disable FPEs while loading the module, then immediately re-enable them
  auto fpes = ekat::get_enabled_fpes();
  ekat::disable_all_fpes();
  try {
    m = pybind11::module::import(module_name.c_str());
  } catch (const pybind11::error_already_set& e) {
    std::cout << "[PySession::safe_import] Error! Python module import failed.\n"
                 " - module name: " + module_name + "\n"
                 " - pybind11 error: " + std::string(e.what()) + "\n"
                 "Did you forget to call PySession::add_path to add the module location to sys.path?\n";
    throw e;
  }
  ekat::enable_fpes(fpes);

  EKAT_REQUIRE_MSG (not m.is_none(),
      "Error! Could not import module '" + module_name + "'.\n");

  return m;
}

} // namespace scream
