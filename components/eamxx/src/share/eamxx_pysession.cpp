#include <eamxx_pysession.hpp>

#include <ekat/ekat_assert.hpp>
#include <Python.h>

#include <filesystem>

namespace scream {

void PySession::initialize () {
  if (num_customers==0) {
    if (not Py_IsInitialized()) {
      Py_Initialize();
      should_finalize = true;
    }
  }
  ++num_customers;
}
void PySession::finalize () {
  if (num_customers==1 and should_finalize) {
    Py_Finalize();
  }
  EKAT_REQUIRE_MSG (num_customers>0,
      "Error! Invalid number of customers.\n"
      "  Did you call PySession::finalize() without calling PySession::initialize()?\n");
  --num_customers;
}

void PySession::add_path (const std::string& path)
{
  EKAT_REQUIRE_MSG (is_initialized(),
      "Error! Cannot modify python's sys.path, since PySession was not initialized yet.\n");

  // Import the sys module
  PyObject *sysModule = PyImport_ImportModule("sys");
  if (sysModule == nullptr) {
      PyErr_Print();
      EKAT_ERROR_MSG("Could not import sys module. Aborting.\n");
  }

  // Get the sys.path list
  PyObject *sysPath = PyObject_GetAttrString(sysModule, "path");
  if (sysPath == nullptr) {
      PyErr_Print();
      Py_DECREF(sysModule);
      EKAT_ERROR_MSG("Could not access sys.path. Aborting.\n");
  }

  // Convert path to py object
  PyObject *pathItem = PyUnicode_FromString(path.c_str());
  if (pathItem == nullptr) {
      PyErr_Print();
      Py_DECREF(sysPath);
      Py_DECREF(sysModule);
      EKAT_ERROR_MSG("Could not convert path string. Aborting.\n");
  }

  // Append the new path to sys.path
  if (PyList_Append(sysPath, pathItem) != 0) {
      PyErr_Print();
      Py_DECREF(pathItem);
      Py_DECREF(sysPath);
      Py_DECREF(sysModule);
      EKAT_ERROR_MSG("Could not append to sys.path. Aborting.\n");
  }

  // Clean up
  Py_DECREF(pathItem);
  Py_DECREF(sysPath);
  Py_DECREF(sysModule);
}

void PySession::add_curr_path ()
{
  auto curr_path = std::filesystem::current_path();
  add_path(curr_path.string());
}

} // namespace scream
