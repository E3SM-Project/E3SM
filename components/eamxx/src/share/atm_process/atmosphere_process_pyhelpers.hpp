#ifndef SCREAM_ATMOSPHERE_PROCESS_PYHELPERS_HPP
#define SCREAM_ATMOSPHERE_PROCESS_PYHELPERS_HPP

#include "share/atm_process/atmosphere_process.hpp"

#ifdef EAMXX_HAS_PYTHON
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#endif

namespace scream
{

template<typename... Args>
void AtmosphereProcess::
py_module_call (const std::string& name, const Args&... args)
{
  const auto& py_module = std::any_cast<const pybind11::module&>(m_py_module);
  py_module.attr(name.c_str())(args...);
}

inline const pybind11::array& AtmosphereProcess::
get_py_field_impl (const strmap_t<strmap_t<std::any>>& py_fields,
                   const std::string& fname, const std::string& grid) const
{
  auto any_f = py_fields.at(fname).at(grid);
  return std::any_cast<const pybind11::array&>(any_f);
}

inline const pybind11::array& AtmosphereProcess::
get_py_field_host (const std::string& fname, const std::string& grid) const
{
  return get_py_field_impl(m_py_fields_host,fname,grid);
}

inline const pybind11::array& AtmosphereProcess::
get_py_field_dev (const std::string& fname, const std::string& grid) const
{
  return get_py_field_impl(m_py_fields_dev,fname,grid);
}

inline const pybind11::array& AtmosphereProcess::
get_py_field_host (const std::string& fname) const
{
  EKAT_REQUIRE_MSG (m_py_fields_host.at(fname).size()==1,
      "Cannot request pyfield by field name only. Multiple copies exist on multiple grids.\n"
      "  - field name: " + fname + "\n");
  return std::any_cast<const pybind11::array&>(m_py_fields_host.at(fname).begin()->second);
}

inline const pybind11::array& AtmosphereProcess::
get_py_field_dev (const std::string& fname) const
{
  EKAT_REQUIRE_MSG (m_py_fields_dev.at(fname).size()==1,
      "Cannot request pyfield by field name only. Multiple copies exist on multiple grids.\n"
      "  - field name: " + fname + "\n");
  return std::any_cast<const pybind11::array&>(m_py_fields_dev.at(fname).begin()->second);
}

} // namespace scream

#endif // SCREAM_ATMOSPHERE_PROCESS_PYHELPERS_HPP
