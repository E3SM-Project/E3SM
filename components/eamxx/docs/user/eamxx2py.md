# Calling python code from EAMxx atmosphere processes

## Requirements

In order to call python code from an EAMxx atmosphere process,
EAMxx must be build with the CMake option `EAMXX_ENABLE_PYTHON=ON`,
and the CMake variable `Python_EXECUTABLE` must point to a python3
executable, with python version >= 3.9. Additionally, the python package
`pybind11` must be installed (e.g., via pip or conda).

If `EAMXX_ENABLE_PYTHON=OFF`, none of the code that is needed to call
python from EAMxx will be compiled.

## Usage

If python support is enabled, every atmosphere process stores
data structures that can hold python-compatible arrays and modules.
During construction, if the input parameter list contains a non-trivial
entry for the key `py_module_name`, EAMxx will automatically set up
these data structures. In particular, EAMxx will

- create python-compatible arrays for each of the input/output/internal
  fields that are registered in the class. These are stored in two maps:
  `m_py_fields_dev` and `m_py_fields_host`, which store python-compatible
  arrays for the device and host views of the Field, respectively. The maps
  are in fact nested maps, so that the python-compatible array for field X
  on grid Y can be retrieved via `m_py_fields_host[Y][X]`.
- load the python module provided via parameter list, so that its interfaces
  can later be called during init/run phases. The module is then stored in
  the local member `m_py_module`. If the module is in a non-standard path,
  the parameter list entry `py_module_path` can be used to specify its path,
  which will be added to python's search path before loading the module.

Due to implementation details in the pybind11 library, and to avoid compiler warnings,
all the python-compatible data structures are stored wrapped inside `std::any` objects.
As such, the need to be properly casted to the correct underlying type before being used.
In particular, the fields and module can be casted as follows:

```c++
  auto& f = std::any_cast<pybind11::array&>(m_py_fields_host[grid_name][fname]);
  auto& pymod = std::any_cast<pybind11::module&>(m_py_module);
```

Once the module is available, a function from it can be called using the `attr` method of
the module. For instance, if the module had a function `run` that takes 2 arrays and a double
(in this order), it can be invoked via

```c++
 pymod.attr("run")(f1,f2,my_double);
```

where `f1` and `f2` are of type `pybind11::array` (e.g., casted from objects in `m_py_fields_host`).

## Example

We provided an example of how to use this feature in `eamxx_cld_fraction_process_interface.cpp`,
which is a very small and simple atmosphere process. We paste here the code, which shows how
to support both C++ and python implemenation in the same cpp file

```c++
#ifdef EAMXX_HAS_PYTHON
  if (m_py_module.has_value()) {
    // For now, we run Python code only on CPU
    const auto& py_fields = m_py_fields_host.at(m_grid->name());

    const auto& py_qi                = std::any_cast<const py::array&>(py_fields.at("qi"));
    const auto& py_liq_cld_frac      = std::any_cast<const py::array&>(py_fields.at("cldfrac_liq"));
    const auto& py_ice_cld_frac      = std::any_cast<const py::array&>(py_fields.at("cldfrac_ice"));
    const auto& py_tot_cld_frac      = std::any_cast<const py::array&>(py_fields.at("cldfrac_tot"));
    const auto& py_ice_cld_frac_4out = std::any_cast<const py::array&>(py_fields.at("cldfrac_ice_for_analysis"));
    const auto& py_tot_cld_frac_4out = std::any_cast<const py::array&>(py_fields.at("cldfrac_tot_for_analysis"));

    // Sync input to host
    liq_cld_frac.sync_to_host();

    const auto& py_module = std::any_cast<const py::module&>(m_py_module);
    double ice_threshold      = m_params.get<double>("ice_cloud_threshold");
    double ice_4out_threshold = m_params.get<double>("ice_cloud_for_analysis_threshold");
    py_module.attr("main")(ice_threshold,ice_4out_threshold,py_qi,py_liq_cld_frac,py_ice_cld_frac,py_tot_cld_frac,py_ice_cld_frac_4out,py_tot_cld_frac_4out);

    // Sync outputs to dev
    qi.sync_to_dev();
    liq_cld_frac.sync_to_dev();
    ice_cld_frac.sync_to_dev();
    tot_cld_frac.sync_to_dev();
    ice_cld_frac_4out.sync_to_dev();
    tot_cld_frac_4out.sync_to_dev();
  } else
#endif
  {
    auto qi_v                = qi.get_view<const Pack**>();
    auto liq_cld_frac_v      = liq_cld_frac.get_view<const Pack**>();
    auto ice_cld_frac_v      = ice_cld_frac.get_view<Pack**>();
    auto tot_cld_frac_v      = tot_cld_frac.get_view<Pack**>();
    auto ice_cld_frac_4out_v = ice_cld_frac_4out.get_view<Pack**>();
    auto tot_cld_frac_4out_v = tot_cld_frac_4out.get_view<Pack**>();

    CldFractionFunc::main(m_num_cols,m_num_levs,m_icecloud_threshold,m_icecloud_for_analysis_threshold,
      qi_v,liq_cld_frac_v,ice_cld_frac_v,tot_cld_frac_v,ice_cld_frac_4out_v,tot_cld_frac_4out_v);
  }
```

A few observations:

- `m_py_module.has_value()` is a good way to check if the `std::any` object is storing anything or it's empty.
  If empty, it means that the user did not specify the `py_module_name` input option. In this case, we interpret
  this as "proceed with the C++ implementation", but of course, another process may only offer a python
  implementation, in which case it would make sense to error out if the check fails.
- the namespace alias `py = pybind11` was used in this implementation. This is a common (and sometimes
  recommended) practice.
- when casting to pybind11 data structures, we used const references, for both inputs and outputs. The reason
  for using a reference is to avoid copy construction of pybind11 structures (even though they are usually
  lightweight). The const qualifier is not really important, as python has no corresponding concept, and the
  code would have been perfectly fine (and working the same way) without `const`.
- when passing host arrays to python, keep in mind that EAMxx only requires that device views be kept up to
  date by atmosphere processes. Hence, you must take care of syncing to host all inputs before calling the
  python interfaces, as well as syncing to device the outputs upon return.

The python implemenetation of the CldFraction process is provided in `cld_fraction.py`, in the same folder
as the process interface.
