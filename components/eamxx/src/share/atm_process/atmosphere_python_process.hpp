#ifndef EAMXX_ATMOSPHERE_PYTHON_PROCESS_HPP
#define EAMXX_ATMOSPHERE_PYTHON_PROCESS_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "share/eamxx_pysession.hpp"

#include <nanobind/nanobind.h>

namespace scream
{

template<typename PyFramework>
class AtmospherePythonProcess : public AtmosphereProcess
{
public:
  template<typename T>
  using strmap_t = std::map<std::string,T>;
  template<typename Framework>
  using nb_field = nanobind::ndarray<Framework>

  using numpy = nanobind::numpy;
  using torch = nanobind::torch;
  using cupy  = nanobind::cupy;

  static_assert(std::is_same_v<PyFramework,numpy> or
                std::is_same_v<PyFramework,torch> or
                std::is_same_v<PyFramework,cupy>,
    "Pythong framework not supported.\n"
    "  - supported options: nb::numpy, nb::torch, nb::cupy.\n");


  AtmospherePythonProcess (const ekat::Comm& comm, const ekat::ParameterList& params)
   : AtmosphereProcess(comm,params)
  {
    PySession::get().initialize();
    const auto& py_module_name = m_params.get<std::string>("py_module_name");

    // Create a py-compatible string
    PyObject* pName = PyUnicode_FromString(py_module_name.c_str());
    if (pName == nullptr) {
      // Handle error
      PyErr_Print();
      EKAT_ERROR_MSG ("Error! Could not convert std::string to py-compatible string.\n");
    }

    // Import the module
    m_py_module = PyImport_Import(pName);
    Py_DECREF(pName); // Decrement reference count for pName

    if (m_py_module == nullptr) {
      // Handle error
      PyErr_Print();
      EKAT_ERROR_MSG ("Error! Could not import module '" + py_module_name + "'.\n");
    }
  }

  ~AtmospherePythonProcess () {
    if (m_py_module != nullptr) {
      Py_DECREF(m_py_module);
    }
    PySession::get().finalize();
  }

protected:

  void run_py_module (const double dt) {
  }

  void set_required_field_impl (const Field& f) override {
    auto& py_f = m_py_fields_in[f.name()];
    create_py_field<true>(py_f,f);
  }
  void set_computed_field_impl (const Field& f) override {
    auto& py_f = m_py_fields_out[f.name()];
    create_py_field<false>(py_f,f);
  }

  template<bool ReadOnly,typename Framework>
  void create_py_field (nb_field<Framework>& py_f, const Field& f)
  {
    using void_ptr = typename std::conditional<ReadOnly,const void*,void*>::type;

    const auto& fh  = f.get_header();
    const auto& fid = fh.get_identifier();

    // Get array shape and strides (nb requires size_t for shape and int64_t for strides)
    // NOTE: since the field may be padded or a subfield, so strides do not necessarily
    //       match the field dims. Also, the strides must be grabbed from the
    //       actual view, since the layout doesn't know them.
    std::vector<int64_t> strides;
    std::vector<size_t>  shape;
    for (auto d : fid.get_layout().dims()) {
      shape.push_back(d);
    }

    // Also the ptr must be grabbed from the view, to make sure it points to the 1st
    // element of the field. In case f is a subfield, this is not the same as the
    // internal data pointer of the field.
    nb::dlpack::dtype dt;
    void_ptr* data;
    switch (fid.data_type()) {
      case DataType::IntType:
        dt = nb::dtype<int>();
        get_ptr_and_strides<int>(strides,data,f);
        break;
      case DataType::FloatType:
        dt = nb::dtype<float>();
        get_ptr_and_strides<float>(strides,data,f);
        break;
      case DataType::DoubleType:
        dt = nb::dtype<double>();
        get_ptr_and_strides<double>(strides,data,f);
        break;
      default:
        EKAT_ERROR_MSG ("Unrecognized/unsupported data type.\n");
    }

    // NOTE: you MUST set the parent handle, or else you won't have view semantic
    auto this_obj = nb::cast(this);
    py_f = nb_field<Framework>(data, shape_t, shape, nb::handle(this_obj), strides.data(), dt);
  }

  template<typename T, typename Ptr>
  void get_ptr_and_strides (std::vector<int64_t>& strides, Ptr& p, const Field& f)
  {
    constexpr bool is_const = std::is_const<typename std::remove_pointer<T>::type>::value;
    using value_type = typename std::conditional<is_const,const T, T>::type;

    const auto& fl  = f.get_header().get_identifier().get_layout();
    strides.resize(fl.rank());
    switch (fl.rank()) {
      case 1:
        {
          auto v = f.get_view<value_type*>();
          strides[0] = v.stride_0();
          p = v.data();
          break;
        }
      case 2:
        {
          auto v = f.get_view<value_type**>();
          strides[0] = v.stride_0();
          strides[1] = v.stride_1();
          p = v.data();
          break;
        }
        break;
      case 3:
        {
          auto v = f.get_view<value_type***>();
          strides[0] = v.stride_0();
          strides[1] = v.stride_1();
          strides[2] = v.stride_2();
          p = v.data();
          break;
        }
        break;
      default:
        EKAT_ERROR_MSG ("Error! Rank " + std::to_string(fl.rank()) + " not supported.\n");
    }
  }

  strmap_t<nb_field<PyFramework>> m_py_fields_in;
  strmap_t<nb_field<PyFramework>> m_py_fields_out;

  PyObject* m_py_module = nullptr;
};

} // namespace scream

#endif // EAMXX_ATMOSPHERE_PYTHON_PROCESS_HPP
