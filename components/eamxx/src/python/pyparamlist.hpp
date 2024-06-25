#ifndef PYPARAMLIST_HPP
#define PYPARAMLIST_HPP

#include <ekat/ekat_parameter_list.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <functional>

namespace scream {

struct PyParamList {
  ekat::ParameterList pl;
  std::reference_wrapper<ekat::ParameterList> pl_ref;

  PyParamList (ekat::ParameterList& src)
   : pl_ref(src)
  {}

  PyParamList(const pybind11::dict& d)
   : PyParamList(d,"")
  {}

  PyParamList(const pybind11::dict& d, const std::string& name)
   : pl(name)
   , pl_ref(pl)
  {
    parse_dict(d,pl);
  }

  PyParamList sublist (const std::string& name) {
    PyParamList spl(pl.sublist(name));
    return spl;
  }

  int get_int (const std::string& name) const {
    return pl_ref.get().get<int>(name);
  }
  double get_dbl (const std::string& name) const {
    return pl_ref.get().get<double>(name);
  }
  std::string get_str (const std::string& name) const {
    return pl_ref.get().get<std::string>(name);
  }

  template<typename T>
  void set (const std::string& name, T val) {
    pl_ref.get().set(name,val);
  }

  void print () {
    pl_ref.get().print();
  }

private:

  void parse_dict(const pybind11::dict& d, ekat::ParameterList& p) {
    for (auto item : d) {
      const std::string key = pybind11::str(item.first);
      if (pybind11::isinstance<pybind11::str>(item.second)) {
        auto pystr = pybind11::str(item.second);
        p.set<std::string>(key,pystr.cast<std::string>());
      } else if (pybind11::isinstance<pybind11::bool_>(item.second)) {
        auto pyint = pybind11::cast<pybind11::bool_>(item.second);
        p.set(key,pyint.cast<bool>());
      } else if (pybind11::isinstance<pybind11::int_>(item.second)) {
        auto pyint = pybind11::cast<pybind11::int_>(item.second);
        p.set(key,pyint.cast<int>());
      } else if (pybind11::isinstance<pybind11::float_>(item.second)) {
        auto pydouble = pybind11::cast<pybind11::float_>(item.second);
        p.set(key,pydouble.cast<double>());
      } else if (pybind11::isinstance<pybind11::list>(item.second)) {
        auto pylist = pybind11::cast<pybind11::list>(item.second);
        parse_list(pylist,p,key);
      } else if (pybind11::isinstance<pybind11::dict>(item.second)) {
        auto pydict = pybind11::cast<pybind11::dict>(item.second);
        parse_dict(pydict,p.sublist(key));
      } else {
        EKAT_ERROR_MSG ("Unsupported/unrecognized dict entry type.\n");
      }
    }
  }

  void parse_list (const pybind11::list& l, ekat::ParameterList&p, const std::string& key) {
    EKAT_REQUIRE_MSG (pybind11::len(l)>0,
        "Error! Cannot deduce type for dictionary list entry '" + key + "'\n");
    auto first = l[0];
    bool are_ints = pybind11::isinstance<pybind11::int_>(first);
    bool are_floats = pybind11::isinstance<pybind11::float_>(first);
    bool are_strings = pybind11::isinstance<pybind11::str>(first);
    if (are_ints) {
      parse_list_impl<int,pybind11::int_>(l,p,key);
    } else if (are_floats) {
      parse_list_impl<double,pybind11::float_>(l,p,key);
    } else if (are_strings) {
      parse_list_impl<std::string,pybind11::str>(l,p,key);
    } else {
      EKAT_ERROR_MSG ("Unrecognized/unsupported list entry type.\n");
    }
  }

  template<typename Txx, typename Tpy>
  void parse_list_impl(const pybind11::list& l, ekat::ParameterList& p, const std::string& key) {
    std::vector<Txx> vals;
    for (auto item : l) {
      EKAT_REQUIRE_MSG (pybind11::isinstance<Tpy>(item),
          "Error! Inconsistent types in list entries.\n");
      auto item_py = pybind11::cast<Tpy>(item);
      vals.push_back(item_py.template cast<Txx>());
    }
    p.set(key,vals);
  }
};

inline void pybind_pyparamlist (pybind11::module& m)
{
  // Param list
  pybind11::class_<PyParamList>(m,"ParameterList")
    .def(pybind11::init<const pybind11::dict&>())
    .def(pybind11::init<const pybind11::dict&,const std::string&>())
    .def("sublist",&PyParamList::sublist)
    .def("print",&PyParamList::print)
    .def("set",&PyParamList::set<int>)
    .def("set",&PyParamList::set<double>)
    .def("set",&PyParamList::set<std::string>)
    .def("get_int",&PyParamList::get_int)
    .def("get_dbl",&PyParamList::get_dbl)
    .def("get_str",&PyParamList::get_str);
}

} // namespace scream

#endif // PYPARAMLIST_HPP
