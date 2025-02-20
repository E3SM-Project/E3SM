#ifndef PYPARAMLIST_HPP
#define PYPARAMLIST_HPP

#include <ekat/ekat_parameter_list.hpp>

#include <nanobind/nanobind.h>
#include <nanobind/stl/list.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <functional>

namespace nb = nanobind;

namespace scream {

struct PyParamList {
  ekat::ParameterList pl;
  std::reference_wrapper<ekat::ParameterList> pl_ref;

  PyParamList (ekat::ParameterList& src)
   : pl_ref(src)
  {}

  PyParamList(const nb::dict& d)
   : PyParamList(d,"")
  {}

  PyParamList(const nb::dict& d, const std::string& name)
   : pl(name)
   , pl_ref(pl)
  {
    parse_dict(d,pl);
  }

  PyParamList sublist (const std::string& name) {
    PyParamList spl(pl.sublist(name));
    return spl;
  }

  bool get_bool (const std::string& name) const {
    return pl_ref.get().get<bool>(name);
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

  std::vector<int> get_int_vec (const std::string& name) const {
    return pl_ref.get().get<std::vector<int>>(name);
  }
  std::vector<double> get_dbl_vec (const std::string& name) const {
    return pl_ref.get().get<std::vector<double>>(name);
  }
  std::vector<std::string> get_str_vec (const std::string& name) const {
    return pl_ref.get().get<std::vector<std::string>>(name);
  }

  template<typename T>
  void set (const std::string& name, T val) {
    pl_ref.get().set(name,val);
  }

  void print () {
    pl_ref.get().print();
  }

private:

  void parse_dict(const nb::dict& d, ekat::ParameterList& p) {
    for (auto item : d) {
      auto key = nb::cast<std::string>(item.first);
      if (nb::isinstance<nb::str>(item.second)) {
        auto pystr = nb::str(item.second);
        p.set<std::string>(key,nb::cast<std::string>(pystr));
      } else if (nb::isinstance<nb::bool_>(item.second)) {
        auto pyint = nb::cast<nb::bool_>(item.second);
        p.set(key,nb::cast<bool>(pyint));
      } else if (nb::isinstance<nb::int_>(item.second)) {
        auto pyint = nb::cast<nb::int_>(item.second);
        p.set(key,nb::cast<int>(pyint));
      } else if (nb::isinstance<nb::float_>(item.second)) {
        auto pydouble = nb::cast<nb::float_>(item.second);
        p.set(key,nb::cast<double>(pydouble));
      } else if (nb::isinstance<nb::list>(item.second)) {
        auto pylist = nb::cast<nb::list>(item.second);
        parse_list(pylist,p,key);
      } else if (nb::isinstance<nb::dict>(item.second)) {
        auto pydict = nb::cast<nb::dict>(item.second);
        parse_dict(pydict,p.sublist(key));
      } else {
        EKAT_ERROR_MSG ("Unsupported/unrecognized dict entry type.\n");
      }
    }
  }

  void parse_list (const nb::list& l, ekat::ParameterList&p, const std::string& key) {
    EKAT_REQUIRE_MSG (nb::len(l)>0,
        "Error! Cannot deduce type for dictionary list entry '" + key + "'\n");
    auto first = l[0];
    bool are_ints = nb::isinstance<nb::int_>(first);
    bool are_floats = nb::isinstance<nb::float_>(first);
    bool are_strings = nb::isinstance<nb::str>(first);
    if (are_ints) {
      parse_list_impl<int,nb::int_>(l,p,key);
    } else if (are_floats) {
      parse_list_impl<double,nb::float_>(l,p,key);
    } else if (are_strings) {
      parse_list_impl<std::string,nb::str>(l,p,key);
    } else {
      EKAT_ERROR_MSG ("Unrecognized/unsupported list entry type.\n");
    }
  }

  template<typename Txx, typename Tpy>
  void parse_list_impl(const nb::list& l, ekat::ParameterList& p, const std::string& key) {
    std::vector<Txx> vals;
    for (auto item : l) {
      EKAT_REQUIRE_MSG (nb::isinstance<Tpy>(item),
          "Error! Inconsistent types in list entries.\n");
      auto item_py = nb::cast<Tpy>(item);
      vals.push_back(nb::cast<Txx>(item_py));
    }
    p.set(key,vals);
  }
};

inline void nb_pyparamlist (nb::module_& m)
{
  // Param list
  nb::class_<PyParamList>(m,"ParameterList")
    .def(nb::init<const nb::dict&>())
    .def(nb::init<const nb::dict&,const std::string&>())
    .def("sublist",&PyParamList::sublist)
    .def("print",&PyParamList::print)
    .def("set",&PyParamList::set<bool>)
    .def("set",&PyParamList::set<int>)
    .def("set",&PyParamList::set<double>)
    .def("set",&PyParamList::set<std::string>)
    .def("set",&PyParamList::set<std::vector<int>>)
    .def("set",&PyParamList::set<std::vector<double>>)
    .def("set",&PyParamList::set<std::vector<std::string>>)
    .def("get_bool",&PyParamList::get_bool)
    .def("get_int",&PyParamList::get_int)
    .def("get_dbl",&PyParamList::get_dbl)
    .def("get_str",&PyParamList::get_str)
    .def("get_int_vec",&PyParamList::get_int_vec)
    .def("get_dbl_vec",&PyParamList::get_dbl_vec)
    .def("get_str_vec",&PyParamList::get_str_vec);
}

} // namespace scream

#endif // PYPARAMLIST_HPP
