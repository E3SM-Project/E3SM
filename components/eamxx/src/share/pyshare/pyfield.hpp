#ifndef PYFIELD_HPP
#define PYFIELD_HPP

#include "share/field/field.hpp"
#include "share/field/field_utils.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace scream {

struct PyField {
  Field f;

  // Create empty field
  PyField () = default;

  PyField(const FieldIdentifier& fid)
  {
    create(fid);
  }

  // Create field and allocate memory
  PyField(const std::string& n,
          const std::vector<int>& d)
  {
    std::vector<FieldTag> t (d.size(),FieldTag::Component);
    FieldIdentifier fid(n,FieldLayout(t,d),ekat::units::Units::invalid(),"");
    create(fid);
  }

  // Create field from pre-existing python nd array
  PyField(const std::string& n,
          pybind11::array_t<double> arr)
  {
    int rank = arr.ndim();
    std::vector<int> d(rank,-1);
    for (int n=0; n<rank; ++n) {
      d[n] = arr.shape(n);
    }
    std::vector<FieldTag> t (rank,FieldTag::Component);
    FieldIdentifier fid(n,FieldLayout(t,d),ekat::units::Units::invalid(),"");
    create(fid,arr);
  }

  PyField(const FieldIdentifier& fid,
          pybind11::array_t<double> arr)
  {
    create(fid,arr);
  }

  void print() const {
    print_field_hyperslab(f);
  }

  void cleanup () {
    f = Field();
  }

private:
  void create (const FieldIdentifier& fid)
  {
    f = Field(fid);
    f.allocate_view();
  }

  void create (const FieldIdentifier& fid,
               pybind11::array_t<double> arr)
  {
    int rank = arr.ndim();
    EKAT_REQUIRE_MSG (rank==fid.get_layout().rank(),
        "Error! Rank mismatch between input FieldIdentifier and pybind11::array_t.\n"
        "  - field name: " + fid.name() + "\n"
        "  - identifier rank: " + std::to_string(fid.get_layout().rank()) + "\n"
        "  - array_t rank   : " + std::to_string(rank) + "\n");
    switch (rank) {
      case 1:
      {
        Field::view_dev_t<double*> v(arr.mutable_data(0),arr.shape(0));
        f = Field(fid,v);
        break;
      }
      case 2:
      {
        Field::view_dev_t<double**> v(arr.mutable_data(0,0),arr.shape(0),arr.shape(1));
        f = Field(fid,v);
        v(0,0) = -1;
        break;
      }
      default:
        EKAT_ERROR_MSG ("AAARGH!\n");
    }
  }
};

} // namespace scream

#endif // PYFIELD_HPP
