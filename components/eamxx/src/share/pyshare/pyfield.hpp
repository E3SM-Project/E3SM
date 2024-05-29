#include "share/field/field.hpp"
#include "share/field/field_utils.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace scream {

struct PyField {
  Field f;

  // Create field and allocate memory
  PyField(const std::string& n,
          const std::vector<int>& d)
  {
    std::vector<FieldTag> t (d.size(),FieldTag::Component);
    FieldIdentifier fid(n,FieldLayout(t,d),ekat::units::Units::invalid(),"");
    f = Field(fid);
    f.allocate_view();
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
    switch (rank) {
      case 1:
      {
        Field::view_dev_t<double*> v(arr.mutable_data(0),d[0]);
        f = Field(fid,v);
        break;
      }
      case 2:
      {
        Field::view_dev_t<double**> v(arr.mutable_data(0,0),d[0],d[1]);
        f = Field(fid,v);
        v(0,0) = -1;
        break;
      }
      default:
        EKAT_ERROR_MSG ("AAARGH!\n");
    }
  }

  void print() const {
    print_field_hyperslab(f);
  }

  void cleanup () {
    f = Field();
  }
};

} // namespace scream
