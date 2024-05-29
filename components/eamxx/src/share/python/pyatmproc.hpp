#ifndef PYATMPROC_HPP
#define PYATMPROC_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "share/io/scorpio_input.hpp"
#include "pygrid.hpp"
#include "pyfield.hpp"

#include <pybind11/pybind11.h>

namespace scream {

struct PyAtmProc {
  std::shared_ptr<AtmosphereProcess> ap;
  PyGrid pygrid;
  std::map<std::string,PyField> fields;

  PyAtmProc (const PyGrid& pyg)
   : pygrid(pyg)
  {
    // Nothing to do here
  }

  std::map<std::string,PyField> create_fields () {
    std::map<std::string,PyField> pyfields;
    for (const auto& req : ap->get_required_field_requests()) {
      if (pyfields.count(req.fid.name())==0) {
        pyfields.emplace(req.fid.name(),PyField(req.fid,req.pack_size));
      }
    }
    for (const auto& req : ap->get_computed_field_requests()) {
      if (pyfields.count(req.fid.name())==0) {
        pyfields.emplace(req.fid.name(),PyField(req.fid,req.pack_size));
      }
    }

    return pyfields;
  }

  void set_fields(const std::map<std::string,PyField>& pyfields) {
    for (const auto& it : pyfields) {
      const auto& f = it.second.f;
      const auto& fid = f.get_header().get_identifier();
      if (ap->has_required_field(fid)) {
        ap->set_required_field(f.get_const());
        fields[fid.name()] = it.second;   
      }
      if (ap->has_computed_field(fid)) {
        ap->set_computed_field(f);
        fields[fid.name()] = it.second;   
      }
    }
  }

  pybind11::array get_arr(const std::string& name) {
    auto it = fields.find(name);

    auto fnames = [&]() {
      std::string s;
      for (auto it : fields) {
        if (s=="") {
          s += it.first;
        } else {
          s += ", " + it.first;
        }
      }
      return s;
    };
    EKAT_REQUIRE_MSG (it!=fields.end(),
        "Error! Field not found in list of P3 fields.\n"
        "  - field name: " + name + "\n"
        "  - p3 fields: " + fnames() + "\n");

    return it->second.get();
  }

  // If running as part of a process group, call the second function, after
  // manually creating/setting the fields
  void initialize () {
    // Create fields
    set_fields(create_fields());

    // TODO: should we allow setting this?
    util::TimeStamp t0({2000,1,1},{0,0,0});
    ap->initialize(t0,RunType::Initial);
  }
  void initialize_without_creating_fields() {
    // TODO: should we allow setting this?
    util::TimeStamp t0({2000,1,1},{0,0,0});
    ap->initialize(t0,RunType::Initial);
  }
  void initialize (const std::string& ic_filename) {
    // Create fields
    set_fields(create_fields());

    // Get input fields, and read them from file
    std::vector<Field> ic_fields;
    for (auto it : fields) {
      const auto& f = it.second.f;
      if (ap->has_required_field(f.get_header().get_identifier())) {
        ic_fields.push_back(f);
      }
    }
    AtmosphereInput reader (ic_filename,pygrid.grid,ic_fields,true);
    reader.read_variables();

    // TODO: should we allow setting this?
    util::TimeStamp t0({2000,1,1},{0,0,0});
    for (auto& f : ic_fields) {
      f.get_header().get_tracking().update_time_stamp(t0);
    }

    ap->initialize(t0,RunType::Initial);
  }
};

} // namespace scream

#endif // PYATMPROC_HPP
