#ifndef PYATMPROC_HPP
#define PYATMPROC_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "pygrid.hpp"
#include "pyfield.hpp"

#include <pybind11/pybind11.h>

namespace scream {

struct PyAtmProc {
  std::shared_ptr<AtmosphereProcess> ap;
  std::map<std::string,PyField> fields;

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
    set_fields(create_fields());
    util::TimeStamp ts;
    ap->initialize(ts,RunType::Initial);
  }
  void initialize_without_creating_fields() {
    util::TimeStamp ts;
    ap->initialize(ts,RunType::Initial);
  }
};

} // namespace scream

#endif // PYATMPROC_HPP
