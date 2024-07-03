#ifndef PYATMPROC_HPP
#define PYATMPROC_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/io/scream_output_manager.hpp"

#include "pygrid.hpp"
#include "pyfield.hpp"
#include "pyparamlist.hpp"
#include "pyeamxx.hpp"

#include <ekat/io/ekat_yaml.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace scream {

struct PyAtmProc {
  std::shared_ptr<AtmosphereProcess> ap;
  std::map<std::string,PyField> fields;
  util::TimeStamp t0;
  util::TimeStamp time;
  ATMBufferManager buffer;

  std::shared_ptr<OutputManager> output_mgr;

  PyAtmProc (const pybind11::dict& d, const std::string& name)
  {
    PyParamList params(d,name);

    // Get the comm
    const auto& comm = PySession::get().comm;

    // Create the atm proc
    auto& apf = AtmosphereProcessFactory::instance();
    const auto& ap_type = params.pl.isParameter("Type")
                        ? params.pl.get<std::string>("Type")
                        : params.pl.name();
    ap = apf.create(ap_type,comm,params.pl);

    // Create the fields
    auto gm = PySession::get().gm;
    ap->set_grids(gm);
    create_fields();
  }

  // I don't think virtual is needed, but just in case
  virtual ~PyAtmProc () = default;

  PyParamList get_params () const {
    PyParamList pypl(ap->get_params());
    return pypl;
  }

  void create_fields () {
    // Create  fields that are input/output to the atm proc
    for (const auto& req : ap->get_required_field_requests()) {
      const auto& fn = req.fid.name();
      auto it_bool = fields.emplace(fn,PyField(req.fid,req.pack_size));
      ap->set_required_field(it_bool.first->second.f.get_const());
    }
    for (const auto& req : ap->get_computed_field_requests()) {
      const auto& fn = req.fid.name();
      auto it_bool = fields.emplace(fn,PyField(req.fid,req.pack_size));
      ap->set_computed_field(it_bool.first->second.f);
    }
  }

  PyField get_field(const std::string& name) {
    auto it = fields.find(name);

    auto print_key = [](auto it) {
      return it.first;
    };
    EKAT_REQUIRE_MSG (it!=fields.end(),
        "Error! Field not found among this atm proc fields.\n"
        "  - atm proc name: " + ap->name() + "\n"
        "  - field name: " + name + "\n"
        "  - atm proc fields: " + ekat::join(fields,print_key,",") + "\n");

    return it->second;
  }

  // If running as part of a process group, call the second function, after
  // manually creating/setting the fields
  void initialize (const std::string& t0_str) {
    int nbytes = ap->requested_buffer_size_in_bytes ();
    buffer.request_bytes(nbytes);
    buffer.allocate();
    ap->init_buffers(buffer);

    time = t0 = util::str_to_time_stamp(t0_str);
    for (auto it : fields) {
      auto& f = it.second.f;
      if (ap->has_required_field(f.get_header().get_identifier())) {
        f.get_header().get_tracking().update_time_stamp(t0);
      }
    }

    ap->initialize(t0,RunType::Initial);
  }

  pybind11::list read_ic (const std::string& ic_filename) {
    // Get input fields, and read them from file (if present).
    // If field is not in the IC, user is responsible for setting
    // it to an initial value
    std::vector<std::string> missing;
    std::vector<Field> ic_fields;
    scorpio::register_file(ic_filename,scorpio::Read);
    for (auto it : fields) {
      auto& f = it.second.f;
      if (ap->has_required_field(f.get_header().get_identifier())) {
        if (scorpio::has_var(ic_filename,f.name())) {
          ic_fields.push_back(f);
        } else {
          missing.push_back(f.name());
        }
      }
    }
    if (ic_fields.size()>0) {
      const auto& gn = ic_fields[0].get_header().get_identifier().get_grid_name();
      auto gm = PySession::get().gm;
      auto grid = gm->get_grid(gn);
      AtmosphereInput reader (ic_filename,grid,ic_fields,true);
      reader.read_variables();
    }
    scorpio::release_file(ic_filename);

    return pybind11::cast(missing);
  }

  void setup_output (const std::string& yaml_file) {
    auto comm = PySession::get().comm;

    // Load output params
    auto params = ekat::parse_yaml_file(yaml_file);

    // Stuff all fields in a field manager
    auto gm = PySession::get().gm;
    std::map<std::string,std::shared_ptr<FieldManager>> fms;
    for (auto it : gm->get_repo()) {
      fms[it.first] = std::make_shared<FieldManager>(it.second);
      fms[it.first]->registration_begins();
      fms[it.first]->registration_ends();
    }
    for (auto it : fields) {
      const auto& gn = it.second.f.get_header().get_identifier().get_grid_name();
      fms.at(gn)->add_field(it.second.f);
    }

    // Create/setup the output mgr
    output_mgr = std::make_shared<OutputManager>();
    output_mgr->setup(comm,params,fms,gm,t0,t0,false);
    output_mgr->set_logger(ap->get_logger());
  }

  void run (double dt) {
    ap->run(dt);
    time += dt;
    if (output_mgr) {
      output_mgr->run(time);
    }
  }
};

// Register type in the py module
inline void pybind_pyatmproc(pybind11::module& m)
{
  pybind11::class_<PyAtmProc>(m,"AtmProc")
    .def(pybind11::init<const pybind11::dict&,const std::string&>())
    .def("get_field",&PyAtmProc::get_field)
    .def("initialize",&PyAtmProc::initialize)
    .def("get_params",&PyAtmProc::get_params)
    .def("setup_output",&PyAtmProc::setup_output)
    .def("run",&PyAtmProc::run)
    .def("read_ic",&PyAtmProc::read_ic);
}
} // namespace scream

#endif // PYATMPROC_HPP
