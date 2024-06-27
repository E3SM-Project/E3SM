#include "share/io/scream_io_utils.hpp"

#include "share/io/scream_scorpio_interface.hpp"
#include "share/util/scream_utils.hpp"

#include <fstream>

namespace scream {

std::string find_filename_in_rpointer (
    const std::string& filename_prefix,
    const bool model_restart,
    const ekat::Comm& comm,
    const util::TimeStamp& run_t0)
{
  std::string filename;
  bool found = false;
  std::string content;
  std::string suffix = model_restart ? ".r." : ".rhist.";
  if (comm.am_i_root()) {
    std::ifstream rpointer_file;
    std::string line;
    rpointer_file.open("rpointer.atm");

    // If the timestamp is in the filename, then the filename ends with "S.nc",
    // with S being the string representation of the timestamp
    auto ts_len = run_t0.to_string().size();
    auto extract_ts = [&] (const std::string& line) -> util::TimeStamp {
      auto min_size = ts_len+3;
      if (line.size()>=min_size) {
        auto ts_str = line.substr(line.size()-min_size,ts_len);
        auto ts = util::str_to_time_stamp(ts_str);
        return ts;
      } else {
        return util::TimeStamp();
      }
    };

    while ((rpointer_file >> line) and not found) {
      content += line + "\n";

      found = line.find(filename_prefix+suffix) != std::string::npos &&
              extract_ts(line)==run_t0;
      filename = line;
    }
  }

  int ifound = int(found);
  comm.broadcast(&ifound,1,0);
  found = bool(ifound);

  if (not found) {
    broadcast_string(content,comm,comm.root_rank());

    // If the history restart file is not found, it must be because the last
    // model restart step coincided with a model output step, in which case
    // a restart history file is not written.
    // If that's the case, *disable* output restart, by setting
    //   'Restart'->'Perform Restart' = false
    // in the input parameter list
    EKAT_ERROR_MSG (
        "Error! Restart requested, but no restart file found in 'rpointer.atm'.\n"
        "   restart filename prefix: " + filename_prefix + "\n"
        "   restart file type: " + std::string(model_restart ? "model restart" : "history restart") + "\n"
        "   run t0           : " + run_t0.to_string() + "\n"
        "   rpointer content:\n" + content);
  }

  // Have the root rank communicate the nc filename
  broadcast_string(filename,comm,comm.root_rank());

  return filename;
}

void write_timestamp (const std::string& filename, const std::string& ts_name,
                      const util::TimeStamp& ts, const bool write_nsteps)
{
  scorpio::set_attribute(filename,"GLOBAL",ts_name,ts.to_string());
  if (write_nsteps) {
    scorpio::set_attribute(filename,"GLOBAL",ts_name+"_nsteps",ts.get_num_steps());
  }
}

util::TimeStamp read_timestamp (const std::string& filename,
                                const std::string& ts_name,
                                const bool read_nsteps)
{
  auto ts = util::str_to_time_stamp(scorpio::get_attribute<std::string>(filename,"GLOBAL",ts_name));
  if (read_nsteps and scorpio::has_attribute(filename,"GLOBAL",ts_name+"_nsteps")) {
    ts.set_num_steps(scorpio::get_attribute<int>(filename,"GLOBAL",ts_name+"_nsteps"));
  }
  return ts;
}

} // namespace scream
