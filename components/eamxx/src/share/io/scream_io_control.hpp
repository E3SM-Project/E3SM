#ifndef SCREAM_IO_CONTROL_HPP
#define SCREAM_IO_CONTROL_HPP

#include "share/util/scream_time_stamp.hpp"

#include <ekat/ekat_assert.hpp>

#include <string>

namespace scream
{

// Mini struct to hold IO frequency info
struct IOControl {

  // If frequency_units is not "none" or "never", frequency *must* be set to a positive number
  int frequency = -1;
  std::string frequency_units = "none";

  int nsamples_since_last_write = 0;  // Needed when updating output data, such as with the OAT::Average flag

  util::TimeStamp next_write_ts;
  util::TimeStamp last_write_ts;

  bool output_enabled () const {
    return frequency_units!="none" && frequency_units!="never";
  }

  bool is_write_step (const util::TimeStamp& ts) const {
    if (not output_enabled()) return false;
    return frequency_units=="nsteps" ? ts.get_num_steps()==next_write_ts.get_num_steps()
                                     : ts==next_write_ts;
  }

  // Computes next_write_ts from frequency and last_write_ts
  void compute_next_write_ts () {
    EKAT_REQUIRE_MSG (last_write_ts.is_valid(),
        "Error! Cannot compute next_write_ts, since last_write_ts was never set.\n");
    if (frequency_units=="nsteps") {
      // This avoids having an invalid date/time in the above check next time this fcn runs
      next_write_ts = last_write_ts;
      next_write_ts.set_num_steps(last_write_ts.get_num_steps()+frequency);
    } else if (frequency_units=="nsecs") {
      next_write_ts = last_write_ts;
      next_write_ts += frequency;
    } else if (frequency_units=="nmins") {
      next_write_ts = last_write_ts;
      next_write_ts += frequency*60;
    } else if (frequency_units=="nhours") {
      next_write_ts = last_write_ts;
      next_write_ts += frequency*3600;
    } else if (frequency_units=="ndays") {
      next_write_ts = last_write_ts;
      next_write_ts += frequency*86400;
    } else if (frequency_units=="nmonths" or frequency_units=="nyears") {
      auto date = last_write_ts.get_date();
      if (frequency_units=="nmonths") {
        int temp = date[1] + frequency - 1;
        date[1] = temp % 12 + 1;
        date[0] += temp / 12;
      } else {
        date[0] += frequency;
      }

      // Fix day, in case we moved to a month/year where current days. E.g., if last_write
      // was on Mar 31st, and units='nmonths', next write is on Apr 30th. HOWEVER, this
      // means we will *always* write on the 30th of each month after then, since we have
      // no memory of the fact that we were writing on the 31st before.
      date[2] = std::min(date[2],util::days_in_month(date[0],date[1]));
      next_write_ts = util::TimeStamp(date,last_write_ts.get_time());
    } else {
      EKAT_ERROR_MSG ("Error! Unrecognized/unsupported frequency unit '" + frequency_units + "'\n");
    }
  }
};

} // namespace scream
#endif // SCREAM_IO_CONTROL_HPP
