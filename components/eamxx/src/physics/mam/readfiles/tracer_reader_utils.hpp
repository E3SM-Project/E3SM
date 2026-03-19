#ifndef EAMXX_MAM_TRACER_READER_UTILS
#define EAMXX_MAM_TRACER_READER_UTILS

#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/util/eamxx_time_stamp.hpp"
#include "share/algorithm/eamxx_time_interpolation.hpp"

namespace scream::mam_coupling {

// We have a similar version in MAM4xx.
// This version was created because the data view cannot be modified
// inside the parallel_for.
// This struct will be used in init while reading nc files.
// The MAM4xx version will be used instead of parallel_for that loops over cols.
constexpr int MAX_SECTION_NUM_FORCING=4;
struct ForcingHelper {
  // This index is in Fortran format. i.e. starts in 1
  int frc_ndx;
  // does this file have altitude?
  bool file_alt_data;
  // number of sectors per forcing
  int nsectors;
  // data of views
  view_2d fields[MAX_SECTION_NUM_FORCING];
};
inline scream::util::TimeStamp convert_date(const int date) {
  constexpr int ten_thousand = 10000;
  constexpr int one_hundred = 100;

  int year = date / ten_thousand;
  int month = (date - year * ten_thousand) / one_hundred;
  int day = date - year * ten_thousand - month * one_hundred;

  return scream::util::TimeStamp(year, month, day, 0, 0, 0);
}
// FIXME: This function is not implemented in eamxx.
// FIXME: Assumes 365 days/year, 30 days/month;
// NOTE: that this assumption is mainly used for plotting.
// NOTE: This is not a direct port from EAM.
// We only use this routine for chlorine.
inline int compute_number_days_from_zero(const util::TimeStamp &ts) {
  return ts.get_year() * 365 + ts.get_month() * 30 + ts.get_day();
}

inline void create_linoz_chlorine_reader(
    const std::string &linoz_chlorine_file, const util::TimeStamp &model_time,
    const int chlorine_loading_ymd,  // in format YYYYMMDD
    std::vector<Real> &values, std::vector<int> &time_secs) {
  auto time_stamp_beg = convert_date(chlorine_loading_ymd);

  const int offset_time = compute_number_days_from_zero(time_stamp_beg) -
                          compute_number_days_from_zero(model_time);
  scorpio::register_file(linoz_chlorine_file, scorpio::Read);
  const int nlevs_time = scorpio::get_time_len(linoz_chlorine_file);
  for(int itime = 0; itime < nlevs_time; ++itime) {
    int date;
    scorpio::read_var(linoz_chlorine_file, "date", &date, itime);
    if(date >= chlorine_loading_ymd) {
      Real value;
      scorpio::read_var(linoz_chlorine_file, "chlorine_loading", &value, itime);
      values.push_back(value);
      auto time_stamp = convert_date(date);
      time_secs.push_back(compute_number_days_from_zero(time_stamp) -
                          offset_time);
    }
  }  // end itime
  scorpio::release_file(linoz_chlorine_file);
}
inline Real chlorine_loading_advance(const util::TimeStamp &ts,
                                     std::vector<Real> &values,
                                     std::vector<int> &time_secs) {
  const int current_time = compute_number_days_from_zero(ts);
  int index              = 0;
  // update index
  for(int i = 0; i < int(values.size()); i++) {
    if(current_time > time_secs[i]) {
      index = i;
      break;
    }
  }  //

  const Real delt = (current_time - time_secs[index]) /
                    (time_secs[index + 1] - time_secs[index]);
  return values[index] + delt * (values[index + 1] - values[index]);
}

}  // namespace scream::mam_coupling
#endif  // EAMXX_MAM_HELPER_MICRO
