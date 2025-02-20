#ifndef SCREAM_TIME_STAMP_HPP
#define SCREAM_TIME_STAMP_HPP

#include <string>
#include <vector>
#include <cstdint>
#include <limits>

namespace scream {
namespace util {

// Micro-struct, to hold a simulation time stamp
class TimeStamp {
public:

  TimeStamp();
  TimeStamp(const std::vector<int>& date,
            const std::vector<int>& time,
            const int num_steps = 0);
  TimeStamp(const int yy, const int mm, int dd,
            const int h, const int min, const int sec,
            const int num_steps = 0);
  TimeStamp(const TimeStamp&) = default;

  // === Query methods === //

  const std::vector<int>& get_date () const { return m_date; }
  const std::vector<int>& get_time () const { return m_time; }
  int get_year    () const { return m_date[0]; }
  int get_month   () const { return m_date[1]; }
  int get_day     () const { return m_date[2]; }
  int get_hours   () const { return m_time[0]; }
  int get_minutes () const { return m_time[1]; }
  int get_seconds () const { return m_time[2]; }
  int get_num_steps () const { return m_num_steps; }

  bool is_valid    () const;

  int sec_of_day () const { return m_time[0]*3600 + m_time[1]*60 + m_time[2]; }
  std::int64_t seconds_from (const TimeStamp& ts) const;
  double days_from (const TimeStamp& ts) const;

  std::string to_string () const;
  std::string get_date_string () const;
  std::string get_time_string () const;
  double frac_of_year_in_days () const;

  int days_in_curr_month () const;
  int days_in_curr_year () const;

  TimeStamp curr_month_beg () const;

  // === Update method(s) === //

  // Set the counter for the number of steps.
  void set_num_steps (const int num_steps) { m_num_steps  = num_steps; }

  TimeStamp& operator= (const TimeStamp&) = default;

  // This method checks that time shifts forward (i.e. that seconds is positive)
  TimeStamp& operator+= (const double seconds);

  // Clones the stamps and sets num steps to given value. If -1, clones num steps too
  TimeStamp clone (const int num_steps);

protected:

  std::vector<int> m_date;  // [year, month, day]
  std::vector<int> m_time;  // [hour, min, sec]

  int m_num_steps = std::numeric_limits<int>::lowest(); // Number of steps since simulation started
};

// Overload operators for TimeStamp
bool operator== (const TimeStamp& ts1, const TimeStamp& ts2);
bool operator<  (const TimeStamp& ts1, const TimeStamp& ts2);
bool operator<= (const TimeStamp& ts1, const TimeStamp& ts2);
TimeStamp operator+ (const TimeStamp& ts, const int dt);

// Difference (in seconds) between two timestamps
std::int64_t operator- (const TimeStamp& ts1, const TimeStamp& ts2);

// Rewind time by given number of seconds
TimeStamp operator- (const TimeStamp& ts, const int dt);

// If input string is not of the format YYYY-MM-DD-XXXXX, returns an invalid time stamp
TimeStamp str_to_time_stamp (const std::string& s);

// An enum describing two ways to look at timestamps:
//  - Linear: treat them as part of a 1d line
//  - YearlyPeriodic: treat them as part of a yearly periodic orbit
// This is used in the TimeInterval class below to correctly handle time stamps differences
enum class TimeLine {
  YearlyPeriodic,
  Linear
};

/*
 * Small struct to deal with time intervals
 *
 * The struct simply contains timestamps for [begin,end] interval,
 * and allows two things: compute the interval length (in days), and check if
 * a timestamp lies within the interval.
 *
 * When the TimeLine arg to the ctor is YearlyPeriodic, the year part of beg/end
 * time points is ignored. In this case, the length of the time interval is bound
 * to be in the interval [0,365] (in non-leap years)
 */
struct TimeInterval {

  TimeStamp beg;
  TimeStamp end;
  TimeLine  timeline = TimeLine::Linear;
  double length  = -1; // the interval length

  TimeInterval () = default;
  TimeInterval (const util::TimeStamp& b, const util::TimeStamp& e, TimeLine tl, bool do_compute_length = true)
   : beg (b), end (e), timeline (tl)
  {
    if (do_compute_length)
      compute_length ();
  }

  bool contains (const util::TimeStamp& t) const {
    if (timeline==TimeLine::Linear) {
      // Compare the full time stamps
      return beg<=t and t<=end;
    } else {
      // Compare the fraction of year for beg/end and t.
      // Pay extra attention to the case where new year's eve
      // is in [bec,end]
      auto t_frac = t.frac_of_year_in_days();
      auto end_frac = end.frac_of_year_in_days();
      auto beg_frac = beg.frac_of_year_in_days();
      bool across_nye = beg.get_month()>end.get_month();
      if (not across_nye) {
        return beg_frac<=t_frac and t_frac<=end_frac;
      } else {
        // We are either PAST beg or BEFORE end (but not both)
        return t_frac>=beg_frac or t_frac<=end_frac;
      }
    }
  }

  void compute_length () {
    if (timeline==TimeLine::Linear) {
      length = end.days_from(beg);
    } else {
      bool across_nye = beg.get_month()>end.get_month();
      auto frac_beg = beg.frac_of_year_in_days();
      auto frac_end = end.frac_of_year_in_days();
      if (across_nye) {
        double year = end.days_in_curr_year();
        length = frac_end + (year - frac_beg);
      } else {
        length = frac_end - frac_beg;
      }
    }
  }

  // Advance the interval, so that it now starts from the old end
  void advance(const TimeStamp& new_end) {
    beg = end;
    end = new_end;
    compute_length();
  }
};

} // namespace util

} // namespace scream

#endif // SCREAM_TIME_STAMP_HPP
