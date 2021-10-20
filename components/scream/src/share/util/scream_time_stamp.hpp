#ifndef SCREAM_TIME_STAMP_HPP
#define SCREAM_TIME_STAMP_HPP

#include <string>
#include <vector>

namespace scream {
namespace util {

// Micro-struct, to hold a simulation time stamp
class TimeStamp {
public:

  TimeStamp();
  TimeStamp(const std::vector<int>& date,
            const std::vector<int>& time);
  TimeStamp(const int yy, const int mm, int dd,
            const int h, const int min, const int sec);
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

  bool is_valid    () const;

  int sec_of_day () const { return m_time[0]*3600 + m_time[1]*60 + m_time[2]; }
  int seconds_from (const TimeStamp& ts) const;

  std::string to_string () const;
  std::string get_date_string () const;
  std::string get_time_string () const;
  double frac_of_year_in_days () const;

  // === Update method(s) === //

  TimeStamp& operator= (const TimeStamp&) = default;

  // This method checks that time shifts forward (i.e. that seconds is positive)
  TimeStamp& operator+= (const int seconds);

protected:

  std::vector<int> m_date;  // [year, month, day]
  std::vector<int> m_time;  // [hour, min, sec]
};

// Overload operators for TimeStamp
bool operator== (const TimeStamp& ts1, const TimeStamp& ts2);
bool operator<  (const TimeStamp& ts1, const TimeStamp& ts2);
bool operator<= (const TimeStamp& ts1, const TimeStamp& ts2);
TimeStamp operator+ (const TimeStamp& ts, const int dt);

// Difference (in seconds) between two timestamps
long long operator- (const TimeStamp& ts1, const TimeStamp& ts2);

// Time-related free-functions
int days_in_month (const int year, const int month);
bool is_leap_year (const int year);

} // namespace util

} // namespace scream

#endif // SCREAM_TIME_STAMP_HPP
