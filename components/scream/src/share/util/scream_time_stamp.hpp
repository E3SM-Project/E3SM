#ifndef SCREAM_TIME_STAMP_HPP
#define SCREAM_TIME_STAMP_HPP

#include <string>

namespace scream {
namespace util {

// Micro-struct, to hold a simulation time stamp
class TimeStamp {
public:

  TimeStamp();
  // TimeStamp(const int yy, const int dd, const double ss);
  TimeStamp(const int yy, const int mm, int dd, const double ss);
  TimeStamp(const TimeStamp&) = default;

  // === Query methods === //

  int    get_years   () const { return m_yy; }
  int    get_months  () const { return m_mm; }
  int    get_days    () const { return m_dd; }
  double get_seconds () const { return m_ss; }
  bool   is_valid    () const;
  double  test_double () const; // Dummy double precision representation of the timestamp used for operator definitions.

  std::string to_string () const;
  double get_julian_day () const;
  double get_julian_day (const int yy, const int mm, const int dd, const double ss) const;

  // === Update method(s) === //

  TimeStamp& operator= (const TimeStamp&) = default;

  // This method checks that time shifts forward (i.e. that seconds is positive)
  TimeStamp& operator+= (const double seconds);

protected:

  int     m_yy;       // Year
  int     m_mm;       // Month
  int     m_dd;       // Day
  double  m_ss;       // Second (of the day)
  double  m_jd;       // Julian day 1.xx - 365.xx
};

bool operator== (const TimeStamp& ts1, const TimeStamp& ts2);
bool operator<  (const TimeStamp& ts1, const TimeStamp& ts2);
bool operator<= (const TimeStamp& ts1, const TimeStamp& ts2);
TimeStamp operator+ (const TimeStamp& ts, const double dt);
double operator- (const TimeStamp& ts, const TimeStamp& dt);

// Define here instead of inside the class, so we can call op==
inline bool TimeStamp::is_valid () const {
  return !(*this==TimeStamp());
}

// Convert timestamp to a double for operator testing using a simple monotonic function.
inline double TimeStamp::test_double () const {
  return (m_yy+1)*1e4 + (m_mm+1)*1e2 + (m_dd+1) + m_ss/86400;
}

} // namespace util

} // namespace scream

#endif // SCREAM_TIME_STAMP_HPP
