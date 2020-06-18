#ifndef SCREAM_TIME_STAMP_HPP
#define SCREAM_TIME_STAMP_HPP

#include "ekat/scream_types.hpp"

#include <string>

namespace scream {
namespace util {

// Micro-struct, to hold a time stamp
// Note: this is NOT to store the current OS time, but rather
//       the time of the simulation.
class TimeStamp {
public:

  TimeStamp();
  // TimeStamp(const int yy, const int dd, const Real ss);
  TimeStamp(const int yy, const int mm, int dd, const Real ss);
  TimeStamp(const TimeStamp&) = default;

  // === Query methods === //

  int    get_years   () const { return m_yy; }
  int    get_months  () const { return m_mm; }
  int    get_days    () const { return m_dd; }
  Real   get_seconds () const { return m_ss; }
  bool   is_valid    () const;

  std::string to_string () const;

  // === Update method(s) === //

  TimeStamp& operator= (const TimeStamp&) = default;

  // This methods will check that time shifts forward
  TimeStamp& operator+= (const Real seconds);

protected:

  int m_yy;       // Year
  int m_mm;       // Month
  int m_dd;       // Day
  Real m_ss;      // Second (of the day)
};

bool operator== (const TimeStamp& ts1, const TimeStamp& ts2);
bool operator<  (const TimeStamp& ts1, const TimeStamp& ts2);
bool operator<= (const TimeStamp& ts1, const TimeStamp& ts2);
TimeStamp operator+ (const TimeStamp& ts, const Real dt);
Real operator- (const TimeStamp& ts, const TimeStamp& dt);

// Define here instead of inside the class, so we can call op==
inline bool TimeStamp::is_valid () const {
  return !(*this==TimeStamp());
}

} // namespace util

} // namespace scream

#endif // SCREAM_TIME_STAMP_HPP
