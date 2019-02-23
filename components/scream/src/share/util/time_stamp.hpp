#ifndef SCREAM_TIME_STAMP_HPP
#define SCREAM_TIME_STAMP_HPP

namespace scream {
namespace util {

// Micro-struct, to hold a time stamp
// Note: this is NOT to store the current OS time, but rather
//       the time of the simulation.
class TimeStamp {
public:

  // Default constructor/copy is good enough
  TimeStamp() = default;
  TimeStamp(const TimeStamp&) = default;
  TimeStamp& operator= (const TimeStamp&) = default;

  friend bool operator== (const TimeStamp& ts1, const TimeStamp& ts2);
  friend bool operator<  (const TimeStamp& ts1, const TimeStamp& ts2);

protected:

  // We prevent modification, so one cannot mess with time stamps.
  // TODO: you will need to grant someone friend's access to this class,
  //       cause you WILL need to modify time stamps. When you have decided
  //       what class/function should have this privilege, come back and
  //       grant them friendship here.
  void set_time (const int year, const int day, const int second);  

  int m_yy;   // Year
  int m_dd;   // Day (of the year)
  int m_ss;   // Second (of the day)
};

bool operator== (const TimeStamp& ts1, const TimeStamp& ts2);
bool operator<  (const TimeStamp& ts1, const TimeStamp& ts2);

} // namespace util

} // namespace scream

#endif // SCREAM_TIME_STAMP_HPP
