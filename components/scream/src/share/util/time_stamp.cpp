#include "time_stamp.hpp"

namespace scream {
namespace util {

void TimeStamp::set_time(const int year, const int day, const int second) {
  m_yy = year;
  m_dd = day;
  m_ss = second;
}

bool operator== (const TimeStamp& ts1, const TimeStamp& ts2) {
  return ts1.m_ss==ts2.m_ss && ts1.m_dd==ts2.m_dd && ts1.m_yy==ts2.m_yy;
}

bool operator<  (const TimeStamp& ts1, const TimeStamp& ts2) {
  if (ts1.m_yy>ts2.m_yy) {
    return false;
  } else if (ts1.m_yy==ts2.m_yy) {
    if (ts1.m_dd>ts2.m_dd) {
      return false;
    } else if (ts1.m_dd==ts2.m_dd) {
      return ts1.m_ss<ts2.m_ss;
    }
  }
  return false;
}

} // namespace util

} // namespace scream

