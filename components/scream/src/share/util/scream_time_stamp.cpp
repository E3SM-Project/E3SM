#include "share/util/scream_time_stamp.hpp"

#include "ekat/ekat_assert.hpp"
#include "share/util/scream_universal_constants.hpp"

#include <numeric>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>

namespace scream {
namespace util {

namespace {

constexpr int nonleap_days [12] = {31,28,31,30,31,30,31,31,30,31,30,31};
constexpr int leap_days    [12] = {31,29,31,30,31,30,31,31,30,31,30,31};

// Utility functions
bool is_leap (const int yy) {
#ifdef EKAT_HAS_LEAP_YEAR
  if (yy%4==0) {
    // Year is divisible by 4 (minimum requirement)
    if (yy%100 != 0) {
      // Not a centennial year => leap.
      return true;
    } else if ((yy/100)%4==0) {
      // Centennial year, but first 2 digids divisible by 4 => leap
      return true;
    }
  }
#else
  (void)yy;
#endif
  return false;
}

int dpm (const int yy, const int mm) {
  auto& arr = is_leap(yy) ? leap_days : nonleap_days;
  return arr[mm];
}

} // anonymous namespace

TimeStamp::TimeStamp()
 : m_date (3,std::numeric_limits<int>::lowest())
 , m_time (3,std::numeric_limits<int>::lowest())
{
  // Nothing to do here
}

TimeStamp::TimeStamp(const std::vector<int>& date,
                     const std::vector<int>& time)
{
  EKAT_REQUIRE_MSG (date.size()==3, "Error! Date should consist of three ints: [year, month, day].\n");
  EKAT_REQUIRE_MSG (time.size()==3, "Error! Time of day should consist of three ints: [hour, min, sec].\n");

  const auto yy   = date[0];
  const auto mm   = date[1];
  const auto dd   = date[2];
  const auto hour = time[0];
  const auto min  = time[1];
  const auto sec  = time[2];

  // Check the days and seconds numbers are non-negative.
  EKAT_REQUIRE_MSG (mm>0   && mm<=12, "Error! Month out of bounds.\n");
  EKAT_REQUIRE_MSG (dd>0   && dd<=dpm(yy,mm), "Error! Day out of bounds.\n");
  EKAT_REQUIRE_MSG (sec>=0  && sec<60, "Error! Seconds out of bounds.\n");
  EKAT_REQUIRE_MSG (min>=0  && min<60, "Error! Minutes out of bounds.\n");
  EKAT_REQUIRE_MSG (hour>=0 && hour<24, "Error! Hours out of bounds.\n");

  // All good, store
  m_date = date;
  m_time = time;
}

TimeStamp::TimeStamp(const int yy, const int mm, const int dd,
                     const int h, const int min, const int sec)
 : TimeStamp({yy,mm,dd},{h,min,sec})
{
  // Nothing to do here
}

std::string TimeStamp::to_string () const {

  std::ostringstream ymdhms;

  ymdhms << std::setw(4) << std::setfill('0') << m_date[0] << "-";
  ymdhms << std::setw(2) << std::setfill('0') << m_date[1] << "-";
  ymdhms << std::setw(2) << std::setfill('0') << m_date[2] << " ";
  ymdhms << std::setw(2) << std::setfill('0') << m_time[0] << ":";
  ymdhms << std::setw(2) << std::setfill('0') << m_time[1] << ":";
  ymdhms << std::setw(2) << std::setfill('0') << m_time[2];

  return ymdhms.str();
}

double TimeStamp::get_julian_day () const {
  return julian_day(get_years(),get_months(),get_days(),sec_of_day());
}

double julian_day (const int yy, const int mm, const int dd, const int sec_of_day) {
  // Initialize Julian Day
  double julianday = (dd-1) + double(sec_of_day) / 86400; // WARNING: avoid integer division
  for (int m=1;m<mm;m++) {
    julianday += dpm(yy,m-1);
  }
  return julianday;
}

int TimeStamp::get_dpm () const {
  return dpm(get_years(),get_months());
}

TimeStamp& TimeStamp::operator+=(const int seconds) {
  EKAT_REQUIRE_MSG(is_valid(), "Error! The time stamp contains uninitialized values.\n"
                                 "       To use this object, use operator= with a valid rhs first.\n");

  auto& sec  = m_time[2];
  auto& min  = m_time[1];
  auto& hour = m_time[0];
  auto& dd = m_date[2];
  auto& mm = m_date[1];
  auto& yy = m_date[0];

  EKAT_REQUIRE_MSG (seconds>=0, "Error! Cannot rewind time, sorry.\n");
  sec += seconds;

  // Carry over
  int carry;
  carry = sec / 60;
  if (carry==0) {
    return *this;
  }

  sec = sec % 60;
  min += carry;
  carry = min / 60;
  if (carry==0) {
    return *this;
  }
  min = min % 60;
  hour += carry;
  carry = hour / 24;

  if (carry==0) {
    return *this;
  }
  hour = hour % 24;
  dd += carry;

  while (dd>dpm(yy,mm)) {
    dd -= dpm(yy,mm);
    ++mm;
    
    if (mm>12) {
      ++yy;
      mm = 1;
    }
  }

  return *this;
}

bool operator== (const TimeStamp& ts1, const TimeStamp& ts2) {
  return ts1.get_date()==ts2.get_date() && ts1.get_time()==ts2.get_time();
}

bool operator< (const TimeStamp& ts1, const TimeStamp& ts2) {
  if (ts1.get_julian_day()<ts2.get_julian_day()) {
    return true;
  } else if (ts1.get_julian_day()==ts2.get_julian_day()) {
    return ts1.sec_of_day()<ts2.sec_of_day();
  }
  return false;
}

bool operator<= (const TimeStamp& ts1, const TimeStamp& ts2) {
  if (ts1.get_julian_day()<ts2.get_julian_day()) {
    return true;
  } else if (ts1.get_julian_day()==ts2.get_julian_day()) {
    return ts1.sec_of_day()<=ts2.sec_of_day();
  }
  return false;
}

TimeStamp operator+ (const TimeStamp& ts, const int dt) {
  TimeStamp sum = ts;
  sum += dt;
  return sum;
}

} // namespace util

} // namespace scream
