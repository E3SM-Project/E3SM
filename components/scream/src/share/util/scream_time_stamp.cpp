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

int days_in_month (const int yy, const int mm) {
  EKAT_REQUIRE_MSG (mm>=1 && mm<=12,
      "Error! Month out of bounds. Did you call `days_in_month` with yy and mm swapped?\n");
  constexpr int nonleap_days [12] = {31,28,31,30,31,30,31,31,30,31,30,31};
  constexpr int leap_days    [12] = {31,29,31,30,31,30,31,31,30,31,30,31};
  auto& arr = is_leap_year(yy) ? leap_days : nonleap_days;
  return arr[mm-1];
}

bool is_leap_year (const int yy) {
#ifdef SCREAM_HAS_LEAP_YEAR
  if (yy%4==0) {
    // Year is divisible by 4 (minimum requirement)
    if (yy%100 != 0) {
      // Not a centennial year => leap.
      return true;
    } else if ((yy/100)%4==0) {
      // Centennial year, AND first 2 digids divisible by 4 => leap
      return true;
    }
  }
#else
  (void)yy;
#endif
  // Either leap year not enabled, or not a leap year at all
  return false;
}

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
  EKAT_REQUIRE_MSG (dd>0   && dd<=days_in_month(yy,mm), "Error! Day out of bounds.\n");
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

bool TimeStamp::is_valid () const {
  return !(*this==TimeStamp());
}

int TimeStamp::seconds_from (const TimeStamp& ts) const {
  return *this-ts;
}

std::string TimeStamp::to_string () const {

  auto time = get_time_string ();
  time.erase(std::remove( time.begin(), time.end(), ':'), time.end());

  return get_date_string() + "." + time;
}

std::string TimeStamp::get_date_string () const {

  std::ostringstream ymd;

  ymd << std::setw(4) << std::setfill('0') << m_date[0] << "-";
  ymd << std::setw(2) << std::setfill('0') << m_date[1] << "-";
  ymd << std::setw(2) << std::setfill('0') << m_date[2];

  return ymd.str();
}

std::string TimeStamp::get_time_string () const {

  std::ostringstream hms;

  hms << std::setw(2) << std::setfill('0') << m_time[0] << ":";
  hms << std::setw(2) << std::setfill('0') << m_time[1] << ":";
  hms << std::setw(2) << std::setfill('0') << m_time[2];

  return hms.str();
}

double TimeStamp::frac_of_year_in_days () const {
  double doy = (m_date[2]-1) + sec_of_day() / 86400.0; // WARNING: avoid integer division
  for (int m=1; m<m_date[1]; ++m) {
    doy += days_in_month(m_date[0],m);
  }
  return doy;
}

TimeStamp& TimeStamp::operator+=(const int seconds) {
  EKAT_REQUIRE_MSG(is_valid(),
      "Error! The time stamp contains uninitialized values.\n"
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

  while (dd>days_in_month(yy,mm)) {
    dd -= days_in_month(yy,mm);
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
  if (ts1.get_year()<ts2.get_year()) {
    return true;
  } else if (ts1.get_year()==ts2.get_year()) {
    return ts1.frac_of_year_in_days()<ts2.frac_of_year_in_days();
  }
  return false;
}

bool operator<= (const TimeStamp& ts1, const TimeStamp& ts2) {
  if (ts1.get_year()<ts2.get_year()) {
    return true;
  } else if (ts1.get_year()==ts2.get_year()) {
    return ts1.frac_of_year_in_days()<=ts2.frac_of_year_in_days();
  }
  return false;
}

TimeStamp operator+ (const TimeStamp& ts, const int dt) {
  TimeStamp sum = ts;
  sum += dt;
  return sum;
}

long long operator- (const TimeStamp& ts1, const TimeStamp& ts2) {
  if (ts1<ts2) {
    return -(ts2-ts1);
  }

  long long diff = 0;
  const auto y1 = ts1.get_year();
  const auto y2 = ts2.get_year();
  const auto m1 = ts1.get_month();
  const auto m2 = ts2.get_month();
  const auto d1 = ts1.get_day();
  const auto d2 = ts2.get_day();

  // Strategy:
  //  - if y1>y2: add fraction of y1, process whole years in (y1,y2), then add fraction of y2.
  //  - if m1>m2: add fraction of m1, process whole months in (m1,m2), then add fraction of m2.
  //  - if d1>m2: add fraction of d1, process whole days in (d1,d2), then add fraction of d2.
  constexpr auto spd = constants::seconds_per_day;
  if (y1>y2) {
    // Rest of the day, month, and year in ts2's current year
    diff += spd - ts2.sec_of_day();
    diff += spd*(days_in_month(y2,m2)-d2);
    for (int m=m2+1; m<=12; ++m) {
      diff += spd*days_in_month(y2,m);
    }

    // Whole years in (y2,y1)
    for (int y=y2+1; y<y1; ++y) {
      diff += spd*(is_leap_year(y) ? 366 : 365);
    }

    // Chunk of current year in ts1
    for (int m=1; m<m1; ++m) {
      diff += spd*days_in_month(y1,m);
    }
    diff += spd*(d1-1);
    diff += ts1.sec_of_day();
  } else if (m1>m2) {
    // Rest of the day and month in ts2's current month
    diff += spd - ts2.sec_of_day();
    diff += spd*(days_in_month(y2,m2)-d2);

    // Whole months in (m2,m1)
    for (int m=m2+1; m<m1; ++m) {
      diff += spd*days_in_month(y1,m);
    }

    // Chunk of current month in ts1
    diff += spd*(d1-1);
    diff += ts1.sec_of_day();
  } else if (d1>d2) {
    // Rest of the day in ts2's current day
    diff += spd-ts2.sec_of_day();

    // Whole days in (d2,d1)
    diff += (d1-d2-1)*spd;

    // Chunk of day in current day in ts1
    diff += ts1.sec_of_day();
  } else {
    diff += ts1.sec_of_day() - ts2.sec_of_day();
  }

  return diff;
}

} // namespace util

} // namespace scream
