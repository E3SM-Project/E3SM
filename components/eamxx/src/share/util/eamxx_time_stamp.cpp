#include "share/util/eamxx_time_stamp.hpp"

#include "share/util/eamxx_universal_constants.hpp"
#include "share/eamxx_config.hpp"
#include "ekat/ekat_assert.hpp"

#include <limits>
#include <numeric>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>

namespace scream {
namespace util {

bool is_leap_year (const int yy) {
  if (use_leap_year()) {
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
  }
  return false;
}

int days_in_month (const int yy, const int mm) {
  EKAT_REQUIRE_MSG (mm>=1 && mm<=12,
      "Error! Month out of bounds. Did you call `days_in_month` with yy and mm swapped?\n");
  constexpr int nonleap_days [12] = {31,28,31,30,31,30,31,31,30,31,30,31};
  constexpr int leap_days    [12] = {31,29,31,30,31,30,31,31,30,31,30,31};
  auto& arr = is_leap_year(yy) ? leap_days : nonleap_days;
  return arr[mm-1];
}

TimeStamp::TimeStamp()
 : m_date (3,std::numeric_limits<int>::lowest())
 , m_time (3,std::numeric_limits<int>::lowest())
{
  // Nothing to do here
}

TimeStamp::TimeStamp(const std::vector<int>& date,
                     const std::vector<int>& time,
                     const int num_steps)
{
  EKAT_REQUIRE_MSG (date.size()==3, "Error! Date should consist of three ints: [year, month, day].\n");
  EKAT_REQUIRE_MSG (time.size()==3, "Error! Time of day should consist of three ints: [hour, min, sec].\n");
  EKAT_REQUIRE_MSG (num_steps>=0,   "Error! Number of steps should be a non-negative number.\n");

  const auto yy   = date[0];
  const auto mm   = date[1];
  const auto dd   = date[2];
  const auto hour = time[0];
  const auto min  = time[1];
  const auto sec  = time[2];

  // Check the days and seconds numbers are non-negative.
  EKAT_REQUIRE_MSG (yy>=0   && yy<=9999, "Error! Year out of bounds.\n");
  EKAT_REQUIRE_MSG (mm>0    && mm<=12, "Error! Month out of bounds.\n");
  EKAT_REQUIRE_MSG (dd>0    && dd<=days_in_month(yy,mm), "Error! Day out of bounds.\n");
  EKAT_REQUIRE_MSG (sec>=0  && sec<60, "Error! Seconds out of bounds.\n");
  EKAT_REQUIRE_MSG (min>=0  && min<60, "Error! Minutes out of bounds.\n");
  EKAT_REQUIRE_MSG (hour>=0 && hour<24, "Error! Hours out of bounds.\n");

  // All good, store
  m_date = date;
  m_time = time;
  m_num_steps = num_steps;
}

TimeStamp::TimeStamp(const int yy, const int mm, const int dd,
                     const int h, const int min, const int sec,
                     const int num_steps)
 : TimeStamp({yy,mm,dd},{h,min,sec},num_steps)
{
  // Nothing to do here
}

bool TimeStamp::is_valid () const {
  return !(*this==TimeStamp());
}

std::int64_t TimeStamp::seconds_from (const TimeStamp& ts) const {
  return *this-ts;
}

double TimeStamp::days_from (const TimeStamp& ts) const {
  return seconds_from(ts)/86400.0;
}

std::string TimeStamp::to_string () const {

  std::ostringstream tod;
  tod << std::setw(5) << std::setfill('0') << sec_of_day();

  return get_date_string() + "-" + tod.str();
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

int TimeStamp::days_in_curr_month () const
{
  return days_in_month(m_date[0],m_date[1]);
}

int TimeStamp::days_in_curr_year () const
{
  return is_leap_year(m_date[0]) ? 366 : 365;
}

TimeStamp TimeStamp::curr_month_beg () const
{
  auto date = m_date;
  date[2] = 1;
  return TimeStamp (date,{0,0,0});
}

TimeStamp& TimeStamp::operator+=(const double seconds) {
  // Sanity checks
  // Note: (x-int(x)) only works for x small enough that can be stored in an int,
  //       but that should be the case here, for use cases in EAMxx.
  EKAT_REQUIRE_MSG (seconds>=0, "Error! Cannot rewind time.\n");
  EKAT_REQUIRE_MSG ((seconds-round(seconds))<std::numeric_limits<double>::epsilon()*10,
      "Error! Cannot update TimeStamp with non-integral number of seconds " << seconds << "\n");

  EKAT_REQUIRE_MSG(is_valid(),
      "Error! The time stamp contains uninitialized values.\n"
      "       To use this object, use operator= with a valid rhs first.\n");

  auto& sec  = m_time[2];
  auto& min  = m_time[1];
  auto& hour = m_time[0];
  auto& dd = m_date[2];
  auto& mm = m_date[1];
  auto& yy = m_date[0];

  ++m_num_steps;

  constexpr auto spd = constants::seconds_per_day;
  auto add_days = std::floor(seconds / spd);
  auto add_secs = seconds - add_days*spd;

  sec += add_secs;
  dd  += add_days;

  // Carry over
  int carry;
  carry = sec / 60;
  if (carry!=0) {
    sec = sec % 60;
    min += carry;
    carry = min / 60;
    if (carry!=0) {
      min = min % 60;
      hour += carry;
      carry = hour / 24;

      if (carry!=0) {
        hour = hour % 24;
        dd += carry;
      }
    }
  }

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

TimeStamp TimeStamp::clone (const int num_steps) {
  return TimeStamp (get_date(),get_time(),num_steps>=0 ? num_steps : m_num_steps);
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

std::int64_t operator- (const TimeStamp& ts1, const TimeStamp& ts2) {
  if (ts1<ts2) {
    return -(ts2-ts1);
  }

  std::int64_t diff = 0;
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

TimeStamp operator- (const TimeStamp& ts, const int dt) {
  if (dt<0) {
    return ts + (-dt);
  }

  const int ndays = dt / constants::seconds_per_day;
  const int nsecs = dt % constants::seconds_per_day;

  auto rewind_one_day = [] (std::vector<int>& date) {
    // Rewind a day
    --date[2];
    if (date[2]==0) {
      // Whoops, rewind a month
      --date[1];
      if (date[1]==0) {
        // Whoops, rewind a year
        --date[0];
        date[1]=12;
      }
      date[2] = days_in_month(date[0],date[1]);
    }
  };

  auto date = ts.get_date();

  for (int i=0; i<ndays; ++i) {
    rewind_one_day (date);
  }

  int frac_of_day;
  if (nsecs>ts.sec_of_day()) {
    rewind_one_day (date);
    frac_of_day = constants::seconds_per_day - (nsecs - ts.sec_of_day());
  } else {
    frac_of_day = ts.sec_of_day()-nsecs;
  }

  TimeStamp t (date,{0,0,0});
  t += frac_of_day;
  return t;
}

TimeStamp str_to_time_stamp (const std::string& s)
{
  auto is_int = [](const std::string& s) -> bool {
    for (const auto& c : s) {
      if (not std::isdigit(c)) {
        return false;
      }
    }
    return true;
  };
  // A time stamp is a string of the form "YYYY-MM-DD-XXXXX"
  // with XXXXX being the time of day in seconds
  if (s.size()!=16 || s[4]!='-' || s[7]!='-' || s[10]!='-') {
    return util::TimeStamp();
  }
  // Check subtokens are valid ints
  auto YY  = s.substr(0,4);
  auto MM  = s.substr(5,2);
  auto DD  = s.substr(8,2);
  auto tod = s.substr(11,6);
  if (not is_int(YY) || not is_int(MM) || not is_int(DD) || not is_int(tod)) {
    return util::TimeStamp();
  }
  std::vector<int> date{std::stoi(YY),std::stoi(MM),std::stoi(DD)};
  auto sec_of_day = std::stoi(tod);
  auto hh = (sec_of_day / 60) / 60;
  auto mm = (sec_of_day / 60) % 60;
  auto ss =  sec_of_day % 60;
  std::vector<int> time{hh,mm,ss};

  try {
    return util::TimeStamp(date,time);
  } catch (...) {
    return util::TimeStamp();
  }
};

} // namespace util

} // namespace scream
