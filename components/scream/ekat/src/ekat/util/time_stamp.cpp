#include "ekat/util/time_stamp.hpp"

#include "ekat/scream_assert.hpp"
#include "ekat/scream_universal_constants.hpp"

#include <numeric>
#include <cmath>
#include <vector>
#include <algorithm>

namespace scream {
namespace util {

namespace {

// Utility functions
bool is_leap (const int yy) {
#ifdef SCREAM_HAS_LEAP_YEAR
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

int dpy (const int yy) {
  int n = constants::days_per_nonleap_year;
  if (is_leap(yy)) {
    ++n;
  }
  return n;
}

std::pair<int,int> get_month_and_day (const int year, const int day) {
  std::vector<int> months = {31, (is_leap(year) ? 29 : 28), 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  std::vector<int> offsets(1,0);
  std::partial_sum (months.begin(), months.end(), std::back_inserter(offsets));

  auto it = std::upper_bound(offsets.begin(),offsets.end(), day);
  scream_require_msg (it!=offsets.end(), "Error! Something went wrong while retrieving month from day '" + std::to_string(day) + "'.\n");

  int month = std::distance(offsets.begin(),it) - 1;
  return std::make_pair(month,day-offsets[month]);
}

} // anonymous namespace

TimeStamp::TimeStamp()
 : m_yy (std::numeric_limits<int>::lowest())
 , m_dd (std::numeric_limits<int>::lowest())
 , m_ss (std::numeric_limits<Real>::lowest())
{
  // Nothing to do here
}

TimeStamp::TimeStamp(const int yy, const int dd, const Real ss)
 : m_yy(yy)
 , m_dd(dd)
 , m_ss(ss)
{
  // Check the days and seconds numbers are non-negative.
  scream_require_msg (dd>=0, "Error! Days are negative.\n");
  scream_require_msg (ss>=0, "Error! Seconds are negative.\n");

  // Adjust if input numbers are too large
  int carry;
  carry = static_cast<int>(std::floor(m_ss / constants::seconds_per_day));
  m_ss  = std::fmod(m_ss, constants::seconds_per_day);
  m_dd += carry;

  while (m_dd>dpy(m_yy)) {
    m_dd -= dpy(m_yy);
    ++m_yy;    
  }
}

std::string TimeStamp::to_string () const {
  auto md = get_month_and_day (m_yy,m_dd);
  const int ss = static_cast<int>(m_ss);
  const int h =  ss / 3600;
  const int m = (ss % 3600) / 60;
  const int s = (ss % 3600) % 60;
  const std::string zero = "00";
  // For h:m:s, check if 0, and if so, use "00" rather than to_string, which returns "0"
  return std::to_string(md.first+1) + "-" + std::to_string(md.second+1) + "-" + std::to_string(m_yy) + " " + 
         (h==0 ? zero : std::to_string(h)) + ":" + (m==0 ? zero : std::to_string(m)) + ":" + (s==0 ? zero : std::to_string(s));
}

TimeStamp& TimeStamp::operator+=(const Real seconds) {
  return (*this += TimeStamp(0,0,seconds));
}

TimeStamp& TimeStamp::operator+=(const TimeStamp& dt) {
  scream_require_msg(is_valid(), "Error! The time stamp contains uninitialized values.\n"
                                 "       To use this object, use operator= with a valid rhs first.\n");
  scream_require_msg(dt.is_valid(), "Error! The input time step contains uninitialized values.\n");
  scream_require_msg(TimeStamp(0,0,0)<dt,
                     "Error! Time can only shift forward during a simulation.\n");

  const int ndays = dpy(m_yy);

  int carry;
  m_ss += dt.get_seconds();
  carry = static_cast<int>(std::floor(m_ss / constants::seconds_per_day));
  m_ss  = std::fmod(m_ss, constants::seconds_per_day);

  m_dd += dt.get_days() + carry;
  carry = m_dd / ndays;
  m_dd  = m_dd % ndays;
 
  m_yy += dt.get_years() + carry;

  return *this;
}

bool operator== (const TimeStamp& ts1, const TimeStamp& ts2) {
  return ts1.get_seconds()==ts2.get_seconds() &&
         ts1.get_days()==ts2.get_days() &&
         ts1.get_years()==ts2.get_years();
}

bool operator< (const TimeStamp& ts1, const TimeStamp& ts2) {
  if (ts1.get_years()>ts2.get_years()) {
    return false;
  } else if (ts1.get_years()==ts2.get_years()) {
    if (ts1.get_days()>ts2.get_days()) {
      return false;
    } else if (ts1.get_days()==ts2.get_days()) {
      return ts1.get_seconds()<ts2.get_seconds();
    }
  }
  return true;
}

bool operator<= (const TimeStamp& ts1, const TimeStamp& ts2) {
  if (ts1.get_years()>ts2.get_years()) {
    return false;
  } else if (ts1.get_years()==ts2.get_years()) {
    if (ts1.get_days()>ts2.get_days()) {
      return false;
    } else if (ts1.get_days()==ts2.get_days()) {
      return ts1.get_seconds()<=ts2.get_seconds();
    }
  }
  return true;
}

TimeStamp operator- (const TimeStamp& ts1, const TimeStamp& ts2) {
  scream_require_msg (ts2<ts1, "Error! We only allow to compute positive time stamp differences");

  Real ss1 = ts1.get_seconds();
  int dd1 = ts1.get_days();
  int yy1 = ts1.get_years();


  if (ss1 < ts2.get_seconds()) {
    ss1 += constants::seconds_per_day;
    --dd1;
  }
  if (dd1 < ts2.get_days()) {
    dd1 += dpy(yy1);
    --yy1;
  }

  return TimeStamp (yy1-ts2.get_years(),dd1-ts2.get_days(),ss1-ts2.get_seconds());
}

TimeStamp operator+ (const TimeStamp& ts, const TimeStamp& dt) {
  TimeStamp sum = ts;
  sum += dt;
  return sum;
}

TimeStamp operator+ (const TimeStamp& ts, const Real dt) {
  TimeStamp sum = ts;
  sum += dt;
  return sum;
}

} // namespace util

} // namespace scream

