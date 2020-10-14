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

int dpy (const int yy) {
  int n = constants::days_per_nonleap_year;
  if (is_leap(yy)) {
    ++n;
  }
  return n;
}

int dpm (const int yy, const int mm) {
  auto& arr = is_leap(yy) ? leap_days : nonleap_days;
  return arr[mm];
}

} // anonymous namespace

TimeStamp::TimeStamp()
 : m_yy (std::numeric_limits<int>::lowest())
 , m_mm (std::numeric_limits<int>::lowest())
 , m_dd (std::numeric_limits<int>::lowest())
 , m_ss (std::numeric_limits<double>::lowest())
{
  // Nothing to do here
}

TimeStamp::TimeStamp(const int yy,
                     const int mm,
                     const int dd,
                     const double ss)
 : m_yy(yy)
 , m_mm(mm)
 , m_dd(dd)
 , m_ss(ss)
{
  // Check the days and seconds numbers are non-negative.
  EKAT_REQUIRE_MSG (mm>=0, "Error! Month is negative.\n");
  EKAT_REQUIRE_MSG (mm<=11, "Error! Month is too large.\n");
  EKAT_REQUIRE_MSG (dd>=0, "Error! Day is negative.\n");
  EKAT_REQUIRE_MSG (dd<(is_leap(m_yy) ? leap_days[mm] : nonleap_days[mm]), "Error! Day is too large.\n");
  EKAT_REQUIRE_MSG (ss>=0, "Error! Seconds are negative.\n");
  EKAT_REQUIRE_MSG (ss<constants::seconds_per_day, "Error! Seconds are too large.\n");

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
  const int ss = static_cast<int>(m_ss);
  const int h =  ss / 3600;
  const int m = (ss % 3600) / 60;
  const int s = (ss % 3600) % 60;

  std::ostringstream ymdhms;
  ymdhms << std::setw(4) << std::setfill('0') << m_yy << "-";
  ymdhms << std::setw(2) << std::setfill('0') << m_mm+1 << "-";
  ymdhms << std::setw(2) << std::setfill('0') << m_dd+1 << " ";
  ymdhms << std::setw(2) << std::setfill('0') << h << ":";
  ymdhms << std::setw(2) << std::setfill('0') << m << ":";
  ymdhms << std::setw(2) << std::setfill('0') << s;
  return ymdhms.str();

}

TimeStamp& TimeStamp::operator+=(const double seconds) {
  EKAT_REQUIRE_MSG(is_valid(), "Error! The time stamp contains uninitialized values.\n"
                                 "       To use this object, use operator= with a valid rhs first.\n");

  int carry;
  m_ss += seconds;
  carry = static_cast<int>(std::floor(m_ss / constants::seconds_per_day));
  m_ss  = std::fmod(m_ss, constants::seconds_per_day);

  if (carry>0) {

    m_dd += carry;
    while (m_dd>dpm(m_yy,m_mm)) {
      m_dd -= dpm(m_yy,m_mm);
      ++m_mm;
      m_yy += m_mm / 12;
      m_mm  = m_mm % 12;
    }
  }

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

double operator- (const TimeStamp& ts1, const TimeStamp& ts2) {
  if (ts1<ts1) {
    return -(ts2-ts1);
  }

  constexpr int spd = constants::seconds_per_day;
  double dt = 0;

  const int y1 = ts1.get_years();
  const int y2 = ts2.get_years();
  const int m1 = ts1.get_months();
  const int m2 = ts2.get_months();
  const int d1 = ts1.get_days();
  const int d2 = ts2.get_days();
  const int s1 = ts1.get_seconds();
  const int s2 = ts2.get_seconds();

  if (y1!=y2) {
    // 1. Process from ts2 till end of year

    // Rest of day
    dt += spd-s2;

    // Rest of month
    for (int d=d2+1; d<dpm(y2,m2); ++d) {
      dt += spd;
    }

    // Rest of year
    for (int m=m2+1; m<12; ++m) {
      dt += spd*dpm(y2,m);
    }

    // 2. Process whole years till y1-1 (if any)
    for (int y=y2+1; y<y1; ++y) {
      dt += spd*dpy(y);
    }

    // 3. Process chunk of last year in ts1

    // Previous months
    for (int m=0; m<m1; ++m) {
      dt += spd*dpm(y1,m);
    }

    // Previous days
    for (int d=0; d<d1; ++d) {
      dt += spd;
    }

    // Previous seconds
    dt += s1;
  } else if (m1!=m2) {
    // 1. Process from ts2 till end of the month

    // Rest of day
    dt += spd-s2;

    // Rest of month
    for (int d=d2+1; d<dpm(y2,m2); ++d) {
      dt += spd;
    }

    // 2. Process whole months till m1-1 (if any)
    for (int m=m2+1; m<m1; ++m) {
      dt += spd*dpm(y1,m);
    }

    // 3. Process chunk of last month in ts1

    // Previous days
    for (int d=0; d<d1; ++d) {
      dt += spd;
    }

    // Previous seconds
    dt += s1;
  } else if (d1!=d2) {
    // 1. Process from ts2 till end of the day
    dt += spd-s2;

    // 2. Process whole days till d1-1 (if any)
    for (int d=d2+1; d<d1; ++d) {
      dt += spd;
    }

    // 3. Process chunk of last day in ts1
    dt += s1;
  } else {
    dt = s1 - s2;
  }
    
  return dt;
}

TimeStamp operator+ (const TimeStamp& ts, const double dt) {
  TimeStamp sum = ts;
  sum += dt;
  return sum;
}

} // namespace util

} // namespace scream
