#include <catch2/catch.hpp>

#include "share/util/eamxx_universal_constants.hpp"
#include "share/util/eamxx_time_stamp.hpp"
#include "share/core/eamxx_config.hpp"

TEST_CASE ("time_stamp") {
  using namespace scream;
  using TS = util::TimeStamp;

  constexpr auto spd = constants::seconds_per_day;

  TS ts1 (2021,10,12,17,8,30);

  SECTION ("ctor_check") {
    REQUIRE (ts1.get_year()==2021);
    REQUIRE (ts1.get_month()==10);
    REQUIRE (ts1.get_day()==12);
    REQUIRE (ts1.get_hours()==17);
    REQUIRE (ts1.get_minutes()==8);
    REQUIRE (ts1.get_seconds()==30);
  }

  SECTION ("getters_checks") {
    // Julian day = frac_of_year_in_days.fraction_of_day, with frac_of_year_in_days=0 at Jan 1st.
    REQUIRE (ts1.frac_of_year_in_days()==(284 + (17*3600+8*60+30)/86400.0));
    REQUIRE (ts1.get_num_steps()==0);
  }

  SECTION ("formatting") {
    REQUIRE (ts1.get_date_string()=="2021-10-12");
    REQUIRE (ts1.get_time_string()=="17:08:30");
    REQUIRE (ts1.to_string()=="2021-10-12-61710");
    REQUIRE (util::str_to_time_stamp("2021-10-12-61710")==ts1);
  }

  SECTION ("comparisons") {
    REQUIRE (ts1==ts1);

    // Comparisons
    REQUIRE ( TS({2021,12,31},{23,59,59}) < TS({2022,1,1},{0,0,0}));
    REQUIRE ( TS({2022,1,1},{0,0,0}) <= TS({2022,1,1},{0,0,0}));
    REQUIRE ( (TS({2021,12,31},{23,59,59})+1) == TS({2022,1,1},{0,0,0}));
  }

  SECTION ("updates") {
    // Cannot rewind time
    REQUIRE_THROWS (ts1+=-10);

    auto ts2 = ts1 + 1;

    REQUIRE (ts1<ts2);
    REQUIRE (ts2<=ts2);
    REQUIRE ( (ts2-1)==ts1 );

    // Update: check carries
    REQUIRE (ts2.get_seconds()==(ts1.get_seconds()+1));
    REQUIRE (ts2.get_minutes()==ts1.get_minutes());
    REQUIRE (ts2.get_hours()==ts1.get_hours());
    REQUIRE (ts2.get_day()==ts1.get_day());
    REQUIRE (ts2.get_month()==ts1.get_month());
    REQUIRE (ts2.get_year()==ts1.get_year());

    ts2 += 60;
    REQUIRE (ts2.get_seconds()==(ts1.get_seconds()+1));
    REQUIRE (ts2.get_minutes()==(ts1.get_minutes()+1));
    REQUIRE (ts2.get_hours()==ts1.get_hours());
    REQUIRE (ts2.get_day()==ts1.get_day());
    REQUIRE (ts2.get_month()==ts1.get_month());
    REQUIRE (ts2.get_year()==ts1.get_year());

    ts2 += 3600;
    REQUIRE (ts2.get_seconds()==(ts1.get_seconds()+1));
    REQUIRE (ts2.get_minutes()==(ts1.get_minutes()+1));
    REQUIRE (ts2.get_hours()==ts1.get_hours()+1);
    REQUIRE (ts2.get_day()==ts1.get_day());
    REQUIRE (ts2.get_month()==ts1.get_month());
    REQUIRE (ts2.get_year()==ts1.get_year());

    ts2 += spd;
    REQUIRE (ts2.get_seconds()==(ts1.get_seconds()+1));
    REQUIRE (ts2.get_minutes()==(ts1.get_minutes()+1));
    REQUIRE (ts2.get_hours()==(ts1.get_hours()+1));
    REQUIRE (ts2.get_day()==(ts1.get_day()+1));
    REQUIRE (ts2.get_month()==ts1.get_month());
    REQUIRE (ts2.get_year()==ts1.get_year());

    ts2 += spd*20;
    REQUIRE (ts2.get_seconds()==(ts1.get_seconds()+1));
    REQUIRE (ts2.get_minutes()==(ts1.get_minutes()+1));
    REQUIRE (ts2.get_hours()==(ts1.get_hours()+1));
    REQUIRE (ts2.get_day()==(ts1.get_day()+1+20-31)); // Add 20 days, subtract Oct 31 days (carry)
    REQUIRE (ts2.get_month()==(ts1.get_month()+1));
    REQUIRE (ts2.get_year()==ts1.get_year());

    ts2 += spd*365;
    REQUIRE (ts2.get_seconds()==ts1.get_seconds()+1);
    REQUIRE (ts2.get_minutes()==(ts1.get_minutes()+1));
    REQUIRE (ts2.get_hours()==(ts1.get_hours()+1));
    REQUIRE (ts2.get_day()==(ts1.get_day()+1+20-31)); // Add 20 days, subtract Oct 31 days (carry)
    REQUIRE (ts2.get_month()==(ts1.get_month()+1));
    REQUIRE (ts2.get_year()==(ts1.get_year()+1));

    REQUIRE (ts2.get_num_steps()==6);
  }

  SECTION ("fractional_update") {
    // Check update with fractional seconds
    REQUIRE ((ts1+0.999)==ts1);
    REQUIRE ((ts1+0.9999)!=ts1); // When seconds frac is <0.001 or >0.999 we round
  }

  SECTION ("leap_years") {
    // Check leap year correctness
    TS ts2({2000,2,28},{23,59,59});
    TS ts3({2012,2,28},{23,59,59});
    TS ts4({2100,2,28},{23,59,59});

    ts2 += 1;
    ts3 += 1;
    ts4 += 1;
#ifdef SCREAM_HAS_LEAP_YEAR
    REQUIRE (use_leap_year());
    REQUIRE (ts2.get_month()==2);
    REQUIRE (ts3.get_month()==2);
#else
    REQUIRE (not use_leap_year());
    REQUIRE (ts2.get_month()==3);
    REQUIRE (ts3.get_month()==3);
#endif
    // Centennial years with first 2 digits not divisible by 4 are not leap
    REQUIRE (ts4.get_month()==3);
  }

  SECTION ("difference") {
    // Difference
    auto ts2 = ts1 + 3600;
    REQUIRE ( (ts2-ts1)==3600 );
    auto ts3 = ts1 + spd;
    REQUIRE ( (ts3-ts1)==spd );
    auto ts4 = ts1 + spd*10;
    REQUIRE ( (ts4-ts1)==spd*10 );
    auto ts5 = ts1 + spd*100;
    REQUIRE ( (ts5-ts1)==spd*100 );
    auto ts6 = ts1 + spd*1000;
    REQUIRE ( (ts6-ts1)==spd*1000 );
  }
  SECTION ("large_updates") {
    TS base ({1850,1,1},{0,0,0},1);
    for (int i=1; i<=500; ++i) {
      TS curr ({1850+i,1,1},{0,0,0},0);
      auto diff = curr.seconds_from(base);
      TS time = base;
      time += diff;
      REQUIRE (time==curr);
    }
  }
}
