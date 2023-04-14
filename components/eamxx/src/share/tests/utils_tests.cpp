#include <catch2/catch.hpp>

#include "share/util/scream_array_utils.hpp"
#include "share/util/scream_universal_constants.hpp"
#include "share/util/scream_utils.hpp"
#include "share/util/scream_time_stamp.hpp"
#include "share/util/scream_setup_random_test.hpp"

TEST_CASE("contiguous_superset") {
  using namespace scream;

  std::string A = "A";
  std::string B = "B";
  std::string C = "C";
  std::string D = "D";
  std::string E = "E";
  std::string F = "F";
  std::string G = "G";

  using LOLS_type = std::list<std::list<std::string>>;

  // These three lists do not allow a superset from which they can all be
  // contiguously subviewed.
  LOLS_type lol1 = { {A,B}, {B,C}, {A,C} };
  REQUIRE(contiguous_superset(lol1).size()==0);

  // Input inner lists are not sorted
  REQUIRE_THROWS(contiguous_superset(LOLS_type{ {B,A} }));

  // The following should both allow the superset (A,B,C,D,E,F,G)
  // Note: lol3 is simply a shuffled version of lol2
  LOLS_type lol2 = { {A,B,C}, {B,C,D,E}, {C,D}, {C,D,E,F}, {D,E,F,G} };
  LOLS_type lol3 = { {D,E,F,G}, {C,D,E,F}, {A,B,C}, {C,D}, {B,C,D,E} };

  // Flipping a list is still a valid solution, so consider both tgt and its reverse.
  std::list<std::string> tgt = {A,B,C,D,E,F,G};
  std::list<std::string> tgt_rev = tgt;
  tgt_rev.reverse();

  auto superset2 = contiguous_superset(lol2);
  auto superset3 = contiguous_superset(lol3);
  REQUIRE ( (superset2==tgt || superset2==tgt_rev) );
  REQUIRE ( (superset3==tgt || superset3==tgt_rev) );
}

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

    REQUIRE (ts1.get_date_string()=="2021-10-12");
    REQUIRE (ts1.get_time_string()=="17:08:30");
    REQUIRE (ts1.to_string()=="2021-10-12-61710");
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

  SECTION ("leap_years") {
    // Check leap year correctness
    TS ts2({2000,2,28},{23,59,59});
    TS ts3({2012,2,28},{23,59,59});
    TS ts4({2100,2,28},{23,59,59});

    ts2 += 1;
    ts3 += 1;
    ts4 += 1;
#ifdef SCREAM_HAS_LEAP_YEAR
    REQUIRE (ts2.get_month()==2);
    REQUIRE (ts3.get_month()==2);
#else
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
}

TEST_CASE ("array_utils") {
  using namespace scream;

  auto engine = setup_random_test ();
  using IPDF = std::uniform_int_distribution<int>;
  IPDF pdf(1,10);

  auto total_size = [](const std::vector<int>& v) -> int {
    int s = 1;
    for (int i : v) {
      s *= i;
    }
    return s;
  };

  // Adds one to fastest striding, doing carrying (if possible) based on max dims d
  // Note: cannot use recursion with a pure lambda
  std::function<bool(int*,int,int*)> add_one = [&](int* v, int n, int* d) -> bool{
    // Increase fastest striding index
    ++v[n];

    // If we reached d[n], we need to carry
    if (v[n]>=d[n]) {
      if (n>0) {
        // Try to carry
        v[n] = 0;

        bool b = add_one(v,n-1,d);

        if (not b) {
          // There was no room to carry. Reset v[n]=d[n] and return false
          v[n] = d[n];
          return false;
        }
      } else {
        v[0] = d[0];
        return false;
      }
    }

    return true;
  };

  for (int rank : {1,2,3,4,5,6}) {
    std::vector<int> dims(rank);
    for (int d=0; d<rank; ++d) {
      dims[d] = pdf(engine);
    }

    std::vector<int> ind(rank,0);
    auto s = total_size(dims);
    for (int idx_1d=0; idx_1d<s; ++idx_1d) {
      auto idx_nd = unflatten_idx(dims,idx_1d);    

      std::cout << "idx1d: " << idx_1d << "\n";
      std::cout << "  indices:";
      for (auto i : ind) {
        std::cout << " " << i;
      }
      std::cout << "\n  unflatten:";
      for (auto i : idx_nd) {
        std::cout << " " << i;
      }
      std::cout << "\n";
      REQUIRE (idx_nd==ind);
      add_one(ind.data(),rank-1,dims.data());
    }
  }
}
