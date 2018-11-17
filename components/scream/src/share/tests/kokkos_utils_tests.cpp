#include "catch2/catch.hpp"

#include "share/util/scream_kokkos_utils.hpp"

namespace {

TEST_CASE("data_type", "") {
  using namespace scream::util;

  // Check meta-util to get underlying scalar type of a raw MD array
  REQUIRE(std::is_same<ValueType<double**&>::type,double>::value);
  REQUIRE(std::is_same<ValueType<double*[3]>::type,double>::value);
  REQUIRE(std::is_same<ValueType<double[2][3]>::type,double>::value);

  // Check meta-util to get rank and dynamic rank of a raw MD array
  REQUIRE(GetRanks<double[2][3]>::rank==2);
  REQUIRE(GetRanks<double[2][3]>::rank_dynamic==0);
  REQUIRE(GetRanks<double*[2][3]>::rank==3);
  REQUIRE(GetRanks<double*[2][3]>::rank_dynamic==1);
  REQUIRE(GetRanks<double**[2][3]>::rank==4);
  REQUIRE(GetRanks<double**[2][3]>::rank_dynamic==2);

  // Check meta-util that allows to reshape a view
  Kokkos::View<double*> v1d("",100);
  auto v2d = reshape<double[2][50]>(v1d);
  REQUIRE(v2d.size()==100);

  auto v3d = reshape<double*[5][5]>(v2d,4);
  REQUIRE (v3d.size()==100);
}

} // anonymous namespace
