#include "share/core/eamxx_config.hpp"
#include "share/core/eamxx_session.hpp"
#include "share/core/eamxx_types.hpp"

#include "eamxx_version.h"

#include <ekat_kokkos_session.hpp>
#include <ekat_arch.hpp>
#include <ekat_fpe.hpp>
#include <ekat_assert.hpp>
#include <algorithm>
#include <cctype>
#include <cmath>

namespace scream {

std::string eamxx_version () {
  return EAMXX_VERSION;
}

std::string eamxx_git_version () {
  return EAMXX_GIT_VERSION;
}

std::string eamxx_config_string() {
  std::string config = "\n-------- EKAT CONFIGS --------\n\n";
  config += "Active AVX settings: " + ekat::active_avx_string () + "\n";
  config += "Compiler Id: " + ekat::compiler_id_string () + "\n";
  config += ekat::kokkos_config_string();
  config += "\n-------- EAMXX CONFIGS --------\n\n";
  config += " sizeof(Real) = " + std::to_string(sizeof(Real)) + "\n";
  config += " default pack size = " + std::to_string(SCREAM_PACK_SIZE) + "\n";
  config += " default FPE mask: " +
      ( get_default_fpes() == 0 ? "0 (NONE) \n" :
        std::to_string(get_default_fpes()) + " (FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW) \n");
  config += "-------------------------------\n";

  return config;
}

bool& use_leap_year_impl () {
#ifdef SCREAM_HAS_LEAP_YEAR
  static bool use_leap = true;
#else
  static bool use_leap = false;
#endif
  return use_leap;
}

bool use_leap_year () {
  return use_leap_year_impl ();
}

void set_use_leap_year (const bool use_leap) {
  use_leap_year_impl () = use_leap;
}

namespace {
bool noleap_input_enabled = false;
bool gregorian_input_enabled = false;
input_calendar_seconds calendar_seconds = nullptr;
input_calendar_resolver calendar_resolver = nullptr;
}

bool use_noleap_input () {
  return noleap_input_enabled;
}

void set_use_noleap_input (const bool enabled) {
  EKAT_REQUIRE_MSG(not enabled or not gregorian_input_enabled,
      "INPUT_CALENDAR: input calendar policies are mutually exclusive");
  EKAT_REQUIRE_MSG(not enabled or use_leap_year(),
      "NOLEAP_INPUT: fudge_noleap_for_gregorian requires the Gregorian model calendar.");
  noleap_input_enabled = enabled;
}

bool use_gregorian_input () {
  return gregorian_input_enabled;
}

bool use_input_calendar_policy () {
  return noleap_input_enabled or gregorian_input_enabled;
}

void set_use_gregorian_input (const bool enabled) {
  EKAT_REQUIRE_MSG(not enabled or not noleap_input_enabled,
      "INPUT_CALENDAR: input calendar policies are mutually exclusive");
  EKAT_REQUIRE_MSG(not enabled or not use_leap_year(),
      "GREGORIAN_INPUT: fudge_gregorian_to_noleap requires the no-leap model calendar.");
  gregorian_input_enabled = enabled;
}

void set_input_calendar_seconds(input_calendar_seconds callback) {
  calendar_seconds = callback;
}

double gregorian_input_seconds(const int date, const double seconds) {
  if (calendar_seconds) return calendar_seconds(date,seconds);
  EKAT_REQUIRE_MSG(is_scream_standalone(),"GREGORIAN_INPUT: coupled input-time callback is missing.");
  const int year = date/10000;
  const int month = (date/100)%100;
  const int day = date%100;
  const bool leap = year%4==0 and (year%100!=0 or year%400==0);
  const int lengths[] = {31,leap ? 29 : 28,31,30,31,30,31,31,30,31,30,31};
  EKAT_REQUIRE_MSG(date>=0 and month>=1 and month<=12 and day>=1 and
      day<=lengths[month-1] and std::isfinite(seconds) and seconds>=0 and seconds<86400,
      "GREGORIAN_INPUT: invalid source date or time of day.");
  long long days = 365LL*year+(year+3)/4-(year+99)/100+(year+399)/400+day-1;
  for (int previous=1; previous<month; ++previous) days += lengths[previous-1];
  return (days-730485LL)*86400.0+seconds;
}

void set_input_calendar_resolver(input_calendar_resolver resolver) {
  calendar_resolver = resolver;
}

void check_input_calendar_range(const std::string& flavor, int epoch, double seconds,
                               double offset_seconds, const std::string& context) {
  if (flavor!="STANDARD") return;
  const double cutoff = gregorian_input_seconds(15821015,0);
  const double source = gregorian_input_seconds(epoch,seconds)+offset_seconds;
  EKAT_REQUIRE_MSG(epoch>=15821015 and std::isfinite(source) and source>=cutoff,
      "INPUT_CALENDAR: unsupported historical calendar=" + flavor + "; epoch=" + std::to_string(epoch) +
      "; offset seconds=" + std::to_string(offset_seconds) + "; require epoch and range >= 1582-10-15; " + context);
}

std::string resolve_input_calendar(const std::string& filename, const std::string& declared,
                                  bool metadata_present, std::string& provenance, std::string& flavor) {
  if (calendar_resolver) {
    char calendar[64] = {};
    char origin[64] = {};
    char source_flavor[64] = {};
    char diagnostic[1024] = {};
    const int status = calendar_resolver(filename.data(),filename.size(),declared.data(),declared.size(),
        metadata_present ? 1 : 0,calendar,origin,source_flavor,diagnostic);
    EKAT_REQUIRE_MSG(status==0,std::string(diagnostic) + "; file=" + filename);
    provenance = origin;
    flavor = source_flavor;
    return calendar;
  }
  EKAT_REQUIRE_MSG(is_scream_standalone(),"NOLEAP_INPUT: coupled calendar resolver was not registered.");
  EKAT_REQUIRE_MSG(metadata_present,"NOLEAP_INPUT: explicit time calendar metadata is required; file=" + filename);
  auto token = declared;
  std::transform(token.begin(),token.end(),token.begin(),[](unsigned char character) {
    if (character==0 or character=='\t' or character=='\n' or character=='\r') return ' ';
    return static_cast<char>(std::tolower(character));
  });
  const auto first = token.find_first_not_of(' ');
  EKAT_REQUIRE_MSG(first!=std::string::npos,"NOLEAP_INPUT: empty calendar metadata; file=" + filename);
  token = token.substr(first,token.find_last_not_of(' ')-first+1);
  provenance = "metadata";
  flavor = "NO_LEAP";
  if (token=="noleap" or token=="no_leap" or token=="365_day" or token=="365day") return "NO_LEAP";
  flavor = token=="proleptic_gregorian" ? "PROLEPTIC_GREGORIAN" : "STANDARD";
  if (token=="gregorian" or token=="standard" or token=="proleptic_gregorian") return "GREGORIAN";
  EKAT_ERROR_MSG("NOLEAP_INPUT: unsupported calendar=" + token + " file=" + filename);
  return {};
}

} // namespace scream
