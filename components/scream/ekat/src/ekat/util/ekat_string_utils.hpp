#ifndef EKAT_STRING_UTILS_HPP
#define EKAT_STRING_UTILS_HPP

#include <string>
#include <vector>

/*
 * A set of utilities for string manipulation
 *
 * This header (and its corresponding cpp file) contain two sets of things:
 *  - A class for case insensitive string. This class inherits from std::string,
 *    so *all* the functionalities of std::string are available. Additionally,
 *    when comparing two strings (with ==, !=, <, or <=), if one of the two
 *    operands is a CaseInsensitiveString object, then we perform a case
 *    insensitive comparison.
 *    Use this class if you want to allow all possible case styles for some inputs.
 *  - Some utility functions to manipulate std::string objects, such as
 *    removing leading/trailing whitespaces, split a string into substrings
 *    at every occurrence of a given char, and more.
 */

namespace scream {
namespace util {

// Strip character c from input string
void strip (std::string& str, const char c);

// Split a string at every occurrence of given delimiter char
std::vector<std::string> split(const std::string& str, const char del);

// Trim leading/trailing whitespaces.
std::string trim (const std::string& s);

// Small utility that cats a space and an integer to an input string.
std::string strint (const std::string& s, const int i);

// A no-overhead class that inherits from std::string, which we only
// use to trigger different behavior in the ==,!=,<,<= operators.
class CaseInsensitiveString final : public std::string
{
public:
  template<typename... Args>
  CaseInsensitiveString (Args... args)
   : std::string(args...)
  {}

  virtual ~CaseInsensitiveString () = default;
};

// Case-insensitive comparison functions
bool caseInsensitiveEqualString (const std::string& s1, const std::string& s2);
bool caseInsensitiveLessString (const std::string& s1, const std::string& s2);
bool caseInsensitiveLessEqualString (const std::string& s1, const std::string& s2);

// Overloads of comparison operators, which use the routines above if at least one
// of the two inputs is indeed a CaseInsensitiveString
template<typename S1, typename S2>
typename std::enable_if<
      std::is_same<S1,CaseInsensitiveString>::value ||
      std::is_same<S2,CaseInsensitiveString>::value,
      bool>::type
operator== (const S1& s1, const S2& s2) {
  return caseInsensitiveEqualString(s1,s2);
}

template<typename S1, typename S2>
typename std::enable_if<
      std::is_same<S1,CaseInsensitiveString>::value ||
      std::is_same<S2,CaseInsensitiveString>::value,
      bool>::type
operator!= (const S1& s1, const S2& s2) {
  return ! (s1==s2);
}

template<typename S1, typename S2>
typename std::enable_if<
      std::is_same<S1,CaseInsensitiveString>::value ||
      std::is_same<S2,CaseInsensitiveString>::value,
      bool>::type
operator< (const S1& s1, const S2& s2) {
  return caseInsensitiveLessString(s1,s2);
}

template<typename S1, typename S2>
typename std::enable_if<
      std::is_same<S1,CaseInsensitiveString>::value ||
      std::is_same<S2,CaseInsensitiveString>::value,
      bool>::type
operator<= (const S1& s1, const S2& s2) {
  return caseInsensitiveLessEqualString(s1,s2);
}

} // namespace util
} // namespace scream

#endif // EKAT_STRING_UTILS_HPP
