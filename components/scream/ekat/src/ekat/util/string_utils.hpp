#ifndef SCREAM_STRING_UTILS_HPP
#define SCREAM_STRING_UTILS_HPP

#include <string>
#include <sstream>
#include <algorithm>

namespace scream {
namespace util {

// Trim leading/trailing whitespaces
inline std::string trim (const std::string& s) {
  if (s=="") {
    return s;
  }
  int nl = 0;
  while (s[nl] == ' ') {
    ++nl;
  }
  int size = s.size();
  int nt = 0;
  while (nt<size && s[size-nt-1] == ' ') {
    ++nt;
  }

  return s.substr(nl,size-nl-nt);
}

// Small utility that cats a space and an integer to an input string.
inline std::string strint (const std::string& s, const int i) {
  std::stringstream ss;
  ss << s << " " << i;
  return ss.str();
}

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
bool caseInsensitiveLessEqualString (const std::string& s1, const std::string& s2);
bool caseInsensitiveLessEqualString (const std::string& s1, const std::string& s2);

inline bool caseInsensitiveEqualString (const std::string& s1, const std::string& s2) {
  auto charComp = [](const char c1, const char c2)->bool{
    return c1==c2 || std::toupper(c1)==std::toupper(c2);
  };
  return s1.size()==s2.size() &&
         std::equal(s1.begin(),s1.end(),s2.begin(),charComp);
}

inline bool caseInsensitiveLessString (const std::string& s1, const std::string& s2) {
  auto charCompLess = [](const char c1, const char c2)->bool{
    return std::toupper(c1)<std::toupper(c2);
  };
  auto charCompEq = [](const char c1, const char c2)->bool{
    return std::toupper(c1)==std::toupper(c2);
  };
  for (auto it1=s1.begin(),it2=s2.begin(); it1!=s1.end() && it2!=s2.end(); ++it1,++it2) {
    if (charCompLess(*it1,*it2)) {
      return true;
    } else if (!charCompEq(*it1,*it2)) {
      return false;
    }
  }
  return s1.size()<s2.size();
}

inline bool caseInsensitiveLessEqualString (const std::string& s1, const std::string& s2) {
  auto charCompLess = [](const char c1, const char c2)->bool{
    return std::toupper(c1)<std::toupper(c2);
  };
  auto charCompEq = [](const char c1, const char c2)->bool{
    return std::toupper(c1)==std::toupper(c2);
  };
  for (auto it1=s1.begin(),it2=s2.begin(); it1!=s1.end() && it2!=s2.end(); ++it1,++it2) {
    if (charCompLess(*it1,*it2)) {
      return true;
    } else if (!charCompEq(*it1,*it2)) {
      return false;
    }
  }
  return s1.size()<=s2.size();
}

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

#endif // SCREAM_STRING_UTILS_HPP
