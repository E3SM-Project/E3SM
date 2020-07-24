#include <ekat/util/ekat_string_utils.hpp>
#include <algorithm>
#include <sstream>

namespace scream {
namespace util {

void strip (std::string& str, const char c) {
  auto new_end = std::remove(str.begin(),str.end(),c);
  str.erase(new_end,str.end());
}

std::vector<std::string> split(const std::string& str, const char del) {
  std::vector<std::string> blocks;

  auto start = 0;
  auto pos = str.find(del);
  while (pos!=std::string::npos) {
    blocks.push_back(str.substr(start,pos-start));
    start = pos + 1;
    pos = str.find(del,start);
  }

  // Don't forget to add the substring from the last occurrence of 'del' (if any) to the end of str
  blocks.push_back(str.substr(start));
  return blocks;
}

std::string trim (const std::string& s) {
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

std::string strint (const std::string& s, const int i) {
  std::stringstream ss;
  ss << s << " " << i;
  return ss.str();
}

std::string upper_case (const std::string& s) {
  std::string s_up = s;
  std::transform(s_up.begin(), s_up.end(), s_up.begin(),
                 [](unsigned char c)->char { return std::toupper(c); }
                );
  return s_up;
}

bool caseInsensitiveEqualString (const std::string& s1, const std::string& s2) {
  auto charComp = [](const char c1, const char c2)->bool{
    return c1==c2 || std::toupper(c1)==std::toupper(c2);
  };
  return s1.size()==s2.size() &&
         std::equal(s1.begin(),s1.end(),s2.begin(),charComp);
}

bool caseInsensitiveLessString (const std::string& s1, const std::string& s2) {
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

bool caseInsensitiveLessEqualString (const std::string& s1, const std::string& s2) {
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

} // namespace util
} // namespace scream
