#ifndef SCREAM_STD_UTILS_HPP
#define SCREAM_STD_UTILS_HPP

#include <algorithm>
#include <memory>
#include <vector>
#include <string>

namespace scream {
namespace util {

// This function returns an iterator which is of the same type of c.begin()
template<typename ContainerType, typename T>
auto find (const ContainerType& c, T&& value) -> decltype(c.begin()) {
  return std::find(c.begin(),c.end(),value);
}

// Note: in C++20, both std::set and std::map will have the 'contains' method.
template<typename ContainerType, typename T>
bool contains (const ContainerType& c, T&& value) {
  return std::find(c.begin(),c.end(),value) != c.end();
}

template<typename ContainerType, typename T>
bool erase (ContainerType& c, const T& value) {
  auto it = std::find(c.begin(),c.end(),value);
  if (it!=c.end()) {
    c.erase(it);
    return true;
  }
  return false;
}

// Shortcuts to avoid using long iterators syntax (like v.begin() and v.end())
template<typename ContainerType, typename T>
int count (const ContainerType& c, const T& value) {
  return std::count(c.begin(), c.end(), value);
}

template<typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
  const int n = v.size();
  if (n==0) {
    return out;
  }

  for (int i=0; i<n-1; ++i) {
    out << v[i] << std::string(" ");
  }
  out << v.back();

  return out;
}

} // namespace util

// A set of weak_ptr would not compile, due to the lack of operator<.
// To overcome this, one could add a Compare type to the set template
// arguments. Or, more easily, overload operator<...
// NOTE: this cannot be in the util namespace, or else the compiler won't detect
//       it, when T is declared in the scream namespace

template<typename T>
inline bool operator< (const std::weak_ptr<T>& p, const std::weak_ptr<T>& q)
{
  return p.owner_before(q);
}

} // namespace scream

#endif // SCREAM_STD_UTILS_HPP
