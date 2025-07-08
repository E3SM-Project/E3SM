#ifndef SCREAM_PHYSICS_TEST_READ_WRITE_HELPERS_HPP
#define SCREAM_PHYSICS_TEST_READ_WRITE_HELPERS_HPP

#include <fstream>
#include <vector>
#include <type_traits>

namespace scream {
namespace impl {

template <typename T>
void read_scalars(std::ifstream& ifile, T& data)
{
  ifile.read(reinterpret_cast<char*>(&data),sizeof(data));
}

template <typename T>
void read_scalars(std::ifstream& ifile, std::vector<T>& data)
{
  for (auto& entry : data) {
    read_scalars(ifile,entry);
  }
}

template <typename T, int N>
void read_scalars(std::ifstream& ifile, T (&data)[N]) {
  for (int i=0; i<N; ++i) {
    read_scalars(ifile,data[i]);
  }
}

template <typename T, typename I>
std::enable_if_t<std::is_integral<I>::value>
read_scalars(std::ifstream& ifile, T* const data, I N) {
  for (I i=0; i<N; ++i) {
    read_scalars(ifile,data[i]);
  }
}

template <typename T, typename... S>
void read_scalars(std::ifstream& ifile, T& data, S&... tail)
{
  read_scalars(ifile,data);
  read_scalars(ifile,tail...);
}

template <typename T>
void write_scalars(std::ofstream& ofile, const T& data)
{
  ofile.write(reinterpret_cast<const char*>(&data),sizeof(data));
}

template <typename T>
void write_scalars(std::ofstream& ofile, const std::vector<T>& data)
{
  for (const auto& entry : data) {
    write_scalars(ofile,entry);
  }
}

template <typename T, int N>
void write_scalars(std::ofstream& ofile, const T (&data)[N]) {
  for (int i=0; i<N; ++i) {
    write_scalars(ofile,data[i]);
  }
}

template <typename T, typename I>
std::enable_if_t<std::is_integral<I>::value>
write_scalars(std::ofstream& ofile, const T* const data, I N) {
  for (I i=0; i<N; ++i) {
    write_scalars(ofile,data[i]);
  }
}


template <typename T, typename... S>
void write_scalars(std::ofstream& ofile, const T& data, const S&... tail)
{
  write_scalars(ofile,data);
  write_scalars(ofile,tail...);
}

} // namespace impl
} // namespace scream

#endif // SCREAM_PHYSICS_TEST_READ_WRITE_HELPERS_HPP
