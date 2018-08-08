/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_LOOPS_UTILS_HPP
#define HOMMEXX_LOOPS_UTILS_HPP

#include "Types.hpp"

namespace Homme {

// Source: https://stackoverflow.com/a/7185723
// Used for iterating over a range of integers
// With this, you can write
// for(int i : int_range(start, end))
template <typename ordered_iterable> class Loop_Range {

public:
  class iterator {
    friend class Loop_Range;

  public:
    KOKKOS_INLINE_FUNCTION
    ordered_iterable operator*() const { return i_; }

    KOKKOS_INLINE_FUNCTION
    const iterator &operator++() {
      ++i_;
      return *this;
    }

    KOKKOS_INLINE_FUNCTION
    iterator operator++(int) {
      iterator copy(*this);
      ++i_;
      return copy;
    }

    KOKKOS_INLINE_FUNCTION
    bool operator==(const iterator &other) const { return i_ == other.i_; }
    KOKKOS_INLINE_FUNCTION
    bool operator!=(const iterator &other) const { return i_ != other.i_; }

  protected:
    KOKKOS_INLINE_FUNCTION
    constexpr iterator(ordered_iterable start) : i_(start) {}

  private:
    ordered_iterable i_;
  };

  KOKKOS_INLINE_FUNCTION
  constexpr iterator begin() const { return begin_; }

  KOKKOS_INLINE_FUNCTION
  constexpr iterator end() const { return end_; }

  KOKKOS_INLINE_FUNCTION
  constexpr int iterations() const { return *end_ - *begin_; }

  KOKKOS_INLINE_FUNCTION
  constexpr Loop_Range(ordered_iterable begin, ordered_iterable end)
      : begin_(begin), end_(end) {}

private:
  iterator begin_;
  iterator end_;
};

} // namespace Homme

#endif // HOMMEXX_LOOPS_UTILS_HPP
