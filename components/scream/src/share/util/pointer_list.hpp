#ifndef SCREAM_POINTER_LIST
#define SCREAM_POINTER_LIST

#include <vector>

namespace scream {

// A pointer_list is a linearly traversible list of pointers that provides
// iterators that automatically double-dereference their referents. Under the
// covers, a pointer_list is just a vector of the given pointer type.
template <typename PointerType, typename ValueType>
class pointer_list final {

public:

  class iterator final {
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = ValueType;
    using difference_type = std::ptrdiff_t;
    using pointer = PointerType;
    using reference = ValueType&;
    using base_iter_type = typename std::vector<PointerType>::iterator;

    explicit iterator(base_iter_type iter) : m_iter(iter) {}
    iterator(const base_iter_type& iter) : m_iter(iter.m_iter) {}
    iterator& operator=(const iterator& iter) {
      m_iter = iter.m_iter;
      return *this;
    }

    iterator& operator++() { m_iter++; return *this;}
    bool operator==(iterator other) const {return m_iter == other.m_iter;}
    bool operator!=(iterator other) const {return !(*this == other);}
    pointer operator->() const {return *m_iter;}
    const reference& operator*() const {return **m_iter;}
  private:
    base_iter_type m_iter;
  };

  class const_iterator final {
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = const ValueType;
    using difference_type = std::ptrdiff_t;
    using pointer = PointerType;
    using reference = const ValueType&;
    using base_iter_type = typename std::vector<PointerType>::const_iterator;

    explicit const_iterator(base_iter_type iter) : m_iter(iter) {}
    const_iterator(const base_iter_type& iter) : m_iter(iter.m_iter) {}
    const_iterator& operator=(const const_iterator& iter) {
      m_iter = iter.m_iter;
      return *this;
    }

    iterator& operator++() { m_iter++; return *this;}
    bool operator==(iterator other) const {return m_iter == other.m_iter;}
    bool operator!=(iterator other) const {return !(*this == other);}
    pointer operator->() const {return *m_iter;}
    const reference& operator*() const {return **m_iter;}
  private:
    base_iter_type m_iter;
  };

  // These iterators provide access to the list's contents.
  iterator begin();
  const_iterator begin() const;
  iterator end();
  const_iterator end() const;

private:

  std::vector<PointerType> m_list;
};

} // end namespace scream

#endif
