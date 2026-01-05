

# File data\_view.hpp

[**File List**](files.md) **>** [**common**](dir_03a82009d675450cd7b26f7887b718e6.md) **>** [**src**](dir_cf65ed25a3cb4ff7a3950496111c5548.md) **>** [**data\_view.hpp**](data__view_8hpp.md)

[Go to the documentation of this file](data__view_8hpp.md)


```C++


#ifndef EMULATOR_DATA_VIEW_HPP
#define EMULATOR_DATA_VIEW_HPP

#include <cstddef>
#include <vector>

namespace emulator {

enum class DataLayout {
  ROW_MAJOR,    
  COLUMN_MAJOR, 
  UNKNOWN       
};

template <typename T> class DataView {
public:
  DataView() : m_data(nullptr), m_size(0), m_layout(DataLayout::UNKNOWN) {}

  DataView(T *data, std::size_t size, DataLayout layout = DataLayout::ROW_MAJOR)
      : m_data(data), m_size(size), m_layout(layout) {}

  DataView(const T *data, std::size_t size,
           DataLayout layout = DataLayout::ROW_MAJOR)
      : m_data(const_cast<T *>(data)), m_size(size), m_layout(layout) {}

  explicit DataView(std::vector<T> &vec,
                    DataLayout layout = DataLayout::ROW_MAJOR)
      : m_data(vec.data()), m_size(vec.size()), m_layout(layout) {}

  explicit DataView(const std::vector<T> &vec,
                    DataLayout layout = DataLayout::ROW_MAJOR)
      : m_data(const_cast<T *>(vec.data())), m_size(vec.size()),
        m_layout(layout) {}

  T *data() { return m_data; }
  const T *data() const { return m_data; }

  std::size_t size() const { return m_size; }

  bool empty() const { return m_size == 0 || m_data == nullptr; }

  DataLayout layout() const { return m_layout; }

  T &operator[](std::size_t idx) { return m_data[idx]; }
  const T &operator[](std::size_t idx) const { return m_data[idx]; }

  T *begin() { return m_data; }
  T *end() { return m_data + m_size; }
  const T *begin() const { return m_data; }
  const T *end() const { return m_data + m_size; }

private:
  T *m_data;
  std::size_t m_size;
  DataLayout m_layout;
};

using DoubleView = DataView<double>;
using FloatView = DataView<float>;

} // namespace emulator

#endif // EMULATOR_DATA_VIEW_HPP
```


