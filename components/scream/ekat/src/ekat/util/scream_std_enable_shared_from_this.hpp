#ifndef SCREAM_STD_ENABLE_SHARED_FROM_THIS_HPP
#define SCREAM_STD_ENABLE_SHARED_FROM_THIS_HPP


#include <memory>

#include "ekat/scream_assert.hpp"

namespace scream {

namespace util {

// ================= std::enable_shared_from_this =============== //

/*
 * README
 *
 * This class already exists in c++11. However, there's an important
 * feature missing in c++11, which will be available only in c++17:
 * when calling the shared_from_this() method, it is crucial that the
 * object is indeed owned by a shared_ptr. If not, undef behavior will
 * ensue. This check is (implicitly) present in c++17 (but not before),
 * guaranteed by the fact that the class stores a weak_ptr, which is
 * assigned to during the shared_ptr constructor. Since we cannot change the
 * content of std::shared_ptr's constructor, to achieve the desired result,
 * one needs to manually set the self pointer upon shared_ptr construction.
 * I thought about rolling our own make_shared routine, which would do that,
 * but there is a problem. The way to check if a newly constructed shared_ptr
 * points to something that inherits from enable_shared_from_this<T>, is to
 * try to cast to it. But to do so, you need to know what T is.
 * Example: say you have a class A inheriting from
 * enable_shared_from_this<A>. Then you have B inheriting from A. When
 * you try to do make_shared<B>(args...), we have no knowledge of the
 * base classes of B (there's no reflection mechanism yet in c++11 or c++17,
 * for that matter). All we can check is whether class B inherits from
 * enable_shared_from_this<B>, which would fail. Adding an extra template
 * parameter to make_shared would make the signature of the function
 * too different from the one in the std namespace, making the switch
 * from util::make_shared to std::make_shared in the future too complicated.
 * Long story short: use util::enable_shared_from_this, but be diligent,
 * and set it up correctly when you create a shared_ptr (assuming you want
 * to exploit the feature of enable_shared_from_this).
 */

template<typename T>
class enable_shared_from_this
{
public:

  enable_shared_from_this () = default;
  virtual ~enable_shared_from_this () = default;

  std::shared_ptr<T> shared_from_this () {
    return m_self.lock();
  }

  std::weak_ptr<T> weak_from_this () {
    return m_self;
  }

  // Make sure the input pointer is indeed a pointer to self
  void setSelfPointer (const std::shared_ptr<T>& ptr) {
    scream_require_msg(ptr.get()==this, "Error! Cannot set self pointer to something different from 'this'.\n");
    m_self = ptr;
  }

protected:

  std::weak_ptr<T>  m_self;
};

} // namespace util

} // namespace scream

#endif // SCREAM_STD_ENABLE_SHARED_FROM_THIS_HPP
