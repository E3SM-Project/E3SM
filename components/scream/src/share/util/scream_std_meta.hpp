#ifndef SCREAM_STD_META_HPP
#define SCREAM_STD_META_HPP

/*
 *  Emulating c++1z features
 *
 *  This file contains some features that belongs to the std namespace,
 *  and that are likely (if not already scheduled) to be in future standards.
 *  However, we have limited access to c++1z standard on target platforms,
 *  so we can only rely on c++11 features. Everything else that we may need,
 *  we try to emulate here.
 */

#include <memory>

#include <share/error_defs.hpp>

namespace scream {

namespace util {

// =============== type_traits utils ============== //

// <type_traits> has remove_all_extents but lacks the analogous
// remove_all_pointers
template<typename T>
struct remove_all_pointers {
  using type = T;
};
template<typename T>
struct remove_all_pointers<T*> {
  using type = typename remove_all_pointers<T>::type;
};

// std::remove_const does not remove the leading const cv from
// the type <const T*>. Indeed, remove_const can only remove the
// const cv from the pointer, not the pointee.
template<typename T>
struct remove_all_consts : std::remove_const<T> {};

template<typename T>
struct remove_all_consts<T*> {
  using type = typename remove_all_consts<T>::type*;
};

// ================ std::any ================= //

namespace impl {

class holder_base {
public:
  virtual ~holder_base() = default;

  virtual const std::type_info& type () const = 0;
};

template <typename HeldType>
class holder : public holder_base {
public:

  template<typename... Args>
  holder (const Args&... args) {
    m_value = std::make_shared<HeldType>(args...);
  }

  template<typename T>
  holder (typename std::enable_if<std::is_base_of<HeldType,T>::value,std::shared_ptr<T>>::type ptr) {
    m_value = ptr;
  }

  const std::type_info& type () const { return typeid(HeldType); }

  HeldType& value () { return *m_value; }
  std::shared_ptr<HeldType> ptr () const { return m_value; }
private:
  // Note: we store a shared_ptr rather than a HeldType directly,
  //       since we may store multiple copies of the concrete object.
  //       Since we don't know if the actual object is copiable, we need
  //       to store copiable pointers (hence, cannot use unique_ptr).
  std::shared_ptr<HeldType> m_value;
};

} // namespace impl

class any {
public:

  any () = default;

  template<typename T, typename... Args>
  void reset (Args... args) {
    m_content.reset( new impl::holder<T>(args...) );
  }

  impl::holder_base& content () const { 
    error::runtime_check(static_cast<bool>(m_content), "Error! Object not yet initialized.\n", -1);
    return *m_content;
  }

  impl::holder_base* content_ptr () const { 
    return m_content.get();
  }

private:

  std::shared_ptr<impl::holder_base> m_content;
};

template<typename ConcreteType>
bool any_is_type (const any& src) {
  const auto& src_type = src.content().type();
  const auto& req_type = typeid(ConcreteType);

  return src_type==req_type;
}

template<typename ConcreteType>
ConcreteType& any_cast (any& src) {
  const auto& src_type = src.content().type();
  const auto& req_type = typeid(ConcreteType);

  error::runtime_check(src_type==req_type, std::string("Error! Invalid cast requested, from '") + src_type.name() + "' to '" + req_type.name() + "'.\n", -1);

  impl::holder<ConcreteType>* ptr = dynamic_cast<impl::holder<ConcreteType>*>(src.content_ptr());

  error::runtime_check(ptr!=nullptr, "Error! Failed dynamic_cast during any_cast. This is an internal problem, please, contact developers.\n", -1);

  return ptr->value();
}

template<typename ConcreteType>
const ConcreteType& any_cast (const any& src) {
  const auto& src_type = src.content().type();
  const auto& req_type = typeid(ConcreteType);

  error::runtime_check(src_type==req_type, std::string("Error! Invalid cast requested, from '") + src_type.name() + "' to '" + req_type.name() + "'.\n", -1);

  impl::holder<ConcreteType>* ptr = dynamic_cast<impl::holder<ConcreteType>*>(src.content_ptr());

  error::runtime_check(ptr!=nullptr, "Error! Failed dynamic_cast during any_cast. This is an internal problem, please, contact developers.\n", -1);

  return ptr->value();
}

template<typename ConcreteType>
std::shared_ptr<ConcreteType> any_ptr_cast (any& src) {
  const auto& src_type = src.content().type();
  const auto& req_type = typeid(ConcreteType);

  error::runtime_check(src_type==req_type, std::string("Error! Invalid cast requested, from '") + src_type.name() + "' to '" + req_type.name() + "'.\n", -1);

  impl::holder<ConcreteType>* ptr = dynamic_cast<impl::holder<ConcreteType>*>(src.content_ptr());

  error::runtime_check(ptr!=nullptr, "Error! Failed dynamic_cast during any_cast. This is an internal problem, please, contact developers.\n", -1);

  return ptr->ptr();
}

} // namespace util

} // namespace scream

#endif // SCREAM_STD_META_HPP
