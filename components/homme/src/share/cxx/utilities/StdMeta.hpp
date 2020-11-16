#ifndef HOMMEXX_STD_META_HPP
#define HOMMEXX_STD_META_HPP

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

#include "ErrorDefs.hpp"

namespace Homme {

// ================ std::any ================= //

namespace Impl {

class holder_base {
public:
  virtual ~holder_base() = default;

  virtual const std::type_info& type () const = 0;
};

template <typename HeldType>
class holder : public holder_base {
public:

  holder () = default;

  virtual const std::type_info& type () const { return typeid(HeldType); }

  HeldType& value () { return *m_value; }
  std::shared_ptr<HeldType> ptr () const { return m_value; }

  template<typename DerivedFromHeldType>
  void reset_ptr (std::shared_ptr<DerivedFromHeldType> ptr) {
    static_assert(std::is_base_of<HeldType,DerivedFromHeldType>::value,
                  "Error! Attempting to set stored value to a pointer to a non-derived class.\n");
    m_value = ptr;
  }

  static holder<HeldType>* create_empty () { return new holder<HeldType>(); }

  template<typename... Args>
  static holder<HeldType>* create (Args... args) {
    auto ptr = create_empty();
    ptr->m_value = std::make_shared<HeldType>(args...);
    return ptr;
  }

private:
  // Note: we store a shared_ptr rather than a HeldType directly,
  //       since we may store multiple copies of the concrete object.
  //       Since we don't know if the actual object is copiable, we need
  //       to store copiable pointers (hence, cannot use unique_ptr).
  std::shared_ptr<HeldType> m_value;
};

} // namespace Impl

class any {
public:

  any () = default;

  template<typename T, typename... Args>
  void reset (Args&&... args) {
    m_content.reset( Impl::holder<T>::create(args...) );
  }

  template<typename BaseType, typename ConcreteType>
  void reset_ptr (std::shared_ptr<ConcreteType> ptr) {
    auto raw_ptr = Impl::holder<BaseType>::create_empty();
    raw_ptr->reset_ptr(ptr);
    m_content.reset(raw_ptr);
  }

  Impl::holder_base& content () const { 
    Errors::runtime_check(static_cast<bool>(m_content), "Error! Object not yet initialized.\n", -1);
    return *m_content;
  }

  Impl::holder_base* content_ptr () const { 
    return m_content.get();
  }

private:

  std::unique_ptr<Impl::holder_base> m_content;
};

template<typename ConcreteType>
ConcreteType& any_cast (any& src) {
  const auto& src_type = src.content().type();
  const auto& req_type = typeid(ConcreteType);

  Errors::runtime_check(src_type==req_type, std::string("Error! Invalid cast requested, from '") + src_type.name() + "' to '" + req_type.name() + "'.\n", -1);

  Impl::holder<ConcreteType>* ptr = dynamic_cast<Impl::holder<ConcreteType>*>(src.content_ptr());

  Errors::runtime_check(ptr!=nullptr, "Error! Failed dynamic_cast during any_cast. This is an internal problem, please, contact developers.\n", -1);

  return ptr->value();
}

template<typename ConcreteType>
std::shared_ptr<ConcreteType> any_ptr_cast (const any& src) {
  const auto& src_type = src.content().type();
  const auto& req_type = typeid(ConcreteType);

  Errors::runtime_check(src_type==req_type, std::string("Error! Invalid cast requested, from '") + src_type.name() + "' to '" + req_type.name() + "'.\n", -1);

  Impl::holder<ConcreteType>* ptr = dynamic_cast<Impl::holder<ConcreteType>*>(src.content_ptr());

  Errors::runtime_check(ptr!=nullptr, "Error! Failed dynamic_cast during any_cast. This is an internal problem, please, contact developers.\n", -1);

  return ptr->ptr();
}

} // namespace Homme

#endif // HOMMEXX_STD_META_HPP
