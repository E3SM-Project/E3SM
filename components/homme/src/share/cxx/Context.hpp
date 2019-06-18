/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_CONTEXT_HPP
#define HOMMEXX_CONTEXT_HPP

#include <string>
#include <map>
#include <memory>

#include "ErrorDefs.hpp"
#include "utilities/StdMeta.hpp"

namespace Homme {

/* A Context manages resources that are morally singletons. Context is
 * meant to have two roles. First, a Context singleton is the only singleton in
 * the program. Second, a context need not be a singleton, and each Context
 * object can have different Elements, ReferenceElement, etc., objects. (That
 * probably isn't needed, but Context immediately supports it.)
 *
 * Finally, Context has two singleton functions: singleton(), which returns
 * Context&, and finalize_singleton(). The second is called in a unit test exe
 * main before Kokkos::finalize().
 */
class Context {
public:

  // Getters for a managed object.
  template<typename ConcreteType>
  bool has () const;

  // Setters for a managed object.
  template<typename ConcreteType, typename... Args>
  ConcreteType& create (Args&&... args);

  // More relaxed than create, it won't throw if the
  // object already exists, and simply return it.
  // NOTE: this is allows more flexibility in the cxx-f90
  //       interfaces, in that you have more freedom in
  //       the order of certain function calls.
  //       With 'create', you have to either check if an
  //       object was already created (and create it if it wasn't),
  //       or make sure that the function that calls 'create<T>'
  //       is *always* invoked *before* those that call 'get<T>'.
  template<typename ConcreteType, typename... Args>
  ConcreteType& create_if_not_there (Args&&... args);

  // Creates a concrete type and sets a second entry in the map, with
  // BaseType's name as key, so one can access the concrete type either
  // with the base or derived type name.
  template<typename BaseType, typename ConcreteType, typename... Args>
  ConcreteType& create_with_base (Args&&... args) {
    ConcreteType& c = create<ConcreteType>(args...);
    alias<BaseType,ConcreteType>();
    return c;
  }

  // Getters for a managed object.
  template<typename ConcreteType>
  ConcreteType& get () const;

  template<typename ConcreteType>
  std::shared_ptr<ConcreteType> get_ptr () const;

  // Exactly one singleton.
  static Context& singleton();

  static void finalize_singleton();
private:
  // Sets an entry in the map for a BaseType given the DerivedType.
  // This allows to set a specific derived type first, then set an
  // alias of it using the base type, so that other files can access
  // the concrete derived type using the base type.
  template<typename BaseType, typename DerivedType>
  void alias ();

  std::map<std::string,Homme::any> m_members;

  // Clear the objects Context manages.
  void clear();
};

// ==================== IMPLEMENTATION =================== //

template<typename ConcreteType>
bool Context::has () const {
  const std::string& name = typeid(ConcreteType).name();
  auto it = m_members.find(name);
  return it!=m_members.end();
}

template<typename ConcreteType, typename... Args>
ConcreteType& Context::create_if_not_there (Args&&... args) {
  // This is needed for emplacing a type whose constructor takes no arguments.
  // We could do emplace(name,ConcreteType()), but then we would be assuming
  // that ConcreteType *has* a move constructor. This implementation here is
  // probably the most cumbersome, but also the safest.
  auto it_bool = m_members.emplace(typeid(ConcreteType).name(),Homme::any());
  if (it_bool.second) {
    // We created it, so init it.
    it_bool.first->second.reset<ConcreteType>(args...);
  }
  return *any_ptr_cast<ConcreteType>(it_bool.first->second);
}

template<typename ConcreteType, typename... Args>
ConcreteType& Context::create (Args&&... args) {
  Errors::runtime_check(!has<ConcreteType>(),
                        "Error! An object for the concrete type " + std::string(typeid(ConcreteType).name()) +
                        " is already stored. The 'Context' class does not allow overwriting or duplicates.\n");
  // This is needed for emplacing a type whose constructor takes no arguments.
  // We could do emplace(name,ConcreteType()), but then we would be assuming
  // that ConcreteType *has* a move constructor. This implementation here is
  // probably the most cumbersome, but also the safest.
  auto it_bool = m_members.emplace(typeid(ConcreteType).name(),Homme::any());
  Errors::runtime_check(it_bool.second, "Error! Something went wrong when inserting a new element in the context. "
                                        "This is an internal error. Please, contact developers.\n", -1);
  it_bool.first->second.reset<ConcreteType>(args...);

  return *any_ptr_cast<ConcreteType>(it_bool.first->second);
}

template<typename BaseType, typename DerivedType>
void Context::alias () {
  static_assert(std::is_base_of<BaseType,DerivedType>::value, "Error! Trying to use 'alias' with types not compatible (BaseType must be a base type of DerivedType).\n");
  Errors::runtime_check(has<DerivedType>(), "Error! An object with the concrete type " + std::string(typeid(DerivedType).name()) +
                                              " is not yet stored. In order to 'alias' it with a base type, you need to set the derived type first.\n");

  auto it_bool = m_members.emplace(typeid(BaseType).name(),Homme::any());
  Errors::runtime_check(it_bool.second, "Error! Something went wrong when inserting a new element in the context. "
                                        "This is an internal error. Please, contact developers.\n");
  auto derived_ptr = get_ptr<DerivedType>();

  it_bool.first->second.reset_ptr<BaseType>(derived_ptr);
}

template<typename ConcreteType>
ConcreteType& Context::get () const {
  return *get_ptr<ConcreteType>();
}

template<typename ConcreteType>
std::shared_ptr<ConcreteType> Context::get_ptr() const {
  const std::string& name = typeid(ConcreteType).name();
  auto it = m_members.find(name);
  Errors::runtime_check(it!=m_members.end(), "Error! Context member '" + name + "' not found.\n", -1);

  return any_ptr_cast<ConcreteType>(it->second);
}

} // namespace Homme

#endif // HOMMEXX_CONTEXT_HPP
