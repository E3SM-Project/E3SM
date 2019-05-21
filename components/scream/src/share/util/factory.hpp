#ifndef SCREAM_FACTORY_HPP
#define SCREAM_FACTORY_HPP

#include <string>
#include <map>

#include "share/scream_assert.hpp"

namespace scream
{

namespace util
{

template<typename AbstractProduct,
         typename KeyType,
         typename PointerType,
         typename... ConstructorArgs>
class Factory
{
public:

  using factory_type  = Factory<AbstractProduct,KeyType,PointerType,ConstructorArgs...>;
  using key_type      = KeyType;
  using obj_type      = AbstractProduct;
  using obj_ptr_type  = PointerType;
  using creator_type  = obj_ptr_type (*) (const ConstructorArgs... args);
  using register_type = std::map<key_type,creator_type>;

  // Make sure that obj_ptr_type is indeed some sort of pointer to obj_type
  using deref_pointer_type = typename std::remove_reference<decltype(*std::declval<PointerType>())>::type;
  static_assert (std::is_same<AbstractProduct,deref_pointer_type>::value,
                 "[Factory] Error! Template parameter for AbstractProduct and PointerType are not compatible.");

  static factory_type& instance ()
  {
    static factory_type factory;

    return factory;
  }

  // For debugging, inspect if the size of the register in the factory
  size_t register_size () const { return m_register.size(); }

  // Register a creator and returns true if it was successfully registered
  bool register_product (const key_type& key,
                         const creator_type& creator,
                         const bool replace_if_found = false);

  // Creates a concrete object using the proper creator
  obj_ptr_type create (const key_type& label, ConstructorArgs&& ...args);

private:
  std::string print_registered_products () const;

  Factory () = default;
  Factory (const factory_type&) = delete;
  factory_type& operator= (const factory_type&) = delete;
  Factory (factory_type&&) = delete;
  factory_type& operator= (factory_type&&) = delete;

  register_type    m_register;
};

// ========================== IMPLEMENTATION ======================== //

template<typename AbstractProduct,
         typename KeyType,
         typename PointerType,
         typename... ConstructorArgs>
bool Factory<AbstractProduct,KeyType,PointerType,ConstructorArgs...>::
register_product (const key_type& key,
                  const creator_type& creator,
                  const bool replace_if_found)
{
  auto it = m_register.find(key);

  if (it==m_register.end() || replace_if_found )
  {
    m_register[key] = creator;
    return true;
  }
  return false;
}

template<typename AbstractProduct,
         typename KeyType,
         typename PointerType,
         typename... ConstructorArgs>
PointerType Factory<AbstractProduct,KeyType,PointerType,ConstructorArgs...>::
create (const key_type& key,ConstructorArgs&& ...args)
{
  // Check that the factory is not empty.
  // Note: this check is redundant, since, if negative, the next one would fail too.
  //       However, this check is symptomatic of a larger problem (not registering product)
  //       than simply not finding the requested one (perhaps because of spelling).
  scream_require_msg(m_register.size()>0,
                     "[Factory] Error! There are no products registered in the factory.\n"
                     "                 Did you forget to call 'register_product'?\n");

  auto it = m_register.find(key);

  // Check that the requested product is registered
  scream_require_msg(it!=m_register.end(),
                     "[Factory] Error! The key '" + key + "' is not associated to any registered product.\n"
                     "                 The list of registered product is: " + print_registered_products() + "\n"
                     "                 Did you forget to register it?\n");

  return (*it->second)(std::forward<ConstructorArgs>(args)...);
}

template<typename AbstractProduct,
         typename KeyType,
         typename PointerType,
         typename... ConstructorArgs>
std::string Factory<AbstractProduct,KeyType,PointerType,ConstructorArgs...>::
print_registered_products () const {
  // This routine simply puts the products name in a string, as "name1, name2, name3,..., name N"
  std::string str;
  if (m_register.size()==0) {
    return str;
  }
  auto next = m_register.begin();
  auto it = m_register.begin();
  for (++next; next!=m_register.end(); ++it, ++next) {
    str += it->first + ", ";
  }
  str += next->first;

  return str;
}

} // namespace util
} // namespace scream 

#endif // SCREAM_FACTORY_HPP
