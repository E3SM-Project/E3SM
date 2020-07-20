#ifndef SCREAM_FACTORY_HPP
#define SCREAM_FACTORY_HPP

#include <string>
#include <sstream>
#include <map>

#include "ekat/scream_assert.hpp"
#include "ekat/util/scream_utils.hpp"

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

    // We rely on whatever the implementation can give us
    factory.m_factory_name = typeid(factory_type).name();

    return factory;
  }

  // For debugging, inspect if the size of the register in the factory
  size_t register_size () const { return m_register.size(); }

  // Return the name of the factory
  static const std::string& name () { return m_factory_name; }

  // Register a creator and returns true if it was successfully registered
  bool register_product (const key_type& key,
                         const creator_type& creator,
                         const bool replace_if_found = false);

  // Creates a concrete object using the proper creator
  obj_ptr_type create (const key_type& label, ConstructorArgs&& ...args);


  // Cleans up the factory
  // This can be useful during unit tests, where several tests in the same
  // execution will try to register the same product(s).
  void clean_up () { m_register.clear(); }

private:
  std::string print_registered_products () const {
    return print_registered_products_impl<check_overloads::StreamExists<key_type>::value>();
  }

  template<int OpStreamExists>
  typename std::enable_if<OpStreamExists==1,std::string>::type
  print_registered_products_impl () const;

  template<int OpStreamExists>
  typename std::enable_if<OpStreamExists==0,std::string>::type
  print_registered_products_impl () const {
    // We cannot use << to print out the keys in this factory, so we just print a message.
    return " Warning! Cannot print registered products in factory " + m_factory_name + ".\n"
           "          The type " + std::string(typeid(key_type).name()) + " does not overload the '<<' operator.\n"
           "          We have no idea how to print the registered products, so we just print this message. Sorry.\n";
  }

  Factory () = default;
  Factory (const factory_type&) = delete;
  factory_type& operator= (const factory_type&) = delete;
  Factory (factory_type&&) = delete;
  factory_type& operator= (factory_type&&) = delete;

  register_type         m_register;
  static std::string    m_factory_name;
};

template<typename AbstractProduct,
         typename KeyType,
         typename PointerType,
         typename... ConstructorArgs>
std::string Factory<AbstractProduct,KeyType,PointerType,ConstructorArgs...>::m_factory_name;

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
                     "[" + m_factory_name + "] Error!\n"
                     "        There are no products registered in the factory.\n"
                     "        Did you forget to call 'register_product'?\n");

  auto it = m_register.find(key);

  // Check that the requested product is registered
  scream_require_msg(it!=m_register.end(),
                     "[" + m_factory_name + "] Error!\n"
                     "        The key '" + key + "' is not associated to any registered product.\n"
                     "        The list of registered product is: " + print_registered_products() + "\n"
                     "        Did you forget to register it?\n");

  return (*it->second)(std::forward<ConstructorArgs>(args)...);
}

template<typename AbstractProduct,
         typename KeyType,
         typename PointerType,
         typename... ConstructorArgs>
template<int OpStreamExists>
typename std::enable_if<OpStreamExists==1,std::string>::type
Factory<AbstractProduct,KeyType,PointerType,ConstructorArgs...>::
print_registered_products_impl () const {
  // This routine simply puts the products name in a string, as "name1, name2, name3,..., name N"
  std::stringstream ss;

  if (m_register.size()==0) {
    return ss.str();
  }
  auto it = m_register.begin();
  ss << it->first;
  for (++it; it!=m_register.end(); ++it) {
    ss << ", " << it->first;
  }

  return ss.str();
}

} // namespace util
} // namespace scream 

#endif // SCREAM_FACTORY_HPP
