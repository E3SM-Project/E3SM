#ifndef SCREAM_FACTORY_HPP
#define SCREAM_FACTORY_HPP

#include <string>
#include <map>
#include <share/error_defs.hpp>

namespace scream 
{

template<typename AbstractProduct,
         typename KeyType,
         typename... ConstructorArgs>
class Factory
{
public:

  typedef Factory<AbstractProduct,KeyType,ConstructorArgs...> factory_type;

  typedef KeyType           key_type;
  typedef AbstractProduct   obj_type;
  typedef obj_type*         obj_ptr_type;

  typedef obj_ptr_type (*creator_type) (const ConstructorArgs... args);
  typedef std::map<key_type,creator_type> register_type;

  static factory_type& instance ()
  {
    static factory_type factory;

    return factory;
  }

  // Register a creator and returns true if it was successfully registered
  bool register_product (const key_type& key,
                         const creator_type& creator,
                         const bool replace_if_found = false);

  // Creates a concrete object using the proper creator
  obj_ptr_type create (const key_type& label, ConstructorArgs&& ...args);

private:
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
         typename... ConstructorArgs>
bool Factory<AbstractProduct,KeyType,ConstructorArgs...>::
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
         typename... ConstructorArgs>
typename Factory<AbstractProduct,KeyType,ConstructorArgs...>::obj_ptr_type
Factory<AbstractProduct,KeyType,ConstructorArgs...>::
create (const key_type& key,ConstructorArgs&& ...args)
{
  // Check that the factory is not empty.
  // Note: this check is redundant, since, if negative, the next one would fail too.
  //       However, this check is symptomatic of a larger problem (not registering product)
  //       than simply not finding the requested one (perhaps because of spelling).
  error::runtime_check(m_register.size()>0, "[Factory] Error! There are no products registered in the factory. Did you forget to call 'register_product'?\n");

  auto it = m_register.find(key);

  // Check that the requested product is registered
  error::runtime_check(it!=m_register.end(), "[Factory] Error! The key '" + key + "' is not associated to any registered product. Did you forget to register it?\n");

  return (*it->second)(std::forward<ConstructorArgs>(args)...);
}

} // namespace scream 

#endif // SCREAM_FACTORY_HPP
