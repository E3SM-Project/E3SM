#ifndef SCREAM_MAP_KEY_ITERATOR_HPP
#define SCREAM_MAP_KEY_ITERATOR_HPP

#include <map>

namespace scream {

/*
 * An iterator that allows to iterate over the keys of a map only.
 *
 * The reason to have this class is that we can expose a way to
 * access the keys of a map without exposing the actual content
 * of the pair. This is very handy for map<K,std::shared_ptr<V>>,
 * where even if the map is const, iterating over it allows one
 * to modify the pointees.
 * By offering a key-iterator, one can iterate over the keys of
 * the map, without exposing the values in any way.
 * Note: this is similar to python's dictionary 'keys()' method.
 */

template<typename MapType>
class map_key_iterator;
template<typename MapType>
class map_key_const_iterator;

template<typename Key, typename Value>
class map_key_iterator<std::map<Key,Value>> final : public std::map<Key, Value>::iterator
{
public:
  using map_iterator = typename std::map<Key,Value>::iterator;

  map_key_iterator ( ) : map_iterator ( ) { };
  map_key_iterator ( map_iterator it_ ) : map_iterator ( it_ ) { };

  Key *operator -> ( ) { return ( Key * ) &( map_iterator::operator -> ( )->first ); }
  Key operator * ( ) { return map_iterator::operator * ( ).first; }
};

template<typename Key, typename Value>
class map_key_const_iterator<std::map<Key,Value>> final : public std::map<Key, Value>::const_iterator
{
public:
  using map_const_iterator = typename std::map<Key,Value>::const_iterator;

  map_key_const_iterator ( ) : map_const_iterator ( ) { };
  map_key_const_iterator ( map_const_iterator it_ ) : map_const_iterator ( it_ ) { };

  const Key *operator -> ( ) { return ( const Key * ) &( map_const_iterator::operator -> ( )->first ); }
  Key operator * ( ) { return map_const_iterator::operator * ( ).first; }
};

} // namespace scream

#endif // SCREAM_MAP_KEY_ITERATOR_HPP
