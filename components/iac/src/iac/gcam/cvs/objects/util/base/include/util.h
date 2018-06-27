#ifndef _UTIL_H_
#define _UTIL_H_
#if defined(_MSC_VER)
#pragma once
#endif

/*
 * LEGAL NOTICE
 * This computer software was prepared by Battelle Memorial Institute,
 * hereinafter the Contractor, under Contract No. DE-AC05-76RL0 1830
 * with the Department of Energy (DOE). NEITHER THE GOVERNMENT NOR THE
 * CONTRACTOR MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
 * LIABILITY FOR THE USE OF THIS SOFTWARE. This notice including this
 * sentence must appear on any copies of this computer software.
 * 
 * EXPORT CONTROL
 * User agrees that the Software will not be shipped, transferred or
 * exported into any country or used in any manner prohibited by the
 * United States Export Administration Act or any other applicable
 * export laws, restrictions or regulations (collectively the "Export Laws").
 * Export of the Software may require some form of license or other
 * authority from the U.S. Government, and failure to obtain such
 * export control license may result in criminal liability under
 * U.S. laws. In addition, if the Software is identified as export controlled
 * items under the Export Laws, User represents and warrants that User
 * is not a citizen, or otherwise located within, an embargoed nation
 * (including without limitation Iran, Syria, Sudan, Cuba, and North Korea)
 *     and that User is not otherwise prohibited
 * under the Export Laws from receiving the Software.
 * 
 * Copyright 2011 Battelle Memorial Institute.  All Rights Reserved.
 * Distributed as open-source under the terms of the Educational Community 
 * License version 2.0 (ECL 2.0). http://www.opensource.org/licenses/ecl2.php
 * 
 * For further details, see: http://www.globalchange.umd.edu/models/gcam/
 * 
 */


/*!
* \file util.h  
* \ingroup Objects
* \brief A set of commonly used functions.
* \details This is a set of functions which are frequently needed within the
*          program.
 * \note These are static functions within a namespace, not static class
 *       functions. This is because partial template specialization cannot be
 *       done for classes. The functions are in the objects namespace, also the
 *       location of the utility container classes. The util namespace is
 *       aliased to the objects, namespace, which means it is an alternative
 *       name for it.
* \author Josh Lurz
*/

#include "util/base/include/definitions.h"
#include "util/base/include/inamed.h"
#include <boost/static_assert.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <limits>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <sstream>
#include <map>
#include <vector>
#include <cassert>
#include <algorithm>

#ifndef _MSC_VER
#include "util/logger/include/ilogger.h"
#endif

// Boost static asserts do not work when included from multiple namespaces.
// Separate them into their own unique namespace.
namespace conditionsCheck {
    BOOST_STATIC_ASSERT( std::numeric_limits<double>::has_quiet_NaN );
    BOOST_STATIC_ASSERT( std::numeric_limits<double>::has_infinity );
}

namespace objects {
    
    /*! \brief Returns the value within this map associated with a given key. 
    * \details This function takes as its input a map and a key to search for.
    *          It will return the value associated with the key, or the default
    *          value for the class of the object if the key is not found.
    * \note Use this function instead of recoding a map search, as this function
    *       should be more efficient and handle errors more appropriately. 
    * \todo Evaluate returning a const reference.
    * \param currMap The map within which to search for the value.
    * \param key The key to find the value with which it is associated.
    * \return The value in the currMap associated with the key, the default
    *         constructed object otherwise. 
    */
    template <class K, class V>
    const V searchForValue( const std::map<K,V>& currMap, const K& key ){
        typedef typename std::map<K,V>::const_iterator CMapIterator;
        CMapIterator iter = currMap.find( key );
        if( iter != currMap.end() ){
            return iter->second;
        } else {
            return V(); // Returns default constructed value, 0 for doubles and
                        // ints
        }
    }

    /*! \brief Returns a constant iterator to a position in the vector which
    *          contains a pointer to an object with the given name. 
    * \details This function searches linearly from the beginning of the vector
    *          until it finds an object which has a getName() function which
    *          returns the given name. If the name is not found the end iterator
    *          will be returned.
    * \param aVector A vector to search.
    * \param aName A name for which to search.
    * \return An iterator to the position in the vector with the given name, the
    *         end iterator if that is not found.
    */
    template <class U>
    inline typename std::vector<U>::const_iterator searchForValue( const std::vector<U>& aVector, const std::string& aName ) {
        typename std::vector<U>::const_iterator iter = aVector.begin();
        for( ; iter!= aVector.end(); ++iter ){
            if( (*iter)->getName() == aName ){
                break;
            }
        }
        return iter;
    }


    /*! \brief Returns a mutable iterator to a position in the vector which
    *          contains a pointer to an object with the given name. 
    * \details This function searches linearly from the beginning of the vector
    *          until it finds an object which has a getName() function which
    *          returns the given name. If the name is not found the end iterator
    *          will be returned.
    * \param aVector A vector to search.
    * \param aName A name for which to search.
    * \return An iterator to the position in the vector with the given name, the
    *         end iterator if that is not found.
    */
    template <class U>
    inline typename std::vector<U>::iterator searchForValue( std::vector<U>& aVector, const std::string& aName ) {
        typename std::vector<U>::iterator iter = aVector.begin();
        for( ; iter!= aVector.end(); ++iter ){
            if( (*iter)->getName() == aName ){
                break;
            }
        }
        return iter;
    }

    template <class U, class T>
    struct InterfaceGetter {
    public:
        const U* operator()( const T& aObject ){
            return &aObject;
        }
    };

    template <class U,class T>
    struct InterfaceGetter<U,T*> {
    public:
        const U* operator()( const T* aObject ){
            return aObject;
        }
    };

    /*! 
     * \ingroup Objects
     * \brief Helper struct which implements the equals operator for INamed pointers.
     * \todo Is this worth a source file?
     * \author Josh Lurz
     */
    template <class T>
    struct NameEquals: public std::unary_function<const T,bool> {
        NameEquals( const std::string& aName );
        bool operator()( const T& aObject );

        //! Stored name to compare against.
        const std::string mName;
    };

    /*!
     * \brief Constructor for the operator which stores the name to compare against
     *        the name of another object.
     * \param aName Name of the object which is being compared against.
     */
    template <class T>
    inline NameEquals<T>::NameEquals( const std::string& aName ):mName( aName ){
        /*! \pre Object name is not empty. */
        assert( !mName.empty() );
    }

    /*! \brief Equals operator which determines if the given object has the same
     *         name as another object.
     * \param aObject Object to compare against.
     * \return Whether the object is equal to the stored object.
     */
    template <class T>
    inline bool NameEquals<T>::operator()( const T& aObject ){
        /*! \pre Object is non-null. */
        assert( aObject );
        /*! \pre Object name is not empty. */
        assert( !aObject->getName().empty() );

        return InterfaceGetter<INamed,T>()( aObject )->getName() == mName;
    }
    
    /*!
     * \brief A binary functor to compare named components by name so that they can be sorted.
     * \note If objects correctly inherited from the INamed interface this could
     *       just use that as the template argument instead of templating NameComparator.
     */
    template<class T>
    struct NameComparator : public std::binary_function<T*, T*, bool> {
        /*!
         * \brief Determine if the left hand operand is less than the
         *        right.
         * \param lhs The left hand side operand for the sort comparison.
         * \param rhs The right hand side operand for the sort comparison.
         * \return True if lhs is less than rhs, false otherwise.
         */
        bool operator()( const T* lhs, const T* rhs ){
            return lhs->getName() < rhs->getName();
        }
    };

    /*! \brief Reorder a container based on an ordering of names.
     *  \details This function is used to reorder a a container using a supplied
     *          ordering. This function handles errors as follows:
     *          <ul>
     *            <li>
     *              If an object is not specified in the list, the function
     *              will output an error.
     *            </li>
     *            <li>
     *              If an object name in the ordering list is not an
     *              existing sector in the model, a debugging warning will
     *              be issued and the name will be skipped.
     *            </li>
     *          </ul>
     * \param aFirst An InputIterator pointing to the start of the container.
     * \param aLast An InputIterator pointing to the end of the container.
     * \param aOrderList A list of names in the order in which the objects in
     *        the container should be put. 
     * \return Whether all objects had orderings in the passed in order list.
     * \warning This function requires that all the strings in aOrderList are
     *          unique.  It's behavior is undefined if the vector contains
     *          duplicate strings.
     */
    template <class InputIterator, class T>
    bool reorderContainer( InputIterator aFirst,
                           InputIterator aLast,
                           const std::vector<std::string>& aOrderList ){
        bool success = true;

        if( aOrderList.empty() ){
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::NOTICE );
            mainLog << "Skipping object reordering due to an empty order list."
                    << std::endl;
            return false;
        }

        // Get the dependency finder logger.
        ILogger& depFinderLog = ILogger::getLogger( "dependency_finder_log" );

        // Construct a map of strings to bools that represent an object in the
        // container NOT being mapped.  The bool is irrelevant because mappings
        // will be removed as they are ordered.
        std::map<std::string,bool> isOrdered;
        for( InputIterator containerIter = aFirst;
            containerIter != aLast; ++containerIter ){
            isOrdered.insert( 
                std::make_pair( objects::InterfaceGetter<INamed,T>()( *containerIter )->getName(), true ) );
        }

        // This loop functions much like insertion sort.  aFirst initially
        // points to the first location in the container.  The loop iterates
        // over the name list.  During each iteration we attempt to find
        // an iterator that points to object in the container that has the
        // right name. If this iterator is found the values are swapped.  Otherwise
        // a warning is emitted.
        typedef std::vector<std::string>::const_iterator OrderIter;
        for( OrderIter nameIter = aOrderList.begin(); 
             nameIter != aOrderList.end(); ++nameIter ){

            // Note that find_if uses a predicate to find the appropriate item.
            // Consult the NameEquals class and the STL documentation for details.
            InputIterator containerIter = std::find_if( aFirst,
                                                        aLast,
                                                        NameEquals<T>( *nameIter ) );
            if( containerIter != aLast ){
                std::swap( *aFirst, *containerIter );
                ++aFirst;
                // Remove it from the isOrdered vector
                isOrdered.erase( *nameIter );
            }
            else {
                // There was a name in the order that did not correspond to
                // and object in the container.
                depFinderLog.setLevel( ILogger::DEBUG );
                depFinderLog << *nameIter << " is not the name of an "
                             << "object in the container. "
                             << "It will not be included in the ordering."
                             << std::endl;
            }
        }
        // Output an error if there was an item in the container that was not
        // in the ordering.
        for( std::map<std::string,bool>::iterator orderIter = isOrdered.begin();
             orderIter != isOrdered.end(); ++orderIter ){
            ILogger& mainLog = ILogger::getLogger( "main_log");
            mainLog.setLevel( ILogger::ERROR );
            mainLog << orderIter->first
                    << " was not assigned a position in the explicit ordering list." << std::endl;
            success = false;
        }
        return success;
    }

    /*! \brief Returns whether a value with the given key exists.
    * \details This function takes as its input a map and a key for which to
    *          search, and returns whether the key exists. 
    * \param aCurrMap The map within which to search for the value.
    * \param aKey The key of which to check for the existance.
    * \return Whether the key exists. 
    */
    template <class K, class V> bool hasValue( const std::map<K,V>& aCurrMap, const K& aKey ){
        return ( aCurrMap.find( aKey ) != aCurrMap.end() );
    }

    /*! \brief A function to determine the sign of a number.
    * \param number A templated parameter which must be comparable to 0.
    * \return Returns -1 if the number is less than 0, +1 otherwise.
    */
    template <class T>
    int sign( const T number ) {
        return ( number < 0 )?(-1):(1);
    }

    /*! \brief Check the validity of a number.
    * \details Occasionally after calculations numbers are no longer in the
    *          range of real numbers, however C++ will continue to perform
    *          calculations on them. This function can then be used to check if
    *          numbers are not not a number or infinity.
    * \warning Some compiler/platforms do not support these checks. In that case
    *          this function will always return true. 
    * \warning Do not try to perform this check without using this function. It
    *          will fail on some compilers.
    * \param number A templated parameter which must be comparable to 0.
    * \return Returns whether the number is valid.
    */
    template <class T>
    inline typename boost::enable_if<boost::is_convertible<T, double>, bool>::type isValidNumber( const T aNumber ) {

        // Need to check whether the type supports not-a-number and infinity.
        const double doubleValue( aNumber );
        return ( aNumber == aNumber )  // check for NaN, by using the fact that == is always false for NaN
            && ( !std::numeric_limits<double>::has_infinity ||
            ( doubleValue != std::numeric_limits<double>::infinity()
            && std::negate<double>()( doubleValue ) != std::numeric_limits<double>::infinity() ) );
    }

    /*!
     * \brief Specialization of isValidNumber for type not convertible to double.
     * \details Since isValidNumber is used with in template datastructures we
     *          must be able to handle types which are not numbers
     * \param aType A non-number type to check.
     * \return true since this wouldn't be a valid check anyway.
     */
    template <class T>
    inline typename boost::disable_if<boost::is_convertible<T, double>, bool>::type isValidNumber( const T aNumber ){
        return true;
    }

    /*!
     * \brief This is a template function which compares two values. 
     * \details This function very simply uses the == operator of the two
     *          arguments to compare them, and returns the return value of the
     *          == operator. The reason for this function is so that it can be
     *          overridden for doubles to perform special comparison not using
     *          the == operator.
     * \param aFirstValue The first value to compare.
     *  \param aSecondValue The second value to compare.
     * \param aTolerance This parameter is unused and only for compatability
     *        with the double specialization of the function.
     * \return Whether or not the two values are equal.
     */
    template<class T>
    inline bool isEqual( const T aFirstValue,
                         const T aSecondValue,
                         const double aTolerance = 1E-10 )
    {
        return ( aFirstValue == aSecondValue );
    }

    /*!
     * \brief A function to determine if two doubles are equal.
     * \details Due to inaccuracies in machine arithmatic, it is virtually
     *          impossible for two doubles with decimal values that are
     *          calculated with different methods to be exactly equal. This
     *          function checks if the two values are within a very small
     *          threshhold. This an explicit template specialization for doubles
     *          which allows isEqual to act differently for doubles. These
     *          function had to be declared inline to avoid linker errors.
     * \warning Do not compare two doubles using the == operator. Use this
     *          function instead. 
     * \param aFirstValue The first double to compare.
     * \param aSecondValue The second double to compare.
     * \param aTolerance Tolerance to use when comparing the numbers. Defaults
     *        to 1E-10.
     * \return Whether the two doubles are within aTolerance of each other. 
     */
    template<>
    inline bool isEqual<double>( const double aFirstValue,
                                 const double aSecondValue,
                                 const double aTolerance )
    {
        return ( std::fabs( aFirstValue - aSecondValue ) < aTolerance );
    }

    /*
    * \brief Interpolate a Y value based on two points and an X value.
    * \param aX X value for which to find an X value.
    * \param aX1 First X value.
    * \param aX2 Second X value.
    * \param aY1 First Y value.
    * \param aY2 Second Y value.
    * \return Linearly interpolated Y value for aX.
    */
    double linearInterpolateY( const double aX,
                               const double aX1,
                               const double aX2,
                               const double aY1,
                               const double aY2 );

    /*! \brief A function to check if a file was opened successfully.
    *
    * When C++ opens a file, for input or output, it does not check whether it
    * was opened successfully. If the file has not been opened correctly, the
    * program will often not fail as FORTRAN programs do, but instead it will
    * have unexpected behavior. This function will check if the file is open,
    * and if it is not it will print an error message before calling abort. 
    * 
    * \todo This function should be more flexible, for use in cases where the
    *       missing file is a non-fatal problem.
    * \param streamIn A templated parameter that must have the function is_open.
    *        This was done so that one function could 
    * check both input files and output files.
    * \param fName The name of the file streamIn references, so that the error
    *        message can be more informative.
    */
    template <class T>
    inline void checkIsOpen( T& streamIn, const std::string& fName ) {
        if( !streamIn.is_open() ) {
            std::cerr << "Severe Error: File " << fName << " could not be opened." << std::endl;
            abort();
        }
    }
    
    std::string replaceSpaces( const std::string& aString );

    /*! \brief Static function which returns SMALL_NUM. 
    * \details This is a static function which is used to find the value of the
    *          constant SMALL_NUM. This avoids the initialization problems of
    *          static variables. This function should be used instead of
    *          defining this constant in multiple locations in the code.
    * \return The constant SMALL_NUM.
    */
   static inline double getSmallNumber() {
      const double SMALL_NUM = 1e-6;
      return SMALL_NUM;
   }

    /*! \brief Static function which returns VERY_SMALL_NUM. 
    * \details This is a static function which is used to find the value of the
    *          constant VERY_SMALL_NUM. This avoids the initialization problems
    *          of static variables. This function should be used instead of
    *          defining this constant in multiple locations in the code.
    * \return The constant VERY_SMALL_NUM.
    */
   static inline double getVerySmallNumber() {
      const double VERY_SMALL_NUM = 1e-8;
      return VERY_SMALL_NUM;
   }

    /*! \brief Static function which returns EXTREMELY_SMALL_NUM. 
    * \details This is a static function which is used to find the value of the
    *          constant EXTREMELY_SMALL_NUM. This avoids the initialization
    *          problems of static variables. This function should be used
    *          instead of defining this constant in multiple locations in the
    *          code.
    * \return The constant EXTREMELY_SMALL_NUM.
    */
   static inline double getTinyNumber() {
      const double EXTREMELY_SMALL_NUM = 1e-16;
      return EXTREMELY_SMALL_NUM;
   }

    /*! \brief Static function which returns LARGE_NUM. 
    * \details This is a static function which is used to find the value of the
    *          constant LARGE_NUM. This avoids the initialization problems of
    *          static variables. This function should be used instead of
    *          defining this constant in multiple locations in the code.
    * \return The constant LARGE_NUM.
    */
   static inline double getLargeNumber() {
      const double LARGE_NUM = 1e+6;
      return LARGE_NUM;
   }

    /*! \brief Function which returns a vector of keys from a map.
    * \details This function takes a map as an argument and returns a vector of
    *           all the keys of the map. It uses the same order as the map
    *           iterator returns.
    * \param aMap A map to return all keys for.
    * \return A vector of all keys from the map in the same order as the map
    *         iterator returns.
    */
    template<class T, class U> const std::vector<T> getKeys( const std::map<T,U> aMap ) {
        typedef typename std::map<T,U>::const_iterator ConstMapIterator;
        std::vector<T> keys;
        for( ConstMapIterator mapIter = aMap.begin(); mapIter != aMap.end(); mapIter++ ){
            keys.push_back( ( *mapIter ).first );
        }
        return keys;
    }
    
   /* \brief Function which returns a vector of values from a map.
   * \details This function takes a map as an argument and returns a vector of
   *          all the values of the map. It uses the same order as the map
   *          iterator returns.
   * \param aMap A map to return all values for.
   * \return A vector of all values from the map in the same order as the map
   *         iterator returns.
   */
    template<class T, class U>
        const std::vector<U> getValues( const std::map<T,U> aMap )
    {
        typedef typename std::map<T,U>::const_iterator ConstMapIterator;
        std::vector<U> values;
        for( ConstMapIterator mapIter = aMap.begin(); mapIter != aMap.end(); mapIter++ ){
            values.push_back( ( *mapIter ).second );
        }
        return values;
    }

    /*!
     * \brief Round a double to the nearest whole number.
     * \param aValue Value to round.
     * \return Value rounded to the nearest whole number.
     */
    static int round( const double aValue ){
        double intPart;
        if( modf( aValue, &intPart ) > 0.5 ){
            ++intPart;
        }
        return static_cast<int>( intPart );
    }

    /*!
     * \brief Calculate the percent difference between an initial and final
     *        value.
     * \details Calculates the percent difference between the final and initial
     *          values as a percent, using the initial value as the denominator.
     *          If the initial value is zero the absolute difference will be
     *          used. This percentage is not an absolute, and will be signed.
     *          The difference is calculated as final minus initial.
     * \param aInitialValue The initial value. This will be used as the
     *        denominator.
     * \param aFinalValue The final value.
     * \return The percent difference between the initial and final terms.
     */
    static double percentDiff( const double aInitialValue,
                               const double aFinalValue )
    {
        double difference = aFinalValue - aInitialValue;
        if( aInitialValue > getSmallNumber() ){
            difference /= aInitialValue;
            difference *= 100;
        }
        return difference;
    }

    /*! \brief Convert a value to a string using the built in stringstream.
    * \todo Remove this function and use the boost equivalent.
    */
    template<class T> std::string toString( const T& value ){
        static std::stringstream converter;
        converter << value;
        std::string output;
        converter >> output;
        converter.clear();
        return output;
    }

   long createMinicamRunID( const time_t& aTime );
   std::string XMLCreateDate( const time_t& time );

   tm* getGMTime( const time_t& aTime );
   tm* getLocalTime( const time_t& aTime );
   void printTime( const time_t& aTime, std::ostream& aOut );
} // End util namespace.

namespace util = objects;

#endif // _UTIL_H_
