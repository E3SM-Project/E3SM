#ifndef _INFO_H_
#define _INFO_H_
#if defined(_MSC_VER_)
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
* \file info.h
* \ingroup objects
* \brief The Info class header file.
* \author Josh Lurz
*/

#include <string>
#include <iosfwd>
#include <boost/any.hpp>
#include <boost/noncopyable.hpp>
#include "containers/include/iinfo.h"

// Can't forward declare because operations are used in template functions.
#include "util/base/include/hash_map.h"
#include "util/base/include/xml_helper.h"
#include "util/base/include/configuration.h"

class Tabs;

/*!
* \ingroup Objects
* \brief This class contains a set of properties which can be accessed by their
*        unique identifier.
* \author Josh Lurz
* \todo Add longevity to properties.
*/

class Info: public IInfo, boost::noncopyable
{
    friend class InfoFactory;
public:
    ~Info();

    bool setBoolean( const std::string& aStringKey, const bool aValue );

    bool setInteger( const std::string& aStringKey, const int aValue );

    bool setDouble( const std::string& aStringKey, const double aValue );

    bool setString( const std::string& aStringKey, const std::string& aValue );

    bool getBoolean( const std::string& aStringKey, const bool aMustExist ) const;

    int getInteger( const std::string& aStringKey, const bool aMustExist ) const;

    double getDouble( const std::string& aStringKey, const bool aMustExist ) const;

    const std::string& getString( const std::string& aStringKey, const bool aMustExist ) const;

    bool getBooleanHelper( const std::string& aStringKey, bool& aFound ) const;

    int getIntegerHelper( const std::string& aStringKey, bool& aFound ) const;

    double getDoubleHelper( const std::string& aStringKey, bool& aFound ) const;

    const std::string& getStringHelper( const std::string& aStringKey, bool& aFound ) const;

    bool hasValue( const std::string& aStringKey ) const;

    void toDebugXML( const int aPeriod, Tabs* aTabs, std::ostream& aOut ) const;
protected:
    Info( const IInfo* aParentInfo, const std::string& aOwnerName );

private:

    std::string mOwnerName;

    /*!
     * \brief Enum representing possible types of each item.
     */
    enum AnyType {
        //! Boolean
        eBoolean,

        //! Integer
        eInteger,

        //! Double
        eDouble,

        //! String
        eString
    };

    template<class T> bool setItemValueLocal( const std::string& aStringKey,
                                              const AnyType aType,
                                              const T& aValue );

    template<class T> const T& getItemValueLocal( const std::string& aStringKey, bool& aExists ) const;

    size_t getInitialSize() const;

    void printItemNotFoundWarning( const std::string& aStringKey ) const;

    void printBadCastWarning( const std::string& aStringKey, bool aIsUpdate ) const;

    void printShadowWarning( const std::string& aStringKey ) const;

    template<class T> void printItem( const boost::any& aValue,
                                      std::ostream& aOut,
                                      Tabs* aTabs ) const;

    //! Type of the item stored as the value in the storage map.
    typedef std::pair<AnyType, boost::any> ValueType;

    //! Type of the internal storage map.
    typedef HashMap<const std::string, ValueType> InfoMap;

    //! Internal storage mapping item names to item values.
    std::auto_ptr<InfoMap> mInfoMap;

    //! A pointer to the parent of this Info object which can be null.
    const IInfo* mParentInfo;
};

/*! \brief Set a name and value for a piece of information related.
* \details This function will check the info map for the associated key, if it
*          exists it will update the associated value to itemValue, otherwise it
*          will create a new name value pair.
* \param aStringKey The string key to use as the key for this information value.
* \param aType Enum value of the type.
* \param aValue The value to be associated with this key. 
*/
template<class T> bool Info::setItemValueLocal( const std::string& aStringKey,
                                                const AnyType aType,
                                                const T& aValue )
{
    /*! \pre A valid key was passed. */
    assert( !aStringKey.empty() );

    // If debug checking is turned on search for the key in the current object
    // to determine if the type of the existing and new types match. Search in
    // the parent to see if this new item will shadow a variable in the parent.
    const static bool debugChecking = Configuration::getInstance()->getBool( "debugChecking" );
    if( debugChecking ){
        InfoMap::const_iterator curr = mInfoMap->find( aStringKey );
        if( curr != mInfoMap->end() ){
            // Check that the types match.
            try {
                boost::any_cast<T>( curr->second.second );
            }
            catch( boost::bad_any_cast ){
                printBadCastWarning( aStringKey, true );
            }
        }
        // Check if the value exists in a parent info.
        else if( mParentInfo && mParentInfo->hasValue( aStringKey ) ){
            printShadowWarning( aStringKey );
        }
    }

    // Add the value regardless of whether a warning was printed.
    mInfoMap->insert( std::make_pair( aStringKey, std::make_pair( aType, boost::any( aValue ) ) ) );
    return true;
}

/*! \brief Get the value of an item based on a string key and a passed in
*          location.
* \author Josh Lurz
* \param aStringKey The string key for which to find the value.
* \param aMustExist Whether it is an error for the item to be missing.
* \param aExists Return parameter to update with whether the item existed.
* \return The value associated with the key if it exists, the default value
*         otherwise.
*/
template<class T>
const T& Info::getItemValueLocal( const std::string& aStringKey,
                                 bool& aExists ) const
{
    /*! \pre A valid key was passed. */
    assert( !aStringKey.empty() );

    // Check for the value.
    InfoMap::const_iterator curr = mInfoMap->find( aStringKey );
    if( curr != mInfoMap->end() ){
        aExists = true;
        // Attempt to set the return value to the found value. This requires
        // converting the data from the actual type to the requested type.
        try {
            return *boost::any_cast<T>( &curr->second.second );
        }
        // Catch bad data conversions and print an error.
        catch( boost::bad_any_cast ){
            printBadCastWarning( aStringKey, false );
        }
    }

    // Return the default value if a successful return has not already occurred.
    aExists = false;
    static const T defaultValue = T();
    return defaultValue;
}

/*!
 * \brief Print a single any type value to XML.
 * \param aValue Value stored as an any type.
 * \param aOut Output stream.
 * \param aTabs Tabs object.
 */
template<class T>
void Info::printItem( const boost::any& aValue,
                      std::ostream& aOut,
                      Tabs* aTabs ) const
{
    try {
        T value = boost::any_cast<T>( aValue );
        XMLWriteElement( value, "Value", aOut, aTabs );
    }
    // Catch bad data conversions.
    catch( boost::bad_any_cast ){
        // Wrong type stored with the item.
        assert( false );
    }
}

#endif // _INFO_H_
