#ifndef _IINFO_H_
#define _IINFO_H_
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
* \file iinfo.h
* \ingroup objects
* \brief The IInfo interface header file.
* \author Josh Lurz
*/
#include <string>
#include <iosfwd>

class Tabs;

/*!
* \ingroup Objects
* \brief This interface represents a set of properties which can be accessed by
*        their unique identifier.
* \details The IInfo interface represents a set of searchable properties
*          accessed by their string key. The properties may be booleans,
*          integers, double or strings. Operations exist to set or update values
*          for a key, query if a key exists, and get the value for a key.
* \todo Evaluate whether functions to add to a double value, and update an
*       average would be useful as additions to the interface.
* \todo Add longevity to properties.
* \author Josh Lurz
*/

class IInfo
{
public:
    virtual inline ~IInfo();

    /*! \brief Set a boolean value for a given key.
    * \details Updates the value associated with the key if it is already
    *          present, and creates a new one key-value pair if it does not
    *          exist.
    * \param aStringKey The key for which to set or update the value.
    * \param aValue The new value.
    */
    virtual bool setBoolean( const std::string& aStringKey,
                             const bool aValue ) = 0;

    /*! \brief Set an integer value for a given key.
    * \details Updates the value associated with the key if it is already
    *          present, and creates a new key-value pair if it does not exist.
    * \param aStringKey The key for which to set or update the value.
    * \param aValue The new value.
    */
    virtual bool setInteger( const std::string& aStringKey,
                             const int aValue ) = 0;
    
    /*! \brief Set a double value for a given key.
    * \details Updates the value associated with the key if it is already
    *          present, and creates a new key-value pair if it does not exist.
    * \param aStringKey The key for which to set or update the value.
    * \param aValue The new value.
    */
    virtual bool setDouble( const std::string& aStringKey,
                            const double aValue ) = 0;
    
    /*! \brief Set a string value for a given key.
    * \details Updates the value associated with the key if it is already
    *          present, and creates a new key-value pair if it does not exist.
    * \param aStringKey The key for which to set or update the value.
    * \param aValue The new value.
    */
    virtual bool setString( const std::string& aStringKey,
                            const std::string& aValue ) = 0;

    /*! \brief Get a boolean from the IInfo with a specified key.
    * \details Searches the key set for the given key and returns the associated
    *          value. If the value does not exist the default value will be
    *          returned.
    * \param aStringKey The key for which to search the IInfo object.
    * \param aMustExist Whether the value should exist in the IInfo.
    * \return The boolean associated with the key or false if it does not exist.
    */
    virtual bool getBoolean( const std::string& aStringKey,
                                   const bool aMustExist ) const = 0;
    
    /*! \brief Get a integer from the IInfo with a specified key.
    * \details Searches the key set for the given key and returns the associated
    *          value. If the value does not exist the default value will be
    *          returned.
    * \param aStringKey The key for which to search the IInfo object.
    * \param aMustExist Whether the value should exist in the IInfo.
    * \return The integer associated with the key or zero if it does not exist.
    */
    virtual int getInteger( const std::string& aStringKey,
                                  const bool aMustExist ) const = 0;
    
    /*! \brief Get a double value from the IInfo with a specified key.
    * \details Searches the key set for the given key and returns the associated
    *          value. If the value does not exist the default value will be
    *          returned.
    * \param aStringKey The key for which to search the IInfo object.
    * \param aMustExist Whether the value should exist in the IInfo.
    * \return The double associated with the key or zero if it does not exist.
    */
    virtual double getDouble( const std::string& aStringKey,
                                    const bool aMustExist ) const = 0;
    
    /*! \brief Get a string from the IInfo with a specified key.
    * \details Searches the key set for the given key and returns the associated
    *          value. If the value does not exist the default value will be
    *          returned.
    * \param aStringKey The key for which to search the IInfo object.
    * \param aMustExist Whether the value should exist in the IInfo.
    * \return The string(by reference) associated with the key or the empty
    *         string if it does not exist.
    */
    const virtual std::string& getString( const std::string& aStringKey,
                                          const bool aMustExist ) const = 0;

    /*! \brief Get a boolean from the IInfo with a specified key.
    * \details Searches the key set for the given key and returns the associated
    *          value. If the value does not exist the default value will be
    *          returned.
    * \param aStringKey The key for which to search the IInfo object.
    * \param aFound Whether the value is found or not.
    * \return The boolean associated with the key or false if it does not exist.
    */
    virtual bool getBooleanHelper( const std::string& aStringKey, bool& aFound ) const = 0;

    /*! \brief Get a integer from the IInfo with a specified key.
    * \details Searches the key set for the given key and returns the associated
    *          value. If the value does not exist the default value will be
    *          returned.
    * \param aStringKey The key for which to search the IInfo object.
    * \param aFound Whether the value is found or not.
    * \return The integer associated with the key or zero if it does not exist.
    */
    virtual int getIntegerHelper( const std::string& aStringKey, bool& aFound ) const = 0;

    /*! \brief Get a double value from the IInfo with a specified key.
    * \details Searches the key set for the given key and returns the associated
    *          value. If the value does not exist the default value will be
    *          returned.
    * \param aStringKey The key for which to search the IInfo object.
    * \param aFound Whether the value is found or not.
    * \return The double associated with the key or zero if it does not exist.
    */
    virtual double getDoubleHelper( const std::string& aStringKey, bool& aFound ) const = 0;

    /*! \brief Get a string from the IInfo with a specified key.
    * \details Searches the key set for the given key and returns the associated
    *          value. If the value does not exist the default value will be
    *          returned.
    * \param aStringKey The key for which to search the IInfo object.
    * \param aFound Whether the value is found or not.
    * \return The string(by reference) associated with the key or the empty
    *         string if it does not exist.
    */
    virtual const std::string& getStringHelper( const std::string& aStringKey, bool& aFound ) const = 0;

    /*! \brief Return whether a value exists in the IInfo.
    * \details Performs a search of the IInfo object using the same method as
    *          all getter methods of the object. The method performs a full
    *          search and so should be avoided whenever the default value and a
    *          warning from the get method will return enough information.
    * \param aStringKey The key for which to search the IInfo object.
    * \return Whether the key exists in the IInfo.
    */
    virtual bool hasValue( const std::string& aStringKey ) const = 0;

    /*! \brief Write the IInfo object to an output stream as XML.
    * \details Writes the set of keys and values to an output stream as XML.
    * \param aPeriod Model period for which to write debugging information.
    * \param aTabs Tabs manager.
    * \param aOut Output stream.
    */
    virtual void toDebugXML( const int aPeriod,
                             Tabs* aTabs,
                             std::ostream& aOut ) const = 0;
};

//! Empty inline destructor needed so that IInfo objects can be deleted through
//! base class pointers.
IInfo::~IInfo(){
}

#endif // _IINFO_
