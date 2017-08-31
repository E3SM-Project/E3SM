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

/* TValidatorInfo.h
 * Created: 02/02/2007
 * Version: 02/21/2007
 *
 * This software, which is provided in confidence, was prepared by employees
 * of Pacific Northwest National Laboratory operated by Battelle Memorial
 * Institute. Battelle has certain unperfected rights in the software
 * which should not be copied or otherwise disseminated outside your
 * organization without the express written authorization from Battelle.
 * All rights to the software are reserved by Battelle.   Battelle makes no
 * warranty, express or implied, and assumes no liability or responsibility
 * for the use of this software.
 */

#if !defined( __TVALIDATORINFO_H )
#define __TVALIDATORINFO_H    // prevent multiple includes

// include files ***********************************************************

#include <string>

// namespaces **************************************************************

namespace ObjECTS {

// class: TValidatorInfo ***************************************************

/*! \ingroup Objects
 *  \brief TValidatorInfo<T> is a class template that can be used to
 *         manage parameter validation information
 *  \details The class template TValidatorInfo<T> is used to store
 *           the status of a named attribute for validation purposes.
 *           T is the type of the attribute, and the default type is double.
 *  \note TValidatorInfo<T> itself does not perform any validation.
 *        The developer is responsible for setting the validation
 *        in the contructor or the TValidatorInfo<T>::setIsValid() method.
 * 
 *  \author Kevin Walker
 */
template <class T = double>
class TValidatorInfo
{
public :

   typedef std::string string_type;
   typedef T           value_type;

   //! Default constructor
   TValidatorInfo(
      const value_type&  aValue   = value_type(),
      const string_type& aName    = string_type(),
      bool               aIsValid = true )
      : mValue( aValue ), mName( aName ), mIsValid( aIsValid ) {}

   //! Copy constructor
   TValidatorInfo( const TValidatorInfo& other )
      : mValue( other.mValue ),
        mName( other.mName ),
        mIsValid( other.mIsValid ) {}

   //! Destructor
   virtual ~TValidatorInfo( void ) {}

   //! Assignment operator
   TValidatorInfo& operator = ( const TValidatorInfo& other )
   {
      if ( &other != this )
      {
         mValue   = other.mValue;
         mName    = other.mName;
         mIsValid = other.mIsValid;
      }
      return *this;
   }

   //! Get if the attribute is valid
   virtual bool getIsValid( void ) const { return mIsValid; }

   //! Get the name of the attribute to validate
   virtual const string_type& getName( void ) const { return mName; }

   //! Get the value of the attribute to validate
   virtual const value_type& getValue( void ) const { return mValue; }

   //! Set if the attribute is valid
   //! \param aIsValid true for valid, false otherwise
   virtual void setIsValid( bool aIsValid ) { mIsValid = aIsValid; }

   //! Set the name of the attribute to validate
   //! \param aName the attribute name
   virtual void setName( const string_type& aName ) { mName = aName; }

   //! Set the value of the attribute to validate
   //! \param aValue the attribute value
   virtual void setValue( const value_type& aValue ) { mValue = aValue; }

private :

   //! The attribute value
   value_type   mValue;

   //! The attribute name
   string_type mName;

   //! If the attribute is valid
   bool        mIsValid;
};

// getInvalidNames *********************************************************

/*! Count and result the number of invalid validators in the specified
 *  collection
 *  \param first - input iterator to the beginning of the collection
 *  \param last - input iterator to the one past the end of the collection
 */
template <class InputIterator>
inline std::string getInvalidNames(
   InputIterator first,
   InputIterator last )
{
   std::string result;
   for ( InputIterator x = first; x != last; ++x )
   {
      if ( !( *x ).getIsValid() )
      {
         if ( result.length() )
         {
            result += ", ";
         }
         result += ( *x ).getName();
      }
   }

   return result;
}

// testValidators **********************************************************

/*! Count and result the number of invalid validators in the specified
 *  collection
 *  \param first - input iterator to the beginning of the collection
 *  \param last - input iterator to the one past the end of the collection
 */
template <class InputIterator>
inline size_t testValidators(
   InputIterator first,
   InputIterator last )
{
   size_t   result = 0;
   for ( InputIterator x = first; x != last; ++x )
   {
      result += ( *x ).getIsValid() ? 0 : 1;
   }

   return result;
}

} // namespace ObjECTS

#endif   // __TVALIDATORINFO_H

// end of TValidatorInfo.h *************************************************

