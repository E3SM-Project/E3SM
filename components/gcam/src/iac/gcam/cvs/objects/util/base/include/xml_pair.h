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

/* xml_pair.h
 * Created: 03/09/2007
 * Version: 04/05/2007
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

#if !defined( XML_PAIR_H )
#define XML_PAIR_H         // prevent multiple includes

// include files ***********************************************************

#include "util/base/include/xml_helper.h"
#include <xercesc/dom/DOMNode.hpp>
#include <string>

// namespaces **************************************************************

namespace ObjECTS {

// class: XMLPair **********************************************************

/*! \ingroup Objects
 *  \brief XMLPair<T> is a class template for a XML <key, value> pair
 *  \details The class template XMLPair<T> is used to the <key, value> pair
 *           of a named XML element that contains a value.
 * 
 *  \author Kevin Walker
 * \date $ Date $
 * \version $ Revision $
 */
template <class T = std::string>
class XMLPair
{
public :

   typedef T value_type;

   /*! Default constructor
    *  \param aKey the key (optional)
    *  \param aValue the value (optional)
    */
   XMLPair(
      const std::string& aKey = std::string(),
      const value_type&  aValue = value_type() )
      : mKey( aKey ), mValue( aValue ) {}
   /*! Copy constructor
    *  \param other the instance to copy
    */
   XMLPair( const XMLPair<T>& other )
      : mKey( other.mKey ), mValue( other.mValue ) {}
   /*! Copy from a STL pair
    *  \param aPair the pair to copy
    */
   XMLPair( const std::pair<std::string, value_type>& aPair )
      : mKey( aPair.first ), mValue( aPair.second ) {}

   /*! Assignment operator
    *  \param other the instance to copy
    *  \return *this
    */
   XMLPair<T>& operator = ( const XMLPair<T>& other )
   {
      if ( &other != this )
      {
         mKey   = other.mKey;
         mValue = other.mValue;
      }
      return *this;
   }
   /*! Assignment operator
    *  \param aPair the pair to copy
    *  \return *this
    */
   XMLPair<T>& operator = ( const std::pair<std::string, value_type>& aPair )
   {
      mKey   = aPair.first;
      mValue = aPair.second;
      return *this;
   }

   /*! Get the key
    *  \return the key
    */
   virtual const std::string& getKey( void ) const { return mKey; }

   /*! Get the value
    *  \return the value
    */
   virtual const value_type& getValue( void ) const { return mValue; }

    /*! Parse a pair from the specified node
    *  \param apNode the node to parse
    *  \return true on success, false otherwise
    */
   virtual bool parse( const xercesc::DOMNode * apNode );

   /*! Print to the specified output stream
    *  \param out the output stream
    *  \param apTabs the current tab setting
    *  \return the output stream
    */
   virtual std::ostream& print( std::ostream& out, Tabs * apTabs = 0 );

   /*! Set the key
    *  \param aKey the key to set
    */
   virtual void setKey( const std::string& aKey ) { mKey = aKey; }

   /*! Set the value
    *  \param aValue the value to set
    */
   virtual void setValue( const value_type& aValue ) { mValue = aValue; }

private :

   //! The key/name for the pair
   std::string mKey;
   //! The value for the pair
   value_type  mValue;
};

// XMLPair<T>::parse *******************************************************

 /*! Parse a pair from the specified node
 *  \param apNode the node to parse
 *  \return true on success, false otherwise
 */
template <class T>
inline bool XMLPair<T>::parse( const xercesc::DOMNode * apNode )
{
   // Validate the node
   if ( !apNode || apNode->getNodeType() != xercesc::DOMNode::ELEMENT_NODE )
   {
      return false;
   }

   // Get the node name (key)
   setKey( XMLHelper<std::string>::safeTranscode( apNode->getNodeName() ) );

   // Get the node value
   try
   {
      setValue( XMLHelper<T>::getValue( apNode ) );
      return true;
   }
   catch ( ... )
   {
   	return false;
   }
}

// XMLPair<T>::print *******************************************************

/*! Print to the specified output stream
 *  \param out the output stream
 *  \param apTabs the current tab setting
 *  \return the output stream
 */
template <class T>
inline std::ostream& XMLPair<T>::print(
   std::ostream& out,
   Tabs *        apTabs )
{
   XMLWriteElement( getValue(), getKey(), out, apTabs );
   return out;
}

} // namespace ObjECTS

#endif   // XML_PAIR_H

// end of xml_pair.h *******************************************************

