#ifndef _XML_HELPER_H_
#define _XML_HELPER_H_
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
 * \file xml_helper.h
 * \ingroup Objects
 * \brief A set of helper function for reading and writing xml data.
 * \note This file contains two things.
 *       - XMLHelper, A static class that has methods for parsing XML data and
 *         static data members which cache information required by the parser.
 *       - A series of global utility functions for writing XML data.
 * \todo XMLHelper should be converted into a non-static XMLReader class. The
 *       static data members could then be regular data members. The interface
 *       to the class should not use template functions, but the class could use
 *       them as helper methods. There are several functions that are used
 *       to read XML that are not part of XMLHelper. These should be moved in.
 * \todo This file needs refactoring. The XML writing utility functions should
 *       be moved to a non-static XMLWriter class. The class should store the
 *       tabs object and output stream.
 * \warning This class is hacked b/c of poor MSVC template support. This makes
 *          it much uglier.
 * \details This library contains a set of routines for reading xml data and
 *          attribute values. It is a templated library so that it should work
 *          with any data type.
 * \author Josh Lurz
 */

#include "util/base/include/definitions.h"
#include <string>
#include <iostream>
#include <sstream>
#include <cassert>
#include <vector>
#include <map>
#include <memory>

#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMAttr.hpp>
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMException.hpp>
#include <xercesc/dom/DOMNamedNodeMap.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include <xercesc/sax/HandlerBase.hpp>
#include <xercesc/util/PlatformUtils.hpp>

#if( __USE_XML_DB__ )
#include <xqilla/utils/XQillaPlatformUtils.hpp>
#endif

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/trim.hpp>

#include "util/base/include/model_time.h"
#include "util/base/include/util.h"
#include "util/logger/include/ilogger.h"
#include "util/base/include/iparsable.h"
#include "util/base/include/time_vector.h"
#include "util/base/include/value.h"

/*!
 * \ingroup Objects
 * \brief A basic class which is a container for a variable containing the
 *        current level of indention in the xml being written.
 * \todo Replace this class with an integer.
 * \author Josh Lurz
 */

class Tabs {

private:
    //! Current number of tabs to write out in order to properly format xml.
   unsigned short mNumTabs;
   //! Number of spaces making up a tab
   unsigned short mTabWidth;
   //! Use tabs or spaces
   bool mUseTabs;
public:

   enum { DEFAULT_TAB_WIDTH = 3 };

   //! Constructor
   Tabs(
      bool           aUseTabs  = true,
      unsigned short aTabWidth = DEFAULT_TAB_WIDTH )
      : mNumTabs(0), mTabWidth(aTabWidth), mUseTabs(aUseTabs) {}

   //! Increase the current level of indentation.
   void increaseIndent() { ++mNumTabs; }

   //! Decrease the current level of indentation.
   void decreaseIndent() {
      if ( mNumTabs ) {
        --mNumTabs;
      }
   }

   /*! Write out the contained number of tabs to the specified output stream.
    *
    * \param out Stream to which to write the tabs->
    * \return void
    */
   void writeTabs( std::ostream& out ) const {
      if ( mUseTabs ) {
         for ( int i = 0; i != mNumTabs; ++i ) {
            out << "\t";
         }
      }
      else if ( mTabWidth ) {
         std::string::size_type  n = mTabWidth * mNumTabs;
         if ( n ) {
            std::string spaces( n, ' ' );
            out << spaces;
         }
      }
   }
};

/*!
 * \ingroup Objects
 * \brief A class with static functions to parse XML DOM trees.
 * \author Josh Lurz
 */

template<class T>
class XMLHelper {
public:
   static T getValue( const xercesc::DOMNode* node );
   static T getAttr( const xercesc::DOMNode* node, const std::string attrName );
   static std::string safeTranscode( const XMLCh* toTranscode );

   static void insertValueIntoVector( const xercesc::DOMNode* node,
                                      std::vector<T>& insertToVector,
                                      const Modeltime* modeltime );

   static void insertValueIntoVector( const xercesc::DOMNode* aNode,
                                      objects::YearVector<T>& aYearVector );

   static void insertValueIntoVector( const xercesc::DOMNode* aNode,
                                      objects::PeriodVector<T>& aPeriodVector,
                                      const Modeltime* aModeltime );

   static int getNodePeriod ( const xercesc::DOMNode* node, const Modeltime* modeltime );
   static bool parseXML( const std::string& aXMLFile, IParsable* aModelElement );
   static const std::string& text();
   static const std::string& name();
   static void cleanupParser();
   static void printXMLTrace( const xercesc::DOMNode* aNode, std::ostream& aOut );
   static void serializeNode( const xercesc::DOMNode* aNode, std::ostream& aOut, Tabs* aTabs,
                              const bool aDeep );
private:
    static xercesc::XercesDOMParser** getParserPointerInternal();
    static xercesc::ErrorHandler** getErrorHandlerPointerInternal();
    static void initParser();
    static xercesc::XercesDOMParser* getParser();
};


/*! \brief Returns the data value associated with the element node.
* \details This function first finds the child node of this element, a text node which contains the data value.
* It then converts the data using the stringstream library into the requested data type and returns the value.
* The function will throw an error if the input node is not an element, but will return the default value if the element has no data value.
* \warning Since this function has a templated type only for the return argument, it must be called as getXMLValue<type>.
* \param node A pointer to a node for which to return the value.
* \return Value of type T from child of node.
* \sa getAttr
*/
template<class T>
T XMLHelper<T>::getValue( const xercesc::DOMNode* node ){
   // make sure we were passed a valid node reference which is an element.
   assert( node );
   assert( node->getNodeType() == xercesc::DOMNode::ELEMENT_NODE
           || node->getNodeType() == xercesc::DOMNode::ATTRIBUTE_NODE );

   // get the first child, which should contain the value.
   xercesc::DOMNode* curr = node->getFirstChild();

   // make sure that the above returned a TEXT_NODE, otherwise value will not be correct.
   if ( !curr || curr->getNodeType() != xercesc::DOMNode::TEXT_NODE ){
      return T();
   }

   // Attempt to convert the XML string to the appropriate type. This uses the
   // boost library lexical_cast operation, which is similar to the C++
   // static_cast, but allows conversion from a string into its numerical value.
   // This operation will fail if the string cannot be converted into the
   // expected type, in which case this function will return the default value
   // of the type.
   try {
       T returnValue = boost::lexical_cast<T>( safeTranscode( curr->getNodeValue() ) );
       return returnValue;
   }
   catch( boost::bad_lexical_cast& ) {
       try {
           // Cast the value to a string to print a more useful error message.
           // This cast should not fail because the value is read as a string.
           const std::string valueAsString = boost::lexical_cast<std::string>( safeTranscode( curr->getNodeValue() ) );
           std::cout << "Cast of node with value " << valueAsString << " to return type failed." << std::endl;
       }
       catch( boost::bad_lexical_cast& ){
           // The cast to a string should never fail because a string is read
           // in.
           assert( false );
       }
       // Set the return value to the default value.
       return T();
   }
}

/*! Returns the requested attribute of the element node passed to the function.
* \details This function searches for the attribute with name attrName of the argument node.
* It then converts it to type T and returns the value. If the function is not passed an element
* it will throw an error. If the requested attribute is not present, the function will return the default
* constructor for type T. (zero for doubles, or false for boolean)
* \warning It must be called as getXMLValue<type> because it is templated only on the return type.
* \param node A pointer to a node for which to fetch the attribute.
* \param attrName The name of the attribute to fetch.
* \return Value of type T from the attribute with name attrName of the node.
* \sa getValue
*/
template<class T>
T XMLHelper<T>::getAttr( const xercesc::DOMNode* node, const std::string attrName ) {
   /*! \pre Make sure we were passed a valid node reference. */
   assert( node );

   /*! \pre Make sure it is an element before we cast, if function is used correctly it will be. */
   assert( node->getNodeType() == xercesc::DOMNode::ELEMENT_NODE );

   // need to cast node to an element.
   const xercesc::DOMElement* element = static_cast<const xercesc::DOMElement*>( node );

   // get the attribute with the name which was passed in.

   XMLCh* nameChars = xercesc::XMLString::transcode( attrName.c_str() );
   xercesc::DOMAttr* nameAttr = element->getAttributeNode( nameChars );
   xercesc::XMLString::release( &nameChars );

   // Check if the name attribute exists.
   if( !nameAttr ){
      return T();
   }

   // Attempt to convert the XML string to the appropriate type. This uses the
   // boost library lexical_cast operation, which is similar to the C++
   // static_cast, but allows conversion from a string into its numerical value.
   // This operation will fail if the string cannot be converted into the
   // expected type, in which case this function will return the default value
   // of the type.
   try {
       T returnValue = boost::lexical_cast<T>( safeTranscode( nameAttr->getValue() ) );
       return returnValue;
   }
   catch( boost::bad_lexical_cast& ) {
       try {
           // Cast the value to a string to print a more useful error message.
           // This cast should not fail because the value is read as a string.
           const std::string valueAsString = boost::lexical_cast<std::string>( safeTranscode( nameAttr->getValue() ) );
           std::cout << "Cast of node with value " << valueAsString << " to return type failed." << std::endl;
       }
       catch( boost::bad_lexical_cast& ){
           // The cast to a string should never fail because a string is read
           // in.
           assert( false );
       }
       // Set the return value to the default value.
       return T();
   }
}

/*!
* \brief Function which takes a node and inserts its value into the correct position
* in an argument vector based on the year attribute.
* Updated to fill out for the rest of the time period if the fillout attribute is true.
*
* This function when passed a node, vector and modeltime object will first extract the year attribute and lookup
* the corresponding period from the modeltime object. It will then insert the item in that position in the vector.
*
* \warning Make sure the node passed as an argument as a year attribute.
* \param node A pointer to a node from which to extract the data.
* \param insertToVector A vector passed by reference in which to insert the value.
* \param modeltime A pointer to the modeltime object to use to determine the correct period.
*/

template<class T>
void XMLHelper<T>::insertValueIntoVector( const xercesc::DOMNode* node, std::vector<T>& insertToVector, const Modeltime* modeltime ) {

   /*! \pre Make sure we were passed a valid node reference. */
   assert( node );

   const int year = XMLHelper<int>::getAttr( node, "year" );
   // boolean to fill out the readin value to all the periods
   const bool fillout = XMLHelper<bool>::getAttr( node, "fillout" );

   // Check to make sure the year attribute returned non-zero.
   if (  year == 0 ) {
       ILogger& mainLog = ILogger::getLogger( "main_log" );
       mainLog.setLevel( ILogger::ERROR );
       mainLog << "Year value not set for vector input tag >>" << XMLHelper<std::string>::safeTranscode( node->getNodeName() ) << "<<" << std::endl;
       return;
   }

   const int maxperiod = modeltime->getmaxper();
   int period = modeltime->getyr_to_per( year );

   // Check that the period returned correctly.
   // Check to make sure the year attribute returned non-zero.
   if ( !( ( period >= 0 ) && ( period < modeltime->getmaxper() ) ) ) {
       ILogger& mainLog = ILogger::getLogger( "main_log" );
       mainLog.setLevel( ILogger::ERROR );
       mainLog << "Period value is out of bounds for vector input tag >>" << XMLHelper<std::string>::safeTranscode( node->getNodeName() ) << "<<" << std::endl;
       return;
   }

   // Check that the period is less than the size of the vector.
   assert( period < static_cast<int>( insertToVector.size() ) );
   insertToVector[ period ] =  XMLHelper<T>::getValue( node );

   if (fillout) {
      // Check that the max period is equal to the size of the vector.
      assert( maxperiod == static_cast<int>( insertToVector.size()) );
      // will not do if period is already last period or maxperiod
      for ( int i = period + 1; i < maxperiod; ++i ) {
         insertToVector[ i ] =  insertToVector[ period ];
      }
   }
}

/*!
* \brief Function which takes a node and inserts its value into the correct
*        position in a YearVector based on the year XML attribute. If the XML
*        attribute "fillout" is set the data will be copied to the end of the
*        vector.
*
* This function when passed a node and YearVector will first extract the year
* attribute and will then insert the item in that position in the vector.
*
* \warning Make sure the node passed as an argument has a year attribute.
* \param aNode A pointer to a node from which to extract the data.
* \param aYearVector A YearVector passed by reference in which to insert the
*        value.
*/
template<class T>
void XMLHelper<T>::insertValueIntoVector( const xercesc::DOMNode* aNode,
                                          objects::YearVector<T>& aYearVector )
{
   /*! \pre Make sure we were passed a valid node reference. */
   assert( aNode );

   const int year = XMLHelper<int>::getAttr( aNode, "year" );

   // boolean to fill out the readin value to all the periods
   const bool fillout = XMLHelper<bool>::getAttr( aNode, "fillout" );

   // Check to make sure the year attribute returned non-zero.
   if( year == 0 ) {
       ILogger& mainLog = ILogger::getLogger( "main_log" );
       mainLog.setLevel( ILogger::ERROR );
       mainLog << "Year value not set for vector input tag: "
           << XMLHelper<std::string>::safeTranscode( aNode->getNodeName() )
           << "." << std::endl;
       return;
   }
   // Ensure that the year is legal.
   typename objects::YearVector<T>::iterator pos = aYearVector.find( year );
   if( pos == aYearVector.end() ){
       ILogger& mainLog = ILogger::getLogger( "main_log" );
       mainLog.setLevel( ILogger::ERROR );
       mainLog << "Year value is out of bounds for input tag:"
           << XMLHelper<std::string>::safeTranscode( aNode->getNodeName() )
           << "." << std::endl;
       return;
   }

   // Get the value for the node.
   double value = XMLHelper<T>::getValue( aNode );

   // Assign the value into the vector.
   *pos = value;

   if( fillout ) {
       // Set the value for all years after the current point.
       ++pos;
       for( ; pos != aYearVector.end(); ++pos ){
           *pos = value;
       }
   }
}

/*!
* \brief Function which takes a node and inserts its value into the correct
*        position in a PeriodVector based on the year XML attribute. If the XML
*        attribute "fillout" is set the data will be copied to the end of the
*        vector.
* \details This function when passed a node and PeriodVector will first extract
*          the year attribute and will then insert the item in that position in
*          the vector.
* \warning Make sure the node passed as an argument has a year attribute.
* \param aNode A pointer to a node from which to extract the data.
* \param aPeriodVector A PeriodVector passed by reference in which to insert the
*                      value.
* \param aModeltime Modeltime object.
*/
template<class T>
void XMLHelper<T>::insertValueIntoVector( const xercesc::DOMNode* aNode,
                                          objects::PeriodVector<T>& aPeriodVector,
                                          const Modeltime* aModeltime )
{

   /*! \pre Make sure we were passed a valid node reference. */
   assert( aNode );

   const int year = XMLHelper<int>::getAttr( aNode, "year" );
   // boolean to fill out the readin value to all the periods
   const bool fillout = XMLHelper<bool>::getAttr( aNode, "fillout" );

   // Check to make sure the year attribute returned non-zero.
   if( year == 0 ) {
       ILogger& mainLog = ILogger::getLogger( "main_log" );
       mainLog.setLevel( ILogger::ERROR );
       mainLog << "Year value not set for vector input tag >>"
           << XMLHelper<std::string>::safeTranscode( aNode->getNodeName() )
           << "<<" << std::endl;
       return;
   }

   int period = aModeltime->getyr_to_per( year );

   // Check that the period returned correctly.
   // Check to make sure the year attribute returned non-zero.
   if ( !( ( period >= 0 ) && ( period < static_cast<int>( aPeriodVector.size() ) ) ) ) {
       ILogger& mainLog = ILogger::getLogger( "main_log" );
       mainLog.setLevel( ILogger::ERROR );
       mainLog << "Period value is out of bounds for vector input tag >>"
           << XMLHelper<std::string>::safeTranscode( aNode->getNodeName() )
           << "<<" << std::endl;
       return;
   }

   aPeriodVector[ period ] =  XMLHelper<T>::getValue( aNode );

   if( fillout ) {
       // will not do if period is already last period or maxperiod
       for ( unsigned int i = period + 1; i < aPeriodVector.size(); ++i ) {
           aPeriodVector[ i ] =  aPeriodVector[ period ];
       }
   }
}

/*!
* \brief Return the period cooresponding to the year in the node,
* works analogous to insertValueIntoVector, returning the appropriate period
* \warning Make sure the node passed as an argument as a year attribute.
* \param node A pointer to a node from which to extract the data.
* \param modeltime A pointer to the modeltime object to use to determine the correct period.
*/

template<class T>
int XMLHelper<T>::getNodePeriod ( const xercesc::DOMNode* node, const Modeltime* modeltime ) {
    /*! \pre Make sure we were passed a valid node reference. */
    assert( node );

    const int year = XMLHelper<int>::getAttr( node, "year" );

    // Check to make sure the year attribute returned non-zero.
    assert( year != 0 );

    int period = modeltime->getyr_to_per( year );

    // Check that the period returned correctly.
    assert( ( period >= 0 ) && ( period <= modeltime->getmaxper() ) );
    return period;
}

/*! \brief Function which converts XMLCh* to a string without leaking memory.
* \details This function when passed an XMLCh* string will call the XMLString::transcode method to
* convert the string into a dynamically allocated char* buffer. The function will then
* convert the buffer into a string and free the buffer. This function should always be used instead
* of the XMLString::transcode( XMLCh* ).
* \warning Always use this function instead of XMLString::transcode( XMLCh* ) otherwise memory will leak.
* \warning This function replaces XMLString::transcode( XMLCh* ) but not XMLString::transcode( char* ).
* The latter version is used to create an XMLCh* string. This must still be done with XMLString::transcode.
* Be very careful to free memory when doing so.
* \param toTranscode string to be converted to a standard string.
* \return An STL string equivalent to the XMLCh* string passed into the function.
*/

template<class T>
std::string XMLHelper<T>::safeTranscode( const XMLCh* toTranscode ) {
   char* transcoded = xercesc::XMLString::transcode( toTranscode );
   std::string retString = transcoded;
   xercesc::XMLString::release( &transcoded );
   return retString;
}

//! Function to write the argument element to xml in proper format.
/*!
* This function is used to write a single element containing a single value and an optional year to the output stream
* in XML. If the year is not passed in, the function will not print the year attribute.
* \param value Value to print to XML.
* \param elementName Name of the element.
* \param out Stream to print to.
* \param tabs A tabs object responsible for printing the correct number of tabs.
* \param year Optional year value to print as an attribute.
* \param name Optional name value to print as an attribute.
* \param fillout Optional attribute which specifies the value should be applied to all following time periods.
*/
template<class T>
void XMLWriteElement( const T value, const std::string elementName, std::ostream& out, const Tabs* tabs, const int year = 0, const std::string name = "", const bool fillout = false ) {

   tabs->writeTabs( out );

   out << "<" << elementName;

   if ( name != "" ) {
      out << " name=\"" << name << "\"";
   }

   if( year != 0 ){
      out << " year=\"" << year << "\"";
   }
   if( fillout ){
       out << " fillout=\"" << 1 << "\"";
   }
   out << ">";

   out << value;

   out << "</" << elementName << ">" << std::endl;
}
//! Function to write the argument element to xml with a integer attribute in proper format.
/*!
* This function is used to write a single element containing a single value along with an integer attribute to the output stream
* in XML.
* \param value Value to print to XML.
* \param elementName Name of the element.
* \param out Stream to print to.
* \param tabs A tabs object responsible for printing the correct number of tabs.
* \param aAttrs Map of attribute name to attribute value.
*/
template<class T, class U>
void XMLWriteElementWithAttributes( const T value, const std::string elementName,
                                   std::ostream& out, const Tabs* tabs,
                                   const std::map<std::string, U> aAttrs )
{
    tabs->writeTabs( out );
    out << "<" << elementName;
    typedef typename std::map<std::string, U>::const_iterator MapIterator;
    for( MapIterator entry = aAttrs.begin(); entry != aAttrs.end(); ++entry ){
        out << " " << entry->first <<"=\"" << entry->second << "\"";
    }
    out << ">" << value << "</" << elementName << ">" << std::endl;
}

/*! \brief Write an element XML tag.
* \details This function is used to write an XML element tag and an optional name and tag type to the output stream.
* The name and tag type are optional attributes. The name and tag type may be left out, but the tag will not be written
* without both tag type and tag name present. If the arguments are left out, the function will not write the attribute.
* The function increases the indent level after writing the tag so that subsequent elements are correctly indented.
* \param elementName Name of the element.
* \param out Stream to print to.
* \param Tabs The number of tabs to print before the element.
* \param tagType Optional tag type value to print as an attribute.
* \param typeName Optional type name value to print as an attribute.
*/

inline void XMLWriteElementTag( const std::string& elementName, std::ostream& out, Tabs* tabs, const std::string& tagType = "", const std::string& typeName = "")
{
   tabs->writeTabs( out );

   out << "<" << elementName;

   if ( ( typeName != "" ) &&  ( tagType != "") ){
	   out << " " << tagType << "=\"" << typeName << "\"";
		   
   }
   out << ">" << std::endl;
   tabs->increaseIndent();
}

/*!  \brief Write a closing XML tag.
* \details This function is used to write a closing XML element tag. It decreases the indent before writing the tag,
* and adds a newline.
* \note Closing tags cannot have attributes.
* \param elementName Name of the element.
* \param out Stream to print to.
* \param tabs The number of tabs to print before the element.
*/
inline void XMLWriteCloseElementTag( std::string& elementName, std::ostream& out, Tabs* tabs)
{
    tabs->decreaseIndent();
    tabs->writeTabs( out );
    out << "</" << elementName;
    out << ">" << std::endl;
}

/*! \brief Write an opening XML tag.
* \details This function is used to write an opening XML tag and an optional name and year to the output stream.
* The name and year are optional attributes. The name and year may be left out, or only the year may be left out, but
* the year cannot be written without a year unless the empty string is included in the function arguments, due
* to the way that C++ default arguments work. If the arguments are left out, the function will not write the attribute.
* The function increases the indent level after writing the tag so that subsequent elements are correctly indented.
* \param elementName Name of the element.
* \param out Stream to print to.
* \param Tabs The number of tabs to print before the element.
* \param name Optional name value to print as an attribute.
* \param year Optional year value to print as an attribute.
* \param type Optional type name to print as an attribute.
*/
inline void XMLWriteOpeningTag( const std::string& elementName, std::ostream& out, Tabs* tabs, const std::string& name = "", const int year = 0, const std::string& type = "" ) {

    tabs->writeTabs( out );
    out << "<" << elementName;

    if( year ){
        out << " year=\"" << year << "\"";
    }
    if ( name != "" ){
        out << " name=\"" << name << "\"";
    }
    if( type != "" ){
        out << " type=\"" << type << "\"";
    }
    out << ">" << std::endl;
    tabs->increaseIndent();
}

/*!  \brief Write a closing XML tag.
* \details This function is used to write a closing XML tag. It decreases the indent before writing the tag,
* and adds a newline.
* \note Closing tags cannot have attributes.
* \param elementName Name of the element.
* \param out Stream to print to.
* \param tabs The number of tabs to print before the element.
*/
inline void XMLWriteClosingTag( const std::string& elementName, std::ostream& out, Tabs* tabs ) {

    tabs->decreaseIndent();
    tabs->writeTabs( out );
    out << "</" << elementName;
    out << ">" << std::endl;
}

//! Function to write the argument element to xml in proper format if it is not equal to the default value for the element..
/*!
* This function is used to write a single element containing a single value and an optional year to the output stream
* in XML if the value is not equal to the default value.. If the year is not passed in, the function will not print the year attribute.
* \param value Value to print to XML.
* \param elementName Name of the element.
* \param out Stream to print to.
* \param tabs A tabs object responsible for printing the correct number of tabs.
* \param defaultValue Default value to compare the value to.
* \param year Optional year value to print as an attribute.
* \param name Optional name value to print as an attribute.
* \param fillout Optional boolean whether to add the fillout attribute with a true value.
*/
template<class T>
void XMLWriteElementCheckDefault( const T value, const std::string elementName, std::ostream& out, const Tabs* tabs, const T defaultValue = T(), const int year = 0, const std::string name = "", const bool fillout = false ) {
   if( !util::isEqual( value, defaultValue ) ) {
       XMLWriteElement( value, elementName, out, tabs, year, name, fillout );
   }
}

/*!
* \brief Function which writes out the values contained in a vector.
* \details This function is used to write out the values of a vector in XML format, along with their year tag.
* The function will also avoid writing out elements if they all have default values, and will collapse consecutive
* equal values into one element with a fillout attribute.
* \param aOutputVector The vector of values to write out.
* \param aElementName The elementName to write out for each value.
* \param aOut Stream to print to.
* \param aTabs A tabs object responsible for printing the correct number of tabs.
* \param modeltime A pointer to the global modeltime object.
* \param aDefaultValue Default value for items in this vector.
*/
template<class T>
void XMLWriteVector( const std::vector<T>& aOutputVector, const std::string& aElementName, std::ostream& aOut, Tabs* aTabs, const Modeltime* modeltime, const T aDefaultValue = T() ) {

    for( unsigned int i = 0; i < aOutputVector.size(); i++ ){
        // Determine the correct year.
        unsigned int year = modeltime->getper_to_yr( i );

        // Determine if we can use fillout.
        unsigned int canSkip = 0;
        for( unsigned int j = i + 1; j < aOutputVector.size(); j++ ){
            if( util::isEqual( aOutputVector.at( i ), aOutputVector.at( j ) ) ){
                canSkip++;
            }
            else {
                break;
            }
        }
        if( canSkip > 0 ){
            // Need to write out the default value because they could be cleared
            // by an earlier fillout.
            XMLWriteElement( aOutputVector.at( i ), aElementName, aOut, aTabs, year, "", true );
            // If canSkip gets to last element in vector, value of i forces breakout of first FOR loop.
            i += canSkip;
        }
        else {
            // Can't skip at all or can't skip anymore.
            // Write out default if fillout and other values are being used together,
            // or if fillout is not used and all values are unique.
            XMLWriteElement( aOutputVector.at( i ), aElementName, aOut, aTabs, year, "", false );
        }
    }
}

/*!
* \brief Function which writes out the values contained in a YearVector.
* \details This function is used to write out the values of a YearVector in XML
*          format, along with their year tag. The function will also avoid
*          writing out elements if they all have default values, and will collapse
*          consecutive equal values into one element with a fillout attribute.
* \param aOutputVector The YearVector of values to write out.
* \param aElementName The elementName to write out for each value.
* \param aOut Stream to print to.
* \param aTabs A tabs object responsible for printing the correct number of
*        tabs.
* \param aFirstYear The year of the first element in the vector.
* \param aDefaultValue Default value for items in this vector.
*/
template<class T>
void XMLWriteVector( const objects::YearVector<T>& aOutputVector,
                     const std::string& aElementName,
                     std::ostream& aOut,
                     Tabs* aTabs,
                     const unsigned int aFirstYear,
                     const T aDefaultValue = T() )
{
    for( unsigned int i = aFirstYear; i < aFirstYear + aOutputVector.size(); ++i ){
        // Determine if we can use fillout.
        unsigned int canSkip = 0;
        for( unsigned int j = i + 1; j < aFirstYear + aOutputVector.size(); ++j ){
            if( util::isEqual( aOutputVector[ i ], aOutputVector[ j ] ) ){
                ++canSkip;
            }
            else {
                break;
            }
        }
        if( canSkip > 0 ){
            // Need to write out the default value because they could be cleared
            // by an earlier fillout.
            XMLWriteElement( aOutputVector[ i ], aElementName, aOut, aTabs, i, "", true );
            // If canSkip gets to last element in vector, value of i forces breakout of first FOR loop.
            i += canSkip;
        }
        else {
            // Can't skip at all or can't skip anymore.
            // Write out default if fillout and other values are being used together,
            // or if fillout is not used and all values are unique.
            XMLWriteElement( aOutputVector[ i ], aElementName, aOut, aTabs, i, "", false );
        }
    }
}

/*!
* \brief Function which writes out the values contained in a PeriodVector.
* \details This function is used to write out the values of a PeriodVector in XML
*          format, along with their year tag. The function will also avoid
*          writing out elements if they all have default values, and will collapse
*          consecutive equal values into one element with a fillout attribute.
* \param aOutputVector The PeriodVector of values to write out.
* \param aElementName The elementName to write out for each value.
* \param aOut Stream to print to.
* \param aTabs A tabs object responsible for printing the correct number of
*        tabs.
* \param aModeltime The Modeltime object.
* \param aDefaultValue Default value for items in this vector.
*/
template<class T>
void XMLWriteVector( const objects::PeriodVector<T>& aOutputVector,
                     const std::string& aElementName,
                     std::ostream& aOut,
                     Tabs* aTabs,
                     const Modeltime* aModeltime,
                     const T aDefaultValue = T() )
{
    for( unsigned int i = 0; i < aOutputVector.size(); i++ ){
        // Determine the correct year.
        unsigned int year = aModeltime->getper_to_yr( i );

        // Determine if we can use fillout.
        unsigned int canSkip = 0;
        for( unsigned int j = i + 1; j < aOutputVector.size(); j++ ){
            if( util::isEqual( aOutputVector[ i ], aOutputVector[ j ] ) ){
                canSkip++;
            }
            else {
                break;
            }
        }
        if( canSkip > 0 ){
            // Need to write out the default value because they could be cleared
            // by an earlier fillout.
            XMLWriteElement( aOutputVector[ i ], aElementName, aOut, aTabs, year, "", true );
            // If canSkip gets to last element in vector, value of i forces breakout of first FOR loop.
            i += canSkip;
        }
        else {
            // Can't skip at all or can't skip anymore.
            // Write out default if fillout and other values are being used together,
            // or if fillout is not used and all values are unique.
            XMLWriteElement( aOutputVector[ i ], aElementName, aOut, aTabs, year, "", false );
        }
    }
}

/*!
* \brief Function to parse an XML file, returning a pointer to the root.
*
* This is a very simple function which calls the parse function and handles the exceptions which it may throw.
* It also takes care of fetching the document and its root element.
* \param aXMLFile The name of the file to parse.
* \param aModelElement Element to call XMLParse on.
* \return Whether parsing was successful.
*/

template <class T>
bool XMLHelper<T>::parseXML( const std::string& aXMLFile, IParsable* aModelElement ) {
    // Track the number of active parses to avoid destroying a document that causes other
    // documents to be parsed before its own parsing was complete.
    static unsigned int numParses = 0;
    ++numParses;
    xercesc::XercesDOMParser* parser = XMLHelper<T>::getParser();
    try {
        parser->parse( aXMLFile.c_str() );
    } catch ( const xercesc::XMLException& toCatch ) {
        std::string message = XMLHelper<std::string>::safeTranscode( toCatch.getMessage() );
        std::cout << "ERROR: XML Read Exception message is:" << std::endl << message << std::endl;
        return false;
    } catch ( const xercesc::DOMException& toCatch ) {
        std::string message = XMLHelper<std::string>::safeTranscode( toCatch.msg );
        std::cout << "ERROR: XML Read Exception message is:" << std::endl << message << std::endl;
        return false;
    } catch ( const xercesc::SAXException& toCatch ){
        std::string message = XMLHelper<std::string>::safeTranscode( toCatch.getMessage() );
        std::cout << "ERROR: XML Read Exception message is:" << std::endl << message << std::endl;
        return false;
    } catch (...) {
        std::cout << "ERROR:Unexpected XML Read Exception." << std::endl;
        return false;
    }

    bool success = aModelElement->XMLParse( parser->getDocument()->getDocumentElement() );
    // Cleanup parser memory if there are no active parses.
    if( --numParses == 0 ){
        parser->resetDocumentPool();
        parser->resetCachedGrammarPool();
    }
    return success;
}

/*! \brief Function which initializes the XML Platform and creates an instance
* of an error handler and parser.
* \note Logs are not initialized yet so they cannot be used.
* \author Josh Lurz
*/
template<class T>
void XMLHelper<T>::initParser() {
    try {
        // Initialize the Xerces platform.  Use the XQillaPlatformUtils to initialize
        // Xerces if the XMLDB is in use otherwise it will not be properly initialized
        // and will cause a crash.
#if( __USE_XML_DB__ )
        XQillaPlatformUtils::initialize();
#else
        xercesc::XMLPlatformUtils::Initialize();
#endif
    } catch ( const xercesc::XMLException& toCatch ) {
        std::string message = XMLHelper<std::string>::safeTranscode( toCatch.getMessage() );
        std::cout << "Severe error during XML Platform initialization: "<< std::endl << message << std::endl;
        exit( -1 );
    }

    // Initialize the instances of the parser and error handler.
    *getParserPointerInternal() = new xercesc::XercesDOMParser();
    (*getParserPointerInternal())->setValidationScheme( xercesc::XercesDOMParser::Val_Always );
    (*getParserPointerInternal())->setDoNamespaces( false );
    (*getParserPointerInternal())->setDoSchema( true );
    (*getParserPointerInternal())->setCreateCommentNodes( false ); // No comment nodes
    (*getParserPointerInternal())->setIncludeIgnorableWhitespace( false ); // No text nodes

    *getErrorHandlerPointerInternal() = ( (xercesc::ErrorHandler*)new xercesc::HandlerBase() );
    (*getParserPointerInternal())->setErrorHandler( *getErrorHandlerPointerInternal() );
}

/*! \brief Return the text string.
* \author Josh Lurz
* \return The #text string.
*/
template<class T>
const std::string& XMLHelper<T>::text(){
    const static std::string TEXT = "#text";
    return TEXT;
}

/*! \brief Return the name string.
* \author Josh Lurz
* \return The name string.
*/
template<class T>
const std::string& XMLHelper<T>::name(){
    const static std::string NAME = "name";
    return NAME;
}

/*! \brief Function which returns a pointer to a XercesDOMParser*.
* \details This function first checks if the parser has already been initialized.
* If it hasn't, it initializes the parser. It then returns a pointer to the parser.
* \author Josh Lurz
* \warning The user must call cleanupParser after the parser is finished being used
* to prevent a memory leak.
* \return A pointer to a XercesDOMParser.
*/
template<class T>
xercesc::XercesDOMParser* XMLHelper<T>::getParser() {
    // If the parser has not been initialized already, initialize it.
    if( !(*getParserPointerInternal()) ){
        initParser();
    }

    // Return a pointer to the already initialized parser.
    return *getParserPointerInternal();
}

/*! \brief Function which cleans up the memory used by the XML Parser.
* \details This function deletes the parser, errorhandler, and instructs
* the XMLPlatform to free its memory.
* \author Josh Lurz
* \warning This function must be called if getParser is ever called.
*/
template<class T>
void XMLHelper<T>::cleanupParser(){
    delete *getErrorHandlerPointerInternal();
    delete *getParserPointerInternal();
// The XML database will terminate xerces if it is compiled in. Terminating
// twice will cause a crash.
#if( __USE_XML_DB__ )
    XQillaPlatformUtils::terminate();
#else
    xercesc::XMLPlatformUtils::Terminate();
#endif
}

/*! \brief Reset the name to number mapping for a vector to the current names and numbers of the map.
* \details This function is used to reset and update a map to contain the correct name to index mapping for a
* vector of items.
* \note T must support the getName function.
* \author Josh Lurz
*/
template<class T>
static void resetMapIndices( const std::vector<T>& aItems, std::map<std::string, int>& aIndiceMap ){
    aIndiceMap.clear();
    for( int i = 0; i < static_cast<int>( aItems.size() ); ++i ){
        aIndiceMap[ aItems[ i ]->getName() ] = i;
    }
}

/*! \brief Function which positions a node containing children within a vector
*          of its parent, and if successful parses the new node.
* \details This function will look at the name and delete attributes of the node
*          to determine if the model node which corresponds to the input should
*          be added, modified, or deleted. After it determines this it will make
*          this change to the data structure passed in.
* \param aNode The node pointing to the container node in the XML tree.
* \param aContainerSet The vector of objects of the type pointed to by node.
* \param aNewObject An object to use if the xml node is unique. This node is
*        deleted if it is not needed.
* \param aIDAttr Name of the attribute containing the unique identifier for the
*        node. Defaults to "name".
* \return Whether the node was parsed successfully.
*/
template<class T, class U>
bool parseContainerNode( const xercesc::DOMNode* aNode,
                         std::vector<T*>& aContainerSet,
                         U* aNewObject,
                         const std::string& aIDAttr = XMLHelper<std::string>::name() )
{
    assert( aNode );
    // Force U* to implement IParsable. TODO: Activate this check when all
    // elements of the model properly implement the interface. IParsable*
    // parsable = aNewNode;

    // Have an auto_ptr keep the new memory.
    std::auto_ptr<U> newObjWrapper( aNewObject );

    // First determine if the node exists.
    const std::string objName = XMLHelper<std::string>::getAttr( aNode, aIDAttr );

    // Search the insert to vector for an item with the name.
    typename std::vector<T*>::iterator iter = util::searchForValue( aContainerSet, objName );

    // Determine if we should be deleting a node.
    bool shouldDelete = XMLHelper<bool>::getAttr( aNode, "delete" );

    // Check if the node already exists in the model tree.
    if( iter != aContainerSet.end() ){
        // Modify or delete the node based on the contents of the delete attribute.
        if( shouldDelete ) {
            // Perform deletion.
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::DEBUG );
            mainLog << "Deleting node: " << objName << std::endl;

            // Clean up the memory the vector points at.
            delete *iter;

            // Remove the pointer from the vector.
            aContainerSet.erase( iter );
        }
        // Parse the XML data into the existing node to modify it.
        else {
           (*iter)->XMLParse( aNode );
        }
    }
    // The node does not already exist.
    else if( shouldDelete ) {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "Could not delete node " << objName << " as it does not exist." << std::endl;
    }
    else if( XMLHelper<bool>::getAttr( aNode, "nocreate" ) ) {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Did not create node " << objName << " as the nocreate input flag was set." << std::endl;
        XMLHelper<void>::printXMLTrace( aNode, mainLog );
    }
    else {
        aNewObject->XMLParse( aNode );
        aContainerSet.push_back( newObjWrapper.release() );
    }

    // TODO: Setup error checking. This currently won't work because the return
    //       type of XMLParse can be void or bool.
    return true;
}

/*!
 * \brief Parse the name attribute of a node.
 * \param aNode Node to parse.
 * \return The name attribute of the node or the empty string if it does not exist.
 */
static const std::string parseNameAttr( const xercesc::DOMNode* aNode ) {
    return XMLHelper<std::string>::getAttr( aNode,
                                            XMLHelper<std::string>::name() );
}

/*!
 * \brief Type definition for a function pointer that parses the identifier of a
 *        node.
 */
typedef const std::string (*AttrGetter) (const xercesc::DOMNode* aNode);

/*! \brief Function which positions a node containing children within a vector
*          of its parent, and if successful parses the new node.
* \details This function will look at the name and delete attributes of the node
*          to determine if the model node which corresponds to the input should
*          be added, modified, or deleted. After it determines this it will make
*          this change to the data structure passed in.
* \param aNode The node pointing to the container node in the XML tree.
* \param aContainerSet The vector of objects of the type pointed to by node.
* \param aNewObject An object to use if the xml node is unique. This node is
*        deleted if it is not needed.
* \param aAttrGetter Pointer to a function which parses the ID attribute from
*        the current node.
*/
template<class T, class U>
void parseContainerNode( const xercesc::DOMNode* aNode,
                         std::vector<T>& aContainerSet,
                         U& aNewObject,
                         AttrGetter aAttrGetter = &parseNameAttr )
{
    const std::string aObjName = aAttrGetter( aNode );

    assert( aNode );
    // Force U* to implement IParsable. TODO: Activate this check when all
    // elements of the model properly implement the interface. IParsable*
    // parsable = aNewNode;

    // First determine if the node exists.

    // Search the insert to vector for an item with the name.
    typename std::vector<T>::iterator iter = util::searchForValue( aContainerSet, aObjName );

    // Determine if we should be deleting a node.
    bool shouldDelete = XMLHelper<bool>::getAttr( aNode, "delete" );

    // Check if the node already exists in the model tree.
    if( iter != aContainerSet.end() ){
        // Modify or delete the node based on the contents of the delete attribute.
        if( shouldDelete ) {
            // Perform deletion.
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::DEBUG );
            mainLog << "Deleting node: " << aObjName << std::endl;

            // Remove the pointer from the vector.
            aContainerSet.erase( iter );
        }
        // Parse the XML data into the existing node to modify it.
        else {
           (*iter)->XMLParse( aNode );
        }
    }
    // The node does not already exist.
    else {
        if( shouldDelete ) {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::ERROR );
            mainLog << "Could not delete node " << aObjName << " as it does not exist." << std::endl;
        }
        else if( XMLHelper<bool>::getAttr( aNode, "nocreate" ) ) {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::NOTICE );
            mainLog << "Did not create node " << aObjName << " as the nocreate input flag was set." << std::endl;
            XMLHelper<void>::printXMLTrace( aNode, mainLog );
        }
        else {
            aNewObject->XMLParse( aNode );
            aContainerSet.push_back( aNewObject );
        }
    }
}

/*! \brief Function which parses a node containing model-children, such as a region, and determines what to do with it.
* \details This function will look at the name and delete attributes of the node to determine if the model node which
* corresponds to the input should be added, modified, or deleted. After it determines this it will make this change
* to the model tree.
* \param node The node pointing to the container node in the XML tree.
* \param insertToVector The vector of objects of the type pointed to by node.
* \param corrMap The map of node name attributes to locations within insertToVector.
* \param newNode An object to use if the xml node is unique. This node is deleted if it is not needed.
* \param attrName Name of the attribute which contains the name identifier.
* \todo Remove this function and use only the above variant once the entire model is converted.
*/
template<class T, class U>
void parseContainerNode( const xercesc::DOMNode* node, std::vector<U>& insertToVector,
                         std::map<std::string,int>& corrMap, T* newNode, const std::string& attrName = "name" )
{
    assert( node );
    // Have an auto_ptr keep the new memory.
    std::auto_ptr<T> newNodePtr( newNode );

    // First determine if the node exists.
    const std::string objName = XMLHelper<std::string>::getAttr( node, attrName );
    std::map<std::string,int>::const_iterator iter = corrMap.find( objName );

    // Determine if we should be deleting a node.
    bool shouldDelete = XMLHelper<bool>::getAttr( node, "delete" );

    // Check if the node already exists in the model tree.
    if( iter != corrMap.end() ){
        // Modify or delete the node based on the contents of the delete attribute.
        if( shouldDelete ) {
            // Perform deletion.
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::DEBUG );
            mainLog << "Deleting node: " << objName << std::endl;

            // Create an iterator which points at the location which should be deleted.
            typedef typename std::vector<U>::iterator VectorIterator;
            VectorIterator delIter = insertToVector.begin() + iter->second;
            // Clean up the memory the vector points at.
            delete *delIter;
            // Remove the pointer from the vector.
            insertToVector.erase( delIter );

            // Now reset the map. There is probably a more efficient way to do this.
            resetMapIndices( insertToVector, corrMap );
        }
        // Otherwise modify node.
        else {
           insertToVector[ iter->second ]->XMLParse( node );
        }
    }
    // The node does not already exist.
    else {
        if( shouldDelete ) {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::ERROR );
            mainLog << "Could not delete node " << objName << " as it does not exist." << std::endl;
        }
        else if( XMLHelper<bool>::getAttr( node, "nocreate" ) ) {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::NOTICE );
            mainLog << "Did not create node " << objName << " as the nocreate input flag was set." << std::endl;
            XMLHelper<void>::printXMLTrace( node, mainLog );
        }
        else {
            newNode->XMLParse( node );
            insertToVector.push_back( newNodePtr.release() );
            corrMap[ newNode->getName() ] = static_cast<int>( insertToVector.size() ) - 1;
        }
    }
}

/*! \brief Parses a single container node held by an auto_ptr.
* \details Looks at the name and delete attributes of the node to determine if
*          the model node which corresponds to the input should be added,
*          modified, or deleted and makes the corresponding change to the model
*          tree.
* \param aNode The node pointing to the container node in the XML tree.
* \param aContainer The auto_ptr which contains the model's corrsponding object.
* \param aNewNode An object to insert into the model tree if the node is unique.
*        This node is deleted if it is not needed.
* \param aAttrName Name of the attribute which contains the name identifier.
* \sa parseContainerNode
*/
template<class T, class U>
void parseSingleNode( const xercesc::DOMNode* aNode, std::auto_ptr<U>& aContainer,
                      T* aNewNode, const std::string& aAttrName = "name" )
{
    assert( aNode );
    assert( aNewNode );

    // Have an auto_ptr keep the new memory.
    std::auto_ptr<T> newNodePtr( aNewNode );

    // Determine if we should be deleting a node.
    bool shouldDelete = XMLHelper<bool>::getAttr( aNode, "delete" );

    // Check if the container has already been created.
    if( aContainer.get() ){
        // Modify or delete the node based on the contents of the delete attribute.
        if( shouldDelete ) {
            // Perform deletion.
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::DEBUG );
            mainLog << "Deleting node: " << XMLHelper<std::string>::getAttr( aNode, aAttrName ) << std::endl;
            aContainer.reset( 0 );
        }
        // Otherwise modify node.
        else {
           aContainer->XMLParse( aNode );
        }
    }
    // The node does not already exist.
    else {
        if( shouldDelete ) {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::ERROR );
            mainLog << "Could not delete node " << XMLHelper<std::string>::getAttr( aNode, aAttrName )
                    << " as it does not exist." << std::endl;
        }
        else if( XMLHelper<bool>::getAttr( aNode, "nocreate" ) ) {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::NOTICE );
            mainLog << "Did not create node " << XMLHelper<std::string>::getAttr( aNode, aAttrName )
                    << " as the nocreate input flag was set." << std::endl;
            XMLHelper<void>::printXMLTrace( aNode, mainLog );
        }
        else {
            aContainer = newNodePtr;
            aContainer->XMLParse( aNode );
        }
    }
}

template<class T>
xercesc::XercesDOMParser** XMLHelper<T>::getParserPointerInternal(){
    static xercesc::XercesDOMParser* parser;
    return &parser;
}

template<class T>
xercesc::ErrorHandler** XMLHelper<T>::getErrorHandlerPointerInternal(){
    static xercesc::ErrorHandler* errorHandler;
    return &errorHandler;
}


/*!
 * \brief Print a trace of the XML ancestors to the give node.
 * \details For each ancestor the node name and any attributes will be printed
 *          to the given output stream to be able to better identify where the
 *          node comes from.  This information would be useful to go along with
 *          error messages during parsing.  The file from which the node comes
 *          from is also printed.
 * \param aNode The DOM node to print hierarchy information for.
 * \param aOut The output stream to write messages to.
 */
template<class T>
void XMLHelper<T>::printXMLTrace( const xercesc::DOMNode* aNode, std::ostream& aOut ){
    // start with the parent of the current node? should I start at the given node?
    const xercesc::DOMNode* currNode = aNode->getParentNode();
    
    // Process all the way up the tree until we get to the document node.
    while( currNode->getNodeType() != xercesc::DOMNode::DOCUMENT_NODE ) {
        // print the node name and all attributes
        aOut << "\tat " << safeTranscode( currNode->getNodeName() ) << " ";
        xercesc::DOMNamedNodeMap* attrMap = currNode->getAttributes();
        for( int attrIndex = 0; attrIndex < attrMap->getLength(); ++attrIndex ) {
            aOut << safeTranscode( attrMap->item( attrIndex )->getNodeName() )
                 << " = " << safeTranscode( attrMap->item( attrIndex )->getNodeValue() );
            if( attrIndex != attrMap->getLength() - 1 ) {
                aOut << ", ";
            }
        }
        aOut << std::endl;
        currNode = currNode->getParentNode();
    }
    
    // Print the file this document was contained in.
    aOut << "\tin " << safeTranscode( 
      static_cast<const xercesc::DOMDocument*>( currNode )->getDocumentURI() ) << std::endl;
}
                                      

/*!
 * \brief Write the given XML node into the given stream.
 * \details This method does a very basic serialization of nodes, attributes, 
 *          and text data.  Children can be optionally processed when the deep
 *          parameter is set.
 * \param aNode The DOM node to serialize.
 * \param aOut The output stream to serialize to.
 * \param aTabs The tabs object to control the number of tabs to use for indentation.
 * \param aDeep If all children of the given node should be processed as well.
 */
template<class T>
void XMLHelper<T>::serializeNode( const xercesc::DOMNode* aNode, std::ostream& aOut,
                                  Tabs* aTabs, const bool aDeep )
{
    switch( aNode->getNodeType() ) {
        case xercesc::DOMNode::TEXT_NODE: {
            std::string textData = safeTranscode( aNode->getNodeValue() );
            boost::algorithm::trim( textData );
            aOut << textData;
            break;
        }
        case xercesc::DOMNode::ELEMENT_NODE: {
            // write the first tag
            aTabs->writeTabs( aOut );
            const std::string nodeName = safeTranscode( aNode->getNodeName() );
            aOut << '<' << nodeName;
            
            // write any attributes
            xercesc::DOMNamedNodeMap* attrMap = aNode->getAttributes();
            for( int attrIndex = 0; attrIndex < attrMap->getLength(); ++attrIndex ) {
                aOut << ' ' << safeTranscode( attrMap->item( attrIndex )->getNodeName() )
                     << "=\"" << safeTranscode( attrMap->item( attrIndex )->getNodeValue() ) << '"';
            }
            
            if( aDeep ) {
                // Process children as well so go ahead and close this tag then recursively
                // process them.
                aOut << '>';
                bool onlyTextChildren = true;
                xercesc::DOMNodeList* childNodes = aNode->getChildNodes();
                for( int childIndex = 0; childIndex < childNodes->getLength() && onlyTextChildren; ++childIndex ) {
                    if( childNodes->item( childIndex )->getNodeType() != xercesc::DOMNode::TEXT_NODE ) {
                        onlyTextChildren = false;
                    }
                }
                if( !onlyTextChildren ) {
                    aOut << std::endl;
                    aTabs->increaseIndent();
                }
                
                for( int childIndex = 0; childIndex < childNodes->getLength(); ++childIndex ) {
                    serializeNode( childNodes->item( childIndex ), aOut, aTabs, aDeep );
                }
                
                if( !onlyTextChildren ) {
                    // write the closing tag
                    XMLWriteClosingTag( nodeName, aOut, aTabs );
                }
                else {
                    // close the tag on the same line
                    aOut << "</" << nodeName << '>' << std::endl;
                }
            }
            else {
                // just put a close tag on this to make it valid XML
                aOut << " />" << std::endl;
            }
            break;
        }
            
        // otherwise not handled by this serializer
    }
}

#endif // _XML_HELPER_H_
