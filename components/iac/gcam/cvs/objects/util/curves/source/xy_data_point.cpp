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
* \file xy_data_point.cpp
* \ingroup Util
* \brief XYDataPoint class source file.
* \author Josh Lurz
*/

#include "util/base/include/definitions.h"
#include <iostream>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include "util/curves/include/data_point.h"
#include "util/curves/include/xy_data_point.h"
#include "util/base/include/xml_helper.h"

using namespace std;

const string XYDataPoint::XML_NAME = "XYDataPoint";

/*! \brief Constructor
* \details This is the default constructor for the XYDataPoint. The parameters are used
* to specify the initial x and y values. They default to 0.
* \param xIn The initial x value. Defaults to 0.
* \param yIn The initial y value. Defaults to 0.
*/
XYDataPoint::XYDataPoint( const double xIn, const double yIn ): x( xIn ), y( yIn ){
}

//! Destructor
XYDataPoint::~XYDataPoint(){
}

//! Clone
DataPoint* XYDataPoint::clone() const {
	return new XYDataPoint( x, y );
}

//! Static function to return the name of the XML element associated with this object.
const string& XYDataPoint::getXMLNameStatic() {
    return XML_NAME;
}

//! Return the name of the XML element associated with this object.
const string& XYDataPoint::getXMLName() const {
    return XML_NAME;
}

//! Equality
bool XYDataPoint::operator==( const XYDataPoint& rhs ) const {
    return( ( x == rhs.x ) && ( y == rhs.y ) );
}

//! Inequality
bool XYDataPoint::operator!=( const XYDataPoint& rhs ) const {
    return !( *this == rhs );
}

//! Get the x value.
double XYDataPoint::getX() const {
    return x;
}

//! Get the y value
double XYDataPoint::getY() const {
    return y;
}

//! Set the x value 
void XYDataPoint::setX( const double xValue ){
    x = xValue;
}

//! Set the y value
void XYDataPoint::setY( const double yValue ){
    y = yValue;
}

//! Print the data point.
void XYDataPoint::print( ostream& out ) const {
    out << x << "," << y << endl;
}

//! Print the data point to an XML stream.
void XYDataPoint::toInputXML( ostream& out, Tabs* tabs ) const {
    XMLWriteOpeningTag( DataPoint::getXMLNameStatic(), out, tabs, "", 0, getXMLName() );
    XMLWriteElement( x, "x", out, tabs );
    XMLWriteElement( y, "y", out, tabs );
    XMLWriteClosingTag( DataPoint::getXMLNameStatic(), out, tabs );
}

//! Parse the XYDataPoint from an XML DOM tree.
void XYDataPoint::XMLParse( const xercesc::DOMNode* node ) {
    
    xercesc::DOMNode* curr = 0;
    xercesc::DOMNodeList* nodeList; 
    string nodeName;

    // assume node is valid.
    assert( node );

    // get all children of the node.
    nodeList = node->getChildNodes();

    // loop through the children
    for ( int i = 0; i < static_cast<int>( nodeList->getLength() ); i++ ){
        curr = nodeList->item( i );
        nodeName = XMLHelper<void>::safeTranscode( curr->getNodeName() );

        // select the type of node.
        if( nodeName == "#text" ) {
            continue;
        }
        else if ( nodeName == "x" ){
            x = XMLHelper<double>::getValue( curr );
        }
        else if ( nodeName == "y" ){
            y = XMLHelper<double>::getValue( curr );
        } 
        else {
            cout << "Unrecognized text string: " << nodeName << " found while parsing " << getXMLName() << endl;
        }
    }
}

//! Switch the X and Y values.
void XYDataPoint::invertAxises(){
    swap( x, y );
}
