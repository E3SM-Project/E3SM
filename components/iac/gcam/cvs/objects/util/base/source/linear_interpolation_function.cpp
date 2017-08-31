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
 * \file linear_interpolation_function.cpp
 * \ingroup Objects
 * \brief LinearInterpolationFunction class source file.
 * \author Pralit Patel
 */

#include "util/base/include/definitions.h"
#include <string>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>

#include "util/base/include/linear_interpolation_function.h"
#include "util/curves/include/data_point.h"
#include "util/base/include/xml_helper.h"
#include "util/logger/include/ilogger.h"

using namespace std;
using namespace xercesc;

//! Default Constructor
LinearInterpolationFunction::LinearInterpolationFunction() {
}

//! Destructor
LinearInterpolationFunction::~LinearInterpolationFunction() {
}

/*!
 * \brief The value for the xml name attribute which identifies this
 *        interpolation function.
 * \details The approach for using the name attribute value rather than
 *          the element name was taking to make generating the tags
 *          easier.
 * \return The string which identifies this function.
 * \see InterpolationFunctionFactory
 */
const string& LinearInterpolationFunction::getXMLAttrNameStatic() {
    const static string XML_NAME = "linear";
    return XML_NAME;
}

bool LinearInterpolationFunction::XMLParse( const DOMNode* aNode ) {
    // nothing to parse but will still check to make sure there is
    // nothing here and warn the user if any thing was found

    // assume we were passed a valid node.
    assert( aNode );

    // get the children of the node.
    DOMNodeList* nodeList = aNode->getChildNodes();

    // loop through the children
    for ( unsigned int i = 0; i < nodeList->getLength(); ++i ){
        DOMNode* curr = nodeList->item( i );
        string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );

        if( nodeName == "#text" ) {
            continue;
        }
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Unrecognized text string: " << nodeName << " found while parsing "
                << getXMLAttrNameStatic() << "." << endl;
        }
    }
    return true;
}

void LinearInterpolationFunction::toInputXML( ostream& aOut, Tabs* aTabs ) const {
    XMLWriteElement("", IInterpolationFunction::getXMLNameStatic(), aOut, aTabs, 0, getXMLAttrNameStatic() );
}

double LinearInterpolationFunction::interpolate( const DataPoint* aLeftPoint, const DataPoint* aRightPoint,
                                                 const double aXValue ) const
{
    // TODO: error checking
    return ( aXValue - aLeftPoint->getX() )
        * ( aRightPoint->getY() - aLeftPoint->getY() ) / ( aRightPoint->getX() - aLeftPoint->getX() )
        + aLeftPoint->getY();
}
