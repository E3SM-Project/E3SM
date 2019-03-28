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
 * \file s_curve_interpolation_function.cpp
 * \ingroup Objects
 * \brief SCurveInterpolationFunction class source file.
 * \author Pralit Patel
 */

#include "util/base/include/definitions.h"
#include <string>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>

#include "util/base/include/s_curve_interpolation_function.h"
#include "util/curves/include/data_point.h"
#include "util/base/include/xml_helper.h"
#include "util/logger/include/ilogger.h"

using namespace std;
using namespace xercesc;

//! Default Constructor
SCurveInterpolationFunction::SCurveInterpolationFunction():
mSteepness( 10 ),
mMedianXValue( 0 )
{
}

/*!
 * \brief A constructor which allows a user to set the shape parameters for
 *        the s-curve directly rather than parsing them.
 * \param aSteepness The steepness to use.
 * \param aMedianXValue The median x-value to use.
 */
SCurveInterpolationFunction::SCurveInterpolationFunction( const double aSteepness,
                                                          const double aMedianXValue ):
mSteepness( aSteepness ),
mMedianXValue( aMedianXValue )
{
}

//! Destructor
SCurveInterpolationFunction::~SCurveInterpolationFunction() {
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
const string& SCurveInterpolationFunction::getXMLAttrNameStatic() {
    const static string XML_NAME = "s-curve";
    return XML_NAME;
}

bool SCurveInterpolationFunction::XMLParse( const DOMNode* aNode ) {
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
        else if( nodeName == "steepness" ) {
            mSteepness = XMLHelper<double>::getValue( curr );
        }
        else if( nodeName == "median-x-value" ) {
            mMedianXValue = XMLHelper<double>::getValue( curr );
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

void SCurveInterpolationFunction::toInputXML( ostream& aOut, Tabs* aTabs ) const {
    XMLWriteOpeningTag( IInterpolationFunction::getXMLNameStatic(), aOut, aTabs, getXMLAttrNameStatic() );

    XMLWriteElement( mSteepness, "steepness", aOut, aTabs );
    XMLWriteElement( mMedianXValue, "median-x-value", aOut, aTabs );

    XMLWriteClosingTag( IInterpolationFunction::getXMLNameStatic(), aOut, aTabs );
}

double SCurveInterpolationFunction::interpolate( const DataPoint* aLeftPoint, const DataPoint* aRightPoint,
                                                 const double aXValue ) const
{
    // TODO: error checking
    double medianX = mMedianXValue;
    if( mMedianXValue < aLeftPoint->getX() || mMedianXValue > aRightPoint->getX() ) {
        // TODO: warn?
        medianX = ( aRightPoint->getX() - aLeftPoint->getX() ) / 2 + aLeftPoint->getX();
    }
    return aRightPoint->getY() - ( aRightPoint->getY() - aLeftPoint->getY() )
        / ( 1.0 + exp( mSteepness * ( ( aXValue - medianX ) / ( aRightPoint->getX() - aLeftPoint->getX() ) ) ) );
}
