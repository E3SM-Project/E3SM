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
 * \file s_curve_shutdown_decider.cpp
 * \ingroup Objects
 * \brief S_CurveShutdownDecider class source file.
 * \author Patrick Luckow
 */

#include "util/base/include/definitions.h"
#include <cassert>
#include <cmath>

#include "technologies/include/s_curve_shutdown_decider.h"
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include "util/base/include/xml_helper.h"
#include "containers/include/scenario.h"
#include "util/base/include/model_time.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

/*!
 * \brief Constructor
 */
S_CurveShutdownDecider::S_CurveShutdownDecider()
: mSteepness( 0.1 ), mHalfLife( 45 )
{}

S_CurveShutdownDecider* S_CurveShutdownDecider::clone() const {
    return new S_CurveShutdownDecider( *this );
}

bool S_CurveShutdownDecider::isSameType( const string& aType ) const {
    return aType == getXMLNameStatic();
}

const string& S_CurveShutdownDecider::getName() const {
    return getXMLNameStatic();
}

/*! \brief Get the XML node name in static form for comparison when parsing XML.
* \details This public function accesses the private constant string, XML_NAME.
*          This way the tag is always consistent for both read-in and output and
*          can be easily changed. The "==" operator that is used when parsing,
*          required this second function to return static.
* \note A function cannot be static and virtual.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME as a static.
*/
const string& S_CurveShutdownDecider::getXMLNameStatic() {
    const static string XML_NAME = "s-curve-shutdown-decider";
    return XML_NAME;
}

bool S_CurveShutdownDecider::XMLParse( const xercesc::DOMNode* node ){
    // Assume we have a valid node.
    assert( node );

    const xercesc::DOMNodeList* nodeList = node->getChildNodes();
    for( unsigned int i = 0; i < nodeList->getLength(); i++ ) {
        const xercesc::DOMNode* curr = nodeList->item( i );
        if( curr->getNodeType() != xercesc::DOMNode::ELEMENT_NODE ){
            continue;
        }
        const string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );
        if( nodeName == "steepness" ){
            mSteepness = XMLHelper<double>::getValue( curr );
        }
        else if( nodeName == "half-life" ){
            mHalfLife = XMLHelper<double>::getValue( curr );
        }
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::ERROR );
            mainLog << "Unknown tag " << nodeName << " encountered while processing "
                    << getXMLNameStatic() << endl;
        }
    }

    return true;
}

void S_CurveShutdownDecider::toInputXML( ostream& aOut,
                                        Tabs* aTabs ) const
{
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs );
    XMLWriteElementCheckDefault( mSteepness, "steepness", aOut, aTabs, 0.0 );
    XMLWriteElementCheckDefault( mHalfLife, "half-life", aOut, aTabs, 0.0 );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

void S_CurveShutdownDecider::toDebugXML( const int aPeriod,
                                        ostream& aOut,
                                        Tabs* aTabs ) const
{
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs );
    XMLWriteElementCheckDefault( mSteepness, "steepness", aOut, aTabs, 0.0 );
    XMLWriteElementCheckDefault( mHalfLife, "half-life", aOut, aTabs, 0.0 );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

double S_CurveShutdownDecider::calcShutdownCoef( const ProductionFunctionInfo* aFuncInfo,
                                                const double aCalculatedProfitRate,
                                                const string& aRegionName,
                                                const string& aSectorName,
                                                const int aInstallationYear,
                                                const int aPeriod ) const 
{
    const Modeltime* modeltime = scenario->getModeltime();

    // shutdown the production with the function 1/(1+e^(steepness*(years active-halflife))).
    // this prevents a lot of retirement in the early years, and gives the vintage a long tail;
    // All remaining vintage is cut off at the read in lifetime
    return 1 / ( 1 + exp( mSteepness*( (modeltime->getper_to_yr( aPeriod ) - aInstallationYear) - mHalfLife ) ) );
}
