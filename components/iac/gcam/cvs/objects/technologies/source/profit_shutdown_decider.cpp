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
 * \file profit_shutdown_decider.cpp
 * \ingroup Objects
 * \brief ProfitShutdownDecider class source file.
 * \author Josh Lurz
 */

#include "util/base/include/definitions.h"
#include <cassert>
#include "technologies/include/profit_shutdown_decider.h"
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include "functions/include/function_utils.h"
#include "functions/include/ifunction.h"
#include "functions/include/inested_input.h"
#include "util/base/include/xml_helper.h"

using namespace std;
using namespace xercesc;

//! Constructor
ProfitShutdownDecider::ProfitShutdownDecider():
    mMaxShutdown(1.0),
    mSteepness(6.0),
    mMedianShutdownPoint(-0.1)
{
}

ProfitShutdownDecider* ProfitShutdownDecider::clone() const {
    return new ProfitShutdownDecider( *this );
}

bool ProfitShutdownDecider::isSameType( const std::string& aType ) const {
    return aType == getXMLNameStatic();
}

const string& ProfitShutdownDecider::getName() const {
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
const string& ProfitShutdownDecider::getXMLNameStatic() {
    const static string XML_NAME = "profit-shutdown-decider";
    return XML_NAME;
}

bool ProfitShutdownDecider::XMLParse( const xercesc::DOMNode* node ){
    
    // Assume we have a valid node.
    assert( node );

    const xercesc::DOMNodeList* nodeList = node->getChildNodes();
    for( unsigned int i = 0; i < nodeList->getLength(); i++ ) {
        const xercesc::DOMNode* curr = nodeList->item( i );
        if( curr->getNodeType() != xercesc::DOMNode::ELEMENT_NODE ){
            continue;
        }
        const string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );
        if( nodeName == "max-shutdown" ){
            mMaxShutdown = XMLHelper<double>::getValue( curr );
        }
        else if( nodeName == "steepness" ) {
            mSteepness = XMLHelper<double>::getValue( curr );
        }
        else if( nodeName == "median-shutdown-point" ) {
            mMedianShutdownPoint = XMLHelper<double>::getValue( curr );
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

void ProfitShutdownDecider::toInputXML( ostream& aOut, Tabs* aTabs ) const {
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs );
    XMLWriteElementCheckDefault( mMaxShutdown, "max-shutdown", aOut, aTabs, 1.0 );
    XMLWriteElementCheckDefault( mSteepness, "steepness", aOut, aTabs, 6.0 );
    XMLWriteElementCheckDefault( mMedianShutdownPoint, "median-shutdown-point", aOut, aTabs, -0.1 );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

void ProfitShutdownDecider::toDebugXML( const int aPeriod, ostream& aOut, Tabs* aTabs ) const {
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs );
    XMLWriteElement( mMaxShutdown, "max-shutdown", aOut, aTabs );
    XMLWriteElement( mSteepness, "steepness", aOut, aTabs );
    XMLWriteElement( mMedianShutdownPoint, "median-shutdown-point", aOut, aTabs );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

double ProfitShutdownDecider::calcShutdownCoef( const ProductionFunctionInfo* aFuncInfo,
                                                const double aCalculatedProfitRate,
                                                const string& aRegionName,
                                                const string& aSectorName,
                                                const int aInitialTechYear,
                                                const int aPeriod ) const 
{
    // Default scale factor is not to scale.
    double scaleFactor = 1;
    // There is no shutdown decision in the base period.
    if( aPeriod > 0 ){
        double profitRate;
        // Calculate the profit rate dynamically.
        if( aCalculatedProfitRate == getUncalculatedProfitRateConstant() ){
            // Calculate the profit rate.
            assert( aFuncInfo );

            double priceReceived = FunctionUtils::getPriceReceived( aRegionName,
                aSectorName, aPeriod );

            // note that the variable costs should have already been calculated and set
            // so we can get it by just calling getLevelizedCost.
            // TODO: we can not compuate it dynamically here because it is a const method and
            // calc variable costs need to set intermediate price when calculating it
            double variableCost = aFuncInfo->mNestedInputRoot->getLevelizedCost( aRegionName,
                                                                                 aSectorName,
                                                                                 aPeriod );

            // we make this profit rate relative to the priceReceived rather than
            // absolute difference since that would not be compatible with the
            // numeraire test
            profitRate = ( priceReceived - variableCost ) / priceReceived;
        }
        else {
            // Use the passed in profit rate.
            profitRate = aCalculatedProfitRate;
        }
       
        // Compute Shutdown factor using logistic S-curve.  ScaleFactor that is returned
        // is actually the fraction not shut down, so it is 1.0 - the shutdown fraction.
        const double midPointToSteepness = pow( mMedianShutdownPoint + 1, mSteepness );
        scaleFactor = 1.0 - mMaxShutdown * ( midPointToSteepness / 
                      ( midPointToSteepness + pow( profitRate + 1, mSteepness ) ) );
    }

    // Scale factor is between 0 and 1.
    assert( scaleFactor >= 0 && scaleFactor <= 1 );
    return scaleFactor;
}
