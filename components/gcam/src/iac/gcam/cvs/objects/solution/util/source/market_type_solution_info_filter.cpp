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
 * \file market_type_solution_info_filter.cpp
 * \ingroup Objects
 * \brief MarketTypeSolutionInfoFilter class source file.
 * \author Pralit Patel
 */

#include "util/base/include/definitions.h"
#include <string>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>

#include "solution/util/include/market_type_solution_info_filter.h"
#include "solution/util/include/solution_info.h"
#include "util/base/include/xml_helper.h"
#include "util/logger/include/ilogger.h"
#include "marketplace/include/market.h"

using namespace std;
using namespace xercesc;

MarketTypeSolutionInfoFilter::MarketTypeSolutionInfoFilter()
:mAcceptMarketType( IMarketType::END )
{
}

MarketTypeSolutionInfoFilter::~MarketTypeSolutionInfoFilter() {
}

const string& MarketTypeSolutionInfoFilter::getXMLNameStatic() {
    const static string XML_NAME = "market-type-solution-info-filter";
    return XML_NAME;
}

bool MarketTypeSolutionInfoFilter::XMLParse( const DOMNode* aNode ) {
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
        else if( nodeName == "market-type" ) {
            mAcceptMarketType = getMarketTypeFromString( XMLHelper<string>::getValue( curr ) );
        }
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Unrecognized text string: " << nodeName << " found while parsing "
                << getXMLNameStatic() << "." << endl;
        }
    }
    return true;
}

bool MarketTypeSolutionInfoFilter::acceptSolutionInfo( const SolutionInfo& aSolutionInfo ) const {
    /*!
     * \pre mAcceptMarketType was read in and set to a valid IMarketType::Type.
     */
    assert( mAcceptMarketType != IMarketType::END );
    
    return mAcceptMarketType == aSolutionInfo.getType();
}

/*!
 * \brief Convert a market type string to it's corresponding enumerated type.
 * \details Uses Market::convert_type_to_string to linearly compare strings
 *          to do the conversion.  If none of the strings matched IMarketType::END
 *          will be returned.
 * \param aMarketType The string type which should be converted to an enum.
 * \return The IMarketType::Type which corresponds to the given string, or IMarketType::END
 *         if no matches were found.
 * \see Market::convert_type_to_string
 */
IMarketType::Type MarketTypeSolutionInfoFilter::getMarketTypeFromString( const string& aMarketType ) const {
    // would it be worth it to convert this to a static method with a pre-populated map of string to enum?
    
    // iterate over each possible enum and have the Market convert it to a string
    for( IMarketType::Type currType = IMarketType::NORMAL; currType < IMarketType::END; 
            currType = IMarketType::Type( currType + 1 ) ) {
        // if the market converted string matches the given string then this is the enum
        // we are looking for
        if( aMarketType == Market::convert_type_to_string( currType ) ) {
            return currType;
        }
    }
    
    // did not find the appropriate enum warn the user
    ILogger& mainLog = ILogger::getLogger( "main_log" );
    mainLog.setLevel( ILogger::WARNING );
    mainLog << "Unrecognized market type: " << aMarketType << " found while parsing "
        << getXMLNameStatic() << "." << endl;
    return IMarketType::END;
}
