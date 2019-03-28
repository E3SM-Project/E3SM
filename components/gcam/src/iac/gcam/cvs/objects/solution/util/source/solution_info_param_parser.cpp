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
 * \file solution_info_param_parser.cpp
 * \ingroup Objects
 * \brief SolutionInfoParamParser class source file.
 * \author Pralit Patel
 */

#include "util/base/include/definitions.h"
#include <string>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include <xercesc/dom/DOMNamedNodeMap.hpp>

#include "solution/util/include/solution_info_param_parser.h"
#include "marketplace/include/marketplace.h"
#include "util/base/include/xml_helper.h"
#include "util/logger/include/ilogger.h"
#include "containers/include/iinfo.h"

using namespace std;
using namespace xercesc;

typedef map<pair<string, string>, SolutionInfoParamParser::SolutionInfoValues>::const_iterator
    CMarketIterator;
typedef vector<SolutionInfoParamParser::SolutionInfoValues*>::iterator ValueFilloutIterator;

SolutionInfoParamParser::SolutionInfoParamParser() {
}

SolutionInfoParamParser::SolutionInfoValues::SolutionInfoValues()
:mSolutionTolerance( 0 ),
mSolutionFloor( 0 ),
mBracketInterval( 0 ),
mMaxNRPriceJump( 0 ),
mDeltaPrice( 0 )
{
}

SolutionInfoParamParser::~SolutionInfoParamParser() {
}

const string& SolutionInfoParamParser::getXMLNameStatic() {
    const static string XML_NAME = "solution-info-param-parser";
    return XML_NAME;
}

bool SolutionInfoParamParser::XMLParse( const DOMNode* aNode ) {
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
        else if( nodeName == "solution-tolerance" ) {
            double value = XMLHelper<double>::getValue( curr );
            vector<SolutionInfoValues*> valuesToSet = getSolutionInfoValuesFromAttrs( curr );
            
            for( ValueFilloutIterator valueIter = valuesToSet.begin(); valueIter != valuesToSet.end(); ++valueIter ) {
                (*valueIter)->mSolutionTolerance = value;
            }
        }
        else if( nodeName == "solution-floor" ) {
            double value = XMLHelper<double>::getValue( curr );
            vector<SolutionInfoValues*> valuesToSet = getSolutionInfoValuesFromAttrs( curr );
            
            for( ValueFilloutIterator valueIter = valuesToSet.begin(); valueIter != valuesToSet.end(); ++valueIter ) {
                (*valueIter)->mSolutionFloor = value;
            }
        }
        else if( nodeName == "bracket-interval" ) {
            double value = XMLHelper<double>::getValue( curr );
            vector<SolutionInfoValues*> valuesToSet = getSolutionInfoValuesFromAttrs( curr );
            
            for( ValueFilloutIterator valueIter = valuesToSet.begin(); valueIter != valuesToSet.end(); ++valueIter ) {
                (*valueIter)->mBracketInterval = value;
            }
        }
        else if( nodeName == "max-price-change" ) {
            double value = XMLHelper<double>::getValue( curr );
            vector<SolutionInfoValues*> valuesToSet = getSolutionInfoValuesFromAttrs( curr );
            
            for( ValueFilloutIterator valueIter = valuesToSet.begin(); valueIter != valuesToSet.end(); ++valueIter ) {
                (*valueIter)->mMaxNRPriceJump = value;
            }
        }
        else if( nodeName == "delta-price" ) {
            double value = XMLHelper<double>::getValue( curr );
            vector<SolutionInfoValues*> valuesToSet = getSolutionInfoValuesFromAttrs( curr );
            
            for( ValueFilloutIterator valueIter = valuesToSet.begin(); valueIter != valuesToSet.end(); ++valueIter ) {
                (*valueIter)->mDeltaPrice = value;
            }
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

/*!
 * \brief An accessor method so that a solver info can get direct access to their
 *        solution info params.
 * \details A SolutionInfo will be responsible for setting the values as well as
 *          checking defaults where appropriate.  The values that will be returned
 *          will be a merged result of values in ascending precedence: matching
 *          market type/global(no region name), market type/region name, and good
 *          name/region name.
 * \param aGoodName The name of the good to get the info for.
 * \param aRegionName The name of the region the good is in.
 * \param aMarketType The market type string for the solution info.
 * \param aPeriod The current period for which we are getting values for.
 * \return A copy of the of a solution info values with merged values should there
 *         be entries by both market type and good name.
 */
SolutionInfoParamParser::SolutionInfoValues SolutionInfoParamParser::getSolutionInfoValuesForMarket( const string& aGoodName,
                                                                                                     const string& aRegionName,
                                                                                                     const string& aMarketType,
                                                                                                     const int aPeriod ) const
{
    // get the map for the current period
    map<pair<string, string>, SolutionInfoValues> currMap = mSolutionInfoParams[ aPeriod ];
    
    // start with an empty SolutionInfoValues
    SolutionInfoValues retMergedValues;
    
    // attempt to find a match for the market type in a global context (empty region name)
    pair<string, string> marketNamePair( aMarketType, "" );
    CMarketIterator valuesIter = currMap.find( marketNamePair );
    if( valuesIter != currMap.end() ) {
        retMergedValues.mergeValues( (*valuesIter).second );
    }
    
    // next try for the market type and the given region name
    marketNamePair.second = aRegionName;
    valuesIter = currMap.find( marketNamePair );
    if( valuesIter != currMap.end() ) {
        retMergedValues.mergeValues( (*valuesIter).second );
    }
    
    // finally look up using the good name
    marketNamePair.first = aGoodName;
    valuesIter = currMap.find( marketNamePair );
    if( valuesIter != currMap.end() ) {
        retMergedValues.mergeValues( (*valuesIter).second );
    }
    
    return retMergedValues;
}
    

/*!
 * \brief Gets the good or market type, region, and period from the give xml node and
 *        returns the SolutionInfoValues for the parsed values.
 * \details The is a convenience method to parse the appropriate attributes and
 *          do the lookup into mSolutionInfoParams and return the reference so that
 *          XMLParse has to do minimal work to set the appropriate parameter value.
 * \param aNode The xml node to parse attributes from.
 * \return A reference into the data structure so that the value can be set.  Multiple
 *         values would be returned if the fillout attribute was set.
 */
vector<SolutionInfoParamParser::SolutionInfoValues*> SolutionInfoParamParser::getSolutionInfoValuesFromAttrs( const DOMNode* aNode ) {
    DOMNamedNodeMap* attributeList = aNode->getAttributes();
    
    string goodName;
    string regionName;
    string marketType;
    int period(0);
    bool fillout(false);
    
    // List of allowed attributes, checking them through the attribute list allows us
    // to warn for any mistyped attribute names.
    for ( unsigned int i = 0; i < attributeList->getLength(); ++i ){
        DOMNode* curr = attributeList->item( i );
        string attrName = XMLHelper<string>::safeTranscode( curr->getNodeName() );
        
        if( attrName == "good" ) {
            goodName = XMLHelper<string>::getValue( curr );
        }
        else if( attrName == "market-type" ) {
            marketType = XMLHelper<string>::getValue( curr );
        }
        else if( attrName == "region" ) {
            regionName = XMLHelper<string>::getValue( curr );
        }
        else if( attrName == "period" ) {
            period = XMLHelper<int>::getValue( curr );
        }
        else if( attrName == "fillout" ) {
            fillout = XMLHelper<bool>::getValue( curr );
        }
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Unrecognized text string: " << attrName << " found while parsing "
                << getXMLNameStatic() << "." << endl;
        }
    }
    
    pair<string, string> marketPair;
    vector<SolutionInfoValues*> retValuesToSet;
    
    // Determine if the user intended this to be by good or market type.
    // Note that an empty region name is only allowed for market type since
    // that implies global.
    if( !goodName.empty() && !regionName.empty() ){
        marketPair.first = goodName;
        marketPair.second = regionName;
    }
    else if( !marketType.empty() && !regionName.empty() ){
        marketPair.first = marketType;
        marketPair.second = regionName;
    }
    else if( !marketType.empty() && regionName.empty() ){
        marketPair.first = marketType;
        // marketPair.second is empty
    }
    else{
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Did not find a good name or market type while parsing "
            << XMLHelper<string>::safeTranscode( aNode->getNodeName() ) << " in " << getXMLNameStatic() << "." << endl;
        // return an empty vector so no values will be set
        return retValuesToSet;
    }

    // get the SolutionInfoValues from the map for the parsed info
    retValuesToSet.push_back( &mSolutionInfoParams[ period ][ marketPair ] );
    
    // if fillout was set we should also include all future periods as well so that those values can also
    // get set
    for( unsigned int fillOutPer = period + 1; fillout && fillOutPer < mSolutionInfoParams.size(); ++fillOutPer ) {
        retValuesToSet.push_back( &mSolutionInfoParams[ fillOutPer ][ marketPair ] );
    }
    
    return retValuesToSet;
}

/*!
 * \brief Merges the values from a given set of values into this set of values.
 * \details If the given SolutionInfoValues has a set value (non-zero) then
 *          that value will override the current value in this SolutionInfoValues.
 * \param aMergeValues The other struct with which this one will merge with.
 */
void SolutionInfoParamParser::SolutionInfoValues::mergeValues( const SolutionInfoValues& aMergeValues ) {
    if( aMergeValues.mSolutionTolerance != 0 ) {
        this->mSolutionTolerance = aMergeValues.mSolutionTolerance;
    }
    if( aMergeValues.mSolutionFloor != 0 ) {
        this->mSolutionFloor = aMergeValues.mSolutionFloor;
    }
    if( aMergeValues.mBracketInterval != 0 ) {
        this->mBracketInterval = aMergeValues.mBracketInterval;
    }
    if( aMergeValues.mMaxNRPriceJump != 0 ) {
        this->mMaxNRPriceJump = aMergeValues.mMaxNRPriceJump;
    }
    if( aMergeValues.mDeltaPrice != 0 ) {
        this->mDeltaPrice = aMergeValues.mDeltaPrice;
    }
}
