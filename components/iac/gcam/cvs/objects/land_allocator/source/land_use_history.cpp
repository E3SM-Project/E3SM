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
 * \file land_use_history.cpp
 * \ingroup Objects
 * \brief LandUseHistory class source file.
 * \author Josh Lurz
 */

#include "util/base/include/definitions.h"
#include <xercesc/dom/DOMNodeList.hpp>
#include "util/base/include/xml_helper.h"
#include "land_allocator/include/land_use_history.h"
#include "util/base/include/ivisitor.h"

using namespace std;
using namespace xercesc;

//! Map type for land allocations by year.
typedef std::map<unsigned int, double> LandMapType;

const string& LandUseHistory::getXMLNameStatic(){
    static const string XML_NAME = "land-use-history";
    return XML_NAME;
}

/*!
 * \brief Constructor.
 */
LandUseHistory::LandUseHistory():
mHistoricAboveGroundCarbonDensity( 0.0 ),
mHistoricBelowGroundCarbonDensity( 0.0 )
{
}

LandUseHistory::LandUseHistory(const LandUseHistory &aLandUseHistory){
	this->mHistoricalLand = aLandUseHistory.mHistoricalLand;
}


bool LandUseHistory::XMLParse( const xercesc::DOMNode* aNode ){

    // assume we are passed a valid node.
    assert( aNode );

    // get all the children.
    DOMNodeList* nodeList = aNode->getChildNodes();
    
    for( unsigned int i = 0; i < nodeList->getLength(); ++i ){
        const DOMNode* curr = nodeList->item( i );
        const string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );

        if( nodeName == "#text" ) {
            continue;
        }
        if( nodeName == "allocation" ){
            unsigned int year = XMLHelper<unsigned int>::getAttr( curr, "year" );
            if( year == 0 ){
                ILogger& mainLog = ILogger::getLogger( "main_log" );
                mainLog.setLevel( ILogger::WARNING );
                mainLog << "Land allocations must have a year attribute." << endl;
            }
            else {
                mHistoricalLand[ year ] = XMLHelper<double>::getValue( curr );
                if ( mHistoricalLand[ year ]  <  0 ) {
                    ILogger& mainLog = ILogger::getLogger( "main_log" );
                    mainLog.setLevel( ILogger::ERROR );
                    mainLog << "Negative land allocation in land-use-history. Resetting to zero." << endl;
                    mHistoricalLand[ year ] = 0;
                }
            }
        }
        else if( nodeName == "above-ground-carbon-density" ) {
            mHistoricAboveGroundCarbonDensity = XMLHelper<double>::getValue( curr );
        }
        else if( nodeName == "below-ground-carbon-density" ) {
            mHistoricBelowGroundCarbonDensity = XMLHelper<double>::getValue( curr );
        }
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Unrecognized text string: " << nodeName << " found while parsing "
                    << getXMLNameStatic() << "." << endl;
        }
    }
    // TODO: Improved error checking.
    return true;
}

const string& LandUseHistory::getName() const {
    return getXMLNameStatic();
}

void LandUseHistory::toInputXML( ostream& aOut,
                                    Tabs* aTabs ) const
{
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs );
    for( LandMapType::const_iterator i = mHistoricalLand.begin();
         i != mHistoricalLand.end(); ++i )
    {
        XMLWriteElement( i->second, "allocation", aOut, aTabs, i->first );
    }
    XMLWriteElementCheckDefault( mHistoricAboveGroundCarbonDensity, "above-ground-carbon-density", aOut, aTabs, 0.0 );
    XMLWriteElementCheckDefault( mHistoricBelowGroundCarbonDensity, "below-ground-carbon-density", aOut, aTabs, 0.0 );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

void LandUseHistory::toDebugXML( const int aPeriod,
                                    ostream& aOut,
                                    Tabs* aTabs ) const
{
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs );
    for( LandMapType::const_iterator i = mHistoricalLand.begin();
         i != mHistoricalLand.end(); ++i )
    {
        XMLWriteElement( i->second, "allocation", aOut, aTabs, i->first );
    }
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

void LandUseHistory::accept( IVisitor* aVisitor,
                             const int aPeriod ) const 
{
    // This is not called because land-use history information is printed out in the leaf instead
}

/*!
 * \brief Get the earliest year with a historical land allocation.
 * \details Gets the earliest year a historical land allocation has been set
 *          for. Returns 0 if there are none.
 * \return The earliest year with a historical land allocation.
 */
unsigned int LandUseHistory::getMinYear() const {
    if( mHistoricalLand.empty() ){
        return 0;
    }
    // The first year is in the begin iterator. This must be dereferencable if
    // the map is not empty.
    return mHistoricalLand.begin()->first;
}

/*!
 * \brief Get the last year with a historical land allocation.
 * \details Gets the last year a historical land allocation has been set for.
 *          Returns 0 if there are none.
 * \return The last year with a historical land allocation.
 */
unsigned int LandUseHistory::getMaxYear() const {
    if( mHistoricalLand.empty() ){
        return 0;
    }
    // The last year is in the reverse begin iterator. This must be
    // dereferencable if the map is not empty.
    return mHistoricalLand.rbegin()->first;
}

/*! 
 * \brief Get the land allocation for a given year.
 * \details Gets the historical land allocation for a given year. If the year is
 *          outside the range of the first and last year this will return the
 *          closest known year. If the year was read-in that value will be
 *          returned. Otherwise a linear interpolation will be performed using
 *          the closest values above and below the year.
 * \param aYear Year for which to get the historical land allocation.
 * \return Historical land allocation for the year.
 */
double LandUseHistory::getAllocation( const unsigned int aYear ) const {
    // If there are no values in the map always return zero.
    if( mHistoricalLand.empty() ){
        return 0;
    }

    // If the year is above the maximum year, assume
    // constant allocations after the last year.
    if( aYear >= getMaxYear() ){
        return mHistoricalLand.rbegin()->second;
    }

    // If the year is below the minimum year, assume
    // constant allocations before the first year.
    if( aYear <= getMinYear() ){
        return mHistoricalLand.begin()->second;
    }

    // Find the first element with a key either equal or greater than the year.
    // Note: The STL function is named strangely. This is the upper bound.
    LandMapType::const_iterator upperBound = mHistoricalLand.lower_bound( aYear );

    // Check if the upper bound is one of the read-in years.
    if( upperBound->first == aYear ){
        return upperBound->second;
    }
    
    // Find the previous element from the upper bound.
    LandMapType::const_iterator lowerBound = upperBound;
    --lowerBound;

    // Perform a linear interpolation to find the correct value.
    return util::linearInterpolateY( aYear, lowerBound->first, upperBound->first,
                                     lowerBound->second, upperBound->second );

}

/*! 
 * \brief Get the above ground carbon density.
 * \details Gets the historical carbon density
 * \return Historical above ground carbon density
 */
double LandUseHistory::getHistoricAboveGroundCarbonDensity( ) const {
    return mHistoricAboveGroundCarbonDensity;
}

/*! 
 * \brief Get the above ground carbon density.
 * \details Gets the historical carbon density
 * \return Historical above ground carbon density
 */
double LandUseHistory::getHistoricBelowGroundCarbonDensity( ) const {
    return mHistoricBelowGroundCarbonDensity;
}

const LandMapType LandUseHistory::getHistoricalLand() const {
	return mHistoricalLand;
}

void LandUseHistory::printHistory() const {
	 for( LandMapType::const_iterator i = mHistoricalLand.begin();
         i != mHistoricalLand.end(); ++i )
    {
		std::cout<<i->second<<" "<<i->first<<endl;
    }
}
