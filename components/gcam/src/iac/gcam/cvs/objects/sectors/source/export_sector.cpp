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
* \file export_sector.cpp
* \ingroup Objects
* \brief ExportSector class source file.
* \author Josh Lurz
*/

#include "util/base/include/definitions.h"
#include <string>
#include <cassert>

// xml headers
#include <xercesc/dom/DOMNode.hpp>

#include "util/base/include/xml_helper.h"
#include "sectors/include/export_sector.h"
#include "sectors/include/subsector.h"
#include "marketplace/include/marketplace.h"
#include "util/logger/include/ilogger.h"
#include "marketplace/include/imarket_type.h"
#include "util/base/include/configuration.h"
#include "containers/include/dependency_finder.h"
#include "containers/include/scenario.h" // for marketplace
#include "containers/include/iinfo.h"
#include "sectors/include/sector_utils.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario; // for marketplace

/* \brief Constructor
* \param aRegionName The region name.
*/
ExportSector::ExportSector ( const string& aRegionName ) : SupplySector ( aRegionName ) {
    mFixedPrices.resize( scenario->getModeltime()->getmaxper() );
}

/*! \brief Return the export sector price.
* \details Export sectors currently have fixed prices, so return the read-in
*          fixed price.
* \param aPeriod Model period.
* \return Price.
*/
double ExportSector::getPrice( const GDP* aGDP, const int aPeriod ) const {
    return mFixedPrices[ aPeriod ];
}

/*! \brief Calculate the final supply price for the ExportSector, which will
*          leave the international price unchanged.
* \details Currently this function does not calculate or set a price into the
*          marketplace, so that the read-in price is preserved.
* \param aGDP The regional GDP container.
* \param aPeriod The period in which to calculate the final supply price.
*/
void ExportSector::calcFinalSupplyPrice( const GDP* aGDP, const int aPeriod ){
    // Calculate the costs for all subsectors.
    calcCosts( aPeriod );
}

/*! \brief Set the ExportSector output.
* \details Export sectors currently have only fixed outputs so this function
*          simply calculates the entire amount of fixed output and distributes
*          that amount to the subsectors.
* \author Josh Lurz
* \param aGDP GDP object uses to calculate various types of GDPs.
* \param aPeriod Model period
*/
void ExportSector::supply( const GDP* aGDP, const int aPeriod ) {
    // Do not scale supplies because this is a global market so there
    // is no way to determine a scaling factor.

    // Instruct technologies to only produce fixed output and no variable
    // output. Shares do not need to be calculated because there is
    // no variable output.
    for( unsigned int i = 0; i < subsec.size(); ++i ){
        subsec[ i ]->setOutput( 0, 1, aGDP, aPeriod );
    }

    // This is a global market so demands do not have to equal supplies locally.
    // Fixed output must equal output.
    assert( util::isEqual( getFixedOutput( aPeriod ), getOutput( aPeriod ) ) );
}

/*! \brief Get the XML node name for output to XML using a virtual method so the
*          correct name is always returned.
* \return The name of the object's associated XML element name.
*/
const std::string& ExportSector::getXMLName() const {
	return getXMLNameStatic();
}

/*! \brief Get the XML node name in static form for comparison when parsing XML.
* \return The name of the object's associated XML element name.
*/
const std::string& ExportSector::getXMLNameStatic() {
	const static string XML_NAME = "export-sector";
	return XML_NAME;
}

/*! \brief Parses any child nodes specific to derived classes
* \details Method parses any input data from child nodes that are specific to
*          the classes derived from this class. Since Sector is the generic base
*          class, there are no values here.
* \param aNodeName name of current node
* \param aCurr pointer to the current node in the XML input tree
*/
bool ExportSector::XMLDerivedClassParse( const string& aNodeName, const DOMNode* aCurr ) {
    if( aNodeName == "market" ){
		mMarketName = XMLHelper<string>::getValue( aCurr );
	}
    else if( aNodeName == "sectorprice" ){
        XMLHelper<double>::insertValueIntoVector( aCurr, mFixedPrices, scenario->getModeltime() );
    }
	else {
		return false;
	}
	return true;
}

/*! \brief Write out derived class specific class members as an input file.
* \param aOut The output stream to which to write.
* \param aTabs Object responsible for writing tabs.
*/
void ExportSector::toInputXMLDerived( ostream& aOut, Tabs* aTabs ) const {  
    // write out the market string.
    XMLWriteElement( mMarketName, "market", aOut, aTabs );
    XMLWriteVector( mFixedPrices, "sectorprice", aOut, aTabs, scenario->getModeltime(), 0.0 );
}	

/*! \brief Write out derived class specific class members for debugging.
* \param aPeriod The period in which to write debugging information.
* \param aOut The output stream to which to write.
* \param aTabs Object responsible for writing tabs.
*/
void ExportSector::toDebugXMLDerived( const int aPeriod, ostream& aOut, Tabs* aTabs ) const {
    // write out the market string.
    XMLWriteElement( mMarketName, "market", aOut, aTabs );
    XMLWriteElement( mFixedPrices[ aPeriod ], "sectorprice", aOut, aTabs );
}

/*! \brief Create new market for the ExportSector.
* \details Sets up the appropriate market within the marketplace for this
*          ExportSector. The market is not solved currently and prices are fixed
*          to read-in values.
*/
void ExportSector::setMarket() {
	// Check if the market name was not read-in or is the same as the region name.
	if( mMarketName.empty() || mMarketName == regionName ){
		ILogger& mainLog = ILogger::getLogger( "main_log" );
		mainLog.setLevel( ILogger::WARNING );
		mainLog << "Market name for export sector was unset or set to the region name." << endl;
		// Is there anything logical to do here?
	}

    Marketplace* marketplace = scenario->getMarketplace();
	// Creates a regional market. MiniCAM supply sectors are not independent and 
	// cannot be members of multi-region markets.
    if( marketplace->createMarket( regionName, mMarketName, name, IMarketType::NORMAL ) ) {
        // Set price and output units for period 0 market info
        IInfo* marketInfo = marketplace->getMarketInfo( name, regionName, 0, true );
        marketInfo->setString( "price-unit", mPriceUnit );
        marketInfo->setString( "output-unit", mOutputUnit );

        // Set the base year price which the sector reads in, into the mFixedPrices vector.
        // TODO: Separate SupplySector so this is not needed.
        mFixedPrices[ 0 ] = mBasePrice;

		// Initializes prices with any values that are read-in. 
        marketplace->setPriceVector( name, regionName, mFixedPrices );
    }
}
