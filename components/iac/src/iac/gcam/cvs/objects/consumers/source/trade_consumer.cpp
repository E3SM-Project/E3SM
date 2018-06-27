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
* \file trade_consumer.cpp
* \ingroup Objects
* \brief The TradeConsumer class source file.
*
* \author Sonny Kim
* \author Pralit Patel
*/
#include "util/base/include/definitions.h"
#include <xercesc/dom/DOMNode.hpp>

#include "consumers/include/trade_consumer.h"
#include "util/base/include/xml_helper.h"
#include "containers/include/national_account.h"
#include "functions/include/ifunction.h"
#include "containers/include/scenario.h"
#include "util/base/include/ivisitor.h"
#include "technologies/include/ioutput.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

//! Default constructor
TradeConsumer::TradeConsumer(){
}

void TradeConsumer::copyParam( const BaseTechnology* aTradeConsumer,
                               const int aPeriod ) {
    BaseTechnology::copyParam( aTradeConsumer, aPeriod );
    aTradeConsumer->copyParamsInto( *this, aPeriod );
}

void TradeConsumer::copyParamsInto( TradeConsumer& aTradeConsumer,
                                    const int aPeriod ) const {
}

TradeConsumer* TradeConsumer::clone() const {
	return new TradeConsumer( *this );
}

//! Parse xml file for data
bool TradeConsumer::XMLDerivedClassParse( const string& aNodeName, const DOMNode* aCurr ) {
	return false;
}

//! For derived classes to output XML data
void TradeConsumer::toInputXMLDerived( ostream& aOut, Tabs* aTabs ) const {
}

//! Output debug info for derived class
void TradeConsumer::toDebugXMLDerived( const int aPeriod, ostream& aOut, Tabs* aTabs ) const {
}

//! Complete the initialization
void TradeConsumer::completeInit( const string& aRegionName,
                                  const string& aSectorName,
                                  const string& aSubsectorName )
{
	BaseTechnology::completeInit( aRegionName, aSectorName, aSubsectorName );
}

//! initialize anything that won't change during the calculation
void TradeConsumer::initCalc( const MoreSectorInfo* aMoreSectorInfo, const string& aRegionName, 
                              const string& aSectorName, NationalAccount& aNationalAccount, 
                              const Demographic* aDemographics, const double aCapitalStock, 
                              const int aPeriod ) 
{
    Consumer::initCalc( aMoreSectorInfo, aRegionName, aSectorName,
                        aNationalAccount, aDemographics, aCapitalStock,
                        aPeriod );

    // TODO: again why is this here?
    if( year == scenario->getModeltime()->getper_to_yr( aPeriod ) ) {
        // calculate Price Paid
        //BaseTechnology::calcPricePaid( aMoreSectorInfo, aRegionName, aSectorName, aPeriod );
    }
}

//! calculate income
void TradeConsumer::calcIncome( NationalAccount& aNationalAccount, const Demographic* aDemographics, 
                                const string& aRegionName, const string& aSectorName, int aPeriod ) 
{
    expenditures[ aPeriod ].reset();
    double netExport = mOutputs[ 0 ]->getCurrencyOutput( aPeriod );
	expenditures[ aPeriod ].setType( Expenditure::TOTAL_IMPORTS, netExport );

	// set National Accounts Consumption for GNP calculation
	aNationalAccount.addToAccount( NationalAccount::GNP_NOMINAL, netExport );
	aNationalAccount.addToAccount( NationalAccount::NET_EXPORT_NOMINAL, netExport );
}

//! calculate demand
void TradeConsumer::operate( NationalAccount& aNationalAccount, const Demographic* aDemographics, 
                             const MoreSectorInfo* aMoreSectorInfo, const string& aRegionName, 
                             const string& aSectorName, const bool aIsNewVintageMode, int aPeriod ) 
{
    const Modeltime* modelTime = scenario->getModeltime();
	if( year == modelTime->getper_to_yr( aPeriod ) ) {
		// calculate prices paid for consumer inputs
		BaseTechnology::calcPricePaid( aMoreSectorInfo, aRegionName, aSectorName, aPeriod,
            modelTime->gettimestep( aPeriod ) );
		calcInputDemand( expenditures[ aPeriod ].getValue( Expenditure::CONSUMPTION ), aRegionName, aSectorName, aPeriod );
		calcIncome( aNationalAccount, aDemographics, aRegionName, aSectorName, aPeriod );
        // Trade consumer does not have emissions so it does not calculate them here.
        // calculate the real amount consumed
        // TODO: this could currently just go in post calc
        aNationalAccount.addToAccount( NationalAccount::NET_EXPORT_REAL,
            calcRealGNP( aNationalAccount, aRegionName, aPeriod ) );

        mPricePaidCached = false;
	}
}

/*! \brief Get the XML node name for output to XML.
*
* This public function accesses the private constant string, XML_NAME.
* This way the tag is always consistent for both read-in and output and can be easily changed.
* This function may be virtual to be overriden by derived class pointers.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME.
*/
const string& TradeConsumer::getXMLName() const {
	return getXMLNameStatic();
}

/*! \brief Get the XML node name in static form for comparison when parsing XML.
*
* This public function accesses the private constant string, XML_NAME.
* This way the tag is always consistent for both read-in and output and can be easily changed.
* The "==" operator that is used when parsing, required this second function to return static.
* \note A function cannot be static and virtual.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME as a static.
*/
const string& TradeConsumer::getXMLNameStatic() {
    const static string XML_NAME = "tradeConsumer";
	return XML_NAME;
}

//! SGM version of outputing data to a csv file
void TradeConsumer::csvSGMOutputFile( ostream& aFile, const int aPeriod ) const {
	if ( year == scenario->getModeltime()->getper_to_yr( aPeriod ) ) {
		aFile << "***** Trade Sector Results *****" << endl << endl;
		expenditures[ aPeriod ].csvSGMOutputFile( aFile, aPeriod );
		aFile << endl;

		aFile << "Trade Consumer Expenditure" << endl << endl;
		BaseTechnology::csvSGMOutputFile( aFile, aPeriod );
	}
}

void TradeConsumer::accept( IVisitor* aVisitor, const int aPeriod ) const {
    aVisitor->startVisitTradeConsumer( this, aPeriod );
    Consumer::accept( aVisitor, aPeriod );
    aVisitor->endVisitTradeConsumer( this, aPeriod );
}
