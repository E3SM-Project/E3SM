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
* \file invest_consumer.cpp
* \ingroup Objects
* \brief The InvestConsumer class source file.
*
* \author Sonny Kim
* \author Pralit Patel
*/
#include "util/base/include/definitions.h"
#include <xercesc/dom/DOMNode.hpp>

#include "consumers/include/invest_consumer.h"
#include "util/base/include/xml_helper.h"
#include "containers/include/national_account.h"
#include "functions/include/iinput.h"
#include "functions/include/ifunction.h"
#include "containers/include/scenario.h"
#include "util/base/include/ivisitor.h"
#include "marketplace/include/marketplace.h"
#include "technologies/include/ioutput.h"
#include "functions/include/node_input.h"
#include "functions/include/function_utils.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

InvestConsumer::InvestConsumer(){
}

void InvestConsumer::copyParam( const BaseTechnology* aInvestConsumer,
                                const int aPeriod ) {
	BaseTechnology::copyParam( aInvestConsumer, aPeriod );
    aInvestConsumer->copyParamsInto( *this, aPeriod );
}

void InvestConsumer::copyParamsInto( InvestConsumer& aInvestConsumer,
                                     const int aPeriod ) const {
}

InvestConsumer* InvestConsumer::clone() const {
    return new InvestConsumer( *this );
}

//! parses rest of InvestConsumer xml object
bool InvestConsumer::XMLDerivedClassParse( const string& aNodeName, const DOMNode* aCurr ) {
    return false;
}

//! Write out additional data members to XML output stream.
void InvestConsumer::toInputXMLDerived( ostream& aOut, Tabs* aTabs ) const {
}

//! Output debug info for derived class
void InvestConsumer::toDebugXMLDerived( const int aPeriod, ostream& aOut, Tabs* aTabs ) const {
}

//! Complete the initialization.
void InvestConsumer::completeInit( const string& aRegionName,
                                   const string& aSectorName,
                                   const string& aSubsectorName )
{
	BaseTechnology::completeInit( aRegionName, aSectorName, aSubsectorName );
}

//! initialize anything that won't change during the calculation
void InvestConsumer::initCalc( const MoreSectorInfo* aMoreSectorInfo, const string& aRegionName, 
                               const string& aSectorName, NationalAccount& aNationalAccount, 
                               const Demographic* aDemographics, const double aCapitalStock, const int aPeriod ) 
{
    Consumer::initCalc( aMoreSectorInfo, aRegionName, aSectorName,
                        aNationalAccount, aDemographics, aCapitalStock,
                        aPeriod );

    const Modeltime* modelTime = scenario->getModeltime();
    if ( aPeriod == 0 && year == modelTime->getper_to_yr( aPeriod ) ) {
        Marketplace* marketplace = scenario->getMarketplace();
        // calculate Price Paid
        BaseTechnology::calcPricePaid(aMoreSectorInfo, aRegionName, aSectorName, aPeriod, 
            modelTime->gettimestep( aPeriod ) );
        // the initial capital good price should have been set by now
        mNestedInputRoot->setPricePaid( FunctionUtils::getCapitalGoodPrice( aRegionName, aPeriod ),
            aPeriod );
        
        mNestedInputRoot->calcCoefficient( aRegionName, aSectorName, 0 );
    }
    else if ( year == modelTime->getper_to_yr( aPeriod ) ) {
        mNestedInputRoot->changeSigma( aRegionName, aPeriod, 
                mNestedInputRoot->getCoefficient( aPeriod ) );
    }

    // tech change?
}

//! calculate income
void InvestConsumer::calcIncome( NationalAccount& aNationalAccount, const Demographic* aDemographics, 
                                 const string& aRegionName, int aPeriod ) 
{
    // Whats up with these?-JPL
	//allocateTransportationDemand( nationalAccount, regionName, period );
	//allocateDistributionCost( nationalAccount, regionName, period );
    double investment = aNationalAccount.getAccountValue( NationalAccount::ANNUAL_INVESTMENT );
    expenditures[ aPeriod ].setType( Expenditure::INVESTMENT, investment );
    scenario->getMarketplace()->addToDemand( "Capital", aRegionName, investment, aPeriod );
	// set National Accounts Consumption for GNP calculation
	aNationalAccount.addToAccount( NationalAccount::GNP_NOMINAL, investment );
	aNationalAccount.addToAccount( NationalAccount::INVESTMENT_NOMINAL, investment );
}

//! calculate demand
void InvestConsumer::operate( NationalAccount& aNationalAccount, const Demographic* aDemographics, 
                             const MoreSectorInfo* aMoreSectorInfo, const string& aRegionName, 
                             const string& aSectorName, const bool aIsNewVintageMode, int aPeriod ) 
{
    const Modeltime* modelTime = scenario->getModeltime();
	if( year == modelTime->getper_to_yr( aPeriod ) ) {
        expenditures[ aPeriod ].reset();
		// calculate prices paid for consumer inputs
		BaseTechnology::calcPricePaid(aMoreSectorInfo, aRegionName, aSectorName, aPeriod, 
            modelTime->gettimestep( aPeriod ) );
        mNestedInputRoot->calcLevelizedCost( aRegionName, aSectorName, aPeriod, 
                mNestedInputRoot->getCoefficient( aPeriod ) ); // alpha zero is the root's alpha
        
        // store the price for reporting
        mCapitalGoodPrice.set( mNestedInputRoot->getPricePaid( aRegionName, aPeriod ) );
        
        // we need to convert the annual investment from a value to a quantity by dividing
        // by the price of the capital good and use that to drive demands
        const double capitalQuantity = aNationalAccount.getAccountValue( NationalAccount::ANNUAL_INVESTMENT ) 
            / mCapitalGoodPrice;
        
		// calculate consumption demands for each final good or service
		calcInputDemand( capitalQuantity, aRegionName, aSectorName, aPeriod );
        calcIncome( aNationalAccount, aDemographics, aRegionName, aPeriod );
        calcEmissions( aSectorName, aRegionName, aPeriod );
        Expenditure& currExpenditure = expenditures[ aPeriod ];
        for( vector<IInput*>::const_iterator inputIt = mLeafInputs.begin(); inputIt != mLeafInputs.end(); ++inputIt ) {
            (*inputIt)->calcTaxes( aRegionName, &aNationalAccount, &currExpenditure, aPeriod );
        }

        // TODO: what about emissions taxes?
        double investmentTaxes = expenditures[ aPeriod ].getValue( Expenditure::INDIRECT_TAXES );
        scenario->getMarketplace()->addToDemand( "government-taxes", aRegionName, investmentTaxes, aPeriod );

        // calculate the real amount consumed
        // TODO: this could currently just go in post calc
        aNationalAccount.addToAccount( NationalAccount::INVESTMENT_REAL,
            calcRealGNP( aNationalAccount, aRegionName, aPeriod ) );

        mPricePaidCached = false;
        mNestedInputRoot->resetCalcLevelizedCostFlag();
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
const string& InvestConsumer::getXMLName() const {
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
const string& InvestConsumer::getXMLNameStatic() {
    const static string XML_NAME = "investConsumer";
    return XML_NAME;
}

//! SGM version of outputing data members to a csv file
void InvestConsumer::csvSGMOutputFile( ostream& aFile, const int aPeriod ) const {
	if ( year == scenario->getModeltime()->getper_to_yr( aPeriod ) ) {
		aFile << "***** Investment Sector Results *****" << endl << endl;
        aFile << "Capital good price," << mCapitalGoodPrice << endl;
		expenditures[ aPeriod ].csvSGMOutputFile( aFile, aPeriod );
		aFile << endl;

		aFile << "Investment Consumer Expenditure" << endl << endl;
		BaseTechnology::csvSGMOutputFile( aFile, aPeriod );
	}
}

void InvestConsumer::accept( IVisitor* aVisitor, const int aPeriod ) const {
    aVisitor->startVisitInvestConsumer( this, aPeriod );
    Consumer::accept( aVisitor, aPeriod );
    aVisitor->endVisitInvestConsumer( this, aPeriod );
}
