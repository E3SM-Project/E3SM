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
* \file trade_input.cpp
* \ingroup Objects
* \brief The TradeInput class source file.
* \author Pralit Patel
*/

#include "util/base/include/definitions.h"
#include "functions/include/trade_input.h"
#include "util/base/include/xml_helper.h"
#include "util/base/include/ivisitor.h"
#include "emissions/include/aghg.h"
#include "sectors/include/more_sector_info.h"
#include "containers/include/national_account.h"
#include "technologies/include/expenditure.h"
#include "containers/include/scenario.h"
#include "marketplace/include/marketplace.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

/*! \brief Get the XML node name in static form for comparison when parsing XML.
* \details This public function accesses the private constant string, XML_NAME.
*          This way
* the tag is always consistent for both read-in and output and can be easily
* changed. The "==" operator that is used when parsing, required this second
* function to return static.
* \note A function cannot be static and virtual.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME as a static.
*/
const string& TradeInput::getXMLNameStatic() {
    const static string XML_NAME = "trade-input";
    return XML_NAME;
}

/*! \brief Get the XML name for reporting to XML file.
*
* This public function accesses the private constant string, XML_NAME. This way
* the tag is always consistent for reporting outputs and can be easily
* changed.
* \author Sonny Kim
* \return The constant XML_NAME.
*/
const string& TradeInput::getXMLReportingName() const{
    /*
    const static string XML_REPORTING_NAME = "input-trade";
    return XML_REPORTING_NAME;
    */
    // TODO this is a hack to be able to identify the trading
    // partner in the XML DB results because there is not start/end
    // visit trade input
    return /*"input-trade-" +*/ mTradingPartner;
}

/*! \brief Get the XML node name for output to XML.
*
* This public function accesses the private constant string, XML_NAME.
* This way the tag is always consistent for both read-in and output and can be easily changed.
* This function may be virtual to be overridden by derived class pointers.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME.
*/
const string& TradeInput::getXMLName() const {
    return getXMLNameStatic();
}

//! Default Constructor
TradeInput::TradeInput() {
}

TradeInput::TradeInput( const TradeInput& aTradeInput ):SGMInput( aTradeInput ),
mTradingPartner( aTradeInput.mTradingPartner ),
mImportTaxRate( aTradeInput.mImportTaxRate ),
mExportTaxRate( aTradeInput.mExportTaxRate )
{
}

TradeInput::~TradeInput() {}

TradeInput* TradeInput::clone() const {
    return new TradeInput( *this );
}

void TradeInput::copyParam( const IInput* aInput,
                                 const int aPeriod ) {
    SGMInput::copyParam( aInput, aPeriod );
    aInput->copyParamsInto( *this, aPeriod );
}

void TradeInput::copyParamsInto( TradeInput& aTradeInput,
                                      const int aPeriod ) const {
    aTradeInput.mTradingPartner = mTradingPartner;

    // only copy these over if they have not already been read
    aTradeInput.mImportTaxRate.init( mImportTaxRate );
    aTradeInput.mExportTaxRate.init( mExportTaxRate );
}

//! XML parsing for derived class
bool TradeInput::XMLDerivedClassParse( const string& aNodeName, const DOMNode* aCurrNode ) {
    if( aNodeName == "trade-partner" ) {
        mTradingPartner = XMLHelper<string>::getValue( aCurrNode );
    }
    else if( aNodeName == "import-tax-rate" ) {
        mImportTaxRate.init( XMLHelper<double>::getValue( aCurrNode ) );
    }
    else if( aNodeName == "export-tax-rate" ) {
        mExportTaxRate.init( XMLHelper<double>::getValue( aCurrNode ) );
    }
    else {
        return false;
    }
    return true;
}

//! Output XML for derived class
void TradeInput::toInputXMLDerived( ostream& aOut, Tabs* aTabs ) const {
    XMLWriteElement( mTradingPartner, "trade-partner", aOut, aTabs );
    XMLWriteElement( mImportTaxRate, "import-tax-rate", aOut, aTabs );
    XMLWriteElement( mExportTaxRate, "export-tax-rate", aOut, aTabs );
}

//! Output debug info to XML for derived class
void TradeInput::toDebugXMLDerived( const int period, ostream& aOut, Tabs* aTabs ) const {
    toInputXMLDerived( aOut, aTabs );
}

void TradeInput::initCalc( const string& aRegionName,
                         const string& aSectorName,
                         const bool aIsNewInvestmentPeriod,
                         const bool aIsTrade,
                         const int aPeriod )
{
    // do the base init calc with the trading partner as the region name
    SGMInput::initCalc( mTradingPartner, aSectorName, aIsNewInvestmentPeriod,
        aIsTrade, aPeriod );

    // add traded to the type flags
    mTypeFlags |= IInput::TRADED;
}

void TradeInput::setPhysicalDemand( double aPhysicalDemand, const string& aRegionName, const int aPeriod )
{
    // do the base set physical demand with the trading partner as the region name
    SGMInput::setPhysicalDemand( aPhysicalDemand, mTradingPartner, aPeriod );
}

double TradeInput::getPrice( const string& aRegionName, const int aPeriod ) const {
    // do the base get price with the trading partner as the region name
    // and adjust for the export tax
    return SGMInput::getPrice( mTradingPartner, aPeriod ) / ( 1 + mExportTaxRate );
}

void TradeInput::calcPricePaid( const string& aRegionName, const string& aSectorName, const MoreSectorInfo* aMoreSectorInfo,
                              const vector<AGHG*>& aGhgs, const ICaptureComponent* aSequestrationDevice,
                              const int aLifetimeYears, const int aPeriod )
{
    // we shouldn't be able to trade factors anyways
    assert( !hasTypeFlag( IInput::FACTOR ) );

    // Calculate GHG taxes.
    double carbonTax = 0;
    for( vector<AGHG*>::const_iterator ghg = aGhgs.begin(); ghg != aGhgs.end(); ++ghg ){
        // TODO: do we want the sector name here?
        carbonTax += (*ghg)->getGHGValue( this, aRegionName, aSectorName,
            aSequestrationDevice, aPeriod );
    }
    double tempPricePaid = ( getPrice( mTradingPartner, aPeriod ) *
        ( 1 + mSalesTaxRate + mImportTaxRate ) +
        carbonTax ) * getPriceAdjustment();
    setPricePaid( tempPricePaid, aPeriod );
}

double TradeInput::calcTaxes( const string& aRegionName, NationalAccount* aNationalAccount, 
                         Expenditure* aExpenditure, const int aPeriod ) const
{
    double price = getPrice( mTradingPartner, aPeriod );
    double currencyDemand = price * getPhysicalDemand( aPeriod );
    double exportTax = mExportTaxRate * currencyDemand * -1;
    double importTax = mImportTaxRate * currencyDemand;
    // TODO: I should create an import tax account to break this out
    if( aNationalAccount && aExpenditure ) {
        aNationalAccount->addToAccount( NationalAccount::INDIRECT_BUSINESS_TAX, importTax );
        aExpenditure->addToType( Expenditure::INDIRECT_TAXES, importTax );

        Marketplace* marketplace = scenario->getMarketplace();
        // the government pays the tax/subsidy for the foreign sector 
        marketplace->addToDemand( "government-taxes", mTradingPartner,
            exportTax, aPeriod );
    }
    return importTax;
}

double TradeInput::getPriceReceived( const string& aRegionName, const int aPeriod ) const {
    // do the base get price received with the trading partner as the region name
    return SGMInput::getPriceReceived( mTradingPartner, aPeriod );
}

void TradeInput::accept( IVisitor* aVisitor, const int aPeriod ) const {
    // TODO: do we need a startVisitTradeInput?
    SGMInput::accept( aVisitor, aPeriod );
}
