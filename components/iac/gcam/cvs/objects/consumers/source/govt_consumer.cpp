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
* \file govt_consumer.cpp
* \ingroup Objects
* \brief The GovtConsumer class source file.
*
* \author Sonny Kim
* \author Katherine Chung
*/
#include "util/base/include/definitions.h"
#include <cmath>
#include <xercesc/dom/DOMNode.hpp>

#include "consumers/include/govt_consumer.h"
#include "util/base/include/xml_helper.h"
#include "containers/include/national_account.h"
#include "technologies/include/expenditure.h"
#include "functions/include/iinput.h"
#include "functions/include/ifunction.h"
#include "containers/include/scenario.h"
#include "marketplace/include/marketplace.h"
#include "functions/include/function_manager.h"
#include "demographics/include/demographic.h"
#include "util/base/include/ivisitor.h"
#include "functions/include/function_utils.h"
#include "util/logger/include/ilogger.h"
#include "technologies/include/ioutput.h"
#include "util/base/include/configuration.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

//!< Default Constructor
GovtConsumer::GovtConsumer() {
}

GovtConsumer* GovtConsumer::clone() const {
    return new GovtConsumer( *this );
}

/*! \brief Used to merge an existing consumer with the coefficients from the previous periods
 *
 * \author Sonny Kim
 * \warning Does not copy everything, only non calculated values to pass on to the next period
 */
void GovtConsumer::copyParam( const BaseTechnology* baseTech,
                              const int aPeriod ) {
    BaseTechnology::copyParam( baseTech, aPeriod );
    baseTech->copyParamsInto( *this, aPeriod );
}


/*! \brief Merges consumers, this is a trick to get c++ to let us use the BaseTech pointer as a consumer without casting
 *
 * \author Sonny Kim
 */
void GovtConsumer::copyParamsInto( GovtConsumer& aGovtConsumer,
                                   const int aPeriod ) const {
    aGovtConsumer.mBaseTransferPopCoef.init( mBaseTransferPopCoef );
    aGovtConsumer.mBaseDeficit.init( mBaseDeficit );
    aGovtConsumer.mBaseTransfer.init( mBaseTransfer );
    aGovtConsumer.mTaxProportional.init( mTaxProportional );
    aGovtConsumer.mTaxAdditive.init( mTaxAdditive );
    aGovtConsumer.mRho.init( mRho );
}

//! Parse xml file for data
bool GovtConsumer::XMLDerivedClassParse( const string &nodeName, const DOMNode* curr ) {
    if ( nodeName == "deficit" ) {
        mBaseDeficit = XMLHelper<double>::getValue( curr );
    }
    else if ( nodeName == "rho" ){
        mRho = XMLHelper<double>::getValue( curr );
    }
    // base year transfer to household
    else if ( nodeName == "baseTransfer" ){
        mBaseTransfer = XMLHelper<double>::getValue( curr );
    }
    else {
        return false;
    }
    return true;
}

//! For derived classes to output XML data
void GovtConsumer::toInputXMLDerived( ostream& out, Tabs* tabs ) const {
    XMLWriteElement( mBaseDeficit, "deficit", out, tabs );
}

//! Output debug info for derived class
void GovtConsumer::toDebugXMLDerived( const int period, ostream& out, Tabs* tabs ) const {
    XMLWriteElement( mBaseDeficit, "deficit", out, tabs );
    XMLWriteElement( mTaxProportional, "taxProportional", out, tabs );
    XMLWriteElement( mTaxAdditive, "taxAdditive", out, tabs );
    XMLWriteElement( mBaseTransferPopCoef, "baseTransferPopCoef", out, tabs );
}

//! Complete the initializations.
void GovtConsumer::completeInit( const string& aRegionName,
                                 const string& aSectorName,
                                 const string& aSubsectorName )
{
    //prodDmdFnType = "GovtDemandFn";
    BaseTechnology::completeInit( aRegionName, aSectorName, aSubsectorName );

    assert( mRho != 1 );
    mSigma.init( 1 / ( 1 - mRho ) );

    const Modeltime* modeltime = scenario->getModeltime();

    // Only the base year government consumer should setup the markets.
    if( modeltime->getyr_to_per( year ) == modeltime->getBasePeriod() ){

        // Setup a trial value market for total household taxes so that the
        // ordering of the government and household consumer does not matter.
        // This will always be a regional market.
        Marketplace* marketplace = scenario->getMarketplace();
        const static string GOVT_TAX_MARKET_NAME = "government-taxes";
        if( !marketplace->createMarket( aRegionName, aRegionName, GOVT_TAX_MARKET_NAME,
            IMarketType::TRIAL_VALUE ) )
        {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::ERROR );
            mainLog << "Government taxes trial value market already existed." << endl;
        }

        // Set the market to solve. It may have already been set to solve,
        // but that will not cause an error. Note that they are being solved in
        // the base period.
        if( !Configuration::getInstance()->getBool( "CalibrationActive" ) ){
            for( int period = 0; period < scenario->getModeltime()->getmaxper(); ++period ){
                // we need to set the price to 1 since it defaults to 0 and will only
                // pass forward trail values if the price was 1
                marketplace->setPrice( GOVT_TAX_MARKET_NAME, aRegionName, 1, period );
                marketplace->setMarketToSolve( GOVT_TAX_MARKET_NAME, aRegionName, period );
            }
        }
    }
}

//! initialize anything that won't change during the calculation
void GovtConsumer::initCalc( const MoreSectorInfo* aMoreSectorInfo, const string& aRegionName,
                             const string& aSectorName, NationalAccount& nationalAccount,
                             const Demographic* aDemographics, const double aCapitalStock, const int aPeriod )
{
    Consumer::initCalc( aMoreSectorInfo, aRegionName, aSectorName,
                        nationalAccount, aDemographics, aCapitalStock,
                        aPeriod );

    const Modeltime* modelTime = scenario->getModeltime();
    if ( year == modelTime->getper_to_yr( aPeriod ) ) {
        // calculate Price Paid
        // TODO: is this really needed?
        BaseTechnology::calcPricePaid(aMoreSectorInfo, aRegionName, aSectorName, aPeriod, 
            modelTime->gettimestep( aPeriod ) );
        if( aPeriod == 0 ){

            // set trial price (budget) to the solution price for the base year
            const static string GOVT_TAX_MARKET_NAME = "government-taxes";
            Marketplace* marketplace = scenario->getMarketplace();
            marketplace->setPrice( GOVT_TAX_MARKET_NAME, aRegionName, mBaseTransfer, aPeriod, true );

            calcBaseCoef( nationalAccount, aDemographics );
            // call income calculation once in the base year to calculate consumption
            calcIncome( nationalAccount, aDemographics, aRegionName, aPeriod );
            // TODO: I don't know if I want to worry about govt as it should not have any 
            // inputs anyways
            //prodDmdFn->calcCoefficient(input, expenditures[ aPeriod ].getValue( Expenditure::CONSUMPTION ), aRegionName, aSectorName, aPeriod, mSigma );
        }
        calcTransfer( nationalAccount, aDemographics, aRegionName, aPeriod );
        // Apply technical change to input coefficients.
        //prodDmdFn->applyTechnicalChange( input, TechChange(), aRegionName, aSectorName, aPeriod, 0, 0 );
    }
}

//! calculate subsidy
void GovtConsumer::calcSubsidy( NationalAccount& nationalAccount, const string& regionName, int period ) {
    double subsidy = nationalAccount.getAccountValue( NationalAccount::SUBSIDY );
    expenditures[ period ].setType( Expenditure::SUBSIDY, subsidy );
}

//! calculate deficit
void GovtConsumer::calcDeficit( const string& aRegionName, int aPeriod ) {
    double deficit = mBaseDeficit;
    expenditures[ aPeriod ].setType( Expenditure::SAVINGS, -deficit ); // savings is negative of deficit
    assert( util::isValidNumber( deficit ) );
    
    // add deficit into the demand of capital
    Marketplace* marketplace = scenario->getMarketplace();
    marketplace->addToDemand( "Capital", aRegionName, deficit, aPeriod );
}

/*! \brief Calculate total taxes.
* \param aNationalAccount The container of national accouting values.
* \param aRegionName The name of the region.
* \param aPeriod The period in which to calculate taxes.
*/
void GovtConsumer::calcTotalTax( NationalAccount& aNationalAccount, const string& aRegionName,
                                 const int aPeriod )
{
    // The total taxes comes from the trial government taxes market.  A trial market is 
    // necessary due to ordering issues with other consumers.
    // We need to add the taxes from the production side to the trial market since they
    // do not do that themselves.  Consumers are responsible for adding their taxes into
    // this market.
    Marketplace* marketplace = scenario->getMarketplace();
    mTaxCorporate.set( aNationalAccount.getAccountValue( NationalAccount::CORPORATE_INCOME_TAXES ) );
    mTaxIBT.set( aNationalAccount.getAccountValue( NationalAccount::INDIRECT_BUSINESS_TAX ) );
    double totalProductionTaxes = aNationalAccount.getAccountValue( NationalAccount::CORPORATE_INCOME_TAXES )
        + aNationalAccount.getAccountValue( NationalAccount::INDIRECT_BUSINESS_TAX )
        + aNationalAccount.getAccountValue( NationalAccount::SOCIAL_SECURITY_TAX )
        + aNationalAccount.getAccountValue( NationalAccount::CARBON_TAX );
        // TODO: txpro, txadd ...
    marketplace->addToDemand( "government-taxes", aRegionName, totalProductionTaxes, aPeriod );
    
    // get the trial total taxes
    double totalTaxes = marketplace->getPrice( "government-taxes", aRegionName, aPeriod );
    expenditures[ aPeriod ].setType( Expenditure::INCOME, totalTaxes );

    // is this appropriate to put in here since it did not spend its
    // income on carbon tax but rather it was a part of it
    expenditures[ aPeriod ].setType( Expenditure::CARBON_TAX,
        aNationalAccount.getAccountValue( NationalAccount::CARBON_TAX ) );
}

//! calculate transfer
void GovtConsumer::calcTransfer( NationalAccount& nationalAccount, const Demographic* aDemographics,
                                const string& regionName, int period )
{
    double transfer = mBaseTransferPopCoef * aDemographics->getTotal( period );
    expenditures[ period ].setType( Expenditure::TRANSFERS, transfer );
    nationalAccount.setAccount(NationalAccount::TRANSFERS, transfer );
}

//! Calculate government income.
void GovtConsumer::calcIncome( NationalAccount& nationalAccount, const Demographic* aDemographics,
                              const string& regionName, int period )
{
    calcTotalTax( nationalAccount, regionName, period );
    calcTransfer( nationalAccount, aDemographics, regionName, period );
    calcDeficit( regionName, period );
    calcSubsidy( nationalAccount, regionName, period );

    double consumption = expenditures[ period ].getValue( Expenditure::INCOME )
                         - expenditures[ period ].getValue( Expenditure::SUBSIDY )
                         - expenditures[ period ].getValue( Expenditure::TRANSFERS )
                         - expenditures[ period ].getValue( Expenditure::SAVINGS );


    // move all consumption to the households through the transfer
    double transfer = expenditures[ period ].getValue( Expenditure::TRANSFERS ) + consumption;
    expenditures[ period ].setType( Expenditure::TRANSFERS, transfer );
    nationalAccount.setAccount(NationalAccount::TRANSFERS, transfer );
    
    // consumption is now zero since this will all occur by the households
    // add to accounts just for consistency
    consumption = 0;
    expenditures[ period ].setType( Expenditure::CONSUMPTION, consumption );
    // set National Accounts Consumption for GNP calculation
    nationalAccount.addToAccount( NationalAccount::GNP_NOMINAL, consumption );
    nationalAccount.addToAccount( NationalAccount::GOVERNMENT_NOMINAL, consumption );
}

//! calculate demand
void GovtConsumer::operate( NationalAccount& aNationalAccount, const Demographic* aDemographics,
                           const MoreSectorInfo* aMoreSectorInfo, const string& aRegionName,
                           const string& aSectorName, const bool aIsNewVintageMode, int aPeriod )
{
    const Modeltime* modelTime = scenario->getModeltime();
    if( year == modelTime->getper_to_yr( aPeriod ) ){
        expenditures[ aPeriod ].reset();
        // calculate prices paid for consumer inputs
        BaseTechnology::calcPricePaid( aMoreSectorInfo, aRegionName, aSectorName, aPeriod, 
            modelTime->gettimestep( aPeriod ) );
        // calcIncome is already run in the base period by initCalc
        // use setType in expenditure to override and ensure that marketplace supplies
        // and demands are nulled
        calcIncome( aNationalAccount, aDemographics, aRegionName, aPeriod );
        // calculate consumption demands for each final good or service
        // Government consumers don't shutdown.
        /*const double SHUTDOWN_COEF = 1;
        assert( prodDmdFn );

        double primaryOutput = prodDmdFn->calcDemand( input, 
                                                      expenditures[ aPeriod ].getValue( Expenditure::CONSUMPTION ),
                                                      aRegionName, aSectorName, SHUTDOWN_COEF,
                                                      aPeriod, 0, 0, mSigma, 0 );
        mOutputs[ 0 ]->setCurrencyOutput( aRegionName, primaryOutput, aPeriod );*/

        calcEmissions( aSectorName, aRegionName, aPeriod );
        // calculate the real amount consumed
        // TODO: this could currently just go in post calc
        aNationalAccount.addToAccount( NationalAccount::GOVERNMENT_REAL,
            calcRealGNP( aNationalAccount, aRegionName, aPeriod ) );
        mPricePaidCached = false;
    }
}

//! calculate base coefficient
void GovtConsumer::calcBaseCoef( NationalAccount& nationalAccount, const Demographic* aDemographics ){
    mBaseTransferPopCoef.set( mBaseTransfer / aDemographics->getTotal( 0 ) );
}

//! calculate government capital demand
void GovtConsumer::calcGovtCapitalDemand( const std::string& regionName, int period ){
    assert( false );
    const IInput* capInput = FunctionUtils::getCapitalInput( mLeafInputs );
    assert( capInput );

    double tempCapital = capInput->getCurrencyDemand( period );

    assert( tempCapital >= 0 );
    assert( util::isValidNumber( tempCapital ) );
    Marketplace* marketplace = scenario->getMarketplace();
    // Should this use the demand currency?
    marketplace->addToDemand( "Capital", regionName, tempCapital, period );
    // add capital to ETE, not done
}

//! calculate government tax or subsidy
/*! \todo This isn't called and doesn't work. */
void GovtConsumer::calcGovtTaxOrSubsidy( const string& regionName, int period ){
    assert( false );

    // need to read in transportationCost in the future!
    double transportationCost = 0;
    double eximport = 1;

    double taxGov = 0;
    double subsidyGov = 0;

    Marketplace* marketplace = scenario->getMarketplace();

    for(unsigned int i=0; i<mLeafInputs.size(); i++){

        // P - marketplace price
        double temp1 = ( marketplace->getPrice( mLeafInputs[i]->getName(), regionName, period )
            + transportationCost * eximport )
            * ( marketplace->getDemand( mLeafInputs[i]->getName(), regionName, period )
                - marketplace->getSupply( mLeafInputs[i]->getName(), regionName, period ) )
            * ( mTaxProportional - 1 );

        // This does not seem completely right -JPL
        if( mTaxProportional > 1 ){
            taxGov += temp1;
            //add to nationalaccounts
        }
        else if( mTaxProportional  < 1 ){
            subsidyGov -= temp1;
            // add to nationalaccounts
        }

        double temp2 = ( marketplace->getDemand( mLeafInputs[i]->getName(), regionName, period )
            - marketplace->getSupply( mLeafInputs[i]->getName(), regionName, period ) ) * mTaxAdditive;

        if( mTaxAdditive >= 0 ){
            taxGov += temp2;
        }
        else{
            subsidyGov -= temp2;
        }

    }

    // do something with taxGov and subsidyGov ...
}


//! calculate budget
void GovtConsumer::calcBudget() {
    // TODO: Figure out what is really supposed to happen here.
    double budget = 0; // ??????????????????????????????
    const Modeltime* modeltime = scenario->getModeltime();
    expenditures[ modeltime->getyr_to_per( year ) ].setType( Expenditure::BUDGET, budget );
}

/*! \brief Get the XML node name for output to XML.
*
* This public function accesses the private constant string, XML_NAME.
* This way the tag is always consistent for both read-in and output and can be easily changed.
* This function may be virtual to be overriden by derived class pointers.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME.
*/
const string& GovtConsumer::getXMLName() const {
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
const string& GovtConsumer::getXMLNameStatic() {
    const static string XML_NAME = "govtConsumer";
    return XML_NAME;
}

//! SGM version of outputing data to a csv file
void GovtConsumer::csvSGMOutputFile( ostream& aFile, const int period ) const {
    if ( year == scenario->getModeltime()->getper_to_yr( period ) ) {
        aFile << "***** Government Sector Results *****" << endl << endl;
        aFile << "Tax Accounts" << endl;
        aFile << "Proportional Tax" << ',' << mTaxProportional << endl;
        aFile << "Additive Tax" << ',' << mTaxAdditive << endl;
        aFile << "Corporate Income Tax" << ',' << mTaxCorporate << endl;
        aFile << "Indirect Business Tax" << ',' << mTaxIBT << endl;
        expenditures[ period ].csvSGMOutputFile( aFile, period );
        aFile << endl;
        BaseTechnology::csvSGMOutputFile( aFile, period );
    }
}

void GovtConsumer::accept( IVisitor* aVisitor, const int aPeriod ) const {
    aVisitor->startVisitGovtConsumer( this, aPeriod );
    Consumer::accept( aVisitor, aPeriod );
    aVisitor->endVisitGovtConsumer( this, aPeriod );
}
