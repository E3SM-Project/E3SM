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
* \file production_technology.cpp
* \ingroup Objects
* \brief The ProductionTechnology class source file.
* \author Pralit Patel
* \author Sonny Kim
*/

#include "util/base/include/definitions.h"
#include <iostream>
#include <xercesc/dom/DOMNode.hpp>

#include "technologies/include/production_technology.h"
#include "util/base/include/xml_helper.h"
#include "containers/include/national_account.h"
#include "functions/include/iinput.h"
#include "containers/include/scenario.h"
#include "marketplace/include/marketplace.h"
#include "sectors/include/more_sector_info.h"
#include "util/base/include/ivisitor.h"
#include "investment/include/iexpected_profit_calculator.h"
#include "emissions/include/aghg.h"
#include "technologies/include/technology_type.h"
#include "functions/include/function_utils.h"
#include "technologies/include/ishutdown_decider.h"
#include "technologies/include/profit_shutdown_decider.h"
#include "technologies/include/ioutput.h"
#include "containers/include/iinfo.h"
#include "functions/include/node_input.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

typedef vector<AGHG*>::const_iterator CGHGIterator;
typedef vector<AGHG*>::iterator GHGIterator;

//! Default Constructor
ProductionTechnology::ProductionTechnology() {
    const int maxper = scenario->getModeltime()->getmaxper();
    mCostsReporting.resize( maxper );
    mProfits.resize( maxper );

    // Resize the cached value to the number of types in the CacheValue enum.
    mCachedValues.resize( END );

    mExpectedProfitRateReporting = 0;
    mCapitalStock = 0;
    indBusTax = 0;
    lifeTime = 0;
    delayedInvestTime = 0;
    maxLifeTime = 0;
    retrofitLifeTime = 0;
    periodIniInvest = 0;
    periodInvestUnallowed = 0;
    mAnnualInvestment = 0;
    mBasePhysicalOutput = 0;
    mFixedInvestment = -1;
    mParentTechType = 0;
    mValidCachePeriod = -1;
}

//! Destructor
ProductionTechnology::~ProductionTechnology(){
}

void ProductionTechnology::copyParam( const BaseTechnology* baseTech,
                                      const int aPeriod ) {
    BaseTechnology::copyParam( baseTech, aPeriod );
    baseTech->copyParamsInto( *this, aPeriod );
}

void ProductionTechnology::copyParamsInto( ProductionTechnology& prodTechIn,
                                           const int aPeriod ) const {
     prodTechIn.indBusTax = indBusTax;
     prodTechIn.lifeTime = lifeTime;
     prodTechIn.delayedInvestTime = delayedInvestTime;
     prodTechIn.maxLifeTime = maxLifeTime;
     prodTechIn.retrofitLifeTime = retrofitLifeTime;
     prodTechIn.periodIniInvest = periodIniInvest;
     prodTechIn.periodInvestUnallowed = periodInvestUnallowed;
}

ProductionTechnology* ProductionTechnology::clone() const {
    return new ProductionTechnology( *this );
}

const string& ProductionTechnology::getXMLName() const {
    return getXMLNameStatic();
}

const string& ProductionTechnology::getXMLNameStatic() {
    const static string XML_NAME = "productionTechnology";
    return XML_NAME;
}

//! Parse xml file for data
bool ProductionTechnology::XMLDerivedClassParse( const string& nodeName, const DOMNode* curr ) {
    if ( nodeName == "lifeTime" ) {
        lifeTime = XMLHelper<int>::getValue( curr );
    }
    else if (nodeName == "basePhysicalOutput" ) {
        mBasePhysicalOutput = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "fixed-investment" ){
        mFixedInvestment = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "technicalChangeHicks" ){
        mTechChange.mHicksTechChange = XMLHelper<double>::getValue( curr );
    }
    else if (nodeName == "technicalChangeEnergy" ) {
        mTechChange.mEnergyTechChange = XMLHelper<double>::getValue( curr );
    }
    else if (nodeName == "technicalChangeMaterial" ) {
        mTechChange.mMaterialTechChange = XMLHelper<double>::getValue( curr );
    }

    // I don't think these are utilized so we could possibly get rid of them?
    // also I did not mention them in the XML Specification for this class.
    else if (nodeName == "indirectBusinessTax" ) {
        indBusTax = XMLHelper<double>::getValue( curr );
    }
    else if (nodeName == "delayedInvestTime" ) {
        delayedInvestTime = XMLHelper<int>::getValue( curr );
    }
    else if (nodeName == "maxLifeTime" ) {
        maxLifeTime = XMLHelper<int>::getValue( curr );
    }
    else if (nodeName == "retrofitLifeTime" ) {
        retrofitLifeTime = XMLHelper<int>::getValue( curr );
    }
    else if (nodeName == "periodIniInvest" ) {
        periodIniInvest = XMLHelper<int>::getValue( curr );
    }
    else if (nodeName == "periodInvestUnallowed" ) {
        periodInvestUnallowed = XMLHelper<int>::getValue( curr );
    }
    else if (nodeName == "annual-investment" ) {
        mAnnualInvestment = XMLHelper<double>::getValue( curr );
    }
    else {
        return false;
    }
    return true;
}

//! For derived classes to output XML data
void ProductionTechnology::toInputXMLDerived( ostream& out, Tabs* tabs ) const {
    // Note: We should be using checkDefault version of this function.
    XMLWriteElement( lifeTime, "lifeTime", out, tabs );
    XMLWriteElement( mCapitalStock, "capital-stock", out, tabs );
    XMLWriteElement( indBusTax, "indBusTax", out, tabs );
    XMLWriteElement( mBasePhysicalOutput, "basePhysicalOutput", out, tabs );
    XMLWriteElement( delayedInvestTime, "delayedInvestTime", out, tabs );
    XMLWriteElement( maxLifeTime, "maxLifeTime", out, tabs );
    XMLWriteElement( retrofitLifeTime, "retrofitLifeTime", out, tabs );
    XMLWriteElement( periodIniInvest, "periodIniInvest", out, tabs );
    XMLWriteElement( periodInvestUnallowed, "periodInvestUnallowed", out, tabs );
    XMLWriteElement( mAnnualInvestment, "annual-investment", out, tabs );
    XMLWriteElementCheckDefault( mFixedInvestment, "fixed-investment", out, tabs, -1.0 );
    XMLWriteElement( mTechChange.mEnergyTechChange, "technicalChangeEnergy", out, tabs );
    XMLWriteElement( mTechChange.mMaterialTechChange, "technicalChangeMaterial", out, tabs );
    XMLWriteElement( mTechChange.mHicksTechChange, "technicalChangeHicks", out, tabs );
}

//! Output debug info for derived class
void ProductionTechnology::toDebugXMLDerived( const int period, ostream& out, Tabs* tabs ) const {
    XMLWriteElement( lifeTime, "lifeTime", out, tabs );
    XMLWriteElement( mCapitalStock, "capital-stock", out, tabs );
    XMLWriteElement( indBusTax, "indBusTax", out, tabs );
    XMLWriteElement( mBasePhysicalOutput, "basePhysicalOutput", out, tabs );
    XMLWriteElement( delayedInvestTime, "delayedInvestTime", out, tabs );
    XMLWriteElement( maxLifeTime, "maxLifeTime", out, tabs );
    XMLWriteElement( retrofitLifeTime, "retrofitLifeTime", out, tabs );
    XMLWriteElement( periodIniInvest, "periodIniInvest", out, tabs );
    XMLWriteElement( periodInvestUnallowed, "periodInvestUnallowed", out, tabs );
    XMLWriteElement( mAnnualInvestment, "annual-investment", out, tabs );
    XMLWriteElement( mFixedInvestment, "fixed-investment", out, tabs );
    XMLWriteElement( mTechChange.mEnergyTechChange, "technicalChangeEnergy", out, tabs );
    XMLWriteElement( mTechChange.mMaterialTechChange, "technicalChangeMaterial", out, tabs );
    XMLWriteElement( mTechChange.mHicksTechChange, "technicalChangeHicks", out, tabs );
}

void ProductionTechnology::completeInit( const string& aRegionName,
                                         const string& aSectorName,
                                         const string& aSubsectorName )
{
    BaseTechnology::completeInit( aRegionName, aSectorName, aSubsectorName );
}

void ProductionTechnology::updateMarketplace( const string& sectorName, const string& regionName,
                                              const int period )
{
    BaseTechnology::updateMarketplace( sectorName, regionName, period );
}

void ProductionTechnology::initCalc( const MoreSectorInfo* aMoreSectorInfo, const string& aRegionName,
                                    const string& aSectorName, NationalAccount& nationalAccount,
                                    const Demographic* aDemographics, const double aCapitalStock, const int aPeriod )
{
    const int BASE_PERIOD = 0; // for base period only
    // TODO: this is a hack to avoid copying coefs unecessarily
    if( ( isAvailable( aPeriod ) && !isRetired( aPeriod ) ) || ( aPeriod == BASE_PERIOD && isInitialYear() ) ) {
        mNestedInputRoot->initCalc( aRegionName, aSectorName, isNewInvestment( aPeriod ), isTrade(), aPeriod );
    }
            
    BaseTechnology::initCalc( aMoreSectorInfo, aRegionName, aSectorName,
                              nationalAccount, aDemographics, aCapitalStock,
                              aPeriod );
    mLeafInputs = FunctionUtils::getLeafInputs( mNestedInputRoot );

    // Setup the cached values for the period.
    mValidCachePeriod = aPeriod;
    mCachedValues[ AVAILABLE ] = calcIsAvailable( aPeriod );
    mCachedValues[ RETIRED ] = calcIsRetired( aPeriod );
    mCachedValues[ NEW_INVESTMENT ] = calcIsNewInvestment( aPeriod );

    // calculate coefficients in the base period for techs that are the initial year
    if ( aPeriod == BASE_PERIOD && isInitialYear() ) {
        mNestedInputRoot->initialize();
        // TODO: create a function util to properly calculate life time years
        const int timeStep = scenario->getModeltime()->gettimestep( aPeriod );
        const int lifetimeYears = lifeTime * timeStep;
        const int techPeriod = scenario->getModeltime()->getyr_to_per( year );
        BaseTechnology::calcPricePaid( aMoreSectorInfo, aRegionName, aSectorName, aPeriod, lifetimeYears );

        // calibrate levelized cost by dividing the currency output by the physical output
        // if the physical output is available otherwise use the sector's price received
        double levelizedCost;
        if(  mBasePhysicalOutput != 0 ) {
            levelizedCost = ( mNestedInputRoot->getCurrencyDemand( aPeriod ) ) / 
                ( mBasePhysicalOutput );
        }
        else {
            levelizedCost = FunctionUtils::getPriceReceived( aRegionName, aSectorName, aPeriod );
        }
        // if a levelized cost was read in do not overide it
        if( mNestedInputRoot->getPricePaid( aRegionName, aPeriod ) == 0 ) {
            mNestedInputRoot->setPricePaid( 
                levelizedCost, aPeriod );
        }

        mNestedInputRoot->calcCoefficient( aRegionName, aSectorName, techPeriod );

        // calculate the levelized cost to ensure that the base year node prices
        // are also stored since they are necessary for adjusting ceoefficients
        mNestedInputRoot->calcLevelizedCost( aRegionName, aSectorName, aPeriod,
                                             mNestedInputRoot->getCoefficient( aPeriod ) );
    }

    // Only do full initialization for active vintages.
    else if( aPeriod > BASE_PERIOD && isAvailable( aPeriod ) && !isRetired( aPeriod ) ){
        // Apply technical change to all technologies and all vintages after base period
        if( year == scenario->getModeltime()->getper_to_yr(aPeriod) ) {
            // change coefficients incase the sigma has changed over time
            mNestedInputRoot->changeSigma( aRegionName, aPeriod, 
                mNestedInputRoot->getCoefficient( aPeriod ) );
            // Apply technical change to input coefficients and alpha zero scaler.
            mNestedInputRoot->applyTechnicalChange( aRegionName, aSectorName, aPeriod, mTechChange );

            // these values are not necessary for producers so zero them out
            // so they do not get written out to the database
            // TODO: it would be nice if I didn't have to do this
            for( vector<IInput*>::iterator it = mLeafInputs.begin(); it != mLeafInputs.end(); ++it ) {
                (*it)->setPhysicalDemand( 0, aRegionName, BASE_PERIOD );
                //(*it)->setPricePaid( 0, BASE_PERIOD );
            }
        }
        // Convert all vintages, not just last.
        else if( ( year == scenario->getModeltime()->getper_to_yr(aPeriod - 1) ) || 
            ( year < scenario->getModeltime()->getper_to_yr(0) && (aPeriod == 1) ) ) 
        { 
            // transform base and pre-base period vintages in period 1.
            // this will be using stored base year prices
            mNestedInputRoot->changeElasticity( aRegionName, aPeriod, 
                mNestedInputRoot->getCoefficient( aPeriod ) );
        }
    }
}

/*! \brief Return whether a technology is new investment for the current period.
* \param aPeriod The current period.
* \return Whether the technology is new investment in the period.
* \author Josh Lurz
*/
bool ProductionTechnology::isNewInvestment( const int aPeriod ) const {
    // Check if we have a valid cached result for the period.
    if( aPeriod == mValidCachePeriod ){
        return mCachedValues[ NEW_INVESTMENT ];
    }
    // There is no cached result available for this period, so calculate it
    // dynamically.
    return calcIsNewInvestment( aPeriod );
}

/*! \brief Operate the technology.
* \author Josh Lurz
* \param aNationalAccount Regional national accounts container.
* \param aDemographic Regional demographic information.
* \param aMoreSectorInfo Sector information container.
* \param aRegionName Region name.
* \param aSectorName Sector name.
* \param aIsNewVintageMode Whether to operate only old technologies, or new and old technologies.
* \param aPeriod Period to operate in.
* \note sets the output value for the current period. Also sets input values.
*/

void ProductionTechnology::operate( NationalAccount& aNationalAccount, const Demographic* aDemographic, 
                                   const MoreSectorInfo* aMoreSectorInfo, const string& aRegionName, 
                                   const string& aSectorName, const bool aIsNewVintageMode, const int aPeriod )
{
    // Always operate old techs and operate old techs if they are still .
    if( isAvailable( aPeriod ) && !isRetired( aPeriod ) ){
        assert( !mLeafInputs.empty() );
        if( ( aIsNewVintageMode && isNewInvestment( aPeriod ))
            || (!aIsNewVintageMode && !isNewInvestment( aPeriod )) )
        {
            expenditures[ aPeriod ].reset();
            // calculate prices paid for technology inputs
            // TODO: create a function util to properly calculate life time years
            const Modeltime* modelTime = scenario->getModeltime();
            const int timeStep = modelTime->gettimestep( aPeriod );
            const int lifetimeYears = lifeTime * timeStep;
            BaseTechnology::calcPricePaid( aMoreSectorInfo, aRegionName, aSectorName, aPeriod, lifetimeYears );
            
            double primaryOutput;
            if( aIsNewVintageMode ) {
                // TODO: this logic does not need to be here since newVintagePeriod will
                // be equal to aPeriod however it does save a bit of time with fewer uncessary
                // function calls
                primaryOutput = mOutputs[ 0 ]->getPhysicalOutput( aPeriod );
            }
            else {
                double shutdownCoef = calcShutdownCoef( aRegionName, aSectorName, aPeriod );
                // need to shutdown from the maximum capacity which can be found from the
                // output of the new investment year
                int newVintagePeriod = year <= modelTime->getStartYear() ? modelTime->getBasePeriod() :
                    modelTime->getyr_to_per( year );
                primaryOutput = mOutputs[ 0 ]->getPhysicalOutput( newVintagePeriod ) * shutdownCoef;
            }
            assert( primaryOutput >= 0 );
            
            // now that we have the price paid for leaf inputs calculated we can calculate node
            // prices with calcLevelizedCost
            mNestedInputRoot->calcLevelizedCost( aRegionName, aSectorName, aPeriod, 
                mNestedInputRoot->getCoefficient( aPeriod ) ); // alpha zero is the root's alpha

            // finally we have all prices calculated and we have the level of output invested we can go ahead
            // and calculate input demands
            double demand = mNestedInputRoot->calcInputDemand( aRegionName, aSectorName, aPeriod, 
                primaryOutput, 0, mNestedInputRoot->getCoefficient( aPeriod ) ); // alpha zero is the root's alpha

            if( aIsNewVintageMode ) {
                // Determine what the demand was for capital which is the capital stock required
                // to produce that amount of output.  Then calculate the annual investment from that.
                // TODO: make sure the above comment makes sense
                mCapitalStock = 0;
                // TODO: consider merging this loop with the one just below
                for( unsigned int i = 0; i < mLeafInputs.size(); i++ ) {
                    if( mLeafInputs[ i ]->hasTypeFlag( IInput::CAPITAL ) ) {
                        mCapitalStock += mLeafInputs[ i ]->getPhysicalDemand( aPeriod );
                    }
                }
                mAnnualInvestment = ( mCapitalStock / timeStep ) * FunctionUtils::getCapitalGoodPrice( aRegionName, aPeriod );
                aNationalAccount.addToAccount(NationalAccount::ANNUAL_INVESTMENT, mAnnualInvestment );
            }
            else {
                // Need to reset the currency output in old vintages so that they get added
                // into the marketplace supply again.  New vintages get this done by
                // setInvestment.
                mOutputs[ 0 ]->setPhysicalOutput( primaryOutput, aRegionName, mSequestrationDevice.get(), aPeriod );
            }

            // Add wages and land rents to national account
            for( unsigned int i = 0; i < mLeafInputs.size(); i++ ) {
                mLeafInputs[i]->calcTaxes( aRegionName, &aNationalAccount, 
                    &expenditures[ aPeriod ], aPeriod );
                if( mLeafInputs[i]->hasTypeFlag( IInput::LABOR ) ) {
                    double wages = mLeafInputs[i]->getPhysicalDemand( aPeriod ) 
                        * mLeafInputs[ i ]->getPrice( aRegionName, aPeriod );
                    expenditures[ aPeriod ].addToType( Expenditure::WAGES, wages );
                    aNationalAccount.addToAccount(NationalAccount::LABOR_WAGES, wages );
                }
                else if( mLeafInputs[i]->hasTypeFlag( IInput::LAND ) ) {
                    double landRents = mLeafInputs[i]->getPhysicalDemand( aPeriod )
                                     * mLeafInputs[ i ]->getPrice( aRegionName, aPeriod );
                    expenditures[ aPeriod ].addToType(Expenditure::LAND_RENTS, landRents );
                    aNationalAccount.addToAccount(NationalAccount::LAND_RENTS, landRents );
                }
                else if( !mLeafInputs[ i ]->hasTypeFlag( IInput::CAPITAL ) ){
                    expenditures[ aPeriod ].addToType( Expenditure::INTERMEDIATE_INPUTS,
                                           mLeafInputs[i]->getPhysicalDemand( aPeriod )
                                           * mLeafInputs[ i ]->getPrice( aRegionName, aPeriod ) );
                }
            }

            expenditures[ aPeriod ].addToType( Expenditure::SALES, mOutputs[ 0 ]->getCurrencyOutput( aPeriod ) );

            // calculate emissions, taxes, and subsidies
            calcEmissions( aSectorName, aRegionName, aPeriod );
            calcTaxes( aNationalAccount, aMoreSectorInfo, aRegionName, aSectorName, aPeriod );

            // reset hackish flags so that prices get recalculated next model iteration
            mPricePaidCached = false;
            mNestedInputRoot->resetCalcLevelizedCostFlag();
        }
    } // if operational
}

void ProductionTechnology::calcEmissions( const string& aGoodName, const string& aRegionName, const int aPeriod )
{
    // Loop over GHGs and calculate emissions.
    for( GHGIterator ghg = mGhgs.begin(); ghg != mGhgs.end(); ++ghg ){
        (*ghg)->calcEmission( aRegionName, mLeafInputs, mOutputs, 0, mSequestrationDevice.get(), aPeriod );
    }
}

/*! \brief Set a level of new investment for a given period.
* \param aRegionName Region name.
* \param aAnnualInvestment The level of annual investment at the end of the time period.
* \param aTotalInvestment The level of new investment.
* \param aPeriod The period in which to add investment.
* \return The level of new investment actually set.
* \author Josh Lurz
*/
double ProductionTechnology::setInvestment( const string& aRegionName, const double aAnnualInvestment,
                                            const double aTotalInvestment, const int aPeriod )
{
    /*! \pre Annual investment is greater than zero, total investment is greater than zero,
    and the period is greater than zero. */
    assert( aTotalInvestment >= 0 );
    assert( util::isValidNumber( aAnnualInvestment ) );
    assert( util::isValidNumber( aTotalInvestment ) );
    double amountInvested = aTotalInvestment;
    double annualInvested = aAnnualInvestment;
    if( aAnnualInvestment < 0 ) {
        annualInvested = 0;
        amountInvested += aAnnualInvestment;
    }

    // Check to make sure the technology year is the same as the period in which 
    // we are trying to set investment.
    assert( year == scenario->getModeltime()->getper_to_yr( aPeriod ) );

    mOutputs[ 0 ]->setPhysicalOutput( annualInvested, aRegionName, mSequestrationDevice.get(), aPeriod );
    return amountInvested;
}

void ProductionTechnology::calcTaxes( NationalAccount& aNationalAccount, const MoreSectorInfo* aMoreSectorInfo,
                                     const string& aRegionName, const string aSectorName, const int aPeriod )
{
    double investTaxCreditRate = aMoreSectorInfo->getValue(MoreSectorInfo::INVEST_TAX_CREDIT_RATE);
    double investmentTaxCredit = investTaxCreditRate * mAnnualInvestment;
    double indBusTaxRate = aMoreSectorInfo->getValue(MoreSectorInfo::IND_BUS_TAX_RATE);
    Marketplace* marketplace = scenario->getMarketplace();
    indBusTax = indBusTaxRate * mOutputs[ 0 ]->getPhysicalOutput( aPeriod ) 
        * FunctionUtils::getPriceReceived( aRegionName, aSectorName, aPeriod );

    // The shutdown coef would have already been applied to the output so do not
    // apply it again
    double shutdownCoef = 1;

    double outputEmissTaxAdj = 0;
    for( CGHGIterator ghg = mGhgs.begin(); ghg != mGhgs.end(); ++ghg ){
        outputEmissTaxAdj += (*ghg)->getGHGValue( mOutputs[ 0 ], aRegionName, aSectorName,
            mSequestrationDevice.get(), aPeriod );
    }

    double profits = mNestedInputRoot->getFunction()->calcProfits( mLeafInputs, aRegionName, aSectorName, shutdownCoef, aPeriod,
        mOutputs[ 0 ]->getPhysicalOutput( aPeriod ),
        mNestedInputRoot->getCoefficient( aPeriod ), outputEmissTaxAdj );

    // TODO: this is a hack since import sectors should not have profits
    // although the profits are extremely small is still shows up when calculating
    // derivatives in the household and investment which we would like to avoid
    if( aMoreSectorInfo->getValue(MoreSectorInfo::MAX_CORP_RET_EARNINGS_RATE) == 0 ){
        profits = 0;
    }
    assert( util::isValidNumber( profits ) );

    // corporate income tax rate
    // add other value added to rentals, does not include land and labor
    expenditures[ aPeriod ].addToType( Expenditure::RENTALS, profits );
    double corpIncomeTaxRate = aMoreSectorInfo->getValue(MoreSectorInfo::CORP_INCOME_TAX_RATE);
    double corpIncomeTax = corpIncomeTaxRate * profits - investmentTaxCredit;

    double reRate = aMoreSectorInfo->getValue(MoreSectorInfo::MAX_CORP_RET_EARNINGS_RATE)
        *(1 - exp(aMoreSectorInfo->getValue(MoreSectorInfo::RET_EARNINGS_PARAM)
        *marketplace->getPrice("Capital",aRegionName,aPeriod)));

    double dividends = (profits * (1 - corpIncomeTaxRate) + investmentTaxCredit) * (1 - reRate);
    double retainedEarnings = (profits * (1 - corpIncomeTaxRate) + investmentTaxCredit) * reRate;
    
    // All retained earnings to marketplace capital supply
    //assert( retainedEarnings >= 0 );
    assert( util::isValidNumber( retainedEarnings ) );

    marketplace->addToSupply( "Capital", aRegionName, retainedEarnings, aPeriod );
    // Add all taxes and accounts to expenditure
    expenditures[ aPeriod ].setType(Expenditure::INDIRECT_TAXES, indBusTax);
    expenditures[ aPeriod ].setType(Expenditure::DIRECT_TAXES, corpIncomeTax);
    expenditures[ aPeriod ].setType(Expenditure::DIVIDENDS, dividends);
    expenditures[ aPeriod ].setType(Expenditure::RETAINED_EARNINGS, retainedEarnings);
    // Add all taxes to national account
    aNationalAccount.addToAccount(NationalAccount::RETAINED_EARNINGS, retainedEarnings);
    aNationalAccount.addToAccount(NationalAccount::DIVIDENDS, dividends);
    aNationalAccount.addToAccount(NationalAccount::INDIRECT_BUSINESS_TAX, indBusTax);
    aNationalAccount.addToAccount(NationalAccount::CORPORATE_INCOME_TAXES, corpIncomeTax);
    // Add corporate profits to national account
    aNationalAccount.addToAccount(NationalAccount::CORPORATE_PROFITS, profits);

    // calculate taxes for ghgs
    double ghgTax = 0;
    for( GHGIterator ghg = mGhgs.begin(); ghg != mGhgs.end(); ++ghg ){
        double ghgTaxRate = marketplace->getPrice( (*ghg)->getName(), aRegionName, aPeriod, false );
        // the market would not exist if there was no policy
        if( ghgTaxRate == Marketplace::NO_MARKET_PRICE ){
            ghgTaxRate = 0;
        }
        // Adjust greenhouse gas tax with the proportional tax rate.
        // Retrieve proportional tax rate.
        const IInfo* ghgMarketInfo = marketplace->getMarketInfo( (*ghg)->getName(), aRegionName, aPeriod, false );
        // Note: the key includes the region name.
        const double proportionalTaxRate = 
            ( ghgMarketInfo && ghgMarketInfo->hasValue( "proportional-tax-rate" + aRegionName ) ) 
            ? ghgMarketInfo->getDouble( "proportional-tax-rate" + aRegionName, true )
            : 1.0;
        ghgTax += ghgTaxRate * proportionalTaxRate * (*ghg)->getEmission( aPeriod );
    }

    // set the ghg taxes into the accounting structures
    expenditures[ aPeriod ].setType( Expenditure::CARBON_TAX, ghgTax );
    aNationalAccount.addToAccount( NationalAccount::CARBON_TAX, ghgTax );
}

/*! \brief Returns the capital stock of the technology.
* \return The capital stock.
* \author Sonny Kim
*/
double ProductionTechnology::getCapitalStock() const {
    assert( util::isValidNumber( mCapitalStock ) );
    assert( mCapitalStock >= 0 );
    return mCapitalStock;
}

/*! \brief Get the annual investment for this technology.
* \param aPeriod Period to get the annual investment for. If it is -1 return any
*        annual investment.
* \return The annual investment for this technology. Not applicable to all
*         technologies.
* \author Josh Lurz
*/
double ProductionTechnology::getAnnualInvestment( const int aPeriod ) const {
    // Check if the technology had any investment in the period. The -1 flag
    // specifies that the technology should not check that the period matches
    // the period when this vintage is new.
    /*if( aPeriod == -1 || isNewInvestment( aPeriod ) ){
        return mAnnualInvestment;
    }
    return 0;
    */
    return isNewInvestment( aPeriod ) ? mOutputs[ 0 ]->getPhysicalOutput( aPeriod ) : 0;
}

/*! \brief Returns the expected profit rate of the the technology.
* \param aNationalAccount The regional accounting object.
* \param aRegionName The name of the region containing this subsector.
* \param aSectorName The name of the sector containing this subsector.
* \param aExpProfitRateCalc The calculator of expected profit rates.
* \param aInvestmentLogitExp The investment logit exponential.
* \param aIsShareCalc Whether this expected profit rate is being used to
*        calculate shares. Not great.
* \param aIsDistributing Whether this expected profit rate is being used
*        to distribute investment.
* \param aPeriod The period for which to calculate expected profit.
* \return The expected profit rate.
* \author Josh Lurz
*/
double ProductionTechnology::getExpectedProfitRate( const NationalAccount& aNationalAccount,
                                                    const string& aRegionName,
                                                    const string& aSectorName,
                                                    const IExpectedProfitRateCalculator* aExpProfitRateCalc,
                                                    const double aInvestmentLogitExp,
                                                    const bool aIsShareCalc,
                                                    const bool aIsDistributing,
                                                    const int aPeriod ) const
{
    // If its not operational or it is a fixed investment the expected profit
    // rate is 0.
    if( !isNewInvestment( aPeriod ) || ( mFixedInvestment != -1 && aIsDistributing ) ){
        return 0;
    }

    // Calculate the number of years in the lifetime. This would not be correct
    // with variable time steps. This would cause problems with variable
    // time steps.
    const int timeStep = scenario->getModeltime()->gettimestep( aPeriod );
    const int lifetimeYears = lifeTime * timeStep;

    BaseTechnology::calcPricePaid( 0, aRegionName, aSectorName, aPeriod, lifetimeYears );
    mNestedInputRoot->calcLevelizedCost( aRegionName, aSectorName, aPeriod, 
                mNestedInputRoot->getCoefficient( aPeriod ) ); // alpha zero is the root's alpha
    
    // Create the structure of info for the production function.
    ProductionFunctionInfo prodFunc = { mLeafInputs, mNestedInputRoot->getFunction(), 0,
        mNestedInputRoot->getCoefficient( aPeriod ), mCapitalStock, mNestedInputRoot };

    // Use the expected profit visitor to determine the expected profit rate.
    double finalExpectedProfit = aExpProfitRateCalc->calcTechnologyExpectedProfitRate( prodFunc,
                                                                                       aNationalAccount,
                                                                                       aRegionName,
                                                                                       aSectorName,
                                                                                       delayedInvestTime,
                                                                                       lifetimeYears,
                                                                                       timeStep,
                                                                                       aPeriod );
    // TODO: this should go somewhere else
    double outputEmissTaxAdj = 0;
    for( CGHGIterator ghg = mGhgs.begin(); ghg != mGhgs.end(); ++ghg ){
        outputEmissTaxAdj += (*ghg)->getGHGValue( mOutputs[ 0 ], aRegionName, aSectorName,
            mSequestrationDevice.get(), aPeriod );
    }
    finalExpectedProfit -= outputEmissTaxAdj;
    assert( finalExpectedProfit >= 0 );
    return finalExpectedProfit;
}

/*! \brief Get the amount of capital required to produce one unit of output.
* \param aDistributor Investment distributor.
* \param aExpProfitRateCalc Expected profit rate calculator.
* \param aNationalAccount Regional national accounts container.
* \param aRegionName The name of the region.
* \param aSectorName The name of the sector.
* \param aPeriod The period.
* \return The capital output ratio.
* \author Josh Lurz 
*/
double ProductionTechnology::getCapitalOutputRatio( const IDistributor* aDistributor,
                                                    const IExpectedProfitRateCalculator* aExpProfitRateCalc,
                                                    const NationalAccount& aNationalAccount,
                                                    const string& aRegionName,
                                                    const string& aSectorName, 
                                                    const int aPeriod ) const
{
    // Only new investment has a capital to output ratio.
    if( isNewInvestment( aPeriod ) ){
        // Calculate the number of years in the lifetime. This would not be
        // correct with variable time steps.
        return mNestedInputRoot->calcCapitalOutputRatio( aRegionName, aSectorName, aPeriod,
            mNestedInputRoot->getCoefficient( aPeriod ) ); // alpha zero is the root's alpha
    }
    return 0;
}

/*! \brief Get the quantity of fixed investment for the production technology.
* \author Josh Lurz
* \param aPeriod The period for which to get fixed investment.
* \return Fixed investment amount for the vintage.
*/
double ProductionTechnology::getFixedInvestment( const int aPeriod ) const {
    // The -1 flag specified the value wasn't read in, return zero in that case.
    // Use a value class here.
    /*if( isNewInvestment( aPeriod ) && mFixedInvestment != -1 ){
        return mFixedInvestment;
    }
    return 0;
    */
    return ( isNewInvestment/*isAvailable*/( aPeriod ) && mFixedInvestment != -1 ) ? mFixedInvestment : 0;
}

/*! \brief Distribute investment to the vintage.
* \param aDistributor The distributor of new investment.
* \param aNationalAccount The national accounts.
* \param aExpProfitRateCalc The calculator of expected profit rates.
* \param aRegionName The name of the containing region.
* \param aSectorName The name of the containing sector.
* \param aNewInvestment The amount of new investment.
* \param aPeriod The period in which to add investment.
* \return The total amount of investment distributed.
* \author Josh Lurz
*/
double ProductionTechnology::distributeInvestment( const IDistributor* aDistributor,
                                                   NationalAccount& aNationalAccount,
                                                   const IExpectedProfitRateCalculator* aExpProfitRateCalc,
                                                   const string& aRegionName,
                                                   const string& aSectorName,
                                                   const double aNewInvestment,
                                                   const int aPeriod )
{
    assert( mParentTechType );  
    if( !isNewInvestment( aPeriod ) ) {
        return 0;
    }
    // Use the parent technology type to distribute investment. This will
    // eventually call back into this object. Check if investment is fixed.
    if( mFixedInvestment != -1 ){
        if( !util::isEqual( aNewInvestment, 0.0 ) ){
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Vintage " << name
                 << " received investment to distribute but it has fixed investment. " << endl;
        }
        // for fixed set it directly into the output
        mOutputs[ 0 ]->setPhysicalOutput( mFixedInvestment, aRegionName, mSequestrationDevice.get(), aPeriod );
    }
    else if( aNewInvestment > 0 ){
        assert( scenario->getModeltime()->getyr_to_per( year ) == aPeriod );
        return mParentTechType->setTotalInvestment( aRegionName,
                                                    year - scenario->getModeltime()->gettimestep( aPeriod ),
                                                    year, aNewInvestment, aPeriod );
    }
    // Return that no investment was done.
    return 0;
}

/*! \brief Function to finalize objects after a period is solved.
* \details This function is used to calculate and store variables which are only
*          needed after the current period is complete.
* \param aRegionName Region name.
* \param aSectorName Sector name.
* \param aPeriod The period to finalize.
* \todo Finish this function, could move transform here.
* \author Josh Lurz
*/
void ProductionTechnology::postCalc( const string& aRegionName, const string& aSectorName, 
                                     const int aPeriod )
{
    if( !isRetired( aPeriod ) && year <= scenario->getModeltime()->getper_to_yr( aPeriod ) ){
        // Save cost, profit rate, and output. Profit rate is needed to convert
        // elasticities, not just for reporting. Calculate the number of years
        // in the lifetime. This would not be correct with variable time steps.
        double shutdownCoef = 1; //calcShutdownCoef( aRegionName, aSectorName, aPeriod );
        mProfits[ aPeriod ] = mNestedInputRoot->getFunction()->calcProfits( mLeafInputs, aRegionName, aSectorName, shutdownCoef, aPeriod,
            mOutputs[ 0 ]->getPhysicalOutput( aPeriod ),
            mNestedInputRoot->getCoefficient( aPeriod ) );
        mCostsReporting[ aPeriod ] = mNestedInputRoot->getFunction()->calcCosts( mLeafInputs, aRegionName, 
                mNestedInputRoot->getCoefficient( aPeriod ), aPeriod );

        // If this is new investment store the expected price received and
        // expected profit rate.
        if( isNewInvestment( aPeriod ) ){
            // expected profit rate is really the levelized cost
            // TODO: the levelized cost is what we want to report right? not the expected profit rate right?
            mExpectedProfitRateReporting = mNestedInputRoot->getLevelizedCost( aRegionName, aSectorName, aPeriod );
        }

        // Account for exports and import.  We do this in postCalc because it is inconsequential to
        // operation however the national account is not available here.  As a temporary solution we added it
        // into the market info for Capital.
        Marketplace* marketplace = scenario->getMarketplace();
        IInfo* capitalMarketInfo = marketplace->getMarketInfo( "Capital", aRegionName, aPeriod, true );

        // nominal is the current year quantity * the current year prices, real is the current year
        // quantity * base year prices
        double currExportsNominal = capitalMarketInfo->getDouble( "export-nominal", false );
        double currExportsReal = capitalMarketInfo->getDouble( "export-real", false );
        double currImportsNominal = capitalMarketInfo->getDouble( "import-nominal", false );
        double currImportsReal = capitalMarketInfo->getDouble( "import-real", false );

        // calculate resource depletion, and imports/exports
        for( vector<IInput*>::const_iterator it = mLeafInputs.begin(); it != mLeafInputs.end(); ++it ) {
            if( (*it)->hasTypeFlag( IInput::LAND ) ) {
                IInfo* resourceMarketInfo = marketplace->getMarketInfo( (*it)->getName(), aRegionName, aPeriod, false );
                
                // the resource does not necessarily suppurt depletion
                if( resourceMarketInfo ) {
                    double depletionRate = resourceMarketInfo->getDouble( "depletion-rate", true );
                    double depletion = resourceMarketInfo->getDouble( "depleted-resource", true );
                    depletion += mOutputs[ 0 ]->getPhysicalOutput( aPeriod ) * depletionRate;
                    resourceMarketInfo->setDouble( "depleted-resource", depletion );
                }
            }

            if( (*it)->hasTypeFlag( IInput::TRADED ) ) {
                // add to imports, note we need to use get price here
                // because the input will get the price of the appropriate
                // region
                currImportsNominal += (*it)->getPhysicalDemand( aPeriod ) *
                    (*it)->getPrice( aRegionName, aPeriod );
                currImportsReal += (*it)->getPhysicalDemand( aPeriod ) *
                    (*it)->getPrice( aRegionName, 0 );
            } else if( !(*it)->hasTypeFlag( IInput::FACTOR ) ) {
                // subtract from exports because it was consumed domestically
                currExportsNominal -= (*it)->getPhysicalDemand( aPeriod ) *
                    marketplace->getPrice( (*it)->getName(), aRegionName, aPeriod );
                currExportsReal -= (*it)->getPhysicalDemand( aPeriod ) *
                    marketplace->getPrice( (*it)->getName(), aRegionName, 0 );
            }
        }

        // TODO: this is a hack to exclude TPT-Margin, think of a better way
        if( name != "TPT-Margin" ) {

            // We do not account for exports directly so we need to take the difference between domestic
            // production and domestic consumption of domestic goods.  Note that import sectors are
            // on the production side so they will also be added to the export account however all of the
            // consumption of the import sector is done domestically so they will net zero which is fine
            currExportsNominal += mOutputs[ 0 ]->getPhysicalOutput( aPeriod ) * marketplace->getPrice(
                mOutputs[ 0 ]->getName(), aRegionName, aPeriod );
            currExportsReal += mOutputs[ 0 ]->getPhysicalOutput( aPeriod ) * marketplace->getPrice(
                mOutputs[ 0 ]->getName(), aRegionName, 0 );

            capitalMarketInfo->setDouble( "export-nominal", currExportsNominal );
            capitalMarketInfo->setDouble( "export-real", currExportsReal );
            capitalMarketInfo->setDouble( "import-nominal", currImportsNominal );
            capitalMarketInfo->setDouble( "import-real", currImportsReal);
        }
    }
    // this is sort of a hack but we need to clear technologies with years before
    // the initial model year that didn't really operate so that we don't report
    // data for them
    if( aPeriod == 0 && isRetired( aPeriod ) && year <= scenario->getModeltime()->getper_to_yr( aPeriod ) ){
        for( vector<IInput*>::iterator it = mLeafInputs.begin(); it != mLeafInputs.end(); ++it ) {
            (*it)->setPhysicalDemand( 0, aRegionName, aPeriod );
        }
    }
    // we should clear the list of inputs which are really leaves because they
    // should be managed between here and the next initCalc by the nesting structure
    mLeafInputs.clear();
}

void ProductionTechnology::csvSGMOutputFile( ostream& aFile, const int period ) const {
    // print for all operating vintages
    if( isAvailable( period ) && !isRetired( period ) ){
        // BaseTechnology::csvSGMOutputFile( aFile, period );
    }
}

void ProductionTechnology::accept( IVisitor* aVisitor, const int aPeriod ) const
{
    aVisitor->startVisitProductionTechnology( this, aPeriod );
    BaseTechnology::accept( aVisitor, aPeriod );
    aVisitor->endVisitProductionTechnology( this, aPeriod );
}

//! Set the parent technology type helper object. This may change.
void ProductionTechnology::setTypeHelper( TechnologyType* aTechType ){
    assert( aTechType );
    mParentTechType = aTechType;
}

/*! \brief Calculate the coefficient used to shutdown older unprofitable
*          vintages.
* \details This coefficient is always one for new vintages. MORE HERE.
* \param aRegionName Name of the region containing the vintage.
* \param aSectorName Name of the sector containing the vintage.
* \param aPeriod Period
* \return Coefficient which scales down unprofitable vintages.
* \author Josh Lurz
*/
double ProductionTechnology::calcShutdownCoef( const string& aRegionName,
                                               const string& aSectorName,
                                               const int aPeriod ) const
{
    // Never shutdown new investment. This would be inconsistent, and avoiding
    // the calculation should be faster. This should not occur in a solved
    // iteration.
    if( isNewInvestment( aPeriod ) ){
        return 1;
    }
    // need to calculate the variable levelized cost before we try to get the shutdown
    // coef
    mNestedInputRoot->calcVariableLevelizedCost( aRegionName, aSectorName, aPeriod, 
        mNestedInputRoot->getCoefficient( aPeriod ) );

    // TODO: this really needs to be part of the nested input root
    double outputEmissTaxAdj = 0;
    for( CGHGIterator ghg = mGhgs.begin(); ghg != mGhgs.end(); ++ghg ){
        outputEmissTaxAdj += (*ghg)->getGHGValue( mOutputs[ 0 ], aRegionName, aSectorName,
            mSequestrationDevice.get(), aPeriod );
    }
    mNestedInputRoot->setPricePaid( mNestedInputRoot->getPricePaid( aRegionName, aPeriod )
        - outputEmissTaxAdj, aPeriod );
    // Could optimize by storing the shutdown coef.
    // Create the structure of info for the production function.
    // TODO: figure out what I want to do about this
    ProductionFunctionInfo prodFunc = { mLeafInputs, mNestedInputRoot->getFunction(), 0,
        mNestedInputRoot->getCoefficient( aPeriod ), mCapitalStock, mNestedInputRoot };

    // Create a new shutdown decider.
    // TODO: this should be a member of the class so that the user can override the
    // paramaters from XML
    ProfitShutdownDecider shutdownDecider;
    // TODO: these parameters are arbitrary
    // allow 100% shutdown
    shutdownDecider.mMaxShutdown = 1;
    // 50% shutdown when variable costs are above price recieved by 20%
    shutdownDecider.mMedianShutdownPoint = -0.2;
    // tuning parameter for how quickly we shutdown
    shutdownDecider.mSteepness = 10;
    double shutdownCoef = shutdownDecider.calcShutdownCoef( &prodFunc,
                                                            IShutdownDecider::getUncalculatedProfitRateConstant(), 
                                                            aRegionName,
                                                            aSectorName,
                                                            year,
                                                            aPeriod );
    return shutdownCoef;
}

/*! \brief Calculate dynamically whether a technology is new investment for the
*          current period.
* \param aPeriod The current period.
* \return Whether the technology is new investment in the period.
* \warning This function is slower than the cached version and should only be
*          used to setup the cache or by the optimized version of the function
*          when the cached value is not available for a period.
* \sa isNewInvestment
* \author Josh Lurz
*/
bool ProductionTechnology::calcIsNewInvestment( const int aPeriod ) const {
    // Return whether the technology's first year is the same as the period.
    return( year == scenario->getModeltime()->getper_to_yr( aPeriod ) );
}

/*! \brief Calculate dynamically whether a technology has been retired yet.
* \param aPeriod The current period.
* \return Whether the technology has been retired.
* \warning This function is slower than the cached version and should only be
*          used to setup the cache or by the optimized version of the function
*          when the cached value is not available for a period.
* \sa isRetired
* \author Josh Lurz
*/
bool ProductionTechnology::calcIsRetired( const int aPeriod ) const {
    const Modeltime* modeltime = scenario->getModeltime();
    // Return whether the period is past the lifetime of the technology as
    // calculated from when it came on-line.
    return( modeltime->getper_to_yr( aPeriod ) >= ( year + delayedInvestTime 
                                                  + lifeTime * modeltime->gettimestep( aPeriod ) ) );
}

/*! \brief Calculate dynamically whether a technology is available to go on-line.
* \param aPeriod The current period.
* \return Whether the technology has gone on-line.
* \warning This function is slower than the cached version and should only be
*          used to setup the cache or by the optimized version of the function
*          when the cached value is not available for a period.
* \sa isAvailable
* \author Josh Lurz
*/
bool ProductionTechnology::calcIsAvailable( const int aPeriod ) const {
    // Return whether we are past the year when the technology came on-line.
    return ( scenario->getModeltime()->getper_to_yr( aPeriod ) >= year + delayedInvestTime );
}
