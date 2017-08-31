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
* \file production_sector.cpp
* \ingroup Objects
* \brief Sector class source file.
* \author Sonny Kim
*/

#include "util/base/include/definitions.h"
#include <string>
#include <cassert>
#include <xercesc/dom/DOMNode.hpp>
#include "util/base/include/xml_helper.h"
#include "sectors/include/subsector.h"
#include "sectors/include/production_sector.h"
#include "sectors/include/more_sector_info.h"
#include "containers/include/scenario.h"
#include "investment/include/iinvestor.h"
// Need a factory method.
#include "investment/include/accelerator.h"
#include "investment/include/market_based_investment.h"
#include "investment/include/investment_utils.h"
#include "marketplace/include/marketplace.h"
#include "marketplace/include/imarket_type.h"
#include "util/base/include/ivisitor.h"
#include "util/logger/include/ilogger.h"
#include "containers/include/iinfo.h"
#include "functions/include/function_utils.h"
#include "util/base/include/configuration.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

/*! \brief Constructor
* \details Initializes the Sector and initializes all characteristics flags to false.
* \param aRegionName Name of the region containing this sector.
*/
ProductionSector::ProductionSector ( const string& aRegionName ) : Sector ( aRegionName ) {
    mFixedPrices.resize( scenario->getModeltime()->getmaxper() );
    mIsFixedPrice = false;
    mIsEnergyGood = false;
    mIsPrimaryEnergyGood = false;
    mIsSecondaryEnergyGood = false;
}

/*! \brief Default destructor
* \note This is necessary because of the auto_ptr to the IInvestor object.
*/
ProductionSector::~ProductionSector() {
}

/*! \brief Parses any child nodes specific to the derived class.
* \details Method parses any input data from this child class that is specific
*          to this class.
* \author Sonny Kim, Josh Lurz
* \param nodeName name of current node
* \param curr pointer to the current node in the XML input tree
*/
bool ProductionSector::XMLDerivedClassParse( const string& nodeName, const DOMNode* curr ) {
    // Note: This doesn't handle add-ons here. Difficult part is if a different type of 
    // sector investment is requested in an add on.
    if( nodeName == Accelerator::getXMLNameStatic() ){
        mInvestor.reset( new Accelerator() );
        mInvestor->XMLParse( curr );
    }
    else if( nodeName == MarketBasedInvestor::getXMLNameStatic() ){
        mInvestor.reset( new MarketBasedInvestor() );
        mInvestor->XMLParse( curr );
    }
    else if( nodeName == "market-name" ){
        mMarketName = XMLHelper<string>::getValue( curr );
    }
    // Note: This behavior is either on or off, not by period currently.
    else if( nodeName == "FixedPricePath" ){
        mIsFixedPrice = XMLHelper<bool>::getValue( curr );
    }
    else if( nodeName == "numeraire" ){
        mIsNumeraireSector = XMLHelper<bool>::getValue( curr );
    }
    else if( nodeName == "ghgEmissCoef" ){
        ghgEmissCoefMap[ XMLHelper<string>::getAttr( curr, "name" ) ] = XMLHelper<double>::getValue( curr );
    } 
    else if( nodeName == "IsEnergyGood" ){
        mIsEnergyGood = XMLHelper<bool>::getValue( curr );
    }
    else if( nodeName == "IsPrimaryEnergyGood" ){
        mIsPrimaryEnergyGood = XMLHelper<bool>::getValue( curr );
    }
    else if( nodeName == "IsSecondaryEnergyGood" ){
        mIsSecondaryEnergyGood = XMLHelper<bool>::getValue( curr );
    }
    else if( nodeName == "sectorprice" ){
        XMLHelper<double>::insertValueIntoVector( curr, mFixedPrices, scenario->getModeltime() );
    }
    else {
        return false;
    }
    return true;
}

/* \brief Write out ProductionSector specific data to the input XML file.
* \param out Stream into which to write.
* \tabs Object responsible for tabs in the output.
*/
void ProductionSector::toInputXMLDerived( std::ostream& out, Tabs* tabs ) const {
    if( mInvestor.get() ){
        mInvestor->toInputXML( out, tabs );
    }
    // write out the market string.
    XMLWriteElement( mMarketName, "market-name", out, tabs );
    XMLWriteElementCheckDefault( mIsFixedPrice, "FixedPricePath", out, tabs );
    XMLWriteElementCheckDefault( mIsEnergyGood, "IsEnergyGood", out, tabs );
    XMLWriteVector( mFixedPrices, "sectorprice", out, tabs, scenario->getModeltime(), 0.0 );
    for( map<string, double>::const_iterator coef = ghgEmissCoefMap.begin(); coef != ghgEmissCoefMap.end(); ++coef ){
        XMLWriteElement( coef->second, "ghgEmissCoef", out, tabs, 0, coef->first );
    }
}

/* \brief Write out ProductionSector specific data to the debugging XML file.
* \param period Period for which to write information.
* \param out Stream into which to write.
* \tabs Object responsible for tabs in the output.
*/
void ProductionSector::toDebugXMLDerived( const int period, std::ostream& out, Tabs* tabs ) const {
    if( mInvestor.get() ){
        mInvestor->toDebugXML( period, out, tabs );
    }
    
    // write out the market string.
    XMLWriteElement( mMarketName, "market-name", out, tabs );
    XMLWriteElement( mIsFixedPrice, "FixedPricePath", out, tabs, false );
    XMLWriteElement( mIsEnergyGood, "IsEnergyGood", out, tabs );
    XMLWriteElement( mFixedPrices[ period ], "fixed-price", out, tabs );
    for( map<string, double>::const_iterator coef = ghgEmissCoefMap.begin(); coef != ghgEmissCoefMap.end(); ++coef ){
        XMLWriteElement( coef->second, "ghgEmissCoef", out, tabs, 0, coef->first );
    }
}

/*! \brief Complete the initialization of the ProductionSector.
* \details Completes the initialization of the parent Sector and initializes the
*          IInvestor object. The sector investor is initialized to an
*          Accelerator if a specific investor was not read in. The sector
*          investor is then initialized.
* \param aRegionInfo Regional information object.
* \param aDependencyFinder Regional dependency finder.
* \param aLandAllocator Regional land allocator.
*/
void ProductionSector::completeInit( const IInfo* aRegionInfo,
                                     DependencyFinder* aDependencyFinder,
                                     ILandAllocator* aLandAllocator )
{
    // Set the market.
    setMarket();
    // Call parent class complete init.
    Sector::completeInit( aRegionInfo, aDependencyFinder, aLandAllocator );    
    // Initialize the investment object to the default if one has not been read in.
    if( !mInvestor.get() ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::DEBUG );
        mainLog << "Creating default investment type for sector " << name << "." << endl;
        mInvestor.reset( new Accelerator() );
    }

    mInvestor->completeInit( regionName, name );
}

/*! \brief Initialize the market required by this ProductionSector.
* \details Creates the market for the Production sector and initializes it. The
*          market is initialized by setting the optionally read-in prices for
*          each period into the marketplace, setting the market to solve if the
*          ProductionSector is not a fixed-price sector, and setting the flag in
*          the IInfo which tells whether the market is a fixed price
*          market. Fixed price markets are sectors which do not have solved
*          prices, the equilibrium price is always the initial price. This is
*          for specifying price paths exogenously, usually for resource sectors.
*          The read-in prices will be used as the initial prices for the Market,
*          for fixed price path sectors this will also be the final price.
*/
void ProductionSector::setMarket() {
    // Check if the market name was not read-in
    if( mMarketName.empty()  ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::DEBUG );
        mainLog << "Market name for production sector was not read-in. Defaulting to region name." << endl;
        mMarketName = regionName;
    }
    
    // Create the market.
    Marketplace* marketplace = scenario->getMarketplace();
    const int START_PERIOD = scenario->getModeltime()->getBasePeriod();

    // name is Sector name (name of good supplied or demanded)
    // market is the name of the regional market from the input file (i.e., global, region, regional group, etc.)
    if( marketplace->createMarket( regionName, mMarketName, name, IMarketType::NORMAL ) ) {
        // Set the base year price which the sector reads in, into the mFixedPrices vector.
        // TODO: Separate MiniCAM sector so this is not needed.

        mFixedPrices[ 0 ] = mBasePrice;
        marketplace->setPriceVector( name, regionName, mFixedPrices );
        for( int period = 0; period < scenario->getModeltime()->getmaxper(); ++period ){
            // Setup the market information on the sector.
            IInfo* marketInfo = marketplace->getMarketInfo( name, regionName, period, true );
            if( !mIsFixedPrice ){
                if( !Configuration::getInstance()->getBool( "CalibrationActive" ) ){
                    if( period >= START_PERIOD ){
                        marketplace->setMarketToSolve( name, regionName, period );
                    }
                }
            }
            else {
                marketInfo->setBoolean( "IsFixedPrice", true );
            }
            
            // Set whether it is an energy or material good. 
            marketInfo->setBoolean( "IsEnergyGood", mIsEnergyGood );
            marketInfo->setBoolean( "IsPrimaryEnergyGood", mIsPrimaryEnergyGood );
            marketInfo->setBoolean( "IsSecondaryEnergyGood", mIsSecondaryEnergyGood );
            
            // Set the energy to physical conversion factor if the MoreSectorInfo
            // specifies one.
            if( moreSectorInfo.get() ){
                double newConversionFactor = moreSectorInfo->getValue( MoreSectorInfo::ENERGY_CURRENCY_CONVERSION );
                marketInfo->setDouble( "ConversionFactor", newConversionFactor );
            }

            // add ghg gass coefficients to the market info for this sector
            for( map<string,double>::iterator i = ghgEmissCoefMap.begin(); i != ghgEmissCoefMap.end(); ++i ){
                marketInfo->setDouble( i->first + "coefficient", i->second );
            }
        }
    }
}

/*! \brief Initialize the ProductionSector before a period is begun.
* \details TODO

* \param aNationalAccount National accounts container.
* \param aDemographics The demographics object.
* \param aPeriod Period for which to initialize the ProductionSector.
*/
void ProductionSector::initCalc( NationalAccount* aNationalAccount,
                                 const Demographic* aDemographics,
                                 const int aPeriod )
{
    // retained earnings parameter calculation
    // for base year only, better in completeInit but NationalAccount not accessible.
    if ( aPeriod == 0 && moreSectorInfo.get() && aNationalAccount ) {
        Marketplace* marketplace = scenario->getMarketplace();
        double corpIncTaxRate = aNationalAccount->getAccountValue(NationalAccount::CORPORATE_INCOME_TAX_RATE);

        // get the total retained earnings and total profits which were read in from the base dataset
        // so that we back out a retained earnings param we can use for this sector that will reproduce
        // base year data
        double totalRetEarnings = aNationalAccount->getAccountValue(NationalAccount::RETAINED_EARNINGS);
        double totalProfits = aNationalAccount->getAccountValue(NationalAccount::CORPORATE_PROFITS);

        // Set retained earnings to zero by setting MAX_CORP_RET_EARNINGS_RATE
        // to 0. This is for the technology sectors, like the transportation
        // vehicle sectors. All of the profits goes to dividends and there are
        // no retained earnings.
        // SHK  4/21/2005
        double retEarnParam = 0;
        if( moreSectorInfo->getValue(MoreSectorInfo::MAX_CORP_RET_EARNINGS_RATE) != 0 ) {
            // back out retained earnings param
            retEarnParam = log( 1 - (totalRetEarnings/(totalProfits*moreSectorInfo->getValue(MoreSectorInfo::MAX_CORP_RET_EARNINGS_RATE)
                *(1 - corpIncTaxRate)))) / marketplace->getPrice("Capital", regionName, aPeriod );
        }

        // set these params into the sector info so that they may be retrieved by the technologies
        // when calculating taxes and retained earnings
        moreSectorInfo->setType(MoreSectorInfo::RET_EARNINGS_PARAM, retEarnParam);
        moreSectorInfo->setType(MoreSectorInfo::CORP_INCOME_TAX_RATE, corpIncTaxRate);
    }

    // we should calculate price recived now since technologies may need that info when
    // calibrating it's nested input structure
    calcPriceReceived( aPeriod );
    
    // The ITC is being read in at the sector level but set to the national level?
    // TODO: check this, at the moment ITC is not used
    if( moreSectorInfo.get() ){
        aNationalAccount->setAccount( NationalAccount::INVESTMENT_TAX_CREDIT, 
            moreSectorInfo->getValue( MoreSectorInfo::INVEST_TAX_CREDIT_RATE ) );
    }

    Sector::initCalc( aNationalAccount, aDemographics, aPeriod );

    //*************** begin initialize the investment routine ********
    // Calculate and distribute investment to the subsectors of this sector.
    vector<IInvestable*> investableSubsecs = InvestmentUtils::convertToInvestables( subsec );
    // aNationalAccount is a pointer here, but initCalc needs a reference to aNationalAccount
    mInvestor->initCalc( investableSubsecs, *aNationalAccount, aDemographics, aPeriod );
    //*************** end initialize the investment routine ********

}

/*! \brief Returns the output of the ProductionSector.
* \details Dynamically calculates and returns the sum of the output of all
*          subsectors of this production sector. A check for the validity of
*          each Subsector output level is performed before summing the value.
* \author Sonny Kim
* \param aPeriod Model period
* \return Total sector level output.
*/
double ProductionSector::getOutput( const int aPeriod ) const {
    double output = 0;
    for ( unsigned int i = 0; i < subsec.size(); ++i ) {
        double subsecOutput = subsec[ i ]->getOutput( aPeriod );
        // error check.
        if ( !util::isValidNumber( subsecOutput ) ){
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::ERROR );
            mainLog << "Output for subsector " << subsec[ i ]->getName() << " in Sector " << name 
                    << " in region " << regionName <<" is not valid." << endl;
            continue;
        }
        output += subsecOutput;
    }
    return output;
}

/*!
 * \brief Return the price of a production sector.
 * \param aGDP Regional GDP container(null for SGM).
 * \details Production sectors always have markets, so the price of a production
 *          sector is the market price.
 * \param aPeriod Model period.
 * \return Price.
 */
double ProductionSector::getPrice( const GDP* aGDP,
                                   const int aPeriod ) const {
    return scenario->getMarketplace()->getPrice( name, regionName, aPeriod, true );
}

/*! \brief Operate the capital of the sector.
* \details This is the main function for the sector called in each iteration.
*          The sequence of operations it performs is:
* <ol><li>Operate existing vintage capital stock.</li>
* <li>Invest in new vintages.</li>
* <li>Operate newly created capital stock.</li>
* \param aDemographic The demographics object.
* \param aNationalAccount The national accounts container.
* \param aPeriod The period in which to operate capital.
*/
void ProductionSector::operate( NationalAccount& aNationalAccount, const Demographic* aDemographics,
                                const int aPeriod )
{
    calcPriceReceived( aPeriod );
    operateOldCapital( aDemographics, aNationalAccount, aPeriod );
    calcInvestment( aDemographics, aNationalAccount, aPeriod );
    operateNewCapital( aDemographics, aNationalAccount, aPeriod );
}

/*! \brief Calculate new investment for the sector.
* \pre operateOldCapital has already been called to set the level of output from
*      existing capital without any new investment.
* \param aDemographic The demographics object.
* \param aNationalAccount The national accounts container.
* \param aPeriod The period in which to calculate investment.
*/
void ProductionSector::calcInvestment( const Demographic* aDemographic,
                                       NationalAccount& aNationalAccount,
                                       const int aPeriod )
{
    vector<IInvestable*> investableSubsecs = InvestmentUtils::convertToInvestables( subsec );
    // Calculate and distribute investment to the subsectors of this sector.
    mInvestor->calcAndDistributeInvestment( investableSubsecs, aNationalAccount, aDemographic,
        aPeriod );
    // Set efficiency conditions here.
    // Note: efficiency conditions do not change until the next iteration and so setting them 
    // before or after calculating and distributing investments does not matter unless using an
    // output share levelized cost calculator in which case it must be done after the distribution
    // TODO: perhaps we should collapse setEfficiencyConditions back into calcAndDistributeInvestment
    // so that we do not need to worry about ordering issues anymore
    mInvestor->setEfficiencyConditions( investableSubsecs, aNationalAccount, aDemographic, aPeriod );
}

/*! \brief Operate the old capital for the sector.
* \details Operate the old investment for the sector to determine the level of
*          output without any new investment.
* \param aDemographic The demographics object.
* \param aNationalAccount The national accounts container.
* \param aPeriod The period in which to operate.
*/
void ProductionSector::operateOldCapital( const Demographic* aDemographic, NationalAccount& aNationalAccount,
                                          const int period )
{
    for( CSubsectorIterator currSub = subsec.begin(); currSub != subsec.end(); ++currSub ){
        // flag tells the subsector only to operate old capital.
        (*currSub)->operate( aNationalAccount, aDemographic, moreSectorInfo.get(), false, period );
    }
}

/*! \brief Operate the new capital for the sector.
* \details Operate the new investment for the sector to determine the final
*          output level for the sector. This is always called after investment
*          is determined.
* \param aDemographic The demographics object.
* \param aNationalAccount The national accounts container.
* \param aPeriod The period in which to operate.
*/
void ProductionSector::operateNewCapital( const Demographic* aDemographic, NationalAccount& aNationalAccount,
                                          const int aPeriod )
{
    for( CSubsectorIterator currSub = subsec.begin(); currSub != subsec.end(); ++currSub ){
        // flag tells the subsector to operate new capital.
        (*currSub)->operate( aNationalAccount, aDemographic, moreSectorInfo.get(), true, aPeriod );
    }
}

/*! \brief Function to finalize objects after a period is solved.
* \details This function is used to calculate and store variables which are only needed after the current
* period is complete.
* \param aPeriod The period to finalize.
* \todo Finish this function.
* \author Sonny Kim
*/
void ProductionSector::postCalc( const int aPeriod ){

    // Call base class postCalc.
    Sector::postCalc( aPeriod );

    // TODO: is there any thing else that needs to get postCalc, if not we should just
    // get rid of it in ProductionSector
}

/*! \brief Get the XML node name for output to XML.
* \details This public function accesses the private constant string, XML_NAME.
*          This way the tag is always consistent for both read-in and output and
*          can be easily changed. This function may be virtual to be overridden
*          by derived class pointers.
* \author Josh Lurz, Sonny Kim
* \return The constant XML_NAME.
*/
const std::string& ProductionSector::getXMLName() const {
    return getXMLNameStatic();
}

/*! \brief Get the XML node name in static form for comparison when parsing XML.
* \details This public function accesses the private constant string, XML_NAME.
*          This way the tag is always consistent for both read-in and output and
*          can be easily changed. The "==" operator that is used when parsing,
*          required this second function to return static.
* \note A function cannot be static and virtual.
* \author Josh Lurz, Sonny Kim
* \return The constant XML_NAME as a static.
*/
const std::string& ProductionSector::getXMLNameStatic() {
    const static string XML_NAME = "productionSector";
    return XML_NAME;
}

/*! \brief Calculate the price received for the sector good.
* \details Calculates the price received for the output good of this sector by
*          adjusting the market price for transportation costs and taxes. This
*          price received is then set into the IInfo for this sector so it
*          can be accessed from other places in the model.
* \param period Period for which to calculate the price received.
*/
void ProductionSector::calcPriceReceived( const int period ){

    Marketplace* marketplace = scenario->getMarketplace();
    // set price received in market info
    double priceReceived = ( marketplace->getPrice( name, regionName, period ) + 
        ( moreSectorInfo->getValue( MoreSectorInfo::TRANSPORTATION_COST )
        * moreSectorInfo->getValue( MoreSectorInfo::TRAN_COST_MULT ) ) )
        / ( 1 + moreSectorInfo->getValue( MoreSectorInfo::IND_BUS_TAX_RATE ) );
    // set price received in market info
    FunctionUtils::setPriceReceived( regionName, name, period, priceReceived );
}

/*! \brief Update an aVisitor for reporting.
* \details Updates an aVisitor with information specific to the Sector
*          and ProductionSector. This is done by calling a class specific method
*          of the aVisitor so that it can update its information
*          pertaining to the ProductionSector.
* \param aVisitor OutputContainer to update.
* \param aPeriod Period in which to perform the update.
*/
void ProductionSector::accept( IVisitor* aVisitor, const int aPeriod ) const {
    // Update the output container for the derived class.
    aVisitor->startVisitProductionSector( this, aPeriod );
    // Update the base class
    Sector::accept( aVisitor, aPeriod );
    aVisitor->endVisitProductionSector( this, aPeriod );
}
