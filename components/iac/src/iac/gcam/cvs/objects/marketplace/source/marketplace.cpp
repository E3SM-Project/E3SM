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
* \file marketplace.cpp
* \ingroup Objects
* \brief Marketplace class source file.
* \author Sonny Kim
*/

#include "util/base/include/definitions.h"

#include <vector>

#include "marketplace/include/marketplace.h"
#include "marketplace/include/market.h"
#include "marketplace/include/imarket_type.h"
#include "containers/include/scenario.h"
#include "util/base/include/model_time.h"
#include "util/base/include/xml_helper.h"
#include "util/base/include/configuration.h"
#include "util/logger/include/ilogger.h"
#include "marketplace/include/price_market.h"
#include "marketplace/include/market_locator.h"
#include "util/base/include/ivisitor.h"
#include "containers/include/iinfo.h"
#include "marketplace/include/cached_market.h"

using namespace std;

extern Scenario* scenario;
const double Marketplace::NO_MARKET_PRICE = util::getLargeNumber();

/*! \brief Default constructor 
*
* The default constructor for the Marketplace which initializes several datamembers and
* creates an instance of the selected solver.
*
* \todo Marketplace might be better as a singleton.
*/
Marketplace::Marketplace(): mMarketLocator( new MarketLocator() ){
}

/*! \brief Destructor
*
* Destructor for the marketplace which deletes all the markets.
*/
Marketplace::~Marketplace() {
    // Clean up the markets.
    for ( vector<vector<Market*> >::iterator outerIter = markets.begin();
        outerIter != markets.end(); outerIter++ ) {
            for( vector<Market*>::iterator innerIter = outerIter->begin();
                innerIter != outerIter->end(); innerIter++ ) {
                    delete *innerIter;
                }
        }
}

/*! \brief Get the XML node name in static form for outputting XML.
* \author Josh Lurz
* \return The constant XML element name string as a static.
*/
const string& Marketplace::getXMLNameStatic() {
    const static string XML_NAME = "Marketplace";
    return XML_NAME;
}

/*! \brief Write out XML for debugging purposes.
*
* This method is called hierarchically from the main loop to write out 
* the current state of the model at a given time period. All member variables relevant to the running 
* of the model are written out to assist with debugging. It writes to the output stream
* passed as an argument in an XML format. 
* 
* \warning Currently to limit the size of the file which is written this method is only called for the
* US region.
* \param period The period for which to print the debugging information.
* \param out The output stream to which to print.
* \param tabs Tabs object used to track the number of tabs to print.
* \author Josh Lurz
*/
void Marketplace::toDebugXML( const int period, ostream& out, Tabs* tabs ) const {

    XMLWriteOpeningTag( getXMLNameStatic(), out, tabs );

    // write the xml for the class members.
    XMLWriteElement( static_cast<int>( markets.size() ), "numberOfMarkets", out, tabs );

    // Write out the individual markets
    const static string debugRegion =
        Configuration::getInstance()->getString( "debug-region", "USA" );

    for( unsigned int i = 0; i < markets.size(); i++ ){
        // TODO: This isn't quite right. This should search the contained
        // region list.
        if( markets[ i ][ period ]->getRegionName() == debugRegion ||
            markets[ i ][ period ]->getRegionName() == "global" )
        {
            markets[ i ][ period ]->toDebugXML( period, out, tabs );
        }
    }

    // finished writing xml for the class members.
    XMLWriteClosingTag( getXMLNameStatic(), out, tabs );
}

/*! \brief This function creates a market of the specified type for a given market region and good if it does not already exist.
*
* This function first checks if a market exists for a given market name and good name. If it does, the function add the region 
* to the contained regions of the market and add the region and good name to the regionToMarketMap, then finally returns false.
* Otherwise, it creates a market of the specified type, add the region to the contained region list, and adds keys to the marketMap and regionToMarketMap.
* The key to the marketMap is based on the marketName and goodName and the key to the regionToMarketMap is based on the regionName and goodName.
*
* \todo A factory method should be used so all market classes don't have to be included by this class. 
* \warning There is an important distinction here between the region name vs the market name. The key to the market is the goodName + market name.
* \param regionName The region of the sector for which to create a market.
* \param marketName The market region for which to create a market. This varies from the regionName, it can be global, a multi-region market, or the same as the region name.
* \param goodName The good for which to create a market.
* \param aType The type of market to create.
* \return Whether a market was created.
*/
bool Marketplace::createMarket( const string& regionName, const string& marketName, const string& goodName, const IMarketType::Type aType ) {
    /*! \pre Region name, market name, and sector name must be non null. */
    assert( !regionName.empty() && !marketName.empty() && !goodName.empty() );

    // Create the index within the market locator.
    const int uniqueNumber = static_cast<int>( markets.size() );
    int marketNumber = mMarketLocator->addMarket( marketName, regionName, goodName, uniqueNumber );

    // If the market number is the unique number we passed it, the market did not already exist and 
    // we should create the market objects, one per period.
    const bool isNewMarket = ( marketNumber == uniqueNumber );
    if( isNewMarket ){
        vector<Market*> tempVector( scenario->getModeltime()->getmaxper() );
        for( unsigned int i = 0; i < tempVector.size(); i++ ){
            tempVector[ i ] = Market::createMarket( aType, goodName, marketName, i ).release();
        }
        markets.push_back( tempVector );
    }

    // Add the region onto the market.
    for( unsigned int i = 0; i < markets[ marketNumber ].size(); i++ ) {
        markets[ marketNumber ][ i ]->addRegion( regionName );
    }
    // Return whether we were required to create a new market.
    return isNewMarket;
}

/*! \brief Restructures a market to account for simultaneities.
*
* Changes the named market to a price market, which supplies a trial price for a secondary good.
* It also adds a corresponding demand market that provides a trial value for demand.
* Markets are subclassed to allow Market::addToDemand, Market::getDemand, Market::getPrice, etc. to act 
* differently for PriceMarket and DemandMarket so that these changes are transparent to the 
* rest of the code.  Only secondary goods are allowed to restructure their markets;
* any markets set to solve will be ignored.
* 
* \author Steve Smith
* \param goodName The name of the good of the market to be restructured.
* \param regionName The region of the market to be restructured.
*/
void Marketplace::resetToPriceMarket( const string& goodName, const string& regionName ) {
    // If simuls are off, return immediately. 
    if( !Configuration::getInstance()->getBool( "simulActive" ) ){
        return;
    }

    const int marketNumber = mMarketLocator->getMarketNumber( regionName, goodName );

    if( marketNumber == MarketLocator::MARKET_NOT_FOUND ) {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "Cannot reset Market "<< goodName << " to a price market because it does not exist."  << endl;
    }
    else if( markets[ marketNumber][ 0 ]->getType() != IMarketType::NORMAL ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "Cannot reset market type other than normal to a price market." << endl;
    }
    else if( markets[ marketNumber ][ 1 ]->isSolvable() ) {
        // Solved markets do not need to be split.
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Solved markets do not need trial price/demand markets." << endl;
    }
    else {
        // Setup the corresponding demand markets
        string marketName = markets[ marketNumber ][ 0 ]->getRegionName();
        string demandGoodName = goodName + "Demand_int";
        createMarket( regionName, marketName, demandGoodName, IMarketType::DEMAND );
        // Add units of the corresponding NORMAL market to the DEMAND market info object.
        IInfo* marketInfoFrom = getMarketInfo( goodName, regionName, 0, true );
        IInfo* marketInfo = getMarketInfo( demandGoodName, regionName, 0, true );
        marketInfo->setString( "price-unit", marketInfoFrom->getString( "price-unit", true ) );
        marketInfo->setString( "output-unit", marketInfoFrom->getString( "output-unit", true ) );

        int demandMarketNumber = mMarketLocator->getMarketNumber( regionName, demandGoodName );
        assert( demandMarketNumber != MarketLocator::MARKET_NOT_FOUND );
        // loop through time periods            
        for( unsigned int per = 1; per < markets[ marketNumber ].size(); per++ ){

            // Get the pointer of the new demand market.
            Market* newDemandMarket = markets[ demandMarketNumber ][ per ];
            assert( newDemandMarket );
            
            // Check if the old market had set an initial trial demand which we
            // could use as the initial price.
            Market* oldMarket = markets[ marketNumber ][ per ];
            marketInfoFrom = oldMarket->getMarketInfo();
            double initialDemand = marketInfoFrom->getDouble( "initial-trial-demand", false );
            newDemandMarket->setPrice( initialDemand );

            // Create a new price market from the old market.
            Market* newPriceMarket = new PriceMarket( *oldMarket, newDemandMarket );

            // Delete the old market.
            delete oldMarket;

            // Insert the new price market. 
            markets[ marketNumber ][ per ] = newPriceMarket;

            // Set both markets to solve.
            setMarketToSolve( goodName, regionName, per );
            setMarketToSolve( demandGoodName, regionName, per );
            // this assumes that all markets have the same number of periods
        }
    }
}

/*! \brief Set the prices by period of a market from a vector.
*
* This function sets the price for each period of a market according to the corresponding value in 
* the prices vector.  The function checks so that it will not fail if the prices vector is a 
* different size than the markets vector. 
*
* \author Josh Lurz
* \param goodName The goodName of the market for which to set new prices.
* \param regionName The regionName to use for the lookup to determine the correct market.
* \param prices A vector containing prices to set into the market. 
*/
void Marketplace::setPriceVector( const string& goodName, const string& regionName,
                                 const vector<double>& prices ){
    // determine what market the region and good are in.
    const int marketNumber = mMarketLocator->getMarketNumber( regionName, goodName );
    if( marketNumber == MarketLocator::MARKET_NOT_FOUND ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Market price vector cannot be set the market does not exist: " 
            << goodName << " " << regionName << endl;
    }
    else {
        for( unsigned int i = 0; i < markets[ marketNumber ].size() && i < prices.size(); i++ ){
            markets[ marketNumber ][ i ]->setRawPrice( prices[ i ] );
        }
    }
}

/*! \brief Initialize prices for all markets. 
*
* Supply and demand sector prices should always get set somewhere else except for in the first period.
* However, no such guarantee exists for GHG markets. The function initializes prices for all periods 
* in all markets just to be safe if they are not already set. 
* Initialization also occurs for supply and demand markets that have prices read-in via routine:
* setPriceVector in sector as supply and demand markets are created.
* This has no effect for future periods as these prices are overwritten by 
* Marketplace::init_to_last except for calibration markets.
*/
void Marketplace::initPrices(){
    // initialize supply and demand sector market prices to 1.
    for ( unsigned int i = 0; i < markets.size(); ++i ){
        for( unsigned int j = 0; j < markets[ i ].size(); ++j ){
            markets[ i ][ j ]->initPrice();
        }
    }
}

/*! \brief Set the solve flag for this market for the given period, or all periods if per argument 
* is undefined.
*
* This function determines a market from a good and region and sets the market to solve for the 
* period passed to the function.  This solve flag determines whether the market is solved by the
* solution mechanism, except for cases where the market does not pass certain other criteria 
* related to singularities. If this flag is set to false, as is the default, the market will never be solved. 
* \param goodName The name of the good of the market.
* \param regionName The region name of the market.
* \param per The period for which the market should be solved.
*/
void Marketplace::setMarketToSolve ( const string& goodName, const string& regionName, const int per ) {
    const int marketNumber = mMarketLocator->getMarketNumber( regionName, goodName );

    // If the market exists.
    if ( marketNumber != MarketLocator::MARKET_NOT_FOUND ) {
        markets[ marketNumber ][ per ]->setSolveMarket( true );
    }
    else {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Market cannot be set to solve as it does not exist: " << goodName << " " 
            << regionName << endl;
    }
}

/*! \brief Unset the solve flag for this market for the given period, or all periods if per argument
* is undefined.  This function determines a market from a good and region and unsets the market to
* solve for the period passed to the function.  This solve flag determines whether the market is
* solved by the solution mechanism, except for cases where the market does not pass certain
* other criteria related to singularities. If this flag is set to false, as is the default, the
* market will never be solved.  This function also clears supply and demand for the market at the
* same time. 
* \todo This function is only used in one place.
* \param goodName The name of the good of the market.
* \param regionName The region name of the market.
* \param per The period for which the market should not be solved.
*/
void Marketplace::unsetMarketToSolve ( const string& goodName, const string& regionName, const int per ) {

    const int marketNumber = mMarketLocator->getMarketNumber( regionName, goodName );

    // If the market exists.
    if ( marketNumber != MarketLocator::MARKET_NOT_FOUND ) {
        markets[ marketNumber ][ per ]->setSolveMarket( false );
        markets[ marketNumber ][ per ]->nullSupply();
        markets[ marketNumber ][ per ]->nullDemand();
    }
    else {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Market cannot be unset not to solve as it does not exist: " << goodName << " " 
            << regionName << endl;
    }
}

/*! \brief Clear all market supplies and demands for the given period.
* 
* This function iterates through the markets and nulls the supply and demand 
* of each market in the given period.
*
* \param period Period in which to null the supplies and demands. 
*/
void Marketplace::nullSuppliesAndDemands( const int period ) {
    for ( unsigned int i = 0; i < markets.size(); i++ ) {
        markets[ i ][ period ]->nullDemand();
        markets[ i ][ period ]->nullSupply();
    }
}

/*! \brief Set the market price.
*
* This function uses the type dependent method Market::setPrice to set the passed in value into
* the market.
*
* \param goodName The good of the market.
* \param regionName The region setting the price.
* \param value The value to which to set price.
* \param per The period in which to set the price.
* \param aMustExist Whether it is an error for the market not to exist.
*/
void Marketplace::setPrice( const string& goodName, const string& regionName, const double value,
                            const int per, bool aMustExist )
{
    // Print a warning message if the new price is not a finite number.
    if ( !util::isValidNumber( value ) ) {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Error setting price in markeplace for: " << goodName << ", value: " << value << endl;
        return;
    }

    const int marketNumber = mMarketLocator->getMarketNumber( regionName, goodName );
    if ( marketNumber != MarketLocator::MARKET_NOT_FOUND ) {
        markets[ marketNumber ][ per ]->setPrice( value );
    }
    else if( aMustExist ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Cannot set price for market as it does not exist: " << goodName << " " 
            << regionName << endl;
    }
}

/*! \brief Add to the supply for this market.
*
* This function increments the supply for a market determined by the goodName and regionName
* by a given value. This function is used throughout the model to add supply to markets. 
*
* \param goodName Name of the good for which to add supply.
* \param regionName Name of the region in which supply should be added for the market.
* \param value Amount of supply to add.
* \param per Period in which to add supply.
* \param aMustExist Whether it is an error for the market not to exist.
*/
void Marketplace::addToSupply( const string& goodName, const string& regionName, const double value,
                               const int per, bool aMustExist )
{
    // Print a warning message when adding infinity values to the supply.
    if ( !util::isValidNumber( value ) ) {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::NOTICE );
		mainLog << "Error adding to supply in markeplace for: " << goodName << ", region: " << regionName << ", value: " << value << endl;
        return;
    }

    const int marketNumber = mMarketLocator->getMarketNumber( regionName, goodName );

    if ( marketNumber != MarketLocator::MARKET_NOT_FOUND ) {
        markets[ marketNumber ][ per ]->addToSupply( value );
    }
    else if( aMustExist ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Cannot add to supply for market as it does not exist: " << goodName << " " 
            << regionName << endl;
    }
}

/*! \brief Add to the demand for this market.
*
* This function increments the demand for a market determined by the goodName and regionName
* by a given value. This function is used throughout the model to add demand to markets. 
*
* \param goodName Name of the good for which to add demand.
* \param regionName Name of the region in which demand should be added for the market.
* \param value Amount of demand to add.
* \param per Period in which to add demand.
* \param aMustExist Whether it is an error for the market not to exist.
*/
void Marketplace::addToDemand( const string& goodName, const string& regionName, const double value,
                               const int per, bool aMustExist )
{
    // Print a warning message when adding infinity values to the demand
    if ( !util::isValidNumber( value ) ) {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Error adding to demand in markeplace for: " << goodName << ", region: " << regionName << ", value: " << value << endl;
        return;
    }

    const int marketNumber = mMarketLocator->getMarketNumber( regionName, goodName );
    if ( marketNumber != MarketLocator::MARKET_NOT_FOUND ) {
        markets[ marketNumber ][ per ]->addToDemand( value );
    }
    else if( aMustExist ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Cannot add to demand for market as it does not exist: " << goodName << " " 
            << regionName << endl;
    }
}

/*! \brief Return the market price. 
*
* This function uses a market type dependent function to find the price for a market determined by
* the goodName and regionName.  This price is not always the raw or true price. For non-existant
* markets, this function returns a near infinite price, except for renewable
* markets that do not exist, for which it returns 0.
*
* \param goodName The good for which a price is needed.
* \param regionName The region for which a price is needed.
* \param per The period to return the market price for.
* \param aMustExist Whether it is an error for the market not to exist.
* \return The market price.
*/  
double Marketplace::getPrice( const string& goodName, const string& regionName, const int per,
                             bool aMustExist ) const {
    const int marketNumber = mMarketLocator->getMarketNumber( regionName, goodName );
    
    if( marketNumber != MarketLocator::MARKET_NOT_FOUND ){
        return markets[ marketNumber ][ per ]->getPrice();
    }

    if( aMustExist ) {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Called for price of non-existant market " << goodName << " in region " 
            << regionName << endl;
    }
    return NO_MARKET_PRICE;
}

/*! \brief Return the market supply. 
*
* This function uses a market type dependent function to find the supply for a market determined
* by the goodName and regionName.  This supply is not always the raw or true supply.
* For non-existant markets, this function returns 0.
*
* \param goodName The good for which a supply is needed.
* \param regionName The region for which a supply is needed.
* \param per Period to get the supply for.
* \return The market supply.
*/
double Marketplace::getSupply( const string& goodName, const string& regionName, const int per ) const {
    const int marketNumber = mMarketLocator->getMarketNumber( regionName, goodName );

    if ( marketNumber != MarketLocator::MARKET_NOT_FOUND ) {
        return markets[ marketNumber ][ per ]->getSupply();
    }

    ILogger& mainLog = ILogger::getLogger( "main_log" );
    mainLog.setLevel( ILogger::NOTICE );
    mainLog << "Called for supply of non-existant market " << goodName << " in " << regionName << endl;
    return 0;
}

/*! \brief Return the market demand. 
*
* This function uses a market type dependent function to find the demand for a market determined
* by the goodName and regionName.  This demand is not always the raw or true demand.
* For non-existant markets, this function returns 0.
*
* \param goodName The good for which a demand is needed.
* \param regionName The region for which a demand is needed.
* \param per The period to return the market demand for.
* \return The market demand.
*/
double Marketplace::getDemand(  const string& goodName, const string& regionName, const int per ) const {
    const int marketNumber = mMarketLocator->getMarketNumber( regionName, goodName );

    if ( marketNumber != MarketLocator::MARKET_NOT_FOUND ) {
        return markets[ marketNumber ][ per ]->getDemand();
    }

    ILogger& mainLog = ILogger::getLogger( "main_log" );
    mainLog.setLevel( ILogger::NOTICE );
    mainLog << "Called for demand of non-existant market " << goodName << " in " << regionName << endl;
    return 0;
}

/*! \brief Return the market stored supply. 
*
* This function uses a market type dependent function to find the stored supply for a market
* determined by the goodName and regionName.  This stored supply is not always the raw or true
* supply. For non-existant markets, this function returns 0.
*
* \param goodName The good for which a stored supply is needed.
* \param regionName The region for which a stored supply is needed.
* \param per Period to get the stored supply for.
* \return The stored supply.
*/
double Marketplace::getStoredSupply( const string& goodName, const string& regionName,
                                    const int per ) const {
    const int marketNumber = mMarketLocator->getMarketNumber( regionName, goodName );

    if ( marketNumber != MarketLocator::MARKET_NOT_FOUND ) {
        return markets[ marketNumber ][ per ]->getStoredRawSupply();
    }

    ILogger& mainLog = ILogger::getLogger( "main_log" );
    mainLog.setLevel( ILogger::NOTICE );
    mainLog << "Called for stored supply of non-existant market " << goodName << " in " 
        << regionName << endl;
    return 0;
}

/*! \brief Return the market stored demand. 
*
* This function uses a market type dependent function to find the stored demand for a market
* determined by the goodName and regionName.  This stored demand is not always the raw or true
* demand.  For non-existant markets, this function returns 0.
*
* \param goodName The good for which a stored demand is needed.
* \param regionName The region for which a stored demand is needed.
* \param per The period to return the stored demand for.
* \return The stored demand.
*/
double Marketplace::getStoredDemand( const string& goodName, const string& regionName,
                                    const int per ) const {
    const int marketNumber = mMarketLocator->getMarketNumber( regionName, goodName );

    if ( marketNumber != MarketLocator::MARKET_NOT_FOUND ) {
        return markets[ marketNumber ][ per ]->getStoredRawDemand();
    }

    ILogger& mainLog = ILogger::getLogger( "main_log" );
    mainLog.setLevel( ILogger::NOTICE );
    mainLog << "Called for stored demand of non-existant market " << goodName << " in " 
        << regionName << endl;
    return 0;
}

//! Returns a set of pointers to each market for the period.
vector<Market*> Marketplace::getMarketsToSolve( const int period ) const {
    vector<Market*> toSolve;
    
    // Loop through the markets and add all markets.
    for( unsigned int i = 0; i < markets.size(); ++i ){
        toSolve.push_back( markets[ i ][ period ] );
    }
    return toSolve;
}

/*!
 * \brief Conditionally initializes the market prices to the market prices from
 *        the previous period.
 * \details Resuse parsed prices up to the configuration parameter restart-period
 *          which will allow users to resume the model at that period.  If the
 *          user did not intend on restarting the model we can resuse prices up to
 *          the final calibration year.  After which we would assume the scenario
 *          would change significantly enough that starting from previous period's
 *          prices would be closer to the solution. This only occurs for periods
 *          greater than 0.
 * \author Sonny Kim
 * \param period Period for which to initialize prices.
 */
void Marketplace::init_to_last( const int period ) { 
    // Get the last period to allow using parsed prices, the default is the
    // final calibration period.
    const static int restartPeriod = max( Configuration::getInstance()->getInt(
        "restart-period", scenario->getModeltime()->getFinalCalibrationPeriod() + 1 ),
        scenario->getModeltime()->getFinalCalibrationPeriod() + 1 );
    if ( period > 0 && period < restartPeriod ) {
        for ( unsigned int i = 0; i < markets.size(); i++ ) {
            markets[ i ][ period ]->set_price_to_last_if_default( markets[ i ][ period - 1 ]->getRawPrice() );
        }
    }
    else if( period >= restartPeriod ){
        for ( unsigned int i = 0; i < markets.size(); i++ ) {
            markets[ i ][ period ]->set_price_to_last( markets[ i ][ period - 1 ]->getRawPrice() );
        }
    }
}

/*! \brief Store the demand, supply and price for each market. 
*
* This function called the Market::storenfo function on each market in the marketplace,
* which causes them to set their stored demand, supply, and price to the respective current values
* for those variables.
*
* \author Sonny Kim
* \param period Period for which to store demands, supplies and prices.
*/
void Marketplace::storeinfo( const int period ) {
    for ( unsigned int i = 0; i  < markets.size(); i++ ) {
        markets[ i ][ period ]->storeInfo();
    }
}

/*! \brief Restore the stored demand, supply and price for each market. 
*
* This function called the Market::restoreInfo function on each market in the marketplace,
* which causes them to set their demand, supply, and price to the respective stored values
* for those variables.
*
* \param period Period for which to restore demands, supplies, and prices.
*/
void Marketplace::restoreinfo( const int period) {
    for ( unsigned int i = 0; i < markets.size(); i++ ) {
        markets[ i ][ period ]->restoreInfo();
    }
}

/*! \brief Store market prices for policy cost caluclation.
*
*
* \author Sonny Kim
*/
void Marketplace::store_prices_for_cost_calculation()
{
    for ( unsigned int i = 0; i < markets.size(); ++i ){
        for( unsigned int period = 0; period < scenario->getModeltime()->getmaxper(); ++period ){
            markets[ i ][ period ]->store_original_price();
        }
    }
}

/*! \brief Restore market prices for policy cost caluclation.
*
*
* \author Sonny Kim
*/
void Marketplace::restore_prices_for_cost_calculation()
{
    for ( unsigned int i = 0; i < markets.size(); ++i ){
        for( unsigned int period = 0; period < scenario->getModeltime()->getmaxper(); ++period ){
            markets[ i ][ period ]->restore_original_price();
        }
    }
}

/*! \brief Get the information object for the specified market and period which
*          can then be used to query for specific values.
* \details Returns the internal IInfo object of the specified market and period
*          which represents a set of pairings of information name to value. This
*          specific function returns a constant pointer to the information
*          object, so values can be queried but not added or modified.
* \author Josh Lurz
* \warning Null will be returned if the market does not exist.
* \warning The function is UNRELATED to storeinfo and restoreinfo.
* \param aGoodName The good of the market for which to get the information object.
* \param aRegionName The region used to find the market from which to get the
*        information object.
* \param aPeriod The period to fetch for which the information object. 
* \param aMustExist Whether it is an error for the market not to exist.
* \return A constant pointer to the market information object, null if the
*         market does not exist.
* \note This version of the function is required so that it can be called in
*       constant functions. A second version is available which returns a
*       mutable pointer.
*/
const IInfo* Marketplace::getMarketInfo( const string& aGoodName, const string& aRegionName,
                                         const int aPeriod, const bool aMustExist ) const 
{
    const int marketNumber = mMarketLocator->getMarketNumber( aRegionName, aGoodName );
    const IInfo* info = 0;
    if ( marketNumber != MarketLocator::MARKET_NOT_FOUND ) {
        info = markets[ marketNumber ][ aPeriod ]->getMarketInfo();
        /*! \invariant The market is required to return an information object
        *              that is non-null. 
        */
        assert( info );
    }
    
    // Report the error if requested.
    if( !info && aMustExist ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Market info object cannot be returned because market "
                << aGoodName << " in " << aRegionName << " does not exist." << endl;
    }
    return info;
}

/*! \brief Get the information object for the specified market and period which
*          can then be used to query, add, or modify specific values.
* \details Returns the internal IInfo object of the specified market and period
*          which represents a set of pairings of information name to value. This
*          specific function returns a mutable pointer to the information
*          object, so values can be queried, added or modified.
* \author Josh Lurz
* \warning Null will be returned if the market does not exist.
* \warning The function is UNRELATED to storeinfo and restoreinfo.
* \param aGoodName The good of the market for which to get the information
*        object.
* \param aRegionName The region used to find the market from which to get the
*        information object.
* \param aPeriod The period to fetch for which the information object. 
* \param aMustExist Whether it is an error for the market not to exist.
* \return A mutable pointer to the market information object, null if the market
*         does not exist.
* \note This function returns a mutable pointer to the information object so it
*       cannot be called from constant function.
*/
IInfo* Marketplace::getMarketInfo( const string& aGoodName, const string& aRegionName,
                                   const int aPeriod, const bool aMustExist )
{
    const int marketNumber = mMarketLocator->getMarketNumber( aRegionName, aGoodName );
    IInfo* info = 0;
    if ( marketNumber != MarketLocator::MARKET_NOT_FOUND ) {
        info = markets[ marketNumber ][ aPeriod ]->getMarketInfo();
        /*! \invariant The market is required to return an information object
        *              that is non-null. 
        */
        assert( info );
    }
    
    // Report the error if requested.
    if( !info && aMustExist ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Market info object cannot be returned because market " 
                << aGoodName << " in " << aRegionName << " does not exist." << endl;
    }
    return info;
}

/*!
 * \brief Locate the market for the given good, region, and period.
 * \details Returns a CachedMarket that wraps the found Market object so that
 *          only functionality provided through the marketplace will be available.
 *          Note that it is not an error to locate a market which does not exist
 *          and calls to methods of a CachedMarket which did not exist will have the
 *          same behavior as the equivalent method in this class.
 * \param aGoodName The good of the market to locate.
 * \param aRegionName The region of the market to locate.
 * \param aPeriod The period for which to locate.
 * \return A CachedMarket object which wraps the requested market.  This will always
 *         be a valid object regardless of if the market was not found.
 * \see CachedMarket
 */
auto_ptr<CachedMarket> Marketplace::locateMarket( const string& aGoodName, const string& aRegionName,
                                                  const int aPeriod ) const
{
    const int marketNumber = mMarketLocator->getMarketNumber( aRegionName, aGoodName );
    auto_ptr<CachedMarket> locatedMarket( new CachedMarket( aGoodName, aRegionName, aPeriod,
                                                            marketNumber != MarketLocator::MARKET_NOT_FOUND ?
                                                            markets[ marketNumber ][ aPeriod ] : 0 ) );
    return locatedMarket;
}

/*! \brief Write out the market information to the database.
*
* This function is used to perform a data writeout to the database.
*
* \note This will be replaced by toXMLOutput
*/
void Marketplace::dbOutput() const {

    const Modeltime* modeltime = scenario->getModeltime();

    void dboutput4(string var1name,string var2name,string var3name,string var4name,
        string uname,vector<double> dout);

    const int maxPeriod = modeltime->getmaxper();
    vector<double> temp(maxPeriod);
    int j;
    // write market prices, supply and demand
    for (int i=0;i< static_cast<int>( markets.size() );i++) {
        string tempRegName = markets[i][0]->getRegionName();
        string tempGoodName = markets[i][0]->getGoodName();
        const IInfo* marketInfo = markets[i][0]->getMarketInfo();
 
        for (j=0;j<maxPeriod;j++) {
            temp[j] = markets[i][j]->getRawPrice();
        }
        dboutput4(markets[i][0]->getRegionName(),"Market",markets[i][0]->getGoodName(),"1_price",
            marketInfo->getString( "price-unit", true ), temp);
        for (j=0;j<maxPeriod;j++) {
            temp[j] = markets[i][j]->getRawSupply();
        }
        dboutput4(markets[i][0]->getRegionName(),"Market",markets[i][0]->getGoodName(),"2_supply",
            marketInfo->getString( "output-unit", true ), temp);
        for (j=0;j<maxPeriod;j++) {
            temp[j] = markets[i][j]->getRawDemand();
        }
        dboutput4(markets[i][0]->getRegionName(),"Market",markets[i][0]->getGoodName(),"3_demand",
            marketInfo->getString( "output-unit", true ), temp);
    }
}

/*! \brief Write out market information to a file.
*
* This function is used to perform a data writeout to a plain text file.
*
* \note This will be replaced by toXMLOutput
*/
void Marketplace::csvOutputFile( string marketsToPrint ) const {

    const Modeltime* modeltime = scenario->getModeltime();

    // function protocol
    void fileoutput2(string var1name,string var2name,string var3name,
        string var4name,string var5name,vector<double> dout,string uname);
    void fileoutput3(string var1name,string var2name,string var3name,
        string var4name,string var5name,string uname,vector<double> dout);

    const int maxPeriod = modeltime->getmaxper();
    vector<double> temp(maxPeriod);

    // Function arguments are variable name, double array, db name, and
    // table name.
    // The function writes all years.

    // write market prices and supply (or demand)
    for ( unsigned int i = 0; i < markets.size(); i++ ) {
      if ( marketsToPrint == "" || markets[i][0]->getRegionName() == marketsToPrint ) {
        const IInfo* marketInfo = markets[i][0]->getMarketInfo();
        for ( int j = 0; j < maxPeriod; j++ ) {
            temp[j] = markets[i][j]->getRawPrice();
        }
        fileoutput3(markets[i][0]->getRegionName(),"market",markets[i][0]->getGoodName()," ","price",
            marketInfo->getString( "price-unit", true ), temp);
        for ( int j = 0; j < maxPeriod; j++ ) {
            temp[j] = markets[i][j]->getRawSupply();
        }
        fileoutput3(markets[i][0]->getRegionName(),"market",markets[i][0]->getGoodName()," ","supply",
            marketInfo->getString( "output-unit", true ), temp);
      }
    }
}

/*! \brief For outputing SGM data to a flat csv File
 * 
 * \author Pralit Patel
 * \param period The period which we are outputing for
 */
void Marketplace::csvSGMOutputFile( ostream& aFile, const int period ) const {
    // come back to this
    aFile << "Region" << ',' << "Good" << ',' << "Price" << ','<< "PriceReceived" << ',' << "ED" 
        << ',' << "Demand" << ',' << "Supply" << endl;
    // reset format to default
    aFile.setf(ios_base::fixed, ios_base::floatfield);
    for ( unsigned int i = 0; i < markets.size(); i++ ) {
        aFile << markets[ i ][ period ]->getRegionName() << ',' 
              << markets[ i ][ period ]->getGoodName() << ','
              << markets[ i ][ period ]->getPrice() << ','
              << markets[ i ][ period ]->getMarketInfo()->getDouble( "priceReceived", false ) << ','
              << markets[ i ][ period ]->getDemand() - markets[ i ][ period ]->getSupply()  << ','
              << markets[ i ][ period ]->getDemand() << ',' << markets[ i ][ period ]->getSupply()
              << endl;
    }
    aFile << endl;
}

/*! \brief Update an output container for the Marketplace.
* \param aVisitor Output container to update.
* \param aPeriod Period to update.
*/
void Marketplace::accept( IVisitor* aVisitor, const int aPeriod ) const {
    aVisitor->startVisitMarketplace( this, aPeriod );

    // Update from the markets.
    for( unsigned int i = 0; i < markets.size(); i++ ){
        // If the period is -1 this means to update all periods.
        if( aPeriod == -1 ){
            for( unsigned int j = 0; j < markets[ i ].size(); ++j ){
                markets[ i ][ j ]->accept( aVisitor, aPeriod );
            }
        }
        // Otherwise only update for the current period.
        else {
            markets[ i ][ aPeriod ]->accept( aVisitor, aPeriod );
        }
    }

    aVisitor->endVisitMarketplace( this, aPeriod );
}
