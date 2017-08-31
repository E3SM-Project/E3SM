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
* \file market.cpp
* \ingroup Objects
* \brief Market class source file.
* \author Sonny Kim
*/

#include "util/base/include/definitions.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cassert>

#include "util/base/include/model_time.h" 
#include "util/base/include/xml_helper.h"
#include "marketplace/include/market.h"
#include "marketplace/include/imarket_type.h"
#include "containers/include/scenario.h"
#include "marketplace/include/price_market.h"
#include "marketplace/include/demand_market.h"
#include "marketplace/include/calibration_market.h"
#include "marketplace/include/inverse_calibration_market.h"
#include "marketplace/include/market_tax.h"
#include "marketplace/include/market_subsidy.h"
#include "marketplace/include/market_RES.h"
#include "marketplace/include/normal_market.h"
#include "marketplace/include/trial_value_market.h"
#include "util/base/include/ivisitor.h"
#include "containers/include/info_factory.h"
#include "util/base/include/atom_registry.h"
#include "util/base/include/atom.h"
#include "containers/include/iinfo.h"
#include "util/logger/include/ilogger.h"

using namespace std;
using namespace objects;

extern Scenario* scenario;

/*! \brief Constructor
* \details This is the constructor for the market class. No default constructor
*          exists to prevent the creation of empty markets. 
* \warning The arguments are required to define the good name, region name and
*          model period. These values are invariants.
* \param goodNameIn The good or fuel name for the item in the market.
* \param regionNameIn The region which this market covers. It may include
*        several model regions.
* \param periodIn The period the market exists in.
*/
Market::Market( const string& goodNameIn, const string& regionNameIn, const int periodIn )
: good( goodNameIn ), 
region( regionNameIn ),
solveMarket( false ),
period( periodIn ),
price( 0 ),
storedPrice( 0 ),
demand( 0 ),
storedDemand( 0 ),
supply( 0 ),
storedSupply( 0 ),
mMarketInfo( InfoFactory::constructInfo( 0, regionNameIn+goodNameIn ) )
{
    // Store the market name so that it can be returned without any allocations.
    mName = region + good;
}

//! Destructor. This is needed because of the auto_ptr.
Market::~Market(){
}

/*! \brief Protected copy constructor
* \details This copy constructor is needed because auto_ptr held memory cannot
*          be copied automatically. The copy constructor is protected because it
*          should only be accessed by the PriceMarket derived class.
* \param aMarket The market to copy.
* \author Josh Lurz
*/
Market::Market( const Market& aMarket ): 
good( aMarket.good ),
region( aMarket.region ),
mName( aMarket.mName ),
solveMarket( aMarket.solveMarket ),
period( aMarket.period ),
price( aMarket.price ),
storedPrice( aMarket.storedPrice ),
demand( aMarket.demand ),
storedDemand( aMarket.storedDemand ),
supply( aMarket.supply ),
storedSupply( aMarket.storedSupply ),
mContainedRegions( aMarket.mContainedRegions ),
mMarketInfo( InfoFactory::constructInfo( 0,aMarket.mName ) ){
    // TODO: Cannot currently copy the market info.
}

/*! \brief Static factory method to create a market based on its type.
* \details
* \param aGoodName The good or fuel name for the item in the market to create.
* \param aRegionName The region which the market to create covers.
* \param aPeriod The period the market to create exists in.
* \param aType Type of market to create.
* \return A pointer to the newly allocated market, null if the type did not
*         exist. 
*/
auto_ptr<Market> Market::createMarket( const IMarketType::Type aType,
                                       const string& aGoodName,
                                       const string& aRegionName,
                                       const int aPeriod )
{
    assert( aType < IMarketType::END );
    auto_ptr<Market> rNewMarket;
    if ( aType == IMarketType::NORMAL ){
        rNewMarket.reset( new NormalMarket( aGoodName, aRegionName, aPeriod ) );
    }
    else if ( aType == IMarketType::TAX ) {
        rNewMarket.reset( new MarketTax( aGoodName, aRegionName, aPeriod ) );
    }
    else if ( aType == IMarketType::RES ) {
        rNewMarket.reset( new MarketRES( aGoodName, aRegionName, aPeriod ) );
    }
	else if ( aType == IMarketType::SUBSIDY ) {
        rNewMarket.reset( new MarketSubsidy( aGoodName, aRegionName, aPeriod ) );
    }
    else if ( aType == IMarketType::CALIBRATION ) {
        rNewMarket.reset( new CalibrationMarket( aGoodName, aRegionName, aPeriod ) );
    }
    else if ( aType == IMarketType::INVERSE_CALIBRATION ) {
        rNewMarket.reset( new InverseCalibrationMarket( aGoodName, aRegionName, aPeriod ) );
    }
    else if ( aType == IMarketType::DEMAND ) {
        rNewMarket.reset( new DemandMarket( aGoodName, aRegionName, aPeriod ) );
    }
    else if ( aType == IMarketType::TRIAL_VALUE ) {
        rNewMarket.reset( new TrialValueMarket( aGoodName, aRegionName, aPeriod ) );
    }
    else if ( aType == IMarketType::PRICE ) {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "Price markets are only created internally in the marketplace." << endl;
    }
    else {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Invalid market type: " << aType << endl;
    }
    return rNewMarket;
}

/*! \brief Write out XML for debugging purposes.
* \details This method is called by the Marketplace::toDebugXML method to write
*          out information for each individual market. It prints the current
*          state of all internal variables. It also calls a derived method which
*          prints derived class specific information.
* \param period Model period for which to print information.
* \param out Output stream to print to.
* \param tabs A tabs object responsible for printing the correct number of tabs. 
*/
void Market::toDebugXML( const int period, ostream& out, Tabs* tabs ) const {
    const Modeltime* modeltime = scenario->getModeltime();
    XMLWriteOpeningTag( getXMLNameStatic(), out, tabs, getName(), modeltime->getper_to_yr( period ) , convert_type_to_string( getType() ) );
    XMLWriteElement( solveMarket, "solved_Market_Flag", out, tabs );
    XMLWriteElement( good, "MarketGoodOrFuel", out, tabs );
    XMLWriteElement( region, "MarketRegion", out, tabs );
    XMLWriteElement( price, "price", out, tabs );
    XMLWriteElement( storedPrice, "storedPrice", out, tabs );
    XMLWriteElement( demand, "demand", out, tabs );
    XMLWriteElement( storedDemand, "storedDemand", out, tabs );
    XMLWriteElement( supply, "supply", out, tabs );
    XMLWriteElement( storedSupply, "storedSupply", out, tabs );

    for( vector<const Atom*>::const_iterator i = mContainedRegions.begin(); i != mContainedRegions.end(); ++i ) {
        XMLWriteElement( (*i)->getID(), "ContainedRegion", out, tabs );
    }

    mMarketInfo->toDebugXML( period, tabs, out );

    toDebugXMLDerived( out, tabs );

    // finished writing xml for the class members.
    XMLWriteClosingTag( getXMLNameStatic(), out, tabs );
}

/*! \brief Get the XML node name in static form for comparison when parsing XML.
* \details This public function accesses the private constant string, XML_NAME.
*          This way the tag is always consistent for both read-in and output and
*          can be easily changed. The "==" operator that is used when parsing,
*          required this second function to return static.
* \note A function cannot be static and virtual.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME as a static.
*/
const std::string& Market::getXMLNameStatic() {
    const static string XML_NAME = "market";
    return XML_NAME;
}

/*! \brief Add a region to the list of contained regions.
* \details This function is used to add a region to the list of model regions
*          which are contained in the market. If the region already exists in
*          the list it is not added and a warning is printed.
* \param aRegion The name of the region to add.
*/
void Market::addRegion( const string& aRegion ) {
    // Convert the string to an atom.
    const Atom* regionID = AtomRegistry::getInstance()->findAtom( aRegion );

    /*! \invariant If the region ID is null that means that the world forgot to
    *              create a hashmap of all the regions. 
    */
    assert( regionID );
    /*! \invariant The ID of the found atom is the same as the name of the
    *              region, this ensures the lookup was correct.
    */
    assert( regionID->getID() == aRegion );

    // Check if the region ID does not already exist in the list.
    if( find( mContainedRegions.begin(), mContainedRegions.end(), regionID ) == mContainedRegions.end() ) {
        // Add the region ID.
        mContainedRegions.push_back( regionID );
    }
}

/*! \brief Get the IDs of all regions contained by this market.
* \details Return the list of contained regions implemented as a vector of
*          constant Atoms. This vector consists of the IDs of all regions within
*          this market.
* \return The IDs of all regions contained by this market.
*/
const vector<const Atom*>& Market::getContainedRegions() const {
    /*! \pre There is at least one contained region in the market. */
    assert( !mContainedRegions.empty() );
    return mContainedRegions;
}

/*! \brief Set an initial price for the market.
* \details This function checks if the price of the market is zero
*          in which case it resets the price to 1. This is done
*          because the model needs a non-zero starting price, but should not
*          overwrite read-in prices.
* \warning Prices for periods greater than zero will have their read-in prices
*          overridden by default when prices are initialized from the last
*          period unless that method is overridden.
*/
void Market::initPrice() {
    if ( price == 0 ) {
        price = 1;  
    }
}

/*! \brief Sets the price variable to the value specified.
* \details This method is used when it is necessary to set the price variable
*          to a value regardless of the type of the market. Note that all the
*          functions with "Raw" in the name have this behavior.
* \warning This function is not virtual.
* \author Josh Lurz
* \param priceIn The value to which to set the price member variable.
* \sa setPrice
* \sa setPriceToLast
*/
void Market::setRawPrice( const double priceIn ) {
    price = priceIn;
}

/*! \brief Set the price of the market based on the type.
* \details This method is used throughout the model to set a new price into a
*          market, but this is not used by the solution mechanism.
* \param priceIn The new price to set the market price to.
* \sa setRawPrice
* \sa setPriceToLast
*/
void Market::setPrice( const double priceIn ) {
    price = priceIn;
}

/*! \brief Set the market price using the price from the last period.
* \details This function is used when setting the price for a market to the
*          value from the last period. The reason setRawPrice is not used is so
*          that this method can be overridden to be a no-op. This is because for
*          CalibrationMarket the initial price is read in and should not be set
*          from the last period.
* \todo Markets never override their read-in price now so this function can
*       become non-virtual.
* \warning Use this instead of setRawPrice when setting the price to the price
*          from the last period.
* \param lastPrice Price from the last period to set the market price to.
* \sa setRawPrice
* \sa setPrice
*/
void Market::set_price_to_last_if_default( const double lastPrice ) {
    // Only initialize the price from last period's price if the price is set to
    // the default. This prevents overwriting read-in initial prices.
    if( price == 1 ){
        price = lastPrice;
    }
    // If last period price is null, reset to small number so that solver
    // has a value to start with.
    else if( price == 0 ){
        price = util::getSmallNumber();
    }
}

/*! \brief Set the market price using the price from the last period.
* \details This function is used when setting the price for a market to the
*          value from the last period. The reason setRawPrice is not used is so
*          that this method can be overridden to be a no-op. This is because for
*          CalibrationMarket the initial price is read in and should not be set
*          from the last period.
* \todo Markets never override their read-in price now so this function can
*       become non-virtual.
* \warning Use this instead of setRawPrice when setting the price to the price
*          from the last period.
* \param lastPrice Price from the last period to set the market price to.
* \sa setRawPrice
* \sa setPrice
*/
void Market::set_price_to_last( const double lastPrice ) {
    // Initialize the price from last period's price.
    // This resets all prices to last.
    if( price > 0 ){
        price = lastPrice;
    }
    // If last period price is null, reset to small number so that solver
    // has a value to start with.
    else if( price == 0 ){
        price = util::getSmallNumber();
    }
}

/*! \brief Get the market price. 
* \details This method is used to get the price out of a Market.
* \return The price for the Market.
* \sa getRawPrice
*/
double Market::getPrice() const {
    return price;
}

/*! \brief Get the raw price.
* \details This method is used to get the true value of the price variable in
*          the Market. It is often used in the solution mechanism. Note that all
*          the functions with "Raw" in the name have this behavior.
* \return The true value of the price variable.
* \sa getPrice
*/
double Market::getRawPrice() const {
    return price;
}

/*! \brief Get the stored price.
* \details This method is used to get the value of the storedPrice variable in
*          the Market. It is often used in the solution mechanism. This is used
*          when calculating the derivative of a market, so that the price can be
*          changed and the solution mechanism can determine the difference in
*          price, supply, and demand.
* \return The value of the storedPrice variable.
* \sa getPrice
*/
double Market::getStoredRawPrice() const {
    return storedPrice;
}

/*! \brief Null the demand.
* This function stores the demand and resets demand to zero. 
*/
void Market::nullDemand() {
    storedDemand = demand;
    demand = 0;
}

/*! \brief Set the Raw demand.
* \details This method is used to set the true value of the demand variable in
*          the Market. It is often used in the solution mechanism. Note that all
*          the functions with "Raw" in the name have this behavior.
* \param value The new value to set the demand variable to. 
* \sa setDemand
*/
void Market::setRawDemand( const double value ) {
    demand = value;
}

/*! \brief Add to the the Market an amount of demand in a method based on the
*          Market's type.
* \details This method is used throughout the model to add demand to a market. 
* \param demandIn The new demand to add to the current demand.
* \sa setRawDemand
*/
void Market::addToDemand( const double demandIn ) {
    demand += demandIn;
}

/*! \brief Remove an amount of demand from the raw demand.
* \details This function is used by the solution mechanism to subtract out an
*          amount of demand. This method was needed because addToDemand is
*          virtual, and this function needs to always change the raw demand. 
* \param demandIn Amount of demand to remove.
*/
void Market::removeFromRawDemand( const double demandIn ) {
    demand -= demandIn;
}

/*! \brief Get the raw demand.
* \details This method is used to get the true value of the demand variable in
*          the Market. It is often used in the solution mechanism. Note that all
*          the functions with "Raw" in the name have this behavior.
* \return The true value of the demand variable.
* \sa getDemand
*/
double Market::getRawDemand() const {
    return demand;
}

/*! \brief Get the stored demand.
* \details This method is used to get the value of the storedDemand variable in
*          the Market. It is often used in the solution mechanism. This is used
*          when calculating the derivative of a market, so that the price can be
*          changed and the solution mechanism can determine the difference in
*          price, supply, and demand.
* \return The value of the storedDemand variable.
* \sa getPrice
*/
double Market::getStoredRawDemand() const {
    return storedDemand;
}

/*! \brief Get the demand.
* \details Get the demand out of the market.
* \return Market demand.
*/
double Market::getDemand() const {
    return demand;
}

/*! \brief Null the supply.
* \details This function stores the supply and resets supply to zero. 
*/
void Market::nullSupply() {
    storedSupply = supply;
    supply = 0;
}

/*! \brief Get the raw supply.
* \details This method is used to get the true value of the supply variable in
*          the Market. It is often used in the solution mechanism. Note that all
*          the functions with "Raw" in the name have this behavior.
* \return The true value of the supply variable.
* \sa getSupply
*/
double Market::getRawSupply() const {
    return supply;
}

/*! \brief Get the storedSupply.
* \details This method is used to get the value of the storedSupply variable in
*          the Market. It is often used in the solution mechanism. This is used
*          when calculating the derivative of a market, so that the price can be
*          changed and the solution mechanism can determine the difference in
*          price, supply, and demand.
* \return The value of the storedSupply variable.
* \sa getSupply
*/
double Market::getStoredRawSupply() const {
    return storedSupply;
}

/*! \brief Set the Raw supply.
* \details This method is used to set the true value of the supply variable in
*          the Market. It is often used in the solution mechanism. Note that all
*          the functions with "Raw" in the name have this behavior.
* \param supplyIn The new value to set the supply variable to. 
* \sa setSupply
*/
void Market::setRawSupply( const double supplyIn ) {
    supply = supplyIn;
}

/*! \brief Get the supply.
* \details Get the supply out of the market.
* \return Market supply
*/
double Market::getSupply() const {
    return supply;
}

/*! \brief Add to the the Market an amount of supply in a method based on the
*          Market's type.
* \details This method is used throughout the model to add supply to a market. 
* \param supplyIn The new demand to add to the current demand.
* \sa setRawSupply
*/
void Market::addToSupply( const double supplyIn ) {
    supply += supplyIn;
}

/*! \brief Remove an amount of supply from the raw supply.
* \details This function is used by the solution mechanism to subtract out an
*          amount of supply. This method was needed because addToSupply is
*          virtual, and this function needs to always change the raw supply. 
* \param supplyIn Amount of supply to remove.
* \return void
*/
void Market::removeFromRawSupply( const double supplyIn ) {
    supply -= supplyIn;
}

/*! \brief Return the market name.
* \details This function returns the name of the market, as defined by region
*          name plus good name.
* \return The market name
*/
const string& Market::getName() const {
    return mName;
}

/*! \brief Return the market region.
* \details This method returns the region of the market. This may not be one of
*          the miniCAM regions, as a market region can contain several regions.
* \return The market region.
*/
const string& Market::getRegionName() const {
    return region;
}

/*! \brief Return the market good name.
* \details This function returns the good that the market represents. 
* \return The market good.
*/
const string& Market::getGoodName() const {
    return good;
}

/*! \brief Get the information object for this market which can then be used to
*          query for specific values.
* \details This function returns the internal IInfo object of this market which
*          represents a set of pairings of information name to value. The
*          information object is allocated in the constructor and so cannot be
*          null. This specific function returns a constant pointer to the
*          information object, so values can be queried but not added or
*          modified.
* \note This version of the function is required so that it can be called in
*       constant functions. A second version is available which returns a
*       mutable pointer.
* \return A constant pointer to the market information object.
* \author Josh Lurz
*/
const IInfo* Market::getMarketInfo() const {
    return mMarketInfo.get();
}

/*! \brief Get the information object for this market which can then be used to
*          query, add, or modify specific values.
* \details This function returns the internal IInfo object of this market which
*          represents a set of pairings of information name to value. The
*          information object is allocated in the constructor and so cannot be
*          null. This specific function returns a mutable pointer to the
*          information object, so values can be queried, added and modified.
* \note This function returns a mutable pointer to the information object so it
*       cannot be called from constant function.
* \return A mutable pointer to the market information object.
* \author Josh Lurz
*/
IInfo* Market::getMarketInfo() {
    return mMarketInfo.get();
}

/*! \brief Store the current demand, supply, and price.
* \details This function stores the current values of demand, supply and price
*          into their respective stored variables. 
*/
void Market::storeInfo() {
    storedDemand = demand;
    storedSupply = supply;
    storedPrice = price;
}

/*! \brief Restore the previous demand, supply, and price.
* \details This function sets the market's demand, supply and price to the
*          stored values of those variables. 
*/
void Market::restoreInfo() {
    demand = storedDemand;
    supply = storedSupply;
    price = storedPrice;
}

/*! \brief Store the original price.
*/
void Market::store_original_price() {
    original_price = price;
}

/*! \brief Store the original price.
*/
void Market::restore_original_price() {
    price = original_price;
}

/*! \brief Set that the market should be solved by the solution mechanism.
* \details This function sets a flag within the market telling the solution
*          mechanism whether it should solve it, given that it satifies whatever
*          conditions are set out in the shouldSolve and shouldSolveNR
*          functions.
* \param doSolve A flag representing whether or not to solve the market.
*/
void Market::setSolveMarket( const bool doSolve ) {
    solveMarket = doSolve;
}

/*! \brief Determine if a market should be solved.
* \details This function returns whether a Solver should attempt to solve this market.
* \return Whether to attempt to solve the market. 
*/
bool Market::shouldSolve() const {
    return solveMarket;
}

/*! \brief Determine if a market should be solved for Newton-Rhapson.
* \details This function returns whether or not a Newton-Rhaphon solution
*          mechanism should attempt to solve this market. This function checks
*          that supply and demand are greater than zero.
* \warning This function could cause the solution mechanism to not solve a
*          market that should because the market could be removed from set of
*          markets to solve when its supply or demand temporarily became less
*          than zero. If another market were adjusted however, supply or demand
*          could become positive again.
* \todo This is not the best design right now, this should be contained in the
*       solution mechanism.
* \return Whether or not to solve the market for Newton-Rhaphson.
*/
bool Market::shouldSolveNR() const {
   // Solves all solvable markets with the following conditions
   // including those with null demand.
   return ( solveMarket && price > 0 && supply > 0 );
}

/*! \brief Return whether a market is solved according to market type specific
*          conditions.
* \return Whether the market meets special solution criteria.
* \author Josh Lurz
*/
bool Market::meetsSpecialSolutionCriteria() const {
    // This is a normal market which should not be solved in the base period
    // unless the solve flag is set.
    return ( !solveMarket && period == 0 );
}


/*! \brief Update an output container with information from a Market for a given period.
* \param aVisitor The output container to update.
* \param aPeriod The period for which to update.
*/
void Market::accept( IVisitor* aVisitor, const int aPeriod ) const {
    aVisitor->startVisitMarket( this, aPeriod );
    aVisitor->endVisitMarket( this, aPeriod );
}

/*!
 * \brief Helper method to convert a type into a string for output.
 * \param aType Market type to convert to a string.
 * \return Market type as a string.
 * \todo Move this to a MarketUtils class.
 */
const string& Market::convert_type_to_string( const IMarketType::Type aType ) {
    // Check that the type is legal.
    assert( aType < IMarketType::END );

    // Setup a static array of the types of markets. If you add a type to the
    // IMarketType array make sure to add to this array as well and keep this in
    // the same order as the IMarketType enum.
    static const string types[] = { "Normal", "Calibration", "Inverse-Calibration",
                                    "Tax", "RES", "Subsidy", "Trial-Value",
                                    "Demand", "Price" };

    // Check that the types array is up to date.
    assert( sizeof( types ) / sizeof( types[ 0 ] ) == IMarketType::END );

    return types[ aType ];
}

/*!
 * \brief Returns if the flag which indicate that this market should be solved is set.
 * \details This method should not be used by a solver to determine if it should
 *          attempt to solve this market.  Instead it should use Market::shouldSolve().
 * \return The current value of the solveMarket flag.
 */
bool Market::isSolvable() const {
    return solveMarket;
}
