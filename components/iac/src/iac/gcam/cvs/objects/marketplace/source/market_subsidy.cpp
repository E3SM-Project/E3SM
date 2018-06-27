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
* \file market_subsidy.cpp
* \ingroup Objects
* \brief MarketSubsidy class source file. Originally called MarketPortfolioStandard
* \author Sonny Kim
*/

#include "util/base/include/definitions.h"
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "marketplace/include/market_subsidy.h"
#include "util/base/include/util.h"

using namespace std;

///! Constructor
MarketSubsidy::MarketSubsidy( const string& goodNameIn, const string& regionNameIn, const int periodIn ) :
Market( goodNameIn, regionNameIn, periodIn ) {
}

void MarketSubsidy::toDebugXMLDerived( ostream& out, Tabs* tabs ) const {
}

IMarketType::Type MarketSubsidy::getType() const {
    return IMarketType::SUBSIDY;
}

/* \brief Initialize the MarketSubsidy price.
* \details This method checks first if the price has already been initialized.
* If it has been initialized, the price is left unchanged. Otherwise the method checks the value of the constraint,
* or supply, to determine how to initialize the price. If supply is 0, price is set to 0. If the supply 
* is greater than 0, the price is initialized to a random number.
* \author Josh Lurz
* \todo Some of these changes for the carbon market might be beneficial to the general market, or end up being the same.
*/
void MarketSubsidy::initPrice() {
    const double MIN_PRICE = 0.5;
    // If price is near zero it needs to be initialized.
    if( price < util::getSmallNumber() ){
        // If this market should be solved price should be initialized to a 
        // random number between MIN_PRICE and (1 + MIN_PRICE)
        if( solveMarket ){
            srand( (unsigned)time( NULL ) );
            price = ((double) rand() / (double) RAND_MAX) + MIN_PRICE;
        }
        // The market will not be solved so it should be zero'd out. 
        else {
            price = 0;
        }
    }
}

void MarketSubsidy::setPrice( const double priceIn ) {
    Market::setPrice( priceIn );
}

/* \brief Initialize the MarketSubsidy price from last period's price.
* \details This method first checks if the lastPrice was 0. This would mean that last period's constraint was 
* 0. If it is, then it checks if the constraint in the current period is greater than 0. In this case price is 
* set to a random price as this is the initial constrained period. Otherwise price is set to the previous 
* period's price as is done in the normal market.
* \param lastPrice Previous period's price. This should have already been set in store to last!!
* \author Josh Lurz
*/
void MarketSubsidy::set_price_to_last_if_default( const double lastPrice ) {
    const double MIN_PRICE = 0.5;
    // If the price is zero and the solve flag is set so a constraint exists. 
    if( price < util::getSmallNumber() && solveMarket ){
        // If the last price is 0, we should set the price to a random number.
        // New price is between MIN_PRICE and (1 + MIN_PRICE)
        if( lastPrice < util::getSmallNumber() ){
            srand( (unsigned)time( NULL ) );
            price = ((double) rand() / (double) RAND_MAX) + MIN_PRICE;
        }
        // Otherwise set the price to the previous period's price.
        else {
            price = lastPrice;
        }
    }
    // There is no else here because we do not want to override prices in the case of a fixed tax.
}

void MarketSubsidy::set_price_to_last( const double lastPrice ) {
    Market::set_price_to_last( lastPrice );
}

double MarketSubsidy::getPrice() const {
    return Market::getPrice();
}

void MarketSubsidy::addToDemand( const double demandIn ) {
    Market::addToDemand( demandIn );
}

double MarketSubsidy::getDemand() const {
    return Market::getDemand();
}

//! The demand in MarketSubsidy is the constraint,
//! it should not be removed by calls to nullDemand
void MarketSubsidy::nullDemand() {
    // Virtual function to override Market::nullDemand() and 
    // clearing of demand which is the constraint.
}

void MarketSubsidy::nullSupply() {
   Market::nullSupply();
}

double MarketSubsidy::getSupply() const {
    return Market::getSupply();
}

void MarketSubsidy::addToSupply( const double supplyIn ) {
    Market::addToSupply( supplyIn );
}

/* \brief This method determines whether to solve a MarketSubsidy with the solution mechanism.
* \details This method only returns that the solution mechanism should attempt to solve the market
* if the demand (constraint) is positive. This prevents the solution mechanism from trying to solve 
* unconstrained time periods.
* \return Whether to solve the market.
* \author Sonny Kim
*/
bool MarketSubsidy::shouldSolve() const {
    bool doSolveMarket = false;
    // Check if this market is a type is solved (i.e. resource, policy, etc.)
    // Note: secondary markets are not solved in the MiniCAM
    if ( solveMarket ){
        // if constraint does exist then solve
        if( demand > util::getSmallNumber() ){
            doSolveMarket = true;
            // if constraint exists but not binding with null price then
            // don't solve
            if ( (price <= util::getSmallNumber()) && (supply >= demand) ){
                doSolveMarket = false; 
            }
        }
    }
    return doSolveMarket;
}

/* \brief This method determines whether to solve a MarketSubsidy with the NR solution mechanism.
* \details This method only returns that the NR solution mechanism should attempt to solve the market
* if the demand (constraint) is positive, but the constraint isn't already met with a null price. 
* \todo This needs a more detailed explanation, but the logic is convoluded.
* \return Whether to solve the market in NR.
* \author Sonny Kim
*/
bool MarketSubsidy::shouldSolveNR() const {
    bool doSolveMarket = false;
    // Check if this market is a type is solved (i.e. resource, policy, etc.)
    // Note: secondary markets are not solved in the MiniCAM
    /* Old code: in case subsidey market should be included in the NR solver.
    if ( solveMarket ){
        // if constraint does exist then solve
        if( demand > util::getSmallNumber() ){
            doSolveMarket = true;
            // if constraint exists but not binding with null price then
            // don't solve
            if ( (price <= util::getSmallNumber()) && (supply >= demand) ){
                doSolveMarket = false; 
            }
        }
    }
    */
    // Do not include subsidy market in the NR Solver.
    return doSolveMarket;
}

//! Check whether the market meets market specific solution criteris.
bool MarketSubsidy::meetsSpecialSolutionCriteria() const {
    // If there is no constraint, this market is solved.
    if( !solveMarket ){
        return true;
    }

    // Treat the subsidy market as a normal supply and demand market
    // where supplies respond to prices.
    // InputSubsidy class converts the price into a subsidy.

    // If the subsidy is null and the constraint has been met,
    // this market is solved.
    if( ( price <= util::getSmallNumber() ) && ( supply >= demand ) ){
        return true;
    }
    return false;
}
