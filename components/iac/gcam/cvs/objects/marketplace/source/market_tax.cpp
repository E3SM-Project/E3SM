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
* \file market_tax.cpp
* \ingroup Objects
* \brief MarketTax class source file. (Originally called GHGMarket)
* \author Sonny Kim
*/

#include "util/base/include/definitions.h"
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "marketplace/include/market_tax.h"
#include "util/base/include/util.h"

using namespace std;

///! Constructor
MarketTax::MarketTax( const string& goodNameIn, const string& regionNameIn, const int periodIn ) :
Market( goodNameIn, regionNameIn, periodIn ) {
}

void MarketTax::toDebugXMLDerived( ostream& out, Tabs* tabs ) const {
}

IMarketType::Type MarketTax::getType() const {
    return IMarketType::TAX;
}

/* \brief Initialize the MarketTax price.
* \details This method initializes the price of the tax market to a random number
* or null if not a solved market.
* \author Josh Lurz, Sonny Kim
*/
void MarketTax::initPrice() {
    const double MIN_PRICE = 5;
    // If price is near zero it needs to be initialized.
    if( price < util::getSmallNumber() ){
        // If this market should be solved price should be initialized to a 
        // random number between MIN_PRICE and (1 + MIN_PRICE)
        if( solveMarket ){
            srand( (unsigned)time( NULL ) );
            price = ((double) rand() / (double) RAND_MAX) + MIN_PRICE;
            price = util::getSmallNumber();
        }
        // The market will not be solved so it is set to null. 
        else {
            price = 0;
        }
    }
}

void MarketTax::setPrice( const double priceIn ) {
    Market::setPrice( priceIn );
}

/* \brief Initialize the MarketTax price from last period's price.
* \details This method checks if the lastPrice was 0. This would mean that last period's constraint was 
* 0. If it is, price is set to a random as this is the initial constrained period.
* Otherwise price is set to the previous period's price as is done in the normal market.
* \param lastPrice Previous period's price. This should have already been set in store to last!!
* \author Josh Lurz, Sonny Kim
*/
void MarketTax::set_price_to_last_if_default( const double lastPrice ) {
    const double MIN_PRICE = 5;
    // If the price is zero and the solve flag is set so a constraint exists. 
    if( price <= util::getSmallNumber() && solveMarket ){
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

/* \brief Initialize the MarketTax price from last period's price.
* \details This method checks if the lastPrice was 0. This would mean that last period's constraint was 
* 0. If it is, price is set to a random as this is the initial constrained period.
* Otherwise price is set to the previous period's price as is done in the normal market.
* \param lastPrice Previous period's price. This should have already been set in store to last!!
* \author Josh Lurz, Sonny Kim
*/
void MarketTax::set_price_to_last( const double lastPrice ) {
    const double MIN_PRICE = 5;
    // If the price is zero and the solve flag is set so a constraint exists. 
    if( price <= util::getSmallNumber() && solveMarket ){
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

double MarketTax::getPrice() const {
    return Market::getPrice();
}

void MarketTax::addToDemand( const double demandIn ) {
    Market::addToDemand( demandIn );
}

double MarketTax::getDemand() const {
    return Market::getDemand();
}

//! The supply in MarketTax is the constraint, it should not be removed by calls to nullSupply
void MarketTax::nullSupply() {
}

double MarketTax::getSupply() const {
    return Market::getSupply();
}

void MarketTax::addToSupply( const double supplyIn ) {
    Market::addToSupply( supplyIn );
}

/* \brief This method determines whether to solve a MarketTax with the solution mechanism.
* \details This method returns false for solving the tax market if constraints are not binding
* when the price is essentially null. Otherwise, it returns the default solveMarket boolean
* that is set to true for a valid market by the marketplace object.  The check for binding 
* constraints prevents wasted efforts in the solution mechanism, such as in finding bracket 
* intervals and running bisection.
* \return Whether to solve the market.
* \author Sonny Kim
*/
bool MarketTax::shouldSolve() const {
    bool doSolveMarket = false;
    // Check if this market is a type that is solved (i.e resource, policy, etc.)
    // Note: secondary markets are not solved in miniCAM
    if ( solveMarket) {
        // if constraint does exist then solve
        if( supply > util::getSmallNumber() ) {
            doSolveMarket = true;
            // if constraint exists but not binding with null price
            // don't solve
            if( (price <= util::getSmallNumber()) && (demand <= supply)){
                doSolveMarket = false;
            }
        }
    }
    return doSolveMarket;
}

/* \brief This method determines whether to solve a MarketTax with the NR solution mechanism.
* \details This method only returns that the NR solution mechanism should attempt to solve the market
* if a constraint exists and is not binding with a null price
* \return Whether to solve the market in NR.
* \author Sonny Kim
*/
bool MarketTax::shouldSolveNR() const {
    bool doSolveMarket = false;
    /* Old code: in case tax market should be included in the NR solver.
    // Check if this market is a type that is solved (i.e resource, policy, etc.)
    // Note: secondary markets are not solved in miniCAM
    if ( solveMarket) {
        // if constraint does exist then solve
        if( supply > util::getSmallNumber() ) {
            doSolveMarket = true;
            // if constraint exists but not binding with null price
            // don't solve
            if( (price <= util::getSmallNumber()) && (demand <= supply)){
                doSolveMarket = false;
            }
        }
    }
    */
    // Do not include tax market in the NR Solver.
    return doSolveMarket;
}

/* \brief Check whether the market meets market specific solution criteris.
* \details This method returns true if market is solved or if price is null 
* and constraint is not binding.
* \return Whether special solution criteria is met.
* \author Sonny Kim
*/
bool MarketTax::meetsSpecialSolutionCriteria() const {
    // If there is no constraint, this market is solved.
    if( !solveMarket ){
        return true;
    }

    // If price is zero, demand cannot be driven any higher.
    // The constraint is not binding (greater than the demand), so this market is solved.
    if( ( price <= util::getSmallNumber() ) && ( supply >= demand ) ){
        return true;
    }
    return false;
}
