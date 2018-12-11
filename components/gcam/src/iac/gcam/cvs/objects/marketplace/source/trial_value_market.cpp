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
* \file trial_value_market.cpp
* \ingroup Objects
* \brief The TrialValueMarket class header file.
* \author Steve Smith
*/

#include "util/base/include/definitions.h"
#include "util/base/include/util.h"
#include <string>
#include "marketplace/include/trial_value_market.h"

using namespace std;

//! Constructor
TrialValueMarket::TrialValueMarket( const string& goodNameIn, const string& regionNameIn, const int periodIn ) :
Market( goodNameIn, regionNameIn, periodIn )
{   
    // Initialize to 0.001. Use of previous getSmallNumber() is too small and 
    // takes longer to solve.
    price = 0.001;
}

void TrialValueMarket::toDebugXMLDerived( ostream& out, Tabs* tabs ) const {
}

IMarketType::Type TrialValueMarket::getType() const {
    return IMarketType::TRIAL_VALUE;
}

void TrialValueMarket::initPrice() {
    // Note zero may be a valid price for a trial value.
}

void TrialValueMarket::setPrice( const double priceIn ) {
    Market::setPrice( priceIn );
}

void TrialValueMarket::set_price_to_last_if_default( const double lastPrice ) {
   //Market::set_price_to_last_if_default( lastPrice );
    // Only initialize the price from last period's price if the price is set to
    // the default. This prevents overwriting read-in initial prices.
    if( price == 0.001 ){
        price = lastPrice;
    }
    // Note zero may be a valid price for a trial value.
}

void TrialValueMarket::set_price_to_last( const double lastPrice ) {
    // Initialize the price from last period's price.
    // This resets all prices to last.
    if( price > 0 ){
        price = lastPrice;
    }
    // Note zero may be a valid price for a trial value.
}

double TrialValueMarket::getPrice() const {
    return Market::getPrice();
}

/*! \brief Add to the the Market an amount of demand in a method based on the Market's type.
* This is the only method that is different for the trial market type. 
* Here is where the price variable is copied to supply, thereby setting up the solution mechanism
* to solve for the trial value of this quantity.
*
* \author Steve Smith
*
* \param demandIn The new demand to add to the current demand.
* \sa setRawDemand
*/
void TrialValueMarket::addToDemand( const double demandIn ) {
    Market::addToDemand( demandIn );
    supply = price;
}

double TrialValueMarket::getDemand() const {
    return Market::getDemand();
}

void TrialValueMarket::nullSupply() {
   Market::nullSupply();
}

double TrialValueMarket::getSupply() const {
    return Market::getSupply();
}

void TrialValueMarket::addToSupply( const double supplyIn ) {
    Market::addToSupply( supplyIn );
}

bool TrialValueMarket::meetsSpecialSolutionCriteria() const {
    // Trial value markets must be solved in all periods including the base
    // period.
    return false;
}

bool TrialValueMarket::shouldSolve() const {
    bool doSolveMarket = false;
    // Check if this market is a type that is solved.
    if ( solveMarket ) {
        // Solve all solvable markets.
        doSolveMarket = true;
    }
    return doSolveMarket;
}

bool TrialValueMarket::shouldSolveNR() const {
    bool doSolveMarket = false;
    // Check if this market is a type that is solved.
    if ( solveMarket ) {
        // Solve all solvable markets including those with null demand.
        doSolveMarket = true;
    }
    return doSolveMarket;
}
