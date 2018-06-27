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
* \file supply_demand_curve.cpp
* \ingroup Objects
* \brief SupplyDemandCurve class source file.
* \author Josh Lurz
*/

#include "util/base/include/definitions.h"
#include <iostream>
#include <cassert>
#include <algorithm>
#include <cmath>
#include "util/base/include/supply_demand_curve.h"
#include "marketplace/include/market.h"
#include "util/logger/include/ilogger.h"
#include "containers/include/world.h"
#include "marketplace/include/marketplace.h"

using namespace std;

//! Constructor
SupplyDemandCurve::SupplyDemandCurve( Market* marketIn ) {
    // Make sure the pointer is non-null.
    assert( marketIn );
    market = marketIn;
}

//! Destructor
SupplyDemandCurve::~SupplyDemandCurve() {
    for ( vector<SupplyDemandPoint*>::iterator i = points.begin(); i != points.end(); i++ ) {
        delete *i;
    }
}

/*! \brief Calculate given number of supply and demand points.
*
* This function first determines a series of price ratios to use to determine the prices to 
* create SupplyDemandPoints for. It then saves the original marketplace information, and perturbs the price
* as specified by the price ratios. Using this new perturbed price, it call World::Calc to determine supply and 
* demand for the market. It saves that point for printing later, and continues to perform this process for each price.
* Finally it restores the original market information.
*
* \param numPoints The number of points to calculate.
* \param world The World object to use for World::calc
* \param marketplace The marketplace to use to store and restore information.
* \param period The period to perform the calculations on.
* \todo It would be good if this used a similar logic to SolverLibrary::derivatives to save time. 
* \todo Un-hardcode the prices. 
*/

void SupplyDemandCurve::calculatePoints( const int numPoints, World* world, Marketplace* marketplace, const int period ) {
    vector<double> priceMults;

    // Determine price ratios.
    const int middle = static_cast<int>( floor( double( numPoints ) / double( 2 ) ) );

    for( int pointNumber = 0; pointNumber < numPoints; pointNumber++ ) {

        if( pointNumber < middle ) {
            priceMults.push_back( 1 - double( 1 ) / double( middle - abs( middle - pointNumber ) + 2 ) );
        }
        else if( pointNumber == middle ) {
            priceMults.push_back( 1 );
        }
        else {
            priceMults.push_back( 1 + double( 1 ) / double( middle - abs( middle - pointNumber ) + 2 ) );
        }
    }

    // Store the market info.
    marketplace->storeinfo( period );

    // get the base price of the market.
    // const double basePrice = market->getRawPrice();

    // Save the original point as price 1.
    // ( *iter )->createSDPoint();

    // iterate through the points and determine supply and demand.
    for ( int pointNumber2 = 0; pointNumber2 < numPoints; pointNumber2++ ) {

        // clear demands and supplies.
        marketplace->nullSuppliesAndDemands( period );

        // set the price to the current point.
        // market->setRawPrice( priceMults[ pointNumber2 ] * basePrice );
        market->setRawPrice( pointNumber2 * ( double( 10 ) / double( numPoints - 1 ) ) );

        // calculate the world.
        world->calc( period );

        points.push_back( new SupplyDemandPoint(  market->getRawPrice(), market->getRawDemand(), market->getRawSupply() ) );

    } // Completed iterating through all price points.

    marketplace->restoreinfo( period );
    marketplace->nullSuppliesAndDemands( period );

    // Call world.calc a final time to restore information for summary.
    world->calc( period );
}

/*! \brief Print the supply demand curve.
*
* This function prints the vector of SupplyDemandPoints created during the calculatePoints function.
* It creates a copy of points and sorts them before printing them. This was done so that the printing function
* could remain constant. The points are printed to the Logger passed as an argument.
*
* \param sdLog Logger to print the points to.
*/
void SupplyDemandCurve::print( ILogger& aSDLog ) const {
    aSDLog << "Supply and Demand curves for: " << market->getName() << endl;
    aSDLog << "Price,Demand,Supply," << endl;

    // Create a copy of the points vector so that we can sort it while keeping the print function constant.
    // Since the vector contains pointers to SupplyDemandPoints, this is relatively inexpensive.
    vector<SupplyDemandPoint*> pointsCopy( points );

    // Sort the SupplyDemandPoint object pointers in the pointsCopy vector by increasing price by using the LesserPrice binary operator. 
    sort( pointsCopy.begin(), pointsCopy.end(), SupplyDemandPoint::LesserPrice() );
    for ( vector<SupplyDemandPoint*>::const_iterator i = pointsCopy.begin(); i != pointsCopy.end(); i++ ) {
        ( *i )->print( aSDLog );
    }

    aSDLog << endl;
}

//! Constructor
SupplyDemandCurve::SupplyDemandPoint::SupplyDemandPoint( const double priceIn, const double demandIn, const double supplyIn ) 
: price( priceIn ), demand( demandIn ), supply( supplyIn ){
}

/*! \brief Get the price.
*
* Return the price contained in the point. This is needed for sorting.
*
* \return The price saved within the point.
*/
double SupplyDemandCurve::SupplyDemandPoint::getPrice() const {
    return price;
}

/*! \brief Print the point in a csv format.
*
* Print the point in a csv format to the specified logger.
*
* \param sdLog The Logger to print to.
*/
void SupplyDemandCurve::SupplyDemandPoint::print( ILogger& aSDLog ) const {
    aSDLog << price << "," << demand << "," << supply << endl;
}

