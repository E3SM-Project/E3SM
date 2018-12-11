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
* \file price_market.cpp
* \ingroup Objects
* \brief PriceMarket class source file.
* \author Steve Smith
*/

#include "util/base/include/definitions.h"
#include <string>
#include "marketplace/include/market.h"
#include "marketplace/include/price_market.h"
#include "util/base/include/xml_helper.h"
#include "containers/include/iinfo.h"

using namespace std;

//! Constructor
PriceMarket::PriceMarket( const string& goodNameIn, const string& regionNameIn, const int periodIn, Market* demandMarketIn ) :
Market( goodNameIn, regionNameIn, periodIn ),
mIsPriceNegative( false )
{
    assert( demandMarketIn );
    demandMarketPointer = demandMarketIn;
    mMarketInfo->setBoolean( "has-split-market", true );
}

//! Copy Constructor.
PriceMarket::PriceMarket( const Market& marketIn, Market* demandMarketIn ) : 
Market( marketIn ),
mIsPriceNegative( false )
{
    assert( demandMarketIn );
    demandMarketPointer = demandMarketIn;
    mMarketInfo->setBoolean( "has-split-market", true );
}

void PriceMarket::toDebugXMLDerived( ostream& out, Tabs* tabs ) const {
    XMLWriteElement( demandMarketPointer->getName(), "LinkedDemandMarket", out, tabs );
}

IMarketType::Type PriceMarket::getType() const {
    return IMarketType::PRICE;
}

void PriceMarket::initPrice() {
    Market::initPrice();
}

/*! \brief Set the price of the market based on the type.
*
* This method is used throughout the model to set a new price into a market. 
* But this is not used by the solution mechanism.
* For the price markets the price set directly to this market (priceIn) is the 
* actual price calculated by a supply sector. This is set to the demand side.
* The supply side is set here to equal the price variable (set directly by the solution mechanism)
* This sets up the solution mechanism to always give a trial value for this price.
*
* \author Josh Lurz
* \param priceIn The new price to set the market price to.
* \sa setRawPrice
* \sa setPriceToLast
*/
void PriceMarket::setPrice( const double priceIn ) {
    // If the model thinks the price should be negative store that information
    // and hide the fact from the solver.
    if( priceIn < 0 ) {
        mIsPriceNegative = true;
        demand = abs( priceIn );
    }
    else {
        mIsPriceNegative = false;
        demand = priceIn;
    }
    supply = price;
}

void PriceMarket::set_price_to_last_if_default( const double lastPrice ) {
    Market::set_price_to_last_if_default( lastPrice );
}

void PriceMarket::set_price_to_last( const double lastPrice ) {
    Market::set_price_to_last( lastPrice );
}

double PriceMarket::getPrice() const {
    // This is the model facing price so if it should be negative make it so now.
    return mIsPriceNegative ? -price : price; 
}

void PriceMarket::addToDemand( const double demandIn ) {
    demandMarketPointer->addToDemand( demandIn );
}

double PriceMarket::getDemand() const {
    return demandMarketPointer->getDemand();
}

void PriceMarket::nullSupply() {
    Market::nullSupply();
}

double PriceMarket::getSupply() const {
    return demandMarketPointer->getSupply();
}

void PriceMarket::addToSupply( const double supplyIn ) {
    demandMarketPointer->addToSupply( supplyIn );
}

bool PriceMarket::meetsSpecialSolutionCriteria() const {
    return Market::meetsSpecialSolutionCriteria();
}

bool PriceMarket::shouldSolve() const {
    return Market::shouldSolve();
}

bool PriceMarket::shouldSolveNR() const {
    return Market::shouldSolveNR();
}
