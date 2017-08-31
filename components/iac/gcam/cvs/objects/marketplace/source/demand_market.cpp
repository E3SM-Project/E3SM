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
* \file demand_market.cpp
* \ingroup Objects
* \brief DemandMarket class source file.
* \author Steve Smith
*/

#include "util/base/include/definitions.h"
#include <string>
#include "marketplace/include/demand_market.h"
#include "util/base/include/xml_helper.h"

using namespace std;

//! Constructor
DemandMarket::DemandMarket( const string& goodNameIn, const string& regionNameIn, const int periodIn ) :
Market( goodNameIn, regionNameIn, periodIn ) {
   demMktSupply = 0;
}

void DemandMarket::toDebugXMLDerived( ostream& out, Tabs* tabs ) const {
   XMLWriteElement( demMktSupply, "DemandMarketSupply", out, tabs );
}

IMarketType::Type DemandMarket::getType() const {
    return IMarketType::DEMAND;
}

void DemandMarket::initPrice() {
    Market::initPrice();
}

void DemandMarket::setPrice( const double priceIn ) {
    Market::setPrice( priceIn );
}

void DemandMarket::set_price_to_last_if_default( const double lastPrice ) {
    Market::set_price_to_last_if_default( lastPrice );
}

void DemandMarket::set_price_to_last( const double lastPrice ) {
    Market::set_price_to_last( lastPrice );
}

double DemandMarket::getPrice() const {
    return Market::getPrice();
}

void DemandMarket::addToDemand( const double demandIn ) {
    Market::addToDemand( demandIn );
}

double DemandMarket::getDemand() const {
   return price;
}

void DemandMarket::nullSupply() {
   supply = 0;
   demMktSupply = 0;
}

double DemandMarket::getSupply() const {
    return Market::getSupply();
}

void DemandMarket::addToSupply( const double supplyIn ) {
   supply = price; 
   demMktSupply += supplyIn; // Store Raw supply value to check later
}

bool DemandMarket::meetsSpecialSolutionCriteria() const {
    return Market::meetsSpecialSolutionCriteria();
}

bool DemandMarket::shouldSolve() const {
    return Market::shouldSolve();
}

bool DemandMarket::shouldSolveNR() const {
    return Market::shouldSolveNR();
}
