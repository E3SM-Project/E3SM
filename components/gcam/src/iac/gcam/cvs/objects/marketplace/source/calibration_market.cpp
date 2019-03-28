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
* \file calibration_market.cpp
* \ingroup Objects
* \brief CalibrationMarket class source file.
* \author Josh Lurz
*/

#include "util/base/include/definitions.h"
#include <string>
#include "marketplace/include/calibration_market.h"
#include "util/base/include/util.h"

using namespace std;

//! Constructor
CalibrationMarket::CalibrationMarket( const string& goodNameIn, const string& regionNameIn, const int periodIn ) :
Market( goodNameIn, regionNameIn, periodIn ) {
}

//! Destructor
CalibrationMarket::~CalibrationMarket() {
}

void CalibrationMarket::toDebugXMLDerived( ostream& out, Tabs* tabs ) const {
}

IMarketType::Type CalibrationMarket::getType() const {
    return IMarketType::CALIBRATION;
}

void CalibrationMarket::initPrice() {
    Market::initPrice();
}

void CalibrationMarket::setPrice( const double priceIn ) {
    Market::setPrice( priceIn );
}

void CalibrationMarket::set_price_to_last_if_default( const double lastPrice ) {
   // Do nothing, as in a calibration market the initial values for each period are set from the XML.
}

void CalibrationMarket::set_price_to_last( const double lastPrice ) {
   // Do nothing, as in a calibration market the initial values for each period are set from the XML.
}

double CalibrationMarket::getPrice() const {
    return Market::getPrice();
}

void CalibrationMarket::addToDemand( const double demandIn ) {
    Market::addToDemand( demandIn );
}

double CalibrationMarket::getDemand() const {
    return Market::getDemand();
}

// Do nothing, as constraints should not be cleared.
void CalibrationMarket::nullDemand() {
}

void CalibrationMarket::nullSupply() {
    Market::nullSupply();
}

double CalibrationMarket::getSupply() const {
    return Market::getSupply();
}

void CalibrationMarket::addToSupply( const double supplyIn ) {
    Market::addToSupply( supplyIn );
}

bool CalibrationMarket::shouldSolve() const {
    return Market::shouldSolve();
}

bool CalibrationMarket::shouldSolveNR() const {
    return Market::shouldSolveNR();
}

bool CalibrationMarket::meetsSpecialSolutionCriteria() const {
    // Check if the market was set to solve without a constraint.
    return ( demand < util::getSmallNumber() );
}

void CalibrationMarket::removeFromRawDemand( const double demandIn ) {
    // need to make sure that the demand will never get wiped out
    // because it will not get set again for a calibration market
}
