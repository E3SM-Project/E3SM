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
* \file inverse_calibration_market.cpp
* \ingroup Objects
* \brief InverseCalibrationMarket class source file.
* \author Josh Lurz
*/

#include "util/base/include/definitions.h"
#include <string>
#include "marketplace/include/inverse_calibration_market.h"
#include "util/base/include/util.h"

using namespace std;

//! Constructor
InverseCalibrationMarket::InverseCalibrationMarket( const string& goodNameIn, const string& regionNameIn, const int periodIn ) :
Market( goodNameIn, regionNameIn, periodIn ) {
}

//! Destructor
InverseCalibrationMarket::~InverseCalibrationMarket() {
}

void InverseCalibrationMarket::toDebugXMLDerived( ostream& out, Tabs* tabs ) const {
}

IMarketType::Type InverseCalibrationMarket::getType() const {
    return IMarketType::INVERSE_CALIBRATION;
}

void InverseCalibrationMarket::initPrice() {
    Market::initPrice();
}

void InverseCalibrationMarket::setPrice( const double priceIn ) {
    Market::setPrice( priceIn );
}

void InverseCalibrationMarket::set_price_to_last_if_default( const double lastPrice ) {
   // Do nothing, as in a calibration market the initial values for each period
   // are set from the XML.
}

void InverseCalibrationMarket::set_price_to_last( const double lastPrice ) {
   // Do nothing, as in a calibration market the initial values for each period
   // are set from the XML.
}

double InverseCalibrationMarket::getPrice() const {
    return Market::getPrice();
}

void InverseCalibrationMarket::addToDemand( const double demandIn ) {
    Market::addToDemand( demandIn );
}

double InverseCalibrationMarket::getDemand() const {
    return Market::getDemand();
}

void InverseCalibrationMarket::nullDemand() {
    // Differs from the CalibrationMarket because demand should be cleared each
    // iteration. In the CalibrationMarket, the constraint is stored in the demand
    // variable and so should not be cleared.
    Market::nullDemand();
}

void InverseCalibrationMarket::nullSupply() {
    // Differs from the CalibrationMarket because supply should not be cleared
    // each iteration. In the CalibrationMarket, the constraint is stored in the
    // demand variable and so supply is cleared each iteration.
}

double InverseCalibrationMarket::getSupply() const {
    return Market::getSupply();
}

void InverseCalibrationMarket::addToSupply( const double supplyIn ) {
    Market::addToSupply( supplyIn );
}

bool InverseCalibrationMarket::shouldSolve() const {
    return Market::shouldSolve();
}

bool InverseCalibrationMarket::shouldSolveNR() const {
    return Market::shouldSolveNR();
}

bool InverseCalibrationMarket::meetsSpecialSolutionCriteria() const {
    // Check if the market was set to solve without a constraint. Differs from
    // the CalibrationMarket because the supply variable is checked instead of
    // the demand variable because it is the constraint.
    return ( supply < util::getSmallNumber() );
}
