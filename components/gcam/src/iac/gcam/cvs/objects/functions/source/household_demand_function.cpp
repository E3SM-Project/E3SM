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
* \file household_demand_function.cpp
* \ingroup Objects
* \brief The HouseholdDemandFunction class source file.
*
*  Detailed description.
*
* \author Pralit Patel
* \author Sonny Kim
* \author Josh Lurz
* \date $Date: 2005/06/01 21:23:59 $
* \version $Revision: 1.2 $
*/

#include "util/base/include/definitions.h"
#include <string>
#include <vector>
#include <cmath>
#include <cassert>

#include "functions/include/household_demand_function.h"
#include "functions/include/iinput.h"
#include "containers/include/scenario.h"
#include "marketplace/include/marketplace.h"
#include "util/base/include/util.h"
#include "functions/include/function_utils.h"
#include "util/logger/include/ilogger.h"

using namespace std;

extern Scenario* scenario;

//! Calculate Demand
double HouseholdDemandFunction::calcDemand( InputSet& input, double personalIncome, 
                                            const string& regionName, const string& sectorName,
                                            const double aShutdownCoef,
											int period, double capitalStock, double alphaZero, 
											double sigma, double IBT ) const 
{
    // Find the numeraire to determine the price paid.
    const IInput* numInput = FunctionUtils::getNumeraireInput( input );
    assert( numInput );
    const double pricePaidNumeraire = numInput->getPricePaid( regionName, period );
    assert( pricePaidNumeraire > 0 );

    double totalDemand = 0; // total demand used for scaling
	for ( unsigned int i = 0; i < input.size(); ++i ) {
        if( !input[ i ]->hasTypeFlag( IInput::FACTOR ) ){
            assert( input[i]->getPricePaid( regionName, period ) > 0 );
            assert( input[i]->getCoefficient( period ) > 0 );
            assert( input[i]->getIncomeElasticity() > 0 );

            double demand = input[i]->getCoefficient( period ) * 
                pow( personalIncome / pricePaidNumeraire, input[i]->getIncomeElasticity() ) *
				pow( input[i]->getPricePaid( regionName, period ) / pricePaidNumeraire, input[i]->getPriceElasticity() );

            assert( util::isValidNumber( demand ) );
			input[i]->setPhysicalDemand( demand, regionName, period );
			totalDemand += demand * input[i]->getPricePaid( regionName, period );
		}
	}
	return totalDemand;
}

double HouseholdDemandFunction::calcCoefficient( InputSet& input, double consumption,
                                                 const string& regionName, const string& sectorName,
                                                 int period, double sigma, double IBT,
                                                 double capitalStock ) const
{
	for ( unsigned int i = 0; i < input.size(); ++i ) {
		if( !input[ i ]->hasTypeFlag( IInput::FACTOR ) ){
            // if we use price paid for calcDemand probably need to use it here also
            // really currency demand
			double tempCoefficient = input[i]->getPhysicalDemand( period ) /
                                     pow( consumption, input[i]->getIncomeElasticity() ) *
                                     pow( input[i]->getPrice( regionName, period ),
                                          input[i]->getPriceElasticity() );
			input[i]->setCoefficient( tempCoefficient, period );
		}
	}
	return 0; // return null for CES AlphaZero only
}
