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
* \file aproduction_function.cpp
* \ingroup Objects
* \brief The AProductionFunction class source file.
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

#include "functions/include/aproduction_function.h"
#include "functions/include/iinput.h"
#include "containers/include/scenario.h"
#include "marketplace/include/marketplace.h"
#include "util/base/include/util.h"
#include "functions/include/function_utils.h"
#include "util/logger/include/ilogger.h"

using namespace std;

extern Scenario* scenario;

/*! \brief Calculate Costs
* \pre Demand currencies for all inputs must be known.
* \return The total cost of production.
* \todo Why is alpha zero not used in this function?
*/
double AProductionFunction::calcCosts( const InputSet& input,
                                       const string& regionName,
                                       const double aAlphaZero,
                                       int period ) const
{
	double totalCost = 0;
	for( unsigned int i = 0; i < input.size(); ++i ) {
		// capital input name should be changed to OtherValueAdded
		if( !input[i]->hasTypeFlag( IInput::CAPITAL ) ) {
			totalCost += input[i]->getPhysicalDemand( period )
				         * input[i]->getPricePaid( regionName, period );
		}
        else {
            totalCost += input[i]->calcTaxes( regionName, 0, 0, period );
        }
	}
	return totalCost;
}

//! Calculate the profits using the production scale factor.
double AProductionFunction::calcProfits( InputSet& input, const string& regionName, 
                                         const string& sectorName, const double aShutdownCoef,
                                         int period, 
                                         double capitalStock, double alphaZero, double sigma ) const
{

    double profits = calcUnscaledProfits( input, regionName, sectorName, period,
                                          capitalStock, alphaZero, sigma );
    // Scale the profits by the scaler which determines how much of the vintage
    // was shutdown.
    return aShutdownCoef * profits;
}
