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
* \file government_demand_function.cpp
* \ingroup Objects
* \brief The GovernmentDemandFunction class source file.
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

#include "functions/include/government_demand_function.h"
#include "functions/include/iinput.h"
#include "containers/include/scenario.h"
#include "marketplace/include/marketplace.h"
#include "util/base/include/util.h"
#include "functions/include/function_utils.h"
#include "util/logger/include/ilogger.h"

using namespace std;

extern Scenario* scenario;

double GovernmentDemandFunction::calcDemand( InputSet& input, double consumption,
										     const string& regionName, 
											 const string& sectorName,
											 const double aShutdownCoef,
											 int period, double capitalStock, 
											 double alphaZero, double sigma, double IBT ) const 
{
    double totalUnscaledDemand = 0; // total demand used for scaling
    for( unsigned int i = 0; i < input.size(); ++i ){
        if ( !input[ i ]->hasTypeFlag( IInput::CAPITAL ) ) {
            totalUnscaledDemand += input[ i ]->getCoefficient( period )
				                   * input[ i ]->getPricePaid( regionName, period );
        }
        // input is capital
        else {
            // use numeraire price paid
            const IInput* numInput = FunctionUtils::getNumeraireInput( input );
            assert( numInput );
            const double pricePaidNumeraire = numInput->getPricePaid( regionName, period ); 
            assert( pricePaidNumeraire > 0 );
            totalUnscaledDemand += input[ i ]->getCoefficient( period ) * pricePaidNumeraire;
        }
    }
    double rho = FunctionUtils::getRho( sigma );
	double govtPrefScaler = consumption * pow( totalUnscaledDemand, sigma*(rho-1) );

	double totalDemand = 0; // total final demand
	for( unsigned int i = 0; i < input.size(); ++i ){
		double demand = 0;
		if ( !input[ i ]->hasTypeFlag( IInput::CAPITAL ) ) {
			demand = input[ i ]->getCoefficient( period ) * govtPrefScaler;
			totalDemand += demand;
			input[ i ]->setCurrencyDemand( demand, regionName, period );
		}
	}
	return totalDemand;
}

double GovernmentDemandFunction::calcCoefficient( InputSet& input, double consumption, 
												  const string& regionName, const string& sectorName, 
												  int period, double sigma, double IBT, 
											      double capitalStock ) const 
{
	// get government's demand for capital
	double capitalDemand = 0;
	for( unsigned int i = 0; i < input.size(); ++i ){
		if ( input[ i ]->hasTypeFlag( IInput::CAPITAL ) ) {
			capitalDemand = input[ i ]->getCurrencyDemand( period );
            break;
		}
	}
	// total used for government demand does not include capital
	double tempTotal = FunctionUtils::getCurrencyDemandSum( input, period ) - capitalDemand;
	// government capital demand is calculated separatly
	// this needs to be corrected, legacy SGM probably could not handle
	for( unsigned int i = 0; i<input.size(); ++i ){
		if ( !input[ i ]->hasTypeFlag( IInput::CAPITAL ) ) {
            // not sure which price
			double tempCoefficient = input[ i ]->getCurrencyDemand( period ) 
			                         / tempTotal 
                                     / input[ i ]->getPrice( regionName, period );
			input[ i ]->setCoefficient( tempCoefficient, period );
		}
	}
    // Is this complete or not? -JPL
	//  to be continued ..

	return 0; // fake, for now
}
