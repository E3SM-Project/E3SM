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
* \file trade_demand_function.cpp
* \ingroup Objects
* \brief The TradeDemandFunction class source file.
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

#include "functions/include/trade_demand_function.h"
#include "functions/include/iinput.h"
#include "containers/include/scenario.h"
#include "marketplace/include/marketplace.h"
#include "util/base/include/util.h"
#include "functions/include/function_utils.h"
#include "util/logger/include/ilogger.h"

using namespace std;

extern Scenario* scenario;

double TradeDemandFunction::calcDemand( InputSet& input, double consumption, const string& regionName,
                                        const string& sectorName, const double aShutdownCoef, int period,
                                        double capitalStock, double alphaZero, double sigma, double IBT ) const 
{
	double totalNetExport = 0; // total demand used for scaling
	Marketplace* marketplace = scenario->getMarketplace();
	for( unsigned int i = 0; i<input.size(); ++i ){
        // Trade exists in all comodities except land and labor.
		if( !input[ i ]->hasTypeFlag( IInput::FACTOR ) ){
			// Open Trade
            double netExport;
            if( FunctionUtils::isFixedPrice( regionName, input[ i ]->getName(), period ) ){
				netExport = marketplace->getSupply( input[ i ]->getName(), regionName, period ) -
					marketplace->getDemand( input[ i ]->getName(), regionName, period );
			}
			// Fixed Trade
			else {
				netExport = input[ i ]->getPhysicalDemand( period );
			}
			// set here adds to marketplace demand as well as setting net export in input
			input[ i ]->setPhysicalDemand( netExport, regionName, period );
            totalNetExport += netExport * input[ i ]->getPrice( regionName, period );
		}
		// for capital, add to market demand but not to total net export
		else if( input[ i ]->hasTypeFlag( IInput::CAPITAL ) ){
			double netExport = input[ i ]->getPhysicalDemand( period );
			// TODO: a subtle thing to note here is demand for capital will 
            // not be added to the marketplace directly although I could get
            // this to work by allowing a setCurrencyDemand to add to the marketplace
			input[ i ]->setPhysicalDemand( netExport, regionName, period );
            marketplace->addToSupply( input[ i ]->getName(), regionName, netExport, period );
		}
	}
    assert( util::isValidNumber( totalNetExport ) );
	return totalNetExport;
}

double TradeDemandFunction::calcCoefficient( InputSet& input, double consumption, const string& regionName,
                                             const string& sectorName, int period, double sigma, double IBT,
                                             double capitalStock ) const
{
	return 0; // function not needed for trade consumer class
}
