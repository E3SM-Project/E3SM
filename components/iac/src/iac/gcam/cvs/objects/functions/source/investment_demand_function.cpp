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
* \file investment_demand_function.cpp
* \ingroup Objects
* \brief The InvestmentDemandFunction class source file.
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

#include "functions/include/investment_demand_function.h"
#include "functions/include/iinput.h"
#include "containers/include/scenario.h"
#include "marketplace/include/marketplace.h"
#include "util/base/include/util.h"
#include "functions/include/function_utils.h"
#include "util/logger/include/ilogger.h"

using namespace std;

extern Scenario* scenario;

double InvestmentDemandFunction::calcDemand( InputSet& input, double capitalTotal, const string& regionName,
											 const string& sectorName, const double aShutdownCoef,
											 int period, double capitalStock, double alphaZero,
											 double sigma, double IBT ) const
{
	double capitalSum = 0;
    for( unsigned int i = 0; i < input.size(); ++i ){
		if( !input[ i ]->hasTypeFlag( IInput::FACTOR ) ) {
			double demand = capitalTotal * input[ i ]->getCoefficient( period );
            if( demand < 0 ){
				ILogger& mainLog = ILogger::getLogger( "main_log" );
				mainLog.setLevel( ILogger::WARNING );
                mainLog << "Trying to add negative demand currency for " << input[ i ]->getName()
					    << " in " << sectorName << endl;
            }
			input[ i ]->setCurrencyDemand( demand, regionName, period );
             // looking at the fortran seems as though this should be price paid
            // Yes it should be.
            capitalSum += demand * input[i]->getPrice( regionName, period );
        }
	}
	return capitalSum;
}

// Calculate coefficients for Investment's demand function
double InvestmentDemandFunction::calcCoefficient( InputSet& input, double consumption,
												  const string& regionName, const string& sectorName,
												  int period, double sigma, double IBT,
												  double capitalStock ) const
{
    return 0; // Investment coefficients are read in, not calculated
}
