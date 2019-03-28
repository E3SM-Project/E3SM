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
* \file function_manager.cpp
* \ingroup Objects
* \brief The FunctionManager class source file.
* \author Pralit Patel
* \author Josh Lurz
*/
#include "util/base/include/definitions.h"
#include <string>
#include <map>

#include "util/base/include/util.h"
#include "functions/include/function_manager.h"
#include "functions/include/household_demand_function.h"
#include "functions/include/trade_demand_function.h"
#include "functions/include/government_demand_function.h"
#include "functions/include/investment_demand_function.h"
#include "functions/include/ces_production_function.h"
#include "functions/include/leontief_production_function.h"
#include "functions/include/minicam_leontief_production_function.h"
#include "functions/include/minicam_price_elasticity_function.h"
#include "functions/include/nested_ces_production_function.h"
#include "functions/include/utility_demand_function.h"
#include "functions/include/logit_production_function.h"
#include "util/logger/include/ilogger.h"

using namespace std;

//! Default Constructor
FunctionManager::FunctionManager() {
    // all avaiable functions must be added here
    /*!
     * \todo To avoid spelling issues perhaps we should change the key to
     *       an enumeration and have a static lookup between a name and the
     *       enum for use during XML parsing.
     */
    mFunctions[ "CES" ] = new CESProductionFunction;
    mFunctions[ "Leontief" ] = new LeontiefProductionFunction;
    mFunctions[ "NestedCES" ] = new NestedCESProductionFunction;
    mFunctions[ "minicam-leontief" ] = new MinicamLeontiefProductionFunction;
    mFunctions[ "minicam-price-elasticity" ] = new MinicamPriceElasticityFunction;
    mFunctions[ "HouseholdDemandFn" ] = new HouseholdDemandFunction;
    mFunctions[ "GovtDemandFn" ] = new GovernmentDemandFunction;
    mFunctions[ "TradeDemandFn" ] = new TradeDemandFunction;
    mFunctions[ "InvestDemandFn" ] = new InvestmentDemandFunction;
    mFunctions[ "UtilityDemandFunction" ] = new UtilityDemandFunction;
    mFunctions[ "Logit" ] = new LogitProductionFunction;
}

//! Destructor
FunctionManager::~FunctionManager() {
    // delete pointer to each function
    for( FunctionsIterator funcIter = mFunctions.begin(); funcIter != mFunctions.end(); ++funcIter ) {
        delete funcIter->second;
    }
}

/*! \brief Get the pointer the the function class represented by aFunctionName
* \details This method access the static instance of the FunctionManager and
*          searches its internal function map for the specified function name.
*          The static internal function manager is instantiated on first access
*          to this function and destroyed when the model run is complete.
* \author Josh Lurz, Pralit Patel
* \param aFunctionName The name of the function you want to access.
* \return The pointer to the requested function.
*/
const IFunction* FunctionManager::getFunction( const string& aFunctionName ) {
    // Allocate the static function manager if it does not already exist.
    const static FunctionManager functionManager;

    const IFunction* tempFn = util::searchForValue( functionManager.mFunctions, aFunctionName );
    // checking to see if functionName exists in the map
    if ( !tempFn ) {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "Could not find Production or Demand Function. Check function type: " 
                << aFunctionName << endl;
    }
    return tempFn;
}
