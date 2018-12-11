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
 * \file production_state_factory.cpp
 * \ingroup Objects
 * \brief ProductionStateFactory source file.
 * \author Josh Lurz
 */

#include "util/base/include/definitions.h"
#include <string>
#include "technologies/include/production_state_factory.h"
#include "util/base/include/model_time.h"
#include "containers/include/scenario.h"

// Add new types here.
#include "technologies/include/fixed_production_state.h"
#include "technologies/include/variable_production_state.h"
#include "technologies/include/vintage_production_state.h"
#include "technologies/include/retired_production_state.h"

extern Scenario* scenario;

using namespace std;

/*
* \brief Create the production state for a Technology in a given period.
* \details This function is called at the beginning of each period to create an
*          object that determines how a Technology will determine its output
*          level for the period. All Technologies begin in a retired state, this
*          represents that they have not been built yet. In the initial
*          production period of a Technology, it will set its production state
*          to either fixed or variable. The fixed state is used when a level of
*          fixed investment is specified by the data, the variable state is used
*          otherwise. In the period after the Technology's initial production
*          period, the Technology will change it's production state to either
*          vintaged or retired. Vintaged Technologies will continue in the same
*          state until their lifetime is exceeded at which point they will enter
*          the retired state. Retired Technologies do not produce any output.
* \param aInvestYear The year in which the Technology was new investment.
* \param aLifetimeYears The maximum number of years the Technology may operate.
* \param aFixedOutput The fixed output level, or
*        IProductionState::fixedOutputDefault if there is none.
* \param aInitialOutput The output in the initial operating period of the
*        technology, or IProductionState::fixedOutputDefault if it cannot be
*        calculated.
* \param aPeriod Model period.
*/
auto_ptr<IProductionState> ProductionStateFactory::create( const int aInvestYear,
                                                           const int aLifetimeYears,
                                                           const double aFixedOutput,
                                                           const double aInitialOutput,
                                                           const int aPeriod )
{
    // Initialize the production state.
    auto_ptr<IProductionState> newState;
    int currYear = scenario->getModeltime()->getper_to_yr( aPeriod );

    if( aInvestYear == currYear ){
        // If the new vintage has fixed output use a fixed output production
        // state.
        if( aFixedOutput != IProductionState::fixedOutputDefault() ){
            newState.reset( new FixedProductionState );
            newState->setBaseOutput( aFixedOutput, aInvestYear );
            
        }
        // Otherwise use a variable production state.
        else {
            newState.reset( new VariableProductionState );
        }
    }
    // Check if it is a still operating vintage.
    else if( ( currYear > aInvestYear ) &&
        ( aInvestYear + aLifetimeYears > currYear ) ){
        assert( aPeriod > 0 );
        newState.reset( new VintageProductionState );
        // Set the base level of output to the output in the initial investment
        // year.
        newState->setBaseOutput( aInitialOutput, aInvestYear );
    }
    // Otherwise it is retired. This may occur if the technology has not been
    // created yet as well.
    else {
        newState.reset( new RetiredProductionState );
    }
    return newState;
}
