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
* \file mac_generator_scenario_runner.cpp
* \ingroup Objects
* \brief MACGeneratorScenarioRunner class source file.
* \author Josh Lurz
*/

#include "util/base/include/definitions.h"
#include <cassert>
#include <string>
#include "containers/include/mac_generator_scenario_runner.h"
#include "containers/include/scenario_runner_factory.h"
#include "util/base/include/configuration.h"
#include "util/logger/include/ilogger.h"
#include "containers/include/total_policy_cost_calculator.h"
#include "containers/include/scenario.h"
#include "util/base/include/auto_file.h"
#include "marketplace/include/marketplace.h"

using namespace std;
using namespace xercesc;

extern void closeDB();
extern void createMCvarid();
extern ofstream outFile;

/*! \brief Constructor.
*/
MACGeneratorScenarioRunner::MACGeneratorScenarioRunner(){
    mSingleScenario = ScenarioRunnerFactory::create( "single-scenario-runner" );

    // Check to make sure calibration is off.
    const Configuration* conf = Configuration::getInstance();
    if( conf->getBool( "debugChecking" ) && conf->getBool( "CalibrationActive" ) ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Calibration is incompatible with the generation of marginal abatement curves." << endl;
    }
    else {
        // Create the policy cost calculator.
        mPolicyCostCalculator.reset( new TotalPolicyCostCalculator( mSingleScenario.get() ) );
    }
}

//! Destructor
MACGeneratorScenarioRunner::~MACGeneratorScenarioRunner(){
}

const string& MACGeneratorScenarioRunner::getName() const {
    return getXMLNameStatic();
}

bool MACGeneratorScenarioRunner::XMLParse( const xercesc::DOMNode* aRoot ){
    // No data to parse.
    return true;
}

/*! \brief Setup the Scenario to be run.
* \details This function sets up the contained SingleScenarioRunner.
* \param aTimer The timer used to print out the amount of time spent performing
*        operations.
* \param aName The name to add on to the name read in in the Configuration file.
* \param aScenComponents A list of additional scenario components to read in.
* \return Whether the setup completed successfully.
*/
bool MACGeneratorScenarioRunner::setupScenarios( Timer& aTimer, const string aName, const list<string> aScenComponents ){
    return mSingleScenario->setupScenarios( aTimer, aName, aScenComponents );
}

/*! \brief Function which handles running the scenario and optionally computing
*          a cost curve.
* \details This function wraps around the scenario so that scenario can be
*          called multiple times if necessary to create an abatement cost
*          curve. This function first calls the scenario regularly, outputs all
*          data, and then calls scenario several more times and calculates the
*          abatement cost.
* \param aSinglePeriod This parameter is ignored currently.
* \param aTimer The timer used to print out the amount of time spent performing
*        operations.
* \return Whether all model runs solved successfully.
* \author Josh Lurz
*/
bool MACGeneratorScenarioRunner::runScenarios( const int aSinglePeriod,
                                               const bool aPrintDebugging,
                                               Timer& aTimer )
{
    // TODO: This is only necessary because the cost calculator has trouble solving
    // a zero carbon tax with restart turned on.  This should be unnecessary with
    // a better solver.
    const static bool usingRestartPeriod = Configuration::getInstance()->getInt(
        "restart-period", -1 ) != -1;
    if( usingRestartPeriod ) {
        mSingleScenario->getInternalScenario()->getMarketplace()->store_prices_for_cost_calculation();
    }
    // Run the base scenario. Print debugging for the base scenario run.
    bool success = mSingleScenario->runScenarios( Scenario::RUN_ALL_PERIODS,
                                                  aPrintDebugging, aTimer );

    // Print the output now before it is overwritten.
    mSingleScenario->printOutput( aTimer, false );

    // Now calculate the abatement curves.
    if( mPolicyCostCalculator.get() ){
        success &= mPolicyCostCalculator->calculateAbatementCostCurve();
    }

    // Return whether the initial run and all datapoint calculations completed
    // successfully.
    return success;
}

//! Print the output.
void MACGeneratorScenarioRunner::printOutput( Timer& timer, const bool aCloseDB ) const {
    if( mPolicyCostCalculator.get() ){
        mPolicyCostCalculator->printOutput();
    }
    
    static const bool printDB = Configuration::getInstance()->getBool( "write-access-db", true );
    
    // Close the database.
    if( printDB && aCloseDB ){
        createMCvarid();
        closeDB();
        outFile.close();
    }
}

/*! \brief Get the internal scenario.
* \return The internal scenario.
*/
Scenario* MACGeneratorScenarioRunner::getInternalScenario(){
	return mSingleScenario->getInternalScenario();
}

/*! \brief Get the internal scenario.
* \return Constant pointer to the internal scenario.
*/
const Scenario* MACGeneratorScenarioRunner::getInternalScenario() const {
	return mSingleScenario->getInternalScenario();
}

/*! \brief Get the static class name.
* \return The static class name.
*/
const string& MACGeneratorScenarioRunner::getXMLNameStatic(){
	static const string XML_NAME = "mac-generator-scenario-runner";
	return XML_NAME;
}
