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
* \file main.cpp
* \todo Update this documentation.
* \brief This is the Main program file which controls initialization of run 
* parameters and sets up scenarios to run from input file names in the 
* Configuration file.
*
* The program reads in parameters and file names stored in the 
* configuration file and assigns them to the configuration object.  It creates 
* a scenario object and triggers a read-in of all input data by calling XMLParse.
* If there are scenario components (ScenComponents), then these are read next.
* 
* A switch is checked for the program to run in a mode that just merges files
* into one xml input file rather than running the model (mergeFilesOnly).
* 
* If the model is to be run, it is triggered by the ScenarioRunner class object, which 
* calls runScenarios() to trigger running of the complete model for all periods.
*
* \author Sonny Kim
*/

#include "util/base/include/definitions.h"

// include standard libraries
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <list>

// xerces xml headers
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include "util/base/include/xml_helper.h"

// include custom headers
#include "util/base/include/configuration.h"
#include "containers/include/scenario.h"
#include "containers/include/iscenario_runner.h"
#include "containers/include/scenario_runner_factory.h"
#include "util/logger/include/ilogger.h"
#include "util/logger/include/logger_factory.h"
#include "util/base/include/timer.h"
#include "util/base/include/version.h"

using namespace std;
using namespace xercesc;

// define file (ofstream) objects for outputs, debugging and logs
/* \todo Finish removing globals-JPL */
ofstream outFile;

// Initialize time and set some pointers to null.
// Declared outside Main to make global.
Scenario* scenario; // model scenario info

void parseArgs( unsigned int argc, char* argv[], string& confArg, string& logFacArg );

//! Main program. 
int main( int argc, char *argv[] ) {

    // identify default file names for control input and logging controls
    string configurationArg = "configuration.xml";
    string loggerFactoryArg = "log_conf.xml";
    // Parse any command line arguments.  Can override defaults with command lone args
    parseArgs( argc, argv, configurationArg, loggerFactoryArg );

    // Add OS dependent prefixes to the arguments.
    const string configurationFileName = configurationArg;
    const string loggerFileName = loggerFactoryArg;

    // Initialize the timer.  Create an object of the Timer class.
    Timer timer;
    timer.start();

    // Initialize the LoggerFactory
    LoggerFactoryWrapper loggerFactoryWrapper;
    bool success = XMLHelper<void>::parseXML( loggerFileName, &loggerFactoryWrapper );
    
    // Check if parsing succeeded. Non-zero return codes from main indicate
    // failure.
    if( !success ){
        return 1;
    }


    // Get the main log file.
    ILogger& mainLog = ILogger::getLogger( "main_log" );
    mainLog.setLevel( ILogger::WARNING );

    // print disclaimer
    mainLog << "This computer software was prepared by Battelle Memorial Institute," << endl;
    mainLog << "hereinafter the Contractor, under Contract No. DE-AC05-76RL0 1830 with" << endl;
    mainLog << "the Department of Energy (DOE). NEITHER THE GOVERNMENT NOR THE" << endl;
    mainLog << "CONTRACTOR MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY" << endl;
    mainLog << "LIABILITY FOR THE USE OF THIS SOFTWARE. This notice including this" << endl;
    mainLog << "sentence must appear on any copies of this computer software." << endl << endl;
    
    // print export control notice
    mainLog << "User agrees that the Software will not be shipped, transferred or" << endl;
    mainLog << "exported into any country or used in any manner prohibited by the United" << endl;
    mainLog << "States Export Administration Act or any other applicable export laws," << endl;
    mainLog << "restrictions or regulations (collectively the 'Export Laws'). Export of" << endl;
    mainLog << "the Software may require some form of license or other authority from" << endl;
    mainLog << "the U.S. Government, and failure to obtain such export control license" << endl;
    mainLog << "may result in criminal liability under U.S. laws. In addition, if the" << endl;
    mainLog << "Software is identified as export controlled items under the Export Laws," << endl;
    mainLog << "User represents and warrants that User is not a citizen, or otherwise" << endl;
    mainLog << "located within, an embargoed nation (including without limitation Iran," << endl;
    mainLog << "Syria, Sudan, Cuba, and North Korea) and that User is not otherwise" << endl;
    mainLog << "prohibited under the Export Laws from receiving the Software. Note that" << endl;
    mainLog << "the United States Department of Commerce has determined that this" << endl;
    mainLog << "software has been assigned an export licensing code (ECCN) of EAR99." << endl;
    mainLog << "This ECCN applies to both the ObjECTS modeling framework and the ObjECTS" << endl;
    mainLog << "Model Interface. All rights to use the Software are granted on condition" << endl;
    mainLog << "that such Distributed as open-source under the terms of the Educational Community " << endl;
    mainLog << "this Agreement. User agrees to indemnify, defend and hold harmless" << endl;
    mainLog << "BATTELLE, its officers, agents and employees from all liability" << endl;
    mainLog << "involving the violation of such Export Laws, either directly or" << endl;
    mainLog << "indirectly, by User." << endl << endl;

    mainLog << "Running GCAM model code base version " << __ObjECTS_VER__ << " revision "
        << __REVISION_NUMBER__ << endl << endl;

    // Parse configuration file.
    mainLog.setLevel( ILogger::NOTICE );
    mainLog << "Parsing input files..." << endl;
    Configuration* conf = Configuration::getInstance();
    success = XMLHelper<void>::parseXML( configurationFileName, conf );
    // Check if parsing succeeded. Non-zero return codes from main indicate
    // failure.
    if( !success ){
        return 1;
    }

    // Create an empty exclusion list so that any type of IScenarioRunner can be
    // created.
    list<string> exclusionList;

    // Create an auto_ptr to the scenario runner. This will automatically
    // deallocate memory.
    auto_ptr<IScenarioRunner> runner = ScenarioRunnerFactory::createDefault( exclusionList );
    
    // Setup the scenario.
    success = runner->setupScenarios( timer );
    // Check if setting up the scenario, which often includes parsing,
    // succeeded.
    if( !success ){
        return 1;
    }

    // Run the scenario and print debugging information.
    const int terminalPeriod = conf->getInt( "terminal-period", Scenario::RUN_ALL_PERIODS );
    success = runner->runScenarios( terminalPeriod, true, timer );

    // Print the output.
    runner->printOutput( timer );
    mainLog.setLevel( ILogger::WARNING ); // Increase level so that user will know that model is done
    mainLog << "Model exiting successfully." << endl;
    // Cleanup Xerces. This should be encapsulated with an initializer object to ensure against leakage.
    XMLHelper<void>::cleanupParser();
    
    // Return exit code based on whether the model succeeded(Non-zero is failure by convention).
    return success ? 0 : 1; 
}

/*!
* \brief Function to parse the command line arguments.
* \param argc Number of arguments.
* \param argv List of arguments.
* \param confArg [out] Name of the configuration file.
* \param logFacArg [out] Name of the log configuration file.
* \todo Allow a space between the flags and the file names.
*/
void parseArgs( unsigned int argc, char* argv[], string& confArg, string& logFacArg ) {
    for( unsigned int i = 1; i < argc; i++ ){
        string temp( argv[ i ] );
        if( temp.compare(0,2,"-C" ) == 0 ){
            confArg = temp.substr( 2, temp.length() );
        } 
        else if( temp.compare(0,2,"-L" ) == 0 ){
            logFacArg = temp.substr( 2, temp.length() );
        }
        else {
            cout << "Invalid argument: " << temp << endl;
            cout << "Usage: " << argv[ 0 ] << " [-CconfigurationFileName ][ -LloggerFactoryFileName ]";
        }
    }
}
