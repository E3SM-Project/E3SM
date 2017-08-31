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
* \file ccsm_scenario.cpp
* \ingroup Objects
* \brief Ccsm_Scenario class source file.
* \author Sonny Kim
*/              

#include "util/base/include/definitions.h"
#include <string>
#include <fstream>
#include <cassert>
#include <ctime>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>

#include "containers/include/scenario.h"
#include "containers/include/ccsm_scenario.h"
#include "containers/include/world.h"
#include "util/base/include/xml_helper.h"
#include "util/base/include/configuration.h"
#include "util/logger/include/ilogger.h"
#include "util/curves/include/curve.h"
#include "solution/solvers/include/solver.h"
#include "util/base/include/auto_file.h"
#include "reporting/include/graph_printer.h"
#include "reporting/include/land_allocator_printer.h"
#include "reporting/include/xml_db_outputter.h"
#include "containers/include/output_meta_data.h"
#include "solution/solvers/include/solver_factory.h"
#include "solution/solvers/include/bisection_nr_solver.h"
#include "solution/util/include/solution_info_param_parser.h"

using namespace std;
using namespace xercesc;
using namespace boost;

//! Default constructor
Ccsm_Scenario::Ccsm_Scenario() : Scenario() {
  scenarioTabs = new Tabs;
}

//! Destructor
Ccsm_Scenario::~Ccsm_Scenario() {
}

/*! \brief Run the scenario.
* \param aSinglePeriod Single period to run
* \param aPrintDebugging Whether to print extra debugging files.
* \param aFilenameEnding The string to add to the end of the debug output file
*        for uniqueness.
* \return Whether all model runs solved successfully.
*/
bool Ccsm_Scenario::run( const int aSinglePeriod,
                    const bool aPrintDebugging,
                    const string& aFilenameEnding )
{
    // Avoid accumulating unsolved periods.
    unsolvedPeriods.clear();
    
    // Open the debugging files.
    //    ofstream XMLDebugFile;
    //    ofstream SGMDebugFile;
    //    Tabs tabs;
  
    if (aSinglePeriod == 0) {
      if( aPrintDebugging ){
        openDebuggingFiles( XMLDebugFile, SGMDebugFile, scenarioTabs, aFilenameEnding );
      }

      // Log that a run is beginning.
      logRunBeginning();
    }

    bool success = true;
    ILogger& mainLog = ILogger::getLogger( "main_log" );
    mainLog.setLevel( ILogger::NOTICE );

    // Check if the single period is invalid.
    if( aSinglePeriod >= modeltime->getmaxper() ){
      mainLog.setLevel( ILogger::ERROR );
      mainLog << "Invalid single period " << aSinglePeriod << " passed to run method." << endl;
      success = false;
    } 
    else {
      // Run a period
      mainLog.setLevel( ILogger::NOTICE );
      mainLog << "calculate Period: " << aSinglePeriod << endl;
      success &= calculatePeriod( aSinglePeriod, XMLDebugFile, SGMDebugFile, scenarioTabs, aPrintDebugging );
    }
    
    // Print any unsolved periods.
    // TODO: This should be added to the db.
    mainLog.setLevel( ILogger::WARNING );
    
    // Report if all model periods solved correctly.
    if( getUnsolvedPeriods().empty() ) {
        mainLog << "Model period solved correctly." << endl;
    }
    else {
        // Otherwise print all model periods which did not solve correctly.
        mainLog << "The following model periods did not solve: ";
        for( vector<int>::const_iterator i = getUnsolvedPeriods().begin(); i != getUnsolvedPeriods().end(); i++ ) {
            mainLog << *i << ", ";
        }
        mainLog << endl;
    }

    // Run the climate model.

    if( aSinglePeriod == getModeltime()->getmaxper()-1 ){
      //      mainLog << "running climate model for period " << aSinglePeriod << endl;
      //      world->runClimateModel();

      // Log that the run has finished.
      logRunEnding();

      // Close the debugging files.
      if( aPrintDebugging ){
        closeDebuggingFiles( XMLDebugFile, SGMDebugFile, scenarioTabs );
      }
    
    }
    return success;
}
