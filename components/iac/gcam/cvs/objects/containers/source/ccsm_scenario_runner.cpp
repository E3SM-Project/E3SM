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
 * \file ccsm_scenario_runner.cpp
 * \ingroup CcsmRunner
 * \brief CcsmRunner class source file.
 * \author Josh Lurz
 */

#include "util/base/include/definitions.h"
#include <fstream>
#include <cassert>
#include <xercesc/dom/DOMNode.hpp>
#include "containers/include/ccsm_scenario_runner.h"
#include "containers/include/ccsm_scenario.h"
#include "containers/include/scenario.h"
#include "util/base/include/xml_helper.h"
#include "util/base/include/timer.h"
#include "util/base/include/configuration.h"
#include "util/base/include/auto_file.h"
#include "util/logger/include/ilogger.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;
extern ofstream outFile;

// Function Prototypes. These need a helper class. 
extern void createMCvarid();
extern void closeDB();
extern void openDB();
extern void createDBout();

/*! \brief Constructor */
CcsmRunner::CcsmRunner(){
}

//! Destructor.
CcsmRunner::~CcsmRunner(){
}

const string& CcsmRunner::getName() const {
    return getXMLNameStatic();
}

// IParsable interface
bool CcsmRunner::XMLParse( const xercesc::DOMNode* aRoot ){
    // No data to parse.
    return true;
}

bool CcsmRunner::setupScenarios( Timer& timer,
				 const string aName,
				 const list<string> aScenComponents )
{
    // Ensure that a new scenario is created for each run.
    mScenario.reset( new /*Ccsm_*/Scenario );

    ILogger& mainLog = ILogger::getLogger( "main_log" );

    // Set the global scenario pointer.
    // TODO: Remove global scenario pointer.
    scenario = mScenario.get();

    // Parse the input file.
    const Configuration* conf = Configuration::getInstance();
    bool success;

    bool restart = (conf->getInt( "restart-period", 0 )!=0);

    if (restart) {
      ifstream inFile;
      inFile.open("rpointer.gcam");
      if (!inFile) {
        mainLog.setLevel( ILogger::SEVERE );
	mainLog << "Unable to open restart file: rpointer.gcam" << endl;
        success = false;
      }else{
	string xmlRestartInputFileName;
	int restartPeriod,restartYear;
        inFile >> xmlRestartInputFileName;
        inFile >> restartPeriod;
        inFile >> restartYear;
        success = XMLHelper<void>::parseXML( xmlRestartInputFileName, 
					     mScenario.get() );
      }
    }else{
      success =
        XMLHelper<void>::parseXML( conf->getFile( "xmlInputFileName" ),
                                   mScenario.get() );
    }
    // Check if parsing succeeded.
    if( !success ){
        return false;
    }

    // Fetch the listing of Scenario Components.
    list<string> scenComponents = conf->getScenarioComponents();

    // Add on any scenario components that were passed in.
    for( list<string>::const_iterator curr = aScenComponents.begin();
		curr != aScenComponents.end(); ++curr )
	{
        scenComponents.push_back( *curr );
    }
    
    // Iterate over the vector.
    typedef list<string>::const_iterator ScenCompIter;
    //    ILogger& mainLog = ILogger::getLogger( "main_log" );
    for( ScenCompIter currComp = scenComponents.begin();
		 currComp != scenComponents.end(); ++currComp )
	{
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Parsing " << *currComp << " scenario component." << endl;
        success = XMLHelper<void>::parseXML( *currComp, mScenario.get() );
        
        // Check if parsing succeeded.
        if( !success ){
            return false;
        }
    }
    
    // Override scenario name from data file with that from configuration file
    const string overrideName = conf->getString( "scenarioName" ) + aName;
    if ( !overrideName.empty() ) {
        mScenario->setName( overrideName );
    }

    mainLog.setLevel( ILogger::NOTICE );
    mainLog << "XML parsing complete." << endl;

    // Print data read in time.
    timer.stop();
    mainLog.setLevel( ILogger::DEBUG );
    timer.print( mainLog, "XML Readin Time:" );

    // Finish initialization.
    if( mScenario.get() ){
        mScenario->completeInit();
    }

    return true;
}

bool CcsmRunner::runScenarios( const int aSinglePeriod,
                                        const bool aPrintDebugging,
                                        Timer& aTimer )
{
    ILogger& mainLog = ILogger::getLogger( "main_log" );
    mainLog.setLevel( ILogger::NOTICE );
    mainLog << "Starting a model run. Running ";
    if( aSinglePeriod == -1 ){
        mainLog << "all periods.";
    }
    else {
        mainLog << "period " << aSinglePeriod;
    }
    mainLog << endl;
    
	bool success = false;
	if( mScenario.get() ){
		// Perform the initial run of the scenario.
        success = mScenario->run( aSinglePeriod, aPrintDebugging,
                                  mScenario->getName() );

		// Compute model run time.
	//jt		aTimer.stop();
	//jt		mainLog.setLevel( ILogger::DEBUG );
	//jt		aTimer.print( mainLog, "Data Readin & Initial Model Run Time:" );
	}
	else {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::SEVERE );
        mainLog << "No scenario container was parsed from the input files. Aborting scenario run!" << endl;
	}

    // Return whether the scenario ran correctly. 
    return success;
}

void CcsmRunner::writeRestart(const int aPeriod,const int aYear ) const {
  ILogger& mainLog = ILogger::getLogger( "main_log" );
  mainLog.setLevel( ILogger::NOTICE );
  mainLog << "Writing Restarts" << endl;
  
  // Print restart xml file.
  string restartFileName= "gcam_restart_"+boost::lexical_cast<string>( aYear )+".xml";
  AutoOutputFile rpointer( "rpointer.gcam" );
  rpointer << restartFileName << "\n";
  rpointer << aPeriod << "\n";
  rpointer << aYear << "\n";
  
  AutoOutputFile xmlOut( restartFileName );
  Tabs tabs;
  mScenario->toInputXML( *xmlOut, &tabs );
}


void CcsmRunner::printOutput( Timer& aTimer, const bool aCloseDB ) const {

    ILogger& mainLog = ILogger::getLogger( "main_log" );
    mainLog.setLevel( ILogger::NOTICE );
    mainLog << "Printing output" << endl;

    // Print output xml file.
    AutoOutputFile xmlOut( "xmlOutputFileName", "output.xml" );
    Tabs tabs;
    mScenario->toInputXML( *xmlOut, &tabs );

    // Write csv file output
    mScenario->writeOutputFiles();

    static const bool printDB = Configuration::getInstance()->getBool( "write-access-db", true );
    if( printDB ){
        // Perform the database output. 
	    // Open MS Access database
        openDB();
	    // create main database output table before calling output routines
        createDBout();
        mScenario->dbOutput();

        if( aCloseDB ){
            createMCvarid(); // create MC variable id's     
            // close MS Access database
            closeDB();
        }
    }

    if( aCloseDB ){
        outFile.close();
    }

#if( __USE_XML_DB__ )
    static const bool printXMLDB = Configuration::getInstance()->getBool( "write-xml-db", true );
    if( printXMLDB ){
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Starting output to XML Database." << endl;
        // Print the XML file for the XML database.
        scenario->printOutputXML();
    }
#endif

     // Print the timestamps.
    aTimer.stop();
    mainLog.setLevel( ILogger::DEBUG );
    aTimer.print( mainLog, "Data Readin, Model Run & Write Time:" );
    mainLog.setLevel( ILogger::NOTICE );
    mainLog << "Model run completed." << endl;
}

Scenario* CcsmRunner::getInternalScenario(){
	return mScenario.get();
}

const Scenario* CcsmRunner::getInternalScenario() const {
	return mScenario.get();
}

/*!
 * \brief Get the XML name of the class.
 * \return The XML name of the class.
 */
const string& CcsmRunner::getXMLNameStatic(){
	static const string XML_NAME = "ccsm-scenario-runner";
	return XML_NAME;
}
