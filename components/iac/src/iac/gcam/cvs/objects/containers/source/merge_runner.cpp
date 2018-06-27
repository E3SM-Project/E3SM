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
* \file merge_runner.cpp
* \ingroup objects
* \brief MergeRunner class source file.
* \author Josh Lurz
*/

#include "util/base/include/definitions.h"
#include <xercesc/dom/DOMNode.hpp>
#include "containers/include/merge_runner.h"
#include "util/base/include/timer.h"
#include "util/base/include/xml_helper.h"
#include "util/base/include/configuration.h"
#include "containers/include/scenario.h"
#include "util/logger/include/ilogger.h"

using namespace std;
using namespace xercesc;

//! Constructor
MergeRunner::MergeRunner(){
}

//! Destructor
MergeRunner::~MergeRunner(){
}

const string& MergeRunner::getName() const {
    return getXMLNameStatic();
}

// IParsable interface
bool MergeRunner::XMLParse( const xercesc::DOMNode* aRoot ){
    // No data to parse.
    return true;
}

//! Setup the scenario.
bool MergeRunner::setupScenarios( Timer& timer, const string aName, const list<string> aScenComponents ){
    // Parse the input file.
    const Configuration* conf = Configuration::getInstance();
    
    bool success = XMLHelper<void>::parseXML( conf->getFile( "xmlInputFileName" ), mScenario.get() );

    // Parsing failed.
    if( !success ){
        return false;
    }

    // Need to ensure the Scenario is cleared and set.
    mScenario.reset( new Scenario );
    // Make sure the global scenario points is set.
    scenario = mScenario.get();
    
    // Override scenario name from data file with that from configuration file
    const string overrideName = conf->getString( "scenarioName" ) + aName;
    if ( !overrideName.empty() ) {
        mScenario->setName( overrideName );
    }

    // Fetch the listing of Scenario Components.
    const list<string> scenComponents = conf->getScenarioComponents();

    // Get the main log file.
    ILogger& mainLog = ILogger::getLogger( "main_log" );
    mainLog.setLevel( ILogger::NOTICE );

    // Iterate over the vector.
    typedef list<string>::const_iterator ScenCompIter;
    for( ScenCompIter currComp = scenComponents.begin(); currComp != scenComponents.end(); ++currComp ) {
        mainLog << "Parsing " << *currComp << " scenario component." << endl;
        if( !XMLHelper<void>::parseXML( *currComp, mScenario.get() ) ){
            // Parsing failed.
            return false;
        }
    }
    
    mainLog << "XML parsing complete." << endl;

    // Print data read in time.
    timer.stop();
    mainLog.setLevel( ILogger::DEBUG );
    timer.print( mainLog, "XML Readin Time:" );
    return true;
}
/*! \brief Does nothing, needed for interface.
* \param aSinglePeriod This parameter is ignored.
* \param aTimer This parameter is ignored.
* \return Always returns true.
* \author Josh Lurz
*/
bool MergeRunner::runScenarios( const int aSinglePeriod,
                                const bool aPrintDebugging,
                                Timer& aTimer )
{
    return true;
}

void MergeRunner::printOutput( Timer& timer, const bool aCloseDB ) const {
    // Print output xml file.
    const Configuration* conf = Configuration::getInstance();
    const string xmlOutFileName = conf->getFile( "xmlOutputFileName", "output.xml" );
    ofstream xmlOut;
    xmlOut.open( xmlOutFileName.c_str(), ios::out );
    util::checkIsOpen( xmlOut, xmlOutFileName );

    Tabs tabs;
    mScenario->toInputXML( xmlOut, &tabs );

    // Close the output file. 
    xmlOut.close();
}

/*! \brief Get the internal scenario.
 \return The internal scenario.
*/
Scenario* MergeRunner::getInternalScenario(){
	return mScenario.get();
}

/*! \brief Get the internal scenario.
* \return Constant pointer to the internal scenario.
*/
const Scenario* MergeRunner::getInternalScenario() const {
	return mScenario.get();
}

/*! \brief Get the static name of the class.
* \return The static name of the class.
*/
const string& MergeRunner::getXMLNameStatic(){
	static const string XML_NAME = "merge-runner";
	return XML_NAME;
}
