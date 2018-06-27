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
 * \file batch_runner.cpp
 * \ingroup objects
 * \brief BatchRunner class source file.
 * \author Josh Lurz
 */
#if defined(_MSC_VER)
#pragma warning( disable : 4503 )
#endif 

#include "util/base/include/definitions.h"
#include <string>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include "containers/include/batch_runner.h"
#include "containers/include/scenario_runner_factory.h"
#include "util/base/include/timer.h"
#include "util/base/include/xml_helper.h"
#include "util/base/include/configuration.h"
#include "util/logger/include/ilogger.h"
#include "containers/include/scenario.h"
#include "reporting/include/batch_csv_outputter.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

typedef list<IScenarioRunner*>::iterator RunnerIterator;

/*!
 * \brief Constructor
 */
BatchRunner::BatchRunner() :
mInternalRunner( 0 ){ 
}

//! Destructor
BatchRunner::~BatchRunner(){
}

const string& BatchRunner::getName() const {
    return getXMLNameStatic();
}

bool BatchRunner::setupScenarios( Timer& aTimer, const string aName, const list<string> aScenComponents ){
    // Get the name of the batch file from the Configuration.
    const string batchFileName = Configuration::getInstance()->getFile( "BatchFileName" );

    // Add note so that if XML read fails here user knows what happened
    ILogger& mainLog = ILogger::getLogger( "main_log" );
    mainLog.setLevel( ILogger::NOTICE );
    mainLog << "Reading Batch File " << batchFileName << endl;

    // Parse the batch file.
    bool success = XMLHelper<void>::parseXML( batchFileName, this );

    // Create a default scenario runner if none were read in. This will be used to run all scenarios.
    if( mScenarioRunners.empty() ){
        // Don't allow another BatchRunner to be created.
        list<string> exclusionList;
        exclusionList.push_back( getXMLNameStatic() );

        mScenarioRunners.push_back( ScenarioRunnerFactory::createDefault( exclusionList ).release() );
    }
    return success;
}

bool BatchRunner::runScenarios( const int aSinglePeriod,
                                const bool aPrintDebugging,
                                Timer& aTimer )
{
    // Quick error checking for empty readin.
    ILogger& mainLog = ILogger::getLogger( "main_log" );
    if( mComponentSet.empty() ){
        mainLog.setLevel( ILogger::SEVERE );
        mainLog << "No scenario sets to run!" << endl;
        return false;
    }

    // Initialize each components iterator to the beginning of the vector. 
    for( ComponentSet::iterator currSet = mComponentSet.begin(); currSet != mComponentSet.end(); ++currSet ){
        currSet->mFileSetIterator = currSet->mFileSets.begin();
    }
    
    // The scenarios are created by determining all possible combinations of
    // file sets. The algorithm operates as follows:
    // 1) Set the current file set in each component to the initial position.
    // 2) Run the scenario.
    // 3) Set the current component to the first.
    // 4) Increment the current file set in the current component.
    // 5a) If this is a valid position in the current component and go to 2.
    // 5b) Otherwise, reset the current file set in the current component to
    //     the first position.
    // 6a) If the current component is the last component exit the algorithm.
    // 6b) Otherwise, increment the current component and go to 4.
    //
    // Example: assume there are two components A and B. A has two file sets
    // named 1 and 2, and B has two file sets named 3 and 4. The scenarios would
    // be run in the order: [A1, B1], [A2, B1], [A1, B2], [A2, B2]
    //
    // All generated scenarios are run with each scenario runner in the order in
    // which the scenario runners were read.
    bool shouldExit = false;
    bool success = true;
    BatchCSVOutputter csvOutputter;
    while( !shouldExit ){
        // The data structure containing the current run.
        Component fileSetsToRun;

        // Loop through the ComponentSet to create the current scenario.
        for( ComponentSet::const_iterator currSet = mComponentSet.begin(); currSet != mComponentSet.end(); ++currSet ){
            fileSetsToRun.mFileSets.push_back( *( currSet->mFileSetIterator ) );
            fileSetsToRun.mName += currSet->mFileSetIterator->mName;
        }

        // Run it using each possible type of IScenarioRunner.
        for( RunnerIterator runner = mScenarioRunners.begin(); runner != mScenarioRunners.end(); ++runner ){
            bool scenarioSuccess = runSingleScenario( *runner, fileSetsToRun, aTimer );
            success &= scenarioSuccess;
            (*runner)->getInternalScenario()->accept( &csvOutputter, -1 );
            csvOutputter.writeDidScenarioSolve( scenarioSuccess );
        }

        // Loop forward to find a position to increment.
        for( ComponentSet::iterator outPos = mComponentSet.begin(); outPos != mComponentSet.end(); ++outPos ){
            outPos->mFileSetIterator++;
            if( outPos->mFileSetIterator == outPos->mFileSets.end() ){
                outPos->mFileSetIterator = outPos->mFileSets.begin();

                // If the current iterator is to the end position of the last
                // component set there are no more scenarios to run.
                if( outPos == mComponentSet.end() - 1 ){
                    shouldExit = true;
                }
            }
            else {
                break;
            }
        }
    }
    return success;
}

void BatchRunner::printOutput( Timer& aTimer, const bool aCloseDB ) const {
    // Print out any scenarios that did not solve.
    ILogger& mainLog = ILogger::getLogger( "main_log" );
    mainLog.setLevel( ILogger::NOTICE );
    if( mUnsolvedNames.empty() ){
        mainLog << "All model runs completed successfully." << endl;
    }
    else {
        mainLog << "Model runs that did not solve correctly: " << endl;
        for( list<string>::const_iterator name = mUnsolvedNames.begin(); name != mUnsolvedNames.end(); ++name ){
            mainLog << *name << endl;
        }
    }
}

/*!
 * \brief Run a single scenario created by the BatchRunner.
 * \details Expands the list of FileSets into a list of scenario components
 *          files to parse. The scenario name is created by combining the names
 *          of all the file sets with the name read from the configuration file.
 *          The function selects the appropriate type of ScenarioRunner from the
 *          configuration file, and initializes it with the list of scenario
 *          components. It then runs the Scenario and prints its output.
 * \param aScenarioRunner The scenario runner to use for the scenario.
 * \param aComponent A named list of FileSets which is expanded to create the
 *        list of scenario files to read in.
 * \param aTimer The timer used to print out the amount of time spent performing
 *        operations.
 * \return Whether the model run solved successfully.
 */
bool BatchRunner::runSingleScenario( IScenarioRunner* aScenarioRunner,
                                     const Component& aComponent,
                                     Timer& aTimer )
{
    // Set the current scenario runner.
    mInternalRunner = aScenarioRunner;

    // Expand the file sets into a flat list and a scenario information string.
    list<string> components;
    
    typedef list<BatchRunner::FileSet>::const_iterator CFileSetIterator;
    typedef list<BatchRunner::File>::const_iterator CFileIterator;

    for( CFileSetIterator currFileSet = aComponent.mFileSets.begin(); currFileSet != aComponent.mFileSets.end(); ++currFileSet ){
        for( CFileIterator currFile = currFileSet->mFiles.begin(); currFile != currFileSet->mFiles.end(); ++currFile ){
            components.push_back( currFile->mPath );
        }
    }
    ILogger& mainLog = ILogger::getLogger( "main_log" );
    mainLog.setLevel( ILogger::WARNING );
    mainLog << "Running scenario " << aComponent.mName
            << " with scenario runner " << aScenarioRunner->getName()
            << "." << endl;

    // Setup the scenario.
    const string runName = aComponent.mName;
    bool success = mInternalRunner->setupScenarios( aTimer, runName, components );
    // Check if setting up the scenario, which often includes parsing,
    // succeeded.
    if( !success ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Scenario setup failed. Skipping run." << endl;
        return false;
    }

    // Run the scenario.
    success = mInternalRunner->runScenarios( Scenario::RUN_ALL_PERIODS, false, aTimer );
    
    // Print the output.
    mInternalRunner->printOutput( aTimer );
    
    // If the run failed, add to the list of failed runs. CHECK ME!
    if( !success ){
        mUnsolvedNames.push_back( aComponent.mName );
    }
    return success;
}

bool BatchRunner::XMLParse( const DOMNode* aRoot ){
    // assume we were passed a valid node.
    assert( aRoot );
    
    // get the children of the node.
    DOMNodeList* nodeList = aRoot->getChildNodes();
    
    // loop through the children
    bool success = true;
    for ( unsigned int i = 0; i < nodeList->getLength(); i++ ){
        DOMNode* curr = nodeList->item( i );
        string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );

        if( nodeName == XMLHelper<void>::text() ) {
            continue;
        }
        // This is a three level XMLParse.
        else if ( nodeName == "ComponentSet" ){
            success &= XMLParseComponentSet( curr );
        }
        else if( nodeName == "runner-set" ){
            success &= XMLParseRunnerSet( curr );
        }
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Unrecognized text string: " << nodeName << " found while parsing BatchScenarioRunner." << endl;
            success = false;
        }
    }
    return success;
}

/*!
 * \brief Parse a single scenario set element.
 * \details Parse a single ComponentSet and add it to the BatchRunner's list of
 *          ComponentSets. Dispatch any FileSets found to the XMLParseFileSet
 *          function.
 * \param aNode DOM node corresponding to the current ComponentSet.
 * \return Whether the node was successfully parsed.
 */
bool BatchRunner::XMLParseComponentSet( const DOMNode* aNode ){
    // assume we were passed a valid node.
    assert( aNode );
    
    // Create a new Component
    Component newComponent;

    // Get the name of the component set. 
    newComponent.mName = XMLHelper<string>::getAttr( aNode, XMLHelper<void>::name() );

    // get the children of the node.
    DOMNodeList* nodeList = aNode->getChildNodes();

    // loop through the children
    bool success = true;
    for ( unsigned int i = 0; i < nodeList->getLength(); i++ ){
        DOMNode* curr = nodeList->item( i );
        string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );

        if( nodeName == XMLHelper<void>::text() ) {
            continue;
        }
        else if ( nodeName == "FileSet" ){
            success &= XMLParseFileSet( curr, newComponent );
        }
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Unrecognized text string: " << nodeName << " found while parsing ComponentSet." << endl;
            success = false;
        }
    }
    // Add the new component
    mComponentSet.push_back( newComponent );
    return success;
}

/*!
 * \brief Parse the set of scenario runners.
 * \details Parse the set of scenario runners. Dispatches and XML data below the
 *          ScenarioRunner to the object itself for parsing.
 * \param aNode DOM node corresponding to the runner-set.
 * \return Whether the node was successfully parsed.
 */
bool BatchRunner::XMLParseRunnerSet( const DOMNode* aNode ){
    // assume we were passed a valid node.
    assert( aNode );
    
    // get the children of the node.
    DOMNodeList* nodeList = aNode->getChildNodes();

    // loop through the children
    bool success = true;
    for ( unsigned int i = 0; i < nodeList->getLength(); i++ ){
        const DOMNode* curr = nodeList->item( i );
        const string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );

        if( nodeName == XMLHelper<void>::text() ) {
            continue;
        }
        else if( nodeName == "Value" ){
            // The value should point at a file containing configuration
            // information for a scenario runner.
            ParseHelper parseHelper;
            if( XMLHelper<void>::parseXML( XMLHelper<string>::getValue( curr ), &parseHelper ) ){
                mScenarioRunners.push_back( parseHelper.getParsedScenarioRunner().release() );
            }
        }
        else if( ScenarioRunnerFactory::isOfType( nodeName ) ){
            // This is a shortcut to allow creating a IScenarioRunner directly
            // without creating a file. Most IScenarioRunners only have a tag
            // and no data so this is useful.
            // Ensure that a BatchRunner cannot create another BatchRunner.
            if( nodeName == getXMLNameStatic() ){
                ILogger& mainLog = ILogger::getLogger( "main_log" );
                mainLog.setLevel( ILogger::ERROR );
                mainLog << "Batch scenario runners cannot create Batch scenario runners." << endl;
                success = false;
            }
            else {
                IScenarioRunner* currRunner =
                    ScenarioRunnerFactory::create( nodeName ).release();
                    if( currRunner->XMLParse( curr ) ){
                        mScenarioRunners.push_back( currRunner  );
                    }
            }
        }
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Unrecognized text string: " << nodeName << " found while parsing the runner-set." << endl;
            success = false;
        }
    }
    return success;
}

/*!
 * \brief Parse a single FileSet element.
 * \details This function parses a single FileSet and adds it to the passed in
 *          ComponentSet's list of FileSets. 
 * \param aNode DOM node corresponding to the current FileSet.
 * \param aCurrComponent The ComponentSet to add this FileSet to.
 * \return Whether the node was successfully parsed.
 */
bool BatchRunner::XMLParseFileSet( const DOMNode* aNode, Component& aCurrComponent ){
    // assume we were passed a valid node.
    assert( aNode );
    
    // Create the new file set and set the name.
    FileSet newFileSet;
    newFileSet.mName = XMLHelper<string>::getAttr( aNode, XMLHelper<void>::name() );

    // get the children of the node.
    DOMNodeList* nodeList = aNode->getChildNodes();

    // loop through the children
    bool success = true;
    for ( unsigned int i = 0; i < nodeList->getLength(); i++ ){
        DOMNode* curr = nodeList->item( i );
        string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );

        if( nodeName == XMLHelper<void>::text() ) {
            continue;
        }
        else if ( nodeName == "Value" ){
            // Create the new File
            File newFile;
            // Get the name of the file.
            newFile.mName = XMLHelper<string>::getAttr( curr, XMLHelper<void>::name() );
            // Get the full path of the file.
            newFile.mPath = XMLHelper<string>::getValue( curr );
            // Add the file to the current new file set.
            newFileSet.mFiles.push_back( newFile );
        }
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Unrecognized text string: " << nodeName << " found while parsing FileSet." << endl;
            success = false;
        }
    }
    // Add the new file set to the current component.
    aCurrComponent.mFileSets.push_back( newFileSet );
    return success;
}

const string& BatchRunner::getXMLNameStatic(){
    static const string XML_NAME = "batch-runner";
    return XML_NAME;
}

Scenario* BatchRunner::getInternalScenario(){
    // The internal scenario runner is not set up until runSingleScenario is
    // called.
    if( mInternalRunner ){
        return mInternalRunner->getInternalScenario();
    }
    return 0;
}

const Scenario* BatchRunner::getInternalScenario() const {
    // The internal scenario runner is not set up until runSingleScenario is
    // called.
    if( mInternalRunner ){
        return mInternalRunner->getInternalScenario();
    }
    return 0;
}

// Implementation for the ParseHelper
bool BatchRunner::ParseHelper::XMLParse( const xercesc::DOMNode* aNode ){
    // assume we were passed a valid node.
    assert( aNode );

    // loop through the children
    bool success = true;
    const string nodeName = XMLHelper<string>::safeTranscode( aNode->getNodeName() );
    if( ScenarioRunnerFactory::isOfType( nodeName ) ){
        if( nodeName == getXMLNameStatic() ){
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::ERROR );
            mainLog << "Batch scenario runners cannot create Batch scenario runners." << endl;
            success = false;
        }
        else {
            mScenarioRunner = ScenarioRunnerFactory::create( nodeName );

            // Allow the IScenarioRunner to parse its own data.
            success &= mScenarioRunner->XMLParse( /*curr*/aNode );
        }
    }
    else {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Unrecognized text string: " << nodeName
                << " found while parsing a scenario runner configuration file." << endl;
        success = false;
    }
    return success;
}

/*!
 * \brief Get an auto-pointer to the parsed scenario runner.
 * \note This method transfers ownership.
 * \pre XMLParse has been called.
 * \return The parsed scenario runner.
 */
auto_ptr<IScenarioRunner>& BatchRunner::ParseHelper::getParsedScenarioRunner(){
    return mScenarioRunner;
}
