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
* \file scenario.cpp
* \ingroup Objects
* \brief Scenario class source file.
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
#include "util/base/include/model_time.h"
#include "marketplace/include/marketplace.h"
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
#include "util/base/include/auto_file.h"
#include "solution/solvers/include/solver_factory.h"
#include "solution/solvers/include/bisection_nr_solver.h"
#include "solution/util/include/solution_info_param_parser.h"

using namespace std;
using namespace xercesc;
using namespace boost;

extern ofstream outFile;
time_t gGlobalTime;

//! Default constructor
Scenario::Scenario() {
    // Get time and date before model run.
    time( &gGlobalTime );

    marketplace.reset( new Marketplace() );
}

//! Destructor
Scenario::~Scenario() {
    // model time is really a singleton and so don't
    // try to delete it
    modeltime.release();
}

/*! \brief Get the static XML name of the Scenario.
* \return The XML name of the scenario element.
*/
const string& Scenario::getXMLNameStatic(){
    const static string XML_NAME = "scenario";
    return XML_NAME;
}

//! Return a reference to the modeltime->
const Modeltime* Scenario::getModeltime() const {
    return modeltime.get();
}

//! Return a constant reference to the goods and services marketplace.
const Marketplace* Scenario::getMarketplace() const {
    return marketplace.get();
}

//! Return a mutable reference to the goods and services marketplace.
Marketplace* Scenario::getMarketplace() {
    return marketplace.get();
}

//! Return a constant reference to the world object.
const World* Scenario::getWorld() const {
    return world.get();
}

//! Return a mutable reference to the world object.
World* Scenario::getWorld() {
    return world.get();
}

//! Set data members from XML input.
bool Scenario::XMLParse( const DOMNode* node ){
    // assume we were passed a valid node.
    assert( node );

    // set the scenario name.
    name = XMLHelper<string>::getAttr( node, "name" );

    // get the children of the node.
    DOMNodeList* nodeList = node->getChildNodes();

    // loop through the children
    for ( unsigned int i = 0; i < nodeList->getLength(); ++i ){
        DOMNode* curr = nodeList->item( i );
        string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );

        if( nodeName == "#text" ) {
            continue;
        }
        else if ( nodeName == Modeltime::getXMLNameStatic() ){
            if( !modeltime.get() ) {
                modeltime.reset( Modeltime::getInstance() );
                const_cast<Modeltime*>( modeltime.get() )->XMLParse( curr );
            }
            else {
                ILogger& mainLog = ILogger::getLogger( "main_log" );
                mainLog.setLevel( ILogger::WARNING );
                mainLog << "Modeltime can only be parsed once." << endl;
            }
        }
        else if ( nodeName == World::getXMLNameStatic() ){
            parseSingleNode( curr, world, new World );
        }
        else if( nodeName == OutputMetaData::getXMLNameStatic() ) {
            parseSingleNode( curr, mOutputMetaData, new OutputMetaData );
        }
        else if( nodeName == SolutionInfoParamParser::getXMLNameStatic() ) {
            parseSingleNode( curr, mSolutionInfoParamParser, new SolutionInfoParamParser );
        }
        
        /*!
         * \warning Parsing of solution algorithms are a special case.  They must be 
         *          parsed after world and modeltime have been parsed and they will never
         *          be written out in toInputXML.  They are also not re-parsable for
         *          example trying to overwrite a parameter in a Solver using an add on
         *          file will not work, rather it will completely overwrite the entire
         *          solver for that period with the latest parsed solver.
         */
        else if( SolverFactory::hasSolver( nodeName ) ) {
            /*!
             * \pre World has already been created.
             */
            assert( world.get() );
            
            boost::shared_ptr<Solver> retSolver( SolverFactory::createAndParseSolver( nodeName, marketplace.get(),
                                                                     world.get(), curr ) );
            
            // make sure we don't attempt to set an invalid solver
            if( retSolver.get() ) {
                /*!
                 * \pre Modeltime has already been created.
                 */
                assert( modeltime.get() );
                
                // this must be done here rather than relying on XMLHelper since we require a factory
                // to create our object
                const int period = modeltime->getyr_to_per( XMLHelper<int>::getAttr( curr, "year" ) );
                const bool fillOut = XMLHelper<bool>::getAttr( curr, "fillout" );
                // we may need to resize the mSolvers which we could do now that we know we have a
                // modeltime
                if( mSolvers.size() == 0 ) {
                    mSolvers.resize( modeltime->getmaxper() );
                }
                mSolvers[ period ] = retSolver;
                
                // TODO: I think just using the same object rather than copying should suffice here
                for( int fillOutPeriod = period + 1; fillOut && fillOutPeriod < modeltime->getmaxper(); ++fillOutPeriod ) {
                    mSolvers[ fillOutPeriod ] = retSolver;
                }
            }
        }
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Unrecognized text string: " << nodeName << " found while parsing scenario." << endl;
        }
    } // end for loop
    return true;
}

//! Sets the name of the scenario. 
void Scenario::setName(string newName) {
    // Used to override the read-in scenario name.
    name = newName;
}

//! Finish all initializations needed before the model can run.
void Scenario::completeInit() {
    // Make sure that some name is set.
    if( name.empty() ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "No scenario name was set, using default." << endl;
        name = "NoScenarioName";
    }
    
    // if we didn't parse any solution info parameters we should at least
    // create an empty one
    if( !mSolutionInfoParamParser.get() ) {
        mSolutionInfoParamParser.reset( new SolutionInfoParamParser() );
    }

    // Complete the init of the world object.
    if( world.get() ){
        world->completeInit();

        // initialize solvers
        initSolvers();
    }
    else {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::SEVERE );
        mainLog << "No world container was parsed from the input files." << endl;
    }

    // Set the valid period vector to false.
    mIsValidPeriod.clear();
    mIsValidPeriod.resize( modeltime->getmaxper(), false );
}

//! Write object to xml output stream.
void Scenario::toInputXML( ostream& out, Tabs* tabs ) const {
    out.precision( 20 );
    // write heading for XML input file
    out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;

    out << "<" << getXMLNameStatic() << " name=\"" << name << "\" date=\"" 
        << util::XMLCreateDate( gGlobalTime ) << "\">" << endl;
    
    tabs->increaseIndent();

    // write the xml for the class members.
    modeltime->toInputXML( out, tabs );
    if( mOutputMetaData.get() ){
        mOutputMetaData->toInputXML( out, tabs );
    }

    world->toInputXML( out, tabs );
    // finished writing xml for the class members.

    XMLWriteClosingTag( getXMLNameStatic(), out, tabs );
}

/*! \brief Write out object to output stream for debugging.
* \param aXMLDebugFile XML debugging file.
* \param aTabs Tabs container.
*/
void Scenario::toDebugXMLOpen( ostream& aXMLDebugFile, Tabs* aTabs ) const {
    string dateString = util::XMLCreateDate( gGlobalTime );
    aXMLDebugFile << "<" << getXMLNameStatic() << " name=\"" << name << "\" date=\"" << dateString << "\">" << endl;

    aTabs->increaseIndent();
    XMLWriteElement( "Debugging output", "summary", aXMLDebugFile, aTabs );
}

/*! \brief Write out close scenario tag to output stream for debugging.
* \param aXMLDebugFile XML debugging file.
* \param aTabs Tabs container.
*/
void Scenario::toDebugXMLClose( ostream& aXMLDebugFile, Tabs* aTabs ) const {
    XMLWriteClosingTag( getXMLNameStatic(), aXMLDebugFile, aTabs );
}

//! Return scenario name.
const string& Scenario::getName() const {
    return name; 
}

/*! \brief Run the scenario.
* \param aSinglePeriod Single period to run or RUN_ALL_PERIODS if all periods
*        should be run.
* \param aPrintDebugging Whether to print extra debugging files.
* \param aFilenameEnding The string to add to the end of the debug output file
*        for uniqueness.
* \return Whether all model runs solved successfully.
*/
bool Scenario::run( const int aSinglePeriod,
                    const bool aPrintDebugging,
                    const string& aFilenameEnding )
{
    // Avoid accumulating unsolved periods.
    unsolvedPeriods.clear();
    
    // Open the debugging files.
    ofstream XMLDebugFile;
    ofstream SGMDebugFile;
    Tabs tabs;
    if( aPrintDebugging ){
        openDebuggingFiles( XMLDebugFile, SGMDebugFile, &tabs, aFilenameEnding );
    }
    
    // Log that a run is beginning.
    logRunBeginning();

    bool success = true;

    // If the single period is RUN_ALL_PERIODS that means to calculate all periods. Loop over
    // time steps and operate model.
    if( aSinglePeriod == RUN_ALL_PERIODS ){
        for( int per = 0; per < modeltime->getmaxper(); per++ ){
            success &= calculatePeriod( per, XMLDebugFile, SGMDebugFile, &tabs, aPrintDebugging );
        }
    }
    // Check if the single period is invalid.
    else if( aSinglePeriod >= modeltime->getmaxper() ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "Invalid single period " << aSinglePeriod << " passed to run method." << endl;
        success = false;
    } 
    else {
        // Run all periods up to the single period which are invalid.
        for( int per = 0; per < aSinglePeriod; per++ ){
            if( !mIsValidPeriod[ per ] ){
                success &= calculatePeriod( per, XMLDebugFile, SGMDebugFile, &tabs, aPrintDebugging );
            }
        }
        
        // Invalidate the period about to be run and all periods past it.
        for( int per = aSinglePeriod; per < modeltime->getmaxper(); ++per ){
            mIsValidPeriod[ per ] = false;
        }

        // Now run the requested period. Results past this period will no longer
        // be valid. Do not attempt to use them!
        success &= calculatePeriod( aSinglePeriod, XMLDebugFile, SGMDebugFile, &tabs, aPrintDebugging );
    }
    
    // Print any unsolved periods.
    // TODO: This should be added to the db.
    ILogger& mainLog = ILogger::getLogger( "main_log" );
    mainLog.setLevel( ILogger::WARNING );
    
    // Report if all model periods solved correctly.
    if( unsolvedPeriods.empty() ) {
        mainLog << "All model periods solved correctly." << endl;
    }
    else {
        // Otherwise print all model periods which did not solve correctly.
        mainLog << "The following model periods did not solve: ";
        for( vector<int>::const_iterator i = unsolvedPeriods.begin(); i != unsolvedPeriods.end(); i++ ) {
            mainLog << *i << ", ";
        }
        mainLog << endl;
    }

    // Run the climate model.
    //jtworld->runClimateModel();

    // Close the debugging files.
    if( aPrintDebugging ){
        closeDebuggingFiles( XMLDebugFile, SGMDebugFile, &tabs );
    }
    
    // Log that the run has finished.
    logRunEnding();
    
    return success;
}

/*! \brief Calculate a single period.
* \param aPeriod Period to calculate.
* \param aXMLDebugFile XML debugging file.
* \param aSGMDebugFile SGM debugging file.
* \param aTabs Tabs formatting object.
* \param aPrintDebugging Whether to print debugging information.
* \return Whether the period was calculated successfully.
*/
bool Scenario::calculatePeriod( const int aPeriod,
                                ostream& aXMLDebugFile,
                                ostream& aSGMDebugFile,
                                Tabs* aTabs,
                                bool aPrintDebugging )
{
    logPeriodBeginning( aPeriod );

    // If this is period 0 initialize market price.
    if( aPeriod == 0 ){
        marketplace->initPrices(); // initialize prices
    }

    // Run the iteration of the model.
    marketplace->nullSuppliesAndDemands( aPeriod ); // initialize market demand to null
    marketplace->init_to_last( aPeriod ); // initialize to last period's info
    world->initCalc( aPeriod ); // call to initialize anything that won't change during calc
    
    // SGM Period 0 needs to clear out the supplies and demands put in by initCalc.
    if( aPeriod == 0 ){
        marketplace->nullSuppliesAndDemands( aPeriod );
    }

    world->calc( aPeriod ); // call to calculate initial supply and demand

    bool success = solve( aPeriod ); // solution uses Bisect and NR routine to clear markets

    world->postCalc( aPeriod );
    
    // Output metadata is not required to exist.
    const list<string>& primaryFuelList = mOutputMetaData.get() ? 
                                          mOutputMetaData->getPrimaryFuelList() 
                                          : list<string>();

    world->updateSummary( primaryFuelList, aPeriod ); // call to update summaries for reporting

    // Mark that the period is now valid.
    mIsValidPeriod[ aPeriod ] = true;
    logPeriodEnding( aPeriod );

    // Write out the results for debugging.
    if( aPrintDebugging ){
        writeDebuggingFiles( aXMLDebugFile, aSGMDebugFile, aTabs, aPeriod );
    }
    return success;
}

/*! \brief Perform any logging which should occur when a period begins.
* \param aPeriod Model period.
*/
void Scenario::logPeriodBeginning( const int aPeriod ) const {
    ILogger& calibrationLog = ILogger::getLogger( "calibration_log" );
    calibrationLog.setLevel( ILogger::DEBUG );
    calibrationLog << "Period " << aPeriod <<": "<< modeltime->getper_to_yr( aPeriod ) << endl << endl;

    ILogger& worstMarketLog = ILogger::getLogger( "worst_market_log" );
    worstMarketLog.setLevel( ILogger::DEBUG );
    worstMarketLog << "Period " << aPeriod <<": "<< modeltime->getper_to_yr( aPeriod ) << endl;

    ILogger& mainLog = ILogger::getLogger( "main_log" );
    mainLog.setLevel( ILogger::NOTICE );
    mainLog << "Period " << aPeriod <<": "<< modeltime->getper_to_yr( aPeriod ) << endl;
}

/*! \brief Perform any logging which should occur when a period ends.
* \param aPeriod Model period.
*/
void Scenario::logPeriodEnding( const int aPeriod ) const {
    ILogger& mainLog = ILogger::getLogger( "main_log" );
    mainLog.setLevel( ILogger::NOTICE );
    mainLog << endl;
}

/*! \brief Perform any logging that should occur when the scenario run begins. */
void Scenario::logRunBeginning() const {
    ILogger& mainLog = ILogger::getLogger( "main_log" );
    mainLog.setLevel( ILogger::NOTICE );
    mainLog << "Model run beginning." << endl;
}

/*! \brief Perform any logging that should occur when the scenario run ends.*/
void Scenario::logRunEnding() const {
    ILogger& mainLog = ILogger::getLogger( "main_log" );
    mainLog.setLevel( ILogger::NOTICE );
    mainLog << "Model run completed." << endl;
}

/*! \brief Open debugging files at the beginning of a scenario run.
* \param aXMLDebugFile XML debugging file.
* \param aSGMDebugFile SGM debugging file.
* \param aTabs Tabs formatting object.
* \param aFileNameEnding String to append to file names.
*/
void Scenario::openDebuggingFiles( ofstream& aXMLDebugFile,
                                   ofstream& aSGMDebugFile,
                                   Tabs* aTabs,
                                   const string& aFileNameEnding ) const
{
    // Open the XML debugging file.
    openDebugXMLFile( aXMLDebugFile, aTabs, aFileNameEnding );

    // Write the opening XML tags.
    toDebugXMLOpen( aXMLDebugFile, aTabs );

    // Sgm output file for debugging.
    const string sgmDebugName = Configuration::getInstance()->getFile( "ObjectSGMFileName", "ObjectSGMout.csv" );
    aSGMDebugFile.open( sgmDebugName.c_str(), std::ios::out );
    util::checkIsOpen( aSGMDebugFile, sgmDebugName );
}

/*! \brief Write to the debugging files for a given period.
* \param aXMLDebugFile XML debugging file.
* \param aSGMDebugFile SGM debugging file.
* \param aTabs Tabs formatting object.
* \param aPeriod Model period.
*/
void Scenario::writeDebuggingFiles( ostream& aXMLDebugFile,
                                    ostream& aSGMDebugFile,
                                    Tabs* aTabs,
                                    const int aPeriod ) const
{
    modeltime->toDebugXML( aPeriod, aXMLDebugFile, aTabs );
    world->toDebugXML( aPeriod, aXMLDebugFile, aTabs );
    csvSGMOutputFile( aSGMDebugFile, aPeriod );
}

/*! \brief Close the debugging files.
* \param aXMLDebugFile XML debugging file.
* \param aSGMDebugFile SGM debugging file.
* \param aTabs Tabs formatting object.
*/
void Scenario::closeDebuggingFiles( ofstream& aXMLDebugFile, ofstream& aSGMDebugFile, Tabs* aTabs ) const {
    // Close the xml debugging tag.
    toDebugXMLClose( aXMLDebugFile, aTabs );

    // Close the files.
    aXMLDebugFile.close();
    aSGMDebugFile.close();
}

/*! \brief Update a visitor for the Scenario.
* \param aVisitor Visitor to update.
* \param aPeriod Period to update.
*/
void Scenario::accept( IVisitor* aVisitor, const int aPeriod ) const {
    aVisitor->startVisitScenario( this, aPeriod );
    // Update the meta-data.
    if( mOutputMetaData.get() ){
        mOutputMetaData->accept( aVisitor, aPeriod );
    }
    // Update the world.
    if( world.get() ){
        world->accept( aVisitor, aPeriod );
    }
    aVisitor->endVisitScenario( this, aPeriod );
}

/*! \brief A function which writes output to an XML file so that it can be read by the XML database.
*/
void Scenario::printOutputXML() const {
#if( __USE_XML_DB__ )	
    // Create a graph printer.
    XMLDBOutputter xmlDBOutputter;
    
    // Update the output container with information from the model.
    // -1 flags to update the output container for all periods at once.
    accept( &xmlDBOutputter, -1 );
    
    
    // Print the output.
    xmlDBOutputter.finish();
#endif
}

/*! \brief A function which print dependency graphs showing fuel usage by
*          sector. 
* \details This function creates a filename and stream for printing the graph
*          data in the dot graphing language. The filename is created from the
*          dependencyGraphName configuration attribute concatenated with the
*          period. The function then calls the World::printDependencyGraphs
*          function to perform the printing. Once the data is printed, dot must
*          be called to create the actual graph as follows: dot -Tpng
*          depGraphs_8.dot -o graphs.png where depGraphs_8.dot is the file
*          created by this function and graphs.png is the file you want to
*          create. The output format can be changed, see the dot documentation
*          for further information.
* \param aPeriod The period to print graphs for.
*/
void Scenario::printGraphs( const int aPeriod ) const {
    // Determine which region to print. Default to the US.
    const string regionToGraph = Configuration::getInstance()->getString( "region-to-graph", "USA" );
    
    // Create a unique filename for the period.
    const string fileName = Configuration::getInstance()->getFile( "dependencyGraphName", "graph" ) 
                            + "_" + util::toString( aPeriod ) + ".dot";
    
    // Open the file. It will automatically close.
    AutoOutputFile graphStream( fileName );
    
    // Create a graph printer.
    GraphPrinter graphPrinter( regionToGraph, *graphStream );
    
    // Update the graph printer with information from the model.
    accept( &graphPrinter, aPeriod );
    
    // Print the graph.
    graphPrinter.finish();
}

void Scenario::printLandAllocatorGraph( const int aPeriod, const bool aPrintValues ) const {
    // Determine which region to print.  Default to the US.
    const string regionToGraph = Configuration::getInstance()->getString( "region-to-graph", "USA" );
    
    //Create a unique filename for the period.
    const string laFileName = 
        Configuration::getInstance()->getFile( "landAllocatorGraphName", "LandAllocatorGraph" )
        + "_" + util::toString( aPeriod ) + ".dot";

    // Open the file.  It will automatically close.
    AutoOutputFile landAllocatorStream( laFileName );

    // Create the land allocator printer.
    LandAllocatorPrinter landAllocatorPrinter( regionToGraph, *landAllocatorStream,
                                               aPrintValues, true );

    // Update the land allocator printer with information from the model.
    accept( &landAllocatorPrinter, aPeriod );

    // Print the graph.
    landAllocatorPrinter.finish();
}

/*! \brief Set a tax into all regions.
* \details TODO
* \param aTax Tax to set.
*/
void Scenario::setTax( const GHGPolicy* aTax ){
    world->setTax( aTax );
}

/*! \brief Get the climate model.
* \return The climate model.
*/
const IClimateModel* Scenario::getClimateModel() const {
    return world->getClimateModel();
}

/*! \brief A function to generate a series of ghg emissions quantity curves
*          based on an already performed model run.
* \details This function used the information stored in it to create a series of
*          curves, one for each region, with each datapoint containing a time
*          period and an amount of gas emissions.
* \note The user is responsible for deallocating the memory in the returned
*       Curves.
* \author Josh Lurz
* \param ghgName The name of the ghg to create a set of curves for.
* \return A vector of Curve objects representing ghg emissions quantity by time
*         period by region.
* \todo Use a visitor to remove this method.
*/
const map<const string, const Curve*> Scenario::getEmissionsQuantityCurves( const string& ghgName ) const {
    /*! \pre The run has been completed. */
    return world->getEmissionsQuantityCurves( ghgName );
}

/*! \brief A function to generate a series of ghg emissions price curves based
*          on an already performed model run.
* \details This function used the information stored in it to create a series of
*          curves, one for each period, with each datapoint containing a time
*          period and the price gas emissions. 
* \note The user is responsible for deallocating the memory in the returned
*       Curves.
* \author Josh Lurz
* \param ghgName The name of the ghg to create a set of curves for.
* \return A vector of Curve objects representing the price of ghg emissions by
*         time period by Region.
* \todo Use a visitor to remove this method.
*/
const map<const string,const Curve*> Scenario::getEmissionsPriceCurves( const string& ghgName ) const {
    /*! \pre The run has been completed. */
    return world->getEmissionsPriceCurves( ghgName );
}

/*! \brief Solve the marketplace using the Solver for a given period. 
* \details The solve method calls the solve method of the instance of the Solver
*          object that was created in the constructor. This method then checks
*          for any errors that occurred while solving and reports the errors if
*          it is the last period. 
* \return Whether the model period solved successfully.
* \param period Period of the model to solve.
*/

bool Scenario::solve( const int period ){
    /*! \pre The solver must be instantiated. */
    assert( mSolvers[ period ].get() );

    // Solve the marketplace. If the return code is false than the model did not
    // solve for the period. Add the period to the scenario list of unsolved
    // periods. 
    const bool success = mSolvers[ period ]->solve( period, mSolutionInfoParamParser.get() );
    if( !success ) {
        unsolvedPeriods.push_back( period );
    }
    
    return success;
}

//! Output Scenario members to a CSV file.
// I don't really like this function being hard-coded to an output file, but its very hard-coded.
void Scenario::writeOutputFiles() const {
    // main output file for sgm, general results.
    {
        AutoOutputFile sgmGenFile( "ObjectSGMGenFileName", "ObjectSGMGen.csv" );
        // SGM csv general output, writes for all periods.
        csvSGMGenFile( *sgmGenFile );
    }
    
    // Print out dependency graphs.
    const Configuration* conf = Configuration::getInstance();
    bool printValues = conf->getBool("PrintValuesOnGraphs");
    if( conf->getBool( "PrintDependencyGraphs" ) ) {
        for( int period = 0; period  < getModeltime()->getmaxper(); ++period  ){
            printGraphs( period );
            // We only need to print a graph for each period if we are printing
            // values on the graphs.  Otherwise they are all the same.
            if( period == 0 || printValues ){
                printLandAllocatorGraph( period, printValues );
            }
        }
    }

    // Open the output file.
    const string outFileName = conf->getFile( "outFileName", "outfile.csv" );
    outFile.open( outFileName.c_str(), ios::out );
    util::checkIsOpen( outFile, outFileName ); 
    outFile.precision( 20 );
    // Write results to the output file.
    // Minicam style output.
    outFile << "Region,Sector,Subsector,Technology,Variable,Units,";
    
    for ( int t = 0; t < modeltime->getmaxper(); t++ ) { 
        outFile << modeltime->getper_to_yr( t ) <<",";
    }
    outFile << "Date,Notes" << endl;

    // Write global market info to file
    marketplace->csvOutputFile( "global" );

    // Write world and regional info
    world->csvOutputFile();
}

//! Output Scenario members to the database.
void Scenario::dbOutput() const {
    // Output metadata is not required to exist.
    const list<string>& primaryFuelList = mOutputMetaData.get() ? 
                                          mOutputMetaData->getPrimaryFuelList() 
                                          : list<string>();
    world->dbOutput( primaryFuelList );
    marketplace->dbOutput();
}

/*! \brief Open the debugging XML file with the correct name and check for any
*          errors.
* \param aXMLDebugFile XML debugging file.
* \param aTabs Tabs formatting container.
* \param aFileNameEnding String to append of the name of the filename.
*/
void Scenario::openDebugXMLFile( ofstream& aXMLDebugFile, Tabs* aTabs, const string& aFileNameEnding ) const {
    // Need to insert the filename ending before the file type.
    const Configuration* conf = Configuration::getInstance();
    string debugFileName = conf->getFile( "xmlDebugFileName", "debug.xml" );
    size_t dotPos = debugFileName.find_last_of( "." );
    debugFileName = debugFileName.insert( dotPos, aFileNameEnding );
    aXMLDebugFile.open( debugFileName.c_str(), ios::out );
    util::checkIsOpen( aXMLDebugFile, debugFileName );
    
    // Write in the main log where the debugging information is.
    ILogger& mainLog = ILogger::getLogger( "main_log" );
    mainLog.setLevel( ILogger::DEBUG );
    mainLog << "Debugging information for this run in: " << debugFileName << endl;
}

/*! \brief Write SGM results to csv text file.
* \param aSGMDebugFile SGM debugging file.
* \param aPeriod Model period for which to print debugging information.
*/
void Scenario::csvSGMOutputFile( ostream& aSGMDebugFile, const int aPeriod ) const {
    aSGMDebugFile <<  "**********************" << endl;
    aSGMDebugFile <<  "RESULTS FOR PERIOD:  " << aPeriod << endl;
    aSGMDebugFile <<  "**********************" << endl << endl;
    marketplace->csvSGMOutputFile( aSGMDebugFile, aPeriod );
    world->csvSGMOutputFile( aSGMDebugFile, aPeriod );
}

/*! \brief Write SGM general results for all periods to csv text file.
*/
void Scenario::csvSGMGenFile( ostream& aFile ) const {
    // Write out the file header.
    aFile << "SGM General Output " << endl;
    aFile << "Date & Time: ";
    util::printTime( gGlobalTime, aFile );
    aFile << endl;
    world->csvSGMGenFile( aFile );
}

/*!
 * \brief A convience method to initialize solvers for all periods.
 * \details First look into the configuration file to see if the user
 *          specified a solver-config file.  If so then parse that file
 *          overwritting any previously created solvers.  Then we check
 *          each period to make sure we have a solver for that period, if
 *          not we will set it to the default solver.  The default solver
 *          is currently BisectionNRSolver.  Finally we call init for
 *          each solver.
 */
void Scenario::initSolvers() {
    // check the config file for a solver config file
    const string solverConfigFile = Configuration::getInstance()->getFile( "solver-config", "" );
    
    // parse the solver config if the user specified one
    if( solverConfigFile != "" ) {
        XMLHelper<void>::parseXML( solverConfigFile, this );
    }
    
    // we may need to resize the mSolvers which we could do now that we know we have a
    // modeltime
    if( mSolvers.size() == 0 ) {
        mSolvers.resize( modeltime->getmaxper() );
    }
    
    // unfortunately we have to access this solver directly
    boost::shared_ptr<Solver> defaultSolver( new BisectionNRSolver( marketplace.get(), world.get() ) );
    for(vector<boost::shared_ptr<Solver> >::iterator solverIt = mSolvers.begin(); solverIt != mSolvers.end(); ++solverIt ) {
        if( !(*solverIt).get() ) {
            (*solverIt) = defaultSolver;
        }
        // Complete the init of the solution object.
        (*solverIt)->init();
    }
}

/*!
 * \brief Get the periods that did not solve in the last call to run.
 * \return A vector of model periods that did not solve.
 */
const vector<int>& Scenario::getUnsolvedPeriods() const {
    return unsolvedPeriods;
}
