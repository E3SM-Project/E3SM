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
 * \file user_configurable_solver.cpp
 * \ingroup objects
 * \brief UserConfigurableSolver class source file.
 * \author Pralit Patel
 */

#include "util/base/include/definitions.h"
#include <iostream>
#include <fstream>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>

#include "solution/solvers/include/user_configurable_solver.h"
#include "containers/include/world.h"
#include "solution/solvers/include/solver_component.h"
#include "solution/solvers/include/solver_component_factory.h"
#include "solution/util/include/solution_info_set.h"
#include "util/base/include/configuration.h"
#include "util/logger/include/ilogger.h"
#include "solution/util/include/calc_counter.h"
#include "util/base/include/xml_helper.h"

using namespace std;
using namespace xercesc;

// typedefs
typedef vector<SolverComponent*>::iterator SolverComponentIterator;
typedef vector<SolverComponent*>::const_iterator CSolverComponentIterator;

//! Constructor
UserConfigurableSolver::UserConfigurableSolver( Marketplace* aMarketplace, World* aWorld ):Solver( aMarketplace, aWorld ),
mDefaultSolutionTolerance( 0.001),
mDefaultSolutionFloor( 0.0001 ),
mCalibrationTolerance( 0.01 ),
mMaxModelCalcs( 2000 )
{
    // get the calc counter from the world
    mCalcCounter = world->getCalcCounter();
}

//! Destructor
UserConfigurableSolver::~UserConfigurableSolver() {
    for( CSolverComponentIterator it = mSolverComponents.begin(); it != mSolverComponents.end(); ++it ) {
        delete *it;
    }
}

//! Get the solver name.
const string& UserConfigurableSolver::getXMLNameStatic() {
    const static string SOLVER_NAME = "user-configurable-solver";
    return SOLVER_NAME;
}

bool UserConfigurableSolver::XMLParse( const DOMNode* aNode ) {
    // assume we were passed a valid node.
    assert( aNode );
    
    // get the children of the node.
    DOMNodeList* nodeList = aNode->getChildNodes();
    
    // loop through the children
    for ( unsigned int i = 0; i < nodeList->getLength(); ++i ){
        DOMNode* curr = nodeList->item( i );
        string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );
        
        if( nodeName == "#text" ) {
            continue;
        }
        else if( nodeName == "solution-tolerance" ) {
            mDefaultSolutionTolerance = XMLHelper<double>::getValue( curr );
        }
        else if( nodeName == "solution-floor" ) {
            mDefaultSolutionFloor = XMLHelper<double>::getValue( curr );
        }
        else if( nodeName == "calibration-tolerance" ) {
            mCalibrationTolerance = XMLHelper<double>::getValue( curr );
        }
        else if( nodeName == "max-model-calcs" ) {
            mMaxModelCalcs = XMLHelper<int>::getValue( curr );
        }
        else if( SolverComponentFactory::hasSolverComponent( nodeName ) ) {
            SolverComponent* tempSolverComponent = SolverComponentFactory::createAndParseSolverComponent( nodeName,
                                                                                                          marketplace,
                                                                                                          world,
                                                                                                          mCalcCounter,
                                                                                                          curr );
            
            // only add valid solver components
            if( tempSolverComponent ) {
                mSolverComponents.push_back( tempSolverComponent );
            }
        }
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Unrecognized text string: " << nodeName << " found while parsing "
                << getXMLNameStatic() << "." << endl;
        }
    }
    return true;
}

//! Initialize the solver at the beginning of the model.
void UserConfigurableSolver::init() {
    /*!
     * \pre The user must have parsed at least one solver component
     *      to try to solve with.
     */
    assert( mSolverComponents.size() > 0 );
}

/*!
 * \brief Solution method to solve all markets for one period.
 * \details This solver will rely on the user to provide which solver components to use and
 *          in which order to clear all markets.
 * \param aPeriod The period to solve.
 * \param aSolutionInfoParamParser An object that will be used to set solution info specific parameters
 *                                 when they are initialized.
 * \return Whether the markets are all solved.
 * \author Pralit Patel
 */
bool UserConfigurableSolver::solve( const int aPeriod, const SolutionInfoParamParser* aSolutionInfoParamParser ) {
    const Configuration* conf = Configuration::getInstance();
    
    // Open all log files.
    ILogger& mainLog = ILogger::getLogger( "main_log" );
    ILogger& solverLog = ILogger::getLogger( "solver_log" );
    ILogger& singleLog = ILogger::getLogger( "single_market_log" );
    mainLog.setLevel( ILogger::DEBUG );
    solverLog.setLevel( ILogger::NOTICE );
    singleLog.setLevel( ILogger::DEBUG );
    solverLog << endl << "Solution() Begin. Per " << aPeriod << endl;
    
    // Create and initialize the solution set.
    // This will fetch the markets to solve and update the prices, supplies and demands.
    SolutionInfoSet solution_set( marketplace );
    solution_set.init( aPeriod, mDefaultSolutionTolerance, mDefaultSolutionFloor, aSolutionInfoParamParser ); // determines solvable and unsolvable markets
    
    mainLog << "Starting Solution. Solving for " << solution_set.getNumSolvable()
        << " markets." << endl;
    solution_set.printMarketInfo( "Begin Solve", mCalcCounter->getPeriodCount(), singleLog );
    
    // If no markets to solve, break out of solution.
    if( solution_set.getNumSolvable() == 0 ){
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Model solved with last period's prices." << endl;
        return true;
    }
    else if( aPeriod == 0 ){
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Period zero has solvable markets." << endl;
    }
    
    // we must have at least one solver component otherwise we will never make any progress
    if( mSolverComponents.size() == 0 ) {
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "No solver components were parsed for use in " << getXMLNameStatic() << '.' << endl;
        mainLog << "Model could not solve period " << aPeriod << endl;
        return false;
    }
    
    solverLog << "Number of Markets: " << solution_set.getNumSolvable() << endl;
    solverLog << "Solution Information Initialized: Left and Right values are the same." << endl;
    solverLog.setLevel( ILogger::DEBUG );
    solverLog << solution_set << endl;
    
    // Initialize solver components
    for( SolverComponentIterator it = mSolverComponents.begin(); it != mSolverComponents.end(); ++it ) {
        (*it)->init();
    }
    
    // Loop is done at least once.
    do {
        solverLog.setLevel( ILogger::NOTICE );
        solverLog << "Solution() loop. N: " << mCalcCounter->getPeriodCount() << endl;
        solverLog.setLevel( ILogger::DEBUG );
        
        // try each solver component in the order they were read
        for( SolverComponentIterator it = mSolverComponents.begin(); it != mSolverComponents.end(); ++it ) {
            // Note we are not checking the return code here since even if a solver component was able to
            // solve successfully it is not necessarily working on the entire solution set.
            (*it)->solve( solution_set, aPeriod );
        }
        
        // Determine if the model has solved. 
    } while ( !solution_set.isAllSolved()
             && mCalcCounter->getPeriodCount() < mMaxModelCalcs );
    
    if( conf->getBool( "CalibrationActive" )
            && !world->isAllCalibrated( aPeriod, mCalibrationTolerance, true ) ) {
        mainLog.setLevel( ILogger::WARNING );
        solverLog << "Model did not calibrate successfully in period " << aPeriod << endl;
        mainLog << "Model did not calibrate successfully in period " << aPeriod << endl;
    }
    
    // Determine whether the solver was successful at solving the model.
    if( solution_set.isAllSolved() ){
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Model solved normally. Iterations period "<< aPeriod << ": "
        << mCalcCounter->getPeriodCount() << ". Total iterations: "
        << mCalcCounter->getTotalCount() << endl;
        return true;
    }
    
    mainLog.setLevel( ILogger::ERROR );
    mainLog << "Model did not solve within set iteration " << mCalcCounter->getPeriodCount() << endl;
    solverLog << "Printing solution information after failed attempt to solve." << endl;
    solverLog << solution_set << endl;
    
    // Print unsolved markets.
    solution_set.printUnsolved( mainLog );
    
    if( conf->getBool( "debugFindSD" ) ){
        string logName = conf->getFile( "supplyDemandOutputFileName", "supply_demand_curves" );
        ILogger& sdLog = ILogger::getLogger( logName );
        sdLog.setLevel( ILogger::WARNING );
        sdLog << "Supply and demand curves for markets that did not solve in period: " << aPeriod << endl;
        solution_set.findAndPrintSD( world, marketplace, aPeriod, sdLog );
    }
    return false;
}
