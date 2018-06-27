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
* \file bisect_policy_nr_solver.cpp
* \ingroup objects
* \brief BisectPolicyNRSolver class source file.
* \author Josh Lurz
*/

#include "util/base/include/definitions.h"
#include <memory>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>

#include "solution/solvers/include/bisect_policy_nr_solver.h"
#include "containers/include/world.h"
#include "solution/solvers/include/solver_component.h"
#include "solution/solvers/include/solver_component_factory.h"
#include "solution/util/include/solution_info_set.h"
#include "solution/util/include/calc_counter.h"
#include "util/base/include/configuration.h"
#include "util/logger/include/ilogger.h"
#include "util/base/include/xml_helper.h"

// need to include these so that we can get the xml names
#include "solution/solvers/include/log_newton_raphson.h"
#include "solution/solvers/include/bisect_all.h"
#include "solution/solvers/include/bisect_one.h"
#include "solution/solvers/include/bisect_policy.h"

using namespace std;
using namespace xercesc;

//! Constructor
BisectPolicyNRSolver::BisectPolicyNRSolver( Marketplace* marketplaceIn, World* worldIn ):Solver( marketplaceIn, worldIn ),
mDefaultSolutionTolerance( 0.001),
mDefaultSolutionFloor( 0.0001 ),
mCalibrationTolerance( 0.01 ),
mMaxModelCalcs( 2000 )
{
    // Construct components.
    mCalcCounter = world->getCalcCounter();
    mLogNewtonRaphson.reset( SolverComponentFactory::createAndParseSolverComponent( LogNewtonRaphson::getXMLNameStatic(), marketplace,
                                                                                   world, mCalcCounter, 0 ) );
    mBisectAll.reset( SolverComponentFactory::createAndParseSolverComponent( BisectAll::getXMLNameStatic(), marketplace,
                                                                                   world, mCalcCounter, 0 ) );
    mBisectOne.reset( SolverComponentFactory::createAndParseSolverComponent( BisectOne::getXMLNameStatic(), marketplace,
                                                                                   world, mCalcCounter, 0 ) );
    mBisectPolicy.reset( SolverComponentFactory::createAndParseSolverComponent( BisectPolicy::getXMLNameStatic(), marketplace,
                                                                                   world, mCalcCounter, 0 ) );
}

//! Destructor
BisectPolicyNRSolver::~BisectPolicyNRSolver() {
}

//! Get the name of the SolverComponent
const string& BisectPolicyNRSolver::getXMLNameStatic() {
    const static string SOLVER_NAME = "BisectPolicyNRSolver";
    return SOLVER_NAME;
}

bool BisectPolicyNRSolver::XMLParse( const DOMNode* aNode ){
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
        else if( nodeName == mLogNewtonRaphson->getXMLName() ) {
            mLogNewtonRaphson->XMLParse( curr );
        }
        else if( nodeName == mBisectAll->getXMLName() ) {
            mBisectAll->XMLParse( curr );
        }
        else if( nodeName == mBisectOne->getXMLName() ) {
            mBisectOne->XMLParse( curr );
        }
        else if( nodeName == mBisectPolicy->getXMLName() ) {
            mBisectPolicy->XMLParse( curr );
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
void BisectPolicyNRSolver::init() {
}

/*! \brief Solution method to solve all markets for one period.
* \details This is the main solution function called from within the Marketplace. It is called once for each period to clear all 
* markets which should be solved. This solve method first brackets the markets, then uses several iterations of bisection_all to move 
* the prices into the range of the solution, and then uses Newton-Raphson to clear the markets.
* \param aPeriod The period to solve.
* \param aSolutionInfoParamParser An object that will be used to set solution info specific value
*                                 when they are initialized.
* \return Whether the markets all solved.
*/
bool BisectPolicyNRSolver::solve( const int aPeriod, const SolutionInfoParamParser* aSolutionInfoParamParser ) {
    static unsigned int MAX_BISECT_ONE_MARKETS = 2;
    // Create and initialize the SolutionInfoSet. 
    // This will fetch the markets to solve and update the prices, supplies and demands.
    SolutionInfoSet sol( marketplace );
    sol.init( aPeriod, mDefaultSolutionTolerance, mDefaultSolutionFloor, aSolutionInfoParamParser );
    
    ILogger& mainLog = ILogger::getLogger( "main_log" );
    mainLog.setLevel( ILogger::DEBUG );
    mainLog << "Starting Solution. Solving for " << sol.getNumSolvable() << " markets." << endl;

    ILogger& singleLog = ILogger::getLogger( "single_market_log" );
    mainLog.setLevel( ILogger::DEBUG );
    sol.printMarketInfo( "Begin Solve", mCalcCounter->getPeriodCount(), singleLog );
    
    // if no markets to solve, break out of solution.
    if ( sol.getNumSolvable() == 0 ){
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Model solved with last period's prices." << endl;
        return true;
    }
    else if( aPeriod == 0 ){
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Period zero has solvable markets." << endl;
    }

    // Print out extra debugging information.
    ILogger& solverLog = ILogger::getLogger( "solver_log" );
    solverLog.setLevel( ILogger::NOTICE );
    solverLog << "Solution() Begin. Per " << aPeriod << endl;
    solverLog << "Number of Markets: " << sol.getNumSolvable() << endl;
    solverLog << "Solution Information Initialized: Left and Right values are the same." << endl;
    
    solverLog.setLevel( ILogger::DEBUG );
    solverLog << sol << endl;

    // Loop is done at least once.
    do {
        solverLog.setLevel( ILogger::NOTICE );
        solverLog << "Solution() loop. N: " << mCalcCounter->getPeriodCount() << endl;
        solverLog.setLevel( ILogger::DEBUG );
        solverLog << "Solution before BisectPolicy: " << endl;
        solverLog << sol << endl;
        
        // Bisect the policy market or the worst market if the policy market is non-existant.
        sol.unsetBisectedFlag();
        mBisectPolicy->solve( sol, aPeriod );

        solverLog.setLevel( ILogger::DEBUG );
        solverLog << "Solution before NewtonRaphson: " << endl;
        solverLog << sol << endl;

        // Call mLogNewtonRaphson. Ignore return code because it may have skipped singular markets.
        mLogNewtonRaphson->solve( sol, aPeriod );

        solverLog.setLevel( ILogger::DEBUG );
        solverLog << "After NewtonRaphson " << mCalcCounter->getPeriodCount() << endl;
        solverLog << sol << endl;

        if( !sol.isAllSolved() ){
            unsigned int count = 0;
            while( !sol.isAllSolved() && count < MAX_BISECT_ONE_MARKETS && ( sol.hasSingularUnsolved() 
                || sol.getMaxRelativeExcessDemand() > 1 ) ){
                // Try bisecting a single market.
                mBisectOne->solve( sol, aPeriod );
                ++count;
                }
                cout << "Exiting Bisect-One loop. Count: " << count << " HasSingularUnsolved: " << sol.hasSingularUnsolved() << endl;
        }

        solverLog.setLevel( ILogger::DEBUG );
        solverLog << "Solution before NewtonRaphson: " << endl;
        solverLog << sol << endl;
        if( !sol.isAllSolved() ){
            // Call mLogNewtonRaphson. Ignore return code because it may have skipped singular markets.
            mLogNewtonRaphson->solve( sol, aPeriod );
        }
        solverLog.setLevel( ILogger::DEBUG );
        solverLog << "After NewtonRaphson " << mCalcCounter->getPeriodCount() << endl;
        solverLog << sol << endl;
    // Determine if the model has solved. 
    } while ( !sol.isAllSolved() && mCalcCounter->getPeriodCount() < mMaxModelCalcs );
    
    mainLog.setLevel( ILogger::NOTICE );
    if( sol.isAllSolved() ){
        mainLog << "Model solved normally. Iterations period "<< aPeriod << ": " << mCalcCounter->getPeriodCount() << ". Total iterations: "<< mCalcCounter->getTotalCount() << endl;
        return true;
    }
   
    mainLog << "Model did not solve within set iteration " << mCalcCounter->getPeriodCount() << endl;
    solverLog << "Printing solution information after failed attempt to solve." << endl;
    solverLog << sol << endl;

    // Print unsolved markets.
    sol.printUnsolved( mainLog );

    if( Configuration::getInstance()->getBool( "debugFindSD" ) ){
        string logName = Configuration::getInstance()->getFile( "supplyDemandOutputFileName", "supply_demand_curves" );
        ILogger& sdLog = ILogger::getLogger( logName );
        sdLog.setLevel( ILogger::WARNING );
        sdLog << "Supply and demand curves for markets that did not solve in period: " << aPeriod << endl;
        sol.findAndPrintSD( world, marketplace, aPeriod, sdLog );
    }
    return false;
}
