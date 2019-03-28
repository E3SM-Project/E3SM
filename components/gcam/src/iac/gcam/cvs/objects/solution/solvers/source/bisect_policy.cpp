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
* \file bisect_policy.cpp
* \ingroup objects
* \brief BisectPolicy class source file.
* \author Josh Lurz
*/

#include "util/base/include/definitions.h"
#include <string>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>

#include "solution/solvers/include/solver_component.h"
#include "solution/solvers/include/bisect_policy.h"
#include "solution/util/include/calc_counter.h"
#include "marketplace/include/marketplace.h"
#include "marketplace/include/imarket_type.h"
#include "containers/include/world.h"
#include "solution/util/include/solution_info.h"
#include "solution/util/include/solution_info_set.h"
#include "solution/util/include/solver_library.h"
#include "util/logger/include/ilogger.h"
#include "util/base/include/xml_helper.h"
#include "solution/util/include/solution_info_filter_factory.h"
// TODO: this filter is hard coded here since it is the default, is this ok?
#include "solution/util/include/solvable_solution_info_filter.h"

using namespace std;
using namespace xercesc;

//! Default Constructor. Constructs the base class. 
BisectPolicy::BisectPolicy( Marketplace* marketplaceIn, World* worldIn, CalcCounter* calcCounterIn ):SolverComponent( marketplaceIn, worldIn, calcCounterIn ),
mMaxIterations( 30 ),
mDefaultBracketInterval( 0.4 ),
mMaxBracketIterations( 40 )
{
}

//! Init method. Currently does nothing.
void BisectPolicy::init() {
    if( !mSolutionInfoFilter.get() ) {
        // note we are hard coding this as the default
        mSolutionInfoFilter.reset( new SolvableSolutionInfoFilter() );
    }
}

//! Get the name of the SolverComponent
const string& BisectPolicy::getXMLName() const {
    return getXMLNameStatic();
}

//! Get the name of the SolverComponent
const string& BisectPolicy::getXMLNameStatic() {
    const static string SOLVER_NAME = "bisect-policy-solver-component";
    return SOLVER_NAME;
}

bool BisectPolicy::XMLParse( const DOMNode* aNode ) {
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
        else if( nodeName == "max-iterations" ) {
            mMaxIterations = XMLHelper<unsigned int>::getValue( curr );
        }
        else if( nodeName == "bracket-interval" ) {
            mDefaultBracketInterval = XMLHelper<double>::getValue( curr );
        }
        else if( nodeName == "max-bracket-iterations" ) {
            mMaxBracketIterations = XMLHelper<unsigned int>::getValue( curr );
        }
        else if( nodeName == "solution-info-filter" ) {
            mSolutionInfoFilter.reset(
                SolutionInfoFilterFactory::createSolutionInfoFilterFromString( XMLHelper<string>::getValue( curr ) ) );
        }
        else if( SolutionInfoFilterFactory::hasSolutionInfoFilter( nodeName ) ) {
            mSolutionInfoFilter.reset( SolutionInfoFilterFactory::createAndParseSolutionInfoFilter( nodeName, curr ) );
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

/*! \brief Bisection the policy or worst market.
* \details Bisect the the policy or worst market which passes our soultion info filters.
 * \param aSolutionSet Object which contains a set of objects with information on
 *        each market.
* \param aPeriod Model period.
*/
SolverComponent::ReturnCode BisectPolicy::solve( SolutionInfoSet& aSolutionSet, const int aPeriod ) {
    startMethod();

    // If all markets are solved, then return with success code.
    if( aSolutionSet.isAllSolved() ){
        return SolverComponent::SUCCESS;
    }

    // Constants.
    // TODO: do we really want this? it has been taken out of the other SolverComponents
    const static unsigned int MAX_ITER_NO_IMPROVEMENT = 8; // Maximum number of iterations without improvement.
    // Setup logging.
    ILogger& solverLog = ILogger::getLogger( "solver_log" );
    solverLog.setLevel( ILogger::NOTICE );
    solverLog << "BisectPolicy function called." << endl;
    ILogger& worstMarketLog = ILogger::getLogger( "worst_market_log" );
    ILogger& singleLog = ILogger::getLogger( "single_market_log" );
    singleLog.setLevel( ILogger::DEBUG );

    // Make sure we have all updated information.
    aSolutionSet.updateFromMarkets();
    aSolutionSet.updateSolvable( mSolutionInfoFilter.get() );

    // Select the worst market.
    SolutionInfo* worstSol = aSolutionSet.getPolicyOrWorstSolutionInfo();
    worstSol->setBisectedFlag();
    unsigned int numIterations = 0;
    if( !worstSol->isSolved() ){
        SolverLibrary::bracketOne( marketplace, world, mDefaultBracketInterval, mMaxBracketIterations,
                                   aSolutionSet, worstSol, calcCounter, mSolutionInfoFilter.get(), aPeriod );
        if( !worstSol->isSolved() ){
            do {
                aSolutionSet.printMarketInfo( "BisectPolicy on " + worstSol->getName(), calcCounter->getPeriodCount(), singleLog );

                // Move the right price bracket in if Supply > Demand
                if ( worstSol->getED() < 0 ) {
                    worstSol->moveRightBracketToX();
                }
                // Move the left bracket in if Demand >= Supply
                else {
                    worstSol->moveLeftBracketToX();
                }
                // Set new trial value to center
                worstSol->setPriceToCenter();

                // price=0 and supply>demand. only true for constraint case
                // other markets cannot have supply>demand as price->0
                // Another condition that should be moved. 
                if ( fabs( worstSol->getPrice() ) < util::getSmallNumber() && worstSol->getED() < 0 ) { 
                    worstSol->setPrice( 0 ); 
                } 

                aSolutionSet.updateToMarkets();
                marketplace->nullSuppliesAndDemands( aPeriod );

                world->calc( aPeriod );
                aSolutionSet.updateFromMarkets();
                aSolutionSet.updateSolvable( mSolutionInfoFilter.get() );
                addIteration( worstSol->getName(), worstSol->getRelativeED() );
                worstMarketLog << "BisectPolicy-MaxRelED: "  << *worstSol << endl;
            } // end do loop        
            while ( isImproving( MAX_ITER_NO_IMPROVEMENT ) &&
                ( ++numIterations < mMaxIterations ) &&
                !worstSol->isWithinTolerance() );
        }
    }
    // Report results.
    solverLog.setLevel( ILogger::NOTICE );
    if( numIterations >= mMaxIterations ){
        solverLog << "Exiting BisectPolicy due to reaching max iterations." << endl;
    }
    else if( !isImproving( MAX_ITER_NO_IMPROVEMENT ) ){
        solverLog << "Exiting BisectPolicy due to lack of improvement." << endl;
    }
    else {
        solverLog << "Exiting BisectPolicy because chosen market is solved." << endl;
    }
    return aSolutionSet.isAllSolved() ? SUCCESS: FAILURE_ITER_MAX_REACHED; // WRONG ERROR CODE
}
