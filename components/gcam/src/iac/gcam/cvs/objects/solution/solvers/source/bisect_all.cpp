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
* \file bisect_all.cpp
* \ingroup objects
* \brief BisectAll class source file.
* \author Josh Lurz
*/

#include "util/base/include/definitions.h"
#include <string>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>

#include "solution/solvers/include/solver_component.h"
#include "solution/solvers/include/bisect_all.h"
#include "solution/util/include/calc_counter.h"
#include "marketplace/include/marketplace.h"
#include "containers/include/world.h"
#include "solution/util/include/solution_info.h"
#include "solution/util/include/solution_info_set.h"
#include "solution/util/include/solver_library.h"
#include "util/base/include/util.h"
#include "util/logger/include/ilogger.h"
#include "util/base/include/xml_helper.h"
#include "solution/util/include/solution_info_filter_factory.h"
// TODO: this filter is hard coded here since it is the default, is this ok?
#include "solution/util/include/solvable_solution_info_filter.h"

using namespace std;
using namespace xercesc;

//! Default Constructor. Constructs the base class. 
BisectAll::BisectAll( Marketplace* marketplaceIn, World* worldIn, CalcCounter* calcCounterIn ):SolverComponent( marketplaceIn, worldIn, calcCounterIn ),
mMaxIterations( 30 ),
mDefaultBracketInterval( 0.4 ),
mMaxBracketIterations( 40 )
{
}

//! Init method.
void BisectAll::init() {
    if( !mSolutionInfoFilter.get() ) {
        // note we are hard coding this as the default
        mSolutionInfoFilter.reset( new SolvableSolutionInfoFilter() );
    }
}

//! Get the name of the SolverComponent
const string& BisectAll::getXMLNameStatic() {
    const static string SOLVER_NAME = "bisect-all-solver-component";
    return SOLVER_NAME;
}

//! Get the name of the SolverComponent
const string& BisectAll::getXMLName() const {
    return getXMLNameStatic();
}

bool BisectAll::XMLParse( const DOMNode* aNode ) {
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

/*! \brief Bisection Solution Mechanism (all markets)
* \details This solution mechanism bisects all markets at once. 
* \todo Update this documentation.
* Bisection is always period-formed at least a few times. Bisection stops if the maximum 
* relative ED does not change at a rate larger than BREAK_OUT_THRESHOLD.
* If the maximum relative "ED" is larger than BRACKET_THRESHOLD, then the unsolved markets
* are re-bracketed. Then bisection continues. The bracketing interval is smaller than that
* used initially so as to not perturb trial values too much. If a further re-bracket
* is necessary, the bracketing interval is decreased further. 
*
* Also, the price and demand markets are very prone to move outside their brackets.
* A check for that is period-formed each time and the brackets are adjusted accordingly.
* This check is critical for solution with simultaneously.
*
* Tracking the excess demand is turned on from the logging configuration file.
* \author Sonny Kim, Josh Lurz, Steve Smith
* \warning Unless stated otherwise, ED values are normalized (i.e., that 10 == 10% difference).
* \todo need more general way to reset price and demand market types within bisect
* \todo implement check on price and demand markets within bracket?
* \param aSolutionSet Initial set of objects with information on each market which may be filtered.
* \param aPeriod Model period.
*/
SolverComponent::ReturnCode BisectAll::solve( SolutionInfoSet& aSolutionSet, const int aPeriod ) {
    // If all markets are solved, then return with success code.
    if( aSolutionSet.isAllSolved() ){
        return SolverComponent::SUCCESS;
    }

    // Setup Logging.
    ILogger& solverLog = ILogger::getLogger( "solver_log" );
    solverLog.setLevel( ILogger::NOTICE );
    ILogger& worstMarketLog = ILogger::getLogger( "worst_market_log" );
    worstMarketLog.setLevel( ILogger::NOTICE );
    
    // need to do bracketing first, does this need to be before or after startMethod?
    aSolutionSet.resetBrackets();
    solverLog << "Solution set before Bracket: " << endl << aSolutionSet << endl;
    // Currently attempts to bracket but does not necessarily bracket all markets.
    SolverLibrary::bracket( marketplace, world, mDefaultBracketInterval, mMaxBracketIterations,
                            aSolutionSet, calcCounter, mSolutionInfoFilter.get(), aPeriod );
    
    startMethod();
    ReturnCode code = ORIGINAL_STATE; // code that reports success 1 or failure 0
    
    worstMarketLog << "Policy All, X, XL, XR, ED, EDL, EDR, RED, bracketed, supply, demand" << endl;
    solverLog << "Bisection_all routine starting" << endl; 

    aSolutionSet.updateFromMarkets();
    aSolutionSet.updateSolvable( mSolutionInfoFilter.get() );
    
    // Select the worst market.
    ILogger& singleLog = ILogger::getLogger( "single_market_log" );
    singleLog.setLevel( ILogger::DEBUG );

    unsigned int numIterations = 1; // number of iterations

    do {
        solverLog.setLevel( ILogger::NOTICE );
        solverLog << "BisectionAll " << numIterations << endl;
        aSolutionSet.printMarketInfo( "Bisect All", calcCounter->getPeriodCount(), singleLog );

        // Since bisection is called after bracketing, the current price and ED will be the
        // one of the brackets.
        // Start bisection with mid-point to improve efficiency.
        for ( unsigned int i = 0; i < aSolutionSet.getNumSolvable(); ++i ) {
            SolutionInfo& currSol = aSolutionSet.getSolvable( i );
            // Skip markets that were not bracketed
            if( !currSol.isBracketed() ) {
                continue;
            }
            // If not solved.
            if ( !currSol.isWithinTolerance() ) {
                // Set new trial value to center
                currSol.setPriceToCenter();
            }   
            // price=0 and supply>demand is true only for constraint case.
            // Other markets cannot have supply>demand as price->0.
            // Another condition that should be moved.
            if ( fabs( currSol.getPrice() ) < util::getSmallNumber() && currSol.getED() < 0 ) { 
                currSol.setPrice( 0 ); 
            } 
        }

        aSolutionSet.updateToMarkets();
        marketplace->nullSuppliesAndDemands( aPeriod );
        world->calc( aPeriod );
        aSolutionSet.updateFromMarkets();
        aSolutionSet.updateSolvable( mSolutionInfoFilter.get() );
        // The  price and demand markets are very prone to moving beyond their brackets. 
        // So check and adjust if needed. Lines below check if XL < Demand, or XR > Demand, 
        // and move brackets if necessary. A more general bracket check is below, but this
        // is needed more often and can be done simply.
        aSolutionSet.adjustBrackets();
        // Print solution set information to solver log.
        solverLog << aSolutionSet << endl;

        // Move brackets, both price and ED, after solving mid-point.  This ensures that
        // both price and ED for each bracket is valid and up to date.
        for ( unsigned int i = 0; i < aSolutionSet.getNumSolvable(); ++i ) {
            SolutionInfo& currSol = aSolutionSet.getSolvable( i );
            // If not solved.
            if ( !currSol.isWithinTolerance() && currSol.isBracketed() ) {
                // Move the right price bracket in if Supply > Demand
                if ( currSol.getED() < 0 ) {
                    currSol.moveRightBracketToX();
                }
                // Move the left price bracket in if Demand >= Supply
                else {
                    currSol.moveLeftBracketToX();
                }
            }   
        }

        if( aSolutionSet.getNumSolvable() > 0 ) {
            const SolutionInfo* maxSol = aSolutionSet.getWorstSolutionInfo();
            addIteration( maxSol->getName(), maxSol->getRelativeED() );
            worstMarketLog << "BisectAll-maxRelED: " << *maxSol << endl;
        }
    } // end do loop        
    while ( ++numIterations <= mMaxIterations 
            && !aSolutionSet.isAllSolved() && !areAllBracketsEqual( aSolutionSet ) );

    // Set the return code. 
    code = ( aSolutionSet.isAllSolved() ? SUCCESS : FAILURE_ITER_MAX_REACHED ); // report success, or failure
    
    // Report exit conditions.
    solverLog.setLevel( ILogger::NOTICE );
    if ( numIterations > mMaxIterations ){
        solverLog << "Exiting BisectionAll due to reaching the maximum number of iterations." << endl;
    }
	else if( code != SUCCESS ) {
		solverLog << "Exiting BisectionAll due to reaching the brackets for all markets." << endl;
	}
    else {
        solverLog << "Exiting BisectionAll with model fully solved." << endl;
    }
    return code;
}

/*!
 * \brief Check if all of the solvable solution infos have either left and right brackets
 *        separated by less than the solution tolerance or are solved.
 * \details Since we do not attempt to reset to resize brackets this will provide an early
 *          exit condition.  The left and right brackets be close enough that means
 *          bisection does not have a chance to make progress on that solution info.
 * \return True if the current bracket interval or relative excess demand for all solution
 *         infos are equal with in tolerance otherwise false.
 */
bool BisectAll::areAllBracketsEqual( SolutionInfoSet& aSolutionSet ) const {
	for ( unsigned int i = 0; i < aSolutionSet.getNumSolvable(); ++i ) {
        // TODO: need a get solution tolerance instead of hard coding the .001
		if( !util::isEqual( aSolutionSet.getSolvable( i ).getCurrentBracketInterval(), 0.0, .001 ) 
			&& !aSolutionSet.getSolvable( i ).isWithinTolerance() ) {
			return false;
		}
	}
	return true;
}
