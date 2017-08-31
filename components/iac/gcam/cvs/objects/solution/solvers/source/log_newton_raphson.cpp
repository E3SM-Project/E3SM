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
* \file log_newton_raphson.cpp
* \ingroup objects
* \brief LogNewtonRaphson class source file.
* \author Josh Lurz
*/

#include "util/base/include/definitions.h"
#include <string>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>

#include "solution/solvers/include/solver_component.h"
#include "solution/solvers/include/log_newton_raphson.h"
#include "solution/util/include/calc_counter.h"
#include "marketplace/include/marketplace.h"
#include "containers/include/world.h"
#include "solution/util/include/solution_info_set.h"
#include "solution/util/include/solution_info.h"
#include "solution/util/include/solver_library.h"
#include "util/base/include/util.h"
#include "util/logger/include/ilogger.h"
#include "util/base/include/xml_helper.h"
#include "solution/util/include/solution_info_filter_factory.h"
// TODO: this filter is hard coded here since it is the default, is this ok?
#include "solution/util/include/solvable_nr_solution_info_filter.h"

using namespace std;
using namespace xercesc;

//! Default Constructor. Constructs the base class. 
LogNewtonRaphson::LogNewtonRaphson( Marketplace* aMarketplace,
                                    World* aWorld, CalcCounter* aCalcCounter ):
SolverComponent( aMarketplace, aWorld, aCalcCounter ),
mDefaultDeltaPrice( 1e-6 ),
mMaxIterations( 4 ),
mDefaultMaxPriceChange( 1.2 )
{
}

//! Destructor
LogNewtonRaphson::~LogNewtonRaphson() {
}

//! Init method.
void LogNewtonRaphson::init() {
    if( !mSolutionInfoFilter.get() ) {
        // note we are hard coding this as the default
        mSolutionInfoFilter.reset( new SolvableNRSolutionInfoFilter() );
    }
}

//! Get the name of the SolverComponent
const string& LogNewtonRaphson::getXMLNameStatic() {
    static const string SOLVER_NAME = "log-newton-raphson-solver-component";
    return SOLVER_NAME;
}

//! Get the name of the SolverComponent
const string& LogNewtonRaphson::getXMLName() const {
    return getXMLNameStatic();
}

bool LogNewtonRaphson::XMLParse( const DOMNode* aNode ) {
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
        else if( nodeName == "delta-price" ) {
            mDefaultDeltaPrice = XMLHelper<double>::getValue( curr );
        }
        else if( nodeName == "max-iterations" ) {
            mMaxIterations = XMLHelper<unsigned int>::getValue( curr );
        }
        else if( nodeName == "max-price-change" ) {
            mDefaultMaxPriceChange = XMLHelper<double>::getValue( curr );
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

/*! \brief Ron's version of the Newton Raphson Solution Mechanism (all markets)
* \details Calculates derivatives and then uses the newton-raphson method to calculate
*          new prices in each iterations.  We use log price changes and limit the max price
*          changes produced to minimize the risk of this algorithm diverging.
* \author Sonny Kim, Josh Lurz, Steve Smith
* \warning Unless stated otherwise, ED values are normalized (i.e., that 10 == 10% difference). -- TODO: Is this still true?
* \param aSolutionSet An initial set of SolutionInfo objects representing all markets which can be filtered.
* \param aPeriod Model period.
* \return A status code to indicate if the algorithm was successful or not.
* \see SolverLibrary::calculateNewPricesLogNR
*/
SolverComponent::ReturnCode LogNewtonRaphson::solve( SolutionInfoSet& aSolutionSet, const int aPeriod ) {
    ReturnCode code = SolverComponent::ORIGINAL_STATE;

    // If all markets are solved, then return with success code.
    if( aSolutionSet.isAllSolved() ){
        return code = SolverComponent::SUCCESS;
    }

    startMethod();
    
    // TODO: is the following necessary
    // Update the solution vector for the correct markets to solve.
    aSolutionSet.updateFromMarkets();
    // Need to update solvable status before starting solution (Ignore return code)
    aSolutionSet.updateSolvable( mSolutionInfoFilter.get() );

    ILogger& solverLog = ILogger::getLogger( "solver_log" );
    solverLog.setLevel( ILogger::NOTICE );
    solverLog << "Log-Newton-Raphson Beginning" << endl;
    ILogger& worstMarketLog = ILogger::getLogger( "worst_market_log" );
    worstMarketLog.setLevel( ILogger::DEBUG );
    ILogger& singleLog = ILogger::getLogger( "single_market_log" );
    singleLog.setLevel( ILogger::DEBUG );

    if( aSolutionSet.getNumSolvable() == 0 ){
        solverLog << "Exiting Newton-Raphson early. No non-singular markets." << endl;
        return SUCCESS; // Need a new code here.
    }

    // Setup the solution matrices.
    size_t currSize = aSolutionSet.getNumSolvable();
    Matrix JF( currSize, currSize );
    PermutationMatrix permMatrix( currSize );

    unsigned int number_of_NR_iteration = 1;
    bool success = true;
    do {
        singleLog.setLevel( ILogger::DEBUG );
        aSolutionSet.printMarketInfo( "Begin logNR", calcCounter->getPeriodCount(), singleLog );

        // Calculate derivatives
        code = calculateDerivatives( aSolutionSet, JF, permMatrix, aPeriod );
        if ( code != SUCCESS ) {
            return code;
        }

        // Calculate new prices
        if( SolverLibrary::calculateNewPricesLogNR( aSolutionSet, JF, permMatrix, mDefaultMaxPriceChange ) ){
            // Call world.calc and update supplies and demands. 
            aSolutionSet.updateToMarkets();
            marketplace->nullSuppliesAndDemands( aPeriod );
            solverLog << "Supplies and demands calculated with new prices." << endl;
            world->calc( aPeriod );
            aSolutionSet.updateFromMarkets();

            // Add to the iteration list.
            SolutionInfo* currWorstSol = aSolutionSet.getWorstSolutionInfo();
            addIteration( currWorstSol->getName(), currWorstSol->getRelativeED() );

            worstMarketLog.setLevel( ILogger::NOTICE );
            worstMarketLog << "NR-maxRelED: " << *currWorstSol << endl;
            solverLog.setLevel( ILogger::DEBUG );
            solverLog << "Solution after " << number_of_NR_iteration << " iterations in NewtonRhapson: " << endl;
            solverLog << aSolutionSet << endl;

            if( aSolutionSet.updateSolvable( mSolutionInfoFilter.get() ) != SolutionInfoSet::UNCHANGED ){
                size_t newSize = aSolutionSet.getNumSolvable();
                JF.resize( newSize, newSize );
                permMatrix.resize( newSize );
            }

            aSolutionSet.printMarketInfo( "NR routine ", calcCounter->getPeriodCount(), singleLog );
        }
        else {
            success = false;
        }
        ++number_of_NR_iteration;
    } // end do loop    
    while ( success && number_of_NR_iteration <= mMaxIterations &&
            !aSolutionSet.isAllSolved() );
    
    // reset the derivatives regardless of if we were successful
    resetDerivatives();

    // Update the return code. 
    code = ( aSolutionSet.isAllSolved() ? SUCCESS : FAILURE_ITER_MAX_REACHED );

    // Print if we exited NR because it had solved all the markets.
    solverLog.setLevel( ILogger::NOTICE );
    if( code == SUCCESS && !aSolutionSet.isAllSolved() ){
        solverLog << "Newton-Raphson solved all markets except at least one with a singularity." << endl;
        aSolutionSet.printUnsolved( solverLog );
    }
    else if( code == SUCCESS && aSolutionSet.isAllSolved() ){
        solverLog << "Newton-Raphson solved all markets successfully." << endl;
    }
    else if( number_of_NR_iteration > mMaxIterations ){
        solverLog << "Exiting Newton-Rhapson without solving all markets. Maximum iteration count exceeded: " << mMaxIterations << endl;
    }
    else {
        solverLog << "Exiting Newton-Rhapson due to lack of improvement. " << endl;
    }
    
    solverLog << endl; // new line
    return code;
}

/*!
 * \brief Calculate derivatives and invert them so that they may be used to calculate new prices.
 * \details We will rely on the solver library to calculate derivatives.  The returned JF will not
 *          actually be inverted rather it will use LU factorization which is computationally less
 *          intensive than calculating inverse.  Note that SolverLibrary::calculateNewPricesLogNR
 *          has also been modified to expect the factorized derivative rather than the inverse.
 * \param aSolutionSet The set of solutions infos which will be included in the derivative calculation.
 * \param JF The LU factorized derivative matrix or garbage if for some reason the calculation failed.
 * \param aPermMatrix A matrix to keep track of row operations done while factorizing the derivative
 *                    matrix.  This will be necessary when doing LU back substitution.
 * \param aPeriod The current model period.
 * \return A return code to indicate if the derivative calculation as well as the factorization was
 *         successful.  It is very important that this value be checked since if it failed JF will
 *         have garbage in it.
 * \see SolverLibrary::derivatives
 */
SolverComponent::ReturnCode LogNewtonRaphson::calculateDerivatives( SolutionInfoSet& aSolutionSet, Matrix& JF, PermutationMatrix& aPermMatrix, int aPeriod ) {
    // Calculate derivatives.
    SolverLibrary::derivatives( marketplace, world, aSolutionSet, mDefaultDeltaPrice, aPeriod );
    ILogger& solverLog = ILogger::getLogger( "solver_log" );
    solverLog.setLevel( ILogger::NOTICE );
    solverLog << "Derivatives calculated" << endl;
    size_t currSize = aSolutionSet.getNumSolvable();
    // we need to keep track of these since if luFactorizeMatrix fails we can reconstruct
    // what JF was
    Matrix JFDM( currSize, currSize );
    Matrix JFSM( currSize, currSize );
    
    // Update the JF, JFDM, and JFSM matrices
    SolverLibrary::updateMatrices( aSolutionSet, JFSM, JFDM, JF );
    if( SolverLibrary::luFactorizeMatrix( JF, aPermMatrix ) ) {
        // print the derivative matrix to help the user understand why it was singular
        solverLog.setLevel( ILogger::ERROR );
        solverLog << "Matrix came back as singular, could not invert." << endl;
        solverLog << ',';
        for( int col = 0; col < aSolutionSet.getNumSolvable(); ++col ) {
            solverLog << aSolutionSet.getSolvable( col ).getName() << ',';
        }
        solverLog << endl;
        for( int row = 0; row < aSolutionSet.getNumSolvable(); ++row ) {
            solverLog << aSolutionSet.getSolvable( row ).getName() << ',';
            for( int col = 0; col < aSolutionSet.getNumSolvable(); ++col ) {
                // we need to recalculate JF since luFactorizeMatrix will not necessarily
                // return it with the rows properly ordered when there is a singularity
                solverLog << JFSM( row , col ) - JFDM( row , col ) << ',';
            }
            solverLog << endl;
        }
        solverLog.setLevel( ILogger::NOTICE );
        return FAILURE_SINGULAR_MATRIX;
    }
    return SUCCESS;
}

/*!
 * \brief Resets any derivatives so that the next iteration of Newton-Raphson will
 *        be forced to calculate a new one.
 * \details Nothing is done for LogNewtonRaphson since a new derivative is calculated
 *          for each iteration of Newton-Raphson regardless.
 */
void LogNewtonRaphson::resetDerivatives() {
    // nothing to do since derivatives are recalculated each time
}
