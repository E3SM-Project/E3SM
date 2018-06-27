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
 * \file secanter.cpp
 * \ingroup Objects
 * \brief Secanter class source file.
 * \author Pralit Patel
 */

#include "util/base/include/definitions.h"
#include <cassert>
#include <string>
#include <cmath>
#include "util/logger/include/ilogger.h"
#include "target_finder/include/secanter.h"
#include "target_finder/include/itarget.h"
#include "util/base/include/util.h"

using namespace std;

/*!
 * \brief Construct the Secanter.
 * \param aTarget The policy target.
 * \param aTolerance Solution tolerance.
 * \param aInitialPrice The initial guess.
 * \param aInitialValue The solution status of the initial guess.
 * \param aInitialPriceChange A percentage change from the initial price to come
 *                            up with the second initial guess necessary for the
 *                            secant method.  This may be a percentage increase
 *                            or decrease depending on if the first initial guess
 *                            was too high or too low.
 * \param aYear Year to check the solution status in.
 */
Secanter::Secanter( const ITarget* aTarget,
                    const double aTolerance,
                    const double aInitialPrice,
                    const double aInitialValue,
                    const double aInitialPriceChange,
                    const int aYear ) :
mTarget( aTarget ),
mTolerance( aTolerance ),
mYear( aYear ),
mIterations( 0 ),
mCurrentTrial( ITargetSolver::undefined(), ITargetSolver::undefined() ),
mPrevTrial( aInitialPrice == 0 ? 1 : aInitialPrice, aInitialValue )
{
    if( aInitialPrice == 0 ) {
        mCurrentTrial.first = aInitialPriceChange + 1;
    }
    else if( fabs( aInitialValue ) < mTolerance ) {
        // Already at the solution
        mCurrentTrial.first = aInitialPrice;
    }
    else if( aInitialValue > 0 ) {
        mCurrentTrial.first = aInitialPrice * ( aInitialPriceChange + 1 );
    }
    else {
        mCurrentTrial.first = aInitialPrice / ( aInitialPriceChange + 1 );
    }
    ILogger& targetLog = ILogger::getLogger( "target_finder_log" );
    targetLog.setLevel( ILogger::DEBUG );
    targetLog << "Constructing a Secanter. Initial points: (" << mCurrentTrial.first << ", "
        << mCurrentTrial.second << "); (" << mPrevTrial.first << ", " << mPrevTrial.second << ")" << endl;
}

/*! \brief Get the next trial value and check for solution.
 * \details Checks if the PolicyTarget is solved and returns a pair representing
 *          the next trial value and whether the PolicyTarget is solved.
 * \return A pair representing the next trial value and whether the PolicyTarget
 *         is solved.
 */
pair<double, bool> Secanter::getNextValue() {
    // Get the status of the current trial.    
    double targetValue = mTarget->getStatus( mYear );
    
    ILogger& mainLog = ILogger::getLogger( "target_finder_log" );
    mainLog.setLevel( ILogger::WARNING );
    mainLog << "Current trial status is " << targetValue << endl;
    
    SolvedState state = eUnsolved;
    
    if( mIterations == 0 ) {
        // This is just the first arbitrary guess which has not yet been run
        mainLog << "Returning initial guess." << endl;
    }
    else if( fabs( targetValue ) < mTolerance ) {
        //case ITarget::SOLVED:
        state = eSolved;
        mCurrentTrial.second = targetValue;
    }
    else {
        //default:
        mCurrentTrial.second = targetValue;
        double slope = ( mCurrentTrial.second - mPrevTrial.second ) / 
            ( log(mCurrentTrial.first) - log(mPrevTrial.first) );
        double newPrice = log(mCurrentTrial.first) - mCurrentTrial.second / slope;
        newPrice = exp( newPrice );
        mPrevTrial = mCurrentTrial;
        mCurrentTrial.first = newPrice;
        mCurrentTrial.second = ITargetSolver::undefined();
    }
    
    // Increment the number of trials.
    ++mIterations;
    
    printState( state );
    
    assert( util::isValidNumber( mCurrentTrial.first ) );
    assert( mCurrentTrial.first >= 0 );
    return make_pair( mCurrentTrial.first, state != eUnsolved );
}

/*! \brief Print the current state of a bisection iteration.
 * \param aState State enum.
 */
void Secanter::printState( const SolvedState aState ) const {
    ILogger& targetLog = ILogger::getLogger( "target_finder_log" );
    targetLog.setLevel( ILogger::DEBUG );
    switch( aState ){
        case eSolved:
            targetLog << "Found solution. ";
            break;
        case eUnsolved:
            targetLog << "Attempting to solve target. Iteration: "
            << mIterations << " ";
            break;
    }
    
    targetLog << "The current trial is: (" << mCurrentTrial.first << ", "
        << mCurrentTrial.second << "); (" << mPrevTrial.first 
        << ", " << mPrevTrial.second << ")" << endl;
}

/*! \brief Get the current number of iterations performed.
 * \return The current number of iterations performed.
 */
unsigned int Secanter::getIterations() const {
    return mIterations;
}
