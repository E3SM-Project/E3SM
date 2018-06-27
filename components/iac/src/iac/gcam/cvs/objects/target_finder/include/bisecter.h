#ifndef _BISECTER_H_
#define _BISECTER_H_
#if defined(_MSC_VER)
#pragma once
#endif

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
* \file bisecter.h
* \ingroup Objects
* \brief The Bisecter class header file.
* \author Josh Lurz
*/

#include "target_finder/include/itarget_solver.h"

class ITarget;

/*! \brief Object which performs bisection on a given target until it
*          reaches a tolerance.
*/
class Bisecter : public ITargetSolver {
public:
    Bisecter( const ITarget* aTarget,
              const double aTolerance,
              const double aMinimum,
              const double aMaximum,
              const double aInitialValue,
              const double aMultiple,
              const int aYear );

    // ITargetSolver methods
    std::pair<double, bool> getNextValue();

    unsigned int getIterations() const;
private:
    /*
    * \brief An enumeration of all states the bisection algorithm may be in at
    *        the end of an iteration.
    */
    enum SolvedState {
        //! Solution has been found.
        eSolved,

        //! Solution has not been found.
        eUnsolved,

        //! The difference between the upper and lower bounds is less than the
        //! tolerance.
        eEmptyBracket,

        //! The difference between the upper bound and the current trial is less
        //! than the tolerance.
        eUpperBoundReached,

        //! The difference between the lower bound and the current trial is less
        //! than the tolerance.
        eLowerBoundReached
    };

    //! The target.
    const ITarget* mTarget;

    //! The tolerance of the target.
    const double mTolerance;

    //! The minimum of the search.
    double mMinimum;

    //! The maximum of the search.
    double mMaximum;

    //! The initial trial value.
    double mInitialGuess;
    
    //! The adjustment to make during initital set of iterations until maximum value is found.
    double mMultiple;

    //! The current lower bound.
    double mLowerBound;

    //! The current upper bound.
    double mUpperBound;

    //! The current trial value.
    double mCurrentTrial;
    
    //! The current number of trial values returned.
    unsigned int mIterations;

    //! Year in which the bisecter is operating.
    unsigned int mYear;

    void printState( const SolvedState aState ) const;
};

#endif // _BISECTER_H_
