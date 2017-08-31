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
* \file solver_component.cpp
* \ingroup objects
* \brief SolverComponent class source file.
* \author Josh Lurz
*/

#include "util/base/include/definitions.h"
#include <memory>
#include <string>
#include <iostream>

#include "solution/solvers/include/solver_component.h"
#include "solution/util/include/calc_counter.h"

using namespace std;

/*! \brief Constructor.
* \details This constructor takes as arguments the marketplace, and world which it will be solving, and a pointer to the CalcCounter
* which tracks calls to world.calc(). It also initializes several variables from values in the Configuration object.
* \param marketplaceIn The marketplace which will be used for solving.
* \param worldIn The world which will be used for solving.
* \param calcCounterIn A pointer to the object which tracks calls to world.calc()
*/
SolverComponent::SolverComponent( Marketplace* marketplaceIn, World* worldIn, CalcCounter* calcCounterIn ): marketplace( marketplaceIn ), world( worldIn ), calcCounter( calcCounterIn ){
}

//! Default Destructor.
SolverComponent::~SolverComponent(){
}

//! Struct constructor
SolverComponent::IterationInfo::IterationInfo( const std::string& aName, const double aRED )
:mName( aName ), mRED( aRED ){}

//! Add a solution iteration to the stack.
void SolverComponent::addIteration( const std::string& aSolName, const double aRED ){
    mPastIters.push_back( IterationInfo( aSolName, aRED ) );
}

//! Check for improvement over the last n iterations
bool SolverComponent::isImproving( const unsigned int aNumIter ) const {
    // Check if there are enough iterations to check.
    if( aNumIter >= mPastIters.size()  ){
        return true;
    }

    // Check if there has been improvement
    double currValue = mPastIters.back().mRED;
    // double prevValue = mPastIters.at( mPastIters.size() - aNumIter - 1 ).mRED;

    // return( ( prevValue - currValue ) / currValue > 0.1 ); // This value isnt right.
    unsigned int numBetter = 0;
    // Check how many of the previous are greater than the number.
    for( unsigned int i = 1; i < aNumIter; ++i ){
        double prevValue = mPastIters.at( mPastIters.size() - i - 1 ).mRED;
        if( ( prevValue - currValue ) / currValue > 0.1 ){
            ++numBetter;
        }
    }
    return( static_cast<double>( numBetter ) / ( aNumIter - 1 ) > 0.25 );
}

void SolverComponent::startMethod(){
    // Set the current calculation method.  
    calcCounter->setCurrentMethod( getXMLName() );
    // Clear the stack.
    mPastIters.clear();
}
