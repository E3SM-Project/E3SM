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
* \file calc_counter.cpp
* \ingroup Solution
* \brief CalcCounter class source file.
* \author Josh Lurz
*/

#include "util/base/include/definitions.h"
#include <cassert>
#include <iostream>
#include <cmath>

#include "util/base/include/util.h"
#include "solution/util/include/calc_counter.h"

using namespace std;

//! Constructor
CalcCounter::CalcCounter() {
    totalCount = 0;
    periodCount = 0;
}

/* \brief Return the total number of iterations of world.calc called so far for all periods.
* \return Integer value of the number of calls of world.calc called so far for all periods.
*/
int CalcCounter::getTotalCount() const {
    return convertToInt( totalCount );
}

/* \brief Return the total number of iterations of world.calc called so far for the current periods.
* \return Integer value of the number of calls of world.calc called so far for the current periods.
*/
int CalcCounter::getPeriodCount() const {
    return convertToInt( periodCount );
}

/*! \brief Return the number of iterations of world.calc called so far by a given solution method in the current period.
* \param methodName The name of the method for which to get the number of world.calc calls.
* \return The number of times the given solution method has called world.calc in the current period.
*/
int CalcCounter::getMethodCount( const string methodName ) const {
    return convertToInt( util::searchForValue( methodCounts, methodName ) );
}

/*!\brief Increment the world.calc count by a given amount, 1 by default.
* \details This method increments the total count, period count and count for the current method
* by the amount passed as an argument. 
* \param additional Amount to increment the counts by, 1 is the default.
*/
void CalcCounter::incrementCount( const double additional ){
    totalCount += additional;
    periodCount += additional;
    methodCounts[ currMethodName ] += additional;
}

/*! \brief Set the name of the method currently being used to solve.
* \param methodName The name of the method now being used to solve.
*/
void CalcCounter::setCurrentMethod( const string methodName ){
    currMethodName = methodName;
}

/*! \brief Start a new period. 
* \details Starts a new period by resetting the period based counters.
*/
void CalcCounter::startNewPeriod(){
    periodCount = 0;
    methodCounts.clear();
}

/*! \brief Utility helper function to convert to an integer from the ceiling of a double.
* \param value Double value to convert.
* \return Integer with the value of the ceiling of the double passed in.
*/
int CalcCounter::convertToInt( const double value ) {
    return static_cast<int>( ceil( value ) );
}

/*! \brief Print out the information contained within the CalcCounter.
* \param out output stream to print to.
*/
void CalcCounter::print( ostream& out ) const {
    out << "Period Count: " << periodCount << endl;
    out << "Total Count: " << totalCount << endl;
    out << "Per Method Period Counts: " << endl;

    typedef map<string, double>::const_iterator MethodCountIterator;

    for( MethodCountIterator iter = methodCounts.begin(); iter != methodCounts.end(); ++iter ){
        out << "Method: " << iter->first << " Count: " << iter->second << endl;
    }
    out << endl;
}


