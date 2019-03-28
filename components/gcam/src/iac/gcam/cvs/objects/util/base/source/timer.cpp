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
* \file timer.cpp
* \ingroup Objects
* \brief Timer class source file.
* \author Josh Lurz
*/

#include "util/base/include/definitions.h"
#include <iostream>
#include <ctime>
#include <string>
#include "util/base/include/timer.h"

using namespace std;

//! Constructor
Timer::Timer(){
    mStartTime = 0;
    mStopTime = 0;
}
        
/*! \brief Start the timer. 
* \details This function starts the timer. All times will be relative to this time.
*/     
void Timer::start(){
    mStartTime = clock();
}

/*! \brief Stop the timer.
* \details This functions stops the timer so that the amount of time that has passed
* may be fetched or printed.
*/
void Timer::stop(){
    mStopTime = clock();
}

/*! \brief Get the differential between the start time and stop time.
* \return The difference between the stop and start time.
*/
double Timer::getTimeDifference() const {
    return (double)( mStopTime - mStartTime ) / CLOCKS_PER_SEC;
}
/*! \brief Print the stored time.
* \details This function prints the time between the last call to stop() and the time
* start() was called.
* \param aOut The output stream to print to.
* \param aTitle The label to print in front of the time. Defaults to 'Time: '
*/
void Timer::print( std::ostream& aOut, const string& aLabel ) const {
    aOut << aLabel << " " << getTimeDifference() << " seconds. " << endl; 
}
