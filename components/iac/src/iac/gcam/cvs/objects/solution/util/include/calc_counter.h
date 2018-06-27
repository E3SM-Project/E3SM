#ifndef _CALC_COUNTER_H_
#define _CALC_COUNTER_H_
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
* \file calc_counter.h
* \ingroup Solution
* \brief The header file for the CalcCounter class.
* \author Josh Lurz
*/

#include <string>
#include <iosfwd>
#include <map>

/*!
* \ingroup Solution
* \brief A class used to count iterations of world.calc.
* \details This class tracks calls to world.calc and keeps updated counts of the total 
* number of times world.calc has been called for the model, the total times for the period, 
* and the number of times in the current period by solution mechanism. The world has a pointer
* to this object, so it is updated automatically. Solution mechanisms need to notify it when they
* switch solution methods and when a new period is started. They can then use its accessor functions 
* at any time and are guaranteed to have updated values. 
* \author Josh Lurz
*/
class CalcCounter {
    //! Function which allows the use of the << operator on the CalcCounter object.
     friend std::ostream& operator<<( std::ostream& os, const CalcCounter& calcCounter ){
        calcCounter.print( os );
        return os;
    }
public:
    CalcCounter();
    int getTotalCount() const;
    int getPeriodCount() const;
    int getMethodCount( const std::string methodName ) const;
    void incrementCount( const double additional = 1 );
    void setCurrentMethod( const std::string methodName );
    void startNewPeriod();
private:
    std::string currMethodName;
    std::map<std::string, double> methodCounts;
    double totalCount;
    double periodCount;
    static int convertToInt( double );
    void print( std::ostream& out ) const;
};

#endif // _CALC_COUNTER_H_
