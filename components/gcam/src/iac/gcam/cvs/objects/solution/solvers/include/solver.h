#ifndef _SOLVER_H_
#define _SOLVER_H_
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
* \file solver.h
* \ingroup Objects
* \brief This is the header file for the Solver class.
* \author Josh Lurz
*/
#include <memory>
#include <string>

#include "util/base/include/iparsable.h"

class CalcCounter;
class Marketplace;
class World;
class SolutionInfoParamParser;
/*!
* \ingroup objects
* \brief Solver is an abstract class which defines a very basic interface to an object which solves the 
* marketplace. A Solver object must define an init() method for setup, and a solve market which attempts
* to solve the model and returns whether or not it did. 
* \author Josh Lurz
*/

class Solver : public IParsable {
public:
    Solver( Marketplace* aMarketplace, World* aWorld ):marketplace( aMarketplace ), world( aWorld ){};
    virtual ~Solver(){};
    virtual void init() = 0;
    virtual bool solve( const int aPeriod, const SolutionInfoParamParser* aSolutionInfoParamParser ) = 0;

protected:
    Marketplace* marketplace; //<! The marketplace to solve. 
    World* world; //!< The world to call calc on.
    //! Weak pointer to the object which tracks the number of calls to World.calc.
    CalcCounter* mCalcCounter;
};

#endif // _SOLVER_H_

