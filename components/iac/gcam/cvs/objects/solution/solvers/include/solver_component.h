#ifndef _SOLVER_COMPONENT_H_
#define _SOLVER_COMPONENT_H_
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
* \file solver_component.h
* \ingroup objects
* \brief This is the header file for the SolverComponent class.
*
* \author Josh Lurz
*/
#include <string>
#include <memory>
#include <vector>

#include "util/base/include/iparsable.h"

class CalcCounter; 
class Marketplace;
class SolutionInfoSet;
class World;

/*! \brief An abstract class defining an interface to an independent component
*          of a Solver.
* \details A SolverComponent is an independent object which takes a Marketplace
*           and attempts to clear all markets to a given relative excess demand
*           tolerance within a given number of iterations. SolverComponents may
*           use static functions within the SolverLibrary class, they may not
*           however use other SolverComponents. This is the role of the Solver
*           object. The main method within a SolverComponent is the solve
*           method, which attempts to clear the markets.
* \author Josh Lurz
*/
class SolverComponent : public IParsable {
public:
    //! Return code of the solve method. 
    enum ReturnCode {
        ORIGINAL_STATE,
        SUCCESS,
        FAILURE_ITER_MAX_REACHED,
        FAILURE_WRONG_DIRECTION,
        FAILURE_SOLUTION_SIZE_CHANGED,
        FAILURE_SINGULAR_MATRIX
    };
   SolverComponent( Marketplace* marketplaceIn, World* worldIn, CalcCounter* calcCounterIn );
   virtual ~SolverComponent();
   
   virtual void init() = 0;
   
   virtual ReturnCode solve( SolutionInfoSet& aSolutionSet, const int aPeriod ) = 0;

   virtual const std::string& getXMLName() const = 0;

protected:
   Marketplace* marketplace; //<! The marketplace to solve. 
   World* world; //<! World to call calc on.
   CalcCounter* calcCounter; //<! Tracks the number of calls to world.calc
   
   //! A structure used to track maximum relative excess demand over many
   //! iterations.
   struct IterationInfo {
   public:
       //! Name of the SolverComponent active in the iteration.
       std::string mName;

       //! Maximum relative excess demand in the iteration.
       double mRED;

       IterationInfo( const std::string& aName = std::string(), const double aRED = 0 );
   };

   std::vector<IterationInfo> mPastIters;
   void addIteration( const std::string& aSolName, const double aRED );
   bool isImproving( const unsigned int aNumIter ) const;
   void startMethod();
};

#endif // _SOLVER_COMPONENT_H_
