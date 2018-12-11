#ifndef _BISECT_POLICY_NR_SOLVER_H_
#define _BISECT_POLICY_NR_SOLVER_H_
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


#include <memory>
#include "solution/solvers/include/solver.h"

/*! 
* \file bisect_policy_nr_solver.h
* \ingroup Objects
* \brief The BisectPolicyNRSolver class header file.
*
* \author Josh Lurz
*/

class SolverComponent;
class Marketplace;
class World;
class SolutionInfoParamParser;
/*!
* \ingroup Objects
* \brief A class which defines an An instance of the Solver class which uses
*        bisection first and then Newton-Raphson.
* \author Josh Lurz
*/

class BisectPolicyNRSolver: public Solver {
public:
    BisectPolicyNRSolver( Marketplace* marketplaceIn, World* worldIn );
    virtual ~BisectPolicyNRSolver();
    static const std::string& getXMLNameStatic();
    
    // Solver methods
    virtual void init();
    virtual bool solve( const int aPeriod, const SolutionInfoParamParser* aSolutionInfoParamParser );
    
    // IParsable methods
    virtual bool XMLParse( const xercesc::DOMNode* aNode );

private:
    std::auto_ptr<SolverComponent> mLogNewtonRaphson; //!< LogNewtonRaphson solver component.
    std::auto_ptr<SolverComponent> mBisectAll; //!< BisectAll solver component.
    std::auto_ptr<SolverComponent> mBisectOne; //!< BisectOne solver component.
    std::auto_ptr<SolverComponent> mBisectPolicy; //!< BisectPolicy solver component.
    
    //! Default solution tolerance, this value may be overridden by at the SolutionInfo level
    double mDefaultSolutionTolerance;
    
    //! Default solution floor, this value may be overridden by at the SolutionInfo level
    double mDefaultSolutionFloor;
    
    //! Calibration tolerance
    double mCalibrationTolerance;
    
    //! Max total solution iterations
    int mMaxModelCalcs;
};

#endif // _BISECT_POLICY_NR_SOLVER_H_
