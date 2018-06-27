#ifndef _USER_CONFIGURABLE_SOLVER_H_
#define _USER_CONFIGURABLE_SOLVER_H_
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


#include <vector>
#include <memory>
#include "solution/solvers/include/solver.h"

/*! 
 * \file user_configurable_solver.h
 * \ingroup Objects
 * \brief The UserConfigurableSolver class header file. 
 * \author Pralit Patel
 */

class SolverComponent;
class Marketplace;
class World;
class SolutionInfoParamParser;

/*!
 * \ingroup Objects
 * \brief A class which defines an An instance of the Solver class which will follow
 *        a solving procedure completely determined by the solver xml config file.
 * \details The user will be responsible for setting up the solver components to use
 *          and providing the ordering of operations through parsing of the xml solver
 *          config file.
 *          <b>XML specification for UserConfigurableSolver</b>
 *          - XML name: \c user-configurable-solver
 *          - Contained by: scenario
 *          - Parsing inherited from class: None.
 *          - Elements:
 *              - \c solution-tolerance double UserConfigurableSolver::mDefaultSolutionTolerance
 *                      The default solution tolerance to use when solving.
 *              - \c solution-floor double UserConfigurableSolver::mDefaultSolutionFloor
 *                      The default solution floor to use when solving.
 *              - \c calibration-tolerance double UserConfigurableSolver::mCalibrationTolerance
 *                      The tolerance used for checking calibrated values when calibrating.
 *              - \c max-model-calcs int UserConfigurableSolver::mMaxModelCalcs
 *                      The maximum total number of iterations to try to find a solution.
 *              - \c (any SolverComponent) vector<SolverComponent*> UserConfigurableSolver::mSolverComponents
 *                      Can be any solver component contained in SolverComponentFactory, each one
 *                      being added in order to the list of solver components to use.
 *
 * \author Pralit Patel
 */
class UserConfigurableSolver: public Solver {
public:
    UserConfigurableSolver( Marketplace* aMarketplace, World* aWorld );
    virtual ~UserConfigurableSolver();
    static const std::string& getXMLNameStatic();
    
    // Solver methods
    virtual void init();
    virtual bool solve( const int aPeriod, const SolutionInfoParamParser* aSolutionInfoParamParser );
    
    // IParsable methods
    virtual bool XMLParse( const xercesc::DOMNode* aNode );
    
private:
    //! In order list of solver components to use when trying to solve.
    std::vector<SolverComponent*> mSolverComponents;
    
    //! Default solution tolerance, this value may be overridden at the SolutionInfo level
    double mDefaultSolutionTolerance;
    
    //! Default solution floor, this value may be overridden at the SolutionInfo level
    double mDefaultSolutionFloor;
    
    //! Calibration tolerance
    double mCalibrationTolerance;
    
    //! Max total solution iterations
    int mMaxModelCalcs;
};

#endif // _USER_CONFIGURABLE_SOLVER_H_

