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
* \file log_newton_raphson_sd.cpp
* \ingroup objects
* \brief LogNewtonRaphsonSaveDeriv class source file.
* \author Josh Lurz
*/

#include "util/base/include/definitions.h"
#include <string>

#include "solution/solvers/include/solver_component.h"
#include "solution/solvers/include/log_newton_raphson_sd.h"
#include "solution/util/include/calc_counter.h"
#include "marketplace/include/marketplace.h"
#include "containers/include/world.h"
#include "solution/util/include/solution_info_set.h"
#include "solution/util/include/solution_info.h"
#include "solution/util/include/solver_library.h"
#include "util/base/include/configuration.h"
#include "util/base/include/util.h"
#include "util/logger/include/ilogger.h"

using namespace std;

//! Default Constructor. Need to call constructor of class next up in hierarchy. Constructs the base class. 
LogNewtonRaphsonSaveDeriv::LogNewtonRaphsonSaveDeriv( Marketplace* aMarketplaceIn, World* aWorld,
                                                      CalcCounter* aCalcCounter ):
LogNewtonRaphson( aMarketplaceIn, aWorld, aCalcCounter ),
mSavedMatrixSize( 0 ),
mPermSave( 0 ),
mDerivativesCalculated( false )
{
}

//! Default Destructor. Currently does nothing.
LogNewtonRaphsonSaveDeriv::~LogNewtonRaphsonSaveDeriv(){
}

//! Init method.  
void LogNewtonRaphsonSaveDeriv::init() {
    LogNewtonRaphson::init();
    mDerivativesCalculated = false;
}

//! Get the name of the SolverComponent
const string& LogNewtonRaphsonSaveDeriv::getXMLNameStatic() {
    const static string SOLVER_NAME = "log-newton-raphson-save-deriv-solver-component";
    return SOLVER_NAME;
}

//! Get the name of the SolverComponent
const string& LogNewtonRaphsonSaveDeriv::getXMLName() const {
    return getXMLNameStatic();
}

/*!
 * \brief Calculate derivatives and invert them so that they may be used to calculate new prices.
 * \details Only calculate the derivative if we don't have one that is still valid otherwise just
 *          reuse the one we have.  We will check the mDerivativesCalculated flag to determine if
 *          we have a valid derivative which gets set to false when the number of solvables have
 *          changed or are done trying to solve regardless if it was successful or not.
 * \param aSolutionSet The set of solution infos we are solving.
 * \param JF A factorized derivative matrix which can be used to calculate new prices or garbage
 *           if the return code is not SUCCESS.
 * \param aPermMatrix A matrix to keep track of row operations done while factorizing the derivative
 *                    matrix.  This will be necessary when doing LU back substitution.
 * \param aPeriod The current model period.
 * \see LogNewtonRaphson::calculateDerivatives
 */
SolverComponent::ReturnCode LogNewtonRaphsonSaveDeriv::calculateDerivatives( SolutionInfoSet& aSolutionSet, Matrix& JF,
                                                                             PermutationMatrix& aPermMatrix, int aPeriod )
{
    SolverComponent::ReturnCode code = SUCCESS;
    // Only calculated once.
    if ( !mDerivativesCalculated ) {
        code = LogNewtonRaphson::calculateDerivatives( aSolutionSet, JF, aPermMatrix, aPeriod );
        if( code == SUCCESS ) {
            mDerivativesCalculated = true;
            // Save matricies
            mSavedMatrixSize = aSolutionSet.getNumSolvable();
            JFSave.resize( mSavedMatrixSize, mSavedMatrixSize );
            mPermSave.resize( mSavedMatrixSize );
            JFSave.assign( JF );
            mPermSave.assign( aPermMatrix );
        }
    }
    else {
        // reuse derivatives
        ILogger& solverLog = ILogger::getLogger( "solver_log" );
        solverLog.setLevel( ILogger::NOTICE );
        solverLog << "Using cached derivatives" << endl;
        if ( aSolutionSet.getNumSolvable() != mSavedMatrixSize ) {
            solverLog.setLevel( ILogger::ERROR );
            solverLog << "Matrix sizes changed " << aSolutionSet.getNumSolvable() << ", "<< mSavedMatrixSize << endl;
            
            // force the derivatives to be recalculated next time
            mDerivativesCalculated = false;
            code = FAILURE_SOLUTION_SIZE_CHANGED;
        }
        else {
            JF.assign( JFSave );
            aPermMatrix.assign( mPermSave );
        }
    }
    
    return code;
}

/*!
 * \brief Resets the derivatatives so that the next iteration they will have to be recalculated.
 * \details Derivatives should be reset usually after both a failed or successful attempt to solve
 *          since in both cases the next iteration would likely be better off starting with a fresh
 *          derivative.
 */
void LogNewtonRaphsonSaveDeriv::resetDerivatives() {
    mDerivativesCalculated = false;
}
