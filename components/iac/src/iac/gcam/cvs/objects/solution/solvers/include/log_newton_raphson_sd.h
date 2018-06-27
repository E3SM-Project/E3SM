#ifndef _LOG_NEWTON_RAPHSON_SD_H_
#define _LOG_NEWTON_RAPHSON_SD_H_
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
* \file log_newton_raphson_sd.h
* \ingroup objects
* \brief This is the header file for the LogNewtonRaphsonSaveDeriv solver component class.
*
* \author Josh Lurz
*/

#include "solution/solvers/include/log_newton_raphson.h"

/*! 
* \ingroup Objects
* \brief A SolverComponent based on the Newton-Raphson algorithm using logarithmic values
*        which reuses derivatives.
* \details This solver component will calculate the derivative once and will continue to
*          reuse them until either the number of solvables have changed or resetDerivatives
*          is called.
* \author Josh Lurz
*/
class LogNewtonRaphsonSaveDeriv: public LogNewtonRaphson {
public:
    LogNewtonRaphsonSaveDeriv( Marketplace* aMarketplaceIn, World* aWorld, CalcCounter* aCalcCounter );
    virtual ~LogNewtonRaphsonSaveDeriv();
    static const std::string& getXMLNameStatic();
    void init();
    const std::string& getXMLName() const;
    
protected:
    virtual ReturnCode calculateDerivatives( SolutionInfoSet& aSolutionSet, Matrix& JF, PermutationMatrix& aPermMatrix, int aPeriod );
    
    virtual void resetDerivatives();

    //! Flag used to determine if the saved derivatives are suitable for use.
    bool mDerivativesCalculated;
    
    //! The saved lu-factorized derivative from the last time they were calculated.
    Matrix JFSave;
    
    //! The saved permutation matrix from the lu-factorization from the last time
    //! they were calculated.
    PermutationMatrix mPermSave;
    
    //! Keeps track of the JFSave size so that we can detect if the number of solvables
    //! have changed in which case we could no longer reuse the saved derivative.
    unsigned int mSavedMatrixSize;
};

#endif // _LOG_NEWTON_RAPHSON_SD_H_
