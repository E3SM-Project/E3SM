#ifndef _LOG_NEWTON_RAPHSON_H_
#define _LOG_NEWTON_RAPHSON_H_
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
* \file log_newton_raphson.h
* \ingroup objects
* \brief This is the header file for the LogNewtonRaphson solver component class.
*
* \author Josh Lurz
*/
#include <string>
#include <boost/numeric/ublas/matrix.hpp>
// include lu for permutation_matrix
#include <boost/numeric/ublas/lu.hpp>

typedef boost::numeric::ublas::matrix<double> Matrix;
typedef boost::numeric::ublas::permutation_matrix<std::size_t> PermutationMatrix;

class CalcCounter; 
class Marketplace;
class World;
class SolutionInfoSet;
class ISolutionInfoFilter;

/*! 
* \ingroup Objects
* \brief A SolverComponent based on the Newton-Raphson algorithm using logarithmic values.
* \details Newton-Raphson does a better job with markets are very interdependent compared
*          to bisections.  To capture the interactions between markets it is often a good
*          idea to include solution infos that are solved.  Newton-Raphson does rely on
*          calculating derivatives which is very computationally intensive and can also
*          diverge from the solution when it is not close or has to deal with non continuous
*          (supply - demand) curves.
* \author Josh Lurz
*/
class LogNewtonRaphson: public SolverComponent {
public:
    LogNewtonRaphson( Marketplace* aMarketplace, World* aWorld, CalcCounter* aCalcCount );
    virtual ~LogNewtonRaphson();
    static const std::string& getXMLNameStatic();
    
    // SolverComponent methods
    virtual void init();
    virtual ReturnCode solve( SolutionInfoSet& aSolutionSet, const int aPeriod );
    virtual const std::string& getXMLName() const;
    
    // IParsable methods
    virtual bool XMLParse( const xercesc::DOMNode* aNode );
    
protected:
    //! The default amount to adjust prices by when calculation derivatives.
    double mDefaultDeltaPrice;
    
    //! Max iterations for this solver component
    unsigned int mMaxIterations;
    
    //! The default allowed jump in price for newton raphson, this may be overriden by a SolutionInfo
    double mDefaultMaxPriceChange;
    
    //! A filter which will be used to determine which SolutionInfos with solver component
    //! will work on.
    std::auto_ptr<ISolutionInfoFilter> mSolutionInfoFilter;
    
    virtual ReturnCode calculateDerivatives( SolutionInfoSet& aSolutionSet, Matrix& JF, PermutationMatrix& aPermMatrix, int aPeriod );
    
    virtual void resetDerivatives();
};

#endif // _LOG_NEWTON_RAPHSON_H_
