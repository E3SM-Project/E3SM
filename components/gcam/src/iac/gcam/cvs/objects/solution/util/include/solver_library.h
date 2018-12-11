#ifndef _SOLVER_LIBRARY_H_
#define _SOLVER_LIBRARY_H_
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
* \file solver_library.h
* \ingroup Objects
* \brief A file containing the header for the static SolverLibrary class which
*        contains helper methods used by SolverComponents.
* \author Josh Lurz, Sonny Kim
*/
#include <vector>
#include <map>
#include <string>
#include <functional>
#include <boost/numeric/ublas/matrix.hpp>
// include lu for permutation_matrix
#include <boost/numeric/ublas/lu.hpp>

typedef boost::numeric::ublas::matrix<double> Matrix;
typedef boost::numeric::ublas::permutation_matrix<std::size_t> PermutationMatrix;

class Marketplace;
class World;
class SolutionInfo;
class SolutionInfoSet;
class CalcCounter;
class ISolutionInfoFilter;
namespace objects {
    class Atom;
}
/*!
* \ingroup Objects
* \brief A class with all static functions which are used by the
*        SolverComponents classes and contains common functionality. 
* \author Josh Lurz, Sonny Kim
*/

class SolverLibrary {
public:
    // Some of these still might go.
   static double getRelativeED( const double excessDemand, const double demand, const double excessDemandFloor );

   static bool isWithinTolerance( const double excessDemand, const double demand, const double solutionTolerance,
       const double excessDemandSolutionFloor );

   static void derivatives( Marketplace* marketplace, World* world, SolutionInfoSet& solutionVector,
       const double aDeltaPrice, const int per );

   static bool luFactorizeMatrix( Matrix& aInputMatrix, PermutationMatrix& aPermMatrix );

   static bool bracketOne( Marketplace* aMarketplace, World* aWorld, const double aDefaultBracketInterval,
                           const unsigned int aMaxIterations, SolutionInfoSet& aSolSet, SolutionInfo* aSol,
                           CalcCounter* aCalcCounter, const ISolutionInfoFilter* aSolutionInfoFilter, const int aPeriod );

   static void updateMatrices( SolutionInfoSet& sol, Matrix& JFSM, Matrix& JFDM, Matrix& JF );

   static bool calculateNewPricesLogNR( SolutionInfoSet& aSolutionSet, Matrix& JFLUFactorized,
                                        PermutationMatrix& aPermMatrix, const double aDefaultMaxPriceJump ); 

   static bool bracket( Marketplace* aMarketplace, World* aWorld, const double aDefaultBracketInterval,
                        const unsigned int aMaxIterations, SolutionInfoSet& aSolSet, CalcCounter* aCalcCounter,
                        const ISolutionInfoFilter* aSolutionInfoFilter, const int aPeriod );

private:
    typedef std::map<const objects::Atom*, std::vector<double> > RegionalMarketValues;

    //! A simple struct to link Supplies and Demands.
    struct RegionalSDDifferences {
        RegionalMarketValues supplies;
        RegionalMarketValues demands;
    };

    //! A function object to compare to values and see if they are approximately equal. 
    struct ApproxEqual : public std::unary_function<double, bool> {
        const double compareValue; //!< A value to compare the argument value against.
        const double tolerance; //!< The tolerance within which to return that the values are equal.
        ApproxEqual( double compareValueIn, double toleranceIn ):
        compareValue( compareValueIn ), tolerance( toleranceIn ){}
        bool operator()( const double value ){
            return fabs( value - compareValue ) < tolerance;
        }
    };

    static bool doRegionalValuesSum( const RegionalMarketValues& regionalValues,
        const std::vector<double>& worldTotals, const SolutionInfoSet& aSolutionSet );

    static const RegionalSDDifferences calcRegionalSDDifferences( Marketplace* marketplace, World* world,
        SolutionInfoSet& sol, const int per );

    static std::vector<double> storePrices( const SolutionInfoSet& aSolutionSet );
    static void restorePrices( SolutionInfoSet& aSolutionSet, const std::vector<double>& aPrices );
};

#endif // _SOLVER_LIBRARY_H_
