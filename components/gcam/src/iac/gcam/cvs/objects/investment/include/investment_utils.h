#ifndef _INVESTMENT_UTILS_H_
#define _INVESTMENT_UTILS_H_
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
 * \file investment_utils.h
 * \ingroup Objects
 * \brief The InvestmentUtils class header file.
 * \author Josh Lurz
 */

#include <string>
#include <vector>

class IInvestable;
class BaseTechnology;

/*! 
 * \ingroup Objects
 * \brief A set of utility functions for calculating investment.
 * \author Josh Lurz
 */
class InvestmentUtils
{
public:
    typedef std::vector<IInvestable*>::const_iterator CInvestableIterator;
    typedef std::vector<IInvestable*>::iterator InvestableIterator;
    
    static double interpolateAndSumFlows( const double aPrevInvestment, const double aCurrInvestment,
                                          const int aIntervalYears );

    static double calcBaseCapital( const std::string& aRegionName, 
                                   const double aPrevInvestment,
                                   const double aAggInvFrac,
                                   const int aPeriod );

    static double sumInvestment( const std::vector<IInvestable*>& aInvestables,
                                 const int aPeriod );

    static double normalizeShares( std::vector<double>& aShares );

    static double sumFixedInvestment( const std::vector<IInvestable*>& aInvestables,
                                      const int aPeriod );

    static std::vector<IInvestable*> getTechInvestables( const std::vector<BaseTechnology*>& aAllTechs,
        const int aPeriod );
    
    /*!
     * \brief Template function which converts a vector of IInvestable subtypes
     *          to a vector of Investables.
     * \details Convert the vector of IInvestable subtypes to a vector of
     *          IInvestable objects. This is legal given that the template class
     *          passed in is a subtype of IInvestable. This will fail at
     *          compilation time if this is not true. This must be done
     *          explicitly because even if class B inherits from class A,
     *          vector<B> does not inherit from vector<A>.
     * \param aSubtypeVector A vector of subtypes of IInvestable to upcast.
     * \return A vector of IInvestable pointers corresponding to the subtype
     *         pointer vector passed in.
     * \todo Replace this function with an enumerator interface.
     */
    template <class T> 
        static std::vector<IInvestable*> convertToInvestables( const std::vector<T>& aSubtypeVector ){
            std::vector<IInvestable*> investables( aSubtypeVector.size() );
            std::copy( aSubtypeVector.begin(), aSubtypeVector.end(), investables.begin() );
            return investables;
        }
};

#endif // _INVESTMENT_UTILS_H_
