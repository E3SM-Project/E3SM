#ifndef _OUTPUT_SHARE_LEVELIZED_COST_CALCULATOR_H_
#define _OUTPUT_SHARE_LEVELIZED_COST_CALCULATOR_H_
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
 * \file output_share_levelized_cost_calculator.h
 * \ingroup Objects
 * \brief The OutputShareLevelizedCostCalculator header file.
 * \author Josh Lurz
 */

#include "investment/include/iexpected_profit_calculator.h"
class IInvestable;

/*! 
* \ingroup Objects
* \brief An object responsible for calculating the levelized cost of a
*        technology or sector.
* \details When calculating sector level levelized cost an Investable will
*          use it's levelized cost times it's share of annual investemnet
*          at that particular level.
* \warning This calculator should only be used after investment has been
*          distributed.
* \author Pralit Patel
*/
class OutputShareLevelizedCostCalculator: public IExpectedProfitRateCalculator
{
public:
    OutputShareLevelizedCostCalculator();
    static const std::string& getXMLNameStatic();
    void toInputXML( std::ostream& aOut, Tabs* aTabs ) const;
    void toDebugXML( const int aPeriod, std::ostream& aOut, Tabs* aTabs ) const;

    double calcSectorExpectedProfitRate( const std::vector<IInvestable*>& aInvestables,
                                         const NationalAccount& aNationalAccount,
                                         const std::string& aRegionName,
                                         const std::string& aGoodName,
                                         const double aInvestmentLogitExp,
                                         const bool aIsShareCalc,
                                         const bool aIsDistributing,
                                         const int aPeriod ) const;

    double calcTechnologyExpectedProfitRate( const ProductionFunctionInfo& aTechProdFuncInfo,
                                             const NationalAccount& aNationalAccount,
                                             const std::string& aRegionName,
                                             const std::string& aSectorName,
                                             const double aDelayedInvestmentTime,
                                             const int aLifetime,
                                             const int aTimeStep,
                                             const int aPeriod ) const;
};

#endif // _OUTPUT_SHARE_LEVELIZED_COST_CALCULATOR_H_
