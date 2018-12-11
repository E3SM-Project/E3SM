#ifndef _RATE_LOGIT_DISTRIBUTOR_H_
#define _RATE_LOGIT_DISTRIBUTOR_H_
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
 * \file rate_logit_distributor.h
 * \ingroup Objects
 * \brief The RateLogitDistributor class header file.
 * \author Josh Lurz
 */

class IInvestable;
class IExpectedProfitRateCalc;

#include "investment/include/idistributor.h"
#include <vector>

/*! 
 * \ingroup Objects
 * \brief This object distributes investment using a logit determined by
 *        expected profits.
 * \details The base functionality of the RateLogitDistributor is to calculate a
 *          set of shares which are used to distribute investment. These shares
 *          are based on profits calculated by an IExpectedProfitCalculator
 *          object which is passed into the functions on this interface. Once a
 *          set of shares is calculated, they can be used to either distribute
 *          an amount of investment or calculate an average capital to output
 *          coefficient. Calculating a capital to output coefficient must use
 *          the same shares as distributing investment since the capital to
 *          output ratio is on a per unit basis, which is only valid given the
 *          same distribution of investment.
 * \author Josh Lurz
 */
class RateLogitDistributor: public IDistributor
{
public:
    RateLogitDistributor();

    double distribute( const IExpectedProfitRateCalculator* IExpectedProfitRateCalculator,
                       std::vector<IInvestable*>& aInvestables,
                       NationalAccount& aNationalAccount,
                       const double aInvestmentExp,
                       const std::string& aRegionName,
                       const std::string& aSectorName,
                       const double aAmount,
                       const int aPeriod ) const;
    
    double calcCapitalOutputRatio( const std::vector<IInvestable*>& aInvestables,
                                   const IExpectedProfitRateCalculator* aRateCalc,  
                                   const NationalAccount& aNationalAccount,
                                   const std::string& aRegionName,
                                   const std::string& aSectorName,
                                   const int aPeriod ) const;
private:

    const std::vector<double> calcInvestmentShares( const std::vector<IInvestable*>& aInvestables,
                                                    const IExpectedProfitRateCalculator* aRateCalc,
                                                    const NationalAccount& aNationalAccount,
                                                    const double aInvestmentExp,
                                                    const std::string& aRegionName,
                                                    const std::string& aSectorName,
                                                    const int aPeriod ) const;
};


#endif // _RATE_LOGIT_DISTRIBUTOR_H_
