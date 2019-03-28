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

#ifndef _MARGINAL_PROFIT_CALCULATOR_H_
#define _MARGINAL_PROFIT_CALCULATOR_H_
#if defined(_MSC_VER)
#pragma once
#endif

/*! 
 * \file marginal_profit_calculator.h
 * \ingroup Objects
 * \brief The MarginalProfitCalculator header file.
 * \author Josh Lurz
 */

class Technology;

/*!
 * \brief Calculates the short term marginal profit for Technologies, normalized
 *        to the non-energy cost.
 * \details Determines the short term marginal profit for Technologies. Marginal
 *          revenue is first calculated, this includes the value of both primary
 *          and secondary good, and emissions subsidies and taxes. Marginal cost
 *          is normalized to the non-energy cost. If non-energy costs are zero,
 *          the marginal cost is returned as an absolute. Marginal costs are
 *          calculated as the variable or fuel costs, and do not include
 *          non-energy costs. This profit rate would be greater than the long
 *          term profit rate given that non-energy costs were greater than zero.
 */
class MarginalProfitCalculator
{
public:
    MarginalProfitCalculator( const Technology* aTechnology );
    
    double calcShortTermMarginalProfit( const std::string& aRegionName,
                                        const std::string& aSectorName,
                                        const int aPeriod ) const;
private:
    //! Technology for which to calculate marginal profits.
    const Technology* mTechnology;
};

#endif // _MARGINAL_PROFIT_CALCULATOR_H_
