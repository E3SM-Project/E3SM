#ifndef _CALC_CAPITAL_GOOD_PRICE_VISITOR_H_
#define _CALC_CAPITAL_GOOD_PRICE_VISITOR_H_
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
* \file calc_capital_good_price_visitor.h
* \ingroup Objects
* \brief The CalcCapitalGoodPriceVisitor class header file.
*
* \author Pralit Patel
*/

#include <string>
#include "util/base/include/default_visitor.h"

class InvestConsumer;

/*! 
* \ingroup Objects
* \brief Calculates the price of the capital good
* \details Uses the CES aggregate of the inputs to the investment consumer to calculate
*          the price of the capital good.  This must be calculated and set in the market
*          info before the production sectors operate.
* \author Pralit Patel
*/
class CalcCapitalGoodPriceVisitor : public DefaultVisitor {
public:
    //! Default Constructor
    CalcCapitalGoodPriceVisitor( std::string& aRegionName );

    // IVisitor methods
    virtual void startVisitInvestConsumer( const InvestConsumer* aInvestConsumer, 
        const int aPeriod );
private:
    const std::string mRegionName;
};

#endif // _CALC_CAPITAL_GOOD_PRICE_VISITOR_H_
